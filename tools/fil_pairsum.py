#!/usr/bin/env python3
"""
fil_pairsum.py — sum adjacent time samples of a 2-bit SIGPROC filterbank, losslessly, into 4 bits.

Why: in the CRAFTS 2-bit filterbanks the noise-diode cal toggles every sample (it is mistimed),
putting a strong tone at exactly half the sampling rate (5086 Hz at 98.304 us) in every channel.
PRESTO rfifind then sees a periodic signal in every channel and every interval and zaps
essentially the whole band (1022/1024 channels in a test).  Summing sample pairs removes the
cal (on + off in every output sample) and the zap fraction falls to ~37%.

Two 2-bit values (0..3) sum to 0..6, which fits in 4 bits exactly: no requantisation, no loss,
and the output is the same size as the input (half the samples, twice the bits).

  out[t, c] = in[2t, c] + in[2t+1, c]        tsamp -> 2*tsamp, nbits 2 -> 4
  packing: LSB-first, as in the input (first value in the low bits of each byte)

Usage:
  fil_pairsum.py IN.fil OUT.fil
  fil_pairsum.py --outdir DIR IN1.fil IN2.fil ...      # OUT = DIR/<name>_x2.fil

Every header field is copied; only tsamp and nbits change (and nsamples, if present).
An odd final sample is dropped.  Streams in blocks, so memory use is small.
"""

import argparse
import os
import struct
import sys
import time

import numpy as np

INT_KEYS = {"telescope_id", "machine_id", "data_type", "barycentric", "pulsarcentric", "nbits",
            "nsamples", "nchans", "nifs", "nbeams", "ibeam"}
DBL_KEYS = {"az_start", "za_start", "src_raj", "src_dej", "tstart", "tsamp", "fch1", "foff",
            "refdm", "period"}
STR_KEYS = {"source_name", "rawdatafile"}
BYTE_KEYS = {"signed"}          # 1-byte field written by some FAST tools


def read_header(f):
    def rstr():
        n = struct.unpack("<i", f.read(4))[0]
        return f.read(n).decode()
    if rstr() != "HEADER_START":
        raise ValueError("not a SIGPROC filterbank")
    items = []
    while True:
        k = rstr()
        if k == "HEADER_END":
            return items, f.tell()
        if k in INT_KEYS:
            items.append((k, "i", struct.unpack("<i", f.read(4))[0]))
        elif k in DBL_KEYS:
            items.append((k, "d", struct.unpack("<d", f.read(8))[0]))
        elif k in BYTE_KEYS:
            items.append((k, "b", struct.unpack("<b", f.read(1))[0]))
        elif k in STR_KEYS:
            items.append((k, "s", rstr()))
        else:
            raise ValueError(f"unknown header key {k!r}")


def write_header(o, items):
    def wstr(s):
        b = s.encode()
        o.write(struct.pack("<i", len(b)) + b)
    wstr("HEADER_START")
    for k, t, v in items:
        wstr(k)
        if t == "i":
            o.write(struct.pack("<i", v))
        elif t == "d":
            o.write(struct.pack("<d", v))
        elif t == "b":
            o.write(struct.pack("<b", v))
        else:
            wstr(v)
    wstr("HEADER_END")


def pairsum(src, dst, block_samples=65536):
    t0 = time.time()
    with open(src, "rb") as f:
        items, off = read_header(f)
        h = {k: v for k, _, v in items}
        if h.get("nbits") != 2:
            raise ValueError(f"{src}: nbits={h.get('nbits')}, expected 2")
        if h.get("nifs", 1) != 1:
            raise ValueError(f"{src}: nifs={h.get('nifs')}, expected 1")
        nc = h["nchans"]
        if nc % 2:
            raise ValueError("nchans must be even for 4-bit packing")
        bytes_per_samp = nc // 4
        total = (os.path.getsize(src) - off) // bytes_per_samp
        nout = total // 2

        out_items = []
        for k, t, v in items:
            if k == "tsamp":
                v = v * 2
            elif k == "nbits":
                v = 4
            elif k == "nsamples":
                v = nout
            out_items.append((k, t, v))

        tmp = dst + ".part"
        with open(tmp, "wb") as o:
            write_header(o, out_items)
            block_samples -= block_samples % 2
            done = 0
            while done < nout * 2:
                n = min(block_samples, nout * 2 - done)
                raw = np.frombuffer(f.read(n * bytes_per_samp), np.uint8)
                if raw.size != n * bytes_per_samp:
                    raise IOError("short read")
                # unpack 2-bit LSB-first -> (n, nc) values 0..3
                v = np.empty((raw.size, 4), np.uint8)
                for i in range(4):
                    v[:, i] = (raw >> (2 * i)) & 3
                v = v.reshape(n, nc)
                s = v[0::2] + v[1::2]                      # (n/2, nc) values 0..6
                s = s.reshape(-1, 2)
                o.write((s[:, 0] | (s[:, 1] << 4)).astype(np.uint8).tobytes())   # LSB-first nibbles
                done += n
        os.replace(tmp, dst)
    mb = os.path.getsize(dst) / 1e6
    print(f"{os.path.basename(src)} -> {dst}  {total} -> {nout} samples, "
          f"tsamp {h['tsamp'] * 1e6:.3f} -> {2 * h['tsamp'] * 1e6:.3f} us, 2 -> 4 bit, "
          f"{mb:.0f} MB in {time.time() - t0:.1f} s", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--outdir", help="write DIR/<name>_x2.fil for each input")
    a = ap.parse_args()
    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)
        jobs = [(s, os.path.join(a.outdir, os.path.basename(s)[:-4] + "_x2.fil")) for s in a.files]
    elif len(a.files) == 2 and not os.path.exists(a.files[1]) or (len(a.files) == 2 and a.files[1].endswith("_x2.fil")):
        jobs = [(a.files[0], a.files[1])]
    else:
        ap.error("give IN OUT, or --outdir DIR IN...")
    for s, d in jobs:
        if os.path.abspath(s) == os.path.abspath(d):
            sys.exit("refusing to overwrite the input")
        pairsum(s, d)


if __name__ == "__main__":
    main()
