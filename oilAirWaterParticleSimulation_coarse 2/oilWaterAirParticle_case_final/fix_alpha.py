#!/usr/bin/env python3
import re, glob

pattern = re.compile(rb'[0-9]\.[0-9]+e-[0-9]{3}')  # matches e-100 through e-999

files = (
    glob.glob('processor*/[0-9]*/alpha.air')
    + glob.glob('processor*/[0-9]*/alpha.air_0')
    + glob.glob('processor*/[0-9]*/alpha.oil')
    + glob.glob('processor*/[0-9]*/alpha.oil_0')
)
for f in files:
    with open(f, 'rb') as fh:
        data = fh.read()
    new_data, n = pattern.subn(b'0', data)
    if n:
        with open(f, 'wb') as fh:
            fh.write(new_data)
        print(f"{f}: fixed {n} bad values")
