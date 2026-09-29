#!/usr/bin/env python3

import sys

ref = sys.argv[1]
chr = sys.argv[2]

on = False

with open(ref, "r") as f:
    for line in f:
        if line.startswith(">"):
            if line[1:].strip().split()[0] == chr:
               on = True
            else:
               on = False

        if on:
            print(line, end='', file=sys.stdout)
