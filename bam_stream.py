#!/usr/bin/env python3

import os
import sys
import subprocess

def stream(contig, bed_path, bam_path):
    print("stream start", file=sys.stderr)
    cmd = ["samtools", "view", "-h", "-L", bed_path, bam_path, str(contig)]
    subprocess.run(cmd, check=True)

if __name__ == '__main__':
    stream(*sys.argv[1:])