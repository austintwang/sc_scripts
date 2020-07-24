#!/usr/bin/env python3

import os
import sys
import subprocess

def stream(contig, vcf_path, bam_path):
    cmd = ["samtools", "view", "-h", "-L", vcf_path, bam_path, str(contig)]
    subprocess.run(cmd, check=True, shell=True)

if __name__ == '__main__':
    stream(*sys.argv[1:])