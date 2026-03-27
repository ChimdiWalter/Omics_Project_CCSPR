#!/usr/bin/env python3
from ccspr.datasets.tcga_luad import download_tcga_luad

if __name__ == "__main__":
    out = download_tcga_luad("data/tcga_luad")
    print(out)
