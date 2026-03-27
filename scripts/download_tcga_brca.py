#!/usr/bin/env python3
from ccspr.datasets.tcga_brca import download_tcga_brca

if __name__ == "__main__":
    out = download_tcga_brca("data/tcga_brca")
    print(out)
