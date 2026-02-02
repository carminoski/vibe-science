"""
Download supplementary files from GSE157153.
"""

import urllib.request
import gzip
import shutil
from pathlib import Path

OUT_DIR = Path(__file__).parent

# Supplementary files to download
files = [
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157153/suppl/GSE157153_clinical_metadata.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157153/suppl/GSE157153_endometriosis_associated.txt.gz",
]

for url in files:
    filename = url.split("/")[-1]
    outpath = OUT_DIR / filename
    outpath_unzipped = OUT_DIR / filename.replace(".gz", "")

    print(f"Downloading {filename}...")
    try:
        urllib.request.urlretrieve(url, outpath)
        print(f"  Saved: {outpath}")

        # Decompress
        print(f"  Decompressing...")
        with gzip.open(outpath, 'rb') as f_in:
            with open(outpath_unzipped, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"  Decompressed: {outpath_unzipped}")
    except Exception as e:
        print(f"  ERROR: {e}")

print("\nDone!")
