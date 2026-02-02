"""
Download and parse GSE157153 dataset for UPR analysis.
Dataset: Gene expression profiling in endometriosis and ovarian cancer
Platform: GPL17303 (Ion Torrent Proton RNA-seq)
Samples: 66
"""

import GEOparse
import pandas as pd
import os
from pathlib import Path

# Set output directory
OUT_DIR = Path(__file__).parent
OUT_DIR.mkdir(exist_ok=True)

print("=" * 60)
print("Downloading GSE157153...")
print("=" * 60)

# Download GEO dataset (this caches locally)
gse = GEOparse.get_GEO(geo="GSE157153", destdir=str(OUT_DIR))

print(f"\nDataset title: {gse.metadata['title'][0]}")
print(f"Platform: {gse.metadata['platform_id'][0]}")
print(f"Number of samples: {len(gse.gsms)}")

# Extract sample metadata
print("\n" + "=" * 60)
print("Extracting sample metadata...")
print("=" * 60)

sample_info = []
for gsm_name, gsm in gse.gsms.items():
    info = {
        'sample_id': gsm_name,
        'title': gsm.metadata.get('title', [''])[0],
        'source': gsm.metadata.get('source_name_ch1', [''])[0],
    }
    # Extract characteristics
    chars = gsm.metadata.get('characteristics_ch1', [])
    for char in chars:
        if ':' in char:
            key, val = char.split(':', 1)
            info[key.strip()] = val.strip()
    sample_info.append(info)

samples_df = pd.DataFrame(sample_info)
samples_df.to_csv(OUT_DIR / "sample_metadata.csv", index=False)
print(f"Sample metadata saved: {len(samples_df)} samples")
print("\nSample groups:")
if 'tissue' in samples_df.columns:
    print(samples_df['tissue'].value_counts())
elif 'source' in samples_df.columns:
    print(samples_df['source'].value_counts())

# Check for expression data in soft file
print("\n" + "=" * 60)
print("Checking expression data...")
print("=" * 60)

# Try to get expression table from GSE
if hasattr(gse, 'pivot_samples'):
    try:
        expr_df = gse.pivot_samples('VALUE')
        if expr_df is not None and len(expr_df) > 0:
            expr_df.to_csv(OUT_DIR / "expression_matrix.csv")
            print(f"Expression matrix: {expr_df.shape[0]} genes x {expr_df.shape[1]} samples")
        else:
            print("No expression data in SOFT file - need to check supplementary files")
    except Exception as e:
        print(f"Could not extract expression matrix: {e}")
else:
    print("No pivot_samples method - checking supplementary files")

# List supplementary files
print("\n" + "=" * 60)
print("Supplementary files:")
print("=" * 60)
suppl = gse.metadata.get('supplementary_file', [])
for f in suppl:
    print(f"  - {f}")

print("\n" + "=" * 60)
print("Download complete!")
print("=" * 60)
print(f"\nFiles saved to: {OUT_DIR}")
