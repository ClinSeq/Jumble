# Anonymized Test Data for Package Distribution

## Overview
All test data has been anonymized and copied to `inst/testdata/` for package distribution. Sample names have been removed from filenames and internal data structures.

## Test Data Structure

```
inst/testdata/
├── gene_panel/
│   ├── reference/
│   │   ├── ref_sample_1.counts.RDS
│   │   ├── ref_sample_2.counts.RDS
│   │   ├── ref_sample_3.counts.RDS
│   │   ├── ref_sample_4.counts.RDS
│   │   └── ref_sample_5.counts.RDS
│   ├── samples/
│   │   ├── test_sample_1.counts.RDS
│   │   ├── test_sample_1.vcf.gz (for SNP testing)
│   │   └── test_sample_2.counts.RDS
│   └── reference.RDS (pre-built reference)
└── wgs/
    ├── reference/
    │   ├── wgs_ref_1.counts.RDS
    │   ├── wgs_ref_2.counts.RDS
    │   └── wgs_ref_3.counts.RDS
    ├── samples/
    │   └── wgs_test_1.counts.RDS
    └── reference.RDS (pre-built reference)
```

## Anonymization Details

### Filenames
- Original: `[SAMPLE_NAME].counts.RDS`
- Anonymized: `test_sample_1.counts.RDS`

### Internal Data
- Sample names in count objects changed to generic names (e.g., `ref_sample_1`, `test_sample_1`)
- Reference `target_template$sample` column set to empty string
- Reference `samples` vector contains anonymized names
- All other data (genomic coordinates, counts, etc.) preserved

## Test Independence

All tests now use package test data exclusively:

```r
# Tests automatically find package test data
testdata_dir <- system.file("testdata", package = "Jumble")
if (testdata_dir == "") testdata_dir <- "inst/testdata" # For development
```

**No dependency on original files** - Tests work with installed package or during development.

## Test Results

All tests pass using anonymized data:

| Test | Status | Details |
|------|--------|---------|
| Gene panel | ✅ SUCCESS | 327 segments (22 X), SNP file created |
| WGS | ✅ SUCCESS | 383 segments (11 X) |
| Multi-sample | ✅ SUCCESS | 2/2 samples processed |

## Creating Test Data

To regenerate anonymized test data (requires original files):

```r
source("create_test_data.R")
```

This script:
1. Copies 5 gene panel reference samples
2. Copies 2 gene panel test samples + 1 VCF
3. Copies 3 WGS reference samples
4. Copies 1 WGS test sample
5. Anonymizes all filenames and internal sample names
6. Builds pre-computed reference files

## Package Size

Test data adds approximately:
- Gene panel: ~2 MB (5 ref + 2 test samples + reference)
- WGS: ~5 MB (3 ref + 1 test sample + reference)
- **Total: ~7 MB**

## Original Files

Original files in `original files/` directory are:
- ✅ **Kept for development reference** (as requested)
- ❌ **Not used by tests**
- ❌ **Not included in package distribution**

Tests are completely independent and use only `inst/testdata/`.
