# Jumble Test Suite

This directory contains comprehensive tests for the Jumble package.

## Test Structure

### 1. Unit Tests (`tests/testthat/`)
Formal R package tests using testthat framework:
- `test-pipeline.R`: Full pipeline integration test
- `test-regression.R`: Regression tests against known good outputs
- `test-reference-building.R`: Reference building for gene panel and WGS
- `test-snp-processing.R`: SNP VCF processing and plotting

### 2. Test Data (`test_data/`)
Organized test data for different scenarios:

#### Gene Panel Data
- **Reference samples**: `/Users/markus/.gemini/antigravity/scratch/Jumble/original files/reference_samples/`
  - 36 samples for building gene panel reference
- **Test samples**: `/Users/markus/.gemini/antigravity/scratch/Jumble/original files/test_samples/`
  - 8 samples with VCF files for SNP testing
- **Reference file**: `/Users/markus/.gemini/antigravity/scratch/Jumble/original files/reference_file/pancancer3_customized_targets_twist.bed.reference.RDS`

#### WGS Data
- **Reference samples**: `verification_test/wgs_test/reference/`
  - 3 WGS samples for building WGS reference
- **Test samples**: `verification_test/wgs_test/test/`
  - 1 WGS sample for testing

### 3. Master Test Script (`run_all_tests.R`)
Comprehensive test script that runs all tests:
- Gene panel reference building
- WGS reference building
- Gene panel sample processing (with and without SNPs)
- WGS sample processing
- Compatibility with original reference

## Running Tests

### Quick Test (testthat)
```r
devtools::test()
```

### Comprehensive Test
```r
source("run_all_tests.R")
```

### Individual Components
```r
# Gene panel only
source("tests/testthat/test-regression.R")

# WGS only
Rscript -e "library(Jumble); run_jumble(
  bam_file = 'verification_test/wgs_test/test/SARC-P-MOL465424-T-MOL465424-KH20250129-WG20250129_markdups.bam.counts.RDS',
  reference_file = 'verification_test/comprehensive/wgs_reference.RDS',
  output_dir = tempdir()
)"
```

## Test Coverage

- ✅ Gene panel reference building
- ✅ WGS reference building
- ✅ Gene panel sample processing
- ✅ WGS sample processing
- ✅ SNP VCF processing and plotting
- ✅ X/Y chromosome handling
- ✅ Segmentation
- ✅ Normalization
- ✅ Plotting (with and without SNPs)
- ✅ Output file generation (.cnr, .cns, .seg, .RDS, .png)

## Cleanup

Obsolete directories have been removed. Current structure maintains only:
- `verification_test/wgs_test/` - WGS test data
- `verification_test/comprehensive/` - Built references for testing
