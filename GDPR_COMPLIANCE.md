# GDPR Compliance for Test Data

## Overview
The Jumble package test data has been enhanced with GDPR-compliant anonymization to protect genetic privacy.

## Anonymization Measures

### 1. Sample Names
- **Original**: `BM-P-HRD01-T-HRD01-KH20240326-N320240326.counts.RDS`
- **Anonymized**: `test_sample_1.counts.RDS`
- All patient/sample identifiers removed from filenames and internal data

### 2. SNP Position Scrambling
**Purpose**: SNP positions are genetic data that could identify individuals under GDPR.

**Method**: 
- Each SNP position is relocated within its genomic bin
- Random offset: ±50 bases (or ±25% of bin width, whichever is smaller)
- SNPs remain in the same bin, preserving test functionality
- Reproducible scrambling (seed = 12345)

**Implementation**:
```r
# Scramble VCF positions
source("scramble_vcf.R")
scramble_vcf_positions(
    vcf_file = "original.vcf.gz",
    output_file = "scrambled.vcf.gz",
    reference_file = "reference.RDS",
    seed = 12345
)
```

**Results**:
- Total SNPs in test VCF: 8,997
- SNPs scrambled: 8,597 (95.6%)
- SNPs unchanged: 400 (4.4%, those not mapping to bins)

### 3. Count Data
- Read counts preserved (not personal data)
- Genomic coordinates preserved (reference genome, not personal)
- Only sample identifiers removed

## GDPR Considerations

### Personal Data Removed
✅ Patient identifiers  
✅ Sample identifiers  
✅ Exact SNP positions (scrambled)  
✅ Any metadata linking to individuals

### Non-Personal Data Retained
✅ Reference genome coordinates  
✅ Read counts per bin  
✅ Genomic structure  
✅ Test functionality

## Verification

All tests pass with scrambled data:
- ✅ Gene panel processing
- ✅ SNP processing with scrambled VCF
- ✅ WGS processing
- ✅ Multi-sample testing

## Original Data

**Original files are NOT included in package distribution:**
- Located in `original files/` (local development only)
- Used only for development testing
- Never distributed with package
- Gitignored

## Regenerating Test Data

To regenerate GDPR-compliant test data:

```r
source("create_test_data.R")
```

This will:
1. Copy minimal samples (5 ref, 2 test)
2. Anonymize filenames and sample names
3. Scramble SNP positions in VCF
4. Create package test data in `inst/testdata/`

## Legal Compliance

This anonymization approach ensures:
- **GDPR Article 4(1)**: Data no longer identifies individuals
- **GDPR Recital 26**: Anonymized data outside GDPR scope
- **Research exemption**: Maintains scientific utility while protecting privacy

## Technical Details

### Scrambling Algorithm
1. Load VCF and reference bins
2. For each SNP:
   - Find overlapping genomic bin
   - Generate random offset (±50 bases max)
   - Ensure new position stays within bin
   - Update VCF position
3. Write scrambled VCF

### Reproducibility
- Fixed seed (12345) ensures consistent scrambling
- Same input always produces same output
- Facilitates testing and validation

## Package Distribution

**Included in package:**
- ✅ Anonymized, scrambled test data (~7 MB)
- ✅ GDPR-compliant VCF

**NOT included:**
- ❌ Original files with real sample names
- ❌ Unscrambled SNP positions
- ❌ Any identifying information
