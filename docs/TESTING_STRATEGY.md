# Jumble Testing Strategy

Jumble uses a three-tier testing architecture. For the full description see `.roo/TESTING.md`.

---

## Tier 1 — Formal Package Tests

**Location**: `tests/testthat/`
**Data**: Anonymized data in `inst/testdata/`
**Purpose**: Verify package functionality for end users and developers. Included in package distribution and run by CRAN/Bioconductor checks.

**Run with**:
```r
devtools::test()
```

**Coverage**:
- ✅ Basic pipeline functionality
- ✅ Gene panel processing
- ✅ WGS processing
- ✅ SNP processing with VCF
- ✅ Reference building
- ✅ QC metrics (21 columns)
- ✅ MSI classification
- ✅ GIS/HRD scoring
- ✅ Contamination estimation (RF model loading, feature extraction, edge cases)

---

## Tier 2 — Dev Tests (Real Data)

**Location**: `local_dev/dev_tests/`
**Data**: Real (non-anonymized) data in `local_dev/dev_test_input/` — never committed
**Purpose**: Validate pipeline correctness on real clinical data after code changes.

**Run with**:
```r
source("local_dev/dev_tests/run_dev_tests.R")   # Full suite (hg19 + hg38 + somatic)
source("local_dev/dev_tests/run_hg19_tests.R")  # hg19 only
source("local_dev/dev_tests/run_hg38_test.R")   # hg38 only
```

**Coverage**: hg19/hg38 gene panel, hg19/hg38 WGS, direct BAM input, somatic VCF integration

---

## Tier 3 — Additional / GMS-Solid Tests

**Location**: `local_dev/` root
**Data**: Real data in `local_dev/GMS-Solid/` and `local_dev/additional_tests_input/` — never committed
**Purpose**: Validate specific assay types and edge cases.

**Run with**:
```r
source("local_dev/run_gms_solid_with_snps.R")   # GMS-Solid panel (hg38, 4 samples)
source("local_dev/run_additional_tests.R")       # Additional samples
```

---

## Test Data Locations

```
Jumble/
├── inst/testdata/              # Tier 1: Package tests (anonymized, committed)
│   ├── gene_panel/
│   │   ├── reference/
│   │   ├── samples/
│   │   └── reference.RDS
│   └── wgs/
│       ├── reference/
│       ├── samples/
│       └── reference.RDS
└── local_dev/                  # Tiers 2 & 3: Real data (never committed)
    ├── dev_test_input/
    │   ├── hg19/
    │   └── hg38/
    ├── GMS-Solid/
    └── additional_tests_input/
```

---

## Maintenance

### Adding Package Tests
1. Add new test logic to `tests/testthat/test-*.R`
2. Programmatically load data via `system.file("testdata", package = "Jumble")`
3. Keep the file sizes of any new test artifacts minimal to avoid bloating the R package distribution.

### After Code Changes
Run all three tiers in order:
1. `devtools::test()` — must pass (or failures must be understood/accepted)
2. `source("local_dev/dev_tests/run_dev_tests.R")` — must pass
3. `source("local_dev/run_gms_solid_with_snps.R")` and `source("local_dev/run_additional_tests.R")`
