# Jumble Testing Strategy

## Two-Tier Testing Approach

### 1. Package Tests (Included in Distribution)
**Location**: `tests/testthat/`  
**Data**: Anonymized data in `inst/testdata/`  
**Purpose**: Verify package functionality for end users

**Run with**:
```r
devtools::test()
# or
source("run_all_tests.R")
```

**Coverage**:
- ✅ Basic pipeline functionality
- ✅ Gene panel processing (2 samples)
- ✅ WGS processing (1 sample)
- ✅ SNP processing with VCF
- ✅ Reference building

**Data Size**: ~7 MB (anonymized)

---

### 2. Development Tests (Local Only)
**Location**: `dev_tests/` (gitignored)  
**Data**: Original files with real sample names  
**Purpose**: Validate against original implementation, catch regressions

**Run with**:
```r
source("dev_tests/run_dev_tests.R")
```

**Coverage**:
- ✅ All 8 gene panel test samples with VCFs
- ✅ WGS sample processing
- ✅ Reference building from original samples
- ✅ Exact compatibility with original reference
- ✅ Detailed per-sample validation

**Data Size**: Full original dataset

---

## Comparison

| Aspect | Package Tests | Development Tests |
|--------|--------------|-------------------|
| **Distribution** | ✅ Included | ❌ Local only |
| **Data** | Anonymized | Original |
| **Sample Names** | Generic (test_sample_1) | Real (BM-P-HRD01-...) |
| **Sample Count** | Minimal (2-3) | Complete (8+) |
| **Purpose** | End-user validation | Development validation |
| **Run Frequency** | Every build | Before releases |
| **CI/CD** | Yes | No |

---

## When to Use Each

### Package Tests
- ✅ During development (quick validation)
- ✅ Before commits
- ✅ In CI/CD pipelines
- ✅ For package users

### Development Tests
- ✅ Before major releases
- ✅ After significant refactoring
- ✅ When validating against original implementation
- ✅ For comprehensive regression testing

---

## Test Data Locations

```
Jumble/
├── inst/testdata/              # Package tests (anonymized)
│   ├── gene_panel/
│   │   ├── reference/          # 5 samples
│   │   ├── samples/            # 2 samples + VCF
│   │   └── reference.RDS
│   └── wgs/
│       ├── reference/          # 3 samples
│       ├── samples/            # 1 sample
│       └── reference.RDS
│
├── original files/             # Development tests (original)
│   ├── reference_file/
│   ├── reference_samples/      # 36 samples
│   └── test_samples/           # 8 samples + VCFs
│
└── verification_test/
    └── wgs_test/               # WGS original data
        ├── reference/          # 3 samples
        └── test/               # 1 sample
```

---

## Running All Tests

### Quick (Package Only)
```r
devtools::test()
```

### Comprehensive (Package + Development)
```bash
Rscript run_all_tests.R              # Package tests
Rscript dev_tests/run_dev_tests.R    # Development tests
```

---

## Maintenance

### Adding Package Tests
1. Add test to `tests/testthat/test-*.R`
2. Use `system.file("testdata", package = "Jumble")`
3. Keep minimal (2-3 samples max)

### Adding Development Tests
1. Add test to `dev_tests/run_dev_tests.R`
2. Use original files directly
3. Test all available samples

### Updating Test Data
```r
# Regenerate anonymized package data
source("create_test_data.R")

# Original files: keep as-is for development reference
```
