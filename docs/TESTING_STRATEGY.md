# Jumble Testing Strategy

## Package Tests 

**Location**: `tests/testthat/`  
**Data**: Anonymized data in `inst/testdata/`  
**Purpose**: Verify package functionality for end users and developers.

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

## Test Data Locations

The package includes minimal anonymized testing data to validate pipeline execution.

```
Jumble/
├── inst/testdata/              # Package tests (anonymized)
│   ├── gene_panel/
│   │   ├── reference/          # Reference build counts
│   │   ├── samples/            # Test samples + VCF
│   │   └── reference.RDS
│   └── wgs/
│       ├── reference/          # Reference build counts
│       ├── samples/            # Test sample
│       └── reference.RDS
```

## Maintenance

### Adding Package Tests
1. Add new test logic to `tests/testthat/test-*.R`
2. Programmatically load data via `system.file("testdata", package = "Jumble")`
3. Keep the file sizes of any new test artifacts minimal to avoid bloating the R package distribution.
