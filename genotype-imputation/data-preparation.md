# Data Preparation

Proper data preparation is crucial for successful imputation. This guide covers all the requirements and best practices for preparing your genotype data.

::: warning VCF Format Only
Our service accepts **only VCF format** for input. If your data is in PLINK, 23andMe, or other formats, you must convert it to VCF before uploading. See the conversion instructions below.
:::

## Supported Input Format

### VCF Files (Only Supported Format)

VCF (Variant Call Format) is the only supported input format for our service.

**Requirements:**

- **bgzip-compressed VCF files** (`.vcf.gz`) - files must be compressed with bgzip (not regular gzip)
- Separate VCF file for each chromosome
- Must contain genotype information (GT field)
- Variants must be sorted by genomic position
- Chromosome encoding: 1-22, X (or chr1-chr22, chrX)
- No missing chromosome information

**Example VCF header:**
```
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1
1       752566  rs3094315       G       A       .       PASS    AF=0.5  GT      0|1
```

### Converting from Other Formats

If your data is in PLINK or other formats, you'll need to convert to VCF before uploading:

**From PLINK to VCF:**

```bash
# Convert PLINK binary format to bgzip-compressed VCF
plink --bfile input_data \
      --recode vcf bgz \
      --out converted_data

# This creates converted_data.vcf.gz compressed with bgzip
# Split by chromosome if needed:
for chr in {1..22} X; do
  plink --bfile input_data \
        --chr $chr \
        --recode vcf bgz \
        --out converted_data_chr${chr}
done
```

**From 23andMe to VCF:**
```bash
# First convert to PLINK format, then to VCF
# (See conversion scripts in our Resources section)
```

## Genome Build Requirements

### Build 37 (hg19/GRCh37)
- Most widely supported
- Recommended for most analyses
- Required for HRC and 1000G reference panels

### Build 38 (hg38/GRCh38)
- Supported for TOPMed and newer reference panels
- Automatic liftover available for some panels

::: tip Build Detection
Our system automatically detects your genome build, but you can specify it manually during job submission.
:::

## Quality Control Requirements

### Minimum Requirements

| Metric | Threshold | Description |
|--------|-----------|-------------|
| Sample call rate | ≥ 95% | Percentage of non-missing genotypes per sample |
| Variant call rate | ≥ 95% | Percentage of non-missing genotypes per variant |
| Hardy-Weinberg p-value | ≥ 1e-6 | Test for genotyping errors |
| Minor allele frequency | ≥ 1% | For quality control variants |

### Recommended Filtering

```bash
# Example PLINK commands for pre-processing
plink --bfile input_data \
      --geno 0.05 \        # Remove variants with >5% missing data
      --mind 0.05 \        # Remove samples with >5% missing data
      --maf 0.01 \         # Remove variants with MAF < 1%
      --hwe 1e-6 \         # Remove variants failing HWE test
      --make-bed \
      --out clean_data
```

## Chromosome X Considerations

Special handling is required for chromosome X:

- **Males**: Hemizygous genotypes (0 or 1, not 0/0 or 1/1)
- **Females**: Diploid genotypes (0/0, 0/1, 1/1)
- **Pseudoautosomal regions**: Treated as autosomal
- **Sample sex**: Must be correctly specified

## Population and Ancestry

### Recommended Approach

1. **Know your population**: Choose appropriate reference panel
2. **Mixed ancestry**: Use "Mixed" population or contact support
3. **Population outliers**: Will be flagged during QC

### Population Options

| Code | Description | Recommended Reference |
|------|-------------|----------------------|
| EUR | European | HRC, 1000G EUR |
| AFR | African | 1000G AFR, CAAPA |
| AMR | Admixed American | 1000G AMR |
| EAS | East Asian | 1000G EAS |
| SAS | South Asian | 1000G SAS |
| Mixed | Mixed/Unknown | 1000G ALL |

## Common Issues and Solutions

### File Format Issues

::: danger Invalid VCF Format
**Problem**: VCF file rejected during upload

**Solutions**:

- Ensure file is bgzip-compressed (`.vcf.gz`)
- Check that file is split by chromosome
- Verify variants are sorted by position
- Validate with `bcftools view -h yourfile.vcf.gz`
:::

::: danger Format Not Supported
**Problem**: Uploaded file is not in VCF format

**Solutions**:

- Convert PLINK files using: `plink --bfile data --recode vcf bgz --out data`
- Compress with bgzip: `bgzip yourfile.vcf`
- Ensure file extension is `.vcf.gz` (bgzip compressed)
:::

### Coordinate Issues

::: danger Build Mismatch
**Problem**: Coordinates don't match reference panel

**Solutions**:
- Verify your genome build
- Use our automatic liftover service
- Convert coordinates before upload
:::

### Quality Issues

::: danger Low Call Rate
**Problem**: Many samples/variants fail QC

**Solutions**:
- Apply stricter pre-filtering
- Check for batch effects
- Review genotyping quality
:::

## Pre-processing Checklist

Before uploading your data, verify:

- [ ] File format is supported and valid
- [ ] Genome build is known and consistent
- [ ] Sample and variant call rates meet requirements
- [ ] Chromosome encoding follows our standards
- [ ] Population ancestry is specified correctly
- [ ] File size is within limits (2GB per file)
- [ ] Sample size doesn't exceed 110,000

## Advanced Options

### Phasing

- **Pre-phased data**: Upload if already phased with high-quality method
- **Unphased data**: We'll phase using Eagle2 (recommended)
- **Mixed phasing**: Not supported - choose one approach

### Custom Reference Panels

Contact us if you need:
- Population-specific reference panels
- Custom variant sets
- Non-standard genome builds

## Getting Help

If you encounter data preparation issues:

1. Check our [FAQ](faq.md#data-preparation) for common solutions
2. Use our [validation tools](../resources.md#validation-tools)
3. [Contact support](../contact.md) with example data files

## Next Steps

After preparing your data:
- Review [imputation parameters](getting-started.md#imputation-parameters)
- Submit your job following our [getting started guide](getting-started.md)
- Monitor progress and download results 