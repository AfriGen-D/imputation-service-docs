# Getting Started

This guide will walk you through the process of submitting your first genotype imputation job.

## Overview

Genotype imputation is a statistical technique used to infer genotypes that are not directly assayed in a sample. Our service uses the latest Minimac4 algorithm to provide highly accurate imputation results.

## Step 1: Register and Login

1. Go to the imputation service website
2. Click "Register" to create a new account
3. Verify your email address
4. Login with your credentials

## Step 2: Prepare Your Data

Before uploading, make sure your data meets our requirements:

- **Supported format**: VCF only (uncompressed .vcf files)
- **Maximum sample size**: 110,000 samples
- **Coordinates**: Build 37 (hg19) or Build 38 (hg38)
- **Quality requirements**: See [Data Preparation](data-preparation.md) for details

::: warning Important
Please review our [data preparation guidelines](data-preparation.md) carefully to avoid common issues that can cause job failures.
:::

## Step 3: Create a New Job

1. Click "Run" in the main navigation
2. Select "Genotype Imputation"
3. Choose your reference panel (see [Reference Panels](../reference-panels.md))
4. Upload your genotype files
5. Select imputation parameters

### Imputation Parameters

| Parameter | Description | Options |
|-----------|-------------|---------|
| Reference Panel | Population-specific reference for imputation | HRC, 1000G, TOPMed, etc. |
| Phasing | Whether input data is pre-phased | Eagle2, No phasing |
| Population | Population for quality control | EUR, AFR, AMR, EAS, SAS, Mixed |
| Build | Genome build of input data | hg19/GRCh37, hg38/GRCh38 |

## Step 4: Quality Control

After upload, our system automatically performs quality control:

- **Chromosome check**: Validates chromosome encoding
- **Strand check**: Identifies and corrects strand flips
- **Allele frequency check**: Compares with reference panel
- **Sample/variant filtering**: Removes low-quality data

You will receive a QC report via email with detailed statistics.

## Step 5: Monitor Progress

- Check job status on your dashboard
- Receive email notifications at key stages
- View real-time progress updates

## Step 6: Download Results {#downloading-results}

Once imputation is complete (typically 1-24 hours):

1. Go to your job dashboard
2. Click "Download" next to your completed job
3. Download individual chromosome files or the complete dataset
4. Results are available for 7 days

### Output Files

- **Imputed genotypes**: VCF files with dosage information
- **Quality metrics**: Info scores, R² values
- **Statistics**: Per-variant and per-sample statistics
- **QC report**: Detailed quality control report

## Example Workflow

```bash
# 1. Prepare your VCF file
bcftools view input.vcf.gz > prepared_input.vcf

# 2. Upload via web interface
# (Use the web interface for file upload)

# 3. Download results
wget https://imputationserver.com/results/job_123/chr1.dose.vcf.gz
```

## Next Steps

- Learn about [data preparation best practices](data-preparation.md)
- Understand our [pipeline overview](pipeline-overview.md)
- Check the [FAQ](faq.md) for common questions

## Getting Help

If you encounter issues:

1. Check our [FAQ](faq.md) for common solutions
2. Review the [troubleshooting guide](faq.md#troubleshooting)
3. [Contact support](../contact.md) with your job ID 