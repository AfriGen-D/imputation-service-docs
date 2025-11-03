# Common Errors and Fixes

This page provides solutions to frequently encountered errors when using the AfriGen-D Imputation Service. Use the search function (Ctrl+F / Cmd+F) to quickly find your specific error message.

## File Format Errors

### Error: "Invalid VCF format"

**Cause**: VCF file doesn't meet format specifications

**Solutions**:

1. **Validate your VCF file**:
   ```bash
   # Using vcf-validator
   vcf-validator your_file.vcf.gz

   # Using bcftools
   bcftools view -h your_file.vcf.gz | head -20
   ```

2. **Check VCF header**: Ensure header includes required fields:
   - `##fileformat=VCFv4.2` or `VCFv4.3`
   - `##contig` lines for chromosomes
   - `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT` column headers

3. **Fix common issues**:
   ```bash
   # Remove spaces in sample names
   bcftools reheader --samples samples.txt input.vcf.gz > fixed.vcf.gz

   # Sort variants by position
   bcftools sort input.vcf.gz -Oz -o sorted.vcf.gz
   ```

---

### Error: "File must be bgzip-compressed"

**Cause**: File is uncompressed or compressed with regular gzip instead of bgzip

**Solution**:

```bash
# If file is uncompressed (.vcf)
bgzip input.vcf
# Creates input.vcf.gz

# If file is gzip-compressed (check with: file myfile.vcf.gz)
gunzip input.vcf.gz
bgzip input.vcf

# Verify bgzip compression
file input.vcf.gz
# Should show: "gzip compressed data" (bgzip creates valid gzip format)

# Verify VCF can be indexed (bgzip-specific feature)
tabix -p vcf input.vcf.gz
```

**Installation** (if bgzip not available):
```bash
# Ubuntu/Debian
sudo apt-get install tabix

# macOS
brew install htslib

# Conda
conda install -c bioconda htslib
```

---

### Error: "Missing GT (genotype) field"

**Cause**: VCF file lacks genotype information in FORMAT column

**Solution**:

Check if GT field exists:
```bash
bcftools query -f '%CHROM\t%POS\t%FORMAT\n' input.vcf.gz | head -5
```

If GT is missing, you cannot use this file for imputation (genotype data is required). Regenerate VCF from source data with genotype information.

**Converting from PLINK with GT**:
```bash
plink --bfile input_data \
      --recode vcf bgz \
      --out output_with_gt
```

---

### Error: "One file per chromosome required"

**Cause**: VCF contains multiple chromosomes in a single file

**Solution**:

Split VCF by chromosome:

```bash
# Method 1: Using bcftools
for chr in {1..22} X; do
  bcftools view -r ${chr} input.vcf.gz -Oz -o chr${chr}.vcf.gz
  tabix -p vcf chr${chr}.vcf.gz
done

# Method 2: Using PLINK
for chr in {1..22} X; do
  plink --bfile input_data \
        --chr ${chr} \
        --recode vcf bgz \
        --out chr${chr}
done
```

---

## Quality Control Errors

### Error: "Insufficient variants passing QC filters"

**Cause**: Too few variants meet the service's quality control thresholds

**QC Thresholds**:
- Variant call rate ≥ 95%
- Sample call rate ≥ 95%
- Hardy-Weinberg equilibrium p-value ≥ 1e-6
- Minor allele frequency ≥ 1% (for QC, not imputation)

**Solutions**:

1. **Check your data quality**:
   ```bash
   # Calculate missingness
   plink --bfile input_data --missing --out qc_report

   # Calculate HWE
   plink --bfile input_data --hardy --out qc_report

   # Calculate MAF
   plink --bfile input_data --freq --out qc_report
   ```

2. **Pre-filter your data** (if appropriate):
   ```bash
   # Apply standard QC filters
   plink --bfile input_data \
         --geno 0.05 \
         --mind 0.05 \
         --hwe 1e-6 \
         --maf 0.01 \
         --make-bed \
         --out qc_filtered

   # Convert to VCF
   plink --bfile qc_filtered \
         --recode vcf bgz \
         --out qc_filtered
   ```

3. **Consider your study design**: Very rare variants or highly stratified populations may fail standard QC

---

### Error: "Sample call rate too low"

**Cause**: One or more samples have excessive missing genotypes

**Solution**:

1. **Identify problematic samples**:
   ```bash
   plink --bfile input_data --missing --out sample_qc
   # Check sample_qc.imiss file, F_MISS column
   ```

2. **Remove samples with >5% missingness**:
   ```bash
   plink --bfile input_data \
         --mind 0.05 \
         --make-bed \
         --out filtered_samples
   ```

3. **If missing data is biologically relevant**: Contact support to discuss options

---

### Error: "Variant positions not in reference genome"

**Cause**: Chromosome positions don't match the selected genome build

**Solutions**:

1. **Verify your genome build**:
   ```bash
   # Check maximum position on chromosome 1
   bcftools view -H input.vcf.gz | grep "^1\s" | awk '{print $2}' | sort -n | tail -1

   # Build 37: ~249,250,621
   # Build 38: ~248,956,422
   ```

2. **Lift over to correct build** (if needed):
   ```bash
   # Using Picard LiftoverVcf
   java -jar picard.jar LiftoverVcf \
        I=input_hg19.vcf.gz \
        O=output_hg38.vcf.gz \
        CHAIN=hg19ToHg38.over.chain \
        REJECT=rejected_variants.vcf.gz \
        R=hg38.fasta
   ```

3. **Match reference panel build**:
   - H3Africa v6: GRCh38/hg38 only
   - 1000 Genomes: GRCh37/hg19 or GRCh38/hg38

---

## Chromosome-Specific Errors

### Error: "Invalid chromosome name"

**Cause**: Chromosome naming doesn't match expected format

**Expected formats**:
- `1`, `2`, ..., `22`, `X` (preferred)
- `chr1`, `chr2`, ..., `chr22`, `chrX` (also accepted)

**Solution**:

```bash
# Add "chr" prefix
bcftools annotate --rename-chrs chr_names.txt input.vcf.gz -Oz -o output.vcf.gz

# chr_names.txt content:
# 1 chr1
# 2 chr2
# ...

# Remove "chr" prefix
bcftools annotate --rename-chrs chr_names.txt input.vcf.gz -Oz -o output.vcf.gz

# chr_names.txt content:
# chr1 1
# chr2 2
# ...
```

---

### Error: "Chromosome X genotype coding error"

**Cause**: Male samples have diploid genotypes on non-pseudoautosomal regions of chromosome X

**Expected coding**:
- **Males (non-PAR regions)**: Haploid genotypes (0 or 1, not 0/0 or 0/1)
- **Females**: Diploid genotypes (0/0, 0/1, 1/1)
- **Pseudoautosomal regions (PAR)**: Diploid for both sexes

**Solution**:

```bash
# Using PLINK with sex information
plink --bfile input_data \
      --split-x hg19 \  # or hg38
      --make-bed \
      --out corrected_chrX

# Ensure sex information is in .fam file:
# Column 5: 1=male, 2=female
```

---

## Reference Panel Errors

### Error: "Selected reference panel not available for genome build"

**Cause**: Mismatch between data genome build and reference panel

**Reference Panel Builds**:
- **H3Africa v6**: GRCh38/hg38 only
- **1000 Genomes Phase 3**: GRCh37/hg19 or GRCh38/hg38

**Solution**:

Either:
1. **Lift over your data** to match reference panel build, OR
2. **Select compatible reference panel** for your genome build

---

### Error: "Reference panel selection failed"

**Cause**: Issue selecting or accessing the chosen reference panel

**Solution**:

1. **Verify panel selection**: Ensure you've selected a valid reference panel (H3Africa v6 or 1000 Genomes Phase 3)
2. **Check genome build compatibility**:
   - H3Africa v6: GRCh38/hg38 only
   - 1000 Genomes: GRCh37/hg19 or GRCh38/hg38
3. **Try alternative panel**: If issues persist, try the other available reference panel
4. **Contact support**: If problems continue, contact support with your job ID

---

## Upload and Processing Errors

### Error: "File size exceeds limit (2GB)"

**Cause**: Individual VCF file is larger than 2GB

**Solutions**:

1. **Verify file is split by chromosome**: Each chromosome should be a separate file

2. **Check file is bgzip-compressed** (not uncompressed):
   ```bash
   ls -lh *.vcf.gz
   # bgzipped files are typically 10-20% of uncompressed size
   ```

3. **For very large sample sizes**, consider:
   - Removing extremely rare variants (MAF < 0.001)
   - Splitting into smaller batches
   - Contact support for large-scale study assistance

---

### Error: "Job timeout / Processing time exceeded"

**Cause**: Job took longer than maximum allowed processing time

**Expected processing times**:
- 1,000 samples: 2-4 hours
- 10,000 samples: 6-12 hours
- 100,000 samples: 16-24 hours

**Solutions**:

1. **Reduce complexity**:
   - Filter ultra-rare variants
   - Ensure proper QC was performed
   - Split large batches

2. **Check for data issues**:
   - Extremely high missingness
   - Unusual variant density
   - Corrupted files

3. **Contact support** if issues persist with job ID

---

## Strand and Allele Errors

### Error: "Strand alignment issues detected"

**Cause**: Variant alleles on wrong strand compared to reference

**Solution**:

```bash
# Use strand files for your array
# Available at: https://www.well.ox.ac.uk/~wrayner/strand/

# Example for Illumina arrays
plink --bfile input_data \
      --reference-allele strand_file.txt \
      --make-bed \
      --out strand_aligned
```

---

### Error: "Allele frequency discrepancy"

**Cause**: Large difference between your data and reference panel allele frequencies

**Possible causes**:
1. **Population mismatch**: Wrong reference panel for your population
2. **Strand issues**: Alleles on opposite strand
3. **Build mismatch**: Wrong genome build
4. **Genotyping errors**: Poor quality genotyping data

**Solutions**:

1. **Verify population match**: Use appropriate reference panel
2. **Check strand alignment**: See "Strand alignment issues" above
3. **Confirm genome build**: See "Variant positions not in reference genome" above
4. **Review QC metrics**: Check genotyping quality reports

---

## Results and Download Errors

### Error: "Results expired (7-day limit)"

**Cause**: Results are automatically deleted after 7 days

**Solution**:

- **Prevention**: Download results promptly after job completion
- **Recovery**: Re-run imputation job (if you saved input files)
- **Future**: Set up email notifications for job completion

---

### Error: "Download incomplete or corrupted"

**Cause**: Network interruption during download

**Solutions**:

```bash
# Use wget with resume capability
wget -c https://impute.afrigen-d.org/results/job_123456/results.zip

# Or use curl with resume
curl -C - -O https://impute.afrigen-d.org/results/job_123456/results.zip

# Verify download integrity (if MD5 provided)
md5sum results.zip
```

---

## System and Account Errors

### Error: "Too many concurrent jobs"

**Cause**: Job submission limit reached

**Current limits**: Check service documentation for current limits

**Solution**:
- Wait for running jobs to complete
- Prioritize most important jobs
- Contact support for large-scale study accommodations

---

### Error: "Authentication failed"

**Cause**: Login credentials issue

**Solutions**:
1. Reset password through service interface
2. Clear browser cookies and cache
3. Try different browser
4. Contact support if issues persist

---

## Getting Additional Help

If your error isn't listed here or the solutions don't resolve your issue:

1. **Check the [FAQ](faq.md)** for additional answers
2. **Search the documentation** using the search bar
3. **Contact support** at [support@imputationserver.com](mailto:support@imputationserver.com)

When contacting support, include:
- Complete error message
- Job ID (if applicable)
- Steps already tried
- Data description (sample size, population, array type)

---

::: tip Pro Tip
Save your command history and log files when preparing data. This makes troubleshooting much faster if errors occur!
:::

**Last updated**: January 2025
