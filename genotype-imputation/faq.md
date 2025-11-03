# Frequently Asked Questions

This page answers the most common questions about our genotype imputation service.

## General Questions

### What is genotype imputation?

Genotype imputation is a statistical technique that predicts missing genotypes based on observed genotypes and reference panels of known haplotypes. It allows researchers to analyze genetic variants that were not directly genotyped in their study.

### How accurate is imputation?

Imputation accuracy depends on several factors:
- **Reference panel size and diversity**
- **Population matching between study and reference**
- **Density of genotyped markers**
- **Minor allele frequency of imputed variants**

Generally, common variants (MAF > 5%) achieve R² > 0.8, while rare variants may have lower accuracy.

### Is the service really free?

Yes, our imputation service is completely free for both academic and commercial use. There are no hidden fees or usage restrictions.

### What is the maximum sample size?

Currently, we support up to **110,000 samples** per job. For larger datasets, contact us to discuss options.

## Data Preparation {#data-preparation}

### What file formats are supported?

We support only one format:
- **VCF files**: Uncompressed .vcf files only

If your data is in PLINK or other formats, you must convert to VCF before uploading.

### My VCF file won't upload. What's wrong?

Common issues include:
- **Compressed files**: We only accept uncompressed .vcf files
- **Invalid headers**: Check VCF format compliance
- **Large file size**: Files must be < 2GB each
- **Character encoding**: Use UTF-8 encoding

```bash
# Fix compressed VCF
gunzip yourfile.vcf.gz

# Validate VCF format
vcf-validator yourfile.vcf
```

### How do I know which genome build my data uses?

Check your data documentation or use our automatic detection. Common indicators:
- **Build 37/hg19**: Chromosome 1 ends around position 249,250,621
- **Build 38/hg38**: Chromosome 1 ends around position 248,956,422

### Should I filter my data before upload?

Basic quality filtering is recommended:
```bash
# Example PLINK filtering
plink --bfile input \
      --geno 0.05 \    # Remove variants with >5% missing
      --mind 0.05 \    # Remove samples with >5% missing  
      --maf 0.01 \     # Remove rare variants
      --hwe 1e-6 \     # Remove HWE violations
      --make-bed \
      --out filtered
```

## Reference Panels

### Which reference panel should I choose?

Choose based on your study population:

| Population | Recommended Panel |
|------------|-------------------|
| European | HRC, 1000G EUR |
| African | 1000G AFR, CAAPA |
| East Asian | 1000G EAS |
| South Asian | 1000G SAS |
| Mixed/Unknown | 1000G ALL |

### What's the difference between HRC and 1000G?

- **HRC (Haplotype Reference Consortium)**:
  - 32,470 samples
  - 39 million variants
  - European-focused
  - Higher accuracy for Europeans

- **1000 Genomes Project**:
  - 2,504 samples
  - 84 million variants
  - Global populations
  - Better for non-Europeans

### Can I use my own reference panel?

Currently, we only support our curated reference panels. Contact us if you have a specific need for custom panels.

## Job Submission and Processing

### How long does imputation take?

Processing time depends on dataset size:
- **1,000 samples**: 2-4 hours
- **10,000 samples**: 4-8 hours
- **50,000 samples**: 8-16 hours
- **100,000 samples**: 16-24 hours

### My job failed. What happened?

Check your email for the QC report, which details the failure reason:
- **Format errors**: File format issues
- **QC failures**: Too many variants/samples excluded
- **System errors**: Temporary technical issues

### Can I check job progress?

Yes, monitor progress through:
- Your user dashboard
- Email notifications
- Real-time status updates

### My job is stuck at "Running". Is this normal?

Jobs can appear stuck during intensive processing phases:
- **Phasing**: Can take several hours for large datasets
- **Imputation**: Varies by chromosome and reference panel size
- **Post-processing**: Final quality control and formatting

If stuck for >48 hours, contact support.

## Results and Downloads

### How do I interpret imputation quality scores?

**R² (R-squared)**:
- \> 0.8: High quality, recommended for analysis
- 0.5-0.8: Moderate quality, use with caution
- < 0.5: Low quality, consider excluding

**Info Score**:
- \> 0.7: Good imputation quality
- 0.4-0.7: Moderate quality
- < 0.4: Poor quality

### What's the difference between dosage and hard-call genotypes?

- **Dosage**: Continuous values (0-2) representing expected allele count
- **Hard-calls**: Discrete genotypes (0/0, 0/1, 1/1) based on highest probability

Dosage is generally preferred for association analyses.

### How long are results available?

Results are available for **7 days** after job completion. Download them promptly or request an extension.

### Can I get results in PLINK format?

Our service only outputs VCF format. You can convert the VCF results to PLINK format using:
```bash
plink --vcf imputed_chr1.vcf.gz \
      --make-bed \
      --out imputed_chr1
```

## Quality Control

### Many of my variants were excluded. Why?

Common reasons for variant exclusion:
- **Strand mismatches**: Alleles don't match reference panel
- **Low call rate**: Too much missing data
- **Allele frequency discrepancy**: Large difference from reference
- **Duplicate positions**: Multiple variants at same position

### What does "population outlier" mean?

Samples that don't cluster with the selected population based on principal component analysis. This may indicate:
- **Ancestry mismatch**: Wrong population specified
- **Sample mix-up**: Contamination or labeling errors
- **Admixed ancestry**: Mixed genetic background

### How do I fix strand issues?

1. **Download strand files** for your genotyping array
2. **Apply corrections** using PLINK or similar tools
3. **Re-upload** corrected data

```bash
# Example strand correction
plink --bfile original_data \
      --flip strand_corrections.txt \
      --make-bed \
      --out corrected_data
```

## Technical Issues

### The website is slow/unresponsive

Try these steps:
1. **Clear browser cache** and cookies
2. **Try different browser** (Chrome, Firefox, Safari)
3. **Check internet connection** stability
4. **Disable browser extensions** that might interfere

### Upload keeps failing

Common solutions:
- **Check file size**: Must be < 2GB per file
- **Verify file format**: Ensure proper format
- **Stable connection**: Use wired connection if possible
- **Try smaller batches**: Split large datasets

### I can't download my results

Possible issues:
- **Results expired**: Download within 7 days
- **Browser settings**: Check popup blockers
- **File size**: Large files may need alternative download methods

## Account and Authentication

### I forgot my password

Use the "Forgot Password" link on the login page. Check your email for reset instructions.

### Can I change my email address?

Currently, email addresses cannot be changed. Contact support if you need to update your account information.

### Do I need to create separate accounts for different projects?

No, one account can handle multiple projects. Jobs are organized in your dashboard for easy management.

## Error Messages

### "Invalid VCF format"

**Solutions**:
- Validate VCF with `vcf-validator`
- Check header format and required fields
- Ensure proper chromosome encoding

### "Reference panel mismatch"

**Solutions**:
- Verify genome build (hg19 vs hg38)
- Check population ancestry
- Review allele encoding (forward strand)

### "Insufficient variants for imputation"

**Solutions**:
- Reduce quality filters
- Check variant density
- Consider different reference panel

## Advanced Topics

### Can I submit pre-phased data?

Yes, specify "pre-phased" during job submission. Ensure your data is properly phased and formatted.

### How do I handle chromosome X?

Chromosome X requires special consideration:
- **Male samples**: Hemizygous genotypes (0 or 2, not 1)
- **Female samples**: Diploid genotypes (0, 1, or 2)
- **Sex information**: Must be correctly specified

### What about structural variants?

Currently, we only impute SNPs and small indels. Structural variants are not supported.

## Troubleshooting {#troubleshooting}

### Job fails immediately after submission

1. **Check file format** and size requirements
2. **Validate data integrity** with appropriate tools
3. **Review error messages** in email notifications
4. **Contact support** with job ID and error details

### Low imputation quality across all variants

1. **Verify population match** with reference panel
2. **Check data preprocessing** quality
3. **Review strand orientation** and allele coding
4. **Consider alternative reference panel**

### Missing expected variants in results

1. **Check QC report** for exclusion reasons
2. **Verify variant naming** (rs IDs vs positions)
3. **Review reference panel** coverage
4. **Examine allele frequency** differences

## Getting Additional Help

### When should I contact support?

Contact us when:
- Jobs fail repeatedly with unclear error messages
- You need help with data preparation
- You have questions about specific populations
- You need assistance with large-scale studies

### What information should I include?

When contacting support, provide:
- **Job ID** (if applicable)
- **Error messages** or screenshots
- **Data description** (platform, population, sample size)
- **Specific questions** or issues encountered

### Response times

- **General questions**: 1-2 business days
- **Technical issues**: Within 24 hours
- **Urgent problems**: Same day response when possible

For immediate assistance, check our [resources page](../resources.md) for additional tools and documentation. 