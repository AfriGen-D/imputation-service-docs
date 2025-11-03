# Changes in Allele Swap Handling

::: warning Important Changes
With the release of Imputation Service 2, we have made significant changes to how allele swaps are handled. Please read this page carefully if you are upgrading from the previous version.
:::

## Overview

In the previous version of our imputation service, we automatically corrected strand flips and allele swaps to maximize the number of variants that could be imputed. While this approach was convenient, it could sometimes mask underlying data quality issues and lead to unexpected results.

## What Changed

### Previous Behavior (v1.x)
- Automatic strand flip correction
- Automatic allele swap correction
- Silent corrections without detailed reporting
- Priority on maximizing variant inclusion

### New Behavior (v2.0+)
- **Strict allele matching** with reference panels
- **No automatic allele swaps** during imputation
- **Detailed QC reporting** of all checks performed
- **Clear documentation** of any corrections made

## Impact on Your Data

### Variants That Will Be Excluded

The following types of variants will now be excluded from imputation:

1. **Strand mismatches** that cannot be resolved unambiguously
2. **Allele swaps** where REF/ALT don't match the reference panel
3. **Ambiguous variants** (A/T and G/C SNPs) with frequency discrepancies
4. **Triallelic variants** not present in the reference panel

### Quality Control Reports

You will now receive detailed reports showing:

- Number of variants excluded for each reason
- Chromosome-specific statistics
- Allele frequency comparisons
- Strand flip analysis

## Migration Guide

### For Existing Users

If you have been using our service previously:

1. **Review your data preparation**: Ensure proper strand orientation
2. **Check reference panel compatibility**: Verify your data matches the chosen panel
3. **Update preprocessing scripts**: Include strand flip checks
4. **Expect different results**: Variant counts may be lower but more reliable

### Recommended Preprocessing

```bash
# 1. Check strand alignment with reference panel
plink --bfile your_data \
      --reference-panel reference.bim \
      --strand-check \
      --out strand_check

# 2. Apply strand flips if needed
plink --bfile your_data \
      --flip strand_check.flip \
      --make-bed \
      --out corrected_data

# 3. Remove problematic variants
plink --bfile corrected_data \
      --exclude problematic_variants.txt \
      --make-bed \
      --out final_data
```

## Best Practices

### Before Upload

1. **Use latest reference sequences**: Ensure your data uses the same reference genome
2. **Validate strand orientation**: Check against reference panel
3. **Review allele frequencies**: Compare with population-appropriate references
4. **Document data source**: Know the platform and preprocessing applied

### Data Quality Checks

```bash
# Check allele frequencies against 1000G
plink --bfile your_data \
      --freq \
      --out your_data_freq

# Compare with reference frequencies
python compare_frequencies.py your_data_freq.frq reference_freq.frq
```

### Handling Excluded Variants

If many variants are excluded:

1. **Review QC report**: Understand reasons for exclusion
2. **Check data preparation**: Ensure proper preprocessing
3. **Consider different reference panel**: May have better coverage
4. **Contact support**: For guidance on specific issues

## Tools and Resources

### Strand Check Tools

- **GTOOL**: For Oxford format data
- **PLINK**: Built-in strand checking
- **Will Rayner's Tools**: Comprehensive strand checking scripts

### Reference Files

Download strand files for common arrays:
- Illumina arrays: GSA, MEGA, Omni series
- Affymetrix arrays: Axiom, SNP6
- Custom arrays: Contact us for assistance

## FAQ

### Why did you make this change?

**Answer**: To improve result reliability and transparency. The previous automatic corrections could mask data quality issues and lead to inconsistent results across different datasets.

### Will my imputation results be different?

**Answer**: Yes, you may see:
- Fewer variants passing QC
- Different imputation accuracy
- More consistent results across batches
- Better documentation of QC decisions

### How do I fix strand issues?

**Answer**: 
1. Download strand files for your array
2. Use PLINK or similar tools to flip strands
3. Re-upload corrected data
4. Review QC reports carefully

### Can I still use my old data?

**Answer**: Yes, but you may need to reprocess it according to our new guidelines. We recommend reviewing the [data preparation guide](data-preparation.md).

### What if I have custom array data?

**Answer**: Contact our support team. We can help you prepare custom strand files and provide guidance on data preparation.

## Support

If you need help adapting to these changes:

1. Review our updated [data preparation guide](data-preparation.md)
2. Check the [FAQ](faq.md) for common issues
3. Use our [validation tools](../resources.md#validation-tools)
4. [Contact support](../contact.md) with specific questions

## Timeline

- **September 2024**: New allele swap handling implemented
- **October 2024**: Updated documentation and tools
- **November 2024**: Migration guide and support resources
- **December 2024**: Legacy support ends

::: tip Getting Help
Our team is here to help you through this transition. Don't hesitate to reach out if you need assistance preparing your data for the new system.
::: 