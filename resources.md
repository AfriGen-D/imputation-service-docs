# Resources

This page provides additional tools, documentation, and external resources to help you with genotype imputation.

## Validation Tools {#validation-tools}

### Data Validation

**VCF Validator**
```bash
# Install vcf-validator
conda install -c bioconda vcf-validator

# Validate your VCF file
vcf-validator --input yourfile.vcf --report text
```

**PLINK Data Check**
```bash
# Basic data checks
plink --bfile yourdata --missing --out missing_report
plink --bfile yourdata --freq --out frequency_report
plink --bfile yourdata --hardy --out hwe_report
```

**Genome Build Detection**
```python
# Python script to detect genome build
def detect_build(bim_file):
    """Detect genome build from PLINK .bim file"""
    with open(bim_file) as f:
        for line in f:
            fields = line.strip().split()
            chr_pos = int(fields[3])
            if fields[0] == '1' and chr_pos > 248000000:
                return 'GRCh38/hg38'
            elif fields[0] == '1' and chr_pos > 247000000:
                return 'GRCh37/hg19'
    return 'Unknown'
```

### Quality Control Scripts

**Pre-imputation QC Pipeline**
```bash
#!/bin/bash
# Comprehensive QC pipeline

# Input files
BFILE="input_data"
OUT="qc_data"

# Step 1: Remove variants with high missingness
plink --bfile $BFILE \
      --geno 0.05 \
      --make-bed \
      --out ${OUT}_step1

# Step 2: Remove individuals with high missingness  
plink --bfile ${OUT}_step1 \
      --mind 0.05 \
      --make-bed \
      --out ${OUT}_step2

# Step 3: Remove variants with low MAF
plink --bfile ${OUT}_step2 \
      --maf 0.01 \
      --make-bed \
      --out ${OUT}_step3

# Step 4: Remove variants failing HWE
plink --bfile ${OUT}_step3 \
      --hwe 1e-6 \
      --make-bed \
      --out ${OUT}_final

echo "QC complete. Output: ${OUT}_final"
```

## File Format Converters

::: warning Important
Our service only accepts VCF format. Use these tools to convert from other formats before uploading.
:::

### PLINK to VCF (Required for PLINK users)
```bash
# Convert PLINK to VCF
plink --bfile input_data \
      --recode vcf \
      --out converted_data

# With compression
plink --bfile input_data \
      --recode vcf bgz \
      --out converted_data
```

### 23andMe to VCF (Two-step process)
```python
# Python script to convert 23andMe to VCF via PLINK
import pandas as pd

def convert_23andme_to_plink(input_file, output_prefix):
    """Convert 23andMe format to PLINK, then use PLINK to create VCF"""
    
    # Read 23andMe file
    df = pd.read_csv(input_file, sep='\t', comment='#',
                     names=['rsid', 'chromosome', 'position', 'genotype'])
    
    # Filter out missing genotypes
    df = df[df['genotype'] != '--']
    
    # Create .bim file (variant information)
    bim = df[['chromosome', 'rsid', 'position']].copy()
    bim['morgan'] = 0
    bim['allele1'] = df['genotype'].str[0]
    bim['allele2'] = df['genotype'].str[1]
    bim = bim[['chromosome', 'rsid', 'morgan', 'position', 'allele1', 'allele2']]
    bim.to_csv(f'{output_prefix}.bim', sep='\t', header=False, index=False)
    
    # Create .fam file (sample information)
    fam = pd.DataFrame({
        'family_id': ['SAMPLE1'],
        'individual_id': ['SAMPLE1'],
        'paternal_id': [0],
        'maternal_id': [0],
        'sex': [0],  # 0=unknown, 1=male, 2=female
        'phenotype': [-9]  # -9=missing
    })
    fam.to_csv(f'{output_prefix}.fam', sep='\t', header=False, index=False)
    
    print(f"Step 1 complete. Now convert to VCF:")
    print(f"plink --file {output_prefix} --make-bed --out {output_prefix}")
    print(f"plink --bfile {output_prefix} --recode vcf --out {output_prefix}")

# Usage
convert_23andme_to_plink('genome_data.txt', 'converted_data')
```

## Reference Data

### Genome Reference Files

**Human Reference Genomes:**
- **GRCh37/hg19**: [UCSC Downloads](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/)
- **GRCh38/hg38**: [UCSC Downloads](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/)

**dbSNP Reference:**
- **Build 151**: [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/snp/)
- **Common variants**: [dbSNP Common](https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/)

### Allele Frequency References

**1000 Genomes Allele Frequencies:**
```bash
# Download 1000G allele frequencies
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz

# Extract frequency information
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' \
  ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz > 1000G_frequencies.txt
```

### Strand Files

Download strand files for common genotyping arrays:

**Illumina Arrays:**
- [Global Screening Array](https://www.well.ox.ac.uk/~wrayner/strand/)
- [MEGA Array](https://www.well.ox.ac.uk/~wrayner/strand/)
- [Omni Arrays](https://www.well.ox.ac.uk/~wrayner/strand/)

**Affymetrix Arrays:**
- [Axiom Arrays](https://www.well.ox.ac.uk/~wrayner/strand/)
- [SNP Array 6.0](https://www.well.ox.ac.uk/~wrayner/strand/)

## Software Tools

### Essential Software

**PLINK 1.9/2.0**
```bash
# Install PLINK 1.9
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip
unzip plink_linux_x86_64_20210606.zip

# Install PLINK 2.0 (development version)
wget https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64_20210505.zip
```

**BCFtools/SAMtools**
```bash
# Install via conda
conda install -c bioconda bcftools samtools htslib

# Or compile from source
wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2
tar -xjf bcftools-1.14.tar.bz2
cd bcftools-1.14
make
```

**VCFtools**
```bash
# Install via conda
conda install -c bioconda vcftools

# Basic usage
vcftools --vcf input.vcf --freq --out frequency_output
vcftools --vcf input.vcf --missing-indv --out missing_output
```

### Specialized Tools

**SHAPEIT4 (Phasing)**
```bash
# Download SHAPEIT4
wget https://github.com/odelaneau/shapeit4/releases/download/v4.2.2/shapeit4.2.2.tar.gz

# Basic usage
shapeit4 --input input.vcf.gz \
         --map genetic_map.txt \
         --region 1 \
         --output phased_output.vcf.gz
```

**Eagle2 (Phasing)**
```bash
# Download Eagle
wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz

# Basic usage
eagle --vcf input.vcf.gz \
      --geneticMapFile genetic_map.txt \
      --outPrefix phased_output
```

## Documentation and Tutorials

### Academic Papers

**Imputation Methods:**
1. Das S, et al. "Next-generation genotype imputation service and methods." *Nature Genetics* 48, 1284–1287 (2016).
2. Browning SR, Browning BL. "Rapid and accurate haplotype phasing and missing-data inference for whole-genome association studies by use of localized haplotype clustering." *American Journal of Human Genetics* 81, 1084–1097 (2007).

**Quality Control:**
1. Anderson CA, et al. "Data quality control in genetic case-control association studies." *Nature Protocols* 5, 1564–1573 (2010).
2. Turner S, et al. "Quality control procedures for genome‐wide association studies." *Current Protocols in Human Genetics* 68, 1–19 (2011).

### Online Tutorials

**PLINK Tutorials:**
- [PLINK 1.9 Documentation](https://www.cog-genomics.org/plink/1.9/)
- [PLINK 2.0 Documentation](https://www.cog-genomics.org/plink/2.0/)

**Imputation Tutorials:**
- [Michigan Imputation Server Tutorial](https://imputationserver.readthedocs.io/)
- [Sanger Imputation Service](https://imputation.sanger.ac.uk/)

**Population Genetics:**
- [PLINK/SEQ Tutorial](https://zzz.bwh.harvard.edu/plink/tutorial.shtml)
- [GWAS QC Tutorial](https://www.cog-genomics.org/plink/1.9/tutorials)

## API Documentation

### REST API Reference

**Authentication:**
```bash
# Get access token
curl -X POST https://impute.afrigen-d.org/api/auth/login \
     -H "Content-Type: application/json" \
     -d '{"username": "your_username", "password": "your_password"}'
```

**Submit Imputation Job:**
```bash
# Submit job via API
curl -X POST https://impute.afrigen-d.org/api/jobs/imputation \
     -H "Authorization: Bearer YOUR_TOKEN" \
     -F "genotypes=@data.vcf" \
     -F "reference_panel=hrc" \
     -F "population=eur"
```

**Check Job Status:**
```bash
# Get job status
curl -X GET https://impute.afrigen-d.org/api/jobs/JOB_ID \
     -H "Authorization: Bearer YOUR_TOKEN"
```

### Python SDK

```python
# Install SDK
pip install imputation-server-sdk

# Basic usage
from imputation_server import ImputationClient

client = ImputationClient(api_key='YOUR_API_KEY')

# Submit job
job = client.submit_imputation(
    genotype_file='data.vcf',
    reference_panel='hrc',
    population='eur'
)

# Monitor progress
status = client.get_job_status(job.id)

# Download results
results = client.download_results(job.id)
```

## Community Resources

### User Forums

- **Biostars**: [Imputation Questions](https://www.biostars.org/t/imputation/)
- **Reddit r/genomics**: [Genomics Community](https://www.reddit.com/r/genomics/)
- **Stack Overflow**: [Bioinformatics Tag](https://stackoverflow.com/questions/tagged/bioinformatics)

### Mailing Lists

- **PLINK Users**: plink-users@googlegroups.com
- **Bioconductor**: bioconductor@r-project.org
- **ASHG Forums**: [Human Genetics Forums](https://www.ashg.org/)

### Training Resources

**Online Courses:**
- [Coursera Genomic Data Science](https://www.coursera.org/specializations/genomic-data-science)
- [edX Introduction to Genomics](https://www.edx.org/course/introduction-to-genomics)

**Workshops:**
- **ASHG Annual Meeting**: Annual conference with imputation workshops
- **ESHG Conference**: European human genetics meeting
- **Local Bioinformatics Groups**: Check your institution

## Troubleshooting Guides

### Common Error Messages

**"Chromosome mismatch"**
```bash
# Solution: Check chromosome encoding
grep "^#CHROM" input.vcf  # Should show proper chromosome format
```

**"Build mismatch"**
```bash
# Solution: Liftover coordinates
CrossMap.py vcf hg19ToHg38.over.chain input_hg19.vcf hg38.fa output_hg38.vcf
```

**"Memory error"**
```bash
# Solution: Increase memory allocation
plink --bfile large_dataset --memory 16384 --make-bed --out output
```

### Performance Optimization

**Large Dataset Handling:**
```bash
# Split by chromosome for parallel processing
for chr in {1..22}; do
    plink --bfile large_dataset \
          --chr $chr \
          --make-bed \
          --out chr${chr}_data &
done
wait
```

**Memory Management:**
```bash
# Use streaming for very large files
bcftools view -r 1:1000000-2000000 large_file.vcf.gz | \
  plink --vcf /dev/stdin --make-bed --out region_output
```

## Version History and Updates

### Current Versions
- **Documentation**: v2.0 (December 2024)
- **API**: v2.1 (November 2024)
- **Web Interface**: v3.0 (October 2024)

### Update Notifications
Subscribe to our [mailing list](contact.md) for notifications about:
- New reference panels
- Software updates
- Documentation improvements
- Maintenance schedules

### Contributing
We welcome contributions to our documentation and tools:
1. **Report issues**: [GitHub Issues](https://github.com/imputationserver/docs/issues)
2. **Suggest improvements**: [Feature Requests](https://github.com/imputationserver/docs/discussions)
3. **Submit corrections**: [Pull Requests](https://github.com/imputationserver/docs/pulls)

## Getting Additional Help

Need more assistance? Try these resources in order:

1. **Search this documentation**: Use the search function above
2. **Check the FAQ**: [Imputation FAQ](genotype-imputation/faq.md)
3. **Review examples**: Look for similar use cases in tutorials
4. **Ask the community**: Post on relevant forums
5. **Contact support**: [Get in touch](contact.md) with our team

For urgent issues or institutional partnerships, contact us directly at [support@bioinformaticsinstitute.africa](mailto:support@bioinformaticsinstitute.africa). 