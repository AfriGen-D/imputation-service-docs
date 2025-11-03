# AfriGen-D Imputation Service

## Free Genotype Imputation Platform for African Genomics Research

The AfriGen-D Imputation Service provides free, high-quality genotype imputation optimized for African and African-diaspora populations. Part of the [African Genomics Data Hub (AfriGen-D)](https://afrigen-d.org), our service uses state-of-the-art algorithms (Minimac4) and population-optimized reference panels to support genomics research throughout the African data lifecycle.

You can upload phased or unphased GWAS genotypes and receive phased and imputed genomes in return. For all uploaded datasets, comprehensive quality control is performed.

::: tip Citation
If you use the AfriGen-D Imputation Service, please acknowledge:

> AfriGen-D Imputation Service. African Genomics Data Hub (AfriGen-D).
> Available at: https://impute.afrigen-d.org (2025).
> Supported by NIH Grant Number U24HG012750.

For the underlying imputation methodology, also cite:

> Das S, Forer L, Schönherr S, Sidore C, Locke AE, Kwong A, Vrieze S, Chew EY, Levy S, McGue M, Schlessinger D, Stambolian D, Loh PR, Iacono WG, Swaroop A, Scott LJ, Cucca F, Kronenberg F, Boehnke M, Abecasis GR, Fuchsberger C. Next-generation genotype imputation service and methods. Nature Genetics 48, 1284–1287 (2016).
:::

---

## Quick Start

Ready to impute your genomics data? Follow these steps:

1. **[Prepare your data](genotype-imputation/data-preparation.md)** - Format your genotype data (bgzip-compressed VCF files)
2. **[Choose reference panel](reference-panels.md)** - Select appropriate population-specific panel
3. **[Upload and submit](genotype-imputation/getting-started.md)** - Submit your imputation job
4. **[Download results](genotype-imputation/getting-started.md#downloading-results)** - Get your imputed genotypes with quality metrics

## Key Features

### Comprehensive Imputation Platform

- **Free Service**: No cost for academic and research use worldwide
- **Multiple Reference Panels**: Choose from global and population-specific panels
- **Quality Control**: Comprehensive QC performed on all datasets
- **Secure Processing**: Encrypted data handling and automatic deletion after 7 days
- **API Access**: Programmatic access for large-scale studies

### Optimized for Diverse Populations

- **Population-Specific QC**: Quality control optimized for genetic diversity
- **Multiple Reference Panels**: Global coverage including African populations
- **Ancestry-Aware Processing**: Handles complex admixture patterns

### Part of AfriGen-D Ecosystem

- **Integrated Resources**: Connect with AGMP, AGVD, AfPO, and Data Catalogue
- **Training & Support**: Workshops and tutorials through H3ABioNet
- **Community**: H3Africa Consortium and Pan-African bioinformatics network

## Supported Formats

- **Input**: bgzip-compressed VCF files (`.vcf.gz`), one file per chromosome
- **Output**: Phased, imputed VCF with dosage information and quality scores
- **Chromosomes**: Autosomes (1-22), Chromosome X, and HLA region
- **Sample Size**: Up to 110,000 samples
- **Genome Builds**: GRCh37/hg19 and GRCh38/hg38

## About AfriGen-D

The African Genomics Data Hub (AfriGen-D) is dedicated to enabling innovation in African genomics research. Our imputation service is part of a comprehensive ecosystem:

### Connected Resources

- **[AGMP](https://agmp.afrigen-d.org)** - African Genomic Medicine Portal for disease associations
- **[AGVD](https://agvd.afrigen-d.org)** - African Genome Variation Database
- **[AfPO](https://afpo.afrigen-d.org)** - African Population Ontology for standardized terminology
- **[Data Catalogue](https://catalogue.afrigen-d.org)** - Discovery of African genomics datasets

### Our Mission

We facilitate African genomics research throughout the data lifecycle:

- **Data Generation**: Standards and best practices
- **Discovery**: Catalogue and search African datasets
- **Analysis**: Imputation and variant annotation tools
- **Archiving**: Long-term data preservation
- **Sharing**: Responsible data sharing frameworks

### Support & Collaboration

AfriGen-D is funded by NIH Grant Number U24HG012750 and operates in collaboration with:

- H3ABioNet (Pan-African Bioinformatics Network)
- H3Africa Consortium
- Multiple African research institutions

::: info Learn More
Visit [afrigen-d.org](https://afrigen-d.org) to explore all AfriGen-D resources and services.
:::

## Getting Help

### Documentation & Support

- **Comprehensive [FAQ](genotype-imputation/faq.md)** for common questions
- **[Data Preparation Guide](genotype-imputation/data-preparation.md)** for format requirements
- **[Pipeline Overview](genotype-imputation/pipeline-overview.md)** for technical details
- **[Contact Support](contact.md)** for technical assistance

### Training Opportunities

AfriGen-D offers regular training workshops through H3ABioNet:

- Genotype imputation best practices
- Reference panel selection guidance
- Quality control and results interpretation
- Downstream analysis workflows

### Community Resources

- **H3ABioNet**: Pan-African bioinformatics network
- **H3Africa Consortium**: African genomics research community
- **AfriGen-D Training Portal**: Workshops and webinars

---

**Ready to get started?** Visit [impute.afrigen-d.org](https://impute.afrigen-d.org) to submit your first imputation job! 