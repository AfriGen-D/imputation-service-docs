# Reference Panels

Reference panels are collections of phased haplotypes from sequenced genomes used for genotype imputation. The AfriGen-D Imputation Service provides population-optimized reference panels to ensure high-quality imputation for diverse populations, with particular emphasis on African and African-diaspora populations.

## Available Reference Panels

The AfriGen-D Imputation Service offers two reference panels optimized for imputation:

### H3Africa v6 (Build 38)

**Overview:**
The H3Africa v6 reference panel is specifically designed for African and African-diaspora populations, providing the most comprehensive African genomic diversity available for imputation.

**Specifications:**

- **Samples**: 4,447 individuals (8,894 haplotypes)
- **Variants**: 130,028,596 biallelic SNPs
- **Coverage**: Autosomes (chromosomes 1-22)
- **Genome Build**: GRCh38/hg38 (Build 38)
- **Population Focus**: 50% African ancestry across diverse populations

**Population Composition:**

The panel includes **48 populations** from:
- **22 African nations** (Nigeria: 727 samples, plus populations from across the continent)
- **6 European countries**
- **12 Asian regions**
- **7 American populations**
- **Papua New Guinea**

**Key Features:**

- **High-coverage sequencing** for accurate variant calls
- **Largest African genomic diversity** available in public reference panels
- **Optimized for admixed populations** with African ancestry
- **Contemporary African populations** reflecting modern genetic diversity

**Recommended For:**

- African and African-diaspora populations
- Admixed populations with African ancestry (African American, Afro-Caribbean, Afro-Latino)
- Studies focusing on African genetic variation
- Research requiring maximum African ancestry representation

---

### 1000 Genomes Project Phase 3 (Version 5)

**Overview:**
A global reference panel with representation from 26 populations worldwide, providing broad ancestral diversity for imputation across multiple populations.

**Specifications:**

- **Samples**: 2,504 individuals (5,008 haplotypes)
- **Variants**:
  - Autosomes: 81,027,987 biallelic SNPs
  - Chromosome X: 3,209,655 biallelic SNPs
- **Coverage**: Autosomes (chromosomes 1-22) and Chromosome X
- **Genome Builds**: GRCh37/hg19 and GRCh38/hg38
- **Populations**: 26 global populations

**Population Breakdown:**

| Super Population | Code | Samples | Description |
|------------------|------|---------|-------------|
| African | AFR | 661 | Sub-Saharan African and African-diaspora |
| European | EUR | 503 | European ancestry |
| East Asian | EAS | 504 | East Asian ancestry |
| South Asian | SAS | 489 | South Asian ancestry |
| Admixed American | AMR | 347 | Latino and admixed American populations |

**African Populations in 1000 Genomes:**

| Code | Population | Samples | Region |
|------|-----------|---------|--------|
| YRI | Yoruba in Ibadan, Nigeria | 108 | West Africa |
| LWK | Luhya in Webuye, Kenya | 99 | East Africa |
| GWD | Gambian in Western Divisions | 113 | West Africa |
| MSL | Mende in Sierra Leone | 85 | West Africa |
| ESN | Esan in Nigeria | 99 | West Africa |
| ASW | African Ancestry in Southwest US | 61 | African American |
| ACB | African Caribbean in Barbados | 96 | Afro-Caribbean |

**Total African Samples**: 661 individuals

**Key Features:**

- **Global diversity** across five continental ancestries
- **Chromosome X support** for sex-linked variant imputation
- **Well-characterized populations** with extensive validation
- **Multiple genome builds** for flexibility

**Recommended For:**

- Multi-ancestry and admixed populations
- Chromosome X imputation needs
- Global comparative studies
- When H3Africa panel is not accessible

---

## Choosing the Right Reference Panel

### For African and African-Diaspora Populations

**Recommended**: **H3Africa v6**

The H3Africa v6 panel provides:
- **Nearly 7x more samples** than African samples in 1000 Genomes (4,447 vs 661)
- **More comprehensive African diversity** across 22 nations
- **Higher variant density** optimized for African populations
- **Better imputation accuracy** for African-specific variants

### For Other Populations

**Recommended**: **1000 Genomes Phase 3 (Version 5)**

The 1000 Genomes panel is suitable for:
- Non-African populations (European, Asian, South Asian, Latino)
- Multi-ancestry studies requiring global representation
- Studies needing Chromosome X imputation
- Projects with broad ancestral diversity

### For Admixed Populations

**African Ancestry Component**:
- If **≥50% African ancestry**: Use **H3Africa v6**
- If **<50% African ancestry**: Use **1000 Genomes Phase 3**

**Multiple Ancestry Components**:
- For complex admixture patterns, **1000 Genomes Phase 3** provides balanced representation

---

## Imputation Quality Expectations

### Variant Frequency Categories

Reference panel performance varies by variant frequency:

| Frequency Category | MAF Range | Expected R² (Quality) |
|-------------------|-----------|----------------------|
| Common | >5% | >0.9 (Excellent) |
| Low Frequency | 1-5% | 0.7-0.9 (Good) |
| Rare | 0.5-1% | 0.5-0.7 (Moderate) |
| Very Rare | <0.5% | <0.5 (Poor) |

::: tip Understanding R² (Imputation Quality)
- **R² > 0.8**: High-quality imputation, suitable for association analysis
- **R² 0.5-0.8**: Moderate quality, may be suitable depending on study design
- **R² < 0.5**: Low quality, generally not recommended for analysis
:::

### Factors Affecting Imputation Quality

1. **Sample Size of Reference Panel**: Larger panels improve rare variant imputation
2. **Population Match**: Better match = higher accuracy
3. **Marker Density**: Higher density input data improves results
4. **Genome Build**: Must match between input data and reference panel

---

## Technical Details

### Software and Algorithms

The AfriGen-D Imputation Service uses:
- **Minimac4**: State-of-the-art imputation algorithm
- **Eagle2**: Fast and accurate pre-phasing

### Processing Pipeline

1. **Quality Control**: Variant validation and filtering
2. **Pre-phasing** (if needed): Haplotype estimation using Eagle2
3. **Imputation**: Genotype prediction using Minimac4
4. **Quality Metrics**: R² and Info score calculation

### Output Information

Imputed genotypes include:
- **Dosage Format**: Continuous values (0-2) representing expected allele counts
- **Genotype Probabilities**: Likelihood of each genotype (0/0, 0/1, 1/1)
- **Quality Scores**: R² (accuracy) and Info (information content)

---

## Frequently Asked Questions

### Can I use multiple reference panels?

No, the service processes each job with a single selected reference panel. Choose the panel that best matches your study population.

### What if my population isn't represented?

For populations not well-represented in either panel:
- **African ancestry**: Use H3Africa v6 (broad African diversity)
- **Non-African ancestry**: Use 1000 Genomes Phase 3 (global diversity)

### How do I access the H3Africa v6 panel?

The H3Africa v6 panel is available directly through the AfriGen-D Imputation Service. Select it when submitting your imputation job at [impute.afrigen-d.org](https://impute.afrigen-d.org).

### Can I impute with both panels and compare?

While you cannot run both simultaneously, you can submit separate jobs with different panels and compare results. This can be useful for validation or method comparison studies.

---

## Additional Resources

### Reference Panel Publications

**H3Africa Reference Genome:**
- H3Africa Consortium. "Enabling the genomic revolution in Africa." Science 344.6190 (2014): 1346-1348.

**1000 Genomes Project:**
- 1000 Genomes Project Consortium. "A global reference for human genetic variation." Nature 526.7571 (2015): 68-74.

### H3ABioNet Resources

- **H3ABioNet Website**: [h3abionet.org](https://www.h3abionet.org)
- **Imputation Service**: [impute.h3africa.org](https://impute.h3africa.org)
- **Training Materials**: Available through H3ABioNet workshops

### Related AfriGen-D Resources

- **AGVD** ([agvd.afrigen-d.org](https://agvd.afrigen-d.org)): African Genome Variation Database
- **AGMP** ([agmp.afrigen-d.org](https://agmp.afrigen-d.org)): African Genomic Medicine Portal
- **AfPO** ([afpo.afrigen-d.org](https://afpo.afrigen-d.org)): African Population Ontology

---

::: info Need Help Choosing?
If you're unsure which reference panel to use, consult our [FAQ](genotype-imputation/faq.md) or [contact our support team](contact.md) for guidance.
:::

**Ready to submit your imputation job?** Visit [Getting Started](genotype-imputation/getting-started.md) for step-by-step instructions.
