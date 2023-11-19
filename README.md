## Polygenic Risk Score Analysis (PRS)
### Contents
1. [Acknowledgement/source](#source)
2. [Pre-requisite to follow the tutrorial](#prerequisite)
3. [Download Dataset](#datasets)
4. [What is polygenic Risk Score](#what-is-polygenic-risk-score?)
5. [Terminologies used in Polygenic Risk Score  analysis](#terminology)
6. [QC of base data](#qc-of-base-data)
   - [Heritability Check](#heritability-check)
   - [Effect Allele](#effect-allele)
   - [Genome Build verification](#genome-build)
   - [GWAS related QC](#standard-gwas-qc)
   - [Identify ambiguity in SNPs](#ambigous-snps)
   - [Identify Mismatching SNPs](#mismatching-snps)
   - [Identify and remove duplicate SNPs](#duplicate-snps)
   - [Sex verification of samples](#sex-chromosomes)
   - [Identify overlapping samples between Base and target data](#sample-overlap)
   - [Identify Relatednes between individuals](#relatedness)
7. [QC of Target Data](#qc-of-target-data)
   - [Genome Build](#genome-build-of-target-data)
   - [Standard GWAS QC](#standard-gwas-qc-of-target-data)
   - - [Removing highly correlated SNPs](#step1-removing-highly-correlated-snps)
   - - [Heterozygosity rates](#step2-heterozygosity-rates-can-then-be-computed-using-plink)
   - [Mismatching SNPs](#mismatching-snps-in-target-data)
   - [Duplicate SNPs](#duplicate-snps-in-target-data)
   - [Sex Chromosome](#sex-chromosomes-in-target-data)
   - [Sample Overlap in Target Data](#sample-overlap-in-target-data)
   - [Sample Relatedness](#relatedness-of-individuals-target-data)
8. [Generate final QCed Target file](#generate-final-qc-target-file)
9. [PRS analysis](#calculating-and-analysis-of-prs)
   - [PRS with PLINK](#plink)
   - - [Dataset used in PLINK](#required-dataset)
   - - [Effect size](#effect-size-standardisation)
   - - [Clumping](#clumping)
   - - [Generating PRS](#generating-prs-with-plink)
   - - [Accounting for Population Stratification](#accounting-for-population-stratification)
   - - [Finding "best-fit" PRS](#finding-best-fit-prs)






### Acknowledgement/Source
The code here is based on the following 2 sources
 1. [Choi, S.W., Mak, T.SH. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc 15, 2759–2772 (2020)](https://www.nature.com/articles/s41596-020-0353-1)
 2. https://choishingwan.github.io/PRS-Tutorial/

### Prerequisite
The tutorial require basic knowledge of following 
 1. Basic Linux command (awk)
 2. R programming
 3. [Plink](https://www.cog-genomics.org/plink2)

### Datasets
We will use the dataset listed in the tutorial but the code can be extended for custom dataset.  Here we will have 2 data sets. The nature of these data is explained in terminology section below.
 1. [Base data](https://drive.google.com/file/d/1hUgPewcOs4ogWJIzD7UrCUectD0Aoc1J/view?usp=sharing). 
 2. [Target data](https://drive.google.com/file/d/10fWD4z_hbE3YXIdRPs1mvbGlPcpIfdQK/view?usp=sharing): Simulated data based on the 1000 Genomes Project European samples

### What is Polygenic Risk Score?
The definition you provided succinctly explains what polygenic risk scores (PRSs) are and how they are calculated. Let's break it down:

1. **PRSs (Polygenic Scores):** PRSs are numerical estimates that represent an individual's genetic predisposition or liability to a specific phenotype or trait. These traits can include various characteristics, such as disease susceptibility, height, cognitive abilities, and more.

2. **Calculation:** PRSs are calculated by considering an individual's entire genome. This calculation involves summing the genetic information (genotypes) across a vast number of single nucleotide polymorphisms (SNPs) or genetic variants present in the individual's DNA.

3. **Weighted Sum:** Each genotype is assigned a weight, which is determined by the corresponding genotype effect size estimates. These effect size estimates are derived from large-scale genome-wide association studies (GWAS), which identify genetic variants associated with the phenotype of interest. Effect sizes indicate how strongly each SNP is associated with the trait.

4. **Genome-wide Genotypes:** The PRS calculation takes into account all the genotypes (minor allele frequency > 0.01) biallelic SNPs) across the genome, not just a select few genetic variants. This comprehensive approach allows for a more holistic assessment of an individual's genetic predisposition to the phenotype.

In summary, PRSs provide a single numeric value that reflects an individual's genetic liability for a particular trait or phenotype. They are computed by considering the genetic information from across the entire genome and weighting each genetic variant's contribution based on its effect size estimate derived from GWAS data. PRSs are valuable tools in genetics research and are used to predict disease risk, assess genetic contributions to traits, and explore the genetic basis of complex phenotypes.

### Note:
*PRS analysis is similar to machine learning methods in making predictions and hence methods like feature slecetion (SNPs here), parameter optimisation, shrinkage, cross validation  etc are integral part of it.*


### Terminology
Few Terminologies associated with polygenic risk score (PRS)
 1. **Risk Allele**: The risk allele is the allele at a SNP that is associated with an increased risk of developing a particular disease or trait.
 2. **Effect Size**:Each SNP included in a PRS has an associated effect size, which quantifies how strongly that SNP is associated with the trait or disease of interest. This effect size is typically represented as a beta coefficient (β) for continuous traits (blood pressure, height) or an odds ratio (OR) for binary traits (disease). A larger absolute value of the effect size indicates a stronger association between the SNP and the trait or disease.

    * **Beta Coefficient (β)**: The beta coefficient, often denoted as β, is a statistical measure used to quantify the effect size of a genetic variant (like a SNP) on a continuous trait or outcome. It indicates the change in the mean value of the trait for each additional copy of the risk allele. A positive beta coefficient indicates that each additional risk allele is associated with an increase in the trait value, while a negative beta coefficient suggests a decrease in the trait value for each additional risk allele. For example, if β = 0.2, it means that, on average, the trait value increases by 0.2 units for each additional copy of the risk allele.
    
    * **Odds Ratio (OR)**: The odds ratio, denoted as OR, is a statistical measure used to quantify the effect size of a genetic variant on a binary outcome, such as disease risk (affected vs. unaffected). It represents the odds of developing the disease in individuals with the risk allele compared to those without it. An OR greater than 1 indicates an increased risk associated with the risk allele, while an OR less than 1 indicates a decreased risk.
 3. **Summary statistic**:In Genome-Wide Association Studies (GWAS), summary statistics are key results that provide a compact and informative representation of the associations between genetic variants (typically for each single nucleotide polymorphisms or SNPs) and a particular trait or disease across the entire genome. These statistics summarize the data obtained from genotyping or sequencing a large number of individuals and are used to identify genetic variants associated with the trait or disease of interest. Here are some common summary statistics used in GWAS:

    * **p-value**: The p-value assesses the statistical significance of the association between a genetic variant and the trait or disease. A lower p-value suggests a stronger and more reliable association. In GWAS, researchers often use a threshold (e.g., p < 5 x 10^-8) to identify significant associations while accounting for multiple testing.

    * **Effect Size** Explained above.

    * **Allele Frequency**: This statistic indicates the frequency of the risk allele (the allele associated with increased risk) in the population being studied. Understanding allele frequencies is important for assessing the potential impact of the associated variant on disease risk.

    * **Standard Error**: The standard error provides a measure of the precision of the effect size estimate. Smaller standard errors indicate more precise estimates.

    * **Confidence Intervals (CI)**: Confidence intervals provide a range of effect size values values within which the true effect size is likely to fall with a certain level of confidence. They help in assessing the precision of the effect size estimate.

    * **Minor Allele Frequency (MAF)**: MAF is the frequency of the less common allele in the population. Variants with low MAF may have limited statistical power to detect associations in GWAS.

    * **Linkage Disequilibrium (LD) Metrics**: LD measures, such as r-squared (r²) or D', describe the correlation between the tested SNP and nearby genetic variants. High LD can help identify regions of the genome associated with the trait, even if the directly tested SNP does not show a strong association.

    * **Genomic Control**: Genomic control is a quality control measure used to correct for any inflation in the test statistics due to population structure or other confounding factors. It provides an inflation factor (λ) to adjust the p-values accordingly.

   4. **Base data** The GWAS summary statistics that is used in the PRS analysis.  Base trait will be the phenotype that is studied in the PRS analysis.*"Base data (GWAS)consisting of summary statistics (e.g., betas and P values) of genotype-phenotype associations at genetic variants (hereafter SNPs) genome-wide, typically made available online in text format by the investigators who performed the GWAS"*-Choi et al.

   5. **Target Data** The genotype phenotype data that in combination with base data PRS's are calculated.  Generally the data is in the binary format e.g. PLINK files. *"consisting of genotypes, and usually also phenotype (s), in individuals from a sample to which the researchers performing the PRS analysis have access (often not publicly available), which should be independent of the GWAS sample"*-Choi et al.

  6. **Shrinkage** Shrinkage is a method used to correct the overestimated effect sizes by pulling them closer to a more realistic estimate. It does this by essentially reducing the magnitude of the estimated effects. **Why Shrinkage is Used:** Shrinkage is employed because overestimated effect sizes can lead to overly optimistic or misleading conclusions. By applying shrinkage, researchers aim to obtain more accurate and reliable estimates that better reflect the actual relationships in the broader population, helping to make more robust and generalizable findings.
  7. **SNP heritability**  is a concept used in genetics to estimate the proportion of the variation in a trait or phenotype that can be attributed to common genetic variants, specifically single nucleotide polymorphisms (SNPs). It provides insights into how much of the observed variability in a particular trait within a population can be explained by differences in genetic variation at the level of SNPs.
    
     * **Phenotypic Variance (Variance of the Trait)**: Any observable trait or phenotype, such as height, disease risk, or cognitive ability, can vary among individuals within a population. This variation is referred to as the phenotypic variance, denoted as "Vp."

     * **Genetic Variance (Variance Attributable to Genetic Factors)**: Part of the phenotypic variance can be attributed to genetic factors. This genetic variance is due to the effects of genetic variants, including SNPs, that influence the trait. It is referred to as the genetic variance and is denoted as "Vg."
     
     * **SNP Heritability (h²)**: SNP heritability, represented as "h²," is the proportion of the total phenotypic variance (Vp) that is explained by the additive genetic effects of common SNPs. It is expressed as a ratio, where h² = Vg / Vp.

        If h² = 0.2 (20%), it means that approximately 20% of the observed variability in the trait within the population is due to the genetic effects of common SNPs.

        If h² = 0.6 (60%), it suggests that a larger proportion (60%) of the phenotypic variance can be attributed to genetic factors represented by SNPs.

        * The h² is affected by accuracy in measuring effect size and differences in base and target populations. The predictive power of PRS will be lower than h² if there are errors in effect size measurements and diversity between population in base and taget data. As GWAS sample size increases there is substantial improvement in effect size measurements.

  8. **Overfitting** Overfitting can be a concern when developing polygenic risk scores (PRS), just as it is in other statistical modeling contexts. Overfitting in the context of PRS development can lead to overly optimistic estimates of the predictive power of the score, which may not generalize well to new populations or datasets.  This could be due to too many parameters.
  9. **Variance explained** : Variance explained" is a statistical concept that quantifies the proportion of the total variability or dispersion in a phenotype that can be attributed to a specific factor or predictors (SNP). In PRS predictive modelling the effect is assumed as linear realationship.

<br>

</br>

fig1 : *Workflow for PRS analysis (source:Choi et al 2020* 

<img src="./workflow.png" width="800"/>

<br>

</br>

## QC of Base data

The base data file `HeightHeight.gwas.txt.gz` is compressed and to visualise the contents we will use `gunzip` command the uncompress the file and we will pipe it to head command that display top 10 rows of the file. 

```{shell}
cd data

gunzip HeightHeight.gwas.txt.gz | head -10

# Output
CHR	BP	SNP	A1	A2	N	SE	P	OR	INFO	MAF
1	756604	rs3131962	A	G	388028	0.00301666	0.483171	0.997886915712657	0.890557941364774	0.369389592764921
1	768448	rs12562034	A	G	388028	0.00329472	0.834808	1.00068731609353	0.895893511351165	0.336845754096289
1	779322	rs4040617	G	A	388028	0.00303344	0.42897	0.997603556067569	0.897508290615237	0.377368010940814
1	801536	rs79373928	G	T	388028	0.00841324	0.808999	1.00203569922793	0.908962856432993	0.483212245374095
1	808631	rs11240779	G	A	388028	0.00242821	0.590265	1.00130832511154	0.893212523690488	0.450409558999587
1	809876	rs57181708	G	A	388028	0.00336785	0.71475	1.00123165786833	0.923557624081969	0.499743932656759
1	835499	rs4422948	G	A	388028	0.0023758	0.710884	0.999119752645202	0.906437735120596	0.481016005816168
1	838555	rs4970383	A	C	388028	0.00235773	0.150993	0.996619945289758	0.907716506801574	0.327164029672754
1	840753	rs4970382	C	T	388028	0.00207377	0.199967	0.99734567895614	0.914602590137255	0.498936220426316
```
The header is a bit misaligned.  The header values are explained below.

* **CHR**: Reference chromosome in which the SNP resides
* **BP**: Chromosomal location of the SNP on chromosome.
* **SNP**: SNP ID, usually in the form of rs-ID
* **A1**: Allele that is responsible for causing the effect.
* **A2**:  Non-effect allele of the SNP
* **N**: Number of samples used to obtain the effect size estimate
* **SE**: The standard error (SE) of the effect size esimate
* **P**: The P-value of association between the SNP genotypes and the base phenotype
* **OR**: The effect size estimate of the SNP,OR (Odd Ratio) if the outcome is binary/case-control. If the outcome is continuous or treated as continuous then this will usually be BETA
* **INFO**: The imputation information score
* **MAF**: The minor allele frequency (MAF) of the SNP

### Heritability Check
It is recommended that base data have a chip-heritability estimate h²snp >0.05.  What does that mean?

The term "**chip-heritability estimate**" refers to an estimate of the heritability of a trait or phenotype using genetic data obtained from genotyping or sequencing microarray chips, often referred to simply as "chips." These chips are a common technology used in large-scale genetic studies, such as genome-wide association studies (GWAS) or other genetic analyses. Chip-heritability estimates help quantify the genetic contribution to a trait or disease using data from these genotyping chips.

Here's a more detailed explanation:

   * **Heritability**: Heritability is a concept in genetics that quantifies the proportion of the total variation in a trait or phenotype within a population that can be attributed to genetic factors. It reflects the extent to which genetic variation contributes to differences in the trait.

   * **Chip-Heritability Estimate**: Chip-heritability estimates are calculated by analyzing the genetic data obtained from genotyping chips and using statistical methods to assess the heritability of a specific trait or disease. These estimates provide insights into how much of the trait's variation can be explained by the genetic information captured by the chip.

**How to calculate it?**

   * **Calculating LD Scores:** Calculate linkage disequilibrium (LD) scores for pairs of SNPs in the dataset. LD measures the extent to which two genetic variants are correlated or tend to be inherited together. LD scores help quantify the genetic relationship between SNPs in the population.

   * **Heritability Estimation**: Use LD Score Regression to estimate the heritability of height. This involves regressing the observed GWAS summary statistics (effect sizes, standard errors, and p-values for each SNP) against the LD scores. The output of the regression provides an estimate of the heritability of height, which represents the proportion of height variability in the population that can be attributed to genetic factors.
   * **Partitioning Genetic Variance**:LD Score Regression also allows for the partitioning of the estimated heritability into different components. For example:
        * **Common SNPs**: Determine the proportion of heritability explained by common SNPs that are well-captured by the genotyping chip used in the study.
        * **Rare Variants**: Assess the contribution of rare genetic variants, which may have been less comprehensively captured by the chip.
        * **Polygenic Component**: Examine the heritability explained by the combined effects of many common SNPs.

### Effect Allele
It is important to identify the allele that is contributing to phenotype (A1 in the table).  If this is not identified properly then the effect of the PRS in the target data will be in wrong direction.

### Genome Build
Make sure the genome release used in the base and traget data is same.  If that is not true then update the coordinates to recent version.  There are tools to facilitate thise.g [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).

### Standard GWAS QC
Here we filter SNPs from base data based on INFO and MAF field. Both are descroibed below.
*  **INFO** The imputation information score, often referred to as "INFO" or "info score," is a metric used in genetic imputation to assess the quality and reliability of imputed genotypes. Imputation is a process used to predict or infer genotypes at untyped or missing genetic variants (typically single nucleotide polymorphisms or SNPs) based on the patterns of genetic variation observed in a reference panel of genotyped individuals. The INFO score provides an estimate of how confident we can be in the imputed genotypes for a given SNP. Here's what the INFO score indicates:

    * **Quality of Imputation**:The INFO score reflects the quality of the imputation for a particular SNP. It ranges from 0 to 1, with higher values indicating better quality imputation. An INFO score close to 1 suggests high confidence in the imputed genotypes, meaning that the imputation likely accurately reflects the true genotypes. An INFO score closer to 0 indicates lower confidence and may suggest that the imputed genotypes are less reliable.

    * **Accuracy and Reliability**:An INFO score provides insights into the accuracy and reliability of imputed genotypes. Higher INFO scores are associated with more accurate imputations. Researchers often use a threshold INFO score (e.g., INFO > 0.8) to filter or select imputed SNPs for downstream analyses. This threshold helps ensure that only high-quality imputed data is used.

    * **Factors Affecting INFO Score**:The INFO score is influenced by various factors, including the quality and size of the reference panel, the density of genotyped markers, and the LD (Linkage Disequilibrium) patterns in the study population. Imputation tends to perform better when there is a large and diverse reference panel with genotypes closely matching those in the study population.

    * **Use in Genetic Association Studies**:INFO scores are essential in genome-wide association studies (GWAS) and other genetic analyses. They help researchers assess the reliability of imputed genotypes used in association testing.  Higher INFO scores provide greater confidence that genetic variants (SNPs) have been imputed accurately and can be used effectively in genetic analyses to identify associations with traits or diseases.  

*  **MAF** **M**inor **A**llele **F**requencey: MAF stands for Minor Allele Frequency, and it is a fundamental concept in genetics that refers to the frequency at which the less common allele of a genetic variant (usually a single nucleotide polymorphism or SNP) occurs within a population. The MAF provides important information about the genetic variation and diversity present in a population. Here's what MAF means and how it is calculated:
    
    * **Minor Allele**: 
      * Within a population, a genetic variant may have two alleles: a major allele (the more common one) and a minor allele (the less common one).
      * The minor allele is the allele that occurs at a lower frequency in the population.

    * **Minor Allele Frequency (MAF)**: 
      * MAF is a measure of how often the minor allele occurs in a population relative to the total number of alleles at that locus.
      * It is typically expressed as a proportion or percentage, ranging from 0 (no copies of the minor allele) to 0.5 (half the alleles in the population are the minor allele).
      * For example, if a SNP has a MAF of 0.2, it means that 20% of the alleles at that locus in the population are the less common, or minor, allele.

    * **Significance of MAF**: 
      * MAF provides insights into the genetic diversity of a population. High MAF values suggest that a variant is relatively common in the population, while low MAF values indicate that it is rare.
      * Rare variants (low MAF) may have larger effects on traits or diseases when they are associated with them but are more challenging to detect in genetic studies due to their low frequency.
      * Common variants (high MAF) are easier to detect but often have smaller individual effects on traits.

    * **MAF** 
      * In Genetic Studies:In genome-wide association studies (GWAS) and other genetic analyses, MAF is an important consideration. Researchers may filter or prioritize genetic variants based on their MAF.
      * Variants with low MAF may be excluded from analyses if they are too rare to provide meaningful statistical power, while common variants may be included in association tests.

    * **Population-Specific MAF**
      * MAF can vary between different populations or ethnic groups. A variant that is rare in one population may be common in another.

**Other recommemded GWAS QC step**
*"We recommend the following QC criteria for standard analyses: genotyping rate >0.99, sample missingness <0.02, Hardy-Weinberg Equilibrium P >1 × 10−6, heterozygosity within 3 standard deviations of the mean, minor allele frequency (MAF) >1% (MAF >5% if target sample N <1000) and imputation ‘info score’ >0.8. If both the base and target data are large (e.g., N >50,000), then SNPs with MAF <1% may be included, in which case we recommend a minor allele count >100 in the base and target data to ensure the
integrity of normality assumptions implicit in association testing and LD calculation. Future work will be required to integrate the effects of extremely rare and common variants and to establish whether their joint effects are typically additive36.
PLINK is a useful software for performing these, and other, QC procedures"* Choi et al.

Here we will filter SNPs where  MAF < 1% and INFO < 0.8 (with very large base sample sizes these thresholds could be reduced if sensitivity checks indicate reliable results)
```{shell}
# count number of SNPs in base file
wc -l 

#Output
  529505 Height.gwas.txt

# Filter 
awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' Height.gwas.txt > Height.txt

# NR==1 : Prints the header line in output
# $11 > 0.01 : Coloumn 11 (MAF) must be greater than 0.01
# $10 > 0.8  :  Column 10 (INFO) values be greater than 0.8
# Output is in Height.txt 


# Count number of SNPS after filtering.
wc -l Height.txt

#Output
  529496 Height.txt
```

We started with 529504 (first row is header line) and filtered to get 529495.  This filter removed 9 SNPs.

### Ambigous SNPs:
This happens in the cases when the chromosome stand information is lacking regarding SNP(s).   This happens in alleles that have complemantry bases A/T or G/C.  The issue is explained below.
```
        Base Chip (+ strand)               Target Chip(- Strand)

(+ Starnd) ----C------T------    vs  (+ Starnd) ----G------A------
(- Strand) ----G------A------        (- Strand) ----C------T------
```
If in base chip assay (+strand) is used and in Target chip assay (-strand) is used. Both will report C and T , in absence of strand information it will indicate that both SNPs are same, however if we have strand information then we can conclude they are different, base is C and T whereas Target is G and A.  MAF can help to identify which allele is on which strand as long as the population are same between base and target and MAF is not close to 50%.

If you have noticed we donot have strand information in our base data so it is best to filter out ambigous SNPs.
```{shell}
cd data
awk '!( ($4=="A" && $5=="T") || \
      ($4=="T" && $5=="A") || \
      ($4=="G" && $5=="C") || \
      ($4=="C" && $5=="G")) {print}' Height.txt > height_GWASQC.txt 

wc -l height_GWASQC.txt  # count how many SNPs are retained.
# output
499618 height_GWASQC.txt
```
The GWAS QC has dropped 29886 SNPs from the analysis due to ambiguity.

### Mismatching SNPs
This is the case when when we have strand information for both base and target dataset but the strands used in the assay are different. If the alleles are assayed from different strands between base and the target dataset, in most tools this is taken care of by switching both strand and base to match the strand in both files. This will resolve the mismatch and if it cannot be resolved it is removed from analysis.

### Duplicate SNPs
It is critical to not have duplicate SNPs in base or traget dataset as it will either produce errors or crash the PRS application.  The removal of duplicates can be done in linux, R or python as per choice of the user
lets check for duplicate SNPs

```{shell}
awk '{seen[$3]++; if(seen[$3]==1) { print }}' height_GWASQC.txt > height_GWASQC_dedup.txt 

# Count number of SNPs retained
wc -l height_GWASQC_dedup.txt
#output
 499618 height_GWASQC_dedup.txt

```
There were no duplicates as we didnot filter out any SNP.  To be correct there were 2 duplicate SNPs in `height.txt`` as shown below.
```{shell}
awk '{seen[$3]++; if(seen[$3]!=1) { print }}' height.txt  

#Output
3	182254751	rs7622072	C	G	388028	0.00526681	0.126843	0.991991549879677	0.903220366820217	0.472790254192978
5	139208711	rs113309830	C	G	388028	0.00377597	0.542249	0.997701505592967	0.90708730603453	0.436439940493625

```
 Both the SNPs are ambigous SNPs so they got filtered in the previous filteration step. Hence we donot see them here.


### Sex Chromosomes

In standard quality control procedures for Genome-Wide Association Studies (GWAS), it is common practice to exclude individuals when there is a discrepancy between their reported sex and the sex determined by their sex chromosomes. While such discrepancies can arise from variations in sex and gender identity, they could also be indicative of sample mislabeling or incorrect reporting, making the data potentially unreliable. To perform this sex check, software like PLINK can be employed. In this check, individuals are classified as females if their X chromosome homozygosity estimate (F statistic) is less than 0.2 and as males if the estimate is greater than 0.8.

Furthermore, in cases where the primary goal of an analysis is to model autosomal genetics exclusively, it is advisable to remove all genetic variants located on the X and Y chromosomes from both the base and target datasets. This removal eliminates the possibility of non-autosomal sex-related effects influencing the results.

However, it's worth noting that incorporating information from the sex chromosomes has the potential to offer insights into the underlying causes of traits or diseases and can enhance the predictive power of Polygenic Risk Scores (PRSs). Therefore, some researchers may choose to include these chromosomes in their analyses. Nevertheless, when reporting analyses that incorporate sex chromosomes, it is crucial to explicitly highlight how the modeling assumptions may have affected the results.


### Sample Overlap

Overlap of samples between base and target data can inflate PRS-trait associations. The degree of inflation corresponds to the proportion of overlapping individuals in the target data. Ideally, remove overlapping samples from the base data and recalculate GWAS. This allows PRS calculation in all target individuals and enhances power for association testing. In consortium meta-analysis, a practical solution is leave-one-out meta-analysis GWAS. Each study is excluded one by one, creating independent target data. Alternatively, leave-one-out results can be calculated analytically. In complex cases, future methods may address partial or unknown overlap. To minimize overlap risk, choose target samples unlikely to be part of the base sample (e.g., due to age or location). If overlap remains possible, inflation in results can't be ruled out.

### Relatedness
This passage discusses the potential impact of close genetic relatedness between individuals in the base and target data on the association between Polygenic Risk Scores (PRS) and the target phenotype (trait or disease). Here's a breakdown of the key points:

*  **Relatedness and Inflation**: When individuals in the base and target datasets are closely related, it can lead to an inflation of associations between PRS and the target phenotype. In other words, it may falsely appear that PRS is strongly associated with the trait due to the shared genetic and potentially environmental factors.

*  **Population Structure vs. Close Relatives**: Population structure refers to genetic differences among individuals from different ancestral backgrounds. It's a broader issue and requires a general solution. However, when close relatives are included in both datasets, the problem becomes more pronounced because they not only share genetic similarities but might also share the same household environment.

*  **Removing Close Relatives**: To mitigate this issue, it's advisable to identify and remove closely related individuals (e.g., first- or second-degree relatives) that appear in both the base and target datasets. This step helps eliminate the risk of inflated associations caused by shared genetics and environments.

*  **Selecting Diverse Data**: If it's not possible to access genetic data for the relevant base samples to remove close relatives, researchers should make every effort to choose base and target datasets that are less likely to contain highly related individuals.

*  **Statistical Power Trade-Off**: While removing close relatives can help control the risk of inflation, it may come at a cost to statistical power, particularly if the base and target datasets come from different populations. Analyzing such diverse datasets can be challenging, so there's a balance to strike between controlling relatedness and maintaining statistical power.

*  **Similarity Between Base and Target Samples**: Ideally, the base and target datasets should be as similar as possible in terms of population characteristics, genetics, and environmental factors, without compromising the integrity of the analysis. This ensures that the PRS associations are not confounded by differences between the datasets.

In summary, the passage emphasizes the importance of addressing the issue of close genetic relatedness between individuals in the base and target datasets when conducting PRS analysis. Strategies include removing close relatives if possible and selecting datasets that are similar in population characteristics while carefully considering the trade-off between controlling relatedness and statistical power.

## QC of Target Data

This section will use PLINK for some parts so please install PLINK.
### Genome Build of Target data
Genome build is same between base and target data.

### Standard GWAS QC Target data

In this step we will filter the SNPs based on the QC parameters that are set for target data.  To achieve this we will be using PLINK command line.
#### **Plink flags and filters** 
* `bfile` String that corresponds to prefix of genomtype files.  Here is is  `EUR`.
* `maf`   Minor allele frequency cutoff.  Alleles with frequency less than this are removed.  Set as `0.01`.  This value depends on the sample size.  for small sample size genotyping errors will have influence of SNPs with low maf. For larger sample sizes maf can be set to lower thresholds.
*  `hwe`  The HWE (Hardy-Weinberg Equilibrium) p-value is a statistical test used in genetics to assess whether a population is in equilibrium with respect to the genotype frequencies of a specific genetic variant, typically a single nucleotide polymorphism (SNP). The Hardy-Weinberg Equilibrium is a fundamental principle in population genetics, stating that in the absence of evolutionary forces (such as selection, mutation, migration, or genetic drift), the genotype frequencies of a genetic variant will remain constant from generation to generation. The HWE p-value is calculated by comparing the observed genotype frequencies in a sample population to the expected genotype frequencies under Hardy-Weinberg Equilibrium. The test helps researchers determine whether there are deviations from this equilibrium, which could be indicative of various factors, including selection, genotyping errors, population substructure, or other biological processes. Here the cut off is `1e-6`.  Values larger than this are significantly deviating from HWE so need to be discarded.
* `geno` This parameters excludes SNPs that are missing in a high fraction of samples. Here `0.01`.  Missing in >1% of samples 
* `mind` This will filter the samples where genotypes are missing.  These samples are not informative and need to be dropped. Here `0.01`.
* `make-just-fam` Informs `plink`` to only generate the QC'ed sample name to avoid generating the .bed file.
* `write-snplist` Informs plink to only generate the QC'ed SNP list to avoid generating the .bed file.
* `out` Output file prefixes, here `EUR.QC`

```{shell}



plink \
    --bfile EUR \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC

####PROCESS OUTPUT MESSAGE
16384 MB RAM detected; reserving 8192 MB for main workspace.
551892 variants loaded from .bim file.
503 people (240 males, 263 females) loaded from .fam.
14 people removed due to missing genotype data (--mind).
IDs written to EUR.QC.irem .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 489 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 1806 het. haploid genotypes present (see EUR.QC.hh ); many commands
treat these as missing.
Total genotyping rate in remaining samples is 0.999816.
5353 variants removed due to missing genotype data (--geno).
Warning: --hwe observation counts vary by more than 10%, due to the X
chromosome.  You may want to use a less stringent --hwe p-value threshold for X
chromosome variants.
--hwe: 944 variants removed due to Hardy-Weinberg exact test.
5061 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
540534 variants and 489 people pass filters and QC.
Note: No phenotypes present.
List of variant IDs written to EUR.QC.snplist .
--make-just-fam to EUR.QC.fam ... done.
```
* `14`        people remoced due to missing Genotype.
* `0.999816`  Genotyping rate in remaining samples
* `5353`      Variants remopoved due to misising genotypes
* `944`       variants removed due to Hardy-Weinberg exact test.
* `540534 variants and 489 people pass filters and QC.`

#### **Filtering based on heterozygocity rates**
 Very high heterozygocity is indicator of contamination and very low heterozygocity indicates inbreeding.  There fore it is advisable to remove extreme heterozygocity scores. This will be done in a 2 step process.
##### STEP1: Removing highly correlated SNPs
* `bfile`	**EUR**	Informs plink that the input genotype files should have a prefix of EUR
* `keep`	**EUR.QC.fam**	Informs plink that we only want to use samples in EUR.QC.fam in the analysis
* `extract`	**EUR.QC.snplist**	Informs plink that we only want to use SNPs in EUR.QC.snplist in the analysis
* `indep-pairwise`	**200 50 0.25**	Informs plink that we wish to perform pruning with a window size of 200 variants, sliding across the genome with step size of 50 variants at a time, and filter out any SNPs with LD r² higher than 0.25
* `out`	**EUR.QC**	Informs plink that all output should have a prefix of EUR.QC




```{shell}
plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC


    ## OUTPUT
    16384 MB RAM detected; reserving 8192 MB for main workspace.
551892 variants loaded from .bim file.
503 people (240 males, 263 females) loaded from .fam.
--extract: 540534 variants remaining.
--keep: 489 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 489 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 1273 het. haploid genotypes present (see EUR.QC.hh ); many commands
treat these as missing.
Total genotyping rate in remaining samples is 0.999919.
540534 variants and 489 people pass filters and QC.
Note: No phenotypes present.
Pruned 20214 variants from chromosome 1, leaving 21127.
Pruned 20984 variants from chromosome 2, leaving 20879.
Pruned 17463 variants from chromosome 3, leaving 17812.
Pruned 16344 variants from chromosome 4, leaving 16690.
Pruned 15027 variants from chromosome 5, leaving 15981.
Pruned 20731 variants from chromosome 6, leaving 15475.
Pruned 14046 variants from chromosome 7, leaving 14606.
Pruned 13289 variants from chromosome 8, leaving 13667.
Pruned 10723 variants from chromosome 9, leaving 11951.
Pruned 12722 variants from chromosome 10, leaving 13383.
Pruned 13053 variants from chromosome 11, leaving 12666.
Pruned 11654 variants from chromosome 12, leaving 12500.
Pruned 8429 variants from chromosome 13, leaving 9609.
Pruned 8072 variants from chromosome 14, leaving 8791.
Pruned 7978 variants from chromosome 15, leaving 8470.
Pruned 8808 variants from chromosome 16, leaving 9500.
Pruned 8032 variants from chromosome 17, leaving 8724.
Pruned 7204 variants from chromosome 18, leaving 8560.
Pruned 7014 variants from chromosome 19, leaving 6726.
Pruned 6341 variants from chromosome 20, leaving 7508.
Pruned 3629 variants from chromosome 21, leaving 4280.
Pruned 4085 variants from chromosome 22, leaving 4323.
Pruned 16235 variants from chromosome 23, leaving 5229.
Pruning complete.  272077 of 540534 variants removed.
Marker lists written to EUR.QC.prune.in and EUR.QC.prune.out .
```

This will generate two files 1) EUR.QC.prune.in and 2) EUR.QC.prune.out. All SNPs within EUR.QC.prune.in have a pairwise r²<0.25.

#### STEP2 Heterozygosity rates can then be computed using plink:
This will generate the `EUR.QC.het` file, which contains F coefficient estimates for assessing heterozygosity. We will remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean.

The `F coefficient`, often referred to as the inbreeding coefficient or coefficient of inbreeding, is a measure in genetics that quantifies the probability that two alleles at a specific genetic locus in an individual are identical by descent, meaning they are inherited from a common ancestor. In other words, it measures the level of relatedness between an individual's parents in terms of a specific gene or genetic variant.

The F coefficient estimates typically fall within the range of 0 to 1, with higher values indicating a higher degree of inbreeding or relatedness. Here's how to interpret F coefficient estimates:

* **F = 0**: An F coefficient of 0 indicates that there is no inbreeding, and the alleles at the specified locus are not identical by descent. In other words, the individual's parents are not closely related with respect to that specific genetic variant.

* **F > 0**: A positive F coefficient greater than 0 indicates some level of inbreeding or relatedness. The higher the F value, the more closely related the individual's parents are with respect to the specified locus. This suggests that the alleles at that locus have a higher chance of being inherited from a common ancestor.

* **F = 1**: An F coefficient of 1 indicates complete inbreeding or identical-by-descent alleles at the specified locus. This means that both alleles at the locus are inherited from the same ancestor, such as when an individual has parents who are close relatives or siblings.

The F coefficient is used in various genetic studies, particularly in populations with a history of consanguineous (closely related) mating, to assess the risk of homozygosity for rare recessive genetic variants, which can lead to an increased risk of genetic disorders. It is also used in studies of population genetics to measure the degree of relatedness and inbreeding within populations.

```{shell}
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC

    ### OUTPUT
    16384 MB RAM detected; reserving 8192 MB for main workspace.
551892 variants loaded from .bim file.
503 people (240 males, 263 females) loaded from .fam.
--extract: 268457 variants remaining.
--keep: 489 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 489 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 558 het. haploid genotypes present (see EUR.QC.hh ); many commands
treat these as missing.
Total genotyping rate in remaining samples is 0.999958.
268457 variants and 489 people pass filters and QC.
Note: No phenotypes present.
--het: 263228 variants scanned, report written to EUR.QC.het 

head EUR.QC.het
     FID       IID       O(HOM)       E(HOM)        N(NM)            F
  HG00096   HG00096       215057    2.148e+05       263228     0.005196
  HG00097   HG00097       214695    2.148e+05       263228     -0.00228
  HG00099   HG00099       214846    2.148e+05       263228    0.0008384
  HG00101   HG00101       214569    2.148e+05       263228    -0.004882
  HG00102   HG00102       215187    2.148e+05       263228     0.007881
  HG00103   HG00103       215264    2.148e+05       263228     0.009471
  HG00105   HG00105       214725    2.148e+05       263228     -0.00166
  HG00107   HG00107       214943    2.148e+05       263228     0.002842
  HG00108   HG00108       214855    2.148e+05       263228     0.001024

```
The filteration based on F score can be done in `R` as shown below.

```{R}
library(data.table)

# Read the file as dataframe
df <- fread("EUR.QC.het")

# calulate mean and stdev
m <- mean(df$F)
std <- sd(df$F)

df_filter <- df[df$F < (m+(3*std)) & df$F > (m-(3*std)),]
fwrite(df_filter[,c("FID","IID")], file="EUR.valid.sample", sep="\t") 

dim(df)
# [1] 489   6
dim(df_filter)
# [1] 487   6
```
2 samples were dropped because of heterozygocity issue based on F coefficient estimates.  The samples are

```
       FID     IID O(HOM) E(HOM)  N(NM)       F
1: HG00324 HG00324 217015 214800 263228 0.04563
2: NA20585 NA20585 217236 214800 263228 0.05020

## Use this R code to get the list
df[df$F > (m+(3*std)) | df$F < (m-(3*std)),]
```

### Mismatching SNPs in Target data

In the scenerio where there is allele mismatching between base and target data due to the fact that they were reported from different strands of DNA.  This can be resolved by switching/flipping the strands and reporting the complementary allele of the reported allele in target.  Most PRS application can perform strand flipping in case of allele mismatch, however if a user want to do it indfependently below is the code to show that.  Also we may also need to recode the SNPs in target where the SNPs are not reported in correct format.  This will be in cases where base data have **A1** as effect allele and **A2** as Non-effect allele but the target data has **A1** as Non-effect allele and **A2** as effect allele. This can also be combined with bases in target reported from complementary strands. Basically, recoding + strand_flipping

To illustrate this example is shown below,

```{R}
            Base                  Target
        ______________        _______________
          A1      A2            A1      A2
.....................................................
snp1      G       T             C       A          - Needs Strand flipping Base.A1=comp(Target.A1) and Base.A2=comp(Target.A2) 
snp2      A       C             T       G          - Needs Strand flipping Base.A1=comp(Target.A1) and Base.A2=comp(Target.A2)
snp3      T       C             A       G          - Needs Strand flipping Base.A1=comp(Target.A1) and Base.A2=comp(Target.A2)
snp4      C       A             G       T          - Needs Strand flipping Base.A1=comp(Target.A1) and Base.A2=comp(Target.A2)
.....................................................
snp5      G       T             T       G          - Needs recoding base.A1 = Taget.A2 and Base.A2 = Target.A1
snp6      A       C             C       A          - Needs recoding base.A1 = Taget.A2 and Base.A2 = Target.A1
snp7      T       C             C       T          - Needs recoding base.A1 = Taget.A2 and Base.A2 = Target.A1
.....................................................
snp8      G       T             A       C          - Needs recoding and Strand flipping base.A1 = comp(Taget.A2) and Base.A2 = comp(Target.A1)
snp9      A       C             G       T          - Needs recoding and Strand flipping base.A1 = comp(Taget.A2) and Base.A2 = comp(Target.A1)
snp10     T       C             G       A          - Needs recoding and Strand flipping base.A1 = comp(Taget.A2) and Base.A2 = comp(Target.A1)

comp : Complementary base/allele
Note : Above examples are not an exhaustive list but some variations are shown for understanding.
```
Below is the R code to achieve strand flipping, recoding and recoding combined with strand flipping. In these cases we will only use the SNPs that has passed our QC so far listed in `EUR.QC.snplist`.

**Strand flipping**
```
# Strand flipping code, we need 3 files for this
#file1 : Base summary stat file, columns are describe above. 
#file2 : Target bim file that has SNPs 
#         (columns are chromosome, SNP, DistancebetweenSNPs in Centimorgan Scale, Base position, A1 aand A2 allele)
# snplist : SNPs that pass QC
#
#Lets read all the files in R  
library(data.table)
# list files from data folder
list.files("./data")

#read in file1 : Base summary stat file, columns are describe above. 
base <- fread("./data/height_GWASQC_dedup.txt")

head (base,3)
#CHR     BP        SNP A1 A2      N         SE        P        OR      INFO       MAF
#1:   1 756604  rs3131962  A  G 388028 0.00301666 0.483171 0.9978869 0.8905579 0.3693896
#2:   1 768448 rs12562034  A  G 388028 0.00329472 0.834808 1.0006873 0.8958935 0.3368458
#3:   1 779322  rs4040617  G  A 388028 0.00303344 0.428970 0.9976036 0.8975083 0.3773680

# file2 : Target bim file that has SNPs 
target <- fread("./data/EUR.bim")
head(target,3)
# V1         V2       V3     V4 V5 V6
# 1:  1  rs3131962 0.490722 756604  A  G
# 2:  1 rs12562034 0.495714 768448  0  0
# 3:  1  rs4040617 0.500708 779322  G  A

# Assign proper column names and same as in base file where values are matching
colnames(target) <- c("CHR","SNP","CM","BP","T.A1","T.A2")

head(target,3)
#   CHR        SNP       CM     BP T.A1 T.A2
# 1:   1  rs3131962 0.490722 756604    A    G
# 2:   1 rs12562034 0.495714 768448    0    0
# 3:   1  rs4040617 0.500708 779322    G    A 

# Load the QCed SNP list
snplist <- fread("./data/EUR.QC.snplist", header=FALSE)  

head(snplist,3)
# V1
# 1:  rs3131962
# 2:  rs4040617
# 3:  rs79373928

# In order to be consisten convert alleles to uppercase in both abse and target.
base$A1 <- toupper(base$A1)
base$A2 <- toupper(base$A2)
target$T.A1 <- toupper(target$T.A1)
target$T.A2 <- toupper(target$T.A2)

# also merge base and target objects to get final object that we will process together.
mergd <- merge(target, base, by=c("CHR","SNP", "BP"))

# select only alleles that have passed QC listed in snplist
mergd <- mergd[mergd$SNP %in% snplist$V1,]

head(mergd,3)
#CHR       SNP        BP       CM T.A1 T.A2 A1 A2      N         SE           P        OR      INFO MAF
#1:   1 rs1000283 209894661 232.2630    A    G  A  G 388028 0.00253106 3.06580e-09 0.9851075 0.9101943  0.3150172
#2:   1 rs1000313  15405489  31.3668    G    A  G  A 388028 0.00253577 3.09630e-02 0.9945439 0.9203394  0.3383716
#3:   1 rs1000494  29789167  52.1518    G    T  G  T 388028 0.00205227 2.81824e-02 0.9955059 0.9052882  0.4761347

dim(mergd)
# [1] 489816     14

# Step1 Identify SNPs that are consistent between base and target
mergd.match <- subset(mergd, T.A1==A1 & T.A2==A2)

dim(mergd.match)
# [1] 485884     14

# STEP2 lets identify the SNPs that need strand flipping.  We can identify them as
# base.A1 == complementary(target.A1) and base.A2 == complementary(target.A2) 
complement <- function(x) {
  switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  )
}

# to test we will create new column with complementary bases of T.A1 and T.A2
mergd$C.A1 <- sapply(mergd$T.A1, complement)
mergd$C.A2 <- sapply(mergd$T.A2, complement)

mergd.strand <- subset(mergd,A1 == C.A1 & A2 == C.A2)

# SNPs that need strand switching
complement.snps <- target$SNP %in% mergd.strand$SNP

# Update SNPs in the target object
target[complement.snps,]$T.A1 <-
  sapply(target[complement.snps,]$T.A1, complement)
target[complement.snps,]$T.A2 <-
  sapply(target[complement.snps,]$T.A2, complement)


# STEP3 SNPs need encoding

# identify SNPs that need recoding
mergd.recode <- subset(mergd, A1 == T.A2 & A2 == T.A1)
# Update the recode SNPs
recode.snps <- target$SNP %in% mergd.recode$SNP
tmp <- target[recode.snps,]$T.A1
target[recode.snps,]$T.A1 <- target[recode.snps,]$T.A2
target[recode.snps,]$T.A2 <- tmp

# identify SNPs that need recoding & complement
mergd.crecode <- subset(mergd, A1 == C.A2 & A2 == C.A1)
# Update the recode + strand flip SNPs
comp.snps <- target$SNP %in% mergd.crecode$SNP
tmp <- target[comp.snps,]$T.A1
target[comp.snps,]$T.A1 <- as.character(sapply(target[comp.snps,]$T.A2, complement))
target[comp.snps,]$T.A2 <- as.character(sapply(tmp, complement))

# Output updated target file
write.table(
  target[,c("SNP", "T.A1")],
  "EUR.a1",
  quote = F,
  row.names = F,
  col.names = F,
  sep="\t"
)

# STEP3 Detect SNPs with disparate alleles between the base and target sequences, often arising 
# from disparities in genome builds or insertion-deletion (Indel) events."

mismatch <-
  target$SNP[!(target$SNP %in% mergd.match$SNP |
              target$SNP %in% mergd.strand$SNP | 
              target$SNP %in% mergd.recode$SNP |
              target$SNP %in% mergd.crecode$SNP)]

mergd.mismatch<-mergd[mergd$SNP %in% mismatch,]


write.table(
  mergd.mismatch,
  "EUR.mismatch.txt",
  quote = F,
  row.names = F,
  col.names = F,
  sep="\t"
)
```

### Duplicate SNPs in Target data

Remove duplicate SNPsin target dataset using R or anyother method.  This dataset is simulated so it doesnot have any duplicates.

### Sex Chromosomes in Target data
This step to identify inconsistencies between reported sex and the sex identified based on the genotype.   The Reason behind this match could be many. Here we will use PLINK to test X chromosome homozygosity estimate (F statistic) to determine the sex of individuals in a genetic dataset (target).

* **F Statistic Calculation**: PLINK calculates an F statistic for each individual based on their X chromosome genotypes. This F statistic measures the degree of homozygosity on the X chromosome. A low F statistic suggests that an individual has heterozygous X chromosome genotypes, which is more common in females. A high F statistic suggests homozygosity on the X chromosome, which is more common in males.

* **Thresholds for Sex Determination**:
  * If an individual's F statistic is < 0.2, they are classified as female because this indicates a higher likelihood of X chromosome heterozygosity.
  *   If an individual's F statistic is > 0.8, they are classified as male because this suggests a higher likelihood of X chromosome homozygosity.

* **Intermediate Values**: Individuals with F statistics between 0.2 and 0.8 may be more challenging to classify definitively. In such cases, further investigation or manual inspection may be necessary.

```{Shell}
plink \
     --bfile EUR \
     --extract EUR.QC.prune.in \
     --keep EUR.valid.sample \
     --check-sex \
     --out EUR.QC

# PLINK OUTPUT
16384 MB RAM detected; reserving 8192 MB for main workspace.
551892 variants loaded from .bim file.
503 people (240 males, 263 females) loaded from .fam.
--extract: 268457 variants remaining.
--keep: 487 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 487 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 558 het. haploid genotypes present (see EUR.QC.hh ); many commands
treat these as missing.
Total genotyping rate in remaining samples is 0.999958.
268457 variants and 487 people pass filters and QC.
Note: No phenotypes present.
--check-sex: 5229 Xchr and 0 Ychr variant(s) scanned, 4 problems detected.
Report written to EUR.QC.sexcheck .

```
Lets have a look at the --check-sex file output **EUR.QC.sexcheck** 
```
    FID       IID       PEDSEX       SNPSEX       STATUS            F
  HG00096   HG00096            1            1           OK            1
  HG00097   HG00097            2            2           OK       -0.217
  HG00099   HG00099            2            2           OK      -0.1848
  HG00101   HG00101            1            1           OK            1
  HG00102   HG00102            2            2           OK      -0.3558
  HG00103   HG00103            1            1           OK            1
  HG00105   HG00105            1            1           OK            1
  HG00107   HG00107            1            1           OK            1
  HG00108   HG00108            1            1           OK            1
```
Here we will use the STATUS column to filter the samples that fail the sex check.

```{R}
# Read in EUR.valid.sample and UR.QC.sexcheck

valid <- fread("./data/EUR.valid.sample")
dim(valid)
# [1] 487   2

head(valid,3)
# FID     IID
# 1: HG00096 HG00096
# 2: HG00097 HG00097
# 3: HG00099 HG00099

qc <- fread("./data/EUR.QC.sexcheck")
dim(qc)
# [1] 487   6

head(qc,3)
# FID     IID PEDSEX SNPSEX STATUS       F
# 1: HG00096 HG00096      1      1     OK  1.0000
# 2: HG00097 HG00097      2      2     OK -0.2170
# 3: HG00099 HG00099      2      2     OK -0.1848
table(qc$STATUS)
# OK PROBLEM 
# 483       4 

# So four samples are problemetic, lets remove them from the "EUR.valid.sample" file.
qc.ok <- subset(qc,STATUS=="OK" & FID %in% valid$FID)
dim(qc.ok)
# [1] 483   6

write.table(qc.ok[,c("FID","IID")], "EUR.QC.valid",row.names=F, col.names=F, sep="\t", quote=F)
```

### Sample Overlap in Target data

Sample overlap is a critical consideration when developing and using polygenic risk scores (PRS) in genetic research and clinical applications. Sample overlap refers to the situation where individuals or samples are included in both the discovery (training) dataset and the target (validation or application) dataset used for PRS analysis. The importance of sample overlap in PRS studies can be understood in several key ways:

1. **Generalization and Validity**:
   - Sample overlap can affect the generalizability and validity of PRS models. If the same individuals are used in both the discovery and target datasets, the PRS may perform well on the target dataset but fail to generalize to new, unrelated populations. This can limit the broader utility of the PRS.

2. **Inflation of Performance Metrics**:
   - Sample overlap can lead to inflated performance metrics, such as area under the receiver operating characteristic curve (AUC-ROC) or R-squared values. Since the model has already seen some of the same individuals in the training dataset, it may appear to perform better than it would in a completely independent validation scenario.

3. **Biased Estimation**:
   - Sample overlap can introduce bias into the estimation of effect sizes and risk associations. PRS models derived from overlapping samples may capture genetic effects that are specific to the training population but do not generalize well to other populations. This can result in overestimation or underestimation of risk.

4. **Ethnic and Ancestral Differences**:
   - Sample overlap can be particularly problematic when studying populations with diverse ethnic or ancestral backgrounds. If the discovery dataset primarily consists of one ethnic group and the target dataset includes a different ethnic group, the PRS may perform poorly due to differences in genetic architecture.

5. **Ethical and Privacy Concerns**:
   - Using the same individuals across multiple datasets raises ethical and privacy concerns, especially if the individuals' identities or genetic information are not properly protected. Researchers must ensure that data sharing and handling comply with ethical and legal standards.

To address these issues, it is often recommended to use independent discovery and validation datasets with minimal or no sample overlap when developing and evaluating PRS models. This helps ensure that the PRS is more robust, generalizable, and unbiased in its predictive performance. Additionally, careful consideration of population stratification and the inclusion of diverse populations in the discovery phase can help mitigate some of the challenges associated with sample overlap and improve the applicability of PRS across different groups.

Since the dataset is simulated there is no sample overlap.

### Relatedness of individuals Target data
Closely related individuals in the target data can lead to overfitting when calculating polygenic risk scores (PRS). To address this issue, we're considering pruning closely related individuals by removing those who have a first or second-degree relative in the sample (relatedness coefficient π > 0.125). Here's how you can perform this pruning using PLINK, a commonly used tool in genetic data analysis

```
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.QC

# OUTPUT
16384 MB RAM detected; reserving 8192 MB for main workspace.
551892 variants loaded from .bim file.
503 people (240 males, 263 females) loaded from .fam.
--extract: 268457 variants remaining.
--keep: 483 people remaining.
Using up to 8 threads (change this with --threads).
Before main variant filters, 483 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.999957.
268457 variants and 483 people pass filters and QC (before --rel-cutoff).
Note: No phenotypes present.
Excluding 5229 variants on non-autosomes from relationship matrix calc.
Relationship matrix calculation complete.
0 people excluded by --rel-cutoff.
Remaining sample IDs written to EUR.QC.rel.id .
```

### Generate final QC target file

```
plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist \
    --exclude EUR.mismatch \
    --a1-allele EUR.a1
```

**Parameter**	*Value*	Description
* **bfile**	*EUR*  	Informs plink that the input genotype files should have a prefix of EUR
* **keep** *EUR.QC.rel.id*	Informs plink that we only want to keep samples in EUR.QC.rel.id
* **extract** 	*EUR.QC.snplist*	Informs plink that we only want to use SNPs in EUR.QC.snplist in the analysis
* **exclude**	*EUR.mismatch*	Informs plink that we wish to remove any SNPs in EUR.mismatch
* **a1-allele** *EUR.a1*	Fix all A1 alleles to those specified in EUR.a1
* **out** *EUR.QC*	Informs plink that all output should have a prefix of EUR.QC


## Calculating and analysis of PRS

### PLINK
Calculating a polygenic risk score (PRS) using `plink` involves several steps, including data preprocessing, PRS calculation, and possibly further analysis or interpretation. Here's a basic overview of the steps involved in calculating a PRS using `plink`:

#### Required Dataset

- **File Name** *Description*
- **Height.QC.gz**	*The post-QCed summary statistic*
- **EUR.QC.bed**	*The genotype file after performing some basic filtering*
- **EUR.QC.bim**	*This file contains the SNPs that passed the basic filtering*
- **EUR.QC.fam**	*This file contains the samples that passed the basic filtering*
- **EUR.height**	*This file contains the phenotype of the samples*
- **EUR.cov**	*This file contains the covariates of the samples*

#### Effect Size Standardisation
In polygenic risk score (PRS) calculation, effect size transformation methods are used to standardize the effect sizes (usually beta coefficients) obtained from a genome-wide association study (GWAS) before incorporating them into the PRS calculation. Standardization helps ensure that effect sizes from different SNPs are on a common scale and can be combined effectively. Here are some commonly used effect size transformation methods:

1. **Z-Score Transformation**:
   - The most common approach is to convert beta coefficients (effect sizes) to Z-scores using the standard error (SE) of the beta estimate. The formula is: Z-score = beta / SE(beta).
   - Z-scores are commonly used because they have a mean of 0 and a standard deviation of 1, making them directly interpretable and comparable across different SNPs.

2. **Odds Ratio (OR) Transformation**:
   - If your GWAS reports odds ratios (ORs) instead of beta coefficients, you can transform ORs into beta coefficients by taking the natural logarithm (log(OR)) and as an option then use Z-score transformation.

3. **Rank-Based Transformation**:
   - In some cases, researchers may prefer to use rank-based methods, such as converting the rank of the effect size within the GWAS to a Z-score. This approach is less sensitive to extreme outliers.

4. **Beta Rescaling**:
   - Some researchers rescale beta coefficients by multiplying them by a constant, typically derived from the overall variance explained by the SNPs in the GWAS.

5. **Binary Trait Transformation**:
   - If you are working with binary traits (e.g., disease status), you may need to transform effect sizes to ensure they are appropriate for use in PRS. This transformation can involve logistic regression to estimate the effect size of the SNP on the binary outcome.

6. **Clipping or Winsorization**:
   - In some cases, researchers may choose to clip or Winsorize extreme effect size estimates to mitigate the influence of outliers.

7. **Other Transformations**:
   - Depending on the specific requirements of your analysis or the nature of your data, you may explore other transformation methods to suit your needs.

The choice of effect size transformation method can depend on the characteristics of your GWAS data, the assumptions of your PRS model, and the desired interpretability of the PRS. Z-score transformation is a widely used and robust approach, but the choice of method should be carefully considered in the context of your research.

It's important to document and report the effect size transformation method used in your PRS analysis to ensure transparency and reproducibility of your results. Additionally, the choice of transformation may also impact the performance and interpretation of the PRS, so sensitivity analyses and validation are crucial steps in PRS development.

Here we will do `OR transformation` using R.

```{R}
# Update effect size
library(data.table)

df<- read.table("./data/height_GWASQC_dedup.txt",header = TRUE)
head(df,3)
# CHR     BP        SNP A1 A2      N         SE        P        OR      INFO       MAF
# 1   1 756604  rs3131962  A  G 388028 0.00301666 0.483171 0.9978869 0.8905579 0.3693896
# 2   1 768448 rs12562034  A  G 388028 0.00329472 0.834808 1.0006873 0.8958935 0.3368458
# 3   1 779322  rs4040617  G  A 388028 0.00303344 0.428970 0.9976036 0.8975083 0.3773680

df$BETA <- log(df$OR)
head(df,3)
# CHR     BP        SNP A1 A2      N         SE        P        OR      INFO       MAF        BETA
# 1   1 756604  rs3131962  A  G 388028 0.00301666 0.483171 0.9978869 0.8905579 0.3693896 -0.00211532
# 2   1 768448 rs12562034  A  G 388028 0.00329472 0.834808 1.0006873 0.8958935 0.3368458  0.00068708
# 3   1 779322  rs4040617  G  A 388028 0.00303344 0.428970 0.9976036 0.8975083 0.3773680 -0.00239932

write.table(df,"./data/height_GWASQC_dedup_transformed.txt")
```

#### Clumping
In the context of polygenic risk score (PRS) analysis, "clumping" refers to the process of selecting a subset of single nucleotide polymorphisms (SNPs) based on their linkage disequilibrium (LD) patterns. This process is often used to reduce the number of SNPs used in the PRS calculation while retaining the most informative ones. Clumping helps streamline the analysis by including only a representative set of SNPs that are not in strong LD with each other.

Here's how clumping in LD works in the context of PRS analysis:

1. **Input Data**: You typically start with a set of genome-wide SNPs and their associated summary statistics, which include effect sizes (usually beta coefficients) and p-values from a genome-wide association study (GWAS).

2. **Selecting a Reference SNP**: In PRS analysis, you often have a "target" SNP or a set of SNPs of primary interest (e.g., SNPs associated with a specific trait or disease). One of these SNPs is chosen as the "reference SNP" for clumping.

3. **Defining LD Thresholds**: You set LD thresholds, typically measured using metrics like r-squared (r²) or D' (a measure of the strength of LD between two SNPs). These thresholds determine how strongly SNPs must be correlated to be considered for inclusion in the PRS.

4. **Clumping Process**: Starting with the reference SNP, you iteratively select additional SNPs from your dataset that meet the LD threshold criteria. These selected SNPs form a "clumped" set, which represents the LD structure around the reference SNP.

5. **PRS Calculation**: You then use the clumped set of SNPs to calculate the polygenic risk score for individuals in your target dataset.

The benefits of clumping in PRS analysis include:

- **Dimension Reduction**: Clumping reduces the number of SNPs included in the PRS calculation, making the analysis more computationally efficient and interpretable.

- **Tagging SNPs**: The clumped SNPs are often referred to as "tagging SNPs" because they effectively capture the genetic variation in the region. By studying these tag SNPs, researchers can infer information about the entire LD block.

- **Avoiding Redundancy**: Clumping helps avoid redundancy by including only a subset of SNPs that are representative of LD blocks and not in strong LD with each other.

- **Interpretability**: A smaller set of SNPs is easier to visualize and interpret, making it more feasible to identify potential causal variants or associations with traits.

Clumping in LD is a common step in PRS analysis, especially when working with large-scale genomic data. Different software packages, including PLINK, provide tools to perform LD-based clumping based on user-defined LD thresholds and reference SNPs. The choice of LD thresholds and reference SNPs can impact the performance and interpretability of the resulting PRS, so careful consideration of these factors is essential.

**NOTE : If your target data are small (e.g. N < 500) then you can use the 1000 Genomes Project samples for the LD calculation. Make sure to use the population that most closely reflects represents the base sample.**

Clumping can be performed using the following command in plink.

```{shell}
plink \
    --bfile EUR.QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump height_GWASQC_dedup_transformed.txt \
    --clump-snp-field SNP \
    --clump-field P \
    --out EUR

## OUTPUT
16384 MB RAM detected; reserving 8192 MB for main workspace.
489805 variants loaded from .bim file.
483 people (232 males, 251 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 483 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is exactly 1.
489805 variants and 483 people pass filters and QC.
Note: No phenotypes present.
Warning: 'rs1076829' is missing from the main dataset, and is a top variant.
Warning: 'rs3129818' is missing from the main dataset, and is a top variant.
Warning: 'rs3118359' is missing from the main dataset, and is a top variant.
9809 more top variant IDs missing; see log file.
--clump: 193758 clumps formed from 489805 top variants.
Results written to EUR.clumped .

```

- `--bfile EUR.QC`: Specifies the input dataset in PLINK binary format (BED/BIM/FAM files) with the prefix "EUR.QC." This dataset contains genetic information for a population (presumably of European ancestry) that has already undergone quality control.

- `--clump-p1 1`: Sets the significance threshold for clumping. SNPs with p-values less than or equal to 1 will be considered for clumping. All SNPs are used.

- `--clump-r2 0.1`: Specifies the LD (linkage disequilibrium) threshold. SNPs with an LD correlation (r²) greater than or equal to 0.1 will be considered for clumping.

- `--clump-kb 250`: Defines the physical distance threshold for clumping. SNPs within 250 kilobases (kb) of each other will be considered for clumping together.

- `--clump height_GWASQC_dedup_transformed.txt`: Specifies the name of the GWAS summary statistics file to be used for clumping. In this case, the file is named "height_GWASQC_dedup_transformed.txt."

- `--clump-snp-field SNP`: Specifies the field in the GWAS summary statistics file that contains SNP (Single Nucleotide Polymorphism) IDs.

- `--clump-field P`: Specifies the field in the GWAS summary statistics file that contains p-values associated with each SNP.

- `--out EUR`: Specifies the prefix for the output files generated by the clumping analysis. In this case, the prefix "EUR" will be used for output file names.

This will generate EUR.clumped, containing the index SNPs after clumping is performed. 
```
head -n 3 EUR.clumped 
CHR    F         SNP         BP        P    TOTAL   NSIG    S05    S01   S001  S0001    SP2
   6    1   rs3134762   31210866  4.52e-165      180     20      0      0      0    160 rs2523898(1),rs13210132(1),rs28732080(1),...(trimmed to display here)
   6    1    rs417162   29916505  1.92e-162      229     21      1      1      4    202 rs3129190(1),rs35471702(1),rs3094727(1),...(trimmed to display here)
```


We can extract the index SNP ID by performing the following command:
```{shell}
awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp
```
$3 because the third column contains the SNP ID

#### Generating PRS with PLINK
PLINK provides two options for calculating polygenic risk scores (PRS): --score and --q-score-range. These options allow you to calculate PRS based on external summary statistics or reference panels. Here's an explanation of files needed for it.

1. The base data file : **height_GWASQC_dedup_transformed.txt**
2. A file containing different SNP ID and their corresponding p-values.  This is 3rd and 8th column in **height_GWASQC_dedup_transformed.txt**.  To extract this information run the following `awk` code

```{shell}
awk '{print $3,$8}' height_GWASQC_dedup_transformed.txt > SNP.pvalue
```
3. The `--q-score-range` option is used to specify a range of polygenic risk score (PRS) thresholds for analysis. It allows you to divide your sample into different risk groups based on PRS percentiles or score ranges. A file can be created with ranges of choice either with help of code ot using any text editor.  The "range file" should have range labels in the first column, p-value lower bounds in the second column, and upper bounds in the third column (inclusive). Here calculate PRS corresponding to a few thresholds for illustration purposes.

```
echo "0.001 0 0.001" > range_list 
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list
```

PRS is calculated with following score

```{shell}
plink \
    --bfile EUR.QC \
    --score height_GWASQC_dedup_transformed.txt 3 4 12 header \
    --q-score-range range_list SNP.pvalue \
    --extract EUR.valid.snp \
    --out EUR
#OUTPUT
16384 MB RAM detected; reserving 8192 MB for main workspace.
489805 variants loaded from .bim file.
483 people (232 males, 251 females) loaded from .fam.
--extract: 193758 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 483 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is exactly 1.
193758 variants and 483 people pass filters and QC.
Note: No phenotypes present.
Warning: 305859 lines skipped in --score file (305859 due to variant ID
mismatch, 0 due to allele code mismatch); see EUR.nopred for details.
--score: 193758 valid predictors loaded.
Warning: 305860 lines skipped in --q-score-range data file.
--score: 7 ranges processed.
```
for more detail on option follow the [link](https://www.cog-genomics.org/plink/2.0/score).

The above command and range_list will generate 7 files:

- EUR.0.5.profile
- EUR.0.4.profile
- EUR.0.3.profile
- EUR.0.2.profile
- EUR.0.1.profile
- EUR.0.05.profile
- EUR.0.001.profile

First 3 entries in `EUR.0.001.profile` file.

```{shell}
head -n 3 EUR.0.001.profile

      FID       IID  PHENO    CNT   CNT2    SCORE
  HG00096   HG00096     -9  32678   6868 -5.00082e-05
  HG00097   HG00097     -9  32678   6937 5.42519e-05
  HG00099   HG00099     -9  32678   6753 -3.72876e-06
```
- **CNT** is the number of alleles checked
- **CNT2** is the total number of named alleles actually observed. 

#### Accounting for Population Stratification
Accounting for population stratification is a crucial step in polygenic risk score (PRS) analysis to ensure that the genetic risk scores are not confounded by differences in ancestry or population substructure within your dataset. Failing to address population stratification can lead to spurious associations and incorrect interpretations of PRS results. Here are some strategies to account for population stratification in PRS analysis:

1. **Principal Component Analysis (PCA)**:
   - Conduct PCA on the genotype data to identify and correct for population structure. This involves calculating principal components (PCs) that capture the major axes of genetic variation in the dataset.
   - Include the top PCs as covariates in your regression model when assessing the association between PRS and the trait of interest. This helps control for population stratification effects.

2. **Ancestry Informative Markers (AIMs)**:
   - Use ancestry informative markers, which are SNPs known to have large allele frequency differences between populations, to estimate individual ancestry proportions.
   - Include ancestry proportions as covariates in your analysis to control for population structure.

3. **Sample Matching or Stratification**:
   - If possible, match cases and controls or individuals with different phenotypic groups based on ancestry or population substructure.
   - Alternatively, stratify your analysis by subpopulations if there are distinct subgroups within your dataset.

4. **Genomic Control**:
   - Apply genomic control correction to adjust the test statistics for association between PRS and the trait. This method helps to control for inflation in test statistics caused by population stratification.

5. **Mixed Linear Models (MLM)**:
   - Consider using mixed linear models, which explicitly model the genetic relatedness or kinship between individuals to account for population structure.
   - MLMs can be especially useful for large and diverse datasets.

6. **Validation in External Populations**:
   - Validate PRS associations in external populations that are ethnically or ancestrally diverse to assess the generalizability of your findings.
   - Ensure that your PRS associations hold across different populations.

7. **Stratified Analysis**:
   - Conduct stratified analyses if population stratification is strong or if you have a diverse dataset. Analyze subpopulations separately to avoid confounding.

8. **Use Ethnicity or Ancestry Information**:
   - Collect self-reported ethnicity or ancestry information from participants and include this as a covariate in your analysis.

9. **Genomic Control Inflation Factor (λ)**:
   - Calculate the genomic control inflation factor (λ) to assess the extent of population stratification in your analysis. High values of λ suggest potential stratification issues.

Addressing population stratification in PRS analysis is essential for obtaining reliable and interpretable results. The specific approach you choose may depend on the nature of your data, the available covariate information, and the research question you are addressing. Consulting with a statistician or genetic epidemiologist experienced in population genetics can be valuable in designing and conducting robust PRS analyses while accounting for population stratification.

We will use PCA to account for population stratification.  We will use `plink` to calculate PCA

```{shell}
# First, we need to perform prunning
plink \
    --bfile EUR.QC \
    --indep-pairwise 200 50 0.25 \
    --out EUR
# Then we calculate the first 6 PCs
plink \
    --bfile EUR.QC \
    --extract EUR.prune.in \
    --pca 6 \
    --out EUR

```
 - **--indep-pairwise 200 50 0.25** To generate a pruned subset of SNPs that are in approximate linkage equilibrium with each other. This can be achieved via two commands: --indep which prunes based on the variance inflation factor (VIF), which recursively removes SNPs within a sliding window; second, **--indep-pairwise** which is similar, except it is based only on pairwise genotypic correlation. The first two parameters (200 and 50) are the (window size and step); the third parameter represents the r^2 threshold. Note: this represents the pairwise SNP-SNP metric based on the genotypic correlation, i.e. it does not involve phasing. It remove one of a pair of SNPs if the LD is greater than 0.25,

--------------
- #### **How to decide on number of principal components?**
- Principal Component Analysis (PCA), LDSC (Linkage Disequilibrium Score Regression), and Polygenic Risk Score (PRS) analysis are three distinct but related techniques often used in genetic epidemiology and genomics research. Here's an overview of each and how they can be used together:

- 1. **Principal Component Analysis (PCA)**:
   - PCA is a dimensionality reduction technique used to analyze genetic variation within a dataset. It identifies the major axes of genetic variation (principal components) by linearly combining genetic markers (e.g., SNPs) in a way that maximizes variance.
   - PCA can help address population stratification, which is the presence of subpopulations within a dataset that can confound genetic association studies. By including the top principal components as covariates in a regression analysis, you can control for population structure.
   - PCA can also be used to visualize genetic relatedness among individuals, which can be helpful for quality control and identifying outliers.

- 2. **Linkage Disequilibrium Score Regression (LDSC)**:
   - LDSC is a statistical method used to estimate heritability and partition the genetic correlation between traits or diseases. It leverages information on linkage disequilibrium (LD) between genetic variants to estimate the proportion of phenotypic variance explained by common genetic variants.
   - LDSC can be used to estimate the heritability of a trait or disease, which quantifies the proportion of phenotypic variance that is attributable to genetic factors.
   - LDSC can also estimate genetic correlations between two traits, providing insights into their shared genetic basis. This can be useful for understanding the genetic relationship between different traits or diseases.

- 3. **Polygenic Risk Score (PRS) Analysis**:
   - PRS analysis involves calculating a polygenic risk score for each individual in a dataset based on the cumulative effect of multiple genetic variants (usually SNPs) associated with a trait or disease.
   - PRS is used to predict an individual's genetic predisposition to a trait or disease. It can be applied to various traits, including complex diseases, height, cognitive abilities, and more.
   - PRS can be adjusted for population stratification using PCA components or other methods to ensure that the calculated scores are not confounded by genetic ancestry.

How They Can Be Used Together:
- PCA can be performed initially to identify and correct for population stratification by including the top principal components as covariates in LDSC or PRS analyses.
- LDSC can help estimate heritability and assess the genetic correlation between traits, which can inform the selection of SNPs for PRS calculation and prioritize traits for further analysis.
- PRS analysis can utilize LDSC results to identify relevant SNPs for inclusion in the PRS calculation, and it can also incorporate PCA components to account for population structure.

The combination of these techniques allows researchers to better understand the genetic architecture of traits or diseases, control for population stratification, estimate heritability, and predict genetic risk based on polygenic scores. Proper data preprocessing, quality control, and careful interpretation of results are essential when using these methods together in genetic studies.

Note : If the base and target samples are collected from different worldwide populations then the results from the PRS analysis may be biased

----------

### Finding best-fit PRS
This is from the original tutorial as it is pretty concise.

The P-value threshold that provides the "best-fit" PRS under the C+T method is usually unknown. To approximate the "best-fit" PRS, we can perform a regression between PRS calculated at a range of P-value thresholds and then select the PRS that explains the highest phenotypic variance (please see Section 4.6 of our paper on overfitting issues). This can be achieved using R as follows:
```{R}
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
# Read in the phenotype file 
phenotype <- read.table("EUR.height", header=T)
# Read in the PCs
pcs <- read.table("EUR.eigenvec", header=F)
# The default output from plink does not include a header
# To make things simple, we will add the appropriate headers
# (1:6 because there are 6 PCs)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
# Read in the covariates (here, it is sex)
covariate <- read.table("EUR.cov", header=T)
# Now merge the files
pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
# We can then calculate the null model (model with PRS) using a linear regression 
# (as height is quantitative)
null.model <- lm(Height~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
# And the R2 of the null model is 
null.r2 <- summary(null.model)$r.squared
prs.result <- NULL
for(i in p.threshold){
    # Go through each p-value threshold
    prs <- read.table(paste0("EUR.",i,".profile"), header=T)
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # Now perform a linear regression on Height with PRS and the covariates
    # ignoring the FID and IID from our model
    model <- lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    # model R2 is obtained as 
    model.r2 <- summary(model)$r.squared
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    prs.coef <- summary(model)$coeff["SCORE",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}
# Best result is:
prs.result[which.max(prs.result$R2),]
q() # exit R

```
The best fit PRS is 0.3 and it explains 0.1612372 of phenotypic variation.
