Project Proposal: Methylation Badassays
================
February 14, 2017

### Motivation and background work

DNA methylation (DNAm) - the covalent modification of DNA at locations where cytosine is followed by guanine (CpG sites), resulting in attached methyl groups - is the best understood component of epigenetic machinery \[1\]. After the discovery of DNAm, researchers have linked DNAm to many diseases like cancer and autism \[2-3\]. These marks vary across cell types, temporal development, and can change due to environmental stimuli. These factors are all taken into account when researchers are trying to use gene expression analysis when studying diseases. In recent years, several DNAm studies have suggested that a large portion of DNAm variability is associated with genetic ancestry and is heritable \[4\]. The difference in methylation at specific CpG sites between populations with different ethnic backgrounds are likely due to the presence of between-population differences in single nucleotide polymorphism (SNP) allele frequencies \[7-8\] and allele-specific DNA methylation \[9-10\]. Differentially methylated CpG sites associated with pathology can be confounded by CpGs associated with genetic ancestry causing spurious results \[7-8\]. Therefore, genetic ancestry, as a covariate, needs to be accounted for in any epigenome-wide association study (EWAS).

Although there are some studies investigating the population-specificity of human DNAm \[4, 11-12\], DNAm profiles across tissue types is extremely variable \[13\], and the amount of variability that can be accounted for by ethnicity in placenta samples have not yet been examined. Therefore, in order to investigate how epigenetics affects perinatal health, it is important for us to identify genetic ancestry associated CpGs to figure out true positives.

**Hypothesis:** DNAm variability in the placenta due to ethnicity needs to be accounted for in large scale DNAm studies, or else interpretation of results will be difficult because of confounding effects.

**Aim 1:** In this project, we propose to identify a set of top-ranked CpG sites in placental samples that are associated with self-reported ethnicity (Caucasians and Asians).

**Aim 2:** Using ancestry-associated CpG sites from aim 1, we will determine the genetic ancestry of a second dataset where genetic ancestry was not taken into account.

\*\*Note: It is important to distinguish between genetic ancestry, race and ethnicity. The latter two are social constructs and have no genetic definition. In contrast, genetic ancestry is a continuum which describes the architecture of genome variation between populations \[14\].

1.  Law JA, Jacobsen SE: Establishing, maintaining and modifying DNA methylation patterns in plants and animals. Nat Rev Genet 2010, 11:204-220.

2.  Cicek MS, Koestler DC, Fridley BL, Kalli KR, Armasu SM, Larson MC, Wang C, Winham SJ, Vierkant RA, Rider DN. Epigenome-wide ovarian cancer analysis identifies a methylation profile differentiating clear-cell histology with epigenetic silencing of the HERG K+ channel. Human molecular genetics. 2013; 22(15):3038–47.

3.  Wong CC, Meaburn EL, Ronald A, Price TS, Jeffries AR, Schalkwyk LC, Plomin R, Mill J. Methylomic analysis of monozygotic twins discordant for autism spectrum disorder and related behavioural traits. Molecular psychiatry. 2013

4.  Fraser HB, Lam LL, Neumann SM, Kobor MS. Population-specificity of human DNA methylation. Genome Biol. 2012;13(2):R8.

5.  Heyn H, Moran S, Hernando-Herraez I, Sayols S, Gomez A, Sandoval J, Monk D, Hata K, Marques- Bonet T, Wang L. DNA methylation contributes to natural human variation. Genome research. 2013; 23:1363–1372.

6.  Kwabi-Addo B, Wang S, Chung W, Jelinek J, Patierno SR, Wang BD, Andrawis R, Lee NH, Apprey V, Issa JP. Identification of differentially methylated genes in normal prostate tissues from African American and Caucasian men. Clin Cancer Res. 2010; 16(14):3539?7.

7.  Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D. Principal components analysis corrects for stratification in genome-wide association studies. Nat Genet. 2006; 38(8): 904-9.

8.  Cavalli-Sforza LL, Edwards AW. Phylogenetic analysis. Models and estimation procedures. Am J Hum Genet. 1967; 19(3 Pt 1):233–57.

9.  Boks MP, Derks EM, Weisenberger DJ, Strengman E, Janson E, Sommer IE, Kahn RS, Ophoff RA. The relationship of DNA methylation with age, gender and genotype in twins and healthy controls. PLoS One. 2009; 4(8):e6767.

10. Bell JT, Pai AA, Pickrell JK, Gaffney DJ, Pique-Regi R, Degner JF, Gilad Y, Pritchard JK. DNA methylation patterns associate with genetic and gene expression variation in HapMap cell lines. Genome Biol. 2011; 12(1):R10.

11. Barfield RT, Almli LM, Kilaru V, et al. Accounting for population stratification in DNA methylation studies. Genet Epidemiol. 2014;38(3):231-241.

12. Phillips C, Aradas AF, Kriegel A, et al. Eurasiaplex: A forensic SNP assay for differentiating european and south asian ancestries. Forensic Science International: Genetics. 2013;7(3):359-366.

13. Illingworth R, Kerr A, Desousa D, Jorgensen H, Ellis P, et al. (2008) A novel CpG island set identifies tissue-specific methylation at developmental gene loci. PLoS Biol 6: e22.

14. Yudell M, Roberts D, DeSalle R, Tishkoff S. SCIENCE AND SOCIETY. taking race out of human genetics. Science. 2016;351(6273):564-565

### Division of labour

<table style="width:97%;">
<colgroup>
<col width="19%" />
<col width="19%" />
<col width="19%" />
<col width="19%" />
<col width="19%" />
</colgroup>
<thead>
<tr class="header">
<th><strong>Name</strong></th>
<th><strong>Department/Program</strong></th>
<th><strong>Experties/Interests</strong></th>
<th><strong>GitHub ID</strong></th>
<th><strong>Tasks</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Victor Yuan</td>
<td>Genome Science and Technology</td>
<td>Placental Epigenetics</td>
<td><span class="citation">@wvictor14</span></td>
<td>Proposal, first data peak, preprocessing, data normalization, Validation, Poster</td>
</tr>
<tr class="even">
<td>Michael Yuen</td>
<td>Medical Genetics</td>
<td>Cancer Genomics</td>
<td><span class="citation">@myuen89</span></td>
<td>preprocessing, data normalization, Progress report, Poster</td>
</tr>
<tr class="odd">
<td>Nivretta Thatra</td>
<td>Bioinformatics</td>
<td>Neuroscience</td>
<td>@nivretta</td>
<td>Proposal, first data peak, Statistical analysis (ID differentially methylated CpGs), Progress report, Poster</td>
</tr>
<tr class="even">
<td>Ming Wan</td>
<td>Statistics</td>
<td>Statistical Genetics and Machine Learning</td>
<td><span class="citation">@MingWan10</span></td>
<td>preprocessing, data normalization, Statistical analysis (ID differentially methylated CpGs), Poster</td>
</tr>
<tr class="odd">
<td>Anni Zhang</td>
<td>Genome Science and Technology</td>
<td>Pancreatic Cancer/Cancer metabolism</td>
<td><span class="citation">@annizubc</span></td>
<td>Proposal, Validation, Progress report, GitHub repository maintenance, Poster</td>
</tr>
</tbody>
</table>

### Dataset

**Aim 1:** We are working with human placental tissue from 45 subjects; 20 controls and 25 patients. All subjects’ metadata is contained in a txt file with columns for the sample name, group (control vs. ), sex, ethnicity, . DNA methylation was measured at 450,000 CpG sites in each of these samples using the 450K microarray from Illumina.

[Link to our data subdirectory](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/data)

**Aim 2:** We will have access to placental datasets in which the genetic ancestry is unknown. Victor, can you expand on these data? Add a README.md file in the subdirectory to describe your data in more detail, e.g., how many rows, how many columns, experimental design etc.

### Aims and methodology

**Aim 1:** Identify CpG sites which are specifically associated with Caucasians or Asians.

**Methods:** We will preprocess the data using minfi and/or preprocessNoob, both of which are bioconductor packages to import and analyze Illumina methylation data \[15\]. Quantile normalization will then be used \[16\]. We will use ANOVA or logistic regression (rather than linear regression) to find the effect of genetic ancestry on methylation.

References:

1.  Bioconductor tutorial: [Analysis of 450k data using minfi](https://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#introduction)
2.  STAT540 seminar on [quantile normalization](https://stat540-ubc.github.io/seminars/sm08_methylation.html)

**Aim 2:** Use the identified CpG sites from aim 1 to determine the genetic ancestry of a second dataset whose genetic ancestry is unknown.

**Methods:** To visually assess the data, we will generate a heatmap of correlations between subjects from the second dataset to the subjects from the first dataset. To cluster the subjects from the second dataset we will explore hierarchical clustering (dendrograms) and principal component analysis (combine all samples of second dataset and see if they cluster with the Caucasian or Asian subset of the first dataset). We also will do cross validation - using random sampling - to ensure that our clusters are real. This may be a challenge due to our small sample sizes.

Please provide a link to your project\_proposal.md file in your group repo README.md. This is a group-level deliverable so you’d submit it by opening an issue in your group repository, with the issue named “Project proposal of \[group name\]”, and provide the revision SHA and tag the TAs and instructors who are mentors of your group by Feb. 15th, 2017.