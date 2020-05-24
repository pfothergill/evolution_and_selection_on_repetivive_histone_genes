# Asexual vs Sexual Reproductive Evolution and Selection Occuring on Repetivive Histone Genes
It has long been known that modes of reproduction, e.g. sexual vs asexual, can have profound influences on the efficacy of natural selection and thus can influence the rate of evolution and the fate of a population. Sexual reproduction is not only important for genetic diversity, but it is key in the process of natural selection by allowing beneficial mutations to arise in populations and disposing of deleterious ones. But what happens in asexual populations that consist of self-cloning females? A higher accumulation of nonsynonymous changes is often seen within asexual species. However, our knowledge on the effects of asexual reproduction on different genomic characteristics is far from complete. A model organism for studying sexual vs asexual reproduction is Potamopyrgus antipodarum, a New Zealand freshwater snail with both sexual and asexual reproducing natural populations. Whole genome duplications occurring in the evolutionary history of P. antipodarum have created different ploidy snails leading to different modes of reproduction. This fresh water snail species can be sexual (diploid) or asexual (triploid and tetraploids). Genomic sequencing of P. antipodarum has revealed a higher number of repetitive histone genes in asexual lineages compared to sexual lineages. This raises the following questions: Why do asexual lineages need an exaggerate amount of histone genes? Are they all functional? Is natural selection acting in the same way on these asexual histone genes as it does on sexual ones? 
This research focuses on understanding the evolution and selective pressures occurring on repetitive histone genes in asexual vs sexual lineages. This project focuses on the four sequences of histone genes: H2A, H2B, H3, and H4. A statistical approach is used to calculate the rate of substitutions at synonymous sites compared to the rate of substitutions at nonsynonymous sites between repetitive histone genes within a single individual that represents a lineage. The lack of knowledge in functionality and location of these genes, raises the question as to if these genes are under any selective pressure and if that differs in snails of sexual vs asexual reproduction.

## Data
### master/data/deviate_histones_raw
The data has been generated with DeviaTE software. There are 26 lineages of *Potamopyrgus antipodarum*, and 1 out group, and include the sequences of H2A, H2B, H3, and H4 genes
## Methods
### master/script/
A python script was written to analyze the rate of nonysonmous and sysnonymous changes occuring in repetive histone genes and compare sexual and asexual snail lineages.
### master/analysis/
Rstudio was used to create boxplots and run non parameteric statistical tests.
## Results
### master/output/
The python script procuded NS_results.tx. 

R analysis results indclude 2 tables of compiled non parameteric statical tests and 4 figures. 

## Conclusion 
At the start of this research, the following questions were raised:
-	Why do asexual lineages have more copies of repetitive histone genes? 
-	Are they all being used? If not, would we see a higher accumulation of nonsynonymous changes? 
-	If we see about the same amount of nonsynonymous to synonymous changes occurring between sexuals and asexuals, can we assume similar modes of selection? If not, can we assume different modes of selection are acting? 
-	How do we look at selection among repetitive genes? 
 
These questions lead me to create N/S values using a Python script. While N/S values do give us insight into the nonsynonymous and synonymous polymorphisms that are occurring among the repetitive histone genes, N/S does not directly tell us what selective pressures are occurring. I hoped to identify a significant or insignificant difference in N/S values by comparing ploidy and different modes of reproduction. However, only all the histone genes combined and H2B alone showed enough statistical evidence to reject the null hypothesis for both the Kruskal-Wallis and Wilcoxon tests. The fact that most of the histone genes individually did not show statistical evidence for a difference in N/S distributions contradicts the results when all the histone genes are combined. I believe to fully understand this inconsistency that more sampling is needed. Either more lineages need to be sampled or multiple snails from each lineage are needed. I imagine this would give more statistical power to both these tests. While the methods and results of this project may not be conclusive enough to determine selection occurring among repetitive histone genes in Potamopyrgus antipodarum, they do suggest some possible underlying distinction between asexual and sexual repetitive histone genes and further investigation is appropriate.  
