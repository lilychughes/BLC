# BLC

## Branch length correlation of gene trees and species trees

This is a test to filter contamination in phylogenomic datasets based on [Simion et al. (2017) in *Current Biology*](https://www.sciencedirect.com/science/article/pii/S0960982217301999). The code here was used as part of the quality control pipeline by [Arcila et al. (2021) in *Systematic Biology*](https://academic.oup.com/sysbio/article-abstract/70/6/1123/6204118).

There are four main steps:
1. Infer a Maximum Likelihood (ML) tree for all your concatenated genes
2. Create a set of contraint trees from the concatenated ML topology with only the taxa present in each gene
3. Infer ML gene trees that are constrained to the concatenated ML topology
4. Correlate the terminal branch lengths for each taxon in the constrained gene tree with the concatenated ML tree

Gene sequences that are contaminated, paralogous, misassembled or misaligned, should have excessively long branch lengths in the constrained gene trees. These sequences are then flagged and removed from downstream analyses. 
