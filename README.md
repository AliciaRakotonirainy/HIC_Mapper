Ressources très utiles :

- L'article de Mathieu Carrière qui décrit exactement le travail qu'on est en train de faire ([article](https://www.biorxiv.org/content/10.1101/486936v1.full))
- Bonne ressource qui réexplique bien ce qu'est Mapper et ses paramètres ([blog](https://medium.datadriveninvestor.com/the-mapper-algorithm-d0842f926658))
- Article donnant la définition du SCC ([article](https://pubmed.ncbi.nlm.nih.gov/28855260/))

## Résumé de l'article "HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient"

They define the SCC, a new new similarity measure, that quantifies the similarity between Hi-C interaction matrices (= contact matrices = contact maps).

Why is it useful ? 
- Assess the reproducibility of replicate samples (2 biologcial replicates should have a high similarity)
- Quantify the distance bewteen Hi-C matrices from different cell types

Hi-c data have certain unique characteristics, that cannot be taken into account using only Pearson/Spearman correlation between 2 Hi-C matrices:
- Specific domain structure : Hi-C matrices present contiguous regions in which loci tend to interact more frequently. Those regions are stable across cell types (even if the particular interactions whithin the regions are variable between cell types). But Pearson/Spearman correlation only consider point interactions and do not take region-scale information into account! So Spearman correlation can be driven to low values by the stochastic variation in the point interactions and overlook the similarity in domain structures.
- Distance dependence : chromatin interaction between 2 loci decrease substantially as their genomic distance increases. This pattern is generally thought to result from nonspecific (= silent) interactions, which are more likely to occur between closer genomic distances. This pattern results in high Pearson correlation between any two Hi-C matrices, even when the samples are unrelated. Therefore, the Pearson correlation cannot distinguish real biological replicates from unrelated samples!


Steps to calculate the SCC : 
- Minimizes the effect of noise by smoothing the Hi-C matrix; it makes the domain structures more visible. For smoothing, they apply a 2D mean filter.
- Addresses the distance-dependence effect by stratifying Hi-C data according to their genomic distance ; they cut the matrices into stratas with similar genomic distance
- SCC : Calculate Pearson correlation for each stratum ; and then aggregating all Pearson correlations coeffs using a weighted average. The weights are derived from the generalized Cochran-Mantel-Haenszel (CMH) statistic (the exact expression is given in the Methods part)

## Résumé de notre projet

On veut représenter nos cellules sous forme de Mapper Graph. Le but est de voir que les cellules se regroupent en fonction de leur stade de cycle cellulaire ; on devrait voir apparaître une grande boucle.  

1. Calculer la matrice de distances deux-à-deux pour chaque cellules. ⚠️ Le SCC est un score de similarité. On le convertit donc en une distance avec la formule : $d(X,Y) = \sqrt{SCC(X,X) + SCC(Y,Y) - 2SCC(X,Y)}$ (produit scalaire -> norme)

2. On calcule un Kernel PCA sur cette matrice de distances. Les 2 premières composantes sont notre "lens" = "filter".

3. On calcule le graph Mapper à partir de notre matrice de distance et du filter qu'on a choisi.

### TO DO next 

Hyperparameter tuning : 
- [ ] Try various Mapper and SCC parameters, and interpret their influence on the Mapper shape
- [X] Find a set of parameters fro which the cell cycle can be detected as a big loop in the Mapper 
- [ ] Quantify the cell cycle statistical robustness
- [ ] Compare the results with basic dimensionality reduction on the raw contact maps