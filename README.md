
![Illustration du plan D-optimal](Graphiques/DOPEcal.jpg)

# Plan d'experiences pour la calibration de codes de calcul couteux (à sorite scalire)
Ce repertoire contient des criteres de selection de plan d'experiences physiques pour la calibration d'un code de calculs couteux à sortie scalaire.
Le cadre statistique utilisé est celui de [Kennedy et O'Hagan (2001)](https://www.asc.ohio-state.edu/statistics/comp_exp/jour.club/kennedy01.pdf).
L'objectif est selectionner le plan d'experiences physiques par par optimisation d'un critère d'optimalité

$$X_{opt} \in \text{ argmax } C_{opt}(X) \text{ ou } X_{opt} \in \text{ argmin } C_{opt}(X)$$

# Description
Les critères suivants sont disponible dans ce repos : 
a) Critères bayesiens : deux critères basés sur la densité a posteriori 

(somme des variances a posteriori et erreur quadratique moyenne a posteriori et la divergence de Kullback-Leibler entre la densité a priori et la densité 
a posteriori)
b)  Critères bayesiens lineaires : une version des criteres bayesiens correspondant au cas particulier de codes de calcul lineaire.
c) Critères alphabeitques : le critere ED-optimalité et le critere ET-optimalité.

*NB : Ces criteres sont calculés à l'aide d'un émulateur de processus gausien (*model*) à fournir en entrée.*

# Installation <a name="Install"></a>
Pour 
**``git clone https://github.com/TheseAdama/DOPE.git``**

# Package R : 
Executer ce code pour installer les packages R necessaire : 

# Reference
