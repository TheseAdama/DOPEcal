
![Illustration du plan D-optimal](Graphiques/DOPEcal.jpg)

# Plan d'experiences pour la calibration de codes de calcul couteux (à sorite scalire)
Ce repertoire contient des criteres de selection de plan d'experiences physiques pour la calibration d'un code de calculs couteux à sortie scalaire.
Le cadre statistique utilisé est celui de [Kennedy et O'Hagan (2001)](https://www.asc.ohio-state.edu/statistics/comp_exp/jour.club/kennedy01.pdf).
L'objectif est selectionner le plan d'experiences physiques par par optimisation d'un critère d'optimalité

$$X_{opt} \in \text{ argmax } C_{opt}(X) \text{ ou } X_{opt} \in \text{ argmin } C_{opt}(X)$$

# Description
Les critères suivants sont disponible dans ce repos : 
- Critères bayesiens : ces critères sont basés sur la densité a posteriori. La fonction
  
 CoptBayes <- function(X, model, Dx, Dtheta, sigeps,  dprior, rprior, L=100, K=1000, typeCopt='KL', type='SK',...)
 
  - typeCopt='SOV' : somme des variances a posteriori.
  - typeCopt='MSE' : erreur quadratique moyenne a posteriori.
  - typeCopt='KL'  :  la divergence de Kullback-Leibler entre la densité a priori et la densité a posteriori.
    
- Critères bayesiens lineaires : une version des criteres bayesiens correspondant au cas particulier de codes de calcul lineaire et de densité a priori gaussienne.
  
  CoptBayesLin <- function(X, model, Dx, Dtheta,  thetaprior, Sigmaprior, sigeps, L=1000, typeCopt='KL', type='SK',...)
  
- Critères alphabeitques : le critère ED-optimalité et le critere ET-optimalité.
  
    CoptMFisher <- function(X, model, sigeps, Dx, Dtheta, dftheta, typeCopt="Det", L=1000, type='SK')
  
- Critere glouton : ce critere utilise la variation du code de caluls et la repartition du plan d'experiences dans l'espace experimental.

  CVMm <- function(X, model, Dx, Dtheta, L=100, type="SK", alpha=0.5)
  
*NB : Ces critères sont calculés à l'aide d'un émulateur de processus gausien (*model*) à fournir en entrée.*

# Installation <a name="Install"></a>
Cloner ce depot pour telecharger les fichiers en local.
**``git clone https://github.com/TheseAdama/DOPE.git``**
Ou telecharger directement le fichier ZIP sur github.

# Package R : 
Executer ce code pour installer les packages R necessaire : 

# Reference
**Adama Barry, François Bachoc, Sarah Bouquet, Miguel Munoz Zuniga, Clémentine Prieur.Optimal Design ofPhysical and Numerical Experiments for Computer Code Calibration. 2024.〈hal-04615127v2〉)**
