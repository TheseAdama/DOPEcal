![Illustration du plan D-optimal](Graphiques/DOPEcal.jpg)

# Plan d'expériences pour la calibration de codes de calcul coûteux (à sortie scalaire)

Ce dépôt contient des critères de sélection de plans d'expériences physiques pour la calibration d'un code de calculs coûteux à sortie scalaire. Le cadre statistique utilisé est celui de [Kennedy et O'Hagan (2001)](https://www.asc.ohio-state.edu/statistics/comp_exp/jour.club/kennedy01.pdf). L'objectif est de sélectionner le plan d'expériences physiques par optimisation d'un critère d'optimalité :

$$
X_{\text{opt}} \in \operatorname{argmax} C_{\text{opt}}(X) \quad \text{ou} \quad X_{\text{opt}} \in \operatorname{argmin} C_{\text{opt}}(X)
$$

## Description

Les critères suivants sont disponibles dans ce dépôt :

- **Critères bayésiens** : ces critères, basés sur la densité a posteriori, sont implémentés par la fonction suivante :

    ```r
    CoptBayes <- function(X, model, Dx, Dtheta, sigeps, dprior, rprior, L = 100, K = 1000, typeCopt = 'KL', type = 'SK', ...)
    ```
    
    - `typeCopt = 'SOV'` : somme des variances a posteriori ;
    - `typeCopt = 'MSE'` : erreur quadratique moyenne a posteriori ;
    - `typeCopt = 'KL'`  : divergence de Kullback-Leibler entre la densité a priori et la densité a posteriori.

- **Critères bayésiens linéaires** : une version des critères bayésiens applicable pour les codes de calcul linéaires et les densités a priori gaussiennes.

    ```r
    CoptBayesLin <- function(X, model, Dx, Dtheta, thetaprior, Sigmaprior, sigeps, L = 1000, typeCopt = 'KL', type = 'SK', ...)
    ```

- **Critères alphabétiques** : ce groupe inclut le critère ED-optimalité ainsi que le critère ET-optimalité.

    ```r
    CoptMFisher <- function(X, model, sigeps, Dx, Dtheta, dftheta, typeCopt = "Det", L = 1000, type = 'SK')
    ```

- **Critère glouton** : ce critère exploite la variation du code de calcul et la répartition du plan d'expériences dans l'espace expérimental.

    ```r
    CVMm <- function(X, model, Dx, Dtheta, L = 100, type = "SK", alpha = 0.5)
    ```

*NB : Ces critères sont calculés à l'aide d'un émulateur de processus gaussien (*model*) à fournir en entrée.*

## Installation <a name="Install"></a>

Clonez ce dépôt pour télécharger les fichiers en local :

```bash
git clone https://github.com/TheseAdama/DOPE.git
```
Vous pouvez également télécharger directement le fichier ZIP depuis GitHub.

# Package R : 
Exécutez le code suivant dans R pour installer les packages nécessaires : 
 ```r
install.packages(c("Dicekriging", " ", " "))
 ```

# Reference
**Adama Barry, François Bachoc, Sarah Bouquet, Miguel Munoz Zuniga, Clémentine Prieur.Optimal Design ofPhysical and Numerical Experiments for Computer Code Calibration. 2024.〈hal-04615127v2〉)**
