# DOPE : Design of Physical Experiments for scalar output expensive computer code.

![DOPE](/DOPECal.jpg)

## Table 
1. [Introduction](#Introduction)
    * [Contexte : choix de plan d'experiences physiques](#DOE)
    * [Cadre statistique de KOH](#CadreKOH)
    * [Objectifs](#Goal)
2. [Critères d'optimalités](#Copt)
    * [Critères bayesiens](#Cbayes)
    * [Critères alphabetiques](#Calpha)
    * [Critères bayesiens (cas de code linéaire)](#CbayesLin)
3. [Examples](#Examples)
    * [Exemple en 2D](#2D)
    * [Exemple en 2D](#2D)
4. [Installation](#Install)
5. [References](#References)

## Introduction <a name="Introduction"></a>

### Contexte : choix de plan d'experiences physiques <a name="DOE"></a>
Dans le cadre de la calibration d'un code de calculs couteux, la selection des experiences physiques est important pour avoir un maximum d'information
sur les parametres incertains.

### Cadre statistique de KOH <a name="CadreKOH"></a>
Le cadre statistique utilisé est celui de [Kennedy et O'Hagan (2001)](https://www.asc.ohio-state.edu/statistics/comp_exp/jour.club/kennedy01.pdf) avec ...

### Objectifs <a name="Goal"></a>
Selection le plan d'experiences physiques ($X_{opt}$) par optimisation d'un critere d'optimalité

$$X_{opt} \in \text{ argmax } C_{opt}(X) \text{ ou } X_{opt} \in \text{ argmin } C_{opt}(X)$$

## Critères d'optimalités <a name="Copt"></a>
Trois types de criteres sont codés.
### Critères bayesiens <a name="Cbayes"></a>
DAP
### Critères alphabetiques <a name="Calpha"></a>
Litterature
### Critères bayesiens (cas de code linéaire) <a name="CbayesLin"></a>
Cas particulier de code de calculs lineaire ou linearisable.

## Installation <a name="Install"></a>
Pour 
**``git clone https://github.com/TheseAdama/DOPE.git``**

## References
