# ReliabilityMOG: Stress-strength reliability estimation of Marshall-Olkin G distributions

This is a source code of **Statistical inference on stress-strength reliability for the Marshall-Olkin family of distributions (2024)** by Gülce Cüran and Fatih Kızılaslan.

In this study, we investigate the stress-strength reliability of the Marshall–Olkin model when the stress and strength components are non-identical, following not only the same baseline distribution but also different baseline distributions, such as Burr Type XII, Chen, and Weibull. 

Point and interval estimates of the reliability are obtained using maximum likelihood and Bayesian approaches under various settings, including different distribution and parameter choices. 

We analyze an original hydrological dataset from Istanbul, Türkiye, to demonstrate the practical application of our methods. 
This dataset is retrieved from the open data source of the Istanbul Metropolitan Municipality ([IMM Open Data](https://data.ibb.gov.tr/en/)).


## Citation 

 > Cüran, Gülce and Kızılaslan, Fatih (2024). **"Statistical inference on stress-strength reliability for the Marshall-Olkin family of distribution"**. arXiv DOI: [??](https://arxiv.org/??)


## R

It includes some of the source code necessary to run an example of real data application in Cüran and Kızılaslan (2024).

``functions.R`` includes the related functions in the algorithms.

``real_data_analysis.R`` includes the code for applying the methods to the arranged hydrological dataset.

``combined_full_data.xlsx`` is the combined full data from 
[*GENERAL\textunderscore DAM\textunderscore RESERVED\textunderscore WATER*](https://data.ibb.gov.tr/en/dataset/istanbul-dam-occupany-rates-data/resource/b68cbdb0-9bf5-474c-91c4-9256c07c4bdf) measurement (million cubic meters $m^3$) which represents the amount of water reserve and 
[*Istanbul gunluk tuketim*](https://data.ibb.gov.tr/en/dataset/istanbul-barajlarina-dusen-gunluk-toplam-yagis-miktari/resource/762b802e-c5f9-4175-a5c1-78b892d9764b) measurement ($m^3$) which represents daily consumption of Istanbul.

``real_data.xlsx`` is the arranged data in the real data application.
