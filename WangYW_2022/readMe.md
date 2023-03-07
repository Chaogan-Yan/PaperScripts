### 1. Harmo- Harmonization step  
	Linear regression - R 4.2.2  
	ComBat - Matlab2020a  
	SMA - Matlab2020a  
	ICVAE - Python3.6.5  

### 2. Efficiency(Residual site effect) - Evaluation step   
    mean and variance of residual site effect -TSP3/CoRR/FCP_SiteStats.m  
    covariance of residual site effect - randomforest_predictsite_LOO.m  
	
### 3. Identifiability - Evaluation step  
	To calculate the random chance of identifiability - ChanceLevelCalculation.m, stirling2.m   
	Evaluation - Identifiability.R	  
	
	notice：The Identifiability test for SMA target site choices is also in the   
	Identifiability.R file.  
	
### 4. Test-retest Reliability - Evaluation step  

### 5. Replicability - Evaluation step  
	statistics - FCP_MalevsFemale.m  
	calculate the dice of replicability - DICEall.m  
	
### 6. Misuse - Evaluation step  
	SexinSeparateSites.m  
	
### 7. SMA Target site Choice - Practice Experiments  
	Test-retest reliability of CoRR - CoRR-ChangeTS-TestRetestReliability  
		statistics - TargetSiteChoice.m  
		Calculate dice - Targetsitechoice_Dice.m  
	Check if the sample size or distribution of Z variates affect stability - FCP-Bootstrapping-Stability  
		sample size - bootstrappingExps.m  
		distribution - bootstrappingExps2.m  
#

### Harmonization tools 

#### **ComBat**   

Codes:  
https://github.com/Jfortin1/ComBatHarmonization  

Citations:  
1. Jean-Philippe Fortin, Drew Parker, Birkan Tunc, Takanori Watanabe, Mark A Elliott, Kosha Ruparel, David R Roalf, Theodore D Satterthwaite, Ruben C Gur, Raquel E Gur, Robert T Schultz, Ragini Verma, Russell T Shinohara. Harmonization Of Multi-Site Diffusion Tensor Imaging Data. NeuroImage, 161, 149-170, 2017  
2. Jean-Philippe Fortin, Nicholas Cullen, Yvette I. Sheline, Warren D. Taylor, Irem Aselcioglu, Philip A. Cook, Phil Adams, Crystal Cooper, Maurizio Fava, Patrick J. McGrath, Melvin McInnis, Mary L. Phillips, Madhukar H. Trivedi, Myrna M. Weissman, Russell T. Shinohara. Harmonization of cortical thickness measurements across scanners and sites. NeuroImage, 167, 104-120, 2018  
3. W. Evan Johnson and Cheng Li, Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics, 8(1):118-127, 2007.  

#### **CovBat**   

Codes:  
https://github.com/andy1764/CovBat_Harmonization  
Citation:  
Chen, A. A., Beer, J. C., Tustison, N. J., Cook, P. A., Shinohara, R. T., Shou, H., & The Alzheimer's Disease Neuroimaging Initiative (2022). Mitigating site effects in covariance for machine learning in neuroimaging data. Human Brain Mapping, 43(4), 1179–1195. https://doi.org/10. 1002/hbm.25688

#### **SMA** 

Codes:  
https://github.com/hzhoustat/PNAS2018  
Citation:  
https://www.pnas.org/doi/full/10.1073/pnas.1719747115  
