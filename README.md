# Allometric_marine_reserves

This repository contains the files needed to recreate the results contained in:

Marshall, D. J., Bode, M., Mangel, M., Arlinghaus, R., & Dick, E. J. (2021). Reproductive hyperallometry and managing the world’s fisheries. *Proceedings of the National Academy of Sciences*, **118**(34).

[* Note that not all files are yet uploaded *]

The file labelled "Pike_IsometricEffort_vs_Hyperallometric_effort.m" is a set of matlab functions that calculate the fishing effort at MSY (and the associated maximum sustainable yield) for comparable populations with hyperallometric and isometric reproduction. 

The TIFF image "MSY_comparison.tiff" plots the results of the pike analysis.

The file labelled "Central_command.m" contains matlab code for predicting the combined benefits of MPAs and harvest effort with explicit incorporation of nonlinear allometry (i.e., hyperallometric). This code recreates a series of results and figures relevant to the article. It is a spatial model of a harvested stock where a proportion of the habitat may be protected by a no-take marine reserve, and where the remaining proportion is exposed to effort-controlled harvesting. 

Note, the root directory should include an empty folder called "Figures".
