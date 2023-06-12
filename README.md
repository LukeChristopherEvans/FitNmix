# FitNmix
Luke Christopher Evans 2023

This is an archive of testing code for fitting n-mixture models to repeated count data. It was designed
to be shared with staff in the RSPB and Bat Conservation Trust, so many of the comments will not necessarily seem relevant 
to a new viewer, but I have shared it as it could be a useful starting point for a student wishing to apply n-mixture models to count data using either stan or nimble. It is also only development code, so use appropriate caution when adapting any of the code for your own uses and
watch out for bugs in the comments - I have been moving towards using comments more sparingly as they often contain 'bug's that never get checked,
but this code has a lot of comments. 
Any comments or queries can be addressed to lukechristopher.evans@reading.ac.uk

Much of the data used in the scripts is not available (RSPB, BCT, ECN) as data should be requested directly
from the relevant sources which would be RSPB/BTO/JNCC or BCT. However, the ECN data for bats and birds can be freely downloaded from
https://eidc.ac.uk/ and R code for putting these into the format required for the rest of the code is contained in the EcnTest folder. 
There is also code for simulated datasets, that can be run immediately.

I have organised the code using relative paths, so that once the folder path is set to the folder on your own machine, 
the rest of the code should run (as long as any datasets required are added the data folder with the correct filename). 
