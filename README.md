# Talk on Bayesian Sequential Design given @ MRC CBU 29 July 2020
In this repository you'll find all important files talks. 

The following things can be found here:
* The __presentation__ was created with _presentation.Rmd_ and can be viewed with _presentation.html_. The __figures__ for this presentation were created with _figures.R_ and are stored in the folders _figures_.
* The rslurm jobs were created with the _simulationScript.R_, which were then submitted to the cluster. 
* The __data frames__ that are used in the presentation and for the figures can be found in _simulationResults.RData_. This file was generated with the script _getRslurmResults.R_.
* The __raw data__ from them the simulations can be found in the in the folders starting with \_rslurm\_. Note that the last two were not featured in the talk:
	* Traditional Design (\_rslurm\_traditionalDesign),
	* Sequential Design (\_rslurm\_SequentialDesignWithoutLimit),
	* Sequential Design with upper limit (\_rslurm\_SequentialDesignWithtLimit),
	* Sequential Design with upper limit and batch size of 1 (\_rslurm\_SequentialDesignWithtLimit_batcheSize1),
	* Sequential Design with upper limit and start sample size of 2 (\_rslurm\_SequentialDesignWithtLimit_lowMinN).