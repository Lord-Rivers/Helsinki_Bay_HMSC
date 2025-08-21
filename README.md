# Read Me
## File structure
These files rely on a set file structure by defualt, as such a working example for the full model has been set up. The file structure used is as follows
- Data
- HMSC
	- Hmsc outputs
	  - <Model_Name>
		- Models
		  - Fitted
		  - INITS
		  - Raw_HPC
		  - Temp
		  - Unfitted
		- Results
		  - Preds

The data on the HPC uses a similar structure, this is evident in the HPC scripts folder. The base scripts write there results to a folder called Helsinki_Bay which is set up as follows:
- Helsinki_Bay
	- <Model_Name>
		- Fitted
		- INIT
		- Sampled

## Scripts
The scripts are set up to run the base example, the user must manually change the model names at the start of each script. All outputs apart from predictions have been generated for a model run of 500 samples thinned by 500. The model outputs for 750 by 750 where too large for github, as where the predictions. It should take around 9 minutes to generate the predictions for the PA and Con model, these are saved so that they are not recalculated during the construction of the hurdle model.
