## Make model predictions ##

# Warning, this takes a lot of memory and can take 10-20 minutes to run locally
# depending on the size of the converged model. An HPC compatible version can be
# found in the HPC folder along with the .sh script for calling it

#This script automatically detects the longest model run, the values of
#samples_list and thin_list can be manually set to force it select a different
#model if required. 

### Load packages ####
remove(list=ls())
require(Hmsc)
require(ggplot2)
require(cli)
set.seed(369)

#### Define inputs ####
model_description = "PA_TSP_AT_d_dp_yr_PTCYC"

#Define a list of environmental variables to make predictions for, these must
#match the variables in the model exactly.
#env.list = NULL
env.list = c("Temp","Ptot","Salin","jday") #This list will keep the prediction objects at a manageable size

nChains = 4
#### Set up directories ####
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))

localDir = shortPathName(sprintf("./Hmsc Outputs/%s",model_description))
ModelDir = shortPathName(file.path(localDir, "Models/Fitted"))
UnfittedDir = shortPathName(file.path(localDir, "Models/Unfitted"))
ResultDir = shortPathName(file.path(localDir, "Results"))

#If the Preds folder doesn't exist create it
if(!dir.exists(file.path(ResultDir,"Preds"))){
  dir.create(file.path(ResultDir,"Preds"))
}

### Read in data ####
load(file = file.path(UnfittedDir,"unfitted_models.RData"))

file_list = Sys.glob(file.path(ModelDir,sprintf("/HPC_samples_*_thin_*_chains_%.1d.Rdata", nChains)))
samples_list = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",file_list))
thin_list = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?Rdata)","\\2",file_list))
order = sort(samples_list*thin_list, index.return = TRUE)
file_list = file_list[order$ix]
rm(samples_list,thin_list,order)

filename = dplyr::last(file_list)

load(filename)
m = fitted_model$posteriors
rm(fitted_model)
modelnames = model_description

#### Prepare the env.list variable ####
if(is.null(env.list)){
  env.list = list()
  env.list = 0
}

if(all(env.list==0)){
  if(m$XFormula=="~."){
    covariates = colnames(m$XData)
  } else {
    covariates = all.vars(m$XFormula)
  }
} else {
  covariates = env.list
}

#### Calculate the gradients and predictions ####
if(length(covariates)>0){
  outfile = file.path(ResultDir,sprintf("Preds/Preds_%s_HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",model_description, m$samples, m$thin, nChains))
  cli_h1("Making predictions")
  if(file.exists(file.path(outfile))){
    cli_alert_success("Predictions already calculated")
    load(outfile)
  } else {
    Preds = vector("list", length(covariates))
    ptm_tot = proc.time()
    for(k in 1:(length(covariates))){
      covariate = covariates[[k]]
      ptm = proc.time()
      cli_h2("Calculating predictions for {covariate}")
      cli_alert_info("Starting at: {format(Sys.time(), '%H:%M:%S')}")
      cli_progress_step("Starting to construction Gradient:")
      #Don't set jday or year here since they are linearly calculated based off
      #the focalVariable
      Gradient = constructGradient(m,focalVariable = covariate, ngrid=50)
      cli_progress_step("Starting to construction Gradient2:")
      #Set jday to 211 (30th July) for all gradients where jday is not the
      #focalVariable
      Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = list("jday"=list(3,211),"year"=list(3,1988)), ngrid=50)
      
      cli_progress_step("Making predictions based on Gradient 1:")
      #useSocket should be false on all OS's other than windows
      predY = predict(m, Gradient=Gradient, expected = TRUE, nParallel = 4, useSocket = TRUE)
      cli_progress_step("Making predictions based on Gradient2:")
      predY2 = predict(m, Gradient=Gradient2, expected = TRUE, nParallel = 4, useSocket = TRUE)
      step_time = lubridate::duration((proc.time() - ptm)[3])
      cli_progress_cleanup()
      cli_alert_info("Step time: {.strong {step_time}}")
      Preds[[k]]$predY = predY 
      Preds[[k]]$predY2 = predY2 
      Preds[[k]]$Gradient = Gradient
      Preds[[k]]$Gradient2 = Gradient2
      #Clean up to prevent excessive memory use
      rm(predY,predY2,Gradient,Gradient2)
      gc()
    }
    names(Preds) = covariates
    save(Preds, file = outfile)
    computational.time = lubridate::duration((proc.time() - ptm_tot)[3])
    cli_alert_success("Total Time taken: {.strong {computational.time}}")
  }
}
