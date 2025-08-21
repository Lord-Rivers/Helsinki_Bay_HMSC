remove(list=ls())
require(Hmsc)
require(ggplot2)
require(cli)
set.seed(369)

#The names of the arguments aren't important in the R script, only the order. The names will be defined in the batch script which will ensure that regardless of the order that they are parsed in the terminal they will be parsed to R in the order area, model_description
args = commandArgs(trailingOnly=TRUE)

area = args[1]
model_description = args[2]

if(dir.exists("/home/jdehaast")){
  setwd(sprintf("~/HMSC_Data/%s", area))
  ModelDir = file.path(sprintf("./%s/Fitted",model_description))
  UnfittedDir = file.path(sprintf("./%s/Temp",model_description))
  TempDir = file.path(sprintf("./%s/Temp",model_description))
  ResultDir = file.path(sprintf("./%s/Temp",model_description))
}else{
  cat("Can't find HPC directoy.\nExcution may be slow on personal computers\n")
  setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),"../"))
  
  localDir = sprintf("./HMSC/Hmsc Outputs/%s",model_description)
  ModelDir = file.path(localDir, "Models/Fitted")
  UnfittedDir = file.path(ModelDir, "Unfitted")
  ResultDir = file.path(localDir, "Results")
}

nChains = 4
env.list = NULL
#This is done because the list of focus variables in this code set up must be as
#long as the list of different models. This code is for when you want to ask for
#the same env response variable form each model. Note that it should not be used
#if the different models have different X variables in them. Here 3 is used as 3
#models were defined for each run env.list = list()
#for(x in 1:3){
#  env.list[[x]] = c("Temp","Salin")
#}
if(!dir.exists(file.path(ResultDir,"Preds"))){
  dir.create(file.path(ResultDir,"Preds"))
}

load(file = file.path(UnfittedDir,"unfitted_models.RData"))

file_list = Sys.glob(file.path(ModelDir,sprintf("/HPC_samples_*_thin_*_chains_%.1d.Rdata", nChains)))
samples_list = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",file_list))
thin_list = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?Rdata)","\\2",file_list))
order = sort(samples_list*thin_list, index.return = TRUE)
file_list = file_list[order$ix]
rm(samples_list,thin_list,order)

filename = dplyr::last(file_list)

#if(file.exists(filename)){
load(filename)
m = fitted_model$posteriors
rm(fitted_model)
modelnames = model_description

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

#Change this to save the gradients, check if the file exists and if it does skip calculating them and move straight to plotting
if(length(covariates)>0){
  outfile = file.path(ResultDir,sprintf("Preds/Preds_%s_HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",model_description, m$samples, m$thin, nChains))
  cli_h1("Making predictions")
  if(file.exists(file.path(outfile))){
    #cat("----------------\nPredictions already calculated\n")
    cli_alert_success("Predictions already calculated")
    load(outfile)
  } else {
    Preds = vector("list", length(covariates))
	ptm_tot = proc.time()
    for(k in 1:(length(covariates))){
      covariate = covariates[[k]]
      #cat(sprintf("----------------\nPlotting predictions for Covariate: %s\n", covariate))
      cli_h2("Calculating predictions for {covariate}")
      cli_progress_step("Starting to construction Gradient:")
      #ptm = proc.time()
      Gradient = constructGradient(m,focalVariable = covariate, ngrid=50)
      #computational.time = proc.time() - ptm
      #cat(sprintf("Time taken for this step: %.2f s\n\nStarting on construction Gradient 2:\n", computational.time[3]))
      
      cli_progress_step("Starting to construction Gradient2:")
	  #Set jday to 211 (30th July) for all gradients where jday is not the focalVariable, there is no need for an if statement
	  #Naturally this code has been updated for the Helsinki Bay models now
      Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = list("jday"=list(3,211),"year"=list(3,1988)), ngrid=50)
      #computational.time = proc.time() - ptm
      #cat(sprintf("Time taken for this step: %.2fs\n\nMaking predictions based on Gradient 1\n", computational.time[3]))
      
      cli_progress_step("Making predictions based on Gradient 1:")
      predY = predict(m, Gradient=Gradient, expected = TRUE, nParallel = 8, useSocket = FALSE)
      #computational.time = proc.time() - ptm
      #cat(sprintf("Time taken for this step: %.2fs \n\nMaking predictions based on Gradient2\n",computational.time[3]))
      cli_progress_step("Making predictions based on Gradient2:")
      
      ptm = proc.time()
      predY2 = predict(m, Gradient=Gradient2, expected = TRUE, nParallel = 8, useSocket = FALSE)
      #cli_process_done()
      #computational.time = proc.time() - ptm
      #cat(sprintf("Time taken for this step: %.2fs \n", computational.time[3]))
      Preds[[k]]$predY = predY 
      Preds[[k]]$predY2 = predY2 
      Preds[[k]]$Gradient = Gradient
      Preds[[k]]$Gradient2 = Gradient2
    }
    names(Preds) = covariates
    save(Preds, file = outfile)
    computational.time = proc.time() - ptm_tot
    cat(sprintf("Total Time taken: %.2f s \nCurrent time: %s\n\n", computational.time[3],format(Sys.time(), "%H:%M:%S")))
  }
}
#}
