#################################################
# Here we fit the models defined for the system
################################################
remove(list=ls())
require(Hmsc)
require(jsonify)

model_description = "PA_TSP_AT_d_dp_yr_PTCYC"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

localDir = sprintf("./Hmsc Outputs/%s",model_description)
ModelDir = file.path(localDir, "Models")

### Read in the unfitted models ####
load(file = file.path(ModelDir, "Unfitted/unfitted_models.RData"))

#Set the sample and thin for the models that are being generated
#samples_list = c(100,100,250,500,750)
#thin_list = c(10,100,250,500,750)
samples_list = c(500,750)
thin_list = c(500,750)
nChains = 4
nParallel = 4
Lst = 1
verbose = 1

while(Lst <= length(samples_list)){
  #Seed each different run length with the same starting values for
  #reproducible
  set.seed(369)
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  transient = ceiling(samples*thin/2)
  
  filename = file.path(ModelDir,sprintf("INITS/HPC_INIT_samples_%.4d_thin_%.2d_chains_%.1d.rds",samples,thin,nChains))
  #if you are fitting the model using R and not the HPC, change engine to R,
  #skip HPC merge step You'll have to change the output of this step accordingly
  m = sampleMcmc(models[[1]], samples = samples, thin=thin,
                 transient = transient,
                 nChains = nChains,
                 verbose = verbose, 
                 engine = "HPC") 
  
  saveRDS(to_json(m), filename)
  Lst = Lst + 1
}

