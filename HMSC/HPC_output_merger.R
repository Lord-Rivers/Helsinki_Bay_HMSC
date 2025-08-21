##This script merges the chains output by the hmsc-hpc fitting. ##
#It contains a now unnecessary check to address a bug in a previous version of
#the hmsc-hpc code that led to Alpha not always being output as a matrix.

### Loading packages####
remove(list=ls())
require(cli)
require(jsonify)
require(Hmsc)

### Set up directories ####
model_description = "PA_TSP_AT_d_dp_yr_PTCYC"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
localDir = "./Hmsc Outputs"
ModelDir = file.path(localDir, sprintf("%s/Models",model_description))

### Model name and samples and thining lists ####
#This method globs all the files in the RAW_HPC allowing for merging without having to manually defining the sample and thin values of the target files.
temp = Sys.glob(file.path(ModelDir,"/RAW_HPC/HPC_samples*_thin*_chain_0.rds"))
samples_list = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?rds)","\\2",temp))
thin_list = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?rds)","\\2",temp))
# verbose = 10 #This was for now unnecessary alpha fix outputs
#Loading in the unfitted model
load(file = file.path(ModelDir, "unfitted/unfitted_models.RData"))

for (x in 1:length(samples_list)) {
  samples = samples_list[x]
  thin = thin_list[x]
  filename = sprintf("/Raw_HPC/HPC_samples_%.4d_thin_%.2d_chain_*.rds", samples, thin)
  in_files = Sys.glob(file.path(ModelDir,filename))
  chains = length(in_files)
  out_file = file.path(ModelDir,sprintf("/Fitted/HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata", samples, thin, chains))
  cli_h1("Model {model_description}")
  cli_alert_info("Samples: {samples} \t Thin: {thin}")
  if(file.exists(out_file)){
    cli_alert_success("HPC outputs have already been converted")
    cli_rule()
  } else if(file.exists(in_files[1])){
    check = TRUE
    chainList = vector("list", chains)
    fit_times = vector("list", chains)
    ptm0 = proc.time()
    for(y in 1:chains){
      #There is an "extra chain", this is the part of the output that saves the
      #model run time
      cli_progress_step("Reading in file {y}")
      temp = from_json(readRDS(file = in_files[y])[[1]])
      cli_process_done()
      # #There was an issue with how the HPC code outputs the Alpha values, this
      # #converts the matrix it outputs into the expected vectors, a vector per
      # #random factor (nr) each with elements for each latent variable (nf)
      # if(is.matrix(temp[[1]][[1]]$Alpha)){
      #   cli_alert_info("Alpha is a matrix")
      #   cli_progress_step("Fixing alpha issue")
      #   for(i in 1:samples){
      #     temp_Alpha_Mat = temp[[1]][[i]]$Alpha
      #     temp[[1]][[i]]$Alpha = lapply(seq_len(nrow(temp_Alpha_Mat)), function(p) temp_Alpha_Mat[p,])
      #   }
      #   cli_process_done()
      # } else {
      #   cli_alert_info("Alpha is not a matrix\nNo fix required")
      # }
      chainList[[y]] = temp[[1]]
      fit_times[y] = temp[[2]]
    }
    Reading_time = proc.time() - ptm0
    cli_progress_step("Importing and merging Posteriors\n")
    ptm1 = proc.time()
    transient = ceiling(samples*thin/2)
    posteriors = importPosteriorFromHPC(models[[1]], chainList, samples, thin, transient, alignPost = TRUE)
    filename_in = file.path(ModelDir,sprintf("INITS/HPC_INIT_samples_%.4d_thin_%.2d_chains_%.1d.rds",samples,thin,chains))
    init_model = from_json(readRDS(file = filename_in)[[1]])
    posteriors$adaptNf = init_model$adaptNf
    posteriors$verbose = init_model$verbose
    import_time = proc.time() - ptm1
    cli_process_done()
    total_time = proc.time() - ptm0
    timming = cbind(Reading_time[1:3],import_time[1:3],total_time[1:3])
    colnames(timming) = c("Reading","Import","Total")
    noquote(timming)
    cat("Timing:\n")
    cat(sprintf("\tReading\tImport\tTotal\nUser\t%.3fs\t%.3fs\t%.3fs\nSystem\t%.3fs\t%.3fs\t%.3fs\nTotal\t%.3fs\t%.3fs\t%.3fs\n",Reading_time[1],import_time[1],total_time[1],Reading_time[2],import_time[2],total_time[2],Reading_time[3],import_time[3],total_time[3]))
    cli_progress_step("Saving file")
    fitted_model = list(posteriors,fit_times)
    names(fitted_model) = c("posteriors","fit_times")
    save(fitted_model, file = out_file)
    cli_progress_done()
    cli_rule()
  } else {
    check = FALSE
    cli_alert_warning("RDS files could not be found")
    cli_rule()
  }
}
