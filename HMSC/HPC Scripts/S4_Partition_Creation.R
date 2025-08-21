set.seed(369)
require(Hmsc)
require(jsonify)
require(dplyr)

#The names of the arguments aren't important in the R script, only the order. The names will be defined in the batch script which will ensure that regardless of the order that they are parsed in the terminal they will be parsed to R in the order area, model_description
args = commandArgs(trailingOnly=TRUE)

area = args[1]
model_description = args[2]

if(dir.exists("/home/jdehaast")){
  setwd(sprintf("~/HMSC_Data/%s", area))
  ModelDir = file.path(sprintf("./%s/Fitted",model_description))
  TempDir = file.path(sprintf("./%s/Temp",model_description))
}else{
  cat("Can't find HPC directoy.\nExcution may be slow on personal computers\n")
  setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),"../"))
  
  localDir = "./HMSC/Hmsc Outputs"
  ModelDir = file.path(localDir, sprintf("Models/%s/Fitted",model_description))
  TempDir = file.path(localDir,sprintf("Models/%s/Temp",model_description))
}

#Create the Temp directory if it does not already exsist.
 if(!dir.exists(TempDir)){
	dir.create(TempDir)
 }
 
nChains = 4
nfolds = 2

file_list = Sys.glob(file.path(ModelDir,sprintf("HPC_samples_*_thin_*_chains_%.1d.Rdata",
                                                nChains)))
samples_list = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",file_list))
thin_list = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?Rdata)","\\2",file_list))
order = sort(samples_list*thin_list, index.return = TRUE)
file_list = file_list[order$ix]
rm(samples_list,thin_list,order)

#Only run for the longest model run
filename.in = last(file_list)
samples = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",filename.in))
thin = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?Rdata)","\\2",filename.in))
message("Input file:", filename.in)
filename.out = file.path(ModelDir,sprintf("MF_samples_%.4d_thin_%.2d_chains_%.1d_nfolds_%.1d.Rdata",
                                          samples, thin, nChains,nfolds))

if(file.exists(filename.out)){
  cat("---------------------------\nthin = ",as.character(thin),
      "; samples = ",as.character(samples),
      "\nmodel fit had been computed already\n")
  break
} else {
  cat("---------------------------\nthin = ",
      as.character(thin),"; samples = ",as.character(samples),
      "\n",date(),"\n")
  load(file = filename.in) #models
  hM = fitted_model$posteriors
  if(is.null(hM$dfPi)){
    cat("No study Design defined, setting to NA's due to bug in function that requires it")
    hM$dfPi = data.frame(temp = factor(rep(NA, dim(hM$Y)[1])))
  }
  
  rm(fitted_model)
  
  partition = createPartition(hM, nfolds = nfolds)
  thin = 1
  start = 1
  ## STAGE 1: Basic housekeeping
  parts <- sort(unique(partition))
  nfolds <- length(parts)
  if (thin > 1 || start > 1)
    postN <- sum(sapply(hM$postList, function(z)
      length(seq(from=start, to=length(z), by=thin))))
  else
    postN <- Reduce(sum, lapply(hM$postList, length))
  ## output array
  predArray <- array(NA, c(hM$ny, hM$ns, postN))
  ## STEP 1: define new Hmsc model for each nfolds partition
  ## only implement for the simple case first
  if (is.list(hM$X))
    stop("not yet implemented for a list of model matrices")
  ## Pack setting training model into one function
  setHmsc <- function(k, hM) {
    train <- k != partition
    ## X can be a list of matrices
    XTrain <- if(is.matrix(hM$X))
      hM$X[train, , drop=FALSE]
    else if(is.list(hM$X))
      lapply(hM$X, function(a) a[train, , drop=FALSE])
    m <- Hmsc(Y = hM$Y[train, , drop=FALSE], X = XTrain,
              XRRR = hM$XRRR[train, , drop=FALSE],
              ncRRR = hM$ncRRR, XSelect = hM$XSelect,
              distr = hM$distr,
              studyDesign = droplevels(hM$dfPi[train,, drop=FALSE]),
              Tr = hM$Tr, C = hM$C, ranLevels = hM$rL)
    ## old code calls here setPriors, but that does nothing as its
    ## result is not saved, and it is currently skipped: CHECK
    ## THIS!
    m$YScalePar <- hM$YScalePar
    m$YScaled <- scale(m$Y, m$YScalePar[1,], m$YScalePar[2,])
    m$XInterceptInd <- hM$XInterceptInd
    m$XScalePar <- hM$XScalePar
    m$XScaled <- scale(m$X, m$XScalePar[1,], m$XScalePar[2,])
    m$TrInterceptInd <- hM$TrInterceptInd
    m$TrScalePar <- hM$TrScalePar
    m$TrScaled <- scale(m$Tr, m$TrScalePar[1,], m$TrScalePar[2,])
    m
  }
  ## parallelism would only slow-down (due to start-off time)
  hM1 <- lapply(parts, function(k) setHmsc(k, hM))
  ## STEP 2: sample Hmsc models for each nfolds * nChains case
  chains <- length(hM$postList)
  threads <- nfolds * chains
  idfold <- rep(parts, each = chains)
  seeds <- sample.int(.Machine$integer.max, threads)
  ## to be called in parallel for each chain x fold, and
  ## therefore we set nChains=1, nParallel=1 within
  
  ##Loop over all the modes to make nfolds*nchains sets of training and testing data
  for(i in 1:threads){
    set.seed(seeds[i])
    k <- idfold[i]
    message("Writing temp files for thread ", i, "/", threads)
    m <- sampleMcmc(hM1[[k]], samples = hM$samples, thin = hM$thin,
                    transient = hM$transient, initPar = hM$initPar, 
                    nChains = 1, nParallel = 1, updater = hM$updater, 
                    verbose = hM$verbose, alignPost = FALSE, 
                    engine = "HPC")
    #Call the gibs sampler on each thread, each thread is 1 chain so we can simple take the output and convert it to Rdata
    filename = file.path(TempDir,sprintf("temp_samples_%.4d_thin_%.2d_thread_%.1d.rds", hM$samples, hM$thin,i))
    saveRDS(to_json(m), file = filename)
  }
  save(idfold,threads, hM1, hM, partition, file = file.path(TempDir,sprintf("temp_fold_info_samples_%.4d_thin_%.2d.rdata",hM$samples, hM$thin)))
}
