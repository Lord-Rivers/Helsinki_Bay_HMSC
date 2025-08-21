set.seed(369)
require(Hmsc)
require(jsonify)
c.Hmsc = getS3method("c","Hmsc")

### Set up directories #### 

#The names of the arguments aren't important in the R script, only the order. The names will be defined in the batch script which will ensure that regardless of the order that they are parsed in the terminal they will be parsed to R in the order area, model_description
args = commandArgs(trailingOnly=TRUE)

area = args[1]
model_description = args[2]

print(area)
print(model_description)

if(dir.exists("/home/jdehaast")){
  setwd(sprintf("~/HMSC_Data/%s", area))
  ModelDir = file.path(sprintf("./%s/Fitted",model_description))
  TempDir = file.path(sprintf("./%s/Temp",model_description))
}else{
  cat("Can't find HPC directoy.\nExcution may be slow on personal computers\n")
  setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),"../"))
  
  localDir = sprintf("./HMSC/Hmsc Outputs/%s",model_description)
  ModelDir = file.path(localDir, "Models/Fitted")
  TempDir = file.path(localDir,"Models/Temp")
}

nChains = 4
nfolds = 2

file_list = Sys.glob(file.path(TempDir,"temp_fold_info_samples_*_thin_*.rdata"))
samples_list = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?rdata)","\\2",file_list))
thin_list = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?rdata)","\\2",file_list))
order = sort(samples_list*thin_list, index.return = TRUE)
file_list = file_list[order$ix]
rm(samples_list,thin_list,order)

filename.in = dplyr::last(file_list)
samples = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?rdata)","\\2",filename.in))
thin = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?rdata)","\\2",filename.in))

transient = ceiling(0.5*thin*samples)
filename.out = file.path(ModelDir,sprintf("MF_samples_%.4d_thin_%.2d_chains_%.1d_nfolds_%.1d.Rdata",
                                          samples, thin, nChains,nfolds))
if(!is.null(filename.in)){
  ptm = proc.time()
  message("Reading in Fold information")
  load(file = filename.in)
  postN <- Reduce(sum, lapply(hM$postList, length))
  predArray <- array(NA, c(hM$ny, hM$ns, postN))
  mods = vector("list", threads)
  for(i in 1:threads){
    temp_fitted = file.path(TempDir,sprintf("Sampled_HPC_samples_%.4d_thin_%.2d_thread_%.1d.rds", samples, thin,i))
    temp = from_json(readRDS(file = temp_fitted)[[1]])
    if(is.matrix(temp[[1]][[1]]$Alpha)){
      cat("\tAlpha is a matrix\nFixing alpha issue\n")
      #pb = txtProgressBar(min = 0, max = samples, initial = 0)
      for(z in 1:samples){
        temp_Alpha_Mat = temp[[1]][[z]]$Alpha
        temp[[1]][[z]]$Alpha = lapply(seq_len(nrow(temp_Alpha_Mat)), function(p) temp_Alpha_Mat[p,])
        #setTxtProgressBar(pb,i)
      }
    } else {
      cat("\tAlpha is not a matrix\nNo fix required\n")
    }
    k = idfold[i]
    m = importPosteriorFromHPC(hM1[[k]], temp[1], samples, thin, transient, alignPost = TRUE)
    message("finished thread ", i, "/", threads)
    attr(m, "fold") <- k
    mods[[i]] = m
  }
  rm(hM1)
  #Combine predictions: this is still a loop
  idfold <- sapply(mods, attr, which = "fold")
  parts <- sort(unique(partition))
  #This might be used in the future if I package up this code better, but for now its always set to NULL
  partition.sp = NULL
  Yc = NULL
  expected=TRUE
  for (p in parts) {
    message("predictions for partition ", p)
    val <- partition == p
    m <- do.call(c.Hmsc, mods[which(idfold == p)])
    m <- alignPosterior(m)
    postList <- poolMcmcChains(m$postList, start=1, thin=1)
    dfPi <- droplevels(hM$dfPi[val,, drop=FALSE])
    Xval <- if (is.matrix(hM$X)){
      hM$X[val, , drop=FALSE]
    } else{
      lapply(hM$X, function(a) a[val, , drop=FALSE])
    }
    pred1 <- if (is.null(partition.sp)) {
      predict(m, post=postList, X = Xval,
              XRRR = hM$XRRR[val,, drop=FALSE],
              Yc = Yc[val,, drop=FALSE], studyDesign = dfPi,
              mcmcStep = mcmcStep, expected = expected)
    } else {
      getSpeciesFoldPrediction(hM, m, val, postList, dfPi,
                               partition.sp = partition.sp,
                               mcmcStep = mcmcStep,
                               expected = expected,
                               nParallel = nParallel,
                               useSocket = useSocket)
    }
	temp = gc()
    cat(sprintf("current mermory usage Ncels: %.1f MB Vcels: %.1f MB\n", temp[3], temp[4]))
    cat("High memory use section, writing predictions to array\n")
    predArray[val,,] <- simplify2array(pred1)
	temp = gc()
    cat(sprintf("current mermory usage Ncels: %.1f MB Vcels: %.1f MB\n", temp[3], temp[4]))
    cat("Cleaning up memory\n")
    rm(pred1,val,m,postList,Xval,dfPi)
	temp = gc()
    cat("Memory clean complete\n")
    cat(sprintf("current mermory usage Ncels: %.1f MB Vcels: %.1f MB\n", temp[3], temp[4]))
  }
  preds = computePredictedValues(hM)
  cat("Calculating MF\n")
  MF = evaluateModelFit(hM, predY=preds)
  rm(pred)
  cat("Calculating MFCV\n")
  MFCV = evaluateModelFit(hM, predY=predArray)
  rm(predArray)
  cat("Calculating WAIC\n")
  WAIC = computeWAIC(hM)
  computational.time = proc.time() - ptm
  cat("Time taken:", computational.time[3],"s \n\n")
  temp = gc()
  cat(sprintf("max mermory usage\n\tNcels: %.2f MB\n\tNcels: %.2f MB\n\tTotal: %.2f MB\n", temp[11], temp[12], temp[11] + temp[12]))
  save(MF,MFCV,WAIC,file = filename.out)
  break
} else {
  message("Could not find any input files")
}

