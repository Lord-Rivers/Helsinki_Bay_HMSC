#################################################
# Here we fit the models defined for the system
################################################

#### Load packages ####
remove(list=ls())
set.seed(369)
require(Hmsc)
require(dplyr)
require(tidyr)
require(colorspace)
require(vioplot)
require(jsonify)
require(cli)
require(ggplot2)
require(gridExtra)

#### Define what outputs are required ####
showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
showAlpha = TRUE

maxOmega = 100

ma.beta = NULL
labels.beta = NULL
ma.gamma = NULL
labels.gamma = NULL
ma.omega= NULL
labels.omega=NULL

nChains = 4 #Change based on number of chains used
Lst = 1 #Do not change this

#### Set up directories ####
model_description = "Con_TSP_AT_d_dp_yr_PTCYC"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

localDir = sprintf("Hmsc Outputs/%s",model_description)
ModelDir = file.path(localDir, "Models/Fitted")
ResultDir = file.path(localDir, "Results")

#Automatic glob based method for generating the sample and thin list
if(TRUE){
  file_list = Sys.glob(file.path(ModelDir,sprintf("HPC_samples_*_thin_*_chains_%.1d.Rdata",nChains)))
  samples_list = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",file_list))
  thin_list = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?Rdata)","\\2",file_list))
  order = sort(samples_list*thin_list, index.return = TRUE)
  file_list = file_list[order$ix]
  rm(samples_list,thin_list,order)
}

#Create a text output file
text.file = file.path(ResultDir,"model_fit_details.txt")
cat(sprintf("%s\nThis file contains human readable information regarding the model fit.\nFitting Date:\t%s\n%1$s\n\n",
            strrep("-",80), strptime(Sys.time(),format = "%Y-%m-%d %H:%M")),
    file=text.file)

#### Calculating model fit variables ####
for(x in file_list){
  ptm = proc.time()
  #Extract the number of samples and thin from the file name
  samples = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",x))
  thin = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?Rdata)","\\2",x))
  load(file = x)
  posteriors = fitted_model$posteriors
  fit_times = fitted_model$fit_times
  labels = c("Names","Thin", "Samples", "Chains")
  labels2 = c("Min.   :","1st Qu.:","Median :","Mean   :","3rd Qu.:","Max.   :")
  values = c(model_description,sprintf("%.3d",thin),sprintf("%.3d",samples),sprintf("%.2d",nChains))
  rm(fitted_model)
  cli_h1("Model running")
  cli_text("Model Info:\n Samples:{samples} Thin:{thin}")
  cli_text("Current time: {format(Sys.time())}")
  cli_rule()
  cat(sprintf("%s\n%-35s\t\t\t|\t%-s",strrep("-",80),"Model Description","Fitting Times"),file=text.file,append = TRUE)
  cat(sprintf("\n%-05s\t\t%-30s\t|\tChain %.1d\t\t%.2f s", labels, values, c(1,2,3,4), unlist(fit_times)),file=text.file,append = TRUE)
  cat(sprintf("\n%35s\t\t\t|\t\t\t\t-------\n%1$35s\t\t\t|\tTotal\t\t%.2f s\n","",sum(unlist(fit_times))),file=text.file,append = TRUE)
  #Note that this code can only one one model at a time.
  nm = 1
  mpost = convertToCodaObject(posteriors, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
  nr = posteriors$nr
  cat(sprintf("%s\n%50s\n",strrep("-",80),"Parameter information"),file=text.file,append = TRUE)
  #Empty predefined dataframe for Beta and Gamma parameters
  Beta_Gamma = data.frame(Beta_PE = rep(NA,6), Beta_ES = NA, Gamma_PE=NA)
  if(showBeta){
    cli_progress_step("Calculating Beta")
    eff = effectiveSize(mpost$Beta)
    tmp_eff = summary(eff)
    psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    Beta_Gamma$Beta_PE = summary(psrf[,1])
    Beta_Gamma$Beta_ES = summary(eff)
    if(is.null(ma.beta)){
      ma.beta = psrf[,1]
      labels.beta = paste0(as.character(samples),",",as.character(thin))
    } else {
      ma.beta = cbind(ma.beta,psrf[,1])
      labels.beta = c(labels.beta,paste0(as.character(samples),",",as.character(thin)))
    }
    cli_progress_done()
  }
  if(showGamma){
    cli_progress_step("Calculating Gamma")
    psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
    Beta_Gamma$Gamma_PE = summary(psrf[,1])
    if(is.null(ma.gamma)){
      ma.gamma = psrf[,1]
      labels.gamma = paste0(as.character(samples),",",as.character(thin))
    } else {
      ma.gamma = cbind(ma.gamma,psrf[,1])
      labels.gamma = c(labels.gamma,paste0(as.character(samples),",",as.character(thin)))
    }
    cli_progress_done()
  }
  if(showBeta|showGamma){
    cat(sprintf("%s\n%-35s\t\t\t|\t%-s\n%35s\t\t\t|\n%s\t\t\t%-s\t\t\t|\t%5$s",strrep("-",80),"Beta","Gamma","","Point est:","Effective Size:"),file=text.file,append = TRUE)
    cat(sprintf("\n %s\t%.4f\t\t\t\t%.2f\t\t|\t\t%.4f",labels2,Beta_Gamma$Beta_PE,Beta_Gamma$Beta_ES,Beta_Gamma$Gamma_PE),file=text.file,append = TRUE)
  }
  if(showOmega & nr>0){
    cli_progress_step("Calculating Omega")
    labels2 = c(" Min.   :"," 1st Qu.:"," Median :"," Mean   :"," 3rd Qu.:"," Max.   :")
    line = c("","","","","","","","","")
    for(k in 1:nr){
      tmp = mpost$Omega[[k]]
      z = dim(tmp[[1]])[2]
      #Randomly sample maxOmega species pairs, since all pairs would take too long
      if(z > maxOmega){
        sel = sample(1:z, size = maxOmega)
        for(i in 1:length(tmp)){
          tmp[[i]] = tmp[[i]][,sel]
        }
      }
      psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
      tmp = summary(psrf[,1])
      #This line based text output is because there can be any number of random variables and each on is reported next to the last, so the lines grow in length each iteration
      line[1] = paste0(line[1],sprintf("Factor %i\t\t\t|%3s",k,""))
      line[2] = paste0(line[2],sprintf("%3s\t\t\t\t\t|",""))
      line[3] = paste0(line[3],sprintf("Point est:\t\t\t|%3s",""))
      for(x in 1:6){
        line[x+3] = paste0(line[x+3], sprintf("%s%3s%.2f\t|%2$3s",labels2[x],"",tmp[x]))
      }
      if(is.null(ma.omega)){
        ma.omega = psrf[,1]
        labels.omega = sprintf("%s:\n%i,%i",posteriors$ranLevelsUsed[k],samples,thin)
      } else {
        ma.omega = cbind(ma.omega,psrf[,1])
        labels.omega = c(labels.omega,sprintf("%s:\n%i,%i",posteriors$ranLevelsUsed[k],samples,thin))
      }
    }
    cli_progress_done()
    cat(sprintf("\n%s\nOmega\n",strrep("-",80)),file=text.file,append = TRUE)
    for(x in 1:9){
      cat(line[x],"\n",file=text.file,append = TRUE)
    }
  }
  
  if(showAlpha & nr>0){
    cli_progress_step("Calculating Alpha")
    #Initiate the first 3 black lines of the Alpha text output, a new line for
    #each factor is added as needed. Note that the alignment will be incorrect
    #if any Alpha has more factors then the Alpha before it, since each line
    #starts "" which does not include the tabbing characters
    line = c("","","")
    for(k in 1:nr){
      if(posteriors$ranLevels[[k]]$sDim>0){
        ess = effectiveSize(mpost$Alpha[[k]])
        psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
        line[1] = paste0(line[1],sprintf("Random factor: %10s\t\t\t|%1s",names(posteriors$ranLevels)[k],""))
        line[2] = paste0(line[2],sprintf("%25s\t\t\t|",""))
        line[3] = paste0(line[3],sprintf("%-25s\t\t\t|","Point est:"))
        for(x in 1:length(psrf[,1])){
          line[x+3] = paste0(ifelse(is.na(line[x+3]),"",line[x+3]), sprintf("%s%s%1$3s%.2f%1$3s\t\t\t|"," ",rownames(psrf)[x],psrf[x,1]))
        }
      }
    }
    cat(sprintf("%s\nalpha\n%1$s\n",strrep("-",80)),file=text.file,append=TRUE)
    for(x in 1:length(line)){
      cat(line[x],"\n",file=text.file,append = TRUE)
    }
    cli_progress_done()
  }
  computational.time = proc.time() - ptm
  cat("Time taken:", computational.time[3],"s \nCurrent time:",format(Sys.time(), "%H:%M:%S"),"\n\n")
  
  cat(sprintf("%s\n%50s\n%1$s\n",strrep("=",80),"END OF CURRENT MODEL"),file=text.file,append=TRUE)
  Lst = Lst + 1
  cli_progress_cleanup()
}

pdf=FALSE
#Check that an object created by the first step is present, prevents this
#running if no matching files are found in the first check
if(exists("Beta_Gamma")){
  #This is needed to avoid duplicate names for omega labels and it will need to
  #change based on the number of factors, this code needs to be changed
  #This is not needed if the random level name is added to the omega label 
  #labels.omegaV2 = paste0(labels.omega,letters[1:length(labels.omega)])
  switch(menu(c("Yes", "No"), title="Would you like to save the graphs as a pdf (any input to cancel plotting)?")+1,
         break,
         {cli_h1("Saving graphs as pdf objects")
           out_filename = (file= file.path(ResultDir,"MCMC_convergence_ggplot.pdf"))
           pdf=TRUE
         },
         {cli_h1("Plotting graphs in R")
           pdf=FALSE})
  par(mfrow=c(3,2))
  graphs = list()
  for(x in 1:6){
    cli_progress_step("Working with graph {x}")
    #automativally switch between beta, gamma and omega
    plot_data_temp = switch(ceiling(x/2), list(ma.beta,labels.beta,"beta"), list(ma.gamma,labels.gamma,"Gamma"), list(ma.omega,labels.omega,"Omega"))
    plot_data = plot_data_temp[[1]] %>%
      as.data.frame() %>%
      rename_with(.clos = everything(),.fn = ~plot_data_temp[[2]]) %>%
      pivot_longer(names_to="Model",cols = everything()) %>%
      #Needed so ggplot doesn't order the omega graph alphabetically by factor name
      mutate(Model = factor(Model, levels = plot_data_temp[[2]]))
    
    #Switch between the zoomed in and not zoomed in y-axis
    if(x%%2==0){
      min_y=0.9
      max_y=1.1
    } else {
      min_y = 2-round(max(plot_data_temp[[1]]),2)
      max_y = round(max(plot_data_temp[[1]]),2)
    }
    
    graphs[[x]] = ggplot(plot_data,aes(Model,value)) +
            geom_hline(yintercept = 1,linetype = "dashed", colour="grey") +
      geom_violin(scale="width",adjust = 10,trim=FALSE, fill = "pink") +
      geom_boxplot(width=0.1)+
      #Setting limits with lims or scale_*_* discards data, while setting the
      #limit of a coordinate system coord_* does not (basically zooms the graph)
      coord_cartesian(ylim=c(min_y,max_y)) +
      theme(panel.background = element_rect(fill = NA),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.caption = element_text(hjust=0),
            text = element_text(size=8, family = "sans",colour="black"))
    cli_progress_done()
  }
  plot = marrangeGrob(graphs,layout_matrix=matrix(seq_len(6),byrow = TRUE,nrow=3,ncol=2),top=NULL)
  if(pdf==TRUE){
    ggsave(out_filename,plot, height = 18, width = 18, units=c("cm"))
    dev.off()
  } else {
    plot
  }
}
