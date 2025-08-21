###################################################
# Plot variance partitioning based on species groups
##################################################

#This script requires VP to have been calculated for the models in S7, it will
#throw a warning if it can not find the saved VP files.

## Load packages ####
remove(list=ls())
require(dplyr)
require(tidyr)
require(ggplot2)
require(Hmsc)
require(colorspace)
require(corrplot)
require(writexl)
require(cli)

set.seed(369)

#### Set up directories ####
model_description = "PA_TSP_AT_d_dp_yr_PTCYC"
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

localDir = sprintf("./Hmsc Outputs/%s",model_description)
ModelDir = file.path(localDir, "Models")
UnfittedDir = file.path(ModelDir, "Unfitted")
ResultDir = file.path(localDir, "Results")

calc_VP = TRUE
nChains = 4

#### Read in files ####
load(file = file.path(UnfittedDir, "unfitted_models.RData"))

#Globed file names to always pull the highest number of samples*thin
if(TRUE){
  cli_progress_cleanup()
  cli_h1("Reading in files")
  temp = Sys.glob(file.path(ModelDir,sprintf("Fitted/HPC_samples_*_thin_*_chains_%.1d.Rdata",nChains)))
  samples_list = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",temp))
  thin_list = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?Rdata)","\\2",temp))
  order = sort(samples_list*thin_list, index.return = TRUE)
  samples_list = samples_list[order$ix]
  thin_list = thin_list[order$ix]
  samples = last(samples_list)
  thin = last(thin_list)
  filename = file.path(ModelDir,sprintf("Fitted/HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",samples,thin, nChains))
  cli_alert_success("Check good\nFile: {.file {filename}} exists")
  rm(order,samples_list,thin_list,temp)
}

#### Plot the varaince partioning ####
if(file.exists(filename)){
  load(filename)
  m = fitted_model$posteriors
  rm(fitted_model)
  
  #Extract covariate names----
  if(m$XFormula=="~."){
    covariates = colnames(m$X)[-1]
  } else {
    covariates = attr(terms(m$XFormula),"term.labels")
  }
  
  filename = file.path(ResultDir, sprintf("parameter_estimates_VP_%.4d_samp_%.3d_thin.csv",samples,thin))
  if(file.exists(filename)){
    cli_progress_step("Reading in Variance Partitioning Data")
    VP=list()
    temp = read.csv(filename,row.names = NULL, check.names = FALSE)
    VP$vals = temp[1:(m$nr+length(covariates)),]
    R2 = as.numeric(temp[(m$nr+length(covariates))+1,-1])
    cli_progress_done()
  } else {
    cli_alert_danger("Please calculate VP using script S7 before plotting grouped VP")
  }
  
  for(k in 1:m$ns){
    VP$vals[,k+1] = R2[k]*VP$vals[,k+1]
  }
  
  #This makes sure that the random variables are always the last variables in
  #the graph
  #Remove the last "Variable" since this is the R2 value
  var_levels = c(sort(covariates), sort(sprintf("Random: %s",m$ranLevelsUsed)))
  #This is the part of the script that does the grouping
  filename_inst = "../Data/VP_spp_grouping.csv"
  instructions = read.csv(filename_inst)
  #Grouping
  #Pivot long, then merge with instruction sheet and group by merged name using the desired function
  to_plot = VP$vals %>%
    pivot_longer(!Variable, names_to = "Species") %>%
    left_join(instructions,by="Species") %>%
    group_by(Group,Variable) %>%
    summarise(value = mean(value), Spp = n()) %>%
    arrange(match(Group,unique(instructions$Group))) %>%
    mutate(Group_labels = sprintf("%s (%i)",Group,Spp)) %>%
    mutate(Variable = factor(Variable,levels = var_levels))
  
  VP_raw = 
    ggplot(to_plot,aes(y=value,x=factor(Group_labels,levels=unique(Group_labels)),fill=Variable)) +
    geom_bar(position = position_stack(reverse = TRUE), stat="identity",colour="black",linewidth=0.1,width=1) +
    #be careful with this, it renames the variables 
    scale_fill_manual(values = c("#66A61E", "#D95F02", "#E7298A", "#6F1253",
                                 "#1B9E77", "#7570B3", "#E6AB02", "#A6761D",
                                 "#666666"),
                      labels=c("Day of Year", "Mean April Temp", 
                               "Total Phosphorus", "Salinity", 
                               "Sea Surface Temperature", "Station Depth", 
                               "Year", "Random: Location","Random: Year"))+
    #This remove allows you to change the gap between zero on the y-axis and
    #the x-axis
    scale_y_continuous(expand=c(0,0.01), limits = c(0,1))+
    xlab("") +
    ylab("Proportion of variance")+
    theme(panel.background = element_rect(fill = NA),
          panel.grid = element_blank(),
          #axis.line = element_line(colour = "black"),
          plot.caption = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(size=8, family = "sans",colour="black"),
          axis.text.x = element_text(size=8, family = "sans",colour="black",angle = 45, hjust = 1, vjust = 1),
          text = element_text(size=8, family = "sans",colour="black"))
  ggsave(filename = file.path(ResultDir,sprintf("VP_raw_grouped_graph_%.4d_samp_%.3d_thin.pdf",samples,thin)),
         plot=VP_raw, width = 16.5, height = 10, unit = "cm")
}
