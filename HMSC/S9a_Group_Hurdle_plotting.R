################################################################################
## Plot predictive graphs for the hurdle model based on species groups  ##
################################################################################

## Load packages ####
remove(list=ls())
require(vegan)
require(abind)
require(cli)
require(dplyr)
require(ggplot2)


#### Define inputs ####
#The name of the model without the PA or Con part because both will need to be
#imported to create the Hurdle model
model_code = "TSP_AT_d_dp_yr_PTCYC"
#Define if a transformation was used for the abundance data. 
#Currently the only option is "None" or "log"
Trans_Used = "log"
nChains = 4 #Make sure this is correct

#### Set up directories ####
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))

### Read in data ####

#Note this method will load the longest model run for which both PA and CON
#models are present.
if(TRUE){
  types = c("Con","PA")
  file_list = NULL
  filename = NULL
  order = NULL
  for(x in 1:2){
    ModelDir = shortPathName(file.path(sprintf("./Hmsc Outputs/%s_%s/Results/Preds",types[x],model_code)))
    file_list[[types[x]]] = Sys.glob(file.path(ModelDir,sprintf("/Preds_%s_%s_HPC_samples_*_thin_*_chains_%s.Rdata",types[x],model_code,nChains)))
    samples_list = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",file_list[[types[x]]]))
    thin_list = as.numeric(gsub("(.*?_thin_)([0-9]*)(.*?Rdata)","\\2",file_list[[types[x]]]))
    order[[types[x]]] = sort(samples_list*thin_list, decreasing = TRUE, index.return = TRUE)
  }
  
  matched_lengths = intersect(order[[1]]$x,order[[2]]$x)
  filename[types[1]] = file_list[[types[1]]][order[[types[1]]]$ix[which(order[[types[1]]]$x==max(matched_lengths))]]
  filename[types[2]] = file_list[[types[2]]][order[[types[2]]]$ix[which(order[[types[2]]]$x==max(matched_lengths))]]
  samples = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",filename[types[1]]))
  rm(file_list,samples_list,thin_list,order)
}

#This is probably not the best way of doing this, check if Raw_PRedsCon has been
#loaded if it has verify that it is the right length. Length should be
#samples*chains since the predictions are the merged results of each sample from
#each chain, this only verifies that th file is the expected length and will no
#reload if you have made predictions for the same sampling length and number of
#chains. This manually relies on there being 4 chains, if more or less are used
#change the 4 to number of chains.
if(ifelse(exists("Raw_PredsCon"),length(Raw_PredsCon[[1]][[1]])/samples!=4, TRUE)){
  cli_progress_step("Importing Conditional predictions")
  load(filename["Con"])
  Raw_PredsCon = Preds
  rm(Preds)
  cli_progress_step("Importing precence-absence predictions")
  load(filename["PA"])
  Raw_PredsPA = Preds
  rm(Preds)
  cli_progress_step("Importing fitted model for raw observations")
  load(sprintf("./Hmsc Outputs/Con_%s/Models/Unfitted/unfitted_models.Rdata",model_code))
  Raw_ObsAB = models[[1]]$Y
  Raw_Variables_AB = models[[1]]$XData
  cli_progress_step("Importing fitted model for raw precences")
  load(sprintf("./Hmsc Outputs/PA_%s/Models/Unfitted/unfitted_models.Rdata",model_code))
  Raw_ObsPA = models[[1]]$Y
  Raw_Variables_PA = models[[1]]$XData
  cli_progress_done()
  rm(models)
  check_flag = TRUE
  gc()
}

#Read in species groupings
spp_groups = read.csv("../Data/VP_spp_grouping.csv")

#### Extract predictions for target variable ####
# This input is here because the code above this line does not need to be reran
# to make predictions for different variables from the same models
Traget_variable = "Temp"
#Tot for total effect, Mar for marginal. Marginal holds all other variables at
#there median values (jday and year are set as per the prediction script)
effect = "Mar" 

PredsCon = Raw_PredsCon[[Traget_variable]][[
  switch(effect,Mar="predY2",Tot="predY",cli_alert_danger("Select Mar or Tot for effect"))
]]
PredsPA = Raw_PredsPA[[Traget_variable]][[
  switch(effect,Mar="predY2",Tot="predY",cli_alert_danger("Select Mar or Tot for effect"))
]]
Gradient = Raw_PredsPA[[Traget_variable]][[
  switch(effect,Mar="Gradient2",Tot="Gradient",cli_alert_danger("Select Mar or Tot for effect"))
]][["XDataNew"]][[Traget_variable]] 

#Back transform the PredsCon if a transformation was used, only supports none and log
PredsCon_BT = switch (Trans_Used,
                      log = lapply(PredsCon,exp),
                      None = PredsCon
)
PredsAB = Map("*",PredsCon_BT, PredsPA)

q = c(0.025,0.5,0.975)
n = length(Raw_PredsCon[[1]]$Gradient$XDataNew[[1]])

#Create the hurdle mode directory if it does not exist
hurdle.dir=sprintf("./Hmsc Outputs/Hurdle_%s",model_code)
if(!dir.exists(hurdle.dir)){
  dir.create(hurdle.dir)
  dir.create(file.path(hurdle.dir,"Graphs"))
}

#Clean up unneeded objects to reduce memory load
rm(PredsCon,Raw_ObsAB, Raw_Variables_AB, Raw_ObsPA,Raw_Variables_PA,Raw_PredsCon)
gc()


###Proportional richness----
propotion_rich = NULL
predS = abind(lapply(PredsPA, rowSums),along=2)
for(y in unique(spp_groups$Group)){
  print(y)
  spps = spp_groups[which(spp_groups$Group==y),"Species"]
  #This very annoying ifstatment is because the ifstament does not work in the
  #piped section, so it has to be coded like this.
  if(length(spps)>1){
    Pred_Prop_tem = PredsPA %>%
      lapply(function(x) x[,spps]) %>%
      lapply(rowSums) %>%
      abind(along=2) %>%
      `/`(predS) %>%
      apply(c(1), quantile,probs = q, na.rm=TRUE) %>%
      t() %>%
      as.data.frame() %>%
      cbind(Gradient) %>%
      rename("lc" = "2.5%", "mu" = "50%", "uc" = "97.5%", "Variable" = "Gradient") %>%
      mutate(Group=y)
  } else {
    Pred_Prop_tem = PredsPA %>%
      lapply(function(x) x[,spps]) %>%
      abind(along=2) %>%
      `/`(predS) %>%
      apply(c(1), quantile,probs = q, na.rm=TRUE) %>%
      t() %>%
      as.data.frame() %>%
      cbind(Gradient) %>%
      rename("lc" = "2.5%", "mu" = "50%", "uc" = "97.5%", "Variable" = "Gradient") %>%
      mutate(Group=y)
  }
  propotion_rich = rbind(propotion_rich,Pred_Prop_tem)
}

graph_pa = 
  ggplot(propotion_rich,aes(Variable,mu, colour=Group)) +
  geom_point() +
  scale_colour_manual(values=c("#69ef7b","#c5089e","#66d9cf","#154e56","#c0e087",
                               "#6146ca","#c7c2ec","#732a66","#799d10","#bf83f8",
                               "#096013"))+
  ylim(c(0,1))+
  labs(y="Propotion", x = Traget_variable)+
  theme(panel.background = element_rect(fill = NA), 
        panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        plot.caption = element_text(hjust=0),
        axis.title = element_text(size=10, family = "sans"),
        text=element_text(size=10, family = "sans"))


###Proportional abundance----
propotion_ab = NULL
pred_AB = abind(lapply(PredsAB, rowSums),along=2)
for(y in unique(spp_groups$Group)){
  print(y)
  spps = spp_groups[which(spp_groups$Group==y),"Species"]
  #This very annoying ifstatment is because the ifstament does not work in the
  #piped section, so it has to be coded like this.
  if(length(spps)>1){
    Pred_Prop_tem = PredsAB %>%
      lapply(function(x) x[,spps]) %>%
      lapply(rowSums) %>%
      abind(along=2) %>%
      `/`(pred_AB) %>%
      apply(c(1), quantile,probs = q, na.rm=TRUE) %>%
      t() %>%
      as.data.frame() %>%
      cbind(Gradient) %>%
      rename("lc" = "2.5%", "mu" = "50%", "uc" = "97.5%", "Variable" = "Gradient") %>%
      mutate(Group=y)
  } else {
    Pred_Prop_tem = PredsAB %>%
      lapply(function(x) x[,spps]) %>%
      abind(along=2) %>%
      `/`(pred_AB) %>%
      apply(c(1), quantile,probs = q, na.rm=TRUE) %>%
      t() %>%
      as.data.frame() %>%
      cbind(Gradient) %>%
      rename("lc" = "2.5%", "mu" = "50%", "uc" = "97.5%", "Variable" = "Gradient") %>%
      mutate(Group=y)
  }
  propotion_ab = rbind(propotion_ab,Pred_Prop_tem)
}

graph_ab = 
  ggplot(propotion_ab,aes(Variable,mu, colour=Group)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values=c("#69ef7b","#c5089e","#66d9cf","#154e56","#c0e087",
                               "#6146ca","#c7c2ec","#732a66","#799d10","#bf83f8",
                               "#096013"))+
  ylim(c(0,1))+
  labs(y="Propotion", x = Traget_variable)+
  theme(panel.background = element_rect(fill = NA), 
        legend.position="bottom",
        panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
        axis.line = element_line(colour = "black"),
        plot.caption = element_text(hjust=0),
        axis.title = element_text(size=10, family = "sans"),
        text=element_text(size=10, family = "sans"))

### Combind graph output ####
require(ggpubr)
Combind_graph = 
  ggarrange(graph_pa, graph_ab, labels = c("A.","B."), hjust = 0, vjust = 1.5,
            font.label = list(size = 14, color = "black", face = "bold",
                              family = "serif"),
            align = c("v"), ncol=1, common.legend=TRUE, legend="right")
filename = file.path(hurdle.dir,sprintf("Graphs/Group_Combind_graph_Prop_%s_%s.pdf",Traget_variable,effect))
ggsave(file=filename,plot=Combind_graph, width = 16.5, height = 12.5, unit = "cm")
