################################################################################
## Plot predictive graphs for PA, CON, and/or the hurdle model  ##
################################################################################

#This script assumes that a hurdle model is being created; however, it has
#sections for producing graphs related to the PA and CON components of the
#hurdle model in addition to plotting hurdle model predictions

#NB This script should be ran section by section, as some graphs may not be
#required.

#In the MS only the outputs from the Richness section of this script were used
#for the variable Temp

## Load packages ####
remove(list=ls())
require(vegan)
require(abind)
require(cli)
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
  thin = as.numeric(gsub("(.*?_samples_)([0-9]*)(.*?Rdata)","\\2",filename[types[1]]))
  rm(file_list,samples_list,thin_list,order)
}

#This checks if Raw_PRedsCon has been loaded and if it has verifies that it is
#the correct length. The expected length is samples*chains since the predictions
#are the merged results of each sample from each chain. WARNING: This only
#verifies the file based on length; thus, it WILL NOT reload the file if you
#have made new predictions for the same sampling length and number of chains.
#This saves time when plotting multiple different types of graphs by only
#reading in the predictions if they have not already been read in. Note that
#this script unlike the others does not clear the Environment automatically in
#order to keep this file loaded even when running all lines.
if(ifelse(exists("Raw_PredsCon"),length(Raw_PredsCon[[1]][[1]])/samples!=nChains, TRUE)){
  cli_progress_step("Importing Conditional predictions")
  load(filename["Con"])
  Raw_PredsCon = Preds
  rm(Preds)
  cli_progress_step("Importing precence-absence predictions")
  load(filename["PA"])
  Raw_PredsPA = Preds
  rm(Preds)
  cli_progress_step("Importing fitted model for raw observations")
  load(sprintf("./Hmsc Outputs/Con_%s/Models/Unfitted/unfitted_models.RData",model_code))
  Raw_ObsAB = models[[1]]$Y
  Raw_Variables_AB = models[[1]]$XData
  Raw_study_design = models[[1]]$studyDesign
  cli_progress_step("Importing fitted model for raw precences")
  load(sprintf("./Hmsc Outputs/PA_%s/Models/Unfitted/unfitted_models.RData",model_code))
  Raw_ObsPA = models[[1]]$Y
  Raw_Variables_PA = models[[1]]$XData
  cli_progress_done()
  rm(models)
  check_flag = TRUE
  gc()
}

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

Raw_ObsAB_BT = switch (Trans_Used,
                       log = exp(Raw_ObsAB),
                       None = Raw_ObsAB
                       )

Obs_AB = data.frame(Count = rowSums(Raw_ObsAB_BT,na.rm = TRUE),
                    Variable = Raw_Variables_AB[[Traget_variable]])
Obs_PA = data.frame(Count = rowSums(Raw_ObsPA,na.rm = TRUE),
                    Variable = Raw_Variables_PA[[Traget_variable]])
#Back transform the PredsCon if a transformation was used, only supports none and log
PredsCon_BT = switch (Trans_Used,
                      log = lapply(PredsCon,exp),
                      None = PredsCon
                      )
#This line produces the hurdle model
PredsAB = Map("*",PredsCon_BT, PredsPA)

q = c(0.025,0.5,0.975)
n = length(Raw_PredsCon[[1]]$Gradient$XDataNew[[1]])
#Clean up memory, useful if the model objects are large
gc()

#Create the hurdle mode directory if it does not exist
hurdle.dir=sprintf("./Hmsc Outputs/Hurdle_%s",model_code)
if(!dir.exists(hurdle.dir)){
  dir.create(hurdle.dir)
  dir.create(file.path(hurdle.dir,"Graphs"))
}

###Conditional Abundance Graphs ####
predS = PredsCon_BT |>
  lapply(rowSums) |>
  abind(along=2)
qpred_temp = apply(predS, c(1), quantile,
                   probs = q, na.rm=TRUE)

condition_ab = data.frame(lc=qpred_temp[1,],
                          mu=qpred_temp[2,],
                          uc = qpred_temp[3,],
                          Variable = Gradient)

filename = file.path(hurdle.dir,sprintf("Graphs/condition_ab_%s_%s.pdf",Traget_variable,effect))
if(TRUE){
  check_flag = TRUE
  if(file.exists(filename)){
    cli_alert_warning("The graph you are trying to plot already exsits would you like to replot it")
    switch(menu(c("Yes", "No"), title="Do you want to replot?")+1,
           break, 
           {check_flag = TRUE 
           cli_alert_success("Replotting")} ,
           {check_flag = FALSE 
           cli_alert_success("Not Replotting")})
  }
  if(check_flag == TRUE){
    graph = ggplot(condition_ab,aes(Variable,mu)) +
      geom_point(data = Obs_AB ,aes(Variable,Count), colour = "grey")+
      geom_point() +
      geom_line() +
      geom_ribbon(aes(ymin = lc, ymax = uc),alpha=0.5,fill="blue")+
      labs(y="Conitional aboundance\n[A|(PA=1)]", x = Traget_variable, 
           caption = sprintf("The relationship between conditional abundance and %s", 
                             Traget_variable)) +
      coord_cartesian(ylim=c(0,min(max(condition_ab$mu),max(Obs_AB))))+
      theme(panel.background = element_rect(fill = NA), 
            panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
            axis.line = element_line(colour = "black"),
            plot.caption = element_text(hjust=0),
            axis.title = element_text(size=8, family = "sans"),
            text=element_text(size=8, family = "sans"))
    ggsave(file=filename,plot=graph, width = 10.5, height = 5.5, unit = "cm")
  }
}

###Abundance Hurdled Graph #####
predS = PredsAB |>
  lapply(rowSums) |>
  abind(along=2)
qpred_temp = apply(predS, c(1), quantile,
                   probs = q, na.rm=TRUE)

hurdle_ab = data.frame(lc=qpred_temp[1,],
                       mu=qpred_temp[2,],
                       uc = qpred_temp[3,],
                       Variable = Gradient)

filename = file.path(hurdle.dir,sprintf("Graphs/hurdle_ab_%s_%s.pdf",Traget_variable,effect))
if(TRUE){
  check_flag = TRUE
  if(file.exists(filename)){
    cli_alert_warning("The graph you are trying to plot already exsits would you like to replot it")
    switch(menu(c("Yes", "No"), title="Do you want to replot?")+1,
           break, 
           {check_flag = TRUE 
           cli_alert_success("Replotting")} ,
           {check_flag = FALSE 
           cli_alert_success("Not Replotting")})
  }
  if(check_flag == TRUE){
    graph = ggplot(hurdle_ab,aes(Variable,mu)) +
      geom_point(data = Obs_AB ,aes(Variable,Count), colour = "grey")+
      geom_point() +
      geom_line() +
      geom_ribbon(aes(ymin = lc, ymax = uc),alpha=0.5,fill="blue")+
      labs(y="Hurdled abundance", x = Traget_variable, 
           caption = sprintf("The relationship between Hurdled abundance and %s", 
                             Traget_variable)) +
      coord_cartesian(ylim=c(0,min(max(hurdle_ab$mu),max(Obs_AB))))+
      theme(panel.background = element_rect(fill = NA), 
            panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
            axis.line = element_line(colour = "black"),
            plot.caption = element_text(hjust=0),
            axis.title = element_text(size=8, family = "sans"),
            text=element_text(size=8, family = "sans"))
    ggsave(file=filename,plot=graph, width = 10.5, height = 5.5, unit = "cm")
  }
}

###Richness Graph ####
#This calculates the species richness for each sample across each chain
predS = abind(lapply(PredsPA, rowSums),along=2)
qpred_temp = apply(predS, c(1), quantile,
                   probs = q, na.rm=TRUE)

richness = data.frame(lc=qpred_temp[1,],
                      mu=qpred_temp[2,],
                      uc = qpred_temp[3,],
                      Variable = Gradient)

filename = file.path(hurdle.dir,sprintf("Graphs/Richness_%s_%s.pdf",Traget_variable,effect))
if(TRUE){
  check_flag = TRUE
  if(file.exists(filename)){
    cli_alert_warning("The graph you are trying to plot already exsits would you like to replot it")
    switch(menu(c("Yes", "No"), title="Do you want to replot?")+1,
           break, 
           {check_flag = TRUE 
           cli_alert_success("Replotting")} ,
           {check_flag = FALSE 
           cli_alert_success("Not Replotting")})
  }
  if(check_flag == TRUE){
    graph = 
      ggplot(richness,aes(Variable,mu)) +
      geom_point(data = Obs_PA ,aes(Variable,Count), colour = "grey")+
      geom_point() +
      geom_line() +
      ylim(c(0,max(richness$uc,Obs_PA$Count)))+
      geom_ribbon(aes(ymin = lc, ymax = uc),alpha=0.5,fill="blue")+
      labs(y="Richness", x = Traget_variable)+
      theme(panel.background = element_rect(fill = NA), 
            panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
            axis.line = element_line(colour = "black"),
            plot.caption = element_text(hjust=0),
            axis.title = element_text(size=8, family = "sans"),
            text=element_text(size=8, family = "sans"))
    ggsave(file=filename,plot=graph, width = 16.5, height = 5.5, unit = "cm") #Default size is 10.5 cm textwidth is around 16.5cm
  }
}

##The section below was not used in the MS at all ##
###Graphs for a specific species ####
spp_name =  "Anabaenopsis elenkinii"
spps = colnames(PredsAB[[1]])
species_index = which(spps==spp_name)
#Probability of occurrence
qpred_pa = abind(PredsPA,along=3)|>
  apply(c(1, 2), quantile, probs = q, na.rm=TRUE)
Occurrence = data.frame(lc = qpred_pa[1, ,species_index],
                       mu = qpred_pa[2, ,species_index],
                       uc = qpred_pa[3, ,species_index],
                       Variable = Gradient)

#Conditional abundance
tmp = abind(PredsCon_BT,along=3)
qpred_com = apply(tmp,c(1, 2), quantile, probs = q, na.rm=TRUE)
qpred = qpred_com[, , species_index]
conditional = data.frame(lc = qpred_com[1, ,species_index],
                         mu = qpred_com[2, ,species_index],
                         uc = qpred_com[3, ,species_index],
                         Variable = Gradient)
#Combined
qpred = abind(PredsAB,along=3) |>
  apply(c(1, 2), quantile, probs = q, na.rm=TRUE)
abundance = data.frame(lc = qpred[1, ,species_index],
                       mu = qpred[2, ,species_index],
                       uc = qpred[3, ,species_index],
                       Variable = Gradient)
obs = switch (Trans_Used,
              log = exp(Raw_ObsAB[,species_index]),
              None = Raw_ObsAB[,species_index]
)
observations = data.frame(obs = obs,
                          Variable = Raw_Variables_AB[[Traget_variable]])
observations$pa = ifelse(is.na(observations$obs),0,1)

#Separate graphs
#To output species preictions as seperate graphs change this to true
seperate_graphs = FALSE
if(seperate_graphs){
  
  ggplot(abundance,aes(Variable,mu)) +
    geom_point(data = observations ,aes(y=obs,x=Variable), colour = "grey", na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    geom_ribbon(aes(ymin = lc, ymax = uc),alpha=0.5,fill="blue")+
    coord_cartesian(ylim=c(0,min(max(c(abundance$mu,observations$obs),na.rm = TRUE),quantile(observations$obs,na.rm=TRUE,0.95)*2)))+
    labs(y="Abundance", x = Traget_variable) +
    theme(panel.background = element_rect(fill = NA), 
          panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
          axis.line = element_line(colour = "black"),
          plot.caption = element_text(hjust=0),
          axis.title = element_text(size=8, family = "sans"),
          text=element_text(size=8, family = "sans"))
  ggsave(file.path(hurdle.dir,sprintf("Graphs/abundance_%s_%s_%s.pdf",gsub(" ","_",spp_name),Traget_variable,effect)), width = 10.5, height = 5.5, unit = "cm")
  
  ggplot(conditional,aes(Variable,mu)) +
    geom_point(data = observations ,aes(y=obs,x=Variable), colour = "grey", na.rm = TRUE)+
    geom_point(na.rm = TRUE)+
    geom_ribbon(aes(ymin = lc, ymax = uc),alpha=0.5,fill="blue")+
    #scale_y_continuous(breaks=seq(-6,10,by=2))+
    coord_cartesian(ylim=c(0,min(max(c(conditional$mu,observations$obs),na.rm = TRUE),quantile(observations$obs,na.rm=TRUE,0.95)*2)))+
    #coord_trans(y="log") +
    labs(y="Conditional abundance (A|PA=1)", x = Traget_variable) +
    theme(panel.background = element_rect(fill = NA), 
          panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
          axis.line = element_line(colour = "black"),
          plot.caption = element_text(hjust=0),
          axis.title = element_text(size=8, family = "sans"),
          text=element_text(size=8, family = "sans"))
  ggsave(file.path(hurdle.dir,sprintf("Graphs/conditional_%s_%s_%s.pdf",gsub(" ","_",spp_name),Traget_variable,effect)), width = 10.5, height = 5.5, unit = "cm")
}

ggplot(Occurrence,aes(Variable,mu)) +
  geom_point(data = observations ,aes(y=pa, x=Variable),colour="grey") +
  geom_ribbon(aes(ymin = lc, ymax = uc),alpha=0.5,fill="blue")+
  geom_point() +
  labs(y="Occurrence", x = Traget_variable) +
  theme(panel.background = element_rect(fill = NA), 
        panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
        axis.line = element_line(colour = "black"),
        plot.caption = element_text(hjust=0),
        axis.title = element_text(size=8, family = "sans"),
        text=element_text(size=8, family = "sans"))
ggsave(file.path(hurdle.dir,sprintf("Graphs/Occurrence_%s_%s_%s.pdf",gsub(" ","_",spp_name),Traget_variable,effect)), width = 10.5, height = 5.5, unit = "cm")

#Combined graph
filename = file.path(hurdle.dir,sprintf("Graphs/combined_%s_%s_%s.pdf",gsub(" ","_",spp_name),Traget_variable,effect))
if(TRUE){
  check_flag == TRUE
  if(file.exists(filename)){
    cli_alert_warning("The graph you are trying to plot already exsits would you like to replot it")
    switch(menu(c("Yes", "No"), title="Do you want to replot?")+1,
           break, 
           {check_flag = TRUE 
           cli_alert_success("Replotting")} ,
           {check_flag = FALSE 
           cli_alert_success("Not Replotting")})
  }
  if(check_flag == TRUE){
    graph = ggplot(abundance,aes(Variable,mu)) +
      geom_point(data = observations ,aes(y=obs,x=Variable), colour = "grey", alpha = 0.5,na.rm=TRUE)+
      geom_point(data = conditional,aes(Variable,mu), colour="blue") +
      #This allows the probability of occurrence to be scaled by the maximum value
      #of conditional abundance. This must be changed in the future to align the
      #axes.
      #geom_point(data = occurrence,aes(Variable,mu*max(conditional$mu,abundance$uc)), colour="red")+
      #This is here to ensure that the hurdled abundance is plotted on top of all other points
      geom_point() +
      geom_ribbon(data = abundance, aes(ymin = lc, ymax = uc),alpha=0.5,fill="black")+
      #scale_y_continuous(limits = c(0,max(conditional$mu,abundance$uc)), sec.axis = sec_axis(~./max(conditional$mu,abundance$uc), name="Probablity of Occurrence"))+
      #If you want to include the observations use this instead
      geom_point(data = Occurrence,aes(Variable,mu*max(abundance$uc,observations$obs,na.rm=TRUE)), colour="red")+
      scale_y_continuous(limits = c(0,max(abundance$uc,observations$obs,na.rm=TRUE)), sec.axis = sec_axis(~./max(abundance$uc,observations$obs,na.rm=TRUE), name="Probablity of Occurrence"))+
      labs(y="Abundance", x = Traget_variable, caption = 
             sprintf("The Probablity of Occurrence (Red), Conditional Abundance (Blue) and the predicted Abundance of %s in response to %s",
                     spp_name, Traget_variable)) +
      theme(panel.background = element_rect(fill = NA), 
            panel.grid = element_line(colour="grey", linewidth = 0.25, linetype = 3),
            axis.line = element_line(colour = "black"),
            plot.caption = element_text(hjust=0),
            axis.title = element_text(size=8, family = "sans"),
            text=element_text(size=8, family = "sans"))
    ggsave(file=filename,plot=graph, width = 10.5, height = 5.5, unit = "cm")
  }
}

