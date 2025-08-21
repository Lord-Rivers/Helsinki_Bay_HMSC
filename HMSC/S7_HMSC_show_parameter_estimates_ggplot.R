#################################################
# Calculate the parameter estimate values
################################################

#Note that the VP results are saved. This saves a lot of time if graphs need to
#be replotted.

#If a model with a new name, different number of samples, or thinning is fitted
#new VP results will be calculated and saved.
## HOWEVER ##
#If the model have been refitted with the same number name, number of samples
#and thinning, the saved VP results must be deleted manually.

### Load packages ####
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

#### Define what outputs are required ####
calc_VP = TRUE
support.level.beta = 0.95
support.level.gamma =  0.95
support.level.omega =  0.9
#Options for part.orders
#NULL: Alphabetical
#decreasing: decreasing R2 value,
#phylo: ordered by tree tip labels (relatedness)
var.part.order.explained = "phylo" 
var.part.order.raw = "phylo"
beta_param = "Sign"
show.sp.names.beta = TRUE #Default: species names shown in beta plot if there are at most 30 species and no phylogeny
plotTree = NULL #Default: tree is plotted in Beta plot if the model includes it
gamma_param = "Sign"
omega.order = "phylo" #Default: species shown in the order they are in the model
show.sp.names.omega = TRUE #Default: species names shown in beta plot if there are at most 30 species

nChains = 4 #Change based on number of chains used

#### Set up directories ####
model_description = "Con_TSP_AT_d_dp_yr_PTCYC"
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

localDir = sprintf("./Hmsc Outputs/%s",model_description)
ModelDir = file.path(localDir, "Models")
UnfittedDir = file.path(ModelDir, "Unfitted")
ResultDir = file.path(localDir, "Results")

#Create a text output file
text.file = file.path(ResultDir,"parameter_estimates.txt")
cat(c("This file contains additional information regarding parameter estimates.\n\n"),file=text.file)

#Load in the unfitted model
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

##### Calculate parameter estimates####
if(file.exists(filename)){
  load(filename)
  m = fitted_model$posteriors
  rm(fitted_model)
  cat(sprintf("Model Description\nSamples: %.4d\tThin: %.3d\tChains: %.1d\n\n", samples, thin, nChains),file=text.file,sep="",append=TRUE)
  #Define plotting orders----
  for(x in c("var.part.order.explained","var.part.order.raw","omega.order")){
    if(is.null(get(x))){
      assign(x,list())
    }
  }
  #Extract covariate names----
  if(m$XFormula=="~."){
    covariates = colnames(m$X)[-1]
  } else {
    covariates = attr(terms(m$XFormula),"term.labels")
  }
  #Calculate VP----
  if(m$nr+length(covariates)>1 & m$ns>1 & calc_VP==TRUE){
    cli_h2("Variance partitioning")
    
    #Check if VP has already been calculated
    filename = file.path(ResultDir, sprintf("parameter_estimates_VP_%.4d_samp_%.3d_thin.csv",samples,thin))
    if(file.exists(filename)){
      cli_progress_step("Reading in Variance Partitioning Data")
      VP=list()
      temp = read.csv(filename,row.names = NULL, check.names = FALSE)
      VP$vals = temp[1:(m$nr+length(covariates)),]
      R2 = as.numeric(temp[(m$nr+length(covariates))+1,-1])
      R2_method = temp[(m$nr+length(covariates))+1,1]
      cli_progress_done()
    } else {
      #These are only needed when calculating VP
      cli_progress_step("Computing Predicted Values")
      preds = computePredictedValues(m, nParallel = 4)
      gc() #Clean up
      cli_progress_step("Evaluating Model Fit")
      MF = evaluateModelFit(hM=m, predY=preds)
      rm(preds) #Clean up
      gc() #Clean up
      
      cli_progress_step("Computing Variance Partitioning")
      VP = computeVariancePartitioning(m)
      cli_progress_done()
      vals = VP$vals
      
      cli_progress_step("Saving model fit data")
      R2 = NULL
      for(x in c("TjurR2","R2","SR2")){
        if(!is.null(MF[[x]])){
          vals = rbind(vals,MF[[x]])
          #Correctly rename the last row
          rownames(vals)[(m$nr+length(covariates))+1] = x
          R2 = MF[[x]]
          R2_method = x
        }
      }
      vals = as.data.frame(vals) |>
        tibble::rownames_to_column("Variable")
      
      # Save the results of VP
      write.csv(vals,file=filename,row.names=FALSE)
      if(!is.null(VP$R2T$Beta)){
        filename = file.path(ResultDir,sprintf("parameter_estimates_VP_R2T_Beta_%.4d_samp_%.3d_thin.csv",samples,thin))
        write.csv(VP$R2T$Beta,file=filename)
      }
      if(!is.null(VP$R2T$Y)){
        filename = file.path(ResultDir, sprintf("parameter_estimates_VP_R2T_Y_%.4d_samp_%.3d_thin.csv",samples,thin))
        write.csv(VP$R2T$Y,file=filename)
      }
      #Removing the R2 value for plotting since it is stored as the separate
      #global object R2, it must be included in the dataframe for external
      #saving
      VP$vals = vals[1:(m$nr+length(covariates)),]
      rm(vals)
    }
    
    #Allow for the order of the species plotting to be changed
    if(all(var.part.order.explained==0)){
      c.var.part.order.explained = 1:m$ns
      level_names.explained = colnames(VP$vals[-1])[c.var.part.order.explained]
    } else {
      if(all(var.part.order.explained=="decreasing")){
        c.var.part.order.explained = order(R2, decreasing = TRUE)
        level_names.explained = colnames(VP$vals[-1])[c.var.part.order.explained]
      } else {
        if(all(var.part.order.explained=="phylo") & !is.null(m$phyloTree)) {
          #This prevents errors if the tree is somehow longer then the species list
          c.var.part.order.explained  = m$phyloTree$tip.label[which((m$phyloTree$tip.label) %in% colnames(VP$vals[-1]))]
          level_names.explained = c.var.part.order.explained
        } else {
          c.var.part.order.explained  = var.part.order.explained
          level_names.explained = colnames(VP$vals[-1])[c.var.part.order.explained]
        }
      }
    }
    
    #This makes sure that the random variables are always the last variables in
    #the graph
    var_levels = c(sort(covariates), sort(sprintf("Random: %s",m$ranLevelsUsed)))
    
    #The species names are in the same order as the colnames in vp$vals, so they
    #are given levels equal to there c.var.part.order.explained
    to_plot = VP$vals %>%
      pivot_longer(!Variable, names_to = "species") %>%
      mutate(species = factor(species, levels = level_names.explained))%>%
      mutate(Variable = factor(Variable,levels = var_levels))
    VP_explained = 
      ggplot(to_plot,aes(y=value,x=species,fill=Variable)) +
      geom_bar(position = position_stack(reverse = TRUE), stat="identity",colour="black",linewidth=0.1,width=1) +
      scale_fill_manual(values = c("#66A61E", "#D95F02", "#E7298A", "#6F1253",
                                   "#1B9E77", "#7570B3", "#E6AB02", "#A6761D",
                                   "#666666"))+
      #This remove allows you to change the gap between zero on the y-axis and
      #the x-axis
      scale_y_continuous(expand=c(0,0.01))+
      xlab("") +
      ylab("Proportion of variance")+
      theme(panel.background = element_rect(fill = NA),
            panel.grid = element_blank(),
            plot.caption = element_blank(),
            legend.position = "none",
            axis.text.y = element_text(size=8, family = "sans",colour="black"),
            axis.text.x = element_text(size=4, family = "sans",colour="black",
                                       face="italic",angle = 50, hjust = 1, vjust = 1),
            text = element_text(size=8, family = "sans",colour="black"))
    ggsave(filename = file.path(ResultDir,sprintf("VP_explained_graph_%.4d_samp_%.3d_thin.pdf",samples,thin)),
           plot=VP_explained, width = 18, height = 10, unit = "cm")
    
    if(all(var.part.order.raw==0)){
      c.var.part.order.raw = 1:m$ns
      level_names.raw = colnames(VP$vals[-1])[c.var.part.order.raw]
    } else {
      if(all(var.part.order.raw=="decreasing")){
        c.var.part.order.raw = order(R2, decreasing = TRUE)
        level_names.raw = colnames(VP$vals[-1])[c.var.part.order.raw]
      } else {
        if(all(var.part.order.raw=="phylo") & !is.null(m$phyloTree)) {
          #This prevents errors if the tree is somehow longer then the species list
          c.var.part.order.raw  = m$phyloTree$tip.label[which((m$phyloTree$tip.label) %in% colnames(VP$vals[-1]))]
          level_names.raw = c.var.part.order.raw
        } else {
          c.var.part.order.raw  = var.part.order.explained
          level_names.raw = colnames(VP$vals[-1])[c.var.part.order.raw]
        }
      }
    }
    
    if(!is.null(R2)){
      VPr = VP
      for(k in 1:m$ns){
        VPr$vals[,k+1] = R2[k]*VPr$vals[,k+1]
      }
      cli_progress_step("Plotting graph to pdf")
      to_plot = VPr$vals %>%
        pivot_longer(!Variable, names_to = "species") %>%
        mutate(species = factor(species, levels = level_names.raw))%>%
        mutate(Variable = factor(Variable,levels = var_levels))
      VP_raw = 
        ggplot(to_plot,aes(y=value,x=species,fill=Variable)) +
        geom_bar(position = position_stack(reverse = TRUE), stat="identity",colour="black",linewidth=0.1,width=1)+
        scale_fill_manual(values = c("#66A61E", "#D95F02", "#E7298A", "#6F1253",
                                     "#1B9E77", "#7570B3", "#E6AB02", "#A6761D",
                                     "#666666"))+
        #This remove allows you to change the gap between zero on the y-axis and
        #the x-axis
        scale_y_continuous(expand=c(0,0.01),limits = c(0,1))+
        xlab("") +
        ylab("Proportion of variance")+
        theme(panel.background = element_rect(fill = NA),
              panel.grid = element_blank(),
              #axis.line = element_line(colour = "black"),
              legend.position = "bottom",
              plot.caption = element_blank(),
              axis.text.y = element_text(size=8, family = "sans",colour="black"),
              axis.text.x = element_text(size=4, family = "sans",colour="black",
                                         face="italic",angle = 50, hjust = 1, vjust = 1),
              text = element_text(size=8, family = "sans",colour="black"))
      ggsave(filename = file.path(ResultDir,sprintf("VP_raw_graph_%.4d_samp_%.3d_thin.pdf",samples,thin)),
             plot=VP_raw, width = 18, height = 12, unit = "cm")
      cli_process_done()
      Raw = rowMeans(VPr$vals[-c(1)])*100
    }
    cli_progress_cleanup()
    
    #Save mean variance explained by each variable to the text file
    Propotional = data.frame(VP$vals[1],rowMeans(VP$vals[-c(1)])*100)
    R2_frame = data.frame(colnames(VP$vals[-1]),R2)
    cat(c("\n","Mean porotion of explained variances"),file=text.file,sep="",append=TRUE)
    cat(sprintf("\n%-36s\t%-10s\t%-10s","Variable","Mean Prop","Mean Raw"),file=text.file,sep="",append=TRUE)
    cat(sprintf("\n%-36s\t%4f\t%4f",Propotional[[1]],Propotional[[2]],Raw),file=text.file,sep="",append=TRUE)
    #Since ggplot is used and plotting order is defined by the levels of the x
    #factor, the levels of species name can be set to c.var.part.order.explained
    #before plotting
    cat("\n","\n",R2_method," Values",file=text.file,sep="",append=TRUE)
    cat(sprintf("\n%-36s\t%4f","Mean:",mean(R2)),file=text.file,sep="",append=TRUE)
    cat(sprintf("\n\nSpecies R2 Values:\n%-36s\t%-10s","Species","R2"),file=text.file,sep="",append=TRUE)
    cat(sprintf("\n%-36s\t%4f",R2_frame[[1]],R2_frame[[2]]),file=text.file,sep="",append=TRUE)
    rm(VP)
  }
  
  # Calculate Beta ----
  if(m$nc>1){
    cli_h2("Computing Beta")
    cli_progress_step("Computing posteria estimate for Beta")
    postBeta = getPostEstimate(m, parName="Beta")
    cli_progress_step("Writing Beta values to excel")
    filename = file.path(ResultDir, "parameter_estimates_Beta.xlsx")
    me = as.data.frame(t(postBeta$mean))
    me = cbind(m$spNames,me)
    colnames(me) = c("Species",m$covNames)
    po = as.data.frame(t(postBeta$support))
    po = cbind(m$spNames,po)
    colnames(po) = c("Species",m$covNames)
    ne = as.data.frame(t(postBeta$supportNeg))
    ne = cbind(m$spNames,ne)
    colnames(ne) = c("Species",m$covNames)
    vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
    #If there is a tree present calculate and write rho estimate to the excel file
    if(!is.null(m$phyloTree)){
      mpost = convertToCodaObject(m)
      rhovals = unlist(poolMcmcChains(mpost$Rho))
      temp=as.data.frame(t(c(summary(mpost$Rho)$statistics,summary(mpost$Rho)$quantiles,mean(rhovals>0))))
      colnames(temp)[length(temp)] = "Pr[rho>0]"
      #rownames(temp) = "Rho statistic"
      rm(mpost) #Clean up
      vals = append(vals,list("Rho info" = temp))
    }
    writexl::write_xlsx(vals,path = filename)
    cli_progress_done()
    
    c.plotTree = !is.null(m$phyloTree)
    if(!is.null(plotTree)){
      c.plotTree = c.plotTree & plotTree
    }
    #Make sure c.show.sp.name is set to "FALSE" if a tree is present because and
    #plotting is requested. If no tree is present respect the input value
    #regardless of the number of species
    if(is.null(show.sp.names.beta) | isTRUE(show.sp.names.beta)){
      c.show.sp.names = as.character(!c.plotTree)
    } else {
      c.show.sp.names = as.character(show.sp.names.beta)
    }
    
    cli_progress_step("Plotting Beta graph to pdf")
    
    to_plot = switch(beta_param, 
                     Sign = sign(postBeta$mean)*((postBeta$support>support.level.beta) + (postBeta$support<(1-support.level.beta))>0),
                     Mean = postBeta$mean*((postBeta$support>support.level.beta) + (postBeta$support<(1-support.level.beta))>0),
                     Support = (2*postBeta$support-1)*((postBeta$support>support.level.beta) + (postBeta$support<(1-support.level.beta))>0)
    ) %>%
      as.data.frame() %>%
      mutate(covNames = gsub("(.*\\()(.*?)(,.*?\\))(\\d)","\\2\\4",m$covNames)) %>%
      pivot_longer(!covNames,names_to = "species") %>%
      mutate(value = switch(beta_param, 
                            Sign= factor(value),
                            Mean = value,
                            Support = factor(value)))
    
    #This the actual heatmap graph
    heatmap = ggplot(to_plot,aes(x=covNames,y=species,fill=value)) +
      geom_tile() +
      switch(beta_param, 
             Sign = scale_fill_manual(breaks = c(-1,0,1),values = c("blue","white","red"), labels = c("Negative","None","Postive"), drop=FALSE),
             scale_fill_gradient2(low="blue",mid="white", high = "red", midpoint = 0))+
      labs(x=element_blank(),y=element_blank(),fill = switch(beta_param, Sign = element_blank(), Mean = "Mean", Support = "Support"))+
      coord_fixed(1/4)+
      theme(panel.background = element_rect(fill = NA),
            panel.grid = element_blank(),
            plot.caption = element_blank(),
            axis.text.y = switch(c.show.sp.names,"TRUE" = element_text(size=4, family = "sans",colour="black"),
                                 element_blank()),
            axis.text.x = element_text(size=8, family = "sans",colour="black",
                                       angle = 30, hjust = 1, vjust = 1),
            text = element_text(size=8, family = "sans",colour="black"))
    
    if(c.plotTree){
      require(ggtree)
      require(aplot)
      #This section is for edits to the heatmap before the tree is added
      heatmap = heatmap
      #This is the tree part of the graph
      tree = ggtree(m$phyloTree,size=0.1) +
        #Because of the way ggtree works size is in mm not pt and there is no way to change it 1pt ~ 0.352mm
        geom_tiplab(offset = 0.1, size=1.4, fontface = "italic", family = "sans")+
        hexpand(1)
      #aplot does the heavy lifting of matching the tree and the row names 
      beta_plot = heatmap %>% insert_left(tree,width=0.5)
    } else {
      beta_plot = heatmap
    }
    
    #This has been roughly automated to attempt to make the graph always fit, it
    #has not been tested with additional species but 81 species works well at
    #12cm which would give a ratio of 12/81 cm/species or 6.75 species/cm.
    #Scaling is naturally required if the graph is too large, but this can also
    #be done in post.
    auto_height = m$ns * (12/81)
    auto_width = 0.2*length(m$covNames)+10.7
    ggsave(filename = file.path(ResultDir,sprintf("Beta_graph_%.4d_samp_%.3d_thin.pdf",samples,thin)),
           plot=beta_plot, width = auto_width, height = auto_height, unit = "cm")
    
    cli_progress_done()
    cli_progress_cleanup()
  }
  
  # Calculate Gama ----
  if(m$nt>1 & m$nc>1){
    cli_h2("Computing Gamma")
    cli_progress_step("Computing posteria estimate for Gamma")
    postGamma = getPostEstimate(m, parName="Gamma")
    cli_progress_step("Plotting Gamma graph to pdf")
    to_plot = switch(gamma_param, 
                     Sign = sign(postGamma$mean)*((postGamma$support>support.level.gamma) + (postGamma$support<(1-support.level.gamma))>0),
                     Mean = postGamma$mean*((postGamma$support>support.level.gamma) + (postGamma$support<(1-support.level.gamma))>0),
                     Support = (2*postGamma$support-1)*((postGamma$support>support.level.gamma) + (postGamma$support<(1-support.level.gamma))>0)
    ) %>%
      as.data.frame() %>%
      rename_with(~m$trNames,everything()) %>%
      mutate(covNames = gsub("(.*\\()(.*?)(,.*?\\))(\\d)","\\2\\4",m$covNames)) %>%
      pivot_longer(!covNames,names_to = "trait") %>%
      mutate(value = switch(gamma_param, 
                            Sign= factor(value, levels = c(-1,0,1)),
                            Mean = value,
                            Support = value))
    
    #This the actual heatmap graph
    gamma_plot = ggplot(to_plot,aes(x=covNames,y=trait,fill=value)) +
      geom_tile() +
      switch(gamma_param, 
             Sign = scale_fill_manual(breaks = c(-1,0,1),values = c("blue","white","red"), labels = c("Negative","None","Postive"), drop=FALSE),
             scale_fill_gradient2(low="blue",mid="white", high = "red", midpoint = 0))+
      labs(x=element_blank(),y=element_blank(),fill = switch(gamma_param, Sign = element_blank(), Mean = "Mean", Support = "Support"))+
      coord_fixed(1)+
      theme(panel.background = element_rect(fill = NA),
            panel.grid = element_blank(),
            plot.caption = element_blank(),
            axis.text.y = element_text(size=8, family = "sans",colour="black"),
            axis.text.x = element_text(size=8, family = "sans",colour="black",
                                       angle = 30, hjust = 1, vjust = 1),
            text = element_text(size=8, family = "sans",colour="black"))
    
    auto_height = m$nt/length(m$covNames) * 10
    auto_width = 10
    ggsave(filename = file.path(ResultDir,sprintf("Gamma_graph_%.4d_samp_%.3d_thin.pdf",samples,thin)),
           plot=gamma_plot, width = auto_width, height = auto_height, unit = "cm")
    cli_progress_done()
    cli_progress_cleanup()
  }
  
  # Calculate Omega ----
  if(m$nr>0 & m$ns>1){
    cli_h2("Computing Omega")
    cli_progress_step("Computing posteria estimate for Omega")
    OmegaCor = computeAssociations(m)
    cli_progress_done()
    for (r in 1:m$nr){
      toPlot = ((OmegaCor[[r]]$support>support.level.omega) + (OmegaCor[[r]]$support<(1-support.level.omega))>0)*sign(OmegaCor[[r]]$mean)
      if(is.null(show.sp.names.omega)){
        c.show.sp.names = (m$ns<=30) 
      } else {
        c.show.sp.names = show.sp.names.omega
      }
      if(!c.show.sp.names){
        colnames(toPlot)=rep("",m$ns)
        rownames(toPlot)=rep("",m$ns)
      }
      if(all(omega.order==0)){
        plotOrder = 1:m$ns
      } else {
        if(all(omega.order=="AOE")){
          plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
        } else if(all(omega.order=="phylo") & !is.null(m$phyloTree)) {
          #This prevents errors if the tree is somehow longer then the species list
          plotOrder  = m$phyloTree$tip.label[which((m$phyloTree$tip.label) %in% colnames(OmegaCor[[r]]$mean))]
          #This can be the names rather then numbers because a matrix can be sorted by column/row index or name
          level_names.explained = plotOrder
        } else {
          plotOrder = omega.order
        }
      }
      mymain = paste0("Associations, ",model_description, ": ",names(m$ranLevels)[[r]])
      if(m$ranLevels[[r]]$sDim>0){
        mpost = convertToCodaObject(m)
        summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))
        rm(mpost) #Clean up
      }
      cli_progress_step("Plotting Omega graph to pdf")
      #This is needed here because this graph is not a ggplot graph and thus
      #can't be saved with ggsave
      pdf(file.path(ResultDir,sprintf("Omega_%s_graph_%.4d_samp_%.3d_thin.pdf",names(m$ranLevels)[[r]],samples,thin)),width = 6.30, height = 7)
      corrplot(toPlot[plotOrder,plotOrder], method = "color",
               type = "lower",
               col=colorRampPalette(c("blue","white","red"))(3),
               mar=c(0,0,1,0), #main=mymain,
               cex.main=0.05, tl.cex=0.4, tl.col = "black")
      dev.off()
      cli_progress_step("Writing Omega values to excel")
      me = as.data.frame(OmegaCor[[r]]$mean)
      me = cbind(m$spNames,me)
      colnames(me)[1] = ""
      po = as.data.frame(OmegaCor[[r]]$support)
      po = cbind(m$spNames,po)
      colnames(po)[1] = ""
      ne = as.data.frame(1-OmegaCor[[r]]$support)
      ne = cbind(m$spNames,ne)
      colnames(ne)[1] = ""
      vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
      filename = file.path(ResultDir, paste0("parameter_estimates_Omega_",names(m$ranLevels)[[r]],".xlsx"))
      writexl::write_xlsx(vals,path = filename)
      cli_progress_done()
    }
  }
  cat("\n\nSpecies plot orders:",file=text.file,sep="",append=TRUE)
  cat(sprintf("\n%-36s\t%-36s\t%-36s","VP Explained Order","VP Raw Order","Omega Order"),file=text.file,sep="",append=TRUE)
  cat(sprintf("\n%-36s\t%-36s\t%-36s",c.var.part.order.explained,c.var.part.order.raw,plotOrder),file=text.file,sep="",append=TRUE)
  #Clean up plotting
  while(!is.null(dev.list())){
    dev.off()
  }
}
