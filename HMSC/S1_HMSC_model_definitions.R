################################################
#### Hmsc analyses on Helsinki Data ####
################################################

####Clean up ####
#This cleans everything up, including making sure MASS is not loaded
remove(list=ls())
pkg = names(sessionInfo()$otherPkgs)
if(!is.null(pkg)){
  lapply(paste0("package:",pkg), detach, character.only = TRUE)
}
rm(pkg)
gc()

require(dplyr)
require(tidyr)
require(Hmsc)
require(ape)

#### Set up directories ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

set.seed(369)
localDir = "./Hmsc Outputs"
dataDir = "../Data/Formatted/"
treeDir = "../Data/Raw/"
ModelDir = file.path(localDir, "Models")
MixingDir = file.path(localDir, "Mixing")
MFDir = file.path(localDir, "Model_fit")

#### Read in Data ####
sppMat = read.csv(file.path(dataDir,'rawObs.csv'),check.names=FALSE) %>%
  arrange(sampleID)
envMat = read.csv(file.path(dataDir,'rawEnv.csv'))%>%
  arrange(sampleID)
spring_data = read.csv(file.path(treeDir,'spring_perc_temp.csv'))
taxonomy = read.table(file.path(treeDir,'HelSummerTrends_bscl.txt'), header = TRUE, sep = '\t')
spp.tree = read.tree(file.path(treeDir,'HelSummerTrends_bstree.txt'))

#Clean up the trait import and add mixotrophy as a trophic level, set auto- and
#heterotrophy to 0 for mixotrophs
traits = read.csv(file.path(dataDir,'traits.csv'))|>
  group_by(Species) %>%
  select(Species,nfix,silica,motility,autotrophy,heterotrophy) |>
  summarise(across(everything(),~max(.x,na.rm = TRUE))) |>
  mutate(across(everything(),~if_else(.x==-Inf,NA,.x))) |>
  mutate(mixotrophy = if_else(autotrophy+heterotrophy==2,1,0)) |>
  mutate(autotrophy = if_else(mixotrophy==1,0,autotrophy),
         heterotrophy = if_else(mixotrophy==1,0,heterotrophy)) |>
  tibble::column_to_rownames("Species")

#### Data sanitation #### 

envMat = envMat %>%
  #Remove year 2018, since only 3 off shore sites were sampled.
  #Remove sample 3439, this was a temperature outlier. The temperature value
  #(1.7 c) was regarded as being incorrectly recorded as all temperatures from
  #the same day are between 17.8 and 21 C
  filter(year<2018,sampleID != 3439) %>%
  mutate(jday = lubridate::yday(obstime)) %>%
  left_join(spring_data, by=join_by("year"=="Year"))

#Make sure that the sppMat only contains samples for which there are envormental
#records
sppMat = sppMat %>%
  filter(sampleID %in% envMat$sampleID)
  
envMat$stn = as.factor(as.character(envMat$stn))

#Remove uneeded varaibles
rm(dataDir, treeDir)

#Remove non species level observations, but keep all subspecies, forma and
#varieties
spps = taxonomy |>
  filter(!is.na(Species))|>
  filter(valid_name %in% colnames(sppMat)) |>
  select(valid_name) |>
  arrange(valid_name) |>
  unlist(use.names=FALSE)

rm(taxonomy)

sppMat = sppMat %>%
  select(sampleID, all_of(spps)) %>%
  tibble::column_to_rownames("sampleID")

#### Define study Design ####
studyDesign = data.frame(sampleID = as.factor(envMat$sampleID),
                         stn = envMat$stn,
                         Location = as.factor(envMat$Location),
                         year = as.factor(envMat$year),
                         obstime = as.factor(envMat$obstime))

#xy columns must be ordered so that column 1 is lon and column 2 is lat when
#using latitude and longitude
xy = envMat |> 
  select(stn,lon,lat) |>
  distinct(lon,lat, .keep_all = TRUE) |>
  tibble::column_to_rownames("stn")

Year =  envMat %>%
  distinct(year) %>%
  arrange(year) %>%
  mutate(ID = year) %>%
  tibble::column_to_rownames("ID")

#### Clean up the envovromental matrix ####
#The strange which statement makes this dynamic while allowing for SampleID to
#be removed as column and used as the row names. Dplry throws an error if one of
#the columns you are trying to select or drop is not present this allows the
#code to automatically drop all the columns used in the studyDesign matrix from
#the XData matrix.
envMat = envMat %>%
  tibble::column_to_rownames("sampleID") %>%
  dplyr::select(!names(studyDesign)[which(names(studyDesign)!="sampleID")],
                -c(lat,lon,obstime,decimal_date),year)

#Check input matrix/dataframe structures
head(studyDesign)
head(envMat)
head(xy)
head(Year)

#### Preparing the species matrices ####
pa_matrix = sppMat>0
pa_matrix = apply(pa_matrix,MARGIN = 2,FUN = as.numeric)

#Remove rare species
rare_cutoff = ceiling(dim(pa_matrix)[1]*0.05)
pa_matrix = pa_matrix[,colSums(pa_matrix)>rare_cutoff]
#Make sure no species is seen at every sample
pa_matrix = pa_matrix[,colSums(pa_matrix)<dim(envMat)[1]]

sppMat = sppMat[,colnames(pa_matrix)]
sppMat = as.matrix(sppMat)

#Ensure that the trait matrix only contains the species present in the species
#matrix
traits = traits |>
  filter(rownames(traits) %in% colnames(pa_matrix))

###Phylogeny/Taxonomic tree ####
load("../Data/New_Taxanomy.Rdata")

#Clean up the taxanomy data, namely fix the title of some plantae levels and
#manually rename flos-aquae to flosaquae
new_taxanomy = Taxanomy |>
  mutate(Phylum = if_else(is.na(Phylum), `Phylum (Division)`, Phylum, Kingdom),
         Subphylum  = if_else(is.na(Subphylum), `Subphylum (Subdivision)`, Subphylum),
         valid_name = case_when(!is.na(Variety) ~Variety,
           !is.na(Forma) ~Forma,
           #!is.na(Subspecies) ~Subspecies, 
           .default = Species)) |>
  #Manually fix the flos-aquae vs flosaquae issue
  mutate(valid_name = if_else(valid_name=="Aphanizomenon flos-aquae","Aphanizomenon flosaquae",valid_name)) |>
  filter(valid_name %in% colnames(pa_matrix))

to_tree = new_taxanomy |>
  filter(valid_name %in% colnames(pa_matrix)) |>
  select(Kingdom,Phylum,Class,Order,Family,Genus,valid_name) |>
  #Manually fixing missing higher levels, adds hupothetical families to Genera
  #that are currenlty not placed in a family since the tree requires all levels
  #to be defined
  mutate(Phylum = ifelse(Order == "Ebriales", "Cercozoa", Phylum),
         Phylum = ifelse(Family == "Katablepharidaceae","Cryptista", Phylum),
         Phylum = ifelse(Genus == "Telonema", "Cryptophyta", Phylum),
         Order = ifelse(Family == "Katablepharidaceae", "Katablepharidales",Order),
         Order = ifelse(Genus == "Telonema", "Telonemida", Order),
         Family = ifelse(Genus == "Telonema","Telonema hypothetical family",Family)) |>
  mutate(across(where(is.character),factor))
spp.tree = as.phylo(~Kingdom/Phylum/Order/Family/Genus/valid_name, data = to_tree, collapse = FALSE)
spp.tree$edge.length = rep(1, length(spp.tree$edge[,1]))


#### Check community matrix distibutions ####
S = rowSums(pa_matrix)
P = colMeans(pa_matrix)
A = colSums(sppMat)/sum(sppMat)
range(S)
range(P)
par(mfrow=c(2,2))
hist(S, xlab = "Species richness (S)")
hist(P, xlab = "Species prevalence (P)",xlim=c(0,1))
hist(A, xlab = "Species Abundance (A)")
hist(log(A,10), xlab = "Species Abundance (A)")

#### Defining the models ####

#Define the random effects
#Location
rl.coords = HmscRandomLevel(sData = xy, longlat = TRUE)

#Temporal distance (no logner used)
#rl.year = HmscRandomLevel(sData = as.matrix(Year))

#Time as a catagorical random level 
rl.yearCat = HmscRandomLevel(units = levels(studyDesign$year))

#When using HMSC-hpc nAdapt can not be set, model convergence is great but not
#perfect for PA_TSP_AT_day_depth_PhyCoYrCat_FULL_2017_NoOut at Samp=350,
#thin=350 at longer runs nf year Cat increases to 12 and nfcoords increase to 7
#which results in very poor convergence. To simulate setting nAdapt at 61250
#iterations while increase transient to 80000 the nfMax for both random factors
#are adjusted.
rl.yearCat = setPriors(rl.yearCat, nfMin = 1, nfMax = 11)
#rl.year = setPriors(rl.year, nfMin = 1, nfMax = 11)
rl.coords = setPriors(rl.coords, nfMin = 1, nfMax = 6)

#Check that everything is in the same sample order
all(studyDesign$sampleID==rownames(envMat)&studyDesign$sampleID==rownames(sppMat))

#Create the species matrix for abundance (biomass mg/L) conditional on occurrence 
Sum_sppMat_con = sppMat
Sum_sppMat_con[which(Sum_sppMat_con==0)] = NA
Sum_sppMat_log_con = log(Sum_sppMat_con)

#Note that the when using sData the rownames of the sData must match the units
#in the study Design. For example if using coordinates the names of each row
#must be the stn ID, for time it must be the time unit. This is because the
#dimensions of sData and X are not the same. The link is between rowname and
#studyDesign unit is defined in the Hmsc function section randLevels where the
#list(x=y) part must be x: unit in study design and y random level object. IE
#year = rl.year or stn = xy.

#Remember to manually define the XFormula based on the model you are trying to fit
#Change model type to Con for conditional abundance models
model_type = "Con"
Y = switch(model_type, "Con"=Sum_sppMat_log_con, "PA"=pa_matrix)
dist = switch(model_type,"Con" = "normal", "PA"="probit")
YScale = model_type=="Con"
m.basic = Hmsc(Y = Y, YScale = YScale,
               XData = envMat, XFormula = ~ poly(Salin,2,raw = TRUE) + 
                 poly(Temp,2,raw = TRUE) + poly(Ptot,2,raw=TRUE) +
                 poly(jday,2,raw=TRUE) + poly(Mean_april_temp,2,raw=TRUE) + 
                 stnDepth + year,
               TrData = traits, TrFormula = ~nfix+silica+motility+autotrophy+heterotrophy+mixotrophy,
               phyloTree = spp.tree,
               studyDesign = studyDesign,
               ranLevels = list(year = rl.yearCat, stn = rl.coords),
               distr=dist)

#This determines the name of the save folder for the model. The key that is
#currently being used is:
#Automatically add PA or Con to the start based on the type of model
#TSP means Temperature, Salinity and Phosphor
#d, dp, and/or yr for day, depth and year respectively if they are included as
#fixed effects
#P for phylo data, Tr for traits
#Then random effects C for coordinates, Y for year and Yc if year is not
#distance based
model_description = sprintf("%s_TSP_AT_d_dp_yr_PTCYC",model_type)

models = list(m.basic)

#Check if the model with this description has been created before, if not create
#both the model and results directories and all sub folders
ModelDir = file.path(localDir, sprintf("%s",model_description))
save.dir = file.path(ModelDir, "Models")
resluts.dir = file.path(ModelDir, "Results")

if(!dir.exists(save.dir)){
  dir.create(ModelDir)
  dir.create(save.dir)
  for(x in c("Fitted","Raw_HPC","Unfitted","INITS","Temp")){
    dir.create(file.path(save.dir,x))
  }
  if(!dir.exists(resluts.dir)){
    dir.create(resluts.dir)
  }
} else {
  cat(sprintf("\n%1$s\nWARNING\n%1$s\nThis model has been ran before, unfitted model file being over written\n",strrep("-",10)))
}

save(models, file = file.path(save.dir, "/Unfitted/unfitted_models.RData"))
