
#####################
# SCRIPT FOR RUNNING ON YOUR DATA

######################


install.packages("rioja")
install.packages("vegan")

library("vegan")
library("rioja")

mainDir <- "~/Desktop"
subDir <- "ForRes"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

setwd(file.path(mainDir, subDir))
source("functions.R")

### GET THE DATA STEP
dataFile <- "6_USA_Purple.csv"
core <- read.csv(dataFile, header = TRUE)
ages <- (core[,1])
core <- core[,-1]


strat.plot(core, yvar = ages, scale.percent= TRUE)
coreT <- sqrt(core)
plot(ages[-1], diff(ages))
coreCrop <- coreT[-c(1:22),]
agesCrop <- ages[-c(1:22)]

# build a list containing all the info needed to run the functions
purple <- list(core = coreCrop, ages = agesCrop, dataset = "purple")

# clean up the workstation
rm(coreCrop, agesCrop, core, ages)

# Calculate the BC distances
purpleBC <- calcBC(site= purple)

# Make the BC distance plot
plotBC(purpleBC)

# Choose how many perturbations you want and find out where they are in the record
distEvents <- extractDist(calcBC.obj = purpleBC, nDist = 2)

# Calculate the number of zones from the PCA and extract the ones you are interested in
PCs <-sigPC(site = maua)



# Plot a PC that you are interested in and also plot where the disturbance event is



###########################
###########################

# NOW YOU NEED TO ESTIMATE THE RECOVERY RATE YOURSELF....


###########################
###########################


# Now we need to build a table with all these recovery rates involved





