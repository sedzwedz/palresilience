

#####################
# SCRIPT FOR RUNNING ON YOUR DATA

######################


# install.packages("rioja")
# install.packages("vegan")

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
dev.off()
coreT <- (core)
plot(ages[-1], diff(ages))
coreCrop <- coreT[-c(1:22),]
agesCrop <- ages[-c(1:22)]

# build a list containing all the info needed to run the functions
purple <- list(core = coreCrop, ages = agesCrop, dataset = "purple")
purple <- check.list(purple)


# clean up the workstation
rm(coreCrop, agesCrop, core, ages)

# Calculate the BC distances
BC <- calcBC(site= purple)

# Make the BC distance plot
plotBC(BC)

# Choose how many perturbations you want and find out where they are in the record
distEvents <- extractDist(BCobj = BC, nDist = 3, ord = 1)

# Check the number of significant axes
PCs <-sigPC(site = purple)

# Plot a PC that you are interested in and also plot where the disturbance event is

par(mfrow= c(2,1))
plot.pca.time(site = purple, pcaWant= 1, distEvents, PCs)
plot.pca.time(site = purple, pcaWant= 2, distEvents, PCs)


# define how many PCs you want
nPCs <- 2


###########################
###########################

# NOW YOU NEED TO ESTIMATE THE RECOVERY RATE YOURSELF....


recov <- c(50,100, 250)
recov <- cbind(distEvents, recov)


###########################
###########################



# Now we need to test for the significant number of zones using a cluster analysis
clust <- calcClust(site= purple)

# Decide that there is only one split (i.e. 2 groups)
nGroups <- 2

# Make a data.frame with the ages and the zones that can be used for plotting later
zones <- makeZones(clust, nGroups, site = purple)

# Now we need to build a list with all the parts recovery rates involved
resList <- list( ages = zones$ages, zones = zones, recov = recov, PC = PCs, BC = BC, nPCs = nPCs)


### NOW MAKE FIGURE 1
plotFig1(resList, nPCs, PCs)


# Now make FIgure 3
plotFig3(resList$recov)



