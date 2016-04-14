
#####################
# SCRIPT FOR RUNNING ON YOUR DATA

######################


install.packages("vegan")
install.packages("vegan")

library("vegan")
library("rioja")

mainDir <- "~/Desktop"
subDir <- "ForRes1"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

setwd(file.path(mainDir, subDir))
source("functions.R")

### GET THE DATA STEP
dataFile <- "6_USA_Purple.csv"


### GET THE DATA STEP
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

# Do a DCA 
dca<- decorana(purple$core)
summary(dca)
screeplot(dca, bstick = TRUE)

# Extract DC1
DC1 <- as.data.frame(scores(dca))$DCA1

# plot DC1 against time
with(purple, plot(ages, DC1, type= "o", pch = 20))

