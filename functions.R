
############# FUNCTIONS FOR PALAEORESILIENCE WORKSHOP
#' Checks input data and corrects if necessary
#' @param site list with elements core (species data), ages and dataset (Name to add to figures)

check.list <- function(site){
  if(any(diff(site$ages) < 0)){
    stopifnot(all(diff(site$ages) < 0))
    site$core <- site$core[order(site$ages), ]
    site$ages <- sort(site$ages)
  }
  if(all(colSums(site$core) < 1))#check if percent. Convert otherwise 
    {
    site$core <- site$core * 100
  }
  print(paste("N samples  = ", length(site$ages)))
  print(paste("Median resolution  = ", median(diff(site$ages))))
  print(paste("Median sum percent  = ", median(rowSums(site$core))))
              
  site 
}



###calcBC
#` Calculate Taxonomic distances between samples 
#' @param site list with elements core (species data), ages and dataset (Name to add to figures)
#' @param method Distance metric recognised by vegdist
#' @param makeNullDistances logical make NULL distances

calcBC <- function(site, method = "bray", makeNullDistances = TRUE) {
  require(vegan)
  res <- with(site, {
	
  	# Do the distance metric
  	turn <- as.matrix(vegdist(core, method = method))
  	BC1 <- diag(turn[-1, -ncol(turn)])
  	BC2 <- diag(turn[-(1:2), -ncol(turn)])
    if(is.null(site$counts)){
      counts <-300
    }
    if(makeNullDistances){
      nullDistances <- getNullDistances(spp = core, counts = counts, prob = c(0.95), nrep = 100, method = method)
    }else{
      nullDistances <- NULL
    }
  	list(ages = ages[-length(ages)], BC1 = BC1, BC2 = BC2, dataset = dataset, nullDistances = nullDistances)
  })
	return(res)
}


###plotBC
#` Plot Taxonomic distances between samples 
#' @param BCobject result of plotBC
#' @param print.pdf logical. Print to pdf. If FALSE, prints to screen.

plotBC<- function(BCobject, print.pdf= FALSE)	{
	
	# Make the decision plot
	with(BCobject, {
		if(print.pdf == TRUE){ 
		  pdf(paste("decisionPlot_", dataset, ".pdf", sep ="" ))
		}
		opar <- par(mfrow = c(3,2), mar = c(3,3,1,1), mgp = c(1.5, .5, 0), oma = c(0, 0, 1, 0))
    on.exit(par(opar))
		plot(ages, BC1, type = "h")
    if(!is.null(nullDistances)){
      points(ages[-1], nullDistances[-1], col = "blue", pch = "-")
    }
		plot(ages[-2], BC2, type = "h")
		hist(BC1)
		hist(BC2)
		qqnorm(BC1)
		qqline(BC1)
		qqnorm(BC2)
		qqline(BC2)
		title(main = dataset, outer = TRUE)
		if(print.pdf == TRUE){
		  dev.off()
		}
	})
}

############### IDENTIFY THE DISTURBANCE TIME POINTS 

extractDist <- function(BCobject, nDist, ord = 1){
  if(ord == 1) {
    obj  <- BCobject$BC1
  }
  else{
    obj <- BCobject$BC2
  }
	want <- order(obj, decreasing = TRUE)[1:nDist]
	distEvents <-	data.frame(ages = BCobject$ages[want], magnitude = obj[want])
  distEvents <- distEvents[order(distEvents$ages),]
	return(distEvents)
}


### FIND THE SIGNIFICANT PCS

sigPC <- function(site) {
  core <- site$core

  # Do a PCA and check the screeplot
  pca <- rda(core)
  #  summary(pca)

  bstick <- screeplot(pca, bstick = TRUE)

  # Extract Principal components and return them
  PC <- scores(pca)$sites # only gives the first two PCs
  return(PC)
  
}

###### Plot the PC axes for estimating plots

plot.pca.time <- function(site, pcaWant, distEvents, PCs){

  # Plot PCs against time
  with(site, plot(ages, PCs[,pcaWant], type = "o", pch = 20))
  abline(v = distEvents$ages, lty= 2, col= "red")
   
}


# Calculate the cluster analysis
calcClust <- function( site) {
	opar <-par(mfrow = c(2,2), mar=c(5,5,5,5))
  on.exit(par(opar))
	diss <- vegdist(site$core, method="bray")
	clust  <-  chclust(diss)
	plot(clust, hang = -1, main = "Constrained Cluster Analysis", cex = 0.4)
	bstick(clust, 10)
	return(clust)
}


# Make the zones
makeZones <- function(clust, nGroups, site){
	core <- site$core
	ages <- site$ages
	sigClust <- cutree(clust, k=nGroups)
	zones <- data.frame(ages = ages, zone = sigClust)
	return(zones)
}




######################################
####null distance - based on resampling observed counts

#` Calculate Null distance between samples
#' @param spp data.frame of species data (Each column is one species)
#' @param counts Vector of pollen sums. Defaults to 300.
#' @param prob Probabilites of null distribution required.
#' @param nrep Number of trials for null model. Ideally 1000 but slow.
#' @param method Distance metric
#' @details Null distribution of distances expected from counting errors
#' @result Matrix of distances at probabilities requested
#' @example 
#` data(BCI)
#` getNullDistances(BCI, rowSums(BCI))



#for each level, resample pollen counts, find distance between original and resampled counts, find xx% limit on distr. cf obseved ample to sample differences
getNullDistances <- function(spp, counts = 300, prob = c(0.5, 0.95), nrep = 100, method = "bray"){
  sppNames <- names(spp)
  spp <- as.data.frame(t(spp))
  
  getNullDist <- function(sp, count){
    newCounts <- rmultinom(n = nrep, size = count, prob = sp)
    newCounts <- newCounts / count * 100
    dists <- apply(newCounts, 2, function(nc){
      bothCounts <- rbind(sp, nc)
      vegdist(bothCounts, method = "bray")
    })
    quantile(dists, prob = prob)  
  }
  
  res <- mapply(getNullDist, sp = spp, count = counts)    
  if(length(prob) >1) res <- t(res)
  res
}


####### FIGURE 1 PLOT 


plotFig1 <- function(resList, nPCs, PCs) {
	
	opar <- par(mfrow = c((nPCs+1), 1), mar = c(3, 3, 1, 1))
	on.exit(par(opar))
	with(resList, {
		
		for(i in 1:nPCs){
			plot( ages, PCs[,i], type = "o", pch =20, xlim = c(max(ages), min(ages)), ylab = 	"Age (cal yr BP)", xlab = names(PC)[i])
			abline(v = recov$ages, lty = 2, col = "red")
		
		}
		plot(ages[-1], BC$BC1, xlim = c(max(BC$ages), min(BC$ages)), type = "h")
				
		points(recov$ages, rep(max(BC$BC1), length(recov$ages)), xlab = "Bray-Curtis Dissimilarity", ylab = "Recovery Rate", pch = 11, col = "blue")
		
		zones$null <- max(BC$BC1)
		by(zones, zones$zone, function(x){
		  lines(x$ages, x$null, col = (x$zone+1))
		})
	})
}




###### FIGURE 3 PLOT 

# Table contains the disturbance event, the BC score for that disturbance event and the estimated recovery rate
plotFig3 <- function(recov){
	with(recov, {
       plot(magnitude, recov, xlab = "Bray-Curtis Dissimilarity", ylab = "Recovery Rate", pch = 20, col = "blue")
       })
}


