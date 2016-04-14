
############# FUNCTIONS FOR PALAEORESILIENCE WORKSHOP

calcBC <- function(site) {
	
	core <- site$core
	ages <- site$ages
	dataset <- site$dataset
	
	# Do the B-C metric
	turn <- as.matrix(vegdist(core, method = "bray"))
	BC1 <- diag(turn[-1, -ncol(turn)])
	BC2 <- diag(turn[-(1:2), -ncol(turn)])

	return(list(ages = ages, BC1 = BC1, BC2 = BC2, dataset = dataset))
}


plotBC<- function(BCobject, print.pdf= FALSE)	{
	
	# Make the decision plot
	with(BCobject, {
		if(print.pdf == TRUE){ 
		  pdf(paste("decisionPlot_", dataset, ".pdf", sep ="" ))
		}
		opar <- par(mfrow = c(3,2), mar = c(3,3,1,1), mgp = c(1.5, .5, 0), oma = c(0, 0, 1, 0))
		plot(ages[-1], BC1, type = "h")
		plot(ages[-c(1,2)], BC2, type = "h")
		hist(BC1)
		hist(BC2)
		qqnorm(BC1)
		qqline(BC1)
		qqnorm(BC2)
		qqline(BC2)
		title(main = dataset, outer = TRUE)
		par <- opar
		if(print.pdf == TRUE){
		  dev.off()
		}
	})
}

######################################
####null distance - based on resampling observed counts

#` Calculate Null distance between samples
#' @param spp data.frame of species data (Each column is one species)
#' @param counts Vector of pollen sums. Defaults to 300.
#' @param prob Probabilites of null distribution required.
#' @param nrep Number of trials for null model. Ideally 1000 but slow.
#' @param dist Distance metric
#' @details Null distribution of distances expected from counting errors
#' @result Matrix of distances at probabilities requested
#' @example 
#` data(BCI)
#` getNullDistances(BCI, rowSums(BCI))



#for each level, resample pollen counts, find distance between original and resampled counts, find xx% limit on distr. cf obseved ample to sample differences
getNullDistances <- function(spp, counts = 300, prob = c(0.5, 0.95), nrep = 100, dist = "bray"){
  print("This function is slow")
  spp <- spp[, order(names(spp))]
  sppNames <- names(spp)
  spp <- as.data.frame(t(spp))
  
  getNullDist <- function(sp, count){
    dists <- replicate(nrep, {#browser()
      newCounts <- sample(sppNames, prob = sp, replace = TRUE, size = count)
      newCounts <- table(c(sppNames, newCounts)) - 1
      bothCounts <- rbind(sp, newCounts)##need to check species order
      vegdist(bothCounts, method = "bray")
    })
    quantile(dists, prob = prob)  
  }
  
  res <- mapply(getNullDist, sp = spp, count = counts)    
  t(res)
}



####### FIGURE 1 PLOT (WITH ARGUMENTS)

nPCs <- 2

par(mfrow = c((nPCs+1), 1), mar = c(3,3,1,1))
pcAxes <- list(pc1 = rnorm(10), pc2 = rnorm(10))
ages <- 1:10
bc <- rbinom(9, size= 10, prob = 0.5)/10
agesBC <- ages[-1]
chgpt <- c(4,8)


plotFig1 <- function(resObject) {
	
	with(resObject, {
		for(i in 1:nPCs){
			plot(  ages, pcAxes[[i]], type = "o", pch =20, xlim = c(max(ages), min(ages)), ylab = 	"Age (cal yr BP)", xlab = names(pcAxes)[i])
		abline(v = chgpt, lty = 2, col = "red")
		
		}
		plot( ages[-1], bc, xlim = c(max(ages), min(ages)), type = "h", ylim = c(0,1))
		points(chgpt, rep(0.95, length(chgpt)), xlab = "Bray-Curtis Dissimilarity", ylab = "Recovery Rate", pch = 11, col = "blue")
	})
}





###### FIGURE 3 PLOT (WITH ARGUMENTS)

resTab <- data.frame(time= 1:10, BC = rnorm(10), recov = rnorm (10), threshold = rbinom(10, size=1, prob = 0.5))

# Table contains the disturbance event, the BC score for that disturbance event and the estimated recovery rate

plotRes <- function(resTab){
	
	with(resTab[resTab$threshold == 1,], plot(BC, recov, xlab = "Bray-Curtis Dissimilarity", ylab = "Recovery Rate", pch = 11, col = "blue"))
	with(resTab[resTab$threshold == 0, ], points(BC, recov, col = "red", pch = 20))
		
}



