
############# FUNCTIONS FOR PALAEORESILIENCE WORKSHOP
#` Calculate Taxonomic distances between samples 
#' @param site list with elements core (species data), ages and dataset (Name to add to figures)


calcBC <- function(site, method = "bray") {
	
	core <- site$core
	ages <- site$ages
	dataset <- site$dataset
	
	# Do the B-C metric
	turn <- as.matrix(vegdist(core, method = method))
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


