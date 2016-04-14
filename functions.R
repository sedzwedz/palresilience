
############# FUNCTIONS FOR PALAEORESILIENCE WORKSHOP

calcBC <- function(site) {
	
	core <- site$core
	ages <- site$ages
	dataset <- site$dataset
	
	# Do the B-C metric
	turn <- as.matrix(vegdist(core, method="bray"))
	BC1 <- diag(turn[-1, -ncol(turn)])
	BC2 <- diag(turn[-(1:2), -ncol(turn)])

	return(list(ages= ages, BC1 = BC1, BC2= BC2, dataset = dataset))
}


plotBC<- function(BCobject, print.pdf= FALSE)	{
	
	# Make the decision plot
	with(BCobject, {
		if(print.pdf == TRUE) pdf(paste("decisionPlot_", dataset, ".pdf", sep ="" ))
		par(mfrow = c(3,2), mar = c(3,3,1,1), mgp = c(1.5, .5, 0))
		plot(ages[-1], BC1, type = "h")
		plot(ages[-c(1,2)], BC2, type = "h")
		hist(BC1)
		hist(BC2)
		qqnorm(BC1)
		qqline(BC1)
		qqnorm(BC2)
		qqline(BC2)
		par(mfrow = c(1,1))
		if(print.pdf == TRUE) dev.off()
		
	})
}

