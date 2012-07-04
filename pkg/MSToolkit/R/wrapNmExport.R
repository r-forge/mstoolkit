"wrapRNMExport" <- function(sampList, obj, ctlFile, path = getwd(), digits = 4, 
		fixThetas = FALSE, writeFiles = TRUE, subDirs = "nmRuns") 
{
	# Stop if library not there
	if (!require(RNMExport)) ectdStop("RNMExport library not found")
	ZERO.TOLERANCE <- RNMExport:::.RNMExportEnv$ZERO.TOLERANCE
	# Create directory for files
	if (!file.exists(file.path(path, subDirs))) dir.create(file.path(path, subDirs))
	
	# Needed for NONMEM
	options(scipen = digits + 4) 
	options(digits = digits) 
	
	# Get the raw control file text
	nmp <- as.nmprob(ctl=obj@controlText, debug=TRUE)
	whereIsTheControlFile <- dirname(obj@controlFileInfo$fileName)
	# Extract the information from these
	# pars <- sampList[[1]]; ctl = nmp; digits = digits; fixVals = fixThetas; path=whereIsTheControlFile
	updatedCtls <- lapply(sampList, function(pars, ctl, digits, fixVals, path) {
				
				# Replace thetas
				thb <- thetas(ctl)
				newThetas <- round(pars$THETA, digits)
				# th <- thb[1]
				count <- 0
				#TODO: create vectorisation to make loops cleaner
				for(i in 1:length(thb)) {
					thv <- thetas(thb[i])
					for(j in 1:length(thv)){
						count <- count + 1
						if(count > length(newThetas))
							return(simpleError("Insufficient THETAS generated"))
						initial(thv[j]) <- newThetas[count]
						if (fixVals) {
							lower(thv[j]) <- NA
							upper(thv[j]) <- NA
							fixed(thv[j]) <- "FIXED"
							brackets(thv[j]) <- FALSE
						}
					}
					thetas(thb[i])<- thv
				}
				thetas(ctl) <- thb
				
				# Replace Sigmas
				smb <- sigmas(ctl)
				newSigmas <- round(pars$SIGMA, digits)
				# sm <- smb[1]
				diagonal <- 0
				#TODO: create vectorisation to make loops cleaner
				for(i in 1:length(smb)) {
					switch(as.character(not.set(block(smb[i]))),
							"TRUE"={
								smv <- sigmas(smb[i])
								for(j in 1:length(smv)){
									# moving along the diagonal
									diagonal <- diagonal + 1
									if(diagonal > dim(newSigmas)[2]){
										return(simpleError("Insufficient SIGMAS generated"))
									}
									# j <- 1
									initial(smv[j]) <- max(ZERO.TOLERANCE,newSigmas[diagonal, diagonal])
								}
							},
							{
								smv <- sigmas(smb[i])
								n <- blocksize(smb[i])
								# diagonal <- 0; irow <- icol <- 1
								if(diagonal + n > dim(newSigmas)[2]){
									return(simpleError("Insufficient SIGMAS generated"))
								}
								count <- 0
								for(irow in 1:n){
									# moving along the lower triangular
									diagonal <- diagonal + 1
									i <- (irow+1)%%n + 1
									for (icol in 1:irow){
										count <- count + 1
										if(count > length(smb[i])){
											return(simpleError("too many SIGMAS generated"))
										}
										initial(smv[count]) <- ifelse(irow==icol, max(ZERO.TOLERANCE,newSigmas[irow, icol]), newSigmas[irow, icol])
									}
								}
							},
					)
					sigmas(smb[i])<- smv
				}
				sigmas(ctl) <- smb
				
				# Replace Omegas
				omb <- omegas(ctl)
				newOmegas <- round(pars$OMEGA, digits)
				omv <- omegas(omb)
				
#				switch(block(om)=="BLOCK",
#						BLOCK={
#							n <- length(omv)
#							nRows <- as.integer( - 1 + (1 +sqrt(1 + 8*n))/2)
#						},{
#							nRows <- length(omv)                              
#						}
#				)
# 				Change only those that are NOT fixed
#				!fixed(omv)
				n <- dim(newOmegas)[2]
#				Cannot reset the size using methods such as initial
				count <- 0
				omv <- RNMExport:::NMOMEGA()
				for(irow in 1:n){
					for (icol in 1:irow){
						count <- count + 1
						omv[count] <- RNMExport:::NMOMEGA(
								initial=ifelse(irow==icol, max(ZERO.TOLERANCE,newOmegas[irow, icol]),newOmegas[irow, icol]),
								brackets=FALSE,
								fixed=FALSE)
					}
				}
				omegas(omb) <- omv
				block(omb) <- "BLOCK"
				blocksize(omb) <- n
				omegas(ctl) <- omb
				
				# set the $SIM record
				# Sould return a NULL simblock record!
				# force the line for now!
				sib <- sim(ctl)
				if(length(sib)==0)
					block(sib) <- inputtext(sib)<- 
							paste(" (", format(runif(1)*10e8, 9), ") NSUB=1 ONLYSIMULATION")
				sim(ctl) <- sib
				
				# get the datafile records
				dataFile <- datafile(ctl)
				
				datafile(datafile(dataFile)) <- 
						switch(.Platform$OS.type, 
								windows=gsub("/", "\\\\", shortPathName(file.path(path, datafile(dataFile)))),
								file.path(path, datafile(dataFile))
						)
				datafile(ctl) <- dataFile
				
				# Return updated control file object
				purge(ctl)
			}, ctl = nmp, digits = digits, fixVals = fixThetas, path=whereIsTheControlFile)
	
	whereAretheFiles <- character(0)
	if (writeFiles) {
		# Export each control file
		for (i in 1:length(updatedCtls)) {
			whichDir <- file.path(path, subDirs, paste("run", i, sep=""))
			if (!file.exists(whichDir)) dir.create(whichDir)
			whereAretheFiles[i] <- switch(.Platform$OS.type,
					windows=shortPathName(file.path(whichDir, ctlFile, fsep="\\")),
					file.path(whichDir, ctlFile)
			)
			sink(whereAretheFiles[i])
			show(updatedCtls[[i]])
			sink()
		}
	}
	
	invisible(list(controls=updatedCtls, locations=whereAretheFiles))
}
