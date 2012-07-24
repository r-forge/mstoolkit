"analyzeData" <- function(
  replicates = "*",                  #@ Replicates to perform analysis on 
  analysisCode,                      #@ Function taking a data
  macroCode,                         #@ Macro evaluation code
  interimCode = NULL,                #@ Interim analysis code
  software = "R",                    #@ Software for analysis: R or SAS
  grid = TRUE,                       #@ Split analysis across the grid?
  waitAndCombine = TRUE,             #@ Wait for all analyses to finish, then combine into single file?
  cleanUp = FALSE,                   #@ Delete micro/macro directories on completion?
  removeMissing = TRUE,              #@ Remove Missing rows?
  removeParOmit = TRUE,              #@ Remove Parameter Omit rows?
  removeRespOmit = TRUE,             #@ Remove Response Omit rows?
  seed = .deriveFromMasterSeed(),    #@ Random number seed
  parOmitFlag = getEctdColName("ParOmit"),           #@ Parameter omit flag name
  respOmitFlag = getEctdColName("RespOmit"),         #@ Response omit flag name
  missingFlag = getEctdColName("Missing"),           #@ Missing flag name
  interimCol = getEctdColName("Interim"),            #@ Interim variable name
  doseCol = getEctdColName("Dose"),                  #@ Dose variable name
  sleepTime = 15,                    #@ Number of seconds to sleep between checking for grid jobs
  deleteCurrData = TRUE,             #@ Delete current analysis results before executing
  initialDoses = NULL,					#@ Initial doses to use for "Interim 1"
  stayDropped = TRUE,				#@ Dose dropping flag: if a dose is dropped, should it stay dropped?
  fullAnalysis = TRUE,			#@ Perform a full analysis
  workingPath = getwd(),              #@ Working path containing data
  method = getEctdDataMethod()
)
{
	###############################################################################
	# � Mango Solutions, Chippenham SN14 0SQ 2006
	# analyzeData.R Tue Jul 03 16:24:00 BST 2007 @447 /Internet Time/
	#
	# Author: Richard, Romain
	###############################################################################
	# DESCRIPTION: High level function to analyze simulated trial datasets
	# KEYWORDS: high, analyze                                                
	###############################################################################
	# TESTME
	funCall <- match.call()

	## Check network connectivity  
	macroCode <- .checkFun(macroCode, "data")
	replicates <- .checkReplicates( replicates, workingPath = workingPath, method = method)
	
	# If there is not 'Rlsf', use 'parallel' to do parallel processing
	grid.para <- FALSE
	if (require("doParallel", quietly = TRUE) && require("foreach", quietly = TRUE)) {
		nclusters <- parallel:::detectCores() - 1
		if (nclusters > 1 && nclusters <= length(replicates)) {
			grid.para <- grid
		}
	}
	
	if (grid && !.checkGridAvailable()) grid <- FALSE
	if (length(replicates) == 1) grid <- waitAndCombine <- FALSE

	## Check directories
	if (deleteCurrData) removeDirectories(c("Micro", "Macro"), workingPath = workingPath)
	createDirectories(c("MicroEvaluation", "MacroEvaluation"), workingPath = workingPath)
  
	## Split jobs and call grid
	if (grid) {
		funCall[[1]] <- as.name(".ectdSubmit")              # Call the .ectdSubmit function for LSF split
		funCall$grid <- funCall$waitAndCombine <- funCall$deleteCurrData <- FALSE     # Don't split grid job over grid or compile
		funCall$func <- "analyzeData"                       # Grid function to call is analyzeData
		funCall$debug <- TRUE                               # Set debug flag on the grid system
		funCall$packages <- c("MSToolkit", "MASS")          # Required packages
		if (software == "SAS") funCall$reqSas <- TRUE       # Need to queue on a SAS machine
		repSplit <- .splitGridVector(replicates)            # Split replicates for Grid execution
		gridJobs <- lapply(repSplit, function(i, call) {
			call$replicates <- i
			eval(call)
		}, call=funCall)
		evalTime <- Sys.time()                              # Store time at grid evaluation
	} else if (grid.para) {
		nreps <- ceiling(length(replicates) / nclusters )
		repSplit <- .splitGridVector(replicates, nreps)
		cl <- parallel:::makeCluster(nclusters)
		doParallel:::registerDoParallel(cl)
		`%dopar%` <- foreach:::"%dopar%"
		k <- 0
		tmp <- foreach:::foreach(k = 1:nclusters, .packages = c("MSToolkit", "MASS")) %dopar% {
			for (i in repSplit[[k]]) {
			
				microData <- analyzeRep(replicate = i, analysisCode = analysisCode, 
						interimCode = interimCode, software = software, removeMissing = removeMissing, 
						removeParOmit = removeParOmit, removeRespOmit = removeRespOmit, 
						seed = seed + i, parOmitFlag = parOmitFlag, respOmitFlag = respOmitFlag, 
						missingFlag = missingFlag, interimCol = interimCol, doseCol = doseCol, 
						initialDoses = initialDoses, stayDropped = stayDropped, fullAnalysis = fullAnalysis,
						workingPath = workingPath, method = method)
				
				# Write out data
				if (is.data.frame(microData) && nrow(microData)) {
					
					writeData(microData, i, "Micro", workingPath = workingPath)
					
					macroData <- macroEvaluation(microData, macroCode = macroCode, 
							interimCol = interimCol, doseCol = doseCol)
					
					writeData(macroData, i, "Macro", workingPath = workingPath)
				}
				else ectdWarning(paste("No return output from replicate", i))
			}
		}
		parallel:::stopCluster(cl)
	} else {
		
		# Loop through and analyze replicates
		for (i in replicates) {

			## TODO: Update analyzeRep and performAnalysis with data storage method ..
			microData <- analyzeRep(replicate = i, analysisCode = analysisCode, 
				interimCode = interimCode, software = software, removeMissing = removeMissing, 
				removeParOmit = removeParOmit, removeRespOmit = removeRespOmit, 
				seed = seed + i, parOmitFlag = parOmitFlag, respOmitFlag = respOmitFlag, 
	        	missingFlag = missingFlag, interimCol = interimCol, doseCol = doseCol, 
				initialDoses = initialDoses, stayDropped = stayDropped, fullAnalysis = fullAnalysis,
				workingPath = workingPath, method = method)
	
			# Write out data
			if (is.data.frame(microData) && nrow(microData)) {
	
				writeData(microData, i, "Micro", workingPath = workingPath)
	      
				macroData <- macroEvaluation(microData, macroCode = macroCode, 
					interimCol = interimCol, doseCol = doseCol)
	      
				writeData(macroData, i, "Macro", workingPath = workingPath)
			}
			else ectdWarning(paste("No return output from replicate", i))
		}
	}
	
	if (waitAndCombine) {   
		if (grid) {
			lsf.job.status <- get("lsf.job.status")
	      	gridStatus <- sapply(gridJobs, lsf.job.status)
	      	checkJobs <- gridStatus %in% c("EXIT", "DONE")
	      	iter <- 1
	      	while (any(!checkJobs)) {
	        	if (any(gridStatus == "DONE")) {
	          		compileSummary("Micro", workingPath = workingPath)
	          		compileSummary("Macro", workingPath = workingPath)
				}
		        ## Write log file
		        writeLogFile(gridStatus, evalTime, workingPath = workingPath)
		        iter <- iter + 1
		        if (iter > 1000) ectdStop("Job timed out")
		        Sys.sleep(sleepTime)
		        gridStatus <- sapply(gridJobs, lsf.job.status)
		        checkJobs <- gridStatus %in% c("EXIT", "DONE")
			}
	      	writeLogFile(gridStatus, evalTime, workingPath = workingPath)
	    }
    
    	compileSummary("Micro", workingPath = workingPath)
    	compileSummary("Macro", workingPath = workingPath)      
    
	}
	.cleanup( cleanUp = cleanUp, grid = grid, workingPath = workingPath )
  
	invisible()
}

