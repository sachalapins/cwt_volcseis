# SCRIPT TO READ DATA/HEADER FROM SAC BINARY FILE #

readSAC <- function(filename, header.only=F, trace.as.zoo=T) {

	# CHECK FOR zoo PACKAGE, IF trace.as.zoo = TRUE ####
	if (trace.as.zoo == T) {
		require(zoo)
	}
	
	# OPEN BINARY SAC FILE ####
	to.read <- file(filename, "rb")
	
	# READ HEADER VALUES (float32 and integer) ####
	hd <- readBin(to.read, double(), 70, size=4)
	ihd <- readBin(to.read, integer(), 40)

	# CHECK CORRECT ENDIAN ####
	# ie BY CHECKING HEADER VERSION NUMBER = 6
	endiantype <- .Platform$endian
	if (ihd[7] != 6) {
		if (endiantype == "little") {
			endiantype <- "big"
			to.read <- file(filename, "rb")
			hd <- readBin(to.read, double(), 70, size=4, endian=endiantype)
			ihd <- readBin(to.read, integer(), 40, endian=endiantype)
		} else {
			endiantype <- "little"
			to.read <- file(filename, "rb")
			hd <- readBin(to.read, double(), 70, size=4, endian=endiantype)
			ihd <- readBin(to.read, integer(), 40, endian=endiantype)		
		}
	}

	# CHECK IF HEADER VERSION NUMBER IS VALID AFTER SWITCH
	if (ihd[7] != 6) {
		stop("Error: SAC file header version number is out of date; i.e. not equal to 6")
	}

	# READ CHARACTER HEADER VALUES ####
	# Create empty vector to accomodate values
	chd <- NULL

	# There are 23 character header values to read
	for (i in 1:23) {
		if (i == 2) {
			chd[i] <- trimws(readChar(to.read, 16)) # The second value is length 16, the rest are length 8
		} else {
			chd[i] <- trimws(readChar(to.read, 8))
		}
	}
	
	# READ TRACE IF header.only = FALSE ####
	# (as per Wookey's matlab code - not sure whether x2 and x3 ever come up)
	if (header.only == F) {
		x1 <- readBin(to.read, double(), ihd[10], size=4, endian=endiantype)
		if (any(ihd[16] == c(2,3,4,51))) {
			x2 <- readBin(to.read, double(), ihd[10], size=4, endian=endiantype)
		}
		if (ihd[16] == 51) {
			x3 <- readBin(to.read, double(), ihd[10], size=4, endian=endiantype)
		}
	}
	
	# CLOSE CONNECTION TO SAC BINARY FILE ####
	close(to.read)
	

	# ORGANISE READ-IN VALUES ####
	# Create empty list to accomodate values
	sactrace <- list()
	
	# Assign 'double' header values to list
	for (i in 1:length(hd)) {
		sactrace[[i]] <- hd[i]
	}
	
	# Assign 'integer' header values to list
	for (i in 1:length(ihd)) {
		sactrace[[i+length(hd)]] <- ihd[i]
	}
	
	# Assign 'character' header values to list
	for (i in 1:length(chd)) {
		sactrace[[i+length(hd)+length(ihd)]] <- chd[i]
	}
	
	# Header variable names (in order)
	hdvars <- c('delta', 'depmin', 'depmax', 'scale', 'odelta', 'b', 'e', 'o', 'a', 'internal0', 
				't0', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 
				'f', 'resp0', 'resp1', 'resp2', 'resp3', 'resp4', 'resp5', 'resp6', 'resp7', 'resp8',
				'resp9', 'stla', 'stlo', 'stel', 'stdp', 'evla', 'evlo', 'evel', 'evdp', 'mag',
				'user0', 'user1', 'user2', 'user3', 'user4', 'user5', 'user6', 'user7', 'user8', 'user9',
				'dist', 'az', 'baz', 'gcarc', 'internal1', 'internal2', 'depmen', 'cmpaz', 'cmpinc', 'xminimum',
				'xmaximum', 'yminimum', 'ymaximum', 'unused1', 'unused2', 'unused3', 'unused4', 'unused5', 'unused6', 'unused7', #hd values end here
				'nzyear', 'nzjday', 'nzhour', 'nzmin', 'nzsec', 'nzmsec', 'nvhdr', 'norid', 'nevid', 'npts',
				'internal3', 'nwfid', 'nxsize', 'nysize', 'unused8', 'iftype', 'idep', 'iztype', 'unused9', 'iinst',
				'istreg', 'ievreg', 'ievtyp', 'iqual', 'isynth', 'imagtyp', 'imagsrc', 'unused10', 'unused11', 'unused12',
				'unused13', 'unused14', 'unused15', 'unused16', 'unused17', 'leven', 'lpspol', 'lovrok', 'lcalda', 'unused18', #ihd values end here
				'kstnm', 'kevnm', 'khole', 'ko', 'ka', 'kt0', 'kt1', 'kt2', 'kt3', 'kt4',
				'kt5', 'kt6', 'kt7', 'kt8', 'kt9', 'kf', 'kuser0', 'kuser1', 'kuser2', 'kcmpnm',
				'knetwk', 'kdatrd', 'kinst') #chd values end here

	# Assign header variable names to values
	names(sactrace) <- hdvars
	
	# Add trace data to list
	if (header.only == F) {
		if (trace.as.zoo == T) {
		    
		    sactrace$x1 <- zoo(x = x1, order.by = seq(from = as.POSIXct(paste(sactrace$nzyear, "-", sactrace$nzjday, " ", 
		    							sactrace$nzhour, ":", sactrace$nzmin, ":", sactrace$nzsec, ".", sactrace$nzmsec, sep=""), 
		    							tz="utc", format="%Y-%j %H:%M:%OS") + sactrace$b, by = sactrace$delta, length.out = sactrace$npts))
		
		} else {
		
			sactrace$x1 <- x1
		
		}
		
		if (exists('x2')) {
		
			if (trace.as.zoo == T) {
			    
			    sactrace$x2 <- zoo(x = x2, order.by = seq(from = as.POSIXct(paste(sactrace$nzyear, "-", sactrace$nzjday, " ", 
			                                 sactrace$nzhour, ":", sactrace$nzmin, ":", sactrace$nzsec, ".", sactrace$nzmsec, sep=""), 
			                                 tz="utc", format="%Y-%j %H:%M:%OS") + sactrace$b, by = sactrace$delta, length.out = sactrace$npts))
		
			} else {
		
				sactrace$x2 <- x2
		
			}
		
		}
		
		if (exists('x3')) {
		
			if (trace.as.zoo == T) {
		
				sactrace$x3 <- zoo(x = x3, order.by = seq(from = as.POSIXct(paste(sactrace$nzyear, "-", sactrace$nzjday, " ", 
			                                 sactrace$nzhour, ":", sactrace$nzmin, ":", sactrace$nzsec, ".", sactrace$nzmsec, sep=""), 
			                                 tz="utc", format="%Y-%j %H:%M:%OS") + sactrace$b, by = sactrace$delta, length.out = sactrace$npts))
			    
			} else {
		
				sactrace$x3 <- x3
		
			}

		}
	}
	
	# CLOSE ANY OPEN CONNECTIONS
	closeAllConnections()
	
	# RETURN LIST OF ALL VALUES AS OUTPUT ####
	return(invisible(sactrace))

}
