##########################################################################################################
##########################################################################################################
#####                                                                                                #####
#####                    AESOP: Analysis of Electrostatic Similarities of Proteins                   #####
#####                                Last update: October 06, 2013                                   #####
#####                                                                                                #####
#####                   Chris A. Kieslich, Ronald D. Gorham Jr., and Dimitrios Morikis               #####
#####                 University of California, Riverside; Department of Bioengineering              #####
#####                                                                                                #####
#####          Correspondence should be directed to Prof. Dimitrios Morikis at dmorikis@ucr.edu      #####
##### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #####
#####                                                                                                #####
#####           Copyright (C) 2013 Chris A. Kieslich, Ronald D. Gorham Jr., and Dimitrios Morikis    #####
#####                                                                                                #####
#####           This program is free software: you can redistribute it and/or modify                 #####
#####           it under the terms of the GNU General Public License as published by                 #####
#####           the Free Software Foundation, either version 3 of the License, or                    #####
#####           (at your option) any later version.                                                  #####
#####                                                                                                #####
#####           This program is distributed in the hope that it will be useful,                      #####
#####           but WITHOUT ANY WARRANTY; without even the implied warranty of                       #####
#####           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                        #####
#####           GNU General Public License for more details.                                         #####
#####                                                                                                #####
#####           You should have received a copy of the GNU General Public License                    #####
#####           along with this program.  If not, see <http://www.gnu.org/licenses/>.                #####
#####                                                                                                #####
##########################################################################################################
##########################################################################################################

mut2ala <- function(pqr,rnum,ch,ff){
	# mut2ala mutates a specified amino acid residue to alanine by truncating the side chain. mut2ala 
	# takes as an input a PDB/PQR object, as obtained from read.pdb/pqr (bio3d), as well as a residue  
	# number and chain identifier of the residue to be mutated. The truncation is performed by converting 
	# the CG atom of the residue to be mutated to HB1, shortening the CB-HB1 bond length, and removing 
	# unneeded atoms. mut2ala returns a PDB/PQR object for the mutated protein.
	
	## Input Details ##
	#-----------------#
	# pqr         - PDB/PQR object (bio3d) - Contains protein structure to be mutated
	# rnum        - numeric                - Residue number of residue to be mutated
	# ch          - character              - Single letter chain identifier of residue to be mutated
	# ff		  - force field            - PARSE or CHARMM (for use with ions)
	
	
	res <- atom.select(pqr,resno = rnum,chain = ch,verbose=F) # Selects residue to be mutated (see bio3d)
	
	CB  <- atom.select(pqr,resno = rnum,chain = ch,elety="CB",verbose=F) # Selects the CB atom of residue to be mutated
	CG  <- atom.select(pqr,resno = rnum,chain = ch,elety="CG",verbose=F) # Selects the CG atom of residue to be mutated 
	HB2 <- atom.select(pqr,resno = rnum,chain = ch,elety="HB2",verbose=F) # Selects the HB2 atom of residue to be mutated
	HB3 <- atom.select(pqr,resno = rnum,chain = ch,elety="HB3",verbose=F) # Selects the HB3 atom of residue to be mutated
	
	pqr$atom[res$atom,"resid"] <- "ALA" # Changes residue name to alanine (ALA)
	
	if (ff == "PARSE") {
		pqr$atom[CB$atom, "o"] <- 0.0000 # Sets the charge of the CB atom to 0
		pqr$atom[CG$atom, "o"] <- 0.0000 # Sets the charge of the CG atom to 0
		pqr$atom[CG$atom, "b"] <- 0.0000 # Sets the radius of the CG atom to 0
	} else if (ff == "CHARMM") {				#Added conditions for CHARMM forcefield to account for ions and non-integer charges | RM & RG 02/20/14
		pqr$atom[CB$atom, "o"] <- 0.0000 # Sets the charge of the CB atom to 0
		pqr$atom[CG$atom, "o"] <- 0.0000 # Sets the charge of the CG atom to 0
		pqr$atom[HB2$atom, "o"] <- 0.0000 # Sets the charge of the HB2 atom to 0
		pqr$atom[HB3$atom, "o"] <- 0.0000 # Sets the charge of the HB3 atom to 0
		pqr$atom[CG$atom, "b"] <- 1.3200 # Sets the radius of the CG atom to 1.32
	}
	pqr$atom[CG$atom, "elety"] <- "HB1" # Changes the name of the CG atom to HB1
	
	HB1 <- atom.select(pqr, resno = rnum,chain = ch,elety="HB1",verbose=F) # Selects the HB1 atom of the mutated residue
	pqr$xyz[HB1$xyz] <- .7105*(pqr$xyz[HB1$xyz] - pqr$xyz[CB$xyz]) + pqr$xyz[CB$xyz] # Shortens the CB-HB1 bond 
	
	# Selects the atom of the truncated residue that should be kept 
	keep <- atom.select(pqr,resno = rnum,elety=c("N","H","H2","H3","CA","HA","CB","HB1","HB2","HB3","C","O","OXT"),verbose=F)
	# Identifies the element numbers (eleno) of atoms to removed 
	rem.at <- as.numeric(pqr$atom[res$atom,"eleno"][!(pqr$atom[res$atom,"eleno"] %in% pqr$atom[keep$atom,"eleno"])])
		
	#Removes unneeded atoms
	rem.ind <- atom.select(pqr,eleno=rem.at,verbose=F) # Selects atoms to be removed
	pqr$atom <- pqr$atom[-rem.ind$atom,] # Removes atoms from the ATOM entries of the PDB/PQR
	pqr$xyz <- pqr$xyz[-rem.ind$xyz] # Removes atoms from the XYZ entries
	
	pqr$atom[,"eleno"] <- 1:length(pqr$atom[,"eleno"]) # Renumbers atoms
	return(pqr)
} 
## Modified ala.scan and ala.scan.c to use write.pqr3, rather than write.pqr - Ron (1/15/2014) ##
ala.scan <- function(pdb.file,ff="PARSE"){
	# ala.scan performs an alanine scan of charged residues (D, E, K, R, and H). A PQR file for the parent
	# structure is first generated using PDB2PQR. mut2ala is then utilized to generate alanine mutants by
	# truncating sides in the parent PQR file. This version is intended for a single protein structure.
	# (No complex)  
	
	## Input Details ##
	#-----------------#
	# pdb.file    - character              - Name of parent PDB file
	# ff          - character              - PDB2PQR force field [default = PARSE]
	

	protein.name <- strsplit(pdb.file,split=".pdb")[1]
	
	pqr.dir <- paste(getwd(),"/",paste(protein.name,"pqr",sep="_"),"/",sep="")
	if (file.exists(path = pqr.dir) == 0){
		dir.create(pqr.dir)
	}
	
	pqr.file <- paste(pqr.dir,"/",protein.name,".pqr",sep="")	
	system(paste('python',paste(paths$pdb2pqr,'pdb2pqr.py',sep=""),
		     '--chain ',paste('--ff=',ff,sep=""),pdb.file,pqr.file,sep=" "))
	
	# Read PQR
	write("\n",file=pqr.file,append=T)
	pqr <- read.pqr(pqr.file)
	header <- readLines(pqr.file)
	header <- header[grep("REMARK", header)]
	
	# Names of possible titratable groups (only standard groups available currently)       
	ion.name <- c("H","E","D","K","R")
		
	# Extract amino acid sequence, residue numbers, and chain ids
	aaseq <- pdbseq(pqr) #seq.pdb renamed to pdbseq in bio3d Rohith
	resnos <- as.numeric(pqr$atom[pqr$calpha,"resno"])
	chain <- pqr$atom[pqr$calpha,"chain"]
	
	# Determine which residue out of the entire sequence to mutate
	ion.nos <- resnos[aaseq %in% ion.name]
	ion.ids <- aaseq[aaseq %in% ion.name]
	ion.chain <- chain[aaseq %in% ion.name]

	for(i in 1:length(ion.nos)){
		mut.pqr <- mut2ala(pqr = pqr,rnum = ion.nos[i],ch = chain[i], ff) #Added ff so that it is passed through to mut2ala | RM & RG 02/20/14
		charge <-  formatC(sum(as.numeric(mut.pqr$atom[,"o"])),digits=4,format="f")
		out.name <- paste(pqr.dir,paste(protein.name,ion.nos[i],
				  paste(ion.ids[i],"2A.pqr",sep=""),sep="_"),sep="")
	
		LineToChange <- header[grep("Total", header)]
		X<-strsplit(LineToChange, split=" ")
		X[[1]][length(X[[1]])-1] <- charge
		X<-paste(X[[1]], collapse=" ")
		header[grep("Total", header)]<-X
	
		writeLines(header, con=out.name)
		write.pqr3(mut.pqr, file=out.name, append=TRUE)
	}

	return(list(dirs=pqr.dir,list = paste(ion.chain,ion.nos,ion.ids,sep="-")))

}

ala.scan.c <- function(pdb.file,chain.names,chains=unique(pdb$atom[,"chain"]),ff="PARSE"){
	# ala.scan performs an alanine scan of charged residues (D, E, K, R, and H). A PQR file for the parent
	# structure is first generated using PDB2PQR. mut2ala is then utilized to generate alanine mutants by
	# truncating sides in the parent PQR file. This version is intended for a protein complex.  
	
	## Input Details ##
	#-----------------#
	# pdb.file    - character          - Name of parent PDB file
	# chain.names - character	       - Array of protein chain names to be used in output names
	# chains      - character	       - Array of chain IDs to be mutated [default = all included chains]
	# ff          - character          - PDB2PQR force field [default = PARSE]
	
	protein.name <- paste(chain.names,collapse="_")
	dirs <- c(paste(getwd(),"/",paste(protein.name,"pqr",sep="_"),sep=""),
			  paste(getwd(),"/",paste(chain.names[1],"pqr",sep="_"),sep=""),
			  paste(getwd(),"/",paste(chain.names[2],"pqr",sep="_"),sep=""))
	
	if(file.exists(path = dirs[1]) == 0){
		dir.create(dirs[1])
	}
	
	if(file.exists(path = dirs[2]) == 0){
		dir.create(dirs[2])
	}
	
	
	if(file.exists(path = dirs[3]) == 0){
		dir.create(dirs[3])
	}
	
	pqr.file <- paste(dirs[1],"/",paste(protein.name,".pqr",sep=""),sep="")
	system(paste('python',paste(paths$pdb2pqr,'pdb2pqr.py',sep=""),
		     '--chain ',paste('--ff=',ff,sep=""),pdb.file,pqr.file,sep=" "))
	
	# Read PQR
	write("\n",file=pqr.file,append=T)
	pqr <- read.pqr(pqr.file)
	header <- readLines(pqr.file)
	header <- header[grep("REMARK", header)]
	write.pqr3(pqr,file=pqr.file)
	# Names of possible titratable groups (only standard groups available currently)       
	ion.name <- c("H","E","D","K","R")
		
	# Extract amino acid sequence, residue numbers, and chain ids
	aaseq <- pdbseq(pqr) #seq.pdb renamed to pdbseq in bio3d Rohith
	resnos <- as.numeric(pqr$atom[pqr$calpha,"resno"])
	chain <- pqr$atom[pqr$calpha,"chain"]
	chains <- unique(chain)
	
	write.pqr3(trim.pdb(pqr,inds=atom.select(pqr,chain=chains[1],verbose=F)),
		  file=paste(dirs[2],"/",paste(chain.names[1],".pqr",sep=""),sep=""))
	write.pqr3(trim.pdb(pqr,inds=atom.select(pqr,chain=chains[2],verbose=F)),
		  file=paste(dirs[3],"/",paste(chain.names[2],".pqr",sep=""),sep=""))
	
	# Determine which residue out of the entire sequence to mutate
	ion.nos <- resnos[aaseq %in% ion.name & chain %in% chains]
	ion.ids <- aaseq[aaseq %in% ion.name & chain %in% chains]
	ion.chain <- chain[aaseq %in% ion.name & chain %in% chains]
	 
	for(i in 1:length(ion.nos)){
		mut.pqr <- mut2ala(pqr = pqr,rnum = ion.nos[i],ch = ion.chain[i], ff) #Added ff so that it is passed through to mut2ala | RM & RG 02/20/14
		out.name <- paste(dirs[1],"/",paste(protein.name,paste(ion.nos[i],ion.chain[i],sep=""),
			          paste(ion.ids[i],"2A.pqr",sep=""),sep="_"),sep="")
		charge <-  formatC(sum(as.numeric(mut.pqr$atom[,"o"])),digits=4,format="f")
		LineToChange <- header[grep("Total", header)]
		X<-strsplit(LineToChange, split=" ")
		X[[1]][length(X[[1]])-1] <- charge
		X<-paste(X[[1]], collapse=" ")
		header[grep("Total", header)]<-X
		writeLines(header, con=out.name)
		write.pqr3(mut.pqr, file=out.name, append=TRUE)
		
		f.mut.pqr <- trim.pdb(mut.pqr,ind=atom.select(mut.pqr,chain=ion.chain[i],verbose=F))
		out.name <- paste(dirs[1+which(chains==ion.chain[i])],"/",
				  paste(chain.names[which(chains==ion.chain[i])],
				  paste(ion.nos[i],ion.chain[i],sep=""),
				  paste(ion.ids[i],"2A.pqr",sep=""),sep="_"),sep="")
		charge <-  formatC(sum(as.numeric(f.mut.pqr$atom[,"o"])),digits=4,format="f")
		LineToChange <- header[grep("Total", header)]
		X<-strsplit(LineToChange, split=" ")
		X[[1]][length(X[[1]])-1] <- charge
		X<-paste(X[[1]], collapse=" ")
		header[grep("Total", header)]<-X
		writeLines(header, con=out.name)
		write.pqr3(f.mut.pqr, file=out.name, append=TRUE)
	}
	
	return(list(dirs = dirs, list = paste(ion.chain,ion.nos,ion.ids,sep="-")))
}

batch.pqr <- function(pdb.dir,pqr.dir,ff="PARSE",flags=NULL){
	# batch.pqr generates PQR files for all PDb structures in a directory (pdb.dir) using PDB2PQR. 
	
	## Input Details ##
	#-----------------#
	# pdb.dir     - character              - Name of directory of PDB files
	# pqr.dir     - character              - Name of output directory for PQR files		
	# ff          - character              - PDB2PQR force field [default = PARSE]
	# flags       - character              - Array of PDB2PQR input flags [default = NULL] 	

	files=list.files(path=pdb.dir) # List of PDB files to be 
	if (!file.exists(path = pqr.dir)) { # Generates folder for PQRs
		dir.create(pqr.dir)
	}
	pqr.list <- NULL
	for (file in files)	{
		# PDB files are first cleaned by keeping only ATOM lines in PDB
		pdb.file <- paste(pdb.dir,"/",file,sep="")
		pdb <- readLines(pdb.file) 
		clean.pdb <- pdb[grep("ATOM ",pdb)] 
		clean.pdb <- c(clean.pdb,"END")
		writeLines(clean.pdb,con=pdb.file)

		# Makes name for outputted PQR file and calls PDB2PQR 		
		protein.name <- strsplit(file,split=".pdb")[1]
		pqr.file <- paste(protein.name,".pqr",sep="")
		pqr.list <- c(pqr.list,pqr.file)
		system(paste('python',paste(paths$pdb2pqr,'pdb2pqr.py',sep=""),
			     '--chain',paste(flags,collapse=" "),
		       	     paste('--ff=',ff,sep=""),paste(pdb.dir,"/",file,sep=""),
			     paste(pqr.dir,"/",pqr.file,sep=""),sep=" "))
				 
		# Add in extra white space with write.pqr3 -RH 04/08/2014
		pqr.edit <- read.pqr(paste(pqr.dir,"/",pqr.file,sep=""))
		header <- readLines(paste(pqr.dir,"/",pqr.file,sep=""))
		header <- header[grep("REMARK", header)]
	
		write.pqr3(pqr.edit,file=paste(pqr.dir,"/",pqr.file,sep=""))
	}
	return(list(dirs=pqr.dir,list=pqr.list))
}

### Added 3 new arguments to function (avail.mem, grid.sp, and cfac) - Ron (1/6/2014) ###
apbsin.solv <- function(pdb,cent.on,pqr.dir,avail.mem = 2000,
					    grid.sp   = 1.0,
                                            cfac      = 1.5,
                                            r.keys    = NULL,
					    s.keys    = NULL,
					    s.only    = FALSE,
					    pot       = FALSE,
					    c.path    = NULL){

	# apbsin.solv generates an APBS input object that contains parameters for electrostatics
	# calculations in both a solvated and reference state.  
	
	## Input Details ##
	#-----------------#
	# pdb         - PDB/PQR object (bio3d) - Contains protein structure of parent
	# cent.on     - character              - Name of PQR used to center calculations
	# pqr.dir     - character              - Name of directory containing PQR files		
	# r.keys      - character              - Array of APBS key words to be changed for reference state [default = NULL]
	# s.keys      - character              - Array of APBS key words to be changed for solvated state [default = NULL]
	# s.only      - True/False             - Only perform calculations for solvated state [default = FALSE] 	
	# pot         - True/False             - Write potentials for solvated state [default = FALSE] 	
	# c.path      - character              - Path to PQR file to be used as center (cent.on) [default = NULL] 	

	x <- as.numeric(pdb$atom[,"x"])
	y <- as.numeric(pdb$atom[,"y"])
	z <- as.numeric(pdb$atom[,"z"])
	
	#cg <- fg <- ceiling(c(ceiling(max(x)- min(x))+5,ceiling(max(y)- min(y))+5,ceiling(max(z)- min(z))+5)*2)
	#dimes <- c("161","129","97","65")
	#det_dim <- (fg < 65) + (fg < 97) + (fg < 129) + (fg < 161)
	#gdime <- dimes[det_dim]

	########################################################################
	### This is a new method of grid size determination - Ron (1/8/2014) ###
	cg <- fg <- ceiling(c(ceiling(max(x)- min(x))+5,ceiling(max(y)- min(y))+5,ceiling(max(z)- min(z))+5)*cfac)
	dimes <- 32*(1:100)+1
	det_dim <- ceiling(fg/(32*grid.sp))
	gdime <- dimes[det_dim]
	mem_req <- (250/1024/1024)*gdime[1]*gdime[2]*gdime[3]
	if (mem_req > avail.mem) {
		print("Not enough memory available for calculation. Parallel solve required!!")
		break
	}
	########################################################################	

	k <- NULL
	
	# See APBS documentation for details on keywords		    
#	k[["solv"]] <- c("   mg-auto",
#			    "dime 129 129 129",
#	  		     paste("cglen",paste(cg,collapse=" "),sep=" "),
#			    "cgcent mol 2",
#			     paste("fglen",paste(fg,collapse=" "),sep=" "),
#			    "fgcent mol 2",
#			    "mol 1",
#			    "lpbe",
#			    "bcfl sdh",
#			    "srfm smol",
#			    "chgm spl2",
#			    "ion 1 0.150 2.0",
#			    "ion -1 0.150 2.0",
#			    "pdie  20.0",
#			    "sdie  78.54",
#			    "sdens  10.0",
#			    "srad  0.0",
#			    "swin  0.3",
#			    "temp  298.15",
#			    "calcenergy total")

	##########################################################################
	### Use 'mg-manual' instead of 'mg-auto', since cg=fg - Ron (1/8/2014) ###
	
	k[["solv"]] <- c(   "   mg-manual",
			    paste("dime",paste(gdime,collapse=" "),sep=" "),
	  		    paste("glen",paste(cg,collapse=" "),sep=" "),
			    "gcent mol 2",
			    "mol 1",
			    "lpbe",
			    "bcfl sdh",
			    "srfm smol",
			    "chgm spl2",
			    "ion 1 0.150 2.0",
			    "ion -1 0.150 2.0",
			    "pdie  20.0",
			    "sdie  78.54",
			    "sdens  10.0",
			    "srad  0.0",
			    "swin  0.3",
			    "temp  298.15",
			    "calcenergy total")
	##########################################################################


	if(!is.null(s.keys)){
		for(i in 1:length(s.keys)){
			if (strsplit(s.keys[i],split=" ")[[1]][1] == "ion"){
				k$solv[grep(paste(strsplit(s.keys[i],split=" ")[[1]][1],strsplit(s.keys[i],split=" ")[[1]][2],sep=" "),k$solv)] <- s.keys[i]
			}
			else{
				k$solv[grep(strsplit(s.keys[i],split=" ")[[1]][1],k$solv)] <- s.keys[i]
			}
		}
	}
				
	if(!s.only){
#		k[["ref"]] <- c("    mg-auto",
#				    "dime 129 129 129",
#		  		     paste("cglen",paste(cg,collapse=" "),sep=" "),
#				    "cgcent mol 2",
#				     paste("fglen",paste(fg,collapse=" "),sep=" "),
#				    "fgcent mol 2",
#				    "mol 1",
#				    "lpbe",
#				    "bcfl sdh",
#				    "srfm smol",
#				    "chgm spl2",
#				    "pdie  20.0",
#				    "sdie  20.0",
#				    "sdens  10.0",
#				    "srad  0.0",
#				    "swin  0.3",
#				    "temp  298.15",
#				    "calcenergy total")

	##########################################################################
	### Use 'mg-manual' instead of 'mg-auto', since cg=fg - Ron (1/8/2014) ###
		
		k[["ref"]] <- c(   "   mg-manual",
				   paste("dime ",gdime[1]," ",gdime[2]," ",gdime[3],sep=""),
		  		   paste("glen",paste(cg,collapse=" "),sep=" "),
				   "gcent mol 2",
				   "mol 1",
				   "lpbe",
				   "bcfl sdh",
				   "srfm smol",
				   "chgm spl2",
				   "pdie  20.0",
				   "sdie  20.0",
				   "sdens  10.0",
				   "srad  0.0",
				   "swin  0.3",
				   "temp  298.15",
				   "calcenergy total")
	##########################################################################

		if(!is.null(r.keys)){
			for(i in 1:length(r.keys)){
				k$ref[grep(strsplit(r.keys[i],split=" ")[[1]][1],k$ref)] <- r.keys[i]
			}
		}
	}
	
	apbsin <- list(dir=pqr.dir,mols=cent.on,center=cent.on,center.path=c.path,states=k,fglen=fg,cglen=cg,dime=gdime)
	
	if(pot){
		p.key <- NULL
		p.key[["format"]] <- "dx" 
		p.key[["file"]] <- "pot" 
		apbsin[["pot"]] <- p.key
	}
	
	return(apbsin)	
}


write.apbsin <- function(apbsin,f.name="apbs_solv.in",pot=FALSE){
	# write.apbsin writes APBS input object to file in appropriate APBS format.  
	
	## Input Details ##
	#-----------------#
	# apbsin      - APBS input object      - APBS input object (apbsin.solv)
	# f.name      - character              - Output file name [default = apbs_solv.in]
	# pot         - True/False             - Write potentials for solvated state [default = FALSE] 	

	infile <- file(f.name,"w")  # open an output file connection
	cat("read",sep="\n",file = infile)
	for(mol in apbsin$mols){
		cat(paste("   mol pqr ",apbsin$dir,"/",mol,sep=""),sep="\n",file = infile)
	}
	if(!is.null(apbsin$center.path)){
		cat(paste("   mol pqr ",apbsin$center.path,"/",apbsin$center,sep=""),sep="\n",file = infile)
	} else {
		cat(paste("   mol pqr ",apbsin$dir,"/",apbsin$center,sep=""),sep="\n",file = infile)
	}
	cat("end",sep="\n\n",file = infile)
		
	for(state in names(apbsin$states)){
		apbsin$states[[state]][grep("dime",apbsin$states[[state]])] <- paste("dime",paste(apbsin$dime,collapse=" "),sep=" ")
		apbsin$states[[state]][grep("cglen",apbsin$states[[state]])] <- paste("cglen",paste(apbsin$cglen,collapse=" "),sep=" ")
		apbsin$states[[state]][grep("fglen",apbsin$states[[state]])] <- paste("fglen",paste(apbsin$fglen,collapse=" "),sep=" ")
		cat(paste("elec name",state,sep=" "),sep="\n",file = infile)
		cat(apbsin$states[[state]],sep="\n   ",file = infile)
		
		if(length(apbsin$pot) > 0 & pot & state=="solv"){
			cat(paste("   write pot",apbsin$pot[1],apbsin$pot[2]),sep="\n", file = infile)
		}
		cat("end",sep="\n\n",file = infile)
	}
		
	for(state in names(apbsin$states)){
		if(length(grep("calcenergy total",apbsin$states[[state]])) > 0){
			cat(paste("print elecEnergy",state,"end",sep=" "),sep="\n", file = infile)
		}
	}
	cat("quit",sep="\n\n",file = infile)
	close(infile)
}

calc.pot <- function(apbsin,pot.dir){
	# calc.pot calculates spatial distributions of electostatic poential using APBS. The function outputs
	# potential files to pot.dir for every structure (PQR file) in apbsin$dir. When generating apbsin it 
	# is useful to use the s.only and pot options.   
	
	## Input Details ##
	#-----------------#
	# apbsin      - APBS input object      - APBS input object (apbsin.solv)
	# pot.dir     - character              - Output path for potential files

	files=list.files(path=apbsin$dir,pattern=".pqr")
	if (!file.exists(path = pot.dir)) {
		dir.create(pot.dir)
	}

  	pot.list <- NULL
  	for (file in files)	{
		apbsin$mols <- file
		protein.name <- strsplit(file,split=".pqr")[1]
		apbsin$pot["file"] <- paste(pot.dir,"/",protein.name,sep="")
		
		write.apbsin(apbsin,f.name="apbs_pot.in",pot=TRUE)
		system(paste(paths$apbs,"bin/apbs apbs_pot.in",sep=""))
		pot.list <- c(pot.list,apbsin$pot["file"])
		
	}
	return(list(dir=pot.dir,list=pot.list))
}

calc.dGsolv <- function(apbsin,coul.ep,batch.dir=NULL,pot.dir=NULL){
	# calc.dGsolv calculates solvation free energy using APBS. The function outputs solvation (dGsolv),
	# solution (Gsolu), reference (Gref), and Coulombic (Gcoul) free energies according to the parameters
	# of apbsin. If batch.dir is supplied solvation calulations are performed for every structure (PQR) 
	# in the directory. Additionally, if pot.dir is supplied electrostatic potential files for the solvated
	# state are outputted to pot.dir for every structure (PQR file) in apbsin$dir.  
	
	## Input Details ##
	#-----------------#
	# apbsin      - APBS input object      - APBS input object (apbsin.solv)
	# coul.ep	  - numeric				   - Dielectric coeficient for Coulomb's law calclations
	# batch.dir   - character	   		   - Directory of PQR files for batch operation [default = NULL] 
	# pot.dir     - character              - Output path for potential files [default = NULL]

	if(!is.null(batch.dir)){	
		
		files=list.files(path=batch.dir)
		if(!is.null(pot.dir)){
			if (!file.exists(path = pot.dir)){
				dir.create(pot.dir)
			}
		}
		
		apbsin$dir <- batch.dir
		FE <- NULL
		for(file in files){
			apbsin$mols <- file
			protein.name <- strsplit(file,split=".pqr")[1]
				
			if(!is.null(pot.dir)){
				apbsin$pot["file"] <- paste(pot.dir,"/",protein.name,sep="")
			}
			write.apbsin(apbsin,f.name="apbs_solv.in", pot= !is.null(pot.dir))
				
			s.fe.lines <- system(paste(paths$apbs,"bin/apbs apbs_solv.in",sep=""),intern=T)
			## Modified Coulomb path (10/14/13) -- Ron
			cou.fe.lines <- system(paste(paste(paths$coulomb,"coulomb ",sep=""),apbsin$dir,"/",apbsin$mols,sep=""),intern=T)
			s.fe <- as.numeric(unlist(strsplit(grep("Global\ net\ ELEC energy",s.fe.lines,value=T),split=" "))[c(8,17)])
			cou.fe <- as.numeric(unlist(strsplit(grep("Total\ energy\ =",cou.fe.lines,value=T),split=" "))[4])/coul.ep
			## Modified dGsolv and Gsolu and Gref indices (10/16/13) -- Ron
			FE <- rbind(FE,data.frame(dGsolv = s.fe[1] - s.fe[2],Gsolu = s.fe[1], Gref = s.fe[2],Gcoul = cou.fe,row.names=protein.name))
		}
				
	}else {
		if(!is.null(pot.dir)){
			protein.name <- strsplit(file,split=".pqr")[1]
			apbsin$pot["file"] <- paste(pot.dir,"/",protein.name,sep="")
		}
		write.apbsin(apbsin,f.name="apbs_solv.in", pot= !is.null(pot.dir))
		s.fe.lines <- system(paste(paths$apbs,"bin/apbs apbs_solv.in",sep=""),intern=T)
		## Modified Coulomb path (10/14/13) -- Ron
		cou.fe.lines <- system(paste(paste(paths$coulomb,"coulomb ",sep=""),apbsin$dir,"/",apbsin$mols,sep=""),intern=T)
		s.fe <-as.numeric(unlist(strsplit(grep("Global\ net\ ELEC energy",s.fe.lines,value=T),split=" "))[c(8,17)])
		cou.fe <-as.numeric(unlist(strsplit(grep("Total\ energy\ =",cou.fe.lines,value=T),split=" "))[4])/coul.ep
		
		## Modified dGsolv and Gsolu and Gref indices (10/16/13) -- Ron
		FE <- data.frame(dGsolv = s.fe[1] - s.fe[2],Gsolu = s.fe[1],Gref = s.fe[2],Gcoul = cou.fe)
	}
	return(FE)
}

calc.ddGsolv <- function(ab.solv.fe,a.solv.fe,b.solv.fe){
	# calc.ddGsolv calculates solvation free energies of association based on solvation energies for 
	# the complex and two individual components, as calculated by calc.dGsolv.   
	
	## Input Details ##
	#-----------------#
	# ab.solv.fe  - FE table               - Solvation free energies for complex (calc.dGsolv) 
	# a.solv.fe   - FE table               - Solvation free energies for component A (calc.dGsolv) 
	# b.solv.fe   - FE table               - Solvation free energies for component B (calc.dGsolv) 

	tmp <- sapply(1:length(row.names(a.solv.fe)), function(i) unlist(strsplit(row.names(a.solv.fe)[i],split="_")))
	a.comp.names <- sapply(1:length(row.names(a.solv.fe)), function(i) paste(tmp[[i]][1],chain.names[2],tmp[[i]][2],tmp[[i]][3],sep="_"))
	a.comp.names <- unlist(strsplit(a.comp.names,split="_NA_NA"))
	b.comp.names <- paste(chain.names[1],row.names(b.solv.fe),sep="_")
	ab.a.fe <- ab.solv.fe[row.names(ab.solv.fe) %in% a.comp.names,]
	ab.b.fe <- ab.solv.fe[row.names(ab.solv.fe) %in% b.comp.names,]
	
	a.ddGsolv <- ab.a.fe$dGsolv - a.solv.fe$dGsolv - b.solv.fe[chain.names[2],"dGsolv"]
	names(a.ddGsolv) <- row.names(a.solv.fe) 
	b.ddGsolv <- ab.b.fe$dGsolv - a.solv.fe[chain.names[1],"dGsolv"] - b.solv.fe$dGsolv
	names(b.ddGsolv) <- row.names(b.solv.fe) 
	
	return(list(a = a.ddGsolv, b = b.ddGsolv))
}

calc.assocFE <- function(ab.solv.fe,a.solv.fe,b.solv.fe){
	# calc.assocFE calculates free energies of association based on energies for the complex 
	# and two individual components, as calculated by calc.dGsolv. A list object containing 
	# a table of asocciation free energies for mutations of each chain. Calculates solvation 
	# free energy of association, solution association free energy, reference  assciations free 
	# energy, and Coulombic assocation free energy. 
	
	## Input Details ##
	#-----------------#
	# ab.solv.fe  - FE table               - Solvation free energies for complex (calc.dGsolv) 
	# a.solv.fe   - FE table               - Solvation free energies for component A (calc.dGsolv) 
	# b.solv.fe   - FE table               - Solvation free energies for component B (calc.dGsolv) 

	tmp <- sapply(1:length(row.names(a.solv.fe)), function(i) unlist(strsplit(row.names(a.solv.fe)[i],split="_")))
	a.comp.names <- sapply(1:length(row.names(a.solv.fe)), function(i) paste(tmp[[i]][1],chain.names[2],tmp[[i]][2],tmp[[i]][3],sep="_"))
	a.comp.names <- unlist(strsplit(a.comp.names,split="_NA_NA"))
	b.comp.names <- paste(chain.names[1],row.names(b.solv.fe),sep="_")
	ab.a.fe <- ab.solv.fe[row.names(ab.solv.fe) %in% a.comp.names,]
	ab.b.fe <- ab.solv.fe[row.names(ab.solv.fe) %in% b.comp.names,]
	
	a.ddGsolv <- NULL
	b.ddGsolv <- NULL
	for(i in 1:length(ab.solv.fe[1,])){
		a.ddGsolv <- cbind(a.ddGsolv,ab.a.fe[,i] - a.solv.fe[,i] - b.solv.fe[chain.names[2],i])
		b.ddGsolv <- cbind(b.ddGsolv,ab.b.fe[,i] - a.solv.fe[chain.names[1],i] - b.solv.fe[,i])
	}
	
	rownames(a.ddGsolv) <- rownames(a.solv.fe)
	rownames(b.ddGsolv) <- rownames(b.solv.fe) 
	colnames(a.ddGsolv) <- colnames(a.solv.fe)
	colnames(b.ddGsolv) <- colnames(b.solv.fe) 
	
	return(list(a = a.ddGsolv, b = b.ddGsolv))
}

esd.dist <- function(pot.dir,esd='awd'){
	# esd.dist generates a distance matrix of electrostatic dissimilarity values given a directory of 
	# spatial distributions of electrostatic potential (DX files). Four possibile distance measures can be 
	# used in generating the distance matrix see the following reference for more details: (AWD=LD; CBD=CDP)
	# Gorham, R., Kieslich, C. A., and Morikis, D. (2011) Electrostatic clustering and free energy calculations 
	# provide a foundation for protein design and optimization. Ann. Biomed Eng, 39(4): 1252-63. 
	
	## Input Details ##
	#-----------------#
	# pot.dir     - character              - Input path for directory of potential files	
	# esd         - character              - Electrostatic distance measure to be used [default = 'awd']
		  
	files <- list.files(path = pot.dir)
	header <- readLines(con=paste(pot.dir,"/",files[1],sep = ""),n=11)
	grid.dim <- as.numeric(strsplit(strsplit(header[grep("object 1",header)], split="counts ")[[1]][2],split=" ")[[1]])
	num.grid.pts <- grid.dim[1]*grid.dim[2]*grid.dim[3]
	
	dist <- matrix(0,length(files),length(files),dimnames=list(strsplit(files,".pqr.dx"),strsplit(files,".pqr.dx")));
	
	if(esd == 'awd'){
		esd.str <- 'sum(abs((potA+1e-9)-(potB+1e-9))/pmax(abs(potA+1e-9),abs(potB+1e-9)))/num.grid.pts'
	}
	else if(esd == 'cbd'){
		esd.str <- 'sqrt(1 - sum((potA+1e-9)*(potB+1e-9))/(sqrt(sum((potA+1e-9)^2))*sqrt(sum((potB+1e-9)^2))))'
	}
	else if(esd == 'dp'){
		esd.str <- 'sqrt(1 - 2*sum((potA+1e-9)*(potB+1e-9))/(sum((potA+1e-9)^2)+sum((potB+1e-9)^2)))'
	}
	else if(esd == 'ldp'){
		esd.str <- 'sqrt(1 - sum(2*(potA+1e-9)*(potB+1e-9)/((potA+1e-9)^2+(potB+1e-9)^2)))'
	}
 
	for (i in 1:(length(files)-1)) {
	
		# Opens DX file
		potA <- read.delim(paste(pot.dir,"/",files[i],sep = ""), skip=11, header=F, sep="", nrows=num.grid.pts/3, colClasses = c("numeric","numeric","numeric")); 
	 	potA <- as.matrix(potA);
		
		for (k in (i+1):length(files)) {
			potB <- read.delim(paste(pot.dir,"/",files[k],sep = ""), skip=11, header=F, sep="", nrows=num.grid.pts/3, colClasses = c("numeric","numeric","numeric")); 
	        	potB <- as.matrix(potB);
			
			dist[i,k]=dist[k,i] <- eval(parse(text=esd.str))
		}
	}

	return(dist)
}

esi.distr <- function(pot.dir,parent,outname="esi_distr.dx",esi="wd"){
	# esd.distr calculates the spatial distribution of electrostatic similarity for the comparison of a set 
	# of potenitals versus a reference potential (parent). esd.distr takes as input a directory of 
	# spatial distributions of electrostatic potential (DX files). See the following reference for more details: 
	# Kieslich, C.A. and Morikis, D. (2012) The two sides of complement C3d: Evolution of electrostatics in 
	# a link between innate and adaptive immunity. PLoS Comp. Bio

	
	## Input Details ##
	#-----------------#
	# pot.dir     - character              - Input path for directory of potential files	
	# parent      - character              - Name of reference potential (parent)
	# outname     - character              - Output name for generated DX file [default = 'esi_distr.dx']
	# esi         - character              - Electrostatic similarity measure to be used [default = 'wd']
	
	
	files <- list.files(path = pot.dir)
	mut.files <- files[-grep(parent,files)]
	
	header <- readLines(con=paste(pot.dir,"/",files[grep(parent,files)],sep = ""),n=11)
	grid.dim <- as.numeric(strsplit(strsplit(header[grep("object 1",header)], split="counts ")[[1]][2],split=" ")[[1]])
	num.grid.pts <- grid.dim[1]*grid.dim[2]*grid.dim[3] 
	
	if(esi == 'wd'){
		esi.str <- '1-abs((potA+1e-9)-(potB+1e-9))/pmax(abs(potA+1e-9),abs(potB+1e-9))'
	}
	else if(esi == 'dp'){
		esi.str <- '(potA+1e-9)*(potB+1e-9)/((potA+1e-9)^2+(potB+1e-9)^2)'
	}
 
	# Opens parent DX file
	potA <- read.delim(paste(pot.dir,"/",files[grep(parent,files)],sep = ""), skip=11, header=F, sep="", nrows=num.grid.pts/3, colClasses = c("numeric","numeric","numeric")); 
	potA <- as.matrix(potA);
	
	for (k in 1:length(mut.files)) {
		print(paste(k,mut.files[k]))
		potB <- read.delim(paste(pot.dir,"/",mut.files[k],sep = ""), skip=11, header=F, sep="", nrows=num.grid.pts/3, colClasses = c("numeric","numeric","numeric")); 
	    potB <- as.matrix(potB);
		if(k == 1){
			distr <- eval(parse(text=esi.str))
		}
		else {
			distr <- distr + eval(parse(text=esi.str))
		}
	}
	distr <- distr/length(mut.files)
	
	writeLines(header,con=outname)
	f <- file(outname,"a")
	writeLines(sprintf("%7.6e %7.6e %7.6e ",distr[,1],distr[,2],distr[,3]),con=f)
	close(f)
}

mut.aa <- function(tmp.seq,muts,out.pdb,pdb.file=NULL,write.pdb=F, input.chain = inputchain){ 
	# mut.aa is a wrapper function for Dunbrack's Scwrl4 package and is used to perform non-alanine mutations.
	# The function takes as an input the template sequence of the structure to be mutated, along with a data.frame
	# (muts) that contains a column of residue numbers ($res) and a column of target amino acids ($aa).   
	
	#Added changes that bypass scwrl4 and use Chimera instead which uses Dunbrack's rotamer libraries. This change was made due to a couple reasons: 1) zinc (or other metal ions do not work with Scwrl and 2) Scwrl seems to change the atomic coordinates drastically in some cases RM 03/13/14
	
	## Input Details ##
	#-----------------#
	# tmp.seq     - character array        - Array containing single letter amino acid sequence of template
	# muts        - data.frame             - Table of residue numbers and target amino acids to be mutated
	# out.pdb     - character              - Output name for mutant PDB
	# pdb.file    - character              - PDB filename contianing template structure [default = NULL]	
	# write.pdb   - True/False             - Logical indicating whether PDB file is needed to be written [default = F]	
	cat("open scwrl_in.pdb",file="chimera_command.com",sep="\n") #creation of the chimera command script RM 03/13/14
	#lc.seq <- capply(tmp.seq,lower) #this doesn't do anything for Chimera RM 03/13/14
	#lc.seq[muts$res] <- as.character(muts$aa) #this doesn't do anything for Chimera RM 03/13/14
	cat(paste("swapaa ", tolower(aa123(as.character(muts$aa))), " #0:",muts$res,".",input.chain,sep=""),file="chimera_command.com",sep="\n",append=TRUE)
	#Needed to add conditions to account for double/triple simultaneous mutations - RM
	### 04/07/2014 Edits to conditional statement to look for residue number gap (RH)
	#This checks to see that the 2nd res column is not NA - RM 03/04/14
	if(is.na(muts$res2) == FALSE) {
		#lc.seq[muts$res2] <- as.character(muts$aa2) #this doesn't do anything for Chimera RM 03/13/14
		cat(paste("swapaa ", tolower(aa123(as.character(muts$aa2))), " #0:",muts$res2,".",input.chain,sep=""),file="chimera_command.com",sep="\n",append=TRUE)
	}
	#This checks to see that 3rd res column is not NA - RM 03/04/14
	if(is.na(muts$res3) == FALSE) {
		#lc.seq[muts$res3] <- as.character(muts$aa3)#this doesn't do anything for Chimera RM 03/13/14
		cat(paste("swapaa ", tolower(aa123(as.character(muts$aa3))), " #0:",muts$res3,".",input.chain,sep=""),file="chimera_command.com",sep="\n",append=TRUE)
	}
	#cat(lc.seq,sep="",file="mutseq.txt")#this doesn't do anything for Chimera RM 03/13/14
	
	if(write.pdb){
		pdb <- read.pdb(pdb.file)
		write.pdb(pdb,file="scwrl_in.pdb")
	}
	
	#system(paste(paths$scwrl,"Scwrl4 -i scwrl_in.pdb -o scwrl_out.pdb -s mutseq.txt",sep=""))  #This calls scwrl, uncomment if you wish to use scwrl
	cat("write format pdb 0 scwrl_out.pdb",file="chimera_command.com",sep="\n",append=TRUE) 
	cat("stop",file="chimera_command.com",append=TRUE) #exits chimera RM 03/13/14
	system(paste(paths$chimera,"chimera --nogui chimera_command.com",sep="")) #Runs chimera in nogui/commandline mode and runs the command script RM 03/13/14 | added chimera path RM 10/14/14
	rewrite<-read.pdb("scwrl_out.pdb", het2atom = TRUE) #...
	write.pdb(rewrite, file="scwrl_out.pdb") #Need to rewrite the PDB to account for Chimera adding HETATM instead of ATOM RM 03/20/14
	system(paste("cp","scwrl_out.pdb",out.pdb))
}
#Edited this to work with changes made to mut.aa RM 10/14/14
char.rev.c <- function(pdb,chain.names,chains=unique(pdb$atom[,"chain"]),ff="PARSE"){
	# char.rev.c performs a charge reversal scan of charged residues (D, E, K, R, and H), by replacing 
	# acidic residues with lysine and basic residues with glutamic acid. mut.aa is then utilized to 
	# call Swrl4 to generate mutant PDB files. Subsequently the mutant PDBs are used to generate PQR
	# files using batch.pqr() and PDB2PQR. This version is intended for a protein complex.  
	
	## Input Details ##
	#-----------------#
	# pdb         - PDB object (Bio3D) - Parent PDB object as generated by read.pdb
	# chain.names - character	       - Array of protein chain names to be used in output names
	# chains      - character	       - Array of chain IDs to be mutated [default = all included chains]
	# ff          - character          - PDB2PQR force field [default = PARSE]
	inputchain <- chains
	protein.name <- paste(chain.names,collapse="_")
	dirs <- c(paste(getwd(),"/",paste(protein.name,"pdb",sep="_"),sep=""),
			  paste(getwd(),"/",paste(chain.names,"pdb",sep="_"),sep=""))
	
	for(dir in dirs){
		if(file.exists(path = dir) == 0){
			dir.create(dir)
		}
	}
	
	# Write parent files
	chainA <- trim.pdb(pdb,inds=atom.select(pdb,chain=unique(pdb$atom[,"chain"])[1],verbose=F))
	chainB <- trim.pdb(pdb,inds=atom.select(pdb,chain=unique(pdb$atom[,"chain"])[2],verbose=F))
	write.pdb(pdb,file=paste(dirs[1],"/",protein.name,".pdb",sep=""))	
	write.pdb(chainA,file=paste(dirs[2],"/",chain.names[1],".pdb",sep=""))
	write.pdb(chainB,file=paste(dirs[3],"/",chain.names[2],".pdb",sep=""))
	
	# Write tmp file
	write.pdb(pdb,file="scwrl_in.pdb")
	
	# Names of possible titratable groups (only standard groups available currently)       
	ion.name <- c("E","D","K","R")
		
	# Extract amino acid sequence, residue numbers, and chain ids
	aaseq <- pdbseq(pdb) #seq.pdb renamed to pdbseq in bio3d Rohith
	resnos <- as.numeric(pdb$atom[pdb$calpha,"resno"])
	chain <- pdb$atom[pdb$calpha,"chain"]
	seqnos <- 1:length(aaseq)
	
	# Determine which residue out of the entire sequence to mutate
	ion.nos <- resnos[aaseq %in% ion.name & chain %in% chains]
	ion.ids <- aaseq[aaseq %in% ion.name & chain %in% chains]
	ion.chain <- chain[aaseq %in% ion.name & chain %in% chains]
	ion.seq <- seqnos [aaseq %in% ion.name & chain %in% chains]
	originaldf<-data.frame(ion.nos,ion.ids,muts=ifelse(ion.ids == "E" | ion.ids == "D","K","E"))
	muts <- data.frame(res=originaldf$ion.nos,aa=originaldf$muts,res2=NA,aa2=NA,res3=NA,aa3=NA)
	mut.list <- lapply(1:dim(muts)[1],function(x) muts[x,])
	mutlist<- mut.list
	for(i in 1:length(mutlist)){
		#Needed to add conditions to account for double/triple simultaneous mutations - RM 03/04/14
		#hardcoded chain input RM 03/17/14		
		#This is the original condition i.e. Just a single mutation. This checks to see that the 2nd res column is NA (which means the 3rd column will also be NA) - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == TRUE) {
			out_name <- paste(dirs[1],"/",paste(protein.name,paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""),collapse="_"),sep="_"),".pdb",sep="")
		}
		mut.aa(tmp.seq = aaseq,muts = mutlist[[i]],out.pdb = out_name, input.chain = inputchain)
		
		mut.pdb <- read.pdb(out_name, het2atom = TRUE) #het2atom = T This change should be okay as the initial read.pdb in the mut_list script determines whether the PDB will be cleaned up or not by the time it gets to here RM 03/13/14
		mut.chain <- trim.pdb(mut.pdb,inds=atom.select(mut.pdb,chain=toupper(inputchain),verbose=F)) #hardcoded chain input RM 03/17/14
		#Needed to add conditions to account for double/triple simultaneous mutations - RM 03/04/14
		
		#This is the original condition i.e. Just a single mutation. This checks to see that the 2nd res column is NA (which means the 3rd column will also be NA) - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == TRUE) {
			out_name <- paste(dirs[1+which(unique(pdb$atom[,"chain"])==toupper(inputchain))],"/",paste(chain.names[which(unique(pdb$atom[,"chain"])==toupper(inputchain))],paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""),collapse="_"),sep="_"),".pdb",sep="")
		}
		
		write.pdb(mut.chain,file=out_name)
		write.pdb(mut.pdb, file="scwrl_out.pdb") #Need to rewrite the PDB to account for Chimera adding HETATM instead of ATOM RM 03/20/14
	}
	pqr.dirs <- c(paste(getwd(),"/",paste(protein.name,"pqr",sep="_"),sep=""),
			      paste(getwd(),"/",paste(chain.names,"pqr",sep="_"),sep=""))
	
	for(dir in pqr.dirs){
		if(file.exists(path = dir) == 0){
			dir.create(dir)
		}
	}
	batch.pqr(dirs[1],pqr.dirs[1],ff=ff,flags=c("--nodebump","--noopt"))
	batch.pqr(dirs[2],pqr.dirs[2],ff=ff,flags=c("--nodebump","--noopt"))
	batch.pqr(dirs[3],pqr.dirs[3],ff=ff,flags=c("--nodebump","--noopt"))
	
	return(list(dirs = pqr.dirs))
}

mut.list.c <- function(pdb,chain.names,mutlist,ff="PARSE", inputchain){
	# mut.list.c generates mutants specified by mutlist. mutlist is list object, where each element 
	# of the list is a data.frame containing res and aa columns (as is used for input for mut.aa). 
	# mut.aa is utilized to call Swrl4 to generate mutant PDB files. Subsequently the mutant PDBs 
	# are used to generate PQR files using batch.pqr() and PDB2PQR. This version is intended for a 
	# protein complex.  
	
	## Input Details ##
	#-----------------#
	# pdb         - PDB object (Bio3D) - Parent PDB object as generated by read.pdb
	# chain.names - character	       - Array of protein chain names to be used in output names
	# mutlist     - list               - List of data.frames with residue numbers and target amino acids for each mutant
	# ff          - character          - PDB2PQR force field [default = PARSE]
	
	protein.name <- paste(chain.names,collapse="_")
	dirs <- c(paste(getwd(),"/",paste(protein.name,"pdb",sep="_"),sep=""),
			  paste(getwd(),"/",paste(chain.names,"pdb",sep="_"),sep=""))
	
	for(dir in dirs){
		if(file.exists(path = dir) == 0){
			dir.create(dir)
		}
	}
	
	# Write parent files
	chainA <- trim.pdb(pdb,inds=atom.select(pdb,chain=unique(pdb$atom[,"chain"])[1],verbose=F))
	chainB <- trim.pdb(pdb,inds=atom.select(pdb,chain=unique(pdb$atom[,"chain"])[2],verbose=F))
	write.pdb(pdb,file=paste(dirs[1],"/",protein.name,".pdb",sep=""))	
	write.pdb(chainA,file=paste(dirs[2],"/",chain.names[1],".pdb",sep=""))
	write.pdb(chainB,file=paste(dirs[3],"/",chain.names[2],".pdb",sep=""))
	
	# Write tmp file
	write.pdb(pdb,file="scwrl_in.pdb")
	
	aaseq <- pdbseq(pdb) #seq.pdb renamed to pdbseq in bio3d Rohith
	chain <- pdb$atom[pdb$calpha,"chain"]
	resnos <- as.numeric(pdb$atom[pdb$calpha,"resno"])
	
	for(i in 1:length(mutlist)){
		#Needed to add conditions to account for double/triple simultaneous mutations - RM 03/04/14
		#hardcoded chain input RM 03/17/14		
		#This is the original condition i.e. Just a single mutation. This checks to see that the 2nd res column is NA (which means the 3rd column will also be NA) - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == TRUE) {
			out_name <- paste(dirs[1],"/",paste(protein.name,paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""),collapse="_"),sep="_"),".pdb",sep="")
		}
		#This checks to see that the 2nd res column is not NA  and the 3rd column is NA aka a double mutation - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == FALSE & is.na(mutlist[[i]]["res3"]) == TRUE ) {
			out_name <- paste(dirs[1],"/",paste(protein.name,paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""), paste(toupper(inputchain),mutlist[[i]]$res2,mutlist[[i]]$aa2,sep=""),sep="_"),sep="_"),".pdb",sep="")
		}
		#This checks to see that both the 2nd and 3rd res column are not NA aka a triple mutation - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == FALSE & is.na(mutlist[[i]]["res3"]) == FALSE ) {
			out_name <- paste(dirs[1],"/",paste(protein.name,paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""), paste(toupper(inputchain),mutlist[[i]]$res2,mutlist[[i]]$aa2,sep=""), paste(toupper(inputchain),mutlist[[i]]$res3,mutlist[[i]]$aa3,sep=""),sep="_"),sep="_"),".pdb",sep="")
		}
		mut.aa(tmp.seq = aaseq,muts = mutlist[[i]],out.pdb = out_name, input.chain = inputchain)
		
		mut.pdb <- read.pdb(out_name, het2atom = TRUE) #het2atom = T This change should be okay as the initial read.pdb in the mut_list script determines whether the PDB will be cleaned up or not by the time it gets to here RM 03/13/14
		mut.chain <- trim.pdb(mut.pdb,inds=atom.select(mut.pdb,chain=toupper(inputchain),verbose=F)) #hardcoded chain input RM 03/17/14
		#Needed to add conditions to account for double/triple simultaneous mutations - RM 03/04/14
		
		#This is the original condition i.e. Just a single mutation. This checks to see that the 2nd res column is NA (which means the 3rd column will also be NA) - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == TRUE) {
			out_name <- paste(dirs[1+which(unique(pdb$atom[,"chain"])==toupper(inputchain))],"/",paste(chain.names[which(unique(pdb$atom[,"chain"])==toupper(inputchain))],paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""),collapse="_"),sep="_"),".pdb",sep="")
		}
		#This checks to see that the 2nd res column is not NA  and the 3rd column is NA aka a double mutation - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == FALSE & is.na(mutlist[[i]]["res3"]) == TRUE ) {
			out_name <- paste(dirs[1+which(unique(pdb$atom[,"chain"])==toupper(inputchain))],"/",paste(chain.names[which(unique(pdb$atom[,"chain"])==toupper(inputchain))],paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""), paste(toupper(inputchain),mutlist[[i]]$res2,mutlist[[i]]$aa2,sep=""),sep="_"),sep="_"),".pdb",sep="")
		}
		#This checks to see that both the 2nd and 3rd res column are not NA aka a triple mutation - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == FALSE & is.na(mutlist[[i]]["res3"]) == FALSE ) {
			out_name <- paste(dirs[1+which(unique(pdb$atom[,"chain"])==toupper(inputchain))],"/",paste(chain.names[which(unique(pdb$atom[,"chain"])==toupper(inputchain))],paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""), paste(toupper(inputchain),mutlist[[i]]$res2,mutlist[[i]]$aa2,sep=""), paste(toupper(inputchain),mutlist[[i]]$res3,mutlist[[i]]$aa3,sep=""),sep="_"),sep="_"),".pdb",sep="")
		}
		
		write.pdb(mut.chain,file=out_name)
		write.pdb(mut.pdb, file="scwrl_out.pdb") #Need to rewrite the PDB to account for Chimera adding HETATM instead of ATOM RM 03/20/14
	}
	
	pqr.dirs <- c(paste(getwd(),"/",paste(protein.name,"pqr",sep="_"),sep=""),
			      paste(getwd(),"/",paste(chain.names,"pqr",sep="_"),sep=""))
	
	for(dir in pqr.dirs){
		if(file.exists(path = dir) == 0){
			dir.create(dir)
		}
	}
	
	batch.pqr(dirs[1],pqr.dirs[1],ff=ff,flags=c("--nodebump","--noopt"))
	batch.pqr(dirs[2],pqr.dirs[2],ff=ff,flags=c("--nodebump","--noopt"))
	batch.pqr(dirs[3],pqr.dirs[3],ff=ff,flags=c("--nodebump","--noopt"))
	
	return(list(dirs = pqr.dirs))
}

mut.comb.c <- function(pdb,chain.names,muts,mpa,ff="PARSE",inputchain){
	# mut.comb.c generates all possible combinations of mutants given the set of possible mutations,
	# defined in muts (as is used for input for mut.aa), and the specified number of mutations per
	# analog (mpa). mut.aa is utilized to call Swrl4 to generate mutant PDB files. Subsequently the mutant 
	# PDBs are used to generate PQR files using batch.pqr() and PDB2PQR. This version is intended 
	# for a protein complex.  
	
	## Input Details ##
	#-----------------#
	# pdb         - PDB object (Bio3D) - Parent PDB object as generated by read.pdb
	# chain.names - character	       - Array of protein chain names to be used in output names
	# muts        - data.frame         - Table of residue numbers and target amino acids to be considered
	# mpa         - integer            - Number of mutations per analog
	# ff          - character          - PDB2PQR force field [default = PARSE]
	protein.name <- paste(chain.names,collapse="_")
	dirs <- c(paste(getwd(),"/",paste(protein.name,"pdb",sep="_"),sep=""),
			  paste(getwd(),"/",paste(chain.names,"pdb",sep="_"),sep=""))
	
	for(dir in dirs){
		if(file.exists(path = dir) == 0){
			dir.create(dir)
		}
	}
	
	# Write parent files
	chainA <- trim.pdb(pdb,inds=atom.select(pdb,chain=unique(pdb$atom[,"chain"])[1],verbose=F))
	chainB <- trim.pdb(pdb,inds=atom.select(pdb,chain=unique(pdb$atom[,"chain"])[2],verbose=F))
	write.pdb(pdb,file=paste(dirs[1],"/",protein.name,".pdb",sep=""))	
	write.pdb(chainA,file=paste(dirs[2],"/",chain.names[1],".pdb",sep=""))
	write.pdb(chainB,file=paste(dirs[3],"/",chain.names[2],".pdb",sep=""))
	
	# Write tmp file
	write.pdb(pdb,file="scwrl_in.pdb")
	
	aaseq <- pdbseq(pdb) #seq.pdb renamed to pdbseq in bio3d Rohith
	chain <- pdb$atom[pdb$calpha,"chain"]
	resnos <- as.numeric(pdb$atom[pdb$calpha,"resno"])
	
	comb<-combn(muts$res,3)  #generates combinations of 3 RM 10/14/14
	mut.list<-merge(merge(merge(data.frame(res=comb[1,],res2=comb[2,],res3=comb[3,]),muts,by=c("res")),muts,by.x=c("res2"),by.y=c("res"),suffixes=c("","2")),muts,by.x=c("res3"),by.y=c("res"),suffixes=c("","3"))
	mutlist <- lapply(1:dim(mut.list)[1],function(x) mut.list[x,])
	for(i in 1:length(mutlist)){
		#Needed to add conditions to account for double/triple simultaneous mutations - RM 03/04/14
		#hardcoded chain input RM 03/17/14		
		#This is the original condition i.e. Just a single mutation. This checks to see that the 2nd res column is NA (which means the 3rd column will also be NA) - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == TRUE) {
			out_name <- paste(dirs[1],"/",paste(protein.name,paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""),collapse="_"),sep="_"),".pdb",sep="")
		}
		#This checks to see that the 2nd res column is not NA  and the 3rd column is NA aka a double mutation - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == FALSE & is.na(mutlist[[i]]["res3"]) == TRUE ) {
			out_name <- paste(dirs[1],"/",paste(protein.name,paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""), paste(toupper(inputchain),mutlist[[i]]$res2,mutlist[[i]]$aa2,sep=""),sep="_"),sep="_"),".pdb",sep="")
		}
		#This checks to see that both the 2nd and 3rd res column are not NA aka a triple mutation - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == FALSE & is.na(mutlist[[i]]["res3"]) == FALSE ) {
			out_name <- paste(dirs[1],"/",paste(protein.name,paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""), paste(toupper(inputchain),mutlist[[i]]$res2,mutlist[[i]]$aa2,sep=""), paste(toupper(inputchain),mutlist[[i]]$res3,mutlist[[i]]$aa3,sep=""),sep="_"),sep="_"),".pdb",sep="")
		}
		mut.aa(tmp.seq = aaseq,muts = mutlist[[i]],out.pdb = out_name, input.chain = inputchain)
		
		mut.pdb <- read.pdb(out_name, het2atom = TRUE) #het2atom = T This change should be okay as the initial read.pdb in the mut_list script determines whether the PDB will be cleaned up or not by the time it gets to here RM 03/13/14
		mut.chain <- trim.pdb(mut.pdb,inds=atom.select(mut.pdb,chain=toupper(inputchain),verbose=F)) #hardcoded chain input RM 03/17/14
		#Needed to add conditions to account for double/triple simultaneous mutations - RM 03/04/14
		
		#This is the original condition i.e. Just a single mutation. This checks to see that the 2nd res column is NA (which means the 3rd column will also be NA) - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == TRUE) {
			out_name <- paste(dirs[1+which(unique(pdb$atom[,"chain"])==toupper(inputchain))],"/",paste(chain.names[which(unique(pdb$atom[,"chain"])==toupper(inputchain))],paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""),collapse="_"),sep="_"),".pdb",sep="")
		}
		#This checks to see that the 2nd res column is not NA  and the 3rd column is NA aka a double mutation - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == FALSE & is.na(mutlist[[i]]["res3"]) == TRUE ) {
			out_name <- paste(dirs[1+which(unique(pdb$atom[,"chain"])==toupper(inputchain))],"/",paste(chain.names[which(unique(pdb$atom[,"chain"])==toupper(inputchain))],paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""), paste(toupper(inputchain),mutlist[[i]]$res2,mutlist[[i]]$aa2,sep=""),sep="_"),sep="_"),".pdb",sep="")
		}
		#This checks to see that both the 2nd and 3rd res column are not NA aka a triple mutation - RM 03/04/14
		if(is.na(mutlist[[i]]["res2"]) == FALSE & is.na(mutlist[[i]]["res3"]) == FALSE ) {
			out_name <- paste(dirs[1+which(unique(pdb$atom[,"chain"])==toupper(inputchain))],"/",paste(chain.names[which(unique(pdb$atom[,"chain"])==toupper(inputchain))],paste(paste(toupper(inputchain),mutlist[[i]]$res,mutlist[[i]]$aa,sep=""), paste(toupper(inputchain),mutlist[[i]]$res2,mutlist[[i]]$aa2,sep=""), paste(toupper(inputchain),mutlist[[i]]$res3,mutlist[[i]]$aa3,sep=""),sep="_"),sep="_"),".pdb",sep="")
		}
		
		write.pdb(mut.chain,file=out_name)
		write.pdb(mut.pdb, file="scwrl_out.pdb") #Need to rewrite the PDB to account for Chimera adding HETATM instead of ATOM RM 03/20/14
	}
	
	pqr.dirs <- c(paste(getwd(),"/",paste(protein.name,"pqr",sep="_"),sep=""),
			      paste(getwd(),"/",paste(chain.names,"pqr",sep="_"),sep=""))
	
	for(dir in pqr.dirs){
		if(file.exists(path = dir) == 0){
			dir.create(dir)
		}
	}
	
	batch.pqr(dirs[1],pqr.dirs[1],ff=ff,flags=c("--nodebump","--noopt"))
	batch.pqr(dirs[2],pqr.dirs[2],ff=ff,flags=c("--nodebump","--noopt"))
	batch.pqr(dirs[3],pqr.dirs[3],ff=ff,flags=c("--nodebump","--noopt"))
	
	return(list(dirs = pqr.dirs))
}

gpl.license <- function(){

	print(c("","",
	"     GNU GENERAL PUBLIC LICENSE",
	"     Version 3, 29 June 2007",
	"     Preamble",
	"",
	"     The GNU General Public License is a free, copyleft license for software and other kinds of works.",
	"",
	"     The licenses for most software and other practical works are designed to take away your freedom to share ",
	"     and change the works. By contrast, the GNU General Public License is intended to guarantee your freedom ",
	"     to share and change all versions of a program--to make sure it remains free software for all its users. ",
	"     We, the Free Software Foundation, use the GNU General Public License for most of our software; it applies ",
	"     also to any other work released this way by its authors. You can apply it to your programs, too.",
	"",
	"     When we speak of free software, we are referring to freedom, not price. Our General Public Licenses are ",
	"     designed to make sure that you have the freedom to distribute copies of free software (and charge for them", 
	"     if you wish), that you receive source code or can get it if you want it, that you can change the software ",
	"     or use pieces of it in new free programs, and that you know you can do these things.",
	"",
	"     To protect your rights, we need to prevent others from denying you these rights or asking you to surrender ",
	"     the rights. Therefore, you have certain responsibilities if you distribute copies of the software, or if ",
	"     you modify it: responsibilities to respect the freedom of others.",
	"",
	"     For example, if you distribute copies of such a program, whether gratis or for a fee, you must pass on to ",
	"     the recipients the same freedoms that you received. You must make sure that they, too, receive or can get ",
	"     the source code. And you must show them these terms so they know their rights.",
	"",
	"     Developers that use the GNU GPL protect your rights with two steps: (1) assert copyright on the software, ",
	"     and (2) offer you this License giving you legal permission to copy, distribute and/or modify it.",
	"",
	"     For the developers' and authors' protection, the GPL clearly explains that there is no warranty for this ",
	"     free software. For both users' and authors' sake, the GPL requires that modified versions be marked as ",
	"     changed, so that their problems will not be attributed erroneously to authors of previous versions.",
	"",
	"     Some devices are designed to deny users access to install or run modified versions of the software inside ",
	"     them, although the manufacturer can do so. This is fundamentally incompatible with the aim of protecting ",
	"     users' freedom to change the software. The systematic pattern of such abuse occurs in the area of products ",
	"     for individuals to use, which is precisely where it is most unacceptable. Therefore, we have designed this ",
	"     version of the GPL to prohibit the practice for those products. If such problems arise substantially in other", 
	"     domains, we stand ready to extend this provision to those domains in future versions of the GPL, as needed ",
	"     to protect the freedom of users.",
	"",
	"     Finally, every program is threatened constantly by software patents. States should not allow patents to ",
	"     restrict development and use of software on general-purpose computers, but in those that do, we wish to ",
	"     avoid the special danger that patents applied to a free program could make it effectively proprietary. To ",
	"     prevent this, the GPL assures that patents cannot be used to render the program non-free.",
	"",
	"     The precise terms and conditions for copying, distribution and modification follow.",
	"",
	"",
	"     TERMS AND CONDITIONS",
	"     0. Definitions.",
	"     “This License” refers to version 3 of the GNU General Public License.",
	"",
	"     “Copyright” also means copyright-like laws that apply to other kinds of works, such as semiconductor masks.",
	"",
	"     “The Program” refers to any copyrightable work licensed under this License. Each licensee is addressed as ",
	"     “you”. “Licensees” and “recipients” may be individuals or organizations.",
	"",
	"     To “modify” a work means to copy from or adapt all or part of the work in a fashion requiring copyright ",
	"     permission, other than the making of an exact copy. The resulting work is called a “modified version” of ",
	"     the earlier work or a work “based on” the earlier work.",
	"",
	"     A “covered work” means either the unmodified Program or a work based on the Program.",
	"",
	"     To “propagate” a work means to do anything with it that, without permission, would make you directly or ",
	"     secondarily liable for infringement under applicable copyright law, except executing it on a computer or ",
	"     modifying a private copy. Propagation includes copying, distribution (with or without modification), ",
	"     making available to the public, and in some countries other activities as well.",
	"",
	"     To “convey” a work means any kind of propagation that enables other parties to make or receive copies. Mere ",
	"     interaction with a user through a computer network, with no transfer of a copy, is not conveying.",
	"",
	"     An interactive user interface displays “Appropriate Legal Notices” to the extent that it includes a convenient ",
	"     and prominently visible feature that (1) displays an appropriate copyright notice, and (2) tells the user that ",
	"     there is no warranty for the work (except to the extent that warranties are provided), that licensees may convey ",
	"     the work under this License, and how to view a copy of this License. If the interface presents a list of user ",
	"     commands or options, such as a menu, a prominent item in the list meets this criterion.",
	"",
	"     1. Source Code.",
	"     The “source code” for a work means the preferred form of the work for making modifications to it. “Object code” ",
	"     means any non-source form of a work.",
	"",
	"     A “Standard Interface” means an interface that either is an official standard defined by a recognized standards ",
	"     body, or, in the case of interfaces specified for a particular programming language, one that is widely used among", 
	"     developers working in that language.",
	"",
	"     The “System Libraries” of an executable work include anything, other than the work as a whole, that (a) is included ",
	"     in the normal form of packaging a Major Component, but which is not part of that Major Component, and (b) serves only ",
	"     to enable use of the work with that Major Component, or to implement a Standard Interface for which an implementation ",
	"     is available to the public in source code form. A “Major Component”, in this context, means a major essential component ",
	"     (kernel, window system, and so on) of the specific operating system (if any) on which the executable work runs, or a ",
	"     compiler used to produce the work, or an object code interpreter used to run it.",
	"",
	"     The “Corresponding Source” for a work in object code form means all the source code needed to generate, install, and ",
	"     (for an executable work) run the object code and to modify the work, including scripts to control those activities. ",
	"     However, it does not include the work's System Libraries, or general-purpose tools or generally available free programs", 
	"     which are used unmodified in performing those activities but which are not part of the work. For example, Corresponding ",
	"     Source includes interface definition files associated with source files for the work, and the source code for shared ",
	"     libraries and dynamically linked subprograms that the work is specifically designed to require, such as by intimate data", 
	"     communication or control flow between those subprograms and other parts of the work.",
	"",
	"     The Corresponding Source need not include anything that users can regenerate automatically from other parts of the ",
	"     Corresponding Source.",
	"",
	"     The Corresponding Source for a work in source code form is that same work.",
	"",
	"     2. Basic Permissions.",
	"     All rights granted under this License are granted for the term of copyright on the Program, and are irrevocable provided ",
	"     the stated conditions are met. This License explicitly affirms your unlimited permission to run the unmodified Program. ",
	"     The output from running a covered work is covered by this License only if the output, given its content, constitutes a ",
	"     covered work. This License acknowledges your rights of fair use or other equivalent, as provided by copyright law.",
	"",
	"     You may make, run and propagate covered works that you do not convey, without conditions so long as your license otherwise ",
	"     remains in force. You may convey covered works to others for the sole purpose of having them make modifications exclusively ",
	"     for you, or provide you with facilities for running those works, provided that you comply with the terms of this License in ",
	"     conveying all material for which you do not control copyright. Those thus making or running the covered works for you must ",
	"     do so exclusively on your behalf, under your direction and control, on terms that prohibit them from making any copies of ",
	"     your copyrighted material outside their relationship with you.",
	"",
	"     Conveying under any other circumstances is permitted solely under the conditions stated below. Sublicensing is not allowed; ",
	"     section 10 makes it unnecessary.",
	"",
	"     3. Protecting Users' Legal Rights From Anti-Circumvention Law.",
	"     No covered work shall be deemed part of an effective technological measure under any applicable law fulfilling obligations ",
	"     under article 11 of the WIPO copyright treaty adopted on 20 December 1996, or similar laws prohibiting or restricting ",
	"     circumvention of such measures.",
	"",
	"     When you convey a covered work, you waive any legal power to forbid circumvention of technological measures to the extent ",
	"     such circumvention is effected by exercising rights under this License with respect to the covered work, and you disclaim ",
	"     any intention to limit operation or modification of the work as a means of enforcing, against the work's users, your or ",
	"     third parties' legal rights to forbid circumvention of technological measures.",
	"",
	"     4. Conveying Verbatim Copies.",
	"     You may convey verbatim copies of the Program's source code as you receive it, in any medium, provided that you conspicuously ",
	"     and appropriately publish on each copy an appropriate copyright notice; keep intact all notices stating that this License ",
	"     and any non-permissive terms added in accord with section 7 apply to the code; keep intact all notices of the absence of any", 
	"     warranty; and give all recipients a copy of this License along with the Program.",
	"",
	"     You may charge any price or no price for each copy that you convey, and you may offer support or warranty protection for ",
	"     a fee.",
	"",
	"     5. Conveying Modified Source Versions.",
	"     You may convey a work based on the Program, or the modifications to produce it from the Program, in the form of source code ",
	"     under the terms of section 4, provided that you also meet all of these conditions:",
	"",
	"     a) The work must carry prominent notices stating that you modified it, and giving a relevant date.",
	"     b) The work must carry prominent notices stating that it is released under this License and any conditions added under ",
	"     section 7. This requirement modifies the requirement in section 4 to “keep intact all notices”.",
	"     c) You must license the entire work, as a whole, under this License to anyone who comes into possession of a copy. This ",
	"     License will therefore apply, along with any applicable section 7 additional terms, to the whole of the work, and all its ",
	"     parts, regardless of how they are packaged. This License gives no permission to license the work in any other way, but it ",
	"     does not invalidate such permission if you have separately received it.",
	"     d) If the work has interactive user interfaces, each must display Appropriate Legal Notices; however, if the Program has ",
	"     interactive interfaces that do not display Appropriate Legal Notices, your work need not make them do so.",
	"     ",
	"     A compilation of a covered work with other separate and independent works, which are not by their nature extensions of the ",
	"     covered work, and which are not combined with it such as to form a larger program, in or on a volume of a storage or distribution 	",
	"     medium, is called an “aggregate” if the compilation and its resulting copyright are not used to limit the access or legal ",
	"     rights of the compilation's users beyond what the individual works permit. Inclusion of a covered work in an aggregate does ",
	"     not cause this License to apply to the other parts of the aggregate.",
	"",
	"     6. Conveying Non-Source Forms.",
	"     You may convey a covered work in object code form under the terms of sections 4 and 5, provided that you also convey the ",
	"     machine-readable Corresponding Source under the terms of this License, in one of these ways:",
	"",
	"     a) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by ",
	"     the Corresponding Source fixed on a durable physical medium customarily used for software interchange.",
	"     b) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by ",
	"     a written offer, valid for at least three years and valid for as long as you offer spare parts or customer support for that ",
	"     product model, to give anyone who possesses the object code either (1) a copy of the Corresponding Source for all the software", 
	"     in the product that is covered by this License,	on a durable physical medium customarily used for software interchange, for a ",
	"     price no more than your reasonable cost of physically performing this conveying of source, or (2) access to copy the ",
	"     Corresponding Source from a network server at no charge.",
	"     c) Convey individual copies of the object code with a copy of the written offer to provide the Corresponding Source. This ",
	"     alternative is allowed only occasionally and noncommercially, and only if you received the object code with such an offer, ",
	"     in accord with subsection 6b.",
	"     d) Convey the object code by offering access from a designated place (gratis or for a charge), and offer equivalent access ",
	"     to the Corresponding Source in the same way through the same place at no further charge. You need not require recipients to ",
	"     copy the Corresponding Source along with the object code. If the place to copy the object code is a network server, the ",
	"     Corresponding Source may be on a different server (operated by you or a third party) that supports equivalent copying ",
	"     facilities, provided you maintain clear directions next to the object code saying where to find the Corresponding Source. ",
	"     Regardless of what server hosts the Corresponding Source, you remain obligated to ensure that it is available for as long as ",
	"     needed to satisfy these requirements.",
	"     e) Convey the object code using peer-to-peer transmission, provided you inform other peers where the object code and ",
	"     Corresponding Source of the work are being offered to the general public at no charge under subsection 6d.",
	"     ",
	"     A separable portion of the object code, whose source code is excluded from the Corresponding Source as a System Library, ",
	"     need not be included in conveying the object code work.",
	"",
	"     A “User Product” is either (1) a “consumer product”, which means any tangible personal property which is normally used for ",
	"     personal, family, or household purposes, or (2) anything designed or sold for incorporation into a dwelling. In determining ",
	"     whether a product is a consumer product, doubtful cases shall be resolved in favor of coverage. For a particular product ",
	"     received by a particular user, “normally used” refers to a typical or common use of that class of product, regardless of ",
	"     the status of the particular user or of the way in which the particular user actually uses, or expects or is expected to ",
	"     use, the product. A product is a consumer product regardless of whether the product has substantial commercial, industrial ",
	"     or non-consumer uses, unless such uses represent the only significant mode of use of the product.",
	"",
	"     “Installation Information” for a User Product means any methods, procedures, authorization keys, or other information ",
	"     required to install and execute modified versions of a covered work in that User Product from a modified version of its ",
	"     Corresponding Source. The information must suffice to ensure that the continued functioning of the modified object code ",
	"     is in no case prevented or interfered with solely because modification has been made.",
	"",
	"     If you convey an object code work under this section in, or with, or specifically for use in, a User Product, and the ",
	"     conveying occurs as part of a transaction in which the right of possession and use of the User Product is transferred ",
	"     to the recipient in perpetuity or for a fixed term (regardless of how the transaction is characterized), the Corresponding ",
	"     Source conveyed under this section must be accompanied by the Installation Information. But this requirement does not ",
	"     apply if neither you nor any third party retains the ability to install modified object code on the User Product (for ",
	"     example, the work has been installed in ROM).",
	"",
	"     The requirement to provide Installation Information does not include a requirement to continue to provide support service, ",
	"     warranty, or updates for a work that has been modified or installed by the recipient, or for the User Product in which it ",
	"     has been modified or installed. Access to a network may be denied when the modification itself materially and adversely ",
	"     affects the operation of the network or violates the rules and protocols for communication across the network.",
	"",
	"     Corresponding Source conveyed, and Installation Information provided, in accord with this section must be in a format that ",
	"     is publicly documented (and with an implementation available to the public in source code form), and must require no special ",
	"     password or key for unpacking, reading or copying.",
	"",
	"     7. Additional Terms.",
	"     “Additional permissions” are terms that supplement the terms of this License by making exceptions from one or more of its ",
	"     conditions. Additional permissions that are applicable to the entire Program shall be treated as though they were included ",
	"     in this License, to the extent that they are valid under applicable law. If additional permissions apply only to part of the ",
	"     Program, that part may be used separately under those permissions, but the entire Program remains governed by this License ",
	"     without regard to the additional permissions.",
	"",
	"     When you convey a copy of a covered work, you may at your option remove any additional permissions from that copy, or from ",
	"     any part of it. (Additional permissions may be written to require their own removal in certain cases when you modify the work.) ",
	"     You may place additional permissions on material, added by you to a covered work, for which you have or can give appropriate ",
	"     copyright permission.",
	"",
	"     Notwithstanding any other provision of this License, for material you add to a covered work, you may (if authorized by the ",
	"     copyright holders of that material) supplement the terms of this License with terms:",
	"",
	"     a) Disclaiming warranty or limiting liability differently from the terms of sections 15 and 16 of this License; or",
	"     b) Requiring preservation of specified reasonable legal notices or author attributions in that material or in the Appropriate ",
	"     Legal Notices displayed by works containing it; or",
	"     c) Prohibiting misrepresentation of the origin of that material, or requiring that modified versions of such material be ",
	"     marked in reasonable ways as different from the original version; or",
	"     d) Limiting the use for publicity purposes of names of licensors or authors of the material; or",
	"     e) Declining to grant rights under trademark law for use of some trade names, trademarks, or service marks; or",
	"     f) Requiring indemnification of licensors and authors of that material by anyone who conveys the material ",
	"     (or modified versions of it) with contractual assumptions of liability to the recipient, for any liability that these ",
	"     contractual assumptions directly impose on those licensors and authors.",
	"     ",
	"     All other non-permissive additional terms are considered “further restrictions” within the meaning of section 10. If the ",
	"     Program as you received it, or any part of it, contains a notice stating that it is governed by this License along with a ",
	"     term that is a further restriction, you may remove that term. If a license document contains a further restriction but permits ",
	"     relicensing or conveying under this License, you may add to a covered work material governed by the terms of that license ",
	"     document, provided that the further restriction does not survive such relicensing or conveying.",
	"",
	"     If you add terms to a covered work in accord with this section, you must place, in the relevant source files, a statement ",
	"     of the additional terms that apply to those files, or a notice indicating where to find the applicable terms.",
	"",
	"     Additional terms, permissive or non-permissive, may be stated in the form of a separately written license, or stated as ",
	"     exceptions; the above requirements apply either way.",
	"",
	"     8. Termination.",
	"     You may not propagate or modify a covered work except as expressly provided under this License. Any attempt otherwise to ",
	"     propagate or modify it is void, and will automatically terminate your rights under this License (including any patent licenses ",
	"     granted under the third paragraph of section 11).",
	"",
	"     However, if you cease all violation of this License, then your license from a particular copyright holder is reinstated (a) ",
	"     provisionally, unless and until the copyright holder explicitly and finally terminates your license, and (b) permanently, if ",
	"     the copyright holder fails to notify you of the violation by some reasonable means prior to 60 days after the cessation.",
	"",
	"     Moreover, your license from a particular copyright holder is reinstated permanently if the copyright holder notifies you of", 
	"     the violation by some reasonable means, this is the first time you have received notice of violation of this License (for any", 
	"     work) from that copyright holder, and you cure the violation prior to 30 days after your receipt of the notice.",
	"",
	"     Termination of your rights under this section does not terminate the licenses of parties who have received copies or rights ",
	"     from you under this License. If your rights have been terminated and not permanently reinstated, you do not qualify to receive", 
	"     new licenses for the same material under section 10.",
	"",
	"     9. Acceptance Not Required for Having Copies.",
	"     You are not required to accept this License in order to receive or run a copy of the Program. Ancillary propagation of a covered 	",
	"     work occurring solely as a consequence of using peer-to-peer transmission to receive a copy likewise does not require acceptance. 	",
	"     However, nothing other than this License grants you permission to propagate or modify any covered work. These actions infringe ",
	"     copyright if you do not accept this License. Therefore, by modifying or propagating a covered work, you indicate your acceptance 	",
	"     of this License to do so.",
	"",
	"     10. Automatic Licensing of Downstream Recipients.",
	"     Each time you convey a covered work, the recipient automatically receives a license from the original licensors, to run, modify ",
	"     and propagate that work, subject to this License. You are not responsible for enforcing compliance by third parties with this License.",
	"",
	"     An “entity transaction” is a transaction transferring control of an organization, or substantially all assets of one, or ",
	"     subdividing an organization, or merging organizations. If propagation of a covered work results from an entity transaction,", 
	"     each party to that transaction who receives a copy of the work also receives whatever licenses to the work the party's predecessor ",
	"     in interest had or could give under the previous paragraph, plus a right to possession of the Corresponding Source of the work ",
	"     from the predecessor in interest, if the predecessor has it or can get it with reasonable efforts.",
	"",
	"     You may not impose any further restrictions on the exercise of the rights granted or affirmed under this License. For example, ",
	"     you may not impose a license fee, royalty, or other charge for exercise of rights granted under this License, and you may not ",
	"     initiate litigation (including a cross-claim or counterclaim in a lawsuit) alleging that any patent claim is infringed by making,", 
	"     using, selling, offering for sale, or importing the Program or any portion of it.",
	"",
	"     11. Patents.",
	"     A “contributor” is a copyright holder who authorizes use under this License of the Program or a work on which the Program is based. ",
	"     The work thus licensed is called the contributor's “contributor version”.",
	"",
	"     A contributor's “essential patent claims” are all patent claims owned or controlled by the contributor, whether already acquired or ",
	"     hereafter acquired, that would be infringed by some manner, permitted by this License, of making, using, or selling its contributor ",
	"     version, but do not include claims that would be infringed only as a consequence of further modification of the contributor version. ",
	"     For purposes of this definition, “control” includes the right to grant patent sublicenses in a manner consistent with the ",
	"     requirements of this License.",
	"",
	"     Each contributor grants you a non-exclusive, worldwide, royalty-free patent license under the contributor's essential patent claims, ",
	"     to make, use, sell, offer for sale, import and otherwise run, modify and propagate the contents of its contributor version.",
	"",
	"     In the following three paragraphs, a “patent license” is any express agreement or commitment, however denominated, not to enforce a ",
	"     patent (such as an express permission to practice a patent or covenant not to sue for patent infringement). To “grant” such a patent ",
	"     license to a party means to make such an agreement or commitment not to enforce a patent against the party.",
	"",
	"     If you convey a covered work, knowingly relying on a patent license, and the Corresponding Source of the work is not available for ",
	"     anyone to copy, free of charge and under the terms of this License, through a publicly available network server or other readily ",
	"     accessible means, then you must either (1) cause the Corresponding Source to be so available, or (2) arrange to deprive yourself ",
	"     of the benefit of the patent license for this particular work, or (3) arrange, in a manner consistent with the requirements of this", 
	"     License, to extend the patent license to downstream recipients. “Knowingly relying” means you have actual knowledge that, but for the", 
	"     patent license, your conveying the covered work in a country, or your recipient's use of the covered work in a country, would infringe ",
	"     one or more identifiable patents in that country that you have reason to believe are valid.",
	"",
	"     If, pursuant to or in connection with a single transaction or arrangement, you convey, or propagate by procuring conveyance of, a ",
	"     covered work, and grant a patent license to some of the parties receiving the covered work authorizing them to use, propagate, modify ",
	"     or convey a specific copy of the covered work, then the patent license you grant is automatically extended to all recipients of the ",
	"     covered work and works based on it.",
	"",
	"     A patent license is “discriminatory” if it does not include within the scope of its coverage, prohibits the exercise of, or is conditioned 	",
	"     on the non-exercise of one or more of the rights that are specifically granted under this License. You may not convey a covered work if ",
	"     you are a party to an arrangement with a third party that is in the business of distributing software, under which you make payment to ",
	"     the third party based on the extent of your activity of conveying the work, and under which the third party grants, to any of the parties", 
	"     who would receive the covered work from you, a discriminatory patent license (a) in connection with copies of the covered work conveyed by 	",
	"     you (or copies made from those copies), or (b) primarily for and in connection with specific products or compilations that contain the ",
	"     covered work, unless you entered into that arrangement, or that patent license was granted, prior to 28 March 2007.",
	"",
	"     Nothing in this License shall be construed as excluding or limiting any implied license or other defenses to infringement that may ",
	"     otherwise be available to you under applicable patent law.",
	"",
	"     12. No Surrender of Others' Freedom.",
	"     If conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License, ",
	"     they do not excuse you from the conditions of this License. If you cannot convey a covered work so as to satisfy simultaneously your", 
	"     obligations under this License and any other pertinent obligations, then as a consequence you may not convey it at all. For example, ",
	"     if you agree to terms that obligate you to collect a royalty for further conveying from those to whom you convey the Program, the ",
	"     only way you could satisfy both those terms and this License would be to refrain entirely from conveying the Program.",
	"",
	"     13. Use with the GNU Affero General Public License.",
	"     Notwithstanding any other provision of this License, you have permission to link or combine any covered work with a work licensed ",
	"     under version 3 of the GNU Affero General Public License into a single combined work, and to convey the resulting work. The terms ",
	"     of this License will continue to apply to the part which is the covered work, but the special requirements of the GNU Affero ",
	"     General Public License, section 13, concerning interaction through a network will apply to the combination as such.",
	"",
	"     14. Revised Versions of this License.",
	"     The Free Software Foundation may publish revised and/or new versions of the GNU General Public License from time to time. Such ",
	"     new versions will be similar in spirit to the present version, but may differ in detail to address new problems or concerns.",
	"",
	"     Each version is given a distinguishing version number. If the Program specifies that a certain numbered version of the GNU General ",
	"     Public License “or any later version” applies to it, you have the option of following the terms and conditions either of that ",
	"     numbered version or of any later version published by the Free Software Foundation. If the Program does not specify a version ",
	"     number of the GNU General Public License, you may choose any version ever published by the Free Software Foundation.",
	"",
	"     If the Program specifies that a proxy can decide which future versions of the GNU General Public License can be used, that ",
	"     proxy's public statement of acceptance of a version permanently authorizes you to choose that version for the Program.",
	"",
	"     Later license versions may give you additional or different permissions. However, no additional obligations are imposed ",
	"     on any author or copyright holder as a result of your choosing to follow a later version.",
	"",
	"     15. Disclaimer of Warranty.",
	"     THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING ",
	"     THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR ",
	"     IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. ",
	"     THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ",
	"     ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.",
	"",
	"     16. Limitation of Liability.",
	"     IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO ",
	"     MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, ",
	"     INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO", 
	"     LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO ",
	"     OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.",
	"",
	"     17. Interpretation of Sections 15 and 16.",
	"     If the disclaimer of warranty and limitation of liability provided above cannot be given local legal effect according to ",
	"     their terms, reviewing courts shall apply local law that most closely approximates an absolute waiver of all civil liability ",
	"     in connection with the Program, unless a warranty or assumption of liability accompanies a copy of the Program in return for ",
	"     a fee."),quote=F)		
}

# Disclaimer
print(c("","","    AESOP: Analysis of Electrostatic Similarities of Proteins" ,
"",
"    Biomolecular Modeling & Design Lab",
"    Department of Bioengineering",
"    University of California, Riverside", 
"",
"    Copyright (C) 2013 Chris A. Kieslich, Ronald D. Gorham Jr., and Dimitrios Morikis",
"",
"    This program comes with ABSOLUTELY NO WARRANTY; for details type gpl.license().",
"    This is free software, and you are welcome to redistribute it",
"    under certain conditions; type gpl.license() for details.","",""),quote=F)


write.pqr3<-function (pdb = NULL, xyz = pdb$xyz, resno = NULL, resid = NULL, 
    eleno = NULL, elety = NULL, chain = NULL, o = NULL, b = NULL, 
    het = FALSE, append = FALSE, verbose = FALSE, chainter = FALSE, 
    file = "R.pdb") {
    if (is.null(xyz) || !is.numeric(xyz)) 
        stop("write.pqr: please provide a 'pdb' object or numeric 'xyz' vector")
    if (any(is.na(xyz))) 
        stop("write.pqr: 'xyz' coordinates must have no NA's.")
    if (is.vector(xyz)) {
        natom <- length(xyz)/3
        nfile <- 1
    }
    else if (is.matrix(xyz)) {
        stop("write.pqr: no multimodel PQR support")
    }
    else {
        stop("write.pdb: 'xyz' or 'pdb$xyz' must be either a vector or matrix")
    }
    card <- rep("ATOM", natom)
    if (!is.null(pdb)) {
        if (natom == 1) 
            pdb$atom <- t(as.matrix(pdb$atom))
        if (het) 
            card <- c(rep("ATOM", nrow(pdb$atom)), rep("HETATM", 
                nrow(pdb$het)))
        if (is.null(resno)) {
            resno = pdb$atom[, "resno"]
            if (het) {
                resno = c(resno, pdb$het[, "resno"])
            }
        }
        if (is.null(resid)) {
            resid = pdb$atom[, "resid"]
            if (het) {
                resid = c(resid, pdb$het[, "resid"])
            }
        }
        if (is.null(eleno)) {
            eleno = pdb$atom[, "eleno"]
            if (het) {
                eleno = c(eleno, pdb$het[, "eleno"])
            }
        }
        if (is.null(elety)) {
            elety = pdb$atom[, "elety"]
            if (het) {
                elety = c(elety, pdb$het[, "elety"])
            }
        }
        if (is.null(chain)) {
            chain = pdb$atom[, "chain"]
            if (het) {
                chain = c(chain, pdb$het[, "chain"])
            }
        }
        if (is.null(o)) {
            o = pdb$atom[, "o"]
            if (het) {
                o = c(o, pdb$het[, "o"])
            }
        }
        if (is.null(b)) {
            b = pdb$atom[, "b"]
            if (het) {
                b = c(b, pdb$het[, "b"])
            }
        }
        if (any(is.na(o))) {
            o = rep("1.00", natom)
        }
        if (any(is.na(b))) {
            b = rep("0.00", natom)
        }
        chain[is.na(chain)] = " "
    }
    else {
        if (is.null(resno)) 
            resno = c(1:natom)
        if (is.null(resid)) 
            resid = rep("ALA", natom)
        if (is.null(eleno)) 
            eleno = c(1:natom)
        if (is.null(elety)) 
            elety = rep("CA", natom)
        if (is.null(chain)) 
            chain = rep(" ", natom)
        if (is.null(o)) 
            o = rep("1.00", natom)
        if (is.null(b)) 
            b = rep("0.00", natom)
    }
    if (!is.logical(append)) 
        stop("write.pqr: 'append' must be logical TRUE/FALSE")
    if (length(as.vector(xyz))%%3 != 0) {
        stop("write.pqr: 'length(xyz)' must be divisable by 3.")
    }
    check.lengths <- sum(length(resno), length(resid), length(eleno), 
        length(elety), length(o), length(b))
    if (check.lengths%%natom != 0) {
        stop("write.pqr: the lengths of all input vectors != 'length(xyz)/3'.")
    }
    o <- as.numeric(o)
    b <- as.numeric(b)
    eleno <- as.character(eleno)
    resno <- as.character(resno)
    ter.lines <- (which(!duplicated(chain))[-1] - 1)
    if (nfile == 1) {
        coords <- matrix(round(as.numeric(xyz), 3), ncol = 3, 
            byrow = TRUE)
        if (verbose) {
            cat(paste("Writing 1 frame with", natom, "atoms "))
        }
        coords <- matrix(round(as.numeric(xyz), 3), ncol = 3, 
            byrow = TRUE)
        lines <- matrix(, ncol = 1, nrow = natom)
        cases <- matrix(1, ncol = 2, nrow = natom)
        cases[(nchar(eleno) > 5), 1] = 3
        cases[(nchar(elety) < 4), 2] = 0
        cases <- rowSums(cases)
        ind.1 <- which(cases == 1)
        ind.2 <- which(cases == 2)
        ind.3 <- which(cases == 3)
        ind.4 <- which(cases == 4)
        atom.print.1 <- function(card = "ATOM", eleno, elety, 
            alt = "", resid, chain = "", resno, insert = "", 
            x, y, z, o = "1.00", b = "0.00", segid = "") {
            format <- "%-6s%5s  %-3s%1s%-4s%1s%5s%1s%3s%9.3f%9.3f%9.3f%8.4f%7.4f%6s%4s"
            sprintf(format, card, eleno, elety, alt, resid, chain, 
                resno, insert, "", x, y, z, o, b, "", segid)
        }
        atom.print.2 <- function(card = "ATOM", eleno, elety, 
            alt = "", resid, chain = "", resno, insert = "", 
            x, y, z, o = "1.00", b = "0.00", segid = "") {
            format <- "%-6s%5s %-4s%1s%-4s%1s%5s%1s%3s%9.3f%9.3f%9.3f%8.4f%7.4f%6s%4s"
            sprintf(format, card, eleno, elety, alt, resid, chain, 
                resno, insert, "", x, y, z, o, b, "", segid)
        }
        atom.print.3 <- function(card = "ATOM", eleno, elety, 
            alt = "", resid, chain = "", resno, insert = "", 
            x, y, z, o = "1.00", b = "0.00", segid = "") {
            format <- "%-4s%7s  %-3s%1s%-4s%1s%5s%1s%3s%9.3f%9.3f%9.3f%8.4f%7.4f%6s%4s"
            sprintf(format, card, eleno, elety, alt, resid, chain, 
                resno, insert, "", x, y, z, o, b, "", segid)
        }
        atom.print.4 <- function(card = "ATOM", eleno, elety, 
            alt = "", resid, chain = "", resno, insert = "", 
            x, y, z, o = "1.00", b = "0.00", segid = "") {
            format <- "%-4s%7s %-4s%1s%-4s%1s%5s%1s%3s%9.3f%9.3f%9.3f%8.4f%7.4f%6s%4s"
            sprintf(format, card, eleno, elety, alt, resid, chain, 
                resno, insert, "", x, y, z, o, b, "", segid)
        }
        if (length(ind.1) > 0) {
            lines[ind.1, ] <- atom.print.1(card = card[ind.1], 
                eleno = eleno[ind.1], elety = elety[ind.1], resid = resid[ind.1], 
                chain = chain[ind.1], resno = resno[ind.1], x = coords[ind.1, 
                  1], y = coords[ind.1, 2], z = coords[ind.1, 
                  3], o = o[ind.1], b = b[ind.1])
        }
        if (length(ind.2) > 0) {
            lines[ind.2, ] <- atom.print.2(card = card[ind.2], 
                eleno = eleno[ind.2], elety = elety[ind.2], resid = resid[ind.2], 
                chain = chain[ind.2], resno = resno[ind.2], x = coords[ind.2, 
                  1], y = coords[ind.2, 2], z = coords[ind.2, 
                  3], o = o[ind.2], b = b[ind.2])
        }
        if (length(ind.3) > 0) {
            lines[ind.3, ] <- atom.print.3(card = card[ind.3], 
                eleno = eleno[ind.3], elety = elety[ind.3], resid = resid[ind.3], 
                chain = chain[ind.3], resno = resno[ind.3], x = coords[ind.3, 
                  1], y = coords[ind.3, 2], z = coords[ind.3, 
                  3], o = o[ind.3], b = b[ind.3])
        }
        if (length(ind.4) > 0) {
            lines[ind.4, ] <- atom.print.4(card = card[ind.4], 
                eleno = eleno[ind.4], elety = elety[ind.4], resid = resid[ind.4], 
                chain = chain[ind.4], resno = resno[ind.4], x = coords[ind.4, 
                  1], y = coords[ind.4, 2], z = coords[ind.4, 
                  3], o = o[ind.4], b = b[ind.4])
        }
        write.table(lines, file = file, quote = FALSE, row.names = FALSE, 
            col.names = FALSE, append = append)
    }
    else {
        if (verbose) {
            cat(paste("Writing", nfile, "frames with", natom, 
                "atoms"), "\n")
            cat("Frame Progress (x50) ")
        }
        stop("REMOVED code for multimodel PQR as these files dont have much support")
    }
    if (verbose) 
        cat(" DONE", "\n")
}

