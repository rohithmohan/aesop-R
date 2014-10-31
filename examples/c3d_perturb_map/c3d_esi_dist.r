rm(list=ls(all=TRUE))
library(bio3d)
source("../../source/AESOP.r")
source("../../source/AESOP_paths.r")
			   		
ion.strength <- 0.150 # Ionic strength of monovalent ions (NaCl-like) 
pdie <- 20.0 # Protein dielectric constant
sdie <- 78.54 # Solvent dielectric constant		

wkdir <- getwd() # Set working directory
cent.on <- "c3d" # Root name of protein to be used as center (as in PDB directory, without file extension)
pdb.file <- "c3d.pdb"
proj.name <-"c3d_ala" # Root name for created directories and files
pot.dir <- paste(wkdir,"/",proj.name,"_pot",sep="") # Directory generated for DX (electrostatic potential) files

pdb <- read.pdb(file=pdb.file) # Reads parent PDB
pqrs <- ala.scan(pdb.file, ff = "CHARMM") # Genreates PQR files by performing alanine scan

solv.in <- apbsin.solv(pdb,paste(cent.on,".pqr",sep=""),pqrs$dirs,s.only=TRUE, # Generates APBS input
		       pot=TRUE,s.keys=c(paste("ion 1",ion.strength,"2.0"), # Alters parameters for solvated state
					 paste("ion -1",ion.strength,"2.0"),
					 paste("pdie",pdie),
					 paste("sdie",sdie)),
				r.keys=c(paste("pdie",pdie),# Alters parameters for reference state
					 paste("sdie",pdie))) 

pots <- calc.pot(solv.in,pot.dir) # Generates electrostatic potentials based on APBS input (solv.in)
esi.distr(pot.dir,paste(strsplit(cent.on,split=".pqr")[[1]][1],".dx",sep=""), # Computes an ESD distance matrix for the potentials in pot.dir
		  paste(proj.name,"_esi_distr.dx",sep="")) 

save(pqrs,solv.in,pots,file=paste(proj.name,"_esi.dat",sep="")) # Saves generated data in a binary file for later use
