rm(list=ls(all=TRUE))
library(bio3d)
source("../../source/AESOP.r")
source("../../source/AESOP_paths.r")

pH <- 7 # pH used PDB2PQR charge state assignments
ion.strength <- 0.050 # Ionic strength of monovalent ions (NaCl-like) 
pdie <- 20.0 # Protein dielectric constant
sdie <- 78.54 # Solvent dielectric constant		
	   			   
wkdir <- getwd() # Set working directory
pdb.dir <- paste(wkdir,"/LTP_ali_pdbs",sep="") # Directory containing superimposed PDB files
cent.on <- "1MZL" # Root name of protein to be used as center (as in PDB directory, without file extension)
proj.name <-"LTP" # Root name for created directories and files

pdb <- read.pdb(paste(pdb.dir,"/1MZL.pdb",sep=""))

pqr.dir <- paste(wkdir,"/",proj.name,"_pqr",sep="") # Directory generated for PQR files
pot.dir <- paste(wkdir,"/",proj.name,"_pot",sep="") # Directory generated for DX (electrostatic potential) files

pqrs <- batch.pqr(pdb.dir,pqr.dir,ff="CHARMM",flags=paste("--with-ph=",7,sep="")) # Generates PQR files the PDBs in pdb.dir using PDB2PQR
solv.in <- apbsin.solv(pdb,paste(cent.on,".pqr",sep=""),pqr.dir,s.only=TRUE, # Generates APBS input
		       pot=TRUE,s.keys=c(paste("ion 1",ion.strength,"2.0"), # Alters parameters for solvated state
					 paste("ion -1",ion.strength,"2.0"),
					 paste("pdie",pdie),
					 paste("sdie",sdie)),
				r.keys=c(paste("pdie",pdie),# Alters parameters for reference state
					 paste("sdie",pdie))) 

pots <- calc.pot(solv.in,pot.dir) # Generates electrostatic potentials based on APBS input (solv.in)

esd <- esd.dist(pot.dir) # Computes an ESD distance matrix for the potentials in pot.dir
hc <- hclust(as.dist(esd), method = "average", members=NULL) # Performs hierarchical clustering based on the ESD distance matrix
hc$labels <- unlist(strsplit(hc$labels,split=".dx")) # Reformats labels for dendrogram

#  Plots and saves the clustering using a dendrogram in SVG format
svg(paste(proj.name,"_dendro.svg",sep=""))  
#png(paste(proj_name,"_dendro.png",sep="")) # Alternative file types (PNG and PDF). Delete # at beginning to uncomment 
#pdf(paste(proj_name,"_dendro.pdf",sep="")) # Only use one of the three options at a time
par(mar=c(8, 2, 4, 4))
plot(as.dendrogram(hc), edgePar=list(col=3, lwd=4), horiz=F) # Plots clustering as dendrogram
dev.off()

save(pqrs,solv.in,pots,esd,hc,file=paste(proj.name,"_clust.dat")) # Saves generated data in a binary file for later use
