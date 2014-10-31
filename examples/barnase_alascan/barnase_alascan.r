rm(list=ls(all=TRUE))
library(bio3d)
source("../../source/AESOP.r")
source("../../source/AESOP_paths.r")
			   		
ion.strength <- 0.150 # Ionic strength of monovalent ions (NaCl-like) 
pdie <- 20.0 # Protein dielectric constant
sdie <- 78.54 # Solvent dielectric constant		
	   
wkdir <- getwd() # Set working directory
pdb.name <- "barnase.pdb" # Parent PDB file name
cent.on <- "barnase" # Root name of protein to be used as center (parent)
pot.dir <- paste(wkdir,"/",paste(cent.on,"pot",sep="_"),sep="") # Directory generated for DX (electrostatic potential) files

### Example Bio3D code for adding chain names and correcting residue naming ###
AB <- read.pdb(file="bn_bs.pdb") # Read in PDB of barnase-barstar complex
#! Note: Chain names must be used in PDB and can be added as follows !#
AB$atom[atom.select(AB,resno=1:110,verbose=F)$atom,"chain"] <- "A" # Sets residues 1 to 110 to chain A
AB$atom[atom.select(AB,resno=111:199,verbose=F)$atom,"chain"] <- "B" # Sets residues 111 to 199 to chain B
AB$atom[atom.select(AB,resid="HSP",verbose=F)$atom,"resid"] <- "HIS" # Renames histidines from HSP (CHARMM) to HIS
A <- trim.pdb(AB,atom.select(AB,chain="A",verbose=F)) # Extracts chain A (barnase) from complex
write.pdb(A,file=pdb.name) # Writes barnase PDB to file


pqrs <- ala.scan(pdb.name,ff="CHARMM") # Genreates PQR files by performing alanine scan
solv.in <- apbsin.solv(A,paste(cent.on,".pqr",sep=""),pqrs$dirs[1],pot=TRUE, # Generates APBS input
		       s.keys=c(paste("ion 1",ion.strength,"2.0"), # Alters parameters for solvated state
			 			paste("ion -1",ion.strength,"2.0"),
			 			paste("pdie",pdie),
			 			paste("sdie",sdie)),
		       r.keys=c(paste("pdie",pdie),# Alters parameters for reference state
			 			paste("sdie",pdie)),
		       c.path=pqrs$dirs[1]) 

a.pots  <- calc.pot(solv.in,pot.dir) # Generates electrostatic potentials based on APBS input (solv.in)
a.esd <- esd.dist(pot.dir) # Computes an ESD distance matrix for the potentials in pot.dir
a.hc <- hclust(as.dist(a.esd), method = "average", members=NULL) # Performs hierarchical clustering based on the ESD distance matrix

# Reformats labels for dendrogram
a.hc$labels <- unlist(strsplit(a.hc$labels,split=".dx"))
names <- unlist(strsplit(unlist(strsplit(a.hc$labels[grep("_",a.hc$labels)],split=cent.on)),split="_"))
a.hc$labels[grep("_",a.hc$labels)] <- paste(names[seq(2,length(names),3)],names[seq(3,length(names),3)])


svg(paste(cent.on,"_dendro.svg",sep="")) #  Plots and saves the clustering using a dendrogram 
plot(as.dendrogram(a.hc), edgePar=list(col=3, lwd=4), horiz=F) # Plots clustering as dendrogram
dev.off()

save(pqrs,solv.in,a.pots,a.esd,a.hc,file=paste(cent.on,"_clust.dat")) # Saves generated data in a binary file for later use

