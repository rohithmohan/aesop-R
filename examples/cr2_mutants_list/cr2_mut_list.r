rm(list=ls(all=TRUE))
library(bio3d)
library(stringr)
source("../../source/AESOP.r")
source("../../source/AESOP_paths.r")
			   		
ion.strength <- 0.150 # Ionic strength of monovalent ions (NaCl-like) 
pdie <- 20.0 # Protein dielectric constant
sdie <- 78.54 # Solvent dielectric constant		
#dime <- c(65,65,65) # Number of grid points used by APBS (Resolution set low for example)
	   
wkdir <- getwd() # Set working directory
pdb.file <- "c3d_cr2.pdb" # Parent PDB file name
cent.on <- "c3d_cr2" # Root name of protein to be used as center (parent)
chain.names <- c("c3d","cr2") # Names of protein chains
pot.dirs <- paste(wkdir,"/",paste(chain.names,"pot",sep="_"),sep="") # Directories for potential (DX) files

pdb <- read.pdb(pdb.file)
pdb$atom[,"b"] <- "0.00"
write.pdb(pdb,file=pdb.file)

AB <- read.pdb(file=pdb.file) # Reads in parent PDB

## Determines amino acid differences between SUMO4 and SUMO2 and generates mutated structures for each individual mutation ## 
# seqs <- read.fasta("cr2mutants1.fasta") # Reads in SUMO2 sequence
# tmp.seq <- aa321(AB$atom[AB$calpha,"resid"]) # Extracts sequence for parent (SENP2-SUMO4)
# tar.seq <- c(aa321(AB$atom[AB$calpha & AB$atom[,"chain"]=="A","resid"]),seqs$ali) # Generates target sequence (SENP2-SUMO2) |test A -> B |test c(aa321(AB$atom[AB$calpha & AB$atom[,"chain"]=="A","resid"]),seqs$ali)

# mutres <- which(tmp.seq != tar.seq) # Determines which positions differ between SUMO4 and SUMO2 (relative to complex)
# mut2 <- tar.seq[tmp.seq != tar.seq] # Determines which amino acid types each differing position should be mutated to

# Formats list of mutations into correct format
#muts <- data.frame(res=mutres,aa=mut2) 
muts = read.csv("cr2muts.csv")
mut.list <- lapply(1:dim(muts)[1],function(x) muts[x,])

pqrs <- mut.list.c(AB,chain.names,mut.list,ff="CHARMM", inputchain="b") # Generates PQRs for all mutants in mutlist | using CHARMM

solv.in <- apbsin.solv(AB,paste(cent.on,".pqr",sep=""),pqrs$dirs[1],pot=TRUE, # Generates APBS input
		       s.keys=c(paste("ion 1",ion.strength,"2.0"), # Alters parameters for solvated state
			 paste("ion -1",ion.strength,"2.0"),
			 paste("pdie",pdie),
			 paste("sdie",sdie)),
		       r.keys=c(paste("pdie",pdie),# Alters parameters for reference state
			 paste("sdie",pdie)),
		       c.path=pqrs$dirs[1]) 

# Calculates solvation energy (vertical process) for each mutant in complex (no potentials)
ab.solv.fe <- calc.dGsolv(apbsin=solv.in,coul.ep=pdie,batch.dir=pqrs$dirs[1]) 
# Calculates solvation energy for each mutant in chain A and saves electrostatic potentials
a.solv.fe  <- calc.dGsolv(apbsin=solv.in,coul.ep=pdie,batch.dir=pqrs$dirs[2])
# Calculates solvation energy for each mutant in chain B and saves electrostatic potentials
b.solv.fe  <- calc.dGsolv(apbsin=solv.in,coul.ep=pdie,batch.dir=pqrs$dirs[3],pot.dir=pot.dirs[2])

b.esd <- esd.dist(pot.dirs[2]) # Computes an ESD distance matrix for the potentials of chain B mutants

save(list=ls(),file=paste(chain.names[1],"_",chain.names[2],".ses",sep=""))

# Reformats mutant labels
b.comp.names <- paste(chain.names[1],row.names(b.solv.fe),sep="_")
ab.b.fe <- ab.solv.fe[row.names(ab.solv.fe) %in% b.comp.names,]

# Calculates ddGsolv for mutants
b.ddGsolv <- ab.b.fe$dGsolv - a.solv.fe[chain.names[1],"dGsolv"] - b.solv.fe$dGsolv
names(b.ddGsolv) <- row.names(b.solv.fe) 

# Calculates dGsolu for mutants (bottom horizontal)
b.dGsolu <- ab.b.fe$Gsolu - a.solv.fe[chain.names[1],"Gsolu"] - b.solv.fe$Gsolu
names(b.dGsolu) <- row.names(b.solv.fe) 

# Calculates dGref for mutants (top horizontal)
b.dGref <- ab.b.fe$Gref - a.solv.fe[chain.names[1],"Gref"] - b.solv.fe$Gref
names(b.dGref) <- row.names(b.solv.fe)

# Calculates dGcoul for mutants (Coulombic binding energies)
b.dGcoul <- ab.b.fe$Gcoul - a.solv.fe[chain.names[1],"Gcoul"] - b.solv.fe$Gcoul
names(b.dGcoul) <- row.names(b.solv.fe) 

# Calculates dGbind for mutants (ddGsolv + dGcoul)
b.dGbind <- b.ddGsolv + b.dGcoul
names(b.dGbind) <- row.names(b.solv.fe)
b.hc <- hclust(as.dist(b.esd), method = "average", members=NULL) # Performs hierarchical clustering based on the ESD distance matrix
# Reformats chain B labels for dendrogram
b.hc$labels <- unlist(strsplit(b.hc$labels,split=".dx"))
#names <- unlist(strsplit(unlist(strsplit(b.hc$labels[grep("_",b.hc$labels)],split=cent.on)),split="_"))
#b.hc$labels[grep("_",b.hc$labels)] <- names[seq(2,length(names),2)]
#The reformatting has been changed to account for the multiple mutations in the name (eg. A3A_A5A_A9A) - RM 03/04/14
chainsearch <- paste(chain.names[2],"_",sep="")
b.hc$labels[grep("_",b.hc$labels)] <-str_replace(b.hc$labels[grep("_",b.hc$labels)],chainsearch,"")
	
names(b.dGbind) <- b.hc$labels

# Calculates relative dGbind for mutants (dGmut-dGWT)
b.dGbind.rel <- b.dGbind - b.dGbind[chain.names[2]]
names(b.dGbind.rel) <- names(b.dGbind)

# Plots and saves the clustering (using a dendrogram) with corresponding free energy plot for chain B mutants
svg(paste(chain.names[2],"_dendro.svg",sep=""))  
split.screen(c(2,1))
screen(1)
plot(as.dendrogram(b.hc), edgePar=list(col=1, lwd=2), horiz=F) # Plots the clustering using a dendrogram 
screen(2)
barplot(b.dGbind.rel[b.hc$order],xlab="Mutants",ylab="dGbind (kJ/mol)",las=2) # Plots solvation free energies of association
close.screen(all=T)
dev.off()

