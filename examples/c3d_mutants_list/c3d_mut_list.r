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

AB <- read.pdb(file=pdb.file) # Reads in parent PDB | , het2atom=TRUE

## Determines amino acid differences between SUMO4 and SUMO2 and generates mutated structures for each individual mutation ## 
# seqs <- read.fasta("c3dmutants1.fasta") # Reads in SUMO2 sequence
# tmp.seq <- aa321(AB$atom[AB$calpha,"resid"]) # Extracts sequence for parent (SENP2-SUMO4)
# tar.seq <- c(seqs$ali, aa321(AB$atom[AB$calpha & AB$atom[,"chain"]=="B","resid"])) # Generates target sequence (SENP2-SUMO2) |test A -> B |test c(aa321(AB$atom[AB$calpha & AB$atom[,"chain"]=="A","resid"]),seqs$ali)

# mutres <- which(tmp.seq != tar.seq) # Determines which positions differ between SUMO4 and SUMO2 (relative to complex)
# mut2 <- tar.seq[tmp.seq != tar.seq] # Determines which amino acid types each differing position should be mutated to

# Formats list of mutations into correct format
#muts <- data.frame(res=mutres,aa=mut2) 
muts = read.csv("c3dmuts.csv")
mut.list <- lapply(1:dim(muts)[1],function(x) muts[x,])

pqrs <- mut.list.c(AB,chain.names,mut.list,ff="CHARMM", inputchain="a") # Generates PQRs for all mutants in mutlist | using CHARMM

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
a.solv.fe  <- calc.dGsolv(apbsin=solv.in,coul.ep=pdie,batch.dir=pqrs$dirs[2],pot.dir=pot.dirs[1])
# Calculates solvation energy for each mutant in chain B and saves electrostatic potentials
b.solv.fe  <- calc.dGsolv(apbsin=solv.in,coul.ep=pdie,batch.dir=pqrs$dirs[3],pot.dir=pot.dirs[2])

a.esd <- esd.dist(pot.dirs[1]) # Computes an ESD distance matrix for the potentials of chain A mutants

save(list=ls(),file=paste(chain.names[1],"_",chain.names[2],".ses",sep=""))

# Reformats mutant labels
a.comp.names <- str_replace(row.names(a.solv.fe),"c3d","c3d_cr2")
ab.a.fe <- ab.solv.fe[row.names(ab.solv.fe) %in% a.comp.names,]

# Calculates ddGsolv for mutants
a.ddGsolv <- ab.a.fe$dGsolv - a.solv.fe$dGsolv - b.solv.fe[chain.names[2],"dGsolv"]
names(a.ddGsolv) <- row.names(a.solv.fe) 

# Calculates dGsolu for mutants (bottom horizontal)
a.dGsolu <- ab.a.fe$Gsolu - a.solv.fe$Gsolu - b.solv.fe[chain.names[2],"Gsolu"]
names(a.dGsolu) <- row.names(a.solv.fe) 

# Calculates dGref for mutants (top horizontal)
a.dGref <- ab.a.fe$Gref - a.solv.fe$Gref - b.solv.fe[chain.names[2],"Gref"]
names(a.dGref) <- row.names(a.solv.fe) 

# Calculates dGcoul for mutants (Coulombic binding energies)
a.dGcoul <- ab.a.fe$Gcoul - a.solv.fe$Gcoul - b.solv.fe[chain.names[2],"Gcoul"]
names(a.dGcoul) <- row.names(a.solv.fe) 

# Calculates dGbind for mutants (ddGsolv + dGcoul)
a.dGbind <- a.ddGsolv + a.dGcoul
names(a.dGbind) <- row.names(a.solv.fe)


a.hc <- hclust(as.dist(a.esd), method = "average", members=NULL) # Performs hierarchical clustering based on the ESD distance matrix

# Reformats chain B labels for dendrogram
a.hc$labels <- unlist(strsplit(a.hc$labels,split=".dx"))
#names <- unlist(strsplit(unlist(strsplit(a.hc$labels[grep("_",a.hc$labels)],split=cent.on)),split="_"))
#a.hc$labels[grep("_",a.hc$labels)] <- names[seq(2,length(names),2)]
#The reformatting has been changed to account for the multiple mutations in the name (eg. A3A_A5A_A9A) - RM 03/04/14
chainsearch <- paste(chain.names[1],"_",sep="")
a.hc$labels[grep("_",a.hc$labels)] <-str_replace(a.hc$labels[grep("_",a.hc$labels)],chainsearch,"")
save(a.esd,a.hc,file=paste(chain.names[1],chain.names[2],"clust.dat",sep="_"))

names(a.dGbind) <- a.hc$labels

# Calculates relative dGbind for mutants (dGmut-dGWT)
a.dGbind.rel <- a.dGbind - a.dGbind[chain.names[1]]
names(a.dGbind.rel) <- names(a.dGbind)

# Plots and saves the clustering (using a dendrogram) with corresponding free energy plot for chain A mutants
svg(paste(chain.names[1],"_dendro.svg",sep=""))  
split.screen(c(2,1))
screen(1)
plot(as.dendrogram(a.hc), edgePar=list(col=1, lwd=2), horiz=F) # Plots the clustering using a dendrogram 
screen(2)
barplot(a.dGbind.rel[a.hc$order],xlab="Mutants",ylab="dGbind (kJ/mol)",las=2) # Plots solvation free energies of association
close.screen(all=T)
dev.off()
