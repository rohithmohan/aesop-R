rm(list=ls(all=TRUE))
library(bio3d)
source("../../source/AESOP.r")
source("../../source/AESOP_paths.r")
			   		
ion.strength <- 0.150 # Ionic strength of monovalent ions (NaCl-like) 
pdie <- 20.0 # Protein dielectric constant
sdie <- 78.54 # Solvent dielectric constant		
			   			   
wkdir <- getwd() # Set working directory

pdb.name <- "barnase_barstar.pdb" # Parent PDB file name
cent.on <- "barnase_barstar" # Root name of protein to be used as center (parent) 
chain.names <- c("barnase","barstar") # Names of protein chains
pot.dirs <- paste(wkdir,"/",paste(chain.names,"pot",sep="_"),sep="") # Directories for potential (DX) files

pdb <- read.pdb(pdb.name)
pdb$atom[,"b"] <- "0.00"
write.pdb(pdb,file=pdb.name)

AB <- read.pdb(file=pdb.name) # Read in parent PDB
pqrs <- ala.scan.c(pdb.name,chain.names,ff="CHARMM") # Genreates PQR files by performing alanine scan
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
b.esd <- esd.dist(pot.dirs[2]) # Computes an ESD distance matrix for the potentials of chain B mutants

save(list=ls(),file=paste(chain.names[1],"_",chain.names[2],".ses",sep=""))

# Reformats mutant labels
tmp <- sapply(1:length(row.names(a.solv.fe)), function(i) unlist(strsplit(row.names(a.solv.fe)[i],split="_")))
a.comp.names <- sapply(1:length(row.names(a.solv.fe)), function(i) paste(tmp[[i]][1],chain.names[2],tmp[[i]][2],tmp[[i]][3],sep="_"))
a.comp.names <- unlist(strsplit(a.comp.names,split="_NA_NA"))
b.comp.names <- paste(chain.names[1],row.names(b.solv.fe),sep="_")
ab.a.fe <- ab.solv.fe[row.names(ab.solv.fe) %in% a.comp.names,]
ab.b.fe <- ab.solv.fe[row.names(ab.solv.fe) %in% b.comp.names,]

# Calculates ddGsolv for mutants
a.ddGsolv <- ab.a.fe$dGsolv - a.solv.fe$dGsolv - b.solv.fe[chain.names[2],"dGsolv"]
names(a.ddGsolv) <- row.names(a.solv.fe) 
b.ddGsolv <- ab.b.fe$dGsolv - a.solv.fe[chain.names[1],"dGsolv"] - b.solv.fe$dGsolv
names(b.ddGsolv) <- row.names(b.solv.fe) 

# Calculates dGsolu for mutants (bottom horizontal)
a.dGsolu <- ab.a.fe$Gsolu - a.solv.fe$Gsolu - b.solv.fe[chain.names[2],"Gsolu"]
names(a.dGsolu) <- row.names(a.solv.fe) 
b.dGsolu <- ab.b.fe$Gsolu - a.solv.fe[chain.names[1],"Gsolu"] - b.solv.fe$Gsolu
names(b.dGsolu) <- row.names(b.solv.fe) 

# Calculates dGref for mutants (top horizontal)
a.dGref <- ab.a.fe$Gref - a.solv.fe$Gref - b.solv.fe[chain.names[2],"Gref"]
names(a.dGref) <- row.names(a.solv.fe) 
b.dGref <- ab.b.fe$Gref - a.solv.fe[chain.names[1],"Gref"] - b.solv.fe$Gref
names(b.dGref) <- row.names(b.solv.fe)

# Calculates dGcoul for mutants (Coulombic binding energies)
a.dGcoul <- ab.a.fe$Gcoul - a.solv.fe$Gcoul - b.solv.fe[chain.names[2],"Gcoul"]
names(a.dGcoul) <- row.names(a.solv.fe) 
b.dGcoul <- ab.b.fe$Gcoul - a.solv.fe[chain.names[1],"Gcoul"] - b.solv.fe$Gcoul
names(b.dGcoul) <- row.names(b.solv.fe) 

# Calculates dGbind for mutants (ddGsolv + dGcoul)
a.dGbind <- a.ddGsolv + a.dGcoul
names(a.dGbind) <- row.names(a.solv.fe)
b.dGbind <- b.ddGsolv + b.dGcoul
names(b.dGbind) <- row.names(b.solv.fe)

a.hc <- hclust(as.dist(a.esd), method = "average", members=NULL) # Performs hierarchical clustering based on the ESD distance matrix
b.hc <- hclust(as.dist(b.esd), method = "average", members=NULL) # Performs hierarchical clustering based on the ESD distance matrix

# Reformats chain A labels for dendrogram
for (i in 1:length(a.hc$labels)) {
	names <- unlist(strsplit(unlist(strsplit(a.hc$labels[i],split=".dx"))[1],"_"))
	if (length(names) > 1) {
		a.hc$labels[i] <- paste(unlist(strsplit(names[3],""))[1],unlist(strsplit(names[2],"[A-Z]"))[1],unlist(strsplit(names[3],""))[3],sep="")
	} else {
		a.hc$labels[i] <- "WT"
	}
}
	
# Reformats chain B labels for dendrogram
for (i in 1:length(b.hc$labels)) {
	names <- unlist(strsplit(unlist(strsplit(b.hc$labels[i],split=".dx"))[1],"_"))
	if (length(names) > 1) {
		b.hc$labels[i] <- paste(unlist(strsplit(names[3],""))[1],unlist(strsplit(names[2],"[A-Z]"))[1],unlist(strsplit(names[3],""))[3],sep="")
	} else {
		b.hc$labels[i] <- "WT"
	}
}
	
names(a.dGbind) <- a.hc$labels
names(b.dGbind) <- b.hc$labels

# Calculates relative dGbind for mutants (dGmut-dGWT)
a.dGbind.rel <- a.dGbind - a.dGbind["WT"]
names(a.dGbind.rel) <- names(a.dGbind)
b.dGbind.rel <- b.dGbind - b.dGbind["WT"]
names(b.dGbind.rel) <- names(b.dGbind)

# Plots and saves the clustering (using a dendrogram) with corresponding free energy plot for chain A mutants
svg(paste(chain.names[1],"_dendro.svg",sep=""))  
split.screen(c(2,1))
screen(1)
plot(as.dendrogram(a.hc), edgePar=list(col=1, lwd=2), horiz=F) # Plots the clustering using a dendrogram 
screen(2)
barplot(a.dGbind.rel[a.hc$order],xlab="Mutants",ylab="dGbind (kJ/mol)",las=2) # Plots solvation free energies of association
close.screen(all=T)
dev.off()

# Plots and saves the clustering (using a dendrogram) with corresponding free energy plot for chain B mutants
svg(paste(chain.names[2],"_dendro.svg",sep=""))  
split.screen(c(2,1))
screen(1)
plot(as.dendrogram(b.hc), edgePar=list(col=1, lwd=2), horiz=F) # Plots the clustering using a dendrogram 
screen(2)
barplot(b.dGbind.rel[b.hc$order],xlab="Mutants",ylab="dGbind (kJ/mol)",las=2) # Plots solvation free energies of association
close.screen(all=T)
dev.off()
