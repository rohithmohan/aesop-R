rm(list=ls(all=TRUE))
library(bio3d)

pdb.path <- "./LTP_pdb" # Path to PDBs to be superimposed
cent_on <- "1MZL.pdb" # Name of PDB to which all others will be fitted
pdb.out <- "./LTP_ali_pdbs" # Directory where aligned PDBs will be saved
if (!file.exists(path = pdb.out)){
		dir.create(pdb.out)
	}

files  <- list.files(path=pdb.path ,pattern="[.pdb]",full.names=TRUE)

##-- Extract sequences
raw <- NULL
for(i in 1:length(files)) {
	pdb <- read.pdb(files[i])
	raw <- seqbind(raw, aa321(pdb$atom[pdb$calpha,"resid"]))
}

##-- Align sequences
aln <- seqaln(raw, id=basename(files),file="seqaln.fa")

##-- Read Aligned PDBs
pdbs <- read.fasta.pdb(aln, pdb.path=pdb.path, pdbext = "")

gaps <- gap.inspect(pdbs$xyz)

rmsd.before <- mean(rmsd(pdbs$xyz[,gaps$f.inds]))

xyz <- fit.xyz( fixed  = pdbs$xyz[cent_on,],
               mobile = pdbs,
               fixed.inds  = gaps$f.inds,
               mobile.inds = gaps$f.inds,
               pdb.path=pdb.path,
	       	   outpath = paste(pdb.out,"/",sep=""),
               full.pdbs = TRUE )

rmsd.after <- mean(rmsd(xyz[,gaps$f.inds]))

files  <- list.files(path=pdb.out,pattern="[.pdb]",full.names=TRUE)

for(file in files){
	file.rename(file,strsplit(file,split="_flsq.pdb")[[1]][1])
}