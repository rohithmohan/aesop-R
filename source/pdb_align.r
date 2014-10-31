rm(list=ls(all=TRUE))
library(bio3d)

##-- Read a folder/directory of PDB files
pdb.path <- "./cr2_20ns_md_traj_pdb_30sep10" # Folder of PDBs to be aligned
cent_on <- "cr2_0.pdb" # Name of PDB to which all others should be superimposed
pdb.out <- "./cr2_ali_pdbs" # Output directory for aligned PDBs

if (!file.exists(path = pdb.out)){
		dir.create(pdb.out)
}

files  <- list.files(path=pdb.path ,pattern="[.pdb]")

##-- Extract sequences
pdb.xyz <- NULL
for(file in files){
	pdb <- read.pdb(paste(pdb.path,"/",file,sep=""))
	pdb.xyz <- rbind(pdb.xyz,pdb$xyz)
}

sel <- atom.select(pdb,string="//////CA/")

rmsd.before <- mean(rmsd(pdb.xyz[,sel$xyz]))

xyz <- fit.xyz(fixed  = pdb.xyz[grep(cent_on,files),],
               mobile = pdb.xyz,
               fixed.inds  = sel$xyz,
               mobile.inds = sel$xyz)

rmsd.after <- mean(rmsd(xyz[,sel$xyz]))

for(i in 1:length(files)){
	write.pdb(pdb=pdb,xyz=xyz[i,],file=paste(pdb.out,"/",files[i],sep=""))
}
