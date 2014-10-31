rm(list=ls(all=TRUE))

#Name of script to be executed
chimera <- '/Applications/Chimera.app/Contents/MacOS/chimera'
#Directory where this script is located
work_dir <- getwd()
#Directory where input files are located
input_dir <- 'LTP_pot'
#Directory where output files should be written (can be new folder)
output_dir <- 'LTP_images'

#Make a folder to save all the images
if (!file.exists(path = output_dir)){
		dir.create(output_dir)
	}
	
files  <- list.files(path=paste(work_dir,"/",input_dir,sep=""),pattern="[.dx]")

print(files)

for(file in files){
	output <- file(paste(work_dir,'/chimera_images.cmd',sep=""),'w')
	root <- unlist(strsplit(file,split=".dx"))[1]
	cat(c(paste('open',paste(work_dir,"/",input_dir,"/",file,sep="")),
							 "volume #0 show style surface surfaceSmoothing true level -1 color red level 1 color blue",
			                 "set bg_color white",
						     "scale .9", # Can be changed depending on protein/potential size
						     
						     # rotated 0 degrees about vertical 
			                 paste("copy file",paste(work_dir,"/",output_dir,"/",root,"_1.png",sep=""),"png supersample 2"), 
			                 
			                 # rotated 90 degrees about vertical
			                 "turn y 90",
			                 paste("copy file",paste(work_dir,"/",output_dir,"/",root,"_2.png",sep=""),"png supersample 2"), 
			                 
			                 # rotated 180 degrees about vertical
			                 "turn y 90",
			                 paste("copy file",paste(work_dir,"/",output_dir,"/",root,"_3.png",sep=""),"png supersample 2"), 
			                 
			                 # rotated 270 degrees about vertical
			                 "turn y 90",
			                 paste("copy file",paste(work_dir,"/",output_dir,"/",root,"_4.png",sep=""),"png supersample 2"), 
			                 
			                 "stop"),file=output,sep="\n")

	system(paste(chimera,'chimera_images.cmd'))
	close(output)
}		
