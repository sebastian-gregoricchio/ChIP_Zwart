###############
##
##  Script to rename long name .bam (eg 5905_25_wz3909_TGACTTCG_S35.bam) files from GCF core to wzNUMBER.bam files (wz3909.bam).
##  Deploy script from the directory that your bams are in: setwd("/DATA/your.name/your_bams/")
##
##############

# clean up any leftover R stuff laying around
rm(list = ls())

# make a new directory to put your renamed files in.
subDir <- "rename"
dir.create(file.path(subDir),showWarnings = FALSE)

# get list of .bam files in folder (long names)
dn = list.files(".", pattern="\\.bam*|\\.fastq*", full.names = FALSE)

# check if reads are paired-end
p = grep("_R2",dn)

# rename your files with only wz number
newfiles <- list()
tmp <- list()
for(i in 1:length(dn)){
  comp = dn[i]
  x <- unlist(strsplit(dn[i],"_"))#split name
  y <- grep("wz", x,ignore.case=TRUE) # get only 'wz' value out of split
  PE=list("R1","R2")
  r <- which(x %in% PE)
  if(length(p) == 0)
    {
    wzname <- x[[y]] # find for 'wz' name using value from split
    ext <- sub("^[^.]*.", "", dn[i])
    newfiles <- paste(subDir,"/",wzname,".",ext, sep="") # make new names with the folder 'rename/' first
    newfile <- paste(subDir,"/",wzname,"_2.",ext, sep="")
    if(!file.exists(newfiles))
      {
      file.rename(from=comp, to =newfiles) # save your renamed files into the rename folder
      tmp[i]<-print(paste(comp," saved to ", newfiles))
      }
    else if(!file.exists(newfile))
      {
      file.rename(from=comp, to =newfile)
      tmp[i]<-print(paste(newfiles," already exists.  ",comp," saved to ", newfile))
      }
    else
      {
      tmp[i]<-print(paste(newfiles," and ",newfile," already exist. File not renamed. "))
      }
    }
    else if (length(p) != 0)
    {
    if (length(r) > 2)
      {
      tmp[i]<-print("Error. There is more than one read tag (i.e. R1 or R2) in your file name.")
      }
    else if(length(r) == 1)
      {
        wzname <- x[[y]] # find for 'wz' name using value from split
        ext <- sub("^[^.]*.", "", dn[i])
        read <- x[[r]]
        newfiles <- paste(subDir,"/",wzname,"_",read,".",ext, sep="") # make new names with the folder 'rename/' first
        newfile <- paste(subDir,"/",wzname,"_",read,"_2.",ext, sep="")
        if(!file.exists(newfiles))
          {
          file.rename(from=comp, to =newfiles) # save your renamed files into the rename folder
          tmp[i]<-print(paste(comp," saved to ", newfiles))
          }
        else if(!file.exists(newfile))
          {
          file.rename(from=comp, to =newfile)
          tmp[i]<-print(paste(newfiles," already exists.  ",comp," saved to ", newfile))
          }
        else
          {
          tmp[i]<-print(paste(newfiles," and ",newfile," already exist. File not renamed. "))
          }
      }
    }
}
tt <-as.data.frame(do.call(rbind,tmp))
write.table(tt, file="rename_log.txt",row.names = F,col.names = F)
