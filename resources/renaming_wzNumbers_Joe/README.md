----------------------------------
This fastq/bam Zwart-specific renaming procedure has been originally developed by Joseph (Joe) Seifert and can be found at his [*GitHub page*](https://github.com/csijcs/snakepipes#renaming-files) at the [*Renaming files* section](https://github.com/csijcs/snakepipes#renaming-files).

<hr style="border:2px">

<br>

## Renaming files
**NOTE** - all of you're sequencing filenames should contain a wz number (i.e. wz3909). Make sure to submit your samples with a wz number in the name or this script will not work. If there are two samples with the same wz number (i.e. same sample split across two lanes) the second file will be renamed wzNUMBER_2 (i.e. wz3909_2). If there are more than two (not likely, but possible), it will give an error and not rename your additional files. If you do actually have more than two (i.e. same sample split across more than two lanes), seek professional help.

Before starting a pipeline, it's best to rename your files. The files from the core come with a very long filename (i.e. 5905_25_wz3909_TGACTTCG_S35.bam) and we will shorten this to just the wz number (i.e. wz3909.bam).

To accomplish this, we have provided an R script above (rename_files.R). This script can either be run from within R, or from the terminal. To run from within R, set your working directory to the folder contaning your files (bam or fastq):

``setwd("/DATA/first_initial.last_name/your_files/")``


If you copy the script (located at `resources/renaming_wzNumbers_Joe/rename_files.R`) to the same folder as your files, you can run:

``source("rename_files.R")``


Otherwise this can be run from the RStudio script window.

If you prefer, you can also run from terminal by copying the script into the folder containing your files and running:

``Rscript /path/to/rename_files.R``

Either way, this will rename all your files and move them into a folder called `rename`. All of the files should have been moved into this folder, so if there are any remaining then something went wrong and you should seek help.

Once your files are renamed, you are now ready to proceed with the appropriate pipeline below.
