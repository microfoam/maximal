# maximal
maximal is a program to explore the world of maximal homology alignment (MHA), created and written by Albert Erives. 
Additional help options can be seen by running ./maximal without any options specified.

Version v3.43 is the most recent, robust, and easiest version.
The file "maximal" currently corresponds to a MAC OS-compatible executible.
The file "mha_v3.43.c" corresponds to the most recent C code.
Older code versions are saved in the "pointbreak" directory.

_______________________________________________
FILE: surf_wavereport.(log|mha)

This is the standard ouput report to which maximal writes. 
Every run of maximal appends run data to the end of this file and is a log of run-time information and performance.

Other informational report files also begin with the "surf_" base file name.
_______________________________________________
FILES: surfboard-cleanup_set_(glassy|gnarly), surfboard-get_worked(_clean), surfboard-check(mha|log), surfboard-Pop_Up, surfboard-pull_in

All files with the "surfboard" base name are bash script files.
The surfboard-get_worked files test all of the internal sample strings (-a to -h), the solved-tricksy, and Drosophila example corpus. 
The first script runs in default mode (glassy), while the second tests in both default and -x mode (gnarly), which adds extra transition squeezing. 
The surf_wavereport.log file in the top level directory corresponds to the output from the surfboard-cleanup_set_gnarly run.

The surfboard-Pop_Up file is a compile script, while the surfboard-Pull_In file is a sample script to begin to assemble sequences for MSA runs.
_______________________________________________
FILE: waves/solved-##-DESCRP_tricksy.txt ("##-DESCRP" provides a historical index number and a short 6 letter description of the originating issue)

Files named in this format represent short sequence text files with a tricksy microfoam nature. 
Typically, tricsky strings were at one point problem strings that were problematically cinched by maximal. 
After developing the code to conduct MHA successfully on a new tricksy string it is renamed as "solved", given a number, and a description. 
These file names and sequences should never be modified because they are used for historical comparisons of different versions of maximal. 

These files are now stored in waves/ subdirectory.

_______________________________________________
FILE: waves/Dxxxx_ex3.txt, Dxxxx_NEE.txt ("xxxx" is a Drosophila species name)

Files named in this format represent gene sequences from the Drosophila vnd locus. 
Sequences correspond either to a conserved protein-coding region (exon 3) or an intronic enhancer (the neurogenic ectoderm enhancer, or NEE) from several species. 
These file names and sequences should never be modified because they are used for historical comparisons of different developmenal versions of maximal. 

_______________________________________________
FILE: surf_foam_and_chowder.log, waves/surf_foam_and_chowder.log

If this file is present it has record of newly identified "tricksy" strings with problematic cinching as identified by the associated version of maximal. 
Typically, this file is committed after a -XX Fischer-Yates run and is named with additional version and run information.
_______________________________________________

This file is a work in progress.
