# maximal
maximal is a prototype program designed to explore the world of micro-homology alignment (MHA), created and written by Albert Erives (albert-erives@uiowa.edu).
MHA methodology is philosophically-distinct from gapped alignment in embracing local microparalogy as a fundamental feature of biological sequence.

As a prototype implementation of MHA, the current versions of maximal feature a rich set of run options in order to encourage a free, open-minded, heuristics-friendly development cycle.
Program options can be seen by running ./maximal without any options specified.

Version v4.22 is the most recent stable version.
Version v3.64 is the version associated with the most recent preprint (see 'Pointbreak/' sub-directory).

The file "maximal" currently corresponds to a MAC OS-compatible executible.
Older code versions are saved in the "pointbreak" directory.

_______________________________________________
FILE: Surf_wavereport.(log|mha)

This is the standard ouput report to which maximal writes. 
Every run of maximal appends run data to the end of this file and is a log of run-time information and performance.

Other informational report files also begin with the "Surf_" base file name.
_______________________________________________
FILES: surfboard-cleanup_set_(all|glassy|gnarly +/- RC), surfboard-get_worked(_clean), surfboard-check(mha|log), surfboard-Pop_Up, surfboard-pull_in

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

UPDATE FOR RECENT VERSIONS (surf-themed): Solved synthetic sequences in the example corpus are saved in the waves/tubespit/ sub-directory. Solved
natural sequences are saved in the waves/animals sub-directory. Unsolved tricksy sequences are saved in the waves/chowder/ sub-directory.

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

