# lineup
lineup (pre-v5.01 versions were called 'maximal') is a prototype program designed to explore the world of micro-homology alignment (MHA), created and written by Albert Erives (albert-erives@uiowa.edu).
MHA methodology is philosophically-distinct from gapped alignment in embracing local microparalogy as a fundamental feature of biological sequence.

As a prototype, sketch implementation of MHA, the current versions of lineup feature a rich set of run options in order to encourage a free, open-minded, heuristics-friendly development cycle.
Program options can be seen by running ./lineup without any options specified.

Version v4.37 is the current version and the first to be tested without any known issues on both macOS and Linux platforms. 
To compile for either platform, run the script './surfboard-Pop_Up', which will produce an executable named 'lineup'.
For development and testing purposes, two executable binaries are included ('lineup-macOS-binary' and 'lineup-Linux-binary').
To use either of these pre-compiled versions with any of the 'surfboard_lowercase_namespace' script, the appropriate binary can be copied or renamed to 'lineup'.

Version numbers represent program development versions that increase in steps of 0.01.
We are 'pre-release' of the first official release, a little bit of the 'journey is the destination' with public access.
Nonetheless, at this moment and this perspective, we can say that the construction work has moved up from the south-east and south-west and joined here,
carried out by two great armies of workers, and some angels. 

Some of the file descriptions below may be a bit outdated but still convey the spirit of this project.

Note: Most of the scripts in this directory are for program development and testing purposes.

_______________________________________________
FILES: surfboard-*

Files with the "surfboard" base name are bash script files for compiling the program, running the program in different modes on batteries of tests, 
and/or for processing program output.
Many of these surfboard programs are for program development and testing.
These scripts facilitate having routine standard tests that use different run parameters invoked as command-line options.

The surfboard-Pop_Up file is a compile script.

The surfboard-cleanup_set-* files are for running benchmark tests. 
Some of these are set with options to run each test automatically, 
while others are set to run in a "visually-inspect" mode where the user is prompted to hit enter before each run.
The surfboard-powerwash runs all of the benchmark tests.

The surfboard-getworked_rgdoll is the metaphorical washing machine, the turbulent surf in which a downed surfer gets rolled like a ragdoll.
This script has been used to run the program on thousands of randomly shuffled biological sequences. 
This has been useful for identifying trick sequences used in program development.
It has also been useful as a control comparison group with biological sequences.

_______________________________________________
FILE: Surf_wavereport.(log|mha), Corduroy.log

This is the standard ouput report to which lineup writes. 
Every run of lineup appends run data to the end of this file and is a log of run-time information and performance.
Files that end in *.mha are typically temporary files.

Corduroy.log is a summary report based on information extracted from Surf_wavereport.log.
The Corduroy.log report gets new summaries appended to it after each batch of tests ('surfboard-cleanup_set*' scripts).

_______________________________________________
FILE: waves/solved-##-DESCRP_tricksy.txt ("##-DESCRP" provides a historical index number and a short 6 letter description of the originating issue)

Files named in this format represent short sequence text files with a tricksy microfoam nature. 
Typically, tricsky strings were at one point problem strings that were problematically cinched by lineup. 
After developing the code to conduct MHA successfully on a new tricksy string it is renamed as "solved", given a number, and a description. 
These file names and sequences should never be modified because they are used for historical comparisons of different versions of lineup. 
These files are now stored in waves/ subdirectory.

UPDATE #1 FOR RECENT VERSIONS (surf-themed): Solved synthetic sequences in the example corpus are saved in the waves/tubespit/ sub-directory. 
Solved natural sequences are saved in the waves/animals sub-directory. Unsolved tricksy sequences are saved in the waves/chowder/ sub-directory.

UPDATE #2: As hundreds more tricksy sequences have been solved, the sequences in the various sub-directories 
'waves/seals', 'waves/sharks', 'waves/tubespit', and 'waves/chowder'
now merely represent different generations of solved sequence patterns. 
To say a bit more, we have deemed it useful to keep these example corpora in separate sub-directories.

_______________________________________________
FILE: waves/Dxxxx_ex3.txt, Dxxxx_NEE.txt ("xxxx" is a Drosophila species name)

Files named in this format represent gene sequences from the Drosophila vnd locus. 
Sequences correspond either to a conserved protein-coding region (exon 3) or an intronic enhancer (the neurogenic ectoderm enhancer, or NEE) from several species. 
These file names and sequences should never be modified because they are used for historical comparisons of different developmenal versions of lineup. 

UPDATE #1: These files now live either in 'waves/seals' (protein-coding sequences) or 'waves/sharks' (non-protein-coding enhancer sequences). 

_______________________________________________
FILE: waves/foam_and_chowder.log

If this file is present it has record of existing "tricksy" strings with problematic cinching as identified by the associated version of lineup. 
Typically, this file is committed after a -XX Fischer-Yates run and is named with additional version and run information.
_______________________________________________

*Last updated*: 04/14/2024
