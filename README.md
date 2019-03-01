# maximal
maximal is a program to explore the world of maximal homology alignment (MHA), created and written by Albert Erives. Additional help options can be seen by running ./maximal without any options specified.

Version v3.39 is the most recent, robust, and easiest version.
_______________________________________________
FILE: output.mha

This is the standard ouput report to which maximal writes. 
Every run of maximal appends run data to the end of this file.

_______________________________________________
FILE: corpus/solved-##-descrp_tricksy.txt ("##-descrp" provides a historical index number and a short 6 letter description of the originating issue)

Files named in this format represent short sequence text files with a tricksy microfoam nature. Typically, tricsky strings were at one point problem strings that were problematically cinched by maximal. After developing the code to conduct MHA successfully on a new tricksy string it is renamed as "solved", given a number, and a description. These file names and sequences should rarely be modified because they are used for historical comparisons of different versions of maximal. For example, even changing the length of the sequence in order to focus on the particular microfoam knot of interest will likely result in a change in the string's width cinch ratio (WCR) post-cinching.

These files are now stored in corpus/ subdirectory.

_______________________________________________
FILE: corpus/Dxxxx_ex3.txt, Dxxxx_NEE.txt ("xxxx" is a Drosophila species name)

Files named in this format represent vnd exon 3 sequences or NEE sequences from one of several species of flies. These file names and sequences should rarely be modified because they are used for historical comparisons of different versions of maximal. For example, even changing the length of the sequence in order to focus on the particular microfoam knot of interest will likely result in a change in the string's width cinch ratio (WCR) post-cinching.

_______________________________________________
FILES: script-testall, script-xtestall, script-check(mha|log)

These are bash script files that test all of the internal sample strings (a--h), the solved-tricksy, and Drosophila example corpus. The first script runs in default mode, while the second tests in both default and -x mode, which adds extra transition squeezing. The output.mha file in the top level directory corresponds to the output dfrom the script-xtestall run.

Other useful scripts are also included in the parent directory. Except for a compile script named giddyup, all others begin with "script-".
_______________________________________________
FILE: corpus/tricksy-output.mha

If this file is present it has record of newly identified "tricksy" strings with problematic cinching as identified by the associated version of maximal. Typically, this file is committed after a -XX Fischer-Yates run (script-RNDtest) and is named with additional version and run information.
_______________________________________________
