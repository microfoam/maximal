# maximal
maximal is a program for maximal homology alignment (MHA), originally created and written by Albert Erives.
Additional help options can be seen by running ./maximal without any options specified.

Below is a guide to key files in the maximal software project and their naming schemes.
_______________________________________________
FILE: maximal, m272 ("272 is VERSION-SPECIFIC)

This is an Apple Mac executable (for now) of the gcc compiled code. Previous versions of maximal that are useful for comparing runs are named "m#").

_______________________________________________
FILE: mha_v2.72.c  ("2.72" is VERSION-SPECIFIC)

This is the main C programming file. Example name shows this is version 2.72.

_______________________________________________
FILE: output.mha

This is the standard ouput report to which maximal writes. 
Every run of maximal gets added in a line at the bottom of this file.
For each version of maximal, the intent is to run it on a test set of strings.
Thus, this file will reflect the output performance of the latest version of maximal.

_______________________________________________
FILE: solved-##-descrp_tricksy.txt ("##-descrp" is VARIABLE)

Files named in this format represent short sequence text files with a tricksy microfoam nature. Typically, tricsky strings were at one point problem strings that were problematically cinched by maximal. After developing the code to conduct MHA successfully on a new tricksy string it is renamed as "solved", given a number, and a description. These file names and sequences should rarely be modified because they are used for historical comparisons of different versions of maximal. For example, even changing the length of the sequence in order to focus on the particular microfoam knot of interest will likely result in a change in the string's width cinch ratio (WCR) post-cinching.

_______________________________________________
FILE: Dxxxx_ex3.txt, Dxxxx_NEE.txt ("xxxx" is a Drosophila species name)

Files named in this format represent vnd exon 3 sequences or NEE sequences from one of several species of flies. These file names and sequences should rarely be modified because they are used for historical comparisons of different versions of maximal. For example, even changing the length of the sequence in order to focus on the particular microfoam knot of interest will likely result in a change in the string's width cinch ratio (WCR) post-cinching.

_______________________________________________
FILE: script-testall

This is a bash script file designed to test all of the mha internal sample strings (a--h), the solved-tricksy, and Drosophila example corpus in both perfect repeat mode and imperfect repeat mode with transitions ("-x" option). Running this script produces (currently 156 tests) reported in the associated output.mha file.

_______________________________________________
FILE: tricksy-output.mha

If this file is present it has record of newly identified "tricksy" strings with problematic cinching as identified by the associated version of maximal.
_______________________________________________
This readme file was last updated on April 8th, 2018 by AJE.
