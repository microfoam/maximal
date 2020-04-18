The file 'foam_and_chowder.log' is the standard output file for program runs that 
finish with <1000 for either pass quality or recovery of perfect 1-D 
sequence from the 2-D alignment. Most recently, 'maximal' has been set to also
write runs requiring nudging, now that maximal is rarely encountering sequences
it cannot handle. We consider nudging to be a stop-gap tool to recover from errors
in 2-D alignment columns, but not the preferred solution. The preferred solution 
is to not produce the errors in the first place. This is why they are now considered
minor chowder bits to work on.

The file 'foam_and_chowder.ragdoll' is a collation of the unsorted chowder produced
after running the ragdoll script ('surfboard_vnd_getworked_rgdll').

-AJE April 18th, 2020 (COVID-bound in the sun room).


F.A.Q. for the 'maximal/waves/' directory 

********************************************************
QUESTION: What is the purpose of the 'waves/' directory?

ANSWER: The 'waves' directory contains sequences that are used in various
testing routines of the 'maximal' program.  These sequences are essential to 
development of 'maximal', which is an implementation of maximal homology
alignment (see PREPRINT in the 'waves/' directory).

The sequences are stored and organized in three different subdirectories within
'waves/': 'animals/', 'chowder/', and 'tubespit/'. It also contains the 
'foam_and_chowder.log' file, which records information about tricky sequence 
runs.


********************************************************************
QUESTION: What distinguishes the sequences in the different 'waves/'
sub-directories? For example, what is the 'waves/animals/' sub-directory?

ANSWER: The original impetus for the 'maximal' program was the modeling of
replication slippage in evolving regulatory DNA enhancer sequences.  These
sequences are trickier than usual and are represented by 'shark_' files in the
'animals/' sub-directory.  In contrast, there are also control protein-coding
exon sequences, which are under different and more stringent evolutionary
constraints and serve as control comparison sequences to the "shark" sequences.
These are also stored in the 'animals/' sub-directories with the 'seal_' base
file name.

Many surfers have enjoyed watching seals bodysurf the waves, but both surfers
and seals have to look out for the surfers in gray suits!


*****************************************************************************
QUESTION: Okay, but what about the 'chowder/' and 'tubespit/' subdirectories?

ANSWER: The 'chowder/chowder_X.txt' files correspond to sequences that happen
to contain problematic sequences that are not completely cinched by 'maximal'.
These sequences are identified by conducting hundreds of thousands of 
'maximal' runs using the '-XX' Fisher-Yates randomization option and a 
selection of the sequences in 'waves/animals/' according to the
get_worked_rgdoll script. These sequences are used to improve the 'maximal' 
code base. 

Once solved, chowder sequences are sequences are renamed as tubespit sequences
and moved to the 'tubespit/' sub-directory. Tubespit sequences are also used in
the 'cleanup_set_glassy|gnarly' scripts and are thereby used to ensure progress
in future code upgrades.

Surfers do not like to encounter chowder bits while getting worked in the waves
after a fall! A surfer would rather experience some ocean spray coming out of a 
totally tubular wave!

*******************************************************************************
QUESTION: Okay, I am starting to understand how all this works, but why go to 
all the trouble of fitting this surfing theme for a program on maximal homology
alignment?

ANSWER: The first answer is that it helps to make the 'maximal' code development
more fun. Secondly, and relatedly, this is part of gamefying the code 
development. In this code development game, you are a surfer riding sets of
waves and you want to score more points by being able to use your different
surfboards (the various scripts) to successfully ride different waves in
different styles. You get more points by turning waves into exhilirating
tubespit experiences. You get less points by falling and getting worked like 
a ragdoll underwater with all the chowder. Last, the theme helps to organize the
physical and mental organization of the code, the scripts, and the test data
sets. This makes the code development process easier to navigate and grok!

-AJE (March 4th, 2019)

