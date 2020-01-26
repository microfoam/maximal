# maximal/waves/chowder/ README

This 'chowder' sub-directory is used to hold sequence files currently under
attention or slated for attention during development of the 'maximal' program.
Because the benchmark scripts (typically named 'surfboard-cleanup_set_*') use
sequences in the 'waves/tubespit' sub-directory, some chowder sequences will be
located in duplicate forms in both sub-directories with the same file name. In those
cases, the duplicate file in the chowder sub-directory indicates that the
sequence contained in the file is either no longer solved correctly or else is
not solved optimally even if solved correctly. If a sequece has never been
solved correctly it should be located in the 'chowder' sub-directory even if it
is referenced in a 
testing script. 

Often an interesting sequence problem is whittled down to a "knot", a snippet of 
a fuller sequence that isolates the interesting problem under investigation. In 
some of these cases, a duplicate sub-sequence to one maintained in the other 
'waves' sub-directories will be also stored in 'chowder' with some knot desigation
in the file name.


## tubespit-17-symmetric.txt
This is a snippet of `tubespit/seq-017-cycle5.txt` and its reverse complement
separated by four N's, which work as separators to prevent cinching across the
segments. This allows a single sequence to be processed to display how both the
+ and - strands are currently cinched. 
In addtion there are sentinel dinucleotide repeats before and after 
each of the two main repeat knots. These dinculeotides repeats can show whether
cinching continued after the knots whenever this sequence is run through
maximal for either strand. 

This sequence also demonstrates how sequences
can be saved in a pseudo-aligned fashion in the text file. The program will
strip all of this formatting (i.e., spaces, line breaks and non-essential characters)
so it does not interfere with the reading of the file. Therefore, formatting 
can be used in a chowder file to indicate what the target folding pattern
is desired (but not yet acheived). The fasta header indicates this with 
the phrase `target_formatted_in_file'. The purpose of putting this in the fasta
header is so that it is read in as the header and then displayed to the 
user/developer to indicate that one can look at the raw file to see what is desired.

```
 >target_formatted_in_file
  GC
  GCnnGACAC
  ....GACAC
  ......CACnnTG
  ...........TGnnnnCA
  .................CAnnGTG
  .....................GTGTC
  .....................GTGTCnnGC
  ............................GC
```

What currently happens:
```
 >GC/
  GCnnGACAC/
  ....GACAC/
  ......CACnnTG/
  ...........TGnnnnCA/
  .................CAnnGTG/
  .....................GTGTC/
  .......................GT/      <=== 
  .......................GTCnnGC/
  ............................GC>
```
 
*Last updated*: 1/26/2010 AJE
