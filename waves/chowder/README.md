# maximal/waves/chowder/ README
The `chowder` sub-directory contains sequence files currently under
attention or slated for attention during development of the 'maximal' program.
Because the benchmark scripts (typically base named `surfboard-cleanup_set_`) use
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
`/waves` sub-directories will also be stored in `/waves/chowder` sometimes with 
some information in the new file name about how it was derived and/or what it 
intended to be done.

## Chowder to eat

- [ ] seq-016-cyc3_a_knot -Rn
- [x] tubespit-17: strand symmetry (v4.30)
- [x] vnd_NEE_Dsech-snippet_1: skip fractal repeats in first unit so that cinch-t does not flatline from equivalence violation (v4.29)

## seq-016-cyc3_a_knot-Rn.txt
This sequence is folded correctly but not optimal.
The `Rn` designates this was from the first test type in `cleanup_set_all'.

```
>target_formatted_in_file
 ATGGGCTGCA
 ........CAAC
 .........AACAAT
 .........AACAATACG
 ...............ACGA
```

What currently happens:
```
2-D PASS #2: cinch-t 

   >ATGGGCTGCA/
    ........CAA/
    ........CAA/
    ........CAATAACAATACG/
    .........:........ACGA>
    _________|_________|__
             10        20        
```

## tubespit-17-symmetric.txt
This is a snippet of `tubespit/seq-017-cycle5.txt` and its reverse complement
separated by four N's, which work as separators to prevent cinching across the
segments. This combination of tricks allows a single sequence to be processed to display how both the
(+) and (-) strands are cinched. 
In addtion there are sentinel dinucleotide repeats before and after 
each of the two main repeat knots. These dinculeotides repeats can show whether
cinching continued after the knots whenever this sequence is run through
maximal for either strand. 

This sequence also demonstrates how sequences
can be saved in a hand-aligned fashion in the text file. The program will
strip all of this formatting (i.e., spaces, line breaks and non-essential characters)
so it does not interfere with the reading and processing of the file. Therefore, formatting 
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

Solved in v4.30:
```
   >GC/
    GCnnGACAC/
    ....GACAC/
    ......CACnnTG/
    .........:.TGnnnnCA/
    .........:.......CAnnGTG/
    .........:.........:.GTGTC/
    .........:.........:.GTGTCnnGC/
    .........:.........:........GC>
    _________|_________|_________|
             10        20        30     
```

Compare to variant missing the 3-mer overlapping repeat:
```
    1. >GTGC/    4
    2.  ..GCAGCTG/   11
    3.  ........GA/   13
    4.  .......TGAC/   17
    5.  .........AC/   19
    6.  ........GAC/   22
    7.  .........AC/   24
    8.  .........:CAnTATG/   31
    9.  .........:.....TGCAnnTG>   39
        _________|_________|___
                 10        20        
        GTGCAGCTGACAnTATGCAnnTG
```

## vnd_NEE_Dsech-snippet_1.txt 
This addresses an issue resulting from conversion to axiom testing in the cinch-t/cinch-k system,
whereby cinch-t takes repeats in a moving window such that intra-repeat repeats only get taken
in the first unit of a longer repeat with cinch-k taking the remaining intra-TR repeats afterwards.
This is now solved for the first test case but will be kept in the chowder directory as a periodic test. The solution
was to have `mark_tela()` clear small fractal repeats in the first unit of a longer repeat.
This allows the longer repeat (9-mer here) to be cinched without violating equivalence prior to
cinch-k. The 5'-AC dinucleotide repeat works here as a sentinel repeat of cinch-t flatlining. 
Some additional aspects remain to be handled.

```
2-D PASS #2: cinch-t [LONGEST *tANDEM REPEATS, k>1] (width = 39)

   >TAGCTCCTTAAATTAGCCAAGCGCGCA/
    .........:........AAGCGCGCAAGTACAGGAC/
    .........:.........:.........:.....ACAG>
    _________|_________|_________|_________
             10        20        30        

2-D PASS #4: cinch-k [INTRA-REPEAT *k-MERS > 0] (width = 35)

   >TAGCTCCTTAAATTAGCCAAGC/
    .........:.........:GC/
    .........:.........:GCA/
    .........:........AAGC/
    .........:.........:GC/
    .........:.........:GCAAGTACAGGAC/
    .........:.........:.........:.ACAG>
    _________|_________|_________|_____
             10        20        30        
```
 
*Last updated*: 1/27/2010 AJE
