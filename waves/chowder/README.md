# maximal/waves/chowder/ README
The `chowder` sub-directory contains sequence files currently under
attention or slated for attention during development of the *maximal* program.
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

- [ ] seq-64-snippet-Rnx_symmetry.txt 
- [x] seq6-koslip-snippet.txt, evaluation of overlapping repeats
- [x] seq26-snippet.txt: split repetitions of overlapping repeat with a fractal subset.
- [x] seq-016-cyc3_a_knot -Rn, evaluation of overlapping repeats not triggered correctly
- [x] tubespit-17: strand symmetry (v4.30)
- [x] vnd_NEE_Dsech-snippet_1: skip fractal repeats in first unit so that cinch-t does not flatline from equivalence violation (v4.29)

## Chowder chunk: seq-64-snippet-Rnx_symmetry.txt 
A 9-mer repeat is skipped in the process of evaluating earlier overlapping repeats.
This does not occur when the reverse complement is run. This behavior is elicited
under '-x' sequeeze (hence '-Rnx' in the file name). This cinching issue was 
discovered while addressing an unrelated transition calling issue with the same 
overlapping repeats (4-mer and 7-mer).

Top strand comes out like this:
```
2-D PASS #2: cinch-t (width = 46)

   >CAAGTCGACGC/
    .......ACGCTCC/
    .......ACGCCCCCCACGCTCCCAGTATTTTTCCCTAACCAGGGC/
    .........:.........:.........:.........:....GC>
    _________|_________|_________|_________|______
             10        20        30        40        
    CAAGTCGACGCYCCCCACGCTCCCAGTATTTTTCCCTAACCAGGGC

```

But the top strand should actually look like the following, which would match the reverse complement:
```
 >CAAGTCGACGC/
    .......ACGCTCCACGCCCCC/
    .........:...CACGCTCCCAGTATTTTTCCCTAACCAGGGC/
    .........:.........:.......:.........:....GC>
    _________|_________|_________|_________|____
             10        20        30        40    
    CAAGTCGACGCTCCACGCTCCCAGTATTTTTCCCTAACCAGGGC
```

Reverse-complement comes out like this:
```
2-D PASS #2: cinch-t (width = 44)

   >GC/
    GCCCTGGTTAGGGAAAAATACTGGGAGCGTG/
    .........:.........:..GGGGGCGTGGAGCGT/
    .........:.........:.........:...GCGTCGACTTG>
    _________|_________|_________|_________|____
             10        20        30        40        
    GCCCTGGTTAGGGAAAAATACTGGGRGCGTGGAGCGTCGACTTG

```

## Solved: seq6-koslip-snippet.txt (conditional mark_tela break)

Previously cinch-t runs produced an auto-alignment that was correct but not optimal:
```
  >TAGC/    4
   ..GCGTC/    9
   ..GCGTCGTGCACTGAT/   24
   .........:.....ATC>   27
   _________|________
            10        
   TAGCGTCGTGCACTGATC
```

Target fold is the following (pre-cinch-k):
```
  >TAGCG/
   ...CGCGT/
   ...CGCGT/
   .....CGTGCACTGAT/
   ..............ATC>
```

This sequence knot has now been solved by a smarter pre-cinch_t mark_tela() module that conditions the normal
break out of a loop scanning rows starting at a given column *n*.
The mark_tela() code block was exapted from a cinch-t block and like cinch-t, mark_tela() used to break out of 
2-D column *n* at row *m* as soon as it found a repeat, thus favoring higher *k*-mer size. The problem for mark_tela() functionality is that
a fractal repeat might be obscured by a cycling repeat frame of the parent TR containing the fractal. 
This is important because cinch-t needs to skip fractal repeats, which are then addressed by cinch-k.
Now mark_tela() senses whether this column is obscuring a smaller fractal *k*. I didn't know for a few days
how I would eventually solve this, but I am pretty happy with the simplicity of this solution.
Post-cinch-k the 2-D auto-alignment now looks like this:

```
 2-D pass #4: cinch-k (width = 15)

   >TAGCG/
    ...CGT/
    ...CG/
    ...CGT/
    ...CGTGCACTGAT/
    .........:..ATC>
    _________|_____
             10        
    TAGCGTGCACTGATC
```

## Solved: seq26-snippet.txt
This was addressed in a recent v4.30 commit and the solution was pretty cool.
Basically, this one involved splitting out a subset of repeats that was also
being counted as fractal. So the (GA)x2 repeat gets cinched as (GA)x1 in 
cinch-t and another (GA)x1 in cinch-k.

```
2-D PASS #2: cinch-t (width = 56)

   >ACTGATGCAAAACAAATGCAG/ **
    .........:........CAGGCGA/
    .........:.........:...GAGAT/
    .........:.........:...GAGATAAATGGGACGAGCGGTGCATCCG/
    .........:.........:.........:.........:.........CGGGTGC>
    _________|_________|_________|_________|_________|______
             10        20        30        40        50       

 Next: cinch-k for k = 2...

2-D PASS #4: cinch-k (width = 54)

   >ACTGATGCAAAACAAATGCAG/
    .........:........CAGGCGA/
    .........:.........:...GA/
    .........:.........:...GAT/
    .........:.........:...GA/
    .........:.........:...GATAAATGGGACGAGCGGTGCATCCG/
    .........:.........:.........:.........:.......CGGGTGC>
    _________|_________|_________|_________|_________|____
             10        20        30        40        50        
```

## Solved: seq-016-cyc3_a_knot-Rn.txt
This sequence was folded correctly but not optimal (in earlier versions now).
The `Rn` designates this was from the first test type in `cleanup_set_all'.

```
>target_formatted_in_file
 ATGGGCTGCA
 ........CAAC
 .........AACAAT
 .........AACAATACG
 ...............ACGA
```

Versions 4.29, and early 4.30 produced:
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

Later versions of 4.30 produce a 2-D alignment that is nearly at target, short
of one cyclization. This was an issue of not triggering a comparison of
overlapping repeats in cinch-t.

```
2-D PASS #2: cinch-t 

    1. >ATGGGCTGCACA/   12
    2.  .........ACA/   15
    3.  .........ACAAT/   20
    4.  ........AACAATACG/   29
    5.  .........:....ACGA>   33
        _________|________
                 10        
        ATGGGCTGCACAATACGA
        ........A.........
```

This issue is now solved (latest v4.30 version):
```
2-D PASS #2: cinch-t 

    1. >ATGGGCTGCACA/   12
    2.  .........ACA/   15
    3.  .........ACAATA/   21
    4.  .........ACAATACG/   29
    5.  .........:....ACGA>   33
        _________|________
                 10        
        ATGGGCTGCACAATACGA
```


## Solved: tubespit-17-symmetric.txt
This is a snippet of `tubespit/seq-017-cycle5.txt` and its reverse complement
separated by four N's, which work as separators to prevent cinching across the
segments. This combination of tricks allows a single sequence to be processed to display how both the
(+) and (-) strands are cinched. 
In addtion there are sentinel dinucleotide repeats before and after 
each of the two main repeat knots. These dinculeotides repeats can show whether
cinching continued after the knots whenever this sequence is run through
*maximal* for either strand. 

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

## Solved: vnd_NEE_Dsech-snippet_1.txt 
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
 
*Last updated*: 1/31/2020 AJE
