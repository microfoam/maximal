# maximal/waves/chowder/ README
The `chowder`<sup>*</sup> sub-directory contains sequence files currently under
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
some information in the new file name about how it was derived and/or how the
desired solution should be.

<sup>\*</sup> In surfing lingo, chowder is (annoying) flotsam and jetsam. If you eat the chowder, you are crashing into it or 
at the very least tumbling over it.
Outstanding (unsolved) chowder are prefaced with the header leader 'Chowder chunks:', while solved ones will have the header leader 'Solved:'.

## Chowder to eat

- [ ] seq-64-snippet-Rnx_symmetry.txt 
- [x] seq-29-toslp2-snippet.txt: keeping track of different *k*-mers annotated at same position
- [x] seq-146-v344_33-snippet.txt: adjacent non-zero mod k's
- [x] seq6-koslip-snippet.txt, evaluation of overlapping repeats
- [x] seq26-snippet.txt: split repetitions of overlapping repeat with a fractal subset.
- [x] seq-016-cyc3_a_knot -Rn, evaluation of overlapping repeats not triggered correctly
- [x] tubespit-17: strand symmetry (v4.30)
- [x] vnd_NEE_Dsech-snippet_1: skip fractal repeats in first unit so that cinch-t does not flatline from equivalence violation (v4.29)

## Chowder chunk: *seq-64-snippet-Rnx_symmetry.txt* 
A 9-mer repeat is skipped in the process of evaluating earlier overlapping repeats.
This does not occur when the reverse complement is run. This behavior is elicited
under '-x' squeeze (hence '-Rnx' in the file name). This cinching issue was 
discovered while addressing an unrelated transition calling issue with the same 
overlapping repeats (4-mer and 7-mer). As explained below, after working on this 
a bit additional aspects were revealed.

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

Okay, so this one actually got a bit more interesting for a couple of reasons.
First, canceling of non-fractal overlapping repeats is not impossible but
requires making sure to not cancel, or uncancelling previously overlapping
repeats that are no longer overlapping. In addressing this aspect of
overlapping repeats, which are now called in the program, I started to
appreciate that this particular sequence is a special case of a repeat unit
separated by an expanded mono-nucleotide tract. Such tracts are not called
until the cinch-k module. So instead of canceling the 7-mer in this case, I
concluded you actually want to cancel the 9-mer and leave it to cinch-d to 
cinch up the 9-mer as a 6-mer!

Take a look again, this time with some annotation added below the consensus row:

```
 2-D pass #2: cinch-t (width = 46)
    1. >CAAGTCGACGC/   11
    2.  .......ACGCTCC/   18
    3.  .......ACGCCCCCCACGCTCCCAG...
    4.  .........:.........:.........:
        _________|_________|_________|
                 10        20        30
        CAAGTCGACGCYCCCCACGCTCCCAG...
        .......1234567..1234567...... <- the 7-mer, indexed with a downstream third unit
        ..............123456789...... <- frame 1 of the 9-mer
        ...............123456789..... <- frame 2 of the 9-mer

 2-D pass #4: cinch-k (width = 30)
    1. >CA/    2
    2.  .AGTCGACGC/   11
    3.  ......ACGCTC/   17
    4.  .........:.C/   18
    5.  ......ACGCCC/   24
    6.  .........:.C/   25
    7.  .........:.C/   26
    8.  .........:.CACGCTC/   33
    9.  .........:.......C/   34
   10.  .........:.......CAG...
        _________|_________|_________|
                 10        20        30        
        CAGTCGACGCYCACGCTCAG...
        ......123456123456..... <- the 6-mer, indexed

 2-D pass #6: cinch-d (width = 24)
    1. >CA/    2
    2.  .AGTCGACGC/   11
    3.  ......ACGCTC/   17
    4.  .........:.C/   18
    5.  ......ACGCCC/   24
    6.  .........:.C/   25
    7.  .........:.C/   26
    8.  .........:.C/   27
    9.  ......ACGCTC/   33
   10.  .........:.C/   34
   11.  .........:.CAG...
        _________|_________|
                 10        20        
        CAGTCGACGCYCAG...
        ......123456..... <- the 6-mer, all cinched up
```
So *maximal* does this now, but it means that this still counts as chowder because we need to ensure strand symmetry. :thinking:

(Technical point: Incidentally, while working on this sequence I discovered that cinch-d had become hobbled because of an inverted
less-than sign that should have been a greater-than-or-equal sign in evaluating score threshold! I must have missed
it because I also recently increased the size of `PISO` to help prioritize chowder. Like the hobbling of cinch-d,
elevating `PISO` also has the effect of increasing the benchmark average WCR's.)

## Solved: *waves/chowder/seq-29-toslp2-snippet.txt*, keeping track of different *k*-mers annotated at same position

I will document this more fully later b/c it involved installing a more robust functionality for *k*-mers at the same position *n*,
which bears explaining.
For this chowder sequence snippet, position 15 is position *n* for the perfect CA-dinucleotide repeat *and* the imperfect 8-mer repeat, 
5'-CACCAGYA. There was more to it than that b/c of overlapping fractal and parent repeats, some of which were canceled for various
reasons, leaving orphan fractals. Also, involved conditioning one cinch-k loop.

```
 2-D pass #4: cinch-k (width = 11)

>GA/
 .ACAC/
 ..CA/
 ..CAC/
 ....CAG/
 ..CA/
 ..CAC/
 ....CAGTACG>
 _________|_
          10
 GACACAGTACG

```
-AE, 2/6/2020

## Solved: *seq-146-v344_33-snippet.txt*, adjacent non-zero mod *k*'s

I just solved this sequence snippet, which previously necessitated nudging correction. I didn't think this case was possible because 
I had never seen an example of adjacent nucleotides with repeats having non-zero mod *k*'s (*k*<sub>*n*</sub> mod *k* <sub>*n*-1</sub> > 0).

I leave the correct cinch-t solution here as a reverse exercise. Try to spot the non-zero mod *k*-mers. (If you want a hint,
factor in cycling frames.):

```
 2-D pass #2: cinch-t (width = 27)

    1. >GACCT/    5
    2.  ..CCTGACAGC/   14
    3.  .........GCCGCCCG/   22
    4.  .........:CTGCCCGAACAC/   34
    5.  .........:.........CACGAACG>   42
        _________|_________|_______
                 10        20        
        GACCTGACAGCYGCCCGAACACGAACG
```

I am adding a little more exposition here for the curious, but also to collate some notes from my hand-written scratch notebook and commit notebook.

Adjacent *k*-mers for the most part are often the result of cycling frames.
For example, for the sequence 5'-CATGATGAG *maximal* will mark two 3-mers at positions 5 and 6 as follows:
```
n: 123456789
  >CATGATGAG>
k: ....33...
r: ....11...
```

The adjacent 3-mers correspond to two related cycling frames as such:
```
Cycling frame 1:
  >CATG/
   .ATGAG>

Cycling frame 2:
  >CATGA/
   ..TGAG>
```

In my notebooks, I call these cycling islands and think of them as a set. With a higher number of repeats,
you can start to call higher *k* repeats at adjacent positions. So for the sequence
5'-CATGATGATGATGAG you get an annotation as follows:
```
  >CATGATGATGATGAG>
k: ....33366333...
r: ....33211111...
```

So we have adjacent columns with differently sized *k*-mers being called: *k*=3 stepping up to *k*=6, and then later *k*=6 stepping back down to *k*=3.
We can appreciate that this is an artifact of the underlying 3-mers being repeated enough times for higher *k* repeats being called 
that are a multiple of the smaller unit. Accordingly, we can identify these cycling islands by simply requiring
that adjacent *k*-mers have zero-valued mod values (*e.g.*, 6 mod 3 = 0).

Okay, so what about adjacent non-zero mod valued *k*-mers? Well at lower *k*-mer sizes they are not very interesting
because they are the result of mononucleotide tracts, which *maximal* is smart enough to detect.
To explain a bit more, let's imagine an adjacent 3-mer and 4-mer. I will write them with letters representing potentially
different symbols at each column.

***A 4-mer at *n* and a 3-mer at *n*-1:***
```
   m...n...
  >ABCDABCD>
k: ...34...
r: ...11...

This sequence can fold like this:
  >ABCD/
   ABCD>

And it can fold like this:
  >ABC/
   DABCD>
```
This implies that D = C = B = A. So the sequence is really just 5'-AAAAAAAA.

Because of the cinch-t sytem of annotating at *n*, I also wanted to work it out for the reverse order.

***3-mer at *n* and a 4-mer at *n*-1:***

```
   ..m..n..
  >CAABCABC>
k: ....43..
r: ....11..

This sequence can fold like this:
  >CAABC/
   ..ABC>

And it can fold like this:
  >CAAB/
   CABC>
```
This implies that C = B = A. So again the sequence is really just 5'-AAAAAAAA.
Natural *k*-mer tandem repeats have a distribution that is shifted to the smaller *k* sizes, so normally
you don't see "adjacent non-zero mod *k*'s". But they are possible at higher *k* sizes. Here's
an example for a 5-mer and a 3-mer in both orders:

***A 5-mer at *n*, and a 3-mer at *n*-1:***
```
   m....n....
  >ABCDEABCDE>
k: ....35....
r: ....11....

This sequence can fold like this:
  >ABCDE/
   ABCDE>

And it can fold like this:
  >ABCD/
    EABCDE>
```
This implies that C = A and E = D = B. So the sequence is actually 5'-ABABBABABB.
This means we *can* have non-trivial adjacent non-zero mod *k*'s (5 mod 3 = 2). But they are 
incompatible with each other and one or the other but not both could be taken.

Here's the case for the reverse order, which actually is a different sequence.

***A 3-mer at *n* and a 5-mer at *n*-1:***
```
   ...m..n..
  >CABABCABCD>
k: .....53...
r: .....11...

This sequence can fold like this:
  >CABABC/
   ...ABCD>

And it can fold like this:
  >CABAB/
   CABCD>
```
This implies that C = A and D = B. So the sequence is actually 5'-AABABAABAB.

:sweat_smile:

### Abba-Zabba foam
Now here's the reason I have been obsessed with this new fractal way of looking at repeats, 
inspired by studying sequence evolution in relatively unconstrained gene regulatory sequences. 
There is nothing about this particular business that is common sense. For example, having implemented
the capabilities I just described, I thought I had a choice about whether to choose the bigger
or the smaller of the two adjacently-annotated *k*-mers with non-zero mod *k*'s.

In a superficially greedy approach, you might think that you should take the larger of the two *k*-mers,
but what happens is usually messier than choosing the smaller of the two for reasons
I will not articulate here for brevity. Nonetheless, I explored this question by working on the 
two abstract sequences derived above for an adjacent 5-mer and 3-mer in both orders.  For reasons that will become clear, I will call them abba-zabba sequences.

Sequence pattern one is ABABBABABBZ.

Sequence pattern two is AABABAABABZ.

I added the 'Z' at the tail because originally I wanted to work on both of them in one sequence and the 'Z' was a good separator of the A's and B's.
It then turned out I would get different behavior depending on whether I repeated each sequence.
This had to do with whether the internal repeats were solved by cinch-t or cinch-k.
The fractal repeats would be solved by cinch-k only if there were two units. This exercise, 
and I am now certainly glad that I wrote this out for the README, showed me that an adjacent 5-mer and 3-mer
are always constructed of an underlying 2-mer pattern based on 'AB'.
The significance is that rather than choosing either the 3-mer or the 5-mer, one should choose neither *k*-mer! 
Why? Just see for yourself:

***Abba-zabba foam based on pattern one:***
```
2-D pass #4: cinch-k (width = 3)
 >AB/
  AB/
  .B/
  AB/
  AB/
  .BZ/
  AB/
  AB/
  .B/
  AB/
  AB/
  .BZ>
  ___
  ABZ
```

***Abba-zabba foam based on pattern two:***
```
2-D pass #4: cinch-k (width = 3)
 >A/
  AB/
  AB/
  A/
  AB/
  ABZ/
  A/
  AB/
  AB/
  A/
  AB/
  ABZ>
  ___
  ABZ
```

:astonished:


## Solved: *seq6-koslip-snippet.txt*, conditional mark_tela break

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

This sequence knot has now been solved by a smarter pre-cinch_t `mark_tela()` module that conditions the normal
break out of a loop scanning rows starting at a given column *n*.
The `mark_tela()` code block was exapted from a cinch-t block and like cinch-t, `mark_tela()` used to break out of 
2-D column *n* at row *m* as soon as it found a repeat, thus favoring higher *k*-mer size. The problem for `mark_tela()` functionality is that
a fractal repeat might be obscured by a cycling repeat frame of the parent TR containing the fractal. 
This is important because cinch-t needs to skip fractal repeats, which are then addressed by cinch-k.
Now *maximal*  senses whether this column is obscuring a smaller fractal *k*. 

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
I didn't know for a few days how I would eventually solve this, but I am pretty happy with the simplicity of this solution.

:smiley:

## Solved: *seq26-snippet.txt*, fractal splitting
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

:heart_eyes:


## Solved: *seq-016-cyc3_a_knot-Rn.txt*
In earlier program versions, this sequence was folded correctly but not optimal.
The `Rn` designates this was from the first test type in `cleanup_set_all`.

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

:slightly_smiling_face:

## Solved: *tubespit-17-symmetric.txt*, tricks in annotation
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
strip all of this formatting (*i.e.*, spaces, line breaks and non-essential characters)
so it does not interfere with the reading and processing of the file. Therefore, formatting 
can be used in a chowder file to indicate what the target folding pattern
is desired (but not yet acheived). The fasta header indicates this with 
the phrase `target_formatted_in_file`. The purpose of putting this in the fasta
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

:stuck_out_tongue_winking_eye:

## Solved: *vnd_NEE_Dsech-snippet_1.txt*, axiom testing
This addresses an issue resulting from conversion to axiom testing in the cinch-t/cinch-k system,
whereby cinch-t takes repeats in a moving window such that intra-repeat repeats only get taken
in the first unit of a longer repeat with cinch-k taking the remaining intra-TR repeats afterwards.
This is now solved for the first test case but will be kept in the chowder directory as a periodic test. The solution
was to have `mark_tela()` clear small fractal repeats in the first unit of a longer repeat.
This allows the longer repeat (9-mer here) to be cinched without violating equivalence prior to
cinch-k. The 5'-AC dinucleotide repeat works here as a sentinel repeat of cinch-t flatlining. 

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

:nerd_face:
 
*Last updated*: 2/6/2020 AJE
