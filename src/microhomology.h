/******************************************************************/
/* lineup's microhomology.h: header file since mha_v3.97.c       */
/******************************************************************/

#ifndef FILE_LOGY_SEEN
#define FILE_LOGY_SEEN

#define TEMP		2		/* SETS MAGIC MELT TEMP; SETS LOWER BOUND FOR MELTAGE opt_m; UNELECTED, DEFAULT ACTS AS TEMP=1 */
#define CYCMAX     60       /* SEMI-MAGIC NUMBER; SEARCH MAGIC TO FIND OTHER EMBEDDED DECISIONS */
#define MAXROW  20000       /* maximum input line size; TAGGED: <MAGIC> BASED ON WHAT MY 'PUTERS CAN DO AS CURRENTLY WRITTEN */
#define WIDTH      72       /* BANDWIDTH: MAX WIDTH OF HEMIDIAGONAL OF PATHBOX; MAX TR UNIT SIZE */
#define MEMROWS    20       /* NUMBER OF mem[MEMROWS] ROWS IN STRUCT COORD ARRAY TELA */
							/* USE 1: BIT (0/1) VALUES FOR MARK_TELA MARKS ASSOCIATED WITH A SINGLE LOOP OF CLEAR_ALL PRECEDENCE */
							/* USE 2: NUMBER OF AVAILABLE ROWS FOR STORING OVERLAPPING REPEAT FRAMES; MULT. OF 4 - EXTRA */
#define START       0       /* FOR USE WITH 1ine_end() */
#define END         1       /* FOR USE WITH 1ine_end() */
#define RULER       2       /* FOR USE WITH 1ine_end() */
#define SLIPS       3       /* FOR USE WITH 1ine_end() */
#define PATHBOXHEAD 4       /* FOR USE WITH 1ine_end() */
#define BLOCKHEAD   5       /* FOR USE WITH 1ine_end() */
#define SLIPRULER   6       /* FOR USE WITH 1ine_end() */
#define XDIR        0		/* FOR USE WITH push_gPnt() */
#define YDIR        1		/* FOR USE WITH push_gPnt() */
#define SEQHEADS  160		/* Maximum size for sequence headers, FASTA */
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })

struct genPoint {
	short unsigned int rel_xy;	/* Point positions relative to previous point:         */
								/* Point (x,y) can only ever be (1,0) or (0,1).        */
								/* rel_xy = 0 if (1,0); relx_xy = 1 if (0,1).          */
	int topPar;		/* Top-most paralogy position, can be same as previous */
	int prevPar;	/* Previous (closest on the left) paralogy position    */
};

struct coord {
	int x;			/* 2D x-AXIS COORDINATE => COLUMN   */
	int y;			/* 2D y-AXIS COORDINATE => ROW      */
	int X;			/* LOCATIONS OF LEFT-MOST OVERLAPPING TRs */
	char c;			/* CHARACTER LETTER IN SEQUENCE: ASSIGN ONLY ONCE! */
	/*************************************************************************************************/
	int k;			/* k-MERs BY LOCATION; WAS A SLIPLOC_NMER I USED TO KNOW */
	int r;			/* REPEAT NUMBER ALBERT-STYLE: 2nd UNIT OF TR = 1st REPEAT; r+1=TOTAL # OF UNITS */
	int cyc_l;		/* number of frames can cycle through */
	char e;			/* EQUIVALENCE CLASS LETTER: ASSIGN ONLY ONCE! */
	/*************************************************************************************************/
	int cyc_k;		/* k-mer size indicated for each cyclelizable option */
	int cyc_r;		/* repeat number for each frame of cyclelizable option */
	int cyc_P;		/* product of cyc_k * cyc_r for all cycleling options */
	char echoes;	/* PICTOGRAPHIC SLIPLOC_ECHOES, AKA DNA REVERB! */
	/*************************************************************************************************/
	int cyc_S;		/* sum of compatible products for all cycleling options */
	int o;			/* CYCLE LENGTH; FOR CYCLELIZATION FUNCTIONS */
	int Dtr;		/* CUMULATIVE SUM OF DIAGONAL TANDEM REPEAT (DTR) SCORES BY POSITION */
	char t;			/* IUPAC TRANSITIONS IN DNA USUALLY (RY) IN "IMPERFECT" TANDEM REPEATS */
	/*************************************************************************************************/
	int  mem[MEMROWS];	/* Use 1: bits for clearall_tela() in mark_tela(), so that order of precedence can be easily edited */
						/* Use 2: cycling frames; count-off column positions per unit; 32 - 11 = 21 */
						/* one row/frame; row 0 is row # locator; MEMROWS IS BASED ON MEM AL. */
	int ok;			/* ALL SERIES: PRE-CINCH-T: OPERATIONAL k-MER SIZE (prev. all_k)*/
	int or;			/* ALL SERIES: PRE-CINCH-T: REPEAT NUMBER (prev. all_r)       	*/
	int k0;			/* ALL SERIES: PRE-CINCH-T: k-MER; largest k at pos. n, undone	*/
	int k1;			/* ALL SERIES: PRE-CINCH-T: k-MER; largest k at pos. n        	*/
	int k2;			/* ALL SERIES: PRE-CINCH-T: k-MER; smaller k at pos. n        	*/
	int impk;		/* ALL SERIES: PRE-CINCH-T: imperfect k-mer not called for k1	*/
	int all_S;		/* ALL SERIES: PRE-CINCH-T: SUM OF SCORES OVER ALL UNITS      	*/
	int isl;		/* MARKS ARCHIPELAGO OF INDEPENDENT ISLANDS               		*/
	int all_Z;		/* ALL SERIES: PRE-CINCH-T: ALL_S + TIE-BREAKERS			  	*/
	int all_R;		/* ALL SERIES: PRE-CINCH-T: POSITION OF CONFLICTING TR ON RHS 	*/
	int all_L;		/* ALL SERIES: PRE-CINCH-T: POSITION OF CONFLICTING TR ON LHS 	*/
	char stat;		/* ALL SERIES: PRE-CINCH-T: STATUS                            	*/
	char statf;		/* ALL SERIES: PRE-CINCH-T: STATUS - FRACTAL TR               	*/
	char statl;		/* ALL SERIES: PRE-CINCH-T: STATUS - LOW-COMPLEXITY           	*/
	/*************************************************************************************************/
	int cyc_Lf;		/* Left-side overlapping TR */
	int cyc_Rt;		/* Right-side overlapping TR */
	char cyc_o;		/* x => cinched; o => untaken cyclelizable option; !,** => CHECK_TELA VIOLATIONS */
	/*************************************************************************************************/
	struct genPoint gPnt;	/* Will eventually replace absolute points above */
} tela[MAXROW];


/* This unnamed struct type organizes run command options read in main() */
struct {
	int		bit;				/* O-N / O-F-F bit switch   */
	int 	val;				/* incrementable value  	*/
	char	sym;				/* option letter symbol 	*/
	char	description[64];	/* usage description    	*/

}	*Options[53] = {},
                           /*.........|.........|.........|.........|.........|.........|...X*/
	opt_a = {0, 2, 'a', 	"Cinch-k mode 2: all k (default); 1: k=1 only; 0: skip all k.*3 "},
	opt_b = {0, 3, 'b', 	"Change default k-floor; starts transition matching above this.*"},
	opt_c = {0, 0, 'c', 	"Show base 62 single character code used for k-size and number. "},
	opt_d = {1, 1, 'd', 	"Cinch-d ON (1, default), OFF (0), or w/ reversed k-loop (2).*2 "},
	opt_e = {0, 0, 'e', {0}},
	opt_f = {0, 0, 'f', 	"Show foam-free segments (requires allowing default relax-2D.   "},
	opt_g = {0, 0, 'g', 	"Gel-up (counteract base level melt setting in cinch-d).*       "},
	opt_h = {0, 0, 'h', 	"Show help and usage."},
	opt_i = {0, 0, 'i', {0}},
	opt_j = {0, 0, 'j', {0}},
	opt_k = {0, 0, 'k', 	"Show cinch-t k-mer count.                                      "},
	opt_l = {0, 0, 'l', 	"Show cinch-t slip locations in 1-D sequence (in base 62).      "},
	opt_m = {0, 0, 'm', 	"Split, open, and melt for cinch-d (see also -g option).*       "},
	opt_n = {0, 0, 'n', 	"Do not do relax-2D pass.                                       "},
	opt_o = {0, 0, 'o', 	"Show original (pre-formatted) input sequence.                  "},
	opt_p = {0, 0, 'p', 	"Show program run parameters.                                   "},
	opt_q = {0, 0, 'q', {0}},
	opt_r = {0, 0, 'r', 	"Show row numbers.                                              "},
	opt_s = {0, 0, 's', 	"Silence writing to normal output file ('Surf_wavereport.mha'). "},
	opt_t = {0, 0, 't', 	"Skip cinch-t module and proceed to remaining modules.          "},
	opt_u = {0, 0, 'u', 	"Unwrap display of 2-D output (adds 10 columns to wrap).*       "},
	opt_v = {0, 0, 'v', 	"Increment verbosity level.*3                                   "},
	opt_w = {0, 0, 'w', {0}},
	opt_x = {0, 0, 'x', 	"Reduce thresholds for transition matching by one.              "},
	opt_y = {0, 0, 'y', {0}},
	opt_z = {0, 0, 'z', 	"Set mismatch score from default to zero (visible in Pathbox).  "},
	opt_A = {0, 0, 'A', {0}},
	opt_B = {0, 0, 'B', 	"Reset blank character (default '.') to space (' ').*4          "},
	opt_C = {0, 0, 'C', 	"Use the reverse complement (- strand) or reverse as allowed.   "},
	opt_D = {0, 0, 'D', 	"Turn on program development prompts (user-controlled pauses).*7"},
	opt_E = {0, 0, 'E', {0}},
	opt_F = {0, 0, 'F', 	"Fill recessed line ends with spacer marks as a visual aid.     "},
	opt_G = {0, 0, 'G', {0}},
	opt_H = {0, 0, 'H', 	"Show help and usage.                                           "},
	opt_I = {0, 0, 'I', {0}},
	opt_J = {0, 0, 'J', {0}},
	opt_K = {0, 0, 'K', 	"Show consensus row of 2-D alignment.                           "},
	opt_L = {0, 0, 'L', 	"Show 1-D locations at line ends.                               "},
	opt_M = {0,10, 'M', 	"Double the default wrap length for long homomonomer tracts.*   "},
	opt_N = {0, 0, 'N', {0}},
	opt_O = {0, 0, 'O', 	"Append 2-D alignment to consensus file or raw input file.*2    "},
	opt_P = {0, 0, 'P', 	"Show Pathbox grid.                                             "},
	opt_Q = {0, 0, 'Q', {0}},
	opt_R = {0, 0, 'R', 	"Recover sequence from 2-D.                                     "},
	opt_S = {0, 0, 'S', 	"Take a number argument that is added to random seed generator.*"},
	opt_T = {0, 0, 'T', 	"Show diagonal threshold values for transitions scoring.        "},
	opt_U = {0, 0, 'U', {0}},
	opt_V = {0, 0, 'V', {0}},
	opt_W = {0, 0, 'W', {0}},
	opt_X = {0, 0, 'X', 	"Scramble sequence with basic or Fisher-Yates randomization.*2  "},
	opt_Y = {0, 0, 'Y', 	"Set Fisher-Yates length (FY_size) specified by run argument.   "},
	opt_Z = {0, 0, 'Z', {0}};
                           /*.........|.........|.........|.........|.........|.........|...X*/


/* This unnamed struct type organizes deep program settings with pleiotropic effects. */
struct {
	int 	set;				/* program setting       */
	char	units[8];			/* parameter units       */
	char	description[64];	/* parameter description */
	char	runopt[8];			/* related run option    */

}	*Parms[2] = {},
                   /*.......X  .........|.........|.........|.........|.........|.........|...X*/
	par_piso = {  5, "char",    "Sets the floor for transition scoring above this k-mer size.  ", "opt_x"},
	par_wrap = {100, "columns", "Sets screen wrap length for 2-D alignment.                    ", "opt_u"};


/* This unnamed struct type organizes a set of pointers and names for mha's 2-D alignment character set.               */
/* Note that much of lineup's legacy code relies on isalpha(), so many symbols cannot be changed to alpha characters. */
struct {
	char	sym;	/* character symbol      */
	int		cod;	/* decimal unicode value */
}
	ambig    = {'n',110},		/* MHA ambiguous character in DNA              */
	nucl_T   = {'T', 84},		/* MHA nucleobase thymine                      */
	nucl_U   = {'U', 85},		/* MHA nucleobase uracil                       */
	*PyrTU   = &nucl_T,			/* pointer to thymine or uracil                */
	fill_1   = {'.', 46},		/* MHA fill character (default): full stop     */
	fill_0   = {' ', 32},		/* MHA fill character (option): space          */
	fill10   = {':', 58},		/* MHA fill character (10 bp): colon           */ 	/* OTHER POSSIBILITIES: |, ^ */
	*Fill    = &fill_1,			/* pointer to default fill character           */	/* ALWAYS USE POINTER, WHICH CAN BE REDIRECTED */
	margin   = {' ', 32},		/* MHA margin character with pairwise: space   */
	*Margin  = &margin,			/* pointer to zero column margin under start   */
	fastahead= {'>', 62},		/* MHA fasta header to decouple from term.     */
	term00   = {'>', 62},		/* MHA terminator character: greater-than sign */	/* NEEDS TO BE UNIQUE */
	term01   = {'*', 42},		/* MHA terminator character: greater-than sign */	/* NEEDS TO BE UNIQUE */
	*Term    = &term00,			/* pointer to terminator character             */	/* ALWAYS USE POINTER, WHICH CAN BE REDIRECTED */
	monoL    = {'(', 40},		/* MHA long monomeric tract left delimiter     */
	monoR    = {')', 41},		/* MHA long monomeric tract right delimiter    */
	F_str    = {'+', 43},		/* MHA dsDNA strand indicator, forward         */
	R_str    = {'-', 45},		/* MHA dsDNA strand indicator, reverse         */
	*Strand  = &F_str,			/* pointer to forward strand character         */
	gap      = {'-', 45},		/* MHA null character, dash                    */
	slip     = {'/', 47},		/* MHA slip character: forward slash           */
	tick     = {'|',124},		/* MHA ruler tick mark: vertical bar           */
	*Tick    = &tick,			/* pointer to default 10 bp ruler tick mark    */
	*ZTick   = &fill_0,			/* pointer to default column zero tick mark    */
	cyc_take = {'x',120},		/* MHA cinch-t cycle frame evaluations         */
	cyc_skip = {'o',111},		/* MHA cinch-t cycle frame evaluations         */
	st_parent= {'p',112},		/* MHA mark_tela() status mark                 */
	st_fract = {'f',102},		/* MHA mark_tela() status mark, confirmed fr.  */
	st_Fract = {'F', 70},		/* MHA mark_tela() status mark, orphan Fractal */
	st_lowcm = {'L', 76},		/* MHA mark_tela() status mark, low-complexity */
	st_cycle = {'c', 99};		/* MHA mark_tela() status mark                 */
//	st_clash = {'!', 33},		/* MHA mark_tela() status mark                 */
//	st_skip0 = {'_', 95},		/* MHA mark_tela() status mark                 */
//	st_skip1 = {'-', 45},		/* MHA mark_tela() status mark                 */
//	st_skip2 = {'~',126};		/* MHA mark_tela() status mark                 */

struct cinch {
	int pass_Q;					/* Pass quality score                                               */
	int pass_R;					/* Pass runs count; for Current this will hold cumulative bad slips */
	int pass_V;					/* Pass runs value; pass counter                                    */
	int pass_W;					/* Pass 2-D width                                                   */
	int pass_H;					/* Pass 2-D height                                                  */
} *Cinches[9],
	 Start  = {0, 0, 0, 0, 0},	/* Initial pass (stage) to read a raw input sequence */
	 Clean  = {0, 0, 0, 0, 0},	/* Pass to format original input string into acceptable characters and determine sequence type */
	Cinch_T = {0, 0, 0, 0, 0},	/* Pass to cinch Tandem repeats using a special traversal of the Pathbox grid */
	Cinch_L = {0, 0, 0, 0, 0},	/* Pass to cinch Long monomeric tracts that are of size >= 2* mrwrap */
	Cinch_K = {0, 0, 0, 0, 0},	/* Pass to cinch intra-repeat K-mers */
	Cinch_D = {0, 0, 0, 0, 0},	/* Pass to cinch De novo repeat structures (fractal repeats) based on the consensus row */
	 Relax  = {0, 0, 0, 0, 0},	/* Pass to relax monomeric tracts that did not aid cinch-D cinches */
	Recover = {0, 0, 0, 0, 0},	/* Optional pass to recover and check a 1-D sequence from a 2-D alignment */
	Current = {0, 0,-1, 0, 0};	/* Current holds the values from the latest pass, so that generic functions can simply look here */
								/* Current.pass_V is pass counter, initialized to -1 so that it is incremented to 0 for Start pass */
int cyc_count=0;

/* The segment struct type organizes the spine and belly row positions of the 2-D 'snake' for computing MSA tucks */
struct segment {
	int spine;		/* row positons of top-most symbol, including terminators, in each 2D column */
	int belly;		/* row positions of bottom-most buffer row, above which is the bottom most character in that column or adjacent column */
	char topc;		/* character of spine position */
};

struct cindstruct {
	int		mfirst;             /* First row with character */
	int		mlast;	            /* Last row with character */
	int		chcount;			/* Character count per column */
	int		stat;				/* Status information for column */
} *dConsensus = NULL;

char align2D[MAXROW+1][MAXROW] = {{0}};
char * cinch2D = NULL;
char * consensus = NULL;
char file_name[255] = "internal_default";
char dev_notes[32] = " ";      		    /* STRING WRITTEN AS LAST FIELD IN OUTPUT FILE */
short unsigned int cinchled=0;			/* BIT FLAG FOR CINCH-L WRAPS */
char letr_unit[4] = {0};				/* UNIT STRING: "bp" FOR DNA, "nt" FOR RNA, 'aa' FOR PROTEINS, 'ch' FOR ALL OTHER; SET IN MAIN() */
FILE *fp_out;                           /* FILE FOR OUTPUT.LOG */

void 				clear_2D_ar(char wipe_align2D[][MAXROW]);
void 				clear_cinch2D(void);
void 				clear_right(char swipe_align2D[][MAXROW]);
int 				col_isclear(  char  check_array[][MAXROW], unsigned int at_n, int row, short int updown);
int 				col_isclear1D(char *check_array          , unsigned int at_n, int row, short int updown);
unsigned int 		consensus_2D(int n_start, int n_width, short unsigned int print);
int 				count_wrap_blocks(int lcl_width, int lcl_opt_w);	/* lcl_width IS WIDTH OF 2-D MHA ARRAY */
int 				get_1Dz(int x, int y, int ignoreCheck);
void 				line_end(int type, int c, int lcl_width);
char 				mha_base62(int num);
void 				mha_head(int lcl_width);
void 				mha_writeback(char lcl_align2D[][MAXROW], char align2D[][MAXROW]);
void				mha_writeback_1Dto2D(char *cinch2D, char align2D[][MAXROW]);
void 				mha_writeback_2Dto1D(char lcl_align2D[][MAXROW], char *cinch2D);
void 				mha_writeconsensus1D(char *align2D_one, char consensus1D[MAXROW]);
int 				mn1D(int row, int col);
int 				mn1Dbig(int row, int col);
void 				print1D(void);
short unsigned int	print_2Dseq(void);
void 				print_blockhead(int numbl, int totbl);
void				print_protein_waxes(void);
short int 			pushdown(char pusharray[][MAXROW], int push_m, int push_n);
int 				span_ork(int point);
int 				span_rk(int point);
void 				warnhead(char l);
short unsigned int 	recoverlen(void);
short unsigned int	cleanseq(char *s);
void 				mha_randomize1(char *input_seq);
void 				mha_randomize2(char *input_seq, int rsize);
void 				print_base62_table(void);
void 				shuffle(int *array, int n);
void 				usage(char *usage_version);	/* FOR PRINTING UNIFORM USAGE INFO */
void 				free_2D(int **p2D, int lenseq);
struct segment 		makesegment(int top, int bottom);
struct segment * 	makesnake(char *array, int height, int width, int w_plus_rattle, short unsigned int zcol);
short unsigned int 	checkfractals_in_imperfect(int kmer, int n);
void 				print_section_spacer(void);

#include "microhom-devl.h"

/*********************************************************************************************************/
int count_unique_chars(char *string, int window)
{
	int i=0, j=0, l=0;
	char letr;
	short unsigned int seqtype = Clean.pass_V;
	short unsigned int max = 24;

	char unique[24] = {0};

	if (seqtype==1 || seqtype==2)
		max=4;
	else if (seqtype==3)
		max=20;

	for (i=0; i<window; i++) {
		if (j==max || (letr=string[i])==slip.sym || letr==Term->sym || letr=='\0')
			return(j);
		else {
			for (l=0; l<j; l++)	{
				if (letr==unique[l])
					break;
			}
			if (l==j && l<32)
				unique[j++] = letr;
		}
	}
	return(j);	/* SHOULD NEVER GET TO HERE UNDER NORMAL CIRCUMSTANCES */
}

/*********************************************************************************************************/
struct segment * makesnake(char *array, int height, int width, int w_plus_rattle, short unsigned int zcol)
{
	struct segment *snake = NULL;

	snake = (struct segment *)calloc(w_plus_rattle+2, sizeof(struct segment));

	int i, j;
	int maxtop=0;
	int maxbot=0;
	char letr='\0';

	/* POPULATE FIRST COLUMN OF SNAKE STRUCT */
	if (zcol) {
		snake[0].spine = 0;
		snake[0].topc = Term->sym;
		for (i=height; i>0; i--) {
			if (isalpha(*(array + MAXROW*(i-1)))) {
				if (i > maxbot)
					maxbot = i;
				break;
			}
		}
		snake[0].belly = maxbot;
	}

	for (j=0; j<=width; j++) {		/* zcol should be 0 or 1 and refers to which column the 2-D alignment begins */
		for (i=0; i<height; i++) {
			if (isalpha( (letr=*(array + i*MAXROW + j)) ) || letr==slip.sym || letr==Term->sym || letr==monoR.sym) {
				if (i > maxtop)
					maxtop = i;
				break;
			}
		}
		for (i=height; i>0; i--) {
			if (isalpha(*(array + MAXROW*(i-1) + j + 1))) {
				if (i > maxbot)
					maxbot = i;
				break;
			}
		}
		snake[j+zcol].spine = maxtop;
		snake[j+zcol].belly = maxbot;
		snake[j+zcol].topc = letr;
	}
	for (j=width+1; j<=w_plus_rattle; j++) {
		snake[j+zcol].spine = maxtop;
		snake[j+zcol].belly = maxbot;
		snake[j+zcol].topc = Term->sym;
	}

	return(snake);
}



/*****************************************************************************************/
static int random_i(int n)
{
int limit = RAND_MAX - RAND_MAX % n;
int rnd;

	do {
		rnd = rand();
    }
    while (rnd >= limit);

    return rnd % n;
}

/*****************************************************************************************/
int get_1Dz(int x, int y, int ignoreCheck)
{
	int i=0, z=0, count_check=0;
	int lenseq = Clean.pass_W;

	for (i=0; i<lenseq; i++) {
		if (tela[i].x == x && tela[i].y == y) {
			count_check++;
			z = i;
			if (ignoreCheck)
				break;
		}
	}
	if (count_check==1)
		return(z);
	else {
		return(-z);		/* KEEP NEGATIVE TO PRESERVE INFO/ AVOID RETURNED VALUE LETTING A FOR LOOP GO */
	}
}
/*****************************************************************************************/

/*****************************************************************************************/
void clear_2D_ar(char wipe_align2D[][MAXROW])
{
int m=0, n=0;
int lenseq = Clean.pass_W;
int twidth = Cinch_T.pass_W;
int height = Current.pass_H;
char *aptr_start = &wipe_align2D[0][0];
char *aptr = aptr_start;

	for (n=0; n <= lenseq; n++, aptr++)
		*aptr = '\0';
	for (m=1; m <= height; m++) {
		aptr = aptr_start + m*MAXROW;
		for (n=0; n <= twidth; n++, aptr++)
			*aptr = '\0';
	}
}
/*****************************************************************************************/


/** INVOKE PRIOR TO RE-USING CINAR2D - CINCH ARRAY 2D **/
void clear_cinch2D(void)
{
int i=0;
	for (i=0; i<(Cinch_T.pass_W+1)*Clean.pass_W; i++)
		cinch2D[i] = '\0';
}


/*****************************************************************************************/
void clear_right(char swipe_align2D[][MAXROW])
{
int m=0, n=0;
int twidth = Cinch_T.pass_W;
int height = Current.pass_H;
char letr;
char *aptr_start = &swipe_align2D[0][0];
char *aptr = aptr_start;

	/* CLEAR TO THE RIGHT OF ROW TERMINATORS */
	for (m=1; m < height; m++) {
		n = 0;
		while ( (letr=swipe_align2D[m][n]) == Fill->sym || isalpha(letr) || letr==monoL.sym) {
			n++;
		}
		aptr = aptr_start + m*MAXROW + n + 1;
		for (n = n+1; n <= twidth; n++, aptr++)
			*aptr = '\0';
		if (letr == Term->sym) {
			aptr = aptr_start + (m+1)*MAXROW;
			for (m=m+1; m < height; m++) {
				for (n=0; n <=twidth; n++, aptr++) {
					*aptr = '\0';
				}
			}
			break;
		}
	}
}
/*****************************************************************************************/


/*****************************************************************************************/
unsigned int consensus_2D(int n_start, int n_width, short unsigned int print)
{
	int badsites=0, m=0, n=0, n_end, x=1;
	int con_width = Current.pass_W;
	short unsigned int nuctype = Clean.pass_V;			/* FOR SEQ TYPE, DNA=1, RNA=2, OTHER (NOT NA)=0 */
	short unsigned int nuctransit = 0;					/* BIT FLAG FOR HANDLING NUCLEAR TRANSITIONS */
	short unsigned int plustransit=0;					/* BIT FLAG ADDENDUM FOR COUNTING BADSITES AT COL */
	short unsigned int checktransit=0;					/* BIT FLAG FOR CHECKING GOOD TRANSITION MARK */
	char blnk = Fill->sym;
	char letr, ltr2, conletr=blnk;
	char checktransletr=blnk;
	int conrows = 26; 						/* COL   n=0 FOR BIT FLAG */
	int consensus_ar[26][MAXROW] = {{0}};	/* ROW   m=0 FOR COUNTER */
											/* ROW   m=1 FOR CONSENSUS */
											/* ROWS 2-24 FOR VARIANTS STORAGE */
											/* ROW  m=25 FOR CHECKTRANSITS UNICODE */
	n_end = n_start + n_width;

	if (nuctype == 1)	/* IF DNA */
		nuctransit = 1;

	/* FILL CONSENSUS ARRAY WITH FIRST LETTER IN EACH COL */
	for (n=0; n<=con_width; n++) {
		for (m=0; m<Clean.pass_W; m++) {
			if (isalpha(letr=align2D[m][n])) {
				consensus_ar[0][n+1]++;
				if (letr <= 'Z' || (nuctype && letr==ambig.sym))
					consensus_ar[1][n+1] = letr;
				else								/* ELSE MAKE UPPERCASE */
					consensus_ar[1][n+1] = toupper(letr);
				break; /* OUT OF FOR m LOOP */
			}
		}
	}
	consensus_ar[0][0] = 1;			/* COUNTER ROW ON */
	consensus_ar[0][n] = '\0';
	consensus_ar[1][0] = 1;			/* CONSENSUS ROW ON */
	consensus_ar[1][n] = '\0';

	for (n = n_start; n <= n_end; n++) {
		if (nuctransit) {
			checktransit = plustransit = 0;		/* RE-INITIALIZE AT EACH COLUMN */
			if ((conletr=consensus[n])=='R' || conletr=='Y') {
				plustransit = 1;
				if (conletr == 'R')
					checktransit = 18;
				else if (conletr == 'Y')
					checktransit = 25;
				consensus_ar[1][n+1] = conletr;
			}
			else
				conletr = blnk;
		}

		for (m = 0; align2D[m][0] != '\0'; m++) {
			if (isalpha(letr=align2D[m][n])) {
				ltr2=toupper(letr);
				if (nuctransit && ((conletr=='R' && (ltr2=='G' || ltr2=='A')) || (conletr=='Y' && (ltr2=='C' || ltr2=='T'))) ) {
					if (conletr=='R') {
						if (checktransit==18) {
							if (ltr2=='A')
								checktransit -=  1;
							else if (ltr2=='G')
								checktransit -= 17;
						}
						else if (checktransit==17 && ltr2=='G')
							checktransit = 0;
						else if (checktransit== 1 && ltr2=='A')
							checktransit = 0;
					}
					if (conletr=='Y') {
						if (checktransit==25) {
							if (ltr2=='C')
								checktransit -=  3;
							else if (ltr2=='T')
								checktransit -= 22;
						}
						else if (checktransit==22 && ltr2=='T')
							checktransit = 0;
						else if (checktransit== 3 && ltr2=='C')
							checktransit = 0;
					}
				}
				else if (letr!=ambig.sym && letr!=consensus_ar[1][n+1]) {
					if (consensus_ar[0][n+1] == 1) {	/* IF THIS IS THE FIRST CONFLICT NOTED AT THIS POSITION */
						++badsites;						/* COUNT ADDITIONAL BAD SITE COLUMN */
						++consensus_ar[0][n+1];			/* INCREMENT COUNTER FOR DIFFERENT LETTERS AT COLUMN */
						consensus_ar[2][ 0 ] = 1;		/* TURN VARIANT ROW ON */
						consensus_ar[2][n+1] = letr;
					}
					else {		/* ELSE CHECK IF THIS IS A NEW LETTER AT THIS COLUMN */
						for (x = 2; x < conrows-1; x++) {
							if (letr == consensus_ar[x][n+1])
								break;
							else if (consensus_ar[x][n+1] == '\0') {
								consensus_ar[x][ 0 ] = 1;		/* TURN VARIANT ROW ON */
								consensus_ar[x][n+1] = letr;	/* ASSIGN LETTER */
								++consensus_ar[0][n+1];			/* INCREMENT COUNTER FOR # UNIQUE LETTERS AT COLUMN */
								break;
							}
						}
					}
				}
				else if (x >= conrows-1)
					break;
			}
			else if (letr == Term->sym) {
				while (col_isclear(align2D,n,m,-1)>-1)
					n++;
				consensus_ar[1][n+1] = letr;
				Current.pass_H = m+1;		/* STORE HEIGHT IN HEIGHT SLOT */
			}
		} /* END OF m LOOP THROUGH COLUMN n */

		if (checktransit) {		/* IF checktransit IS STILL POSITIVE THEN IT'S SUPPOSED TO NOT HAVE CHECKED OUT */
			checktransletr = consensus_ar[1][n+1];
			plustransit = 0;
			if (Current.pass_V) { /* IF PASS NUMBER */
				consensus_ar[25][0]=1;
				consensus_ar[25][n+1]=checktransletr;
				sprintf(dev_notes, "checktransit=%d at n=%d", checktransit, n);
			}
		}
		else if (plustransit)
			++consensus_ar[0][n+1];
	} /* END OF FOR n LOOP */

	if (print && opt_K.bit) {							/* opt_K SHOW CONSENSUS ROW */
		/* PRINT CONSENSUS ROWS */
		for (m = 1; consensus_ar[m][0] != '\0' && m<conrows-1; m++) {
			line_end(BLOCKHEAD, 9, 9);
			printf(" ");
			for (n = n_start; n < n_end; n++) {
				if (isalpha(letr=consensus_ar[m][n+1]))
					printf("%c", letr);
				else
					printf("%c",blnk);
			}
			printf("\n");
		}
		if (consensus_ar[25][0]) {		/* IF CHECKTRANSITS INCOMPLETE */
			line_end(BLOCKHEAD, 9, 9);
			printf(" ");
			for (n = n_start; n < n_end; n++) {
				if (consensus_ar[25][n]) {
					printf("%c", consensus_ar[25][n]);
					consensus_ar[1][n] = blnk;
				}
				else
					printf(".");
			}
			printf(" <= check transition\n");
		}
		if (badsites) {
			line_end(BLOCKHEAD, 9, 9);
			printf(" ");
			for (n = n_start; n < n_end; n++) {
				letr = mha_base62(consensus_ar[0][n+1]);
				printf("%c", letr);
			}
		    printf("\n");
			warnhead('C');
			if (badsites == 1)
				printf("There is one column with mimatched letters.");
			else if (badsites > 1)
				printf("There are %d columns with mismatched letters.", badsites);
		}
	}

	/* STORE CONSENSUS ROW */
	for (n = n_start; n <= n_end; n++) {
		consensus[n] = consensus_ar[1][n+1];
	}
	consensus[con_width] = '\0';
	if (!print) {
		if (!badsites) 
			Current.pass_Q = 1000;
		else
			Current.pass_Q = round((1000*(con_width-badsites))/con_width);
	}
	return(badsites);
}
/*****************************************************************************************/

/*****************************************************************************************/
int count_wrap_blocks(int lcl_width, int lcl_opt_w)
{
int lcl_blocks2D;
	if ( (lcl_width > lcl_opt_w) && (lcl_width % lcl_opt_w == 0) )
		lcl_blocks2D = lcl_width/lcl_opt_w;
	else if ( (lcl_width > lcl_opt_w) && (lcl_width % lcl_opt_w != 0) )
		lcl_blocks2D = lcl_width/lcl_opt_w + 1;
	else
		lcl_blocks2D = 1;

	return(lcl_blocks2D);
}
/*****************************************************************************************/

/*****************************************************************************************/
unsigned int foam_2D(int n_start, int n_width)
{
int m=0, n=0, n_end;
char letr;
int foam_ar[3][MAXROW] = {{0}};	 	/* ROW m=0 FOR COLUMN LETR COUNTER 	*/
									/* ROW m=1 FOR CONSENSUS 			*/
									/* ROW m=2 FIRST LETTER ROW			*/
	n_end = n_start + n_width;

	for (n = n_start; n <= n_end; n++) {
		for (m = 0; m < MAXROW; m++) {
			if (isalpha(letr=align2D[m][n])) {
				foam_ar[0][n] = 1;
				foam_ar[1][n] = letr;
				foam_ar[2][n] = m;
				break; /* OUT OF FOR m LOOP */
			}
		}
	}
	foam_ar[0][n] = '\0';
	foam_ar[1][n] = '\0';
	foam_ar[2][n] = '\0';

	for (n = n_start; n <= n_end; n++) {
		for (m = foam_ar[2][n]+1; align2D[m][0] != '\0'; m++) {
			if ( isalpha(align2D[m][n]) ) {
				++foam_ar[0][n];		/* INCREMENT COUNTER */
				foam_ar[1][n] = fill_1.sym;
			}
		}
	}

	if (opt_f.bit) {	/* opt_f FOAM-FREE CONSENSUS ROW */
		line_end(BLOCKHEAD, 9, 9);
		printf(" ");
		for (n = n_start; n < n_end; n++) {
			printf("%c", foam_ar[1][n]);
		}
		printf(" = foam-free\n");
	}

	return(0);
}
/*****************************************************************************************/

/*****************************************************************************************/
int next_foamfree(char check_array[][MAXROW], int row, int at_n)
{
	int n = at_n;

	while ( isalpha(check_array[row][n]) ) {
		if (col_isclear(check_array,n,row,-1) == -1 && col_isclear(check_array,n,row,1) == -1)
			return(n);
		else
			n++;
	}

	return(-1);

}
/*****************************************************************************************/

/*****************************************************************************************/
int col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown)
{
int m=0;
	/* INVOKE updown AS -1 OR 1 TO CHECK BELOW OR ABOVE row, RESPECTIVELY */
	for (m = row + updown; m >= 0; m += updown) {
		if (isalpha(check_array[m][at_n]))
			return(m);
		else if (check_array[m][0]=='\0')
			break;
	}
	return(-1);		/* RETURN NEGATIVE VALUE TO INDICATE COLUMN IS CLEAR IN DESIGNATED DIRECTION */
}
/*****************************************************************************************/

/*****************************************************************************************/
int col_isclear1D(char *check_array, unsigned int at_n, int row, short int updown)
{
int m=0;
	/* INVOKE updown AS -1 OR 1 TO CHECK BELOW OR ABOVE row, RESPECTIVELY */
	for (m = row + updown; m >= 0; m += updown) {
		if (isalpha(check_array[mn1D(m,at_n)]))
			return(m);
		else if (check_array[mn1D(m,0)]=='\0')
			break;
	}
	return(-1);		/* RETURN NEGATIVE VALUE TO INDICATE COLUMN IS CLEAR IN DESIGNATED DIRECTION */
}
/*****************************************************************************************/

/*****************************************************************************************/
void line_end(int type, int c, int lcl_width)
{
	char rule1[] = "_________|_________|_________|_________|_________|_________|_________|_________|"
				   "_________|_________|_________|_________|_________|_________|_________|_________|"
				   "_________|_________|_________|_________|_________|_________|_________|_________|"
				   "_________|_________|_________|_________|_________|_________|_________|_________|";
	char rule2[] = "---------|---------|---------|---------|---------|---------|---------|---------|"
				   "---------|---------|---------|---------|---------|---------|---------|---------|"
				   "---------|---------|---------|---------|---------|---------|---------|---------|"
				   "---------|---------|---------|---------|---------|---------|---------|---------|";
	/* START		type=0 */
	/* END			type=1 */
	/* RULER		type=2 */
	/* SLIPS		type=3 */
	/* PATHBOXHEAD	type=4 */
	/* BLOCKHEAD	type=5 */
	/* SLIPRULER	type=6 */

	char *ruler = rule1;	/* USE TO CHANGE RULE STYLE */

	if (type == 1)					/* FORMAT FOR LINE END. c IS CHAR. NUMBER */
		printf("  %d\n", c);

	if (opt_r.bit) {				/* FORMAT WITH LINE NUMBERING. DEFAULT */
		if (type == 0) {			/* FORMAT FOR LINE BEGINNING. c is LINE NUMBER (m+1) */
			if (c == 1)
				printf(" %4d. %c", c, Term->sym);
			else
				printf(" %4d. %c", c, ZTick->sym);
		}
		else if (type == 2) {			/* FORMAT FOR RULER. */
			printf("       %c%.*s\n", ZTick->sym, lcl_width, ruler + c);		/* DOUBLE USE OF c AS O-F-FSET FOR RULER */
		}
		else if (type == 3)			/* FORMAT FOR NUMBERED SLIP LOCATION LINES. */
			printf("        ");
		else if (type == 4)			/* FORMAT FOR PATHBOX HEADERS, #'d */
			printf("      _|");
		else if (type == 5)			/* FORMAT FOR CONSENSUS BLOCKHEAD HEADERS, #'d */
			printf("       ");
		else if (type == 6) {		/* FORMAT FOR SLIP RULER. */
			ruler = rule1;
			printf("        %.*s\n", lcl_width, ruler + c);		/* DOUBLE USE OF c AS O-F-FSET FOR RULER */
		}
	}
	else {		 				 	/* FORMAT WITHOUT LINE NUMBERS. r = REMOVED */
		if (type == 0) {		   	/* FORMAT FOR LINE BEGINNING. c is LINE NUMBER (m+1) */
			if (c == 1)
				printf("   %c", Term->sym);
			else
				printf("   %c", ZTick->sym);
		}
		else if (type == 2) {			/* FORMAT FOR RULER OR SLIPRULER. */
			printf("   %c%.*s\n", ZTick->sym, lcl_width, ruler + c);		/* DOUBLE USE OF c AS O-F-F-SET FOR RULER */
		}
		else if (type == 3)			/* FORMAT FOR UNNUMBERED SLIP LOCATION LINES. */
			printf("    ");
		else if (type == 4)			/* FORMAT FOR PATHBOX HEADERS, STILL #'d */
			printf("      _|");
		else if (type == 5)			/* FORMAT FOR CONSENSUS BLOCKHEAD HEADERS, #'d */
			printf("   ");
		else if (type == 6) {		/* FORMAT FOR SLIPRULER. */
			ruler = rule2;
			printf("   %c%.*s\n", ZTick->sym, lcl_width, ruler + c);		/* DOUBLE USE OF c AS O-F-F-SET FOR RULER */
		}
	}
}
/*****************************************************************************************/

/*****************************************************************************************/
char mha_base62(int num)	/* TO REPRESENT NUMBERS UP TO 62 WITH A SINGLE CHAR. ******/
{
char c;
	if (num < 10)
		c = '0' + num;
	else if (num >=10 && num < 36)
		c = 'A'+ num-10;
	else if (num >=36 && num < 62)
		c = 'A' + num-4;
	else
		c = '!';
	return(c);
}
/*****************************************************************************************/

/*****************************************************************************************/
void mha_head(int lcl_width)
{
/*
char h1[]=	"\n_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//"
			  "_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//"
			  "_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//";
*/
char h1[]=	"\n__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),"
			  "__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),"
			  "__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),__mha__/O)),";
char h2[]=	"\n\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			  "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			  "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_";
char *h_rule = h1;						/* DEFAULT BANNER STYLE */
/*int min_len = 80;*/					/* MINIMUM LENGTH */
int med_len = 12*(lcl_width/10);		/* MEDIUM LENGTH, SCALING */
int max_len = par_wrap.set+8;			/* MAXIMUM LENGTH */
int hr_len = max_len;					/* DEFAULT LENGTH OF HEADER BANNER */

	if (lcl_width > MAXROW && dev_print(LOGY, __LINE__)) {
		printf("Bad news bears: Unexpectedly, lcl_width > MAXROW. Current.pass_V=%d, lcl_width=%d.\n", Current.pass_V, lcl_width);
	}

	if (lcl_width+8 > hr_len && lcl_width+8 < par_wrap.set)
		hr_len = med_len;
	else if (lcl_width+8 >= par_wrap.set)
		hr_len = max_len;

	if (opt_X.bit)
		h_rule = h2;

	if (Current.pass_V == 5) {
		printf("%.*s\n", hr_len, h_rule);
		printf(" 2-D pass #%d: cinch-d run %d (width = %d)\n\n", Current.pass_V, Cinch_D.pass_R, lcl_width);
	}
	else if (Current.pass_V == 4) {
		printf("%.*s\n", hr_len, h_rule);
		printf(" 2-D pass #%d: cinch-k (width = %d)\n\n", Current.pass_V, lcl_width);
	}
	else if (Current.pass_V == 3) {
		printf("%.*s\n", hr_len, h_rule);
		printf(" 2-D pass #%d: cinch-l (width = %d)\n\n", Current.pass_V, lcl_width);
	}
	else if (Current.pass_V == 2) {
		printf("%.*s\n", hr_len, h_rule);
		printf(" 2-D pass #%d: cinch-t (width = %d)\n\n", Current.pass_V, lcl_width);
	}
	else if (Current.pass_V == 1) {
			printf("\n\n 1-D sequence:\n");
	}
	else if (Current.pass_V == 6) {
		printf("%.*s\n", hr_len, h_rule);
		printf(" 2-D pass #%d: relax-2D (width = %d)\n\n", Current.pass_V, lcl_width);
	}
	else if (Current.pass_V == 0) {
		printf("\nOriginal string (length = %d):\n", lcl_width);
	}
	else
		printf("%.*s\n", max_len, h_rule);
}
/*****************************************************************************************/

/*****************************************************************************************/
void mha_randomize1(char *input_seq)
{
    int i, rand_num;
	int seqlen = Clean.pass_W;
	char *random_seq = NULL;

	random_seq = (char *)calloc(seqlen, sizeof(char));

    for (i = 0 ; i < seqlen; i++)
    {
        rand_num = rand() % seqlen;
		random_seq[i] = input_seq[rand_num];
    }

	for (i = 0; i < seqlen; i++) {
		input_seq[i] = random_seq[i];
	}
	input_seq[i] = '\0';

	free(random_seq);
}
/*****************************************************************************************/

/*****************************************************************************************/
void mha_randomize2(char *input_seq, int rsize)
{
	int i, j, tmp;
	int seqlen = Clean.pass_W;
	char *random_seq = NULL;
	int *rnumbers= NULL;

	random_seq 	= (char *)calloc(rsize, sizeof(char));
	rnumbers 	= (int *)calloc(rsize, sizeof(int));

	for (i = 0; i < rsize; i++) {
		rnumbers[i]= i;
	}

	for (i = rsize - 1; i > 0; i--) {
    	j = random_i(i + 1);
        tmp = rnumbers[j];
        rnumbers[j] = rnumbers[i];
        rnumbers[i] = tmp;
	}

	for (i = 0; i < rsize; i++) {
		random_seq[i] = input_seq[(rnumbers[i] % seqlen)];
	}

	for (i = 0; i < rsize; i++) {
		input_seq[i] = random_seq[i];
	}
	input_seq[i] = '\0';	/* B/C SOMETIMES input_seq > rsize */

	free(rnumbers);
	free(random_seq);
}
/*****************************************************************************************/

/*****************************************************************************************/
void mha_writeback(char lcl_align2D[][MAXROW], char align2D[][MAXROW])
{
char letr;
int lenseq = Clean.pass_W;
int i=0, m=0, n=0, widest_n=0;

	/* WRITE BACK TO align2D */
	clear_right(lcl_align2D);
	clear_2D_ar(align2D);
	for (m = 0; lcl_align2D[m][0] != '\0' && m <= lenseq; m++) {
		for (n = 0; (letr=lcl_align2D[m][n]) != slip.sym && letr != Term->sym && letr != monoR.sym; n++) {
			align2D[m][n] = letr;
			if (letr == monoL.sym && lcl_align2D[m][n+opt_M.val+1] == monoR.sym)
				align2D[m][n+opt_M.val+2] = '\0';
		}
		align2D[m][n  ] = letr;	/* MHA-STANDARD TERMINATOR  */
		if (n > widest_n)
			widest_n = n;
		if (letr == Term->sym) {
			Current.pass_W = widest_n;	/* ASSIGN 2-D WIDTH  */
			Current.pass_H = m+1;		/* ASSIGN 2-D HEIGHT */
		}
		for (i=n+1; i<Clean.pass_W; i++)
			align2D[m][i] = '\0';
	}
}

/************ FOR USE WITH cin2Dar *******************************************************/
void mha_writeback_2Dto1D(char lcl_align2D[][MAXROW], char *cinch2D)
{
char letr;
int lenseq = Clean.pass_W;
int i=0, m=0, n=0, widest_n=0;

	/* WRITE BACK TO cinch2D */
	clear_right(lcl_align2D);
	clear_cinch2D();
	for (m = 0; lcl_align2D[m][0] != '\0' && m <= lenseq; m++) {
		for (n = 0; (letr=lcl_align2D[m][n]) != slip.sym && letr != Term->sym && letr != monoR.sym; n++) {
			cinch2D[mn1D(m,n)] = letr;
			if (letr == monoL.sym && lcl_align2D[m][n+opt_M.val+1] == monoR.sym)
				cinch2D[mn1D(m,n+opt_M.val+2)] = '\0';
		}
		cinch2D[mn1D(m,n)] = letr;	/* MHA-STANDARD TERMINATOR  */
		if (n > widest_n)
			widest_n = n;
		if (letr == Term->sym) {
			Current.pass_W = widest_n;	/* ASSIGN 2-D WIDTH  */
			Current.pass_H = m+1;		/* ASSIGN 2-D HEIGHT */
		}
		for (i=n+1; i<Clean.pass_W; i++)
			cinch2D[mn1D(m,i)] = '\0';
	}
}

/************ FOR USE WITH cin2Dar *******************************************************/
void mha_writeback_1Dto2D(char *cinch2D, char align2D[][MAXROW])
{
char letr;
int i=0, m=0, n=0, widest_n=0;
int bound = (Cinch_T.pass_W+1)*Clean.pass_W;

	/* WRITE BACK TO align2D */
	clear_2D_ar(align2D);
	for (i=0; (letr=cinch2D[i])!=Term->sym && i<=bound; i++) {
		if (letr==Fill->sym || isalpha(letr) || letr==monoL.sym)
			align2D[m  ][n++] = letr;
		else if (letr==slip.sym || letr==monoR.sym)	{
			if (n>widest_n)
				widest_n = n;
			align2D[m][n  ] = letr;
			align2D[m][n+1] = '\0';
			m++;
			n = 0;
			i = (Cinch_T.pass_W+1)*m - 1;
		}
		else if (letr == '\0') {
			warnhead('?');
			printf("Missing a line terminator of some sort? prev letr=%c, m=%2d, n=%2d (microhomology.h line %d)\n", cinch2D[i-1], m,n,__LINE__);
		}
	}
	if (letr == Term->sym) {
		if (n>widest_n)
			widest_n = n;
		align2D[m][n  ] = letr;
		align2D[m][n+1] = '\0';
		for (n=0; n<=Cinch_T.pass_W; n++)
			align2D[m+1][n] = '\0';
		Current.pass_W = widest_n;	/* ASSIGN 2-D WIDTH  */
		Current.pass_H = m+1;		/* ASSIGN 2-D HEIGHT */
	}
}

/*****************************************************************************************/
void mha_writeconsensus1D(char *cinch2D, char consensus1D[MAXROW])
{
	int i = 0;
	int constart = Clean.pass_W * (Cinch_T.pass_W+1);

	/* COPY CONSENSUS ROW FROM ARRAY ONE TO 1D CONSENSUS ARRAY */
	for (i=0; i<=Cinch_T.pass_W+1; i++) {
		consensus1D[i] = cinch2D[constart+i];
	}
	consensus1D[i] = '\0';
}


/** RETURN 1D POSITION FOR USE IN CINCH-T 2D ARRAY POINTER **/
int mn1Dbig(int row, int col)
{
	return(row*MAXROW + col);
}

/** RETURN 1D POSITION FOR USE IN POST-CINCH-T CINCH 2D ARRAY POINTER **/
int mn1D(int row, int col)
{
	return(row*(Cinch_T.pass_W+1) + col);
}


/*****************************************************************************************/
short unsigned int print_2Dseq(void)
{
unsigned int foam_2D(int n_start, int n_width);
int all_clear;		/* COUNTER VARIABLE USED FOR CHECKING NEED TO PRINT BOTTOM ROWS */
int blocks2D=0, b=0, c=0, carry_over=0, d=0, fudge=0, g, h, i, j, j_start=0, j_end, m, m_start=0, n;
int mmsites=0, max_n=0;
int linewrap = par_wrap.set;
char letr, next;
char blnk  = Fill->sym;			/* opt_B blank character */
int cinchwidth = Current.pass_W;
int   lenseq = Clean.pass_W;
int head_start;							/* USE TO PASS RULER O-F-F-SET TO line_end() */
int scrimmageline;						/* USE TO INCREMENT AND TEST IF FILLER IS NEEDED, CAN BE OPTION TO DO SO */
char tick = fill10.sym;
short unsigned int lcl_opt_F = 0;

	if (opt_B.val > 1)
		tick = blnk;

	mha_head(cinchwidth);

	blocks2D = count_wrap_blocks(cinchwidth, linewrap);

	for (j = 0; j < blocks2D; j++) {
		if (blocks2D != 1)
			print_blockhead(j+1, blocks2D);

		/* FAST FORWARD TO INFORMATIVE ROW */
		m_start = 0;
		while (align2D[m_start][(j * linewrap)] == '\0') {
			m_start++;
		}

		scrimmageline = 1;
				max_n = 0;

		for (m = m_start; align2D[m][0] != '\0'; ) {

			b = d = 0;		/* VAR b WILL COUNT BLANKS, VAR d WILL COUNT ALPHA CHAR ANEW FOR EACH ROW */

			j_start =   j   * linewrap;
			j_end   = (j+1) * linewrap;

			for (n = j_start; n < j_end && (letr=align2D[m][n])!=slip.sym && letr!=Term->sym && letr!=monoR.sym; n++) {
				if (isalpha(letr)) {
					c++;
					d++;
					if (d+b > scrimmageline)
						scrimmageline++;
				}
				else if (letr == blnk || letr == monoL.sym) {
					b++;
				}
			} /* END OF n SCAN LOOPS */
			if (b+d > max_n)
				max_n = b+d;

			if (b == linewrap) {
				all_clear = 1;
				for (g = m+1; align2D[g][0] != '\0' && g < MAXROW; g++) {
					for (h = 0; h < j_end; h++) {
						if (isalpha(align2D[g][h]))
							all_clear = 0;
					}
				}
				if (all_clear)
					break;
			}

			line_end(START, m+1, 0);

			for (n=j_start; n<j_end && (letr=align2D[m][n])!='\0' && letr!=Term->sym && letr!=monoR.sym; n++) {
				if (d+b == 0) {		/* HANDLES BLOCK CONTINUATION LINES THAT HAVE ZERO CHARACTERS B/C OF ITS TUCK & LENGTH */
					break;
				}
				else if ((letr=align2D[m][n]) == blnk && (n+1) % 10 == 0 && opt_B.val < 2)
					printf("%c", tick);		/* PRINT TICK MARKS AT 10 bp INTERVALS IF NOT BLANK SPACE */
				else if (Fill->sym == 32)		/* opt_B BLANK = SPACE ' ' */
					printf("%c", letr);
				else					/* opt_B BLANK = FULLSTOP '.' */
					printf("%c", letr);

				/* TURN ON OPTION opt_F TO USE BLANK CHAR TO FILL < SCRIMMAGE IF EXCESSIVELY SHORT */
				if (!opt_F.bit && align2D[m][n] == slip.sym && scrimmageline-(b+d-2) > 10)
					lcl_opt_F = 1; 	/* TURN ON local opt_F IN THIS CASE */

				/* OPTION opt_F TO USE BLANK CHAR TO FILL < SCRIMMAGE */
				if ((opt_F.bit || lcl_opt_F) && align2D[m][n] == slip.sym && b+d-2 < scrimmageline) {
					for (i = 0; b+d+i < scrimmageline; i++) {
						if ( (b+d+2+i)%10 == 0)		/* WHY 2? +1 FOR STARTING AT 0, +1 FOR '/' CHAR */
							printf("%c", tick);		/* PRINT TICK MARKS AT 10 bp INTERVALS IF NOT BLANK SPACE */
						else
							printf("%c", blnk);		/* PRINT LINE-END FILLER CHARACTER */
					}
					lcl_opt_F = 0;
				}
			} /* END OF n PRINTING LOOPS */

			next = align2D[m][n];		/* TERMINAL CHARACTERS WILL ALSO PRINT AT END OF BLOCK IF...	*/
										/* ...THEY ARE FIRST CHARACTER OF NEXT BLOCK					*/

			if (next == slip.sym || next == monoR.sym) 		/* TO INDICATE ADJACENCY TO SLIP, WHICH ALSO WILL SHOW IN NEXT BLOCK */
				printf("%c", next);
			else if (next==blnk || (isalpha(next) && d+b != 0))	/* TRUE IF AT EDGE OF opt_w WINDOW */
				printf("=>");			/* TO INDICATE CONTINUATION TO NEXT BLOCK */
			else if (next == '\0' && n == j_start)
				printf("%c", blnk);
			else if (next == Term->sym) {		/* LAST CHARACTER */
				printf("%c", next);
				if (opt_L.bit)	/* opt_L = LINE END NUMBERING (ON) */
					line_end(END, c, 0);
				else
					printf("\n");
				break;	/* BREAK OUT OF m=m_start FOR LOOP */
			}

			if (opt_L.bit)	/* opt_L = LINE END NUMBERING ON */
				line_end(END, c, 0);
			else
				printf("\n");
			m++;

		} /* END OF m=m_start FOR LOOP */

		/* PRINT RULER */
		head_start = (j * linewrap) % 10;
		if (j+1 < blocks2D) {
			line_end(RULER, head_start, linewrap);
		}
		else {
			line_end(RULER, head_start, (linewrap=max_n));
		}

		/* PRINT NUMBERS FOR CONSENSUS RULER */
		if (opt_B.val < 4) {
			line_end(SLIPS, head_start, linewrap);
			fudge = carry_over % 10;
			printf("%.*s", 9-fudge, "          ");
			for (n = 1; n <= linewrap+fudge; n++) {
				if ( (n+carry_over-fudge) % 10 == 0)
					printf("%-10d", n + carry_over-fudge);
			}
			carry_over = carry_over + n - fudge - 1;
			printf("\n");
		}
		/* *********************** */

		/* ADD TO COUNT OF MISMATCHED SITES */
		mmsites = mmsites + consensus_2D(j_start, linewrap, ON);

		if (opt_f.bit && Current.pass_V > 5 && !mmsites) {
			foam_2D(j_start, linewrap);
		}
	} /* END OF FOR j PRINTING LOOP */

	if (c == lenseq) {
		if (!mmsites) {
			Current.pass_Q = 1000;
			print_section_spacer();
			return(0);
		}
		else {
			Current.pass_Q = round((1000*(cinchwidth-mmsites))/cinchwidth);
			print_section_spacer();
			return(Current.pass_Q);
		}
	}
	else if (c < lenseq) {
		Current.pass_Q = round((1000*(cinchwidth-mmsites-(lenseq-c)))/cinchwidth);
		warnhead(gap.sym);
		printf("2-D auto-alignment is missing %d %s(s)!\n\n", lenseq-c, letr_unit);
		print_section_spacer();
		return(0);
	}
	else {
		Current.pass_Q = round((1000*(cinchwidth-mmsites-(c-lenseq)))/cinchwidth);
		warnhead('+');
		printf("2-D auto-alignment contains an extra %d %s(s)!", c-lenseq, letr_unit);
		print_section_spacer();
		return(0);
	}
}

/*****************************************************************************************/
void print_base62_table(void)
{
unsigned int i=0;

	mha_head(80);
	printf("\nMHA base 62 single-letter representation:\n\n Base 10: ");
	for (i = 0; i < 31; i++)
		printf("%2d ", i);
	printf("\n Base 62: ");
	for (i = 0; i < 31; i++)
		printf("%2c ", mha_base62(i));
	printf("\n\n Base 10: ");
	for (i = 31; i < 62; i++)
		printf("%2d ", i);
	printf(">61\n Base 62: ");
	for (i = 31; i < 62; i++)
		printf("%2c ", mha_base62(i));
	printf("  !\n\n");
}

/**********************************************************/
void print_blockhead(int numbl, int totbl)      		/**/
{                                               		/**/
	if (numbl == 1)                           		    /**/
	    printf("   Block %d of %d:\n", numbl, totbl);   /**/
	else                                        		/**/
	    printf("\n   Block %d of %d:\n", numbl, totbl);	/**/
}                                               		/**/
/**********************************************************/

/*****************************************************************************************/
short unsigned int recoverlen(void)
{
int m=0, n=0;
int alpha_count=0;
char  blnk = Fill->sym;
int  width = Current.pass_W;
int height = Current.pass_H;
char letr;

	for (m=0; m<height; m++) {
		for (n=0; n<width; n++) {
			while (align2D[m][n]==blnk) {
				n++;
			}
			if (align2D[m][n]==monoL.sym)
				n++;

			while(isalpha(align2D[m][n])) {
				n++;
				alpha_count++;
			}

			if (m<height-1 && ((letr=align2D[m][n]) != slip.sym && letr!=monoR.sym)) {
				align2D[m][n  ] = slip.sym;
				align2D[m][n+1] = '\0';
				if (dev_print(LOGY,__LINE__)) {
					printf("recoverlen() is adding a missing line terminator at m=%d, n=%d.", m,n);
				}
				break;
			}
			else if (m==height-1 && align2D[m][n] != Term->sym) {
				align2D[m][n  ] = Term->sym;
				align2D[m][n+1] = '\0';
				if (dev_print(LOGY,__LINE__)) {
					printf("recoverlen() is adding a missing final terminator at m=%d, n=%d.", m,n);
				}
				return(alpha_count);
			}
			else if (align2D[m][n] == Term->sym)
				return(alpha_count);
			else
				break;
		}
	}
	return(alpha_count);
}


/*****************************************************************************************/
void usage(char *usage_version)
{
	int i;

	printf("Running lineup version %s, a program for micro-paralogy self-alignment.\n", usage_version);
	printf("\nUsage: ./lineup -[OPTIONS] sequence.txt");
	for (i=1; i<53; i++) {
		if (*Options[i]->description)
			printf("\n\t\t -%c\t%s",Options[i]->sym, Options[i]->description);
	}
	printf("\n\n");
	printf("  *   This option takes an optional number argument with a monotonically-increasing effect.\n");
	printf("  *x  This option takes an optional number argument specifying which mode to use among a set of size x.\n");
	printf("\n Usage examples: ./lineup -v1 sequence_file.txt            [Run in verbose user mode]");
	printf("\n                 ./lineup -Kn sequence_file.txt -X2 -Y800  [Show consensus; skip relax-2D; use Fisher-Yates to generate 800 bp]\n");
	printf("\n Citations:");
	printf("\n  (1) Erives, A. J. (2018) Genetic sequences are two-dimensional. bioRxiv 2018.");
	printf("\n      https://doi.org/10.1101/299867. CC-BY 4.0 International license.");
	printf("\n  (2) Erives, A. J. (2019) Maximal homology alignment: A new method based on two-dimensional homology. bioRxiv 2019.");
	printf("\n      https://doi.org/10.1101/593228. CC-BY 4.0 International license.\n\n");
}

/*****************************************************************************************/
void warnhead(char l)
{
/*	printf("%c", 7); */	/* BELL CHARACTER */
	printf("\n  * Notice (%c): ", l);
}

/*************************************************************************************************/
short int pushdown(char pusharray[][MAXROW], int push_m, int push_n)
{
int i=0, j=0;
int h = Current.pass_H;		/* height slot */
char letr;
char blank = Fill->sym;

	if (h < MAXROW) {
		/* FIRST PUSH DOWN ROWS UNDER m */
		for (i = h; i > push_m+1; i--) {
			for (j = 0; (letr=pusharray[i-1][j]) != '\0'; j++) {
				pusharray[i][j] = letr;
			}
				pusharray[i][j] = '\0';		/* TO ENSURE TERMINATION OF PUSHED DOWN LINE */
		}

		/*  THEN PUSH DOWN ROW m */
		for (j = 0; j < push_n; j++)
			pusharray[push_m+1][j] = blank;
		for (j = push_n; (letr=pusharray[push_m][j]) != '\0'; j++) {
			pusharray[push_m+1][j] = letr;
		}
			pusharray[push_m+1][j] = '\0';			/* TO ENSURE TERMINATION OF PUSHED DOWN LINE */

		return (1);		/* RETURN SUCCESS */
	}
	else
		return (0);		/* RETURN NOT SUCCESSFUL */
}


/*************************************************************************/
void print1D(void)
{
	int i, j, n;
	int lenseq = Clean.pass_W;
	int blocks = count_wrap_blocks(lenseq, par_wrap.set);
	int head_start = 0;			/* USE TO PASS RULER O-F-F-SET TO line_end() */
	char ch;

	mha_head(lenseq);

	for (j = 0; j < blocks; j++) {
		if (Clean.pass_V==1 && opt_v.val) { /* IF DNA AND VERBOSITY */
			line_end(SLIPS, j+1, 0);
   			for (n = j * par_wrap.set; n < (j+1) * par_wrap.set && tela[n].c!=Term->sym && tela[n].c!='\0'; n++) {
				if ((ch=tela[n].t)=='R' || ch=='Y')
					printf("%c", ch);
				else
					printf("_");
			}
			printf("\n");
		}
		line_end(START, j+1, 0);
   		for (n = j * par_wrap.set; n < (j+1) * par_wrap.set && tela[n].c!=Term->sym && tela[n].c!='\0'; n++) {
			printf("%1c", tela[n].c);
		}
		if (tela[n].c == Term->sym) {
			printf("%1c", tela[n].c);  /* PRINTS TERMINAL CHARACTER */
			if (opt_L.bit)
				line_end(END, n, 0);
			else
				printf("\n");
		}
		else {
			printf(" ");
			if (opt_L.bit)
				line_end(END, n, 0);
			else
				printf("\n");
		}

		if (opt_l.bit) {		   /* OPTION opt_l TO SHOW SLIP LOCATIONS */
			/**********************************************/
			line_end(SLIPS, 0, 0);
			for (i = j * par_wrap.set; i < (j+1) * par_wrap.set && i<lenseq; i++) {
				if (!tela[i].k)
					printf(".");
				else
					printf("%c", mha_base62(tela[i].k));
			}
			if (j+1 == blocks)
				printf(" <==== TR unit size (k in base 62)\n");
			else
				printf("\n");
			/**********************************************/
			line_end(SLIPS, 0, 0);
			for (i = j * par_wrap.set; i < (j+1) * par_wrap.set && i<lenseq; i++) {
				if (!tela[i].r)
					printf(".");
				else {
					printf("%c", mha_base62(tela[i].r));
				}
			}
			if (j+1 == blocks)
				printf(" <==== # of TRs > 1 (r in base 62)\n");
			else
				printf("\n");
			/**********************************************/
			line_end(SLIPS, 0, 0);
			for (i = j * par_wrap.set; i < (j+1) * par_wrap.set && i<lenseq; i++) {
				printf("%c", tela[i].echoes);
			}
			if (j+1 == blocks)
				printf(" <==== (((( {  REVERB  } ))))\n");
			else
				printf("\n");
			/**********************************************/
			head_start = (j * par_wrap.set) % 10;
			if (j+1 < blocks) {
				line_end(SLIPRULER, head_start, par_wrap.set);
				printf("\n");
			}
			else
				line_end(SLIPRULER, head_start, lenseq-par_wrap.set*(blocks-1));
			/**********************************************/
		}   /* END OF opt_l PRINT MODULE */

	} /* END OF FOR j LOOP */
	print_section_spacer();
}

/*************************/
int span_ork(int point)
{
	int product = tela[point].or * (tela[point].ok);
	return(product);
}

/**********************/
int span_rk(int point)
{
	int product = tela[point].r * (tela[point].k);
	return(product);
}

/**********************************/
void free_2D(int **p2D, int lenseq)
{
    int row = 0;
    for (row =0; row < lenseq+1; row++)
    {
        free(p2D[row]); // free allocated memory
    }
    free(p2D);
}

/**************************************************************/
short unsigned int checkfractals_in_imperfect(int kmer, int n)
{
	int lenseq = Clean.pass_W;
	int allowed_transits(int k);
	int theta = opt_b.val;			/* threshold for needing to check */

	if (kmer<theta || n<theta-1 || n>lenseq-theta-1)
		return(0);
	else {
		int i, k;
		int transits, max_transits;
		for (k = kmer/2; k>theta; k--) {
			max_transits = allowed_transits(k);
			if (kmer%k==0) {
				int m = n-kmer;
				transits = 0;

				for (i=0; i<k; i++) {
					if (tela[m+i].e != tela[m+k+i].e || tela[n+i].e != tela[n+k+i].e)
						break;
					else if (tela[m+i].c == tela[m+k+i].c && tela[n+i].c == tela[n+k+i].c)
						;
					else if (transits<=max_transits)
						transits++;
					else if (transits > max_transits)
						break;
				}
				if (i==k && transits && transits<=max_transits) {
/*					printf("\n checkfractals_in_imperfect: kmer=%d, n=%d, transits=%d, returning 1 (TRUE).", kmer,n,transits);
*/					return(1);			/* FLAG THAT IT IS THE CASE THAT THE PARENT IMPERFECT TR IS BASED ON SUB-FRACTAL TRs */
				}
			}
		}
		return(0);
	}
}

/**************************************************/
/****  SETS UNIFORM SPACING BETWEEN SECTIONS.  ****/
void print_section_spacer(void) {
	printf("\n");
}

#endif		/* !FILE_LOGY_SEEN */

/*************************************************************************************************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal */
/*************************************************************************************************************/
