/******************************************************************/
/* maximal's microhomology.h: header file since mha_v3.97.c       */
/******************************************************************/

#ifndef FILE_LOGY_SEEN
#define FILE_LOGY_SEEN

#define TEMP		2		/* SETS MAGIC MELT TEMPERATURE; SETS LOWER BOUND FOR MELTAGE opt_m; UNELECTED, DEFAULT ACTS AS TEMP=1 */
#define CYCMAX     60       /* SEMI-MAGIC NUMBER; SEARCH MAGIC TO FIND OTHER EMBEDDED DECISIONS */
#define MAXROW   2000       /* maximum input line size; NOT M-A-G-I-C JUST WHAT MY 'PUTERS CAN DO AS CURRENTLY WRITTEN */
#define WIDTH      72       /* BANDWIDTH: MAX WIDTH OF HEMIDIAGONAL OF PATHBOX; MAX TR UNIT SIZE */ 
#define FRAME_ROWS 25       /* NUMBER OF AVAILABLE ROWS FOR STORING OVERLAPPING REPEAT FRAMES; MULT. OF 4 - EXTRA */
#define START       0       /* FOR USE WITH line_end() */
#define END         1       /* FOR USE WITH line_end() */
#define RULER       2       /* FOR USE WITH line_end() */
#define SLIPS       3       /* FOR USE WITH line_end() */
#define PATHBOXHEAD 4       /* FOR USE WITH line_end() */
#define BLOCKHEAD   5       /* FOR USE WITH line_end() */
#define SLIPRULER   6       /* FOR USE WITH line_end() */
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })

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
	int  cyc_F[FRAME_ROWS];	/* cycling frames; count-off column positions per unit; 32 - 7 = 25 */
							/* one row/frame; row 0 is row # locator; FRAME_ROWS IS BASED ON MEM AL. */
	int all_k;		/* ALL SERIES: PRE-CINCH-T: k-MER SIZE                        */
	int all_r;		/* ALL SERIES: PRE-CINCH-T: REPEAT NUMBER                     */
	int all_S;		/* ALL SERIES: PRE-CINCH-T: SUM OF SCORES OVER ALL UNITS      */
	int all_Z;		/* ALL SERIES: PRE-CINCH-T: ALL_S + TIE-BREAKERS			  */
	int all_R;		/* ALL SERIES: PRE-CINCH-T: POSITION OF CONFLICTING TR ON RHS */
	int all_L;		/* ALL SERIES: PRE-CINCH-T: POSITION OF CONFLICTING TR ON LHS */
	char stat;		/* ALL SERIES: PRE-CINCH-T: STATUS                            */
	/*************************************************************************************************/
	int cyc_Lf;		/* Left-side overlapping TR; 0=lenseq, which is also stored in options[1][1]     */
	int cyc_Rt;		/* Right-side overlapping TR */
	char cyc_o;		/* x => cinched; o => untaken cyclelizable option; !,** => CHECK_TELA VIOLATIONS */
	/*************************************************************************************************/
} tela[MAXROW];

struct stropt {
	int		bit;				/* ON/OFF bit switch    */
	int 	val;				/* incrementable value  */
	char	sym;				/* option letter symbol */
	char	description[64];	/* usage description    */

}	*stroptions[53] = {},
                     /* ".........|.........|.........|.........|.........|.........|...X */	
	opt_a = {0, 0, 'a', {0}},	
	opt_b = {0, 0, 'b', {0}},
	opt_c = {0, 0, 'c', 	"Show base 62 single character code used for k-size and number. "},
	opt_d = {0, 0, 'd',		"Skip cinch-d module. Also automatically skips relax-2D module."},
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
	opt_u = {0, 0, 'u', 	"Do not wrap 2-D output (wrap in one block or other).*          "},
	opt_v = {0, 0, 'v', 	"Increment verbosity level.***                                  "},
	opt_w = {0, 0, 'w', {0}},
	opt_x = {0, 0, 'x', 	"Reduce thresholds for transition matching by one.*             "},
	opt_y = {0, 0, 'y', {0}},
	opt_z = {0, 0, 'z', 	"Set mismatch score from default to zero (visible in pathbox).  "},
	opt_A = {0, 0, 'A', {0}},
	opt_B = {0, 0, 'B', 	"Reset blank character (default '.') to space (' ').****        "},
	opt_C = {0, 0, 'C', 	"Use the reverse complement (- strand) of a DNA sequence.       "},
	opt_D = {0, 0, 'D', 	"Turn on program development prompts (user-controlled pauses).  "},
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
	opt_O = {0, 0, 'O', 	"Append 2-D alignment to consensus file or raw input file.**    "},
	opt_P = {0, 0, 'P', 	"Show pathbox.                                                  "},
	opt_Q = {0, 0, 'Q', {0}},
	opt_R = {0, 0, 'R', 	"Recover sequence from 2-D; toggled on by '-v' and off by '-vR'."},
	opt_S = {0, 0, 'S', {0}},
	opt_T = {0, 0, 'T', 	"Show diagonal threshold values for transitions scoring.        "},
	opt_U = {0, 0, 'U', {0}},
	opt_V = {0, 0, 'V', {0}},
	opt_W = {0, 0, 'W', {0}},
	opt_X = {0, 0, 'X', 	"Scramble sequence with basic or Fisher-Yates randomization.**  "},
	opt_Y = {0, 0, 'Y', 	"Set Fisher-Yates length (FY_size) specified by run argument.   "},
	opt_Z = {0, 0, 'Z', {0}};
                     /* ".........|.........|.........|.........|.........|.........|...X */	

long int options[2][64] = {
/* 0 0 0 0 0 0 0 0 0 0 1  1 1 1 1 1 1 1  1 1 2 2 2 2 2 2  2   2  2 2 3 3 3 3 3  3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6   
   0,1,2,3,4,5,6,7,8,9,0, 1,2,3,4,5,6,7, 8,9,0,1,2,3,4,5, 6,  7, 8,9,0,1,2,3,4, 5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3   
   o>c>t>l>k>n>d>r>R   A  B C D E F G H  I J K L M N O P  Q   R  S T U V W X W  Z a b c d e f g h i j k l m n o p q r s t u v w x y z - - */
 { 0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,  0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, 
 { 0,0,0,0,0,0,0,0,0,0,0,46,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,40, 41,43,4,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, 
}; /*                     |              |                |   |  |              |                                                           
                         46 = "." (FULLSTOP)             "(" ")" |             Zero tick mark                                               
                         32 = " " (SPACE)|               left/right            default = " " (SPACE)                                        
                         blank character |               run delimiters                                                                     
                                                                 |                                                                  
                                                                 Strand characters 43 = '+', 45 = '-'                               
   options[0][  n] = OFF (0) / ON (1) SWITCH TO INDICATE OPTION DETECTED AS ARGUMENT TO MAXIMAL                                      
   options[1][  n] = THIS VALUE IS INCREMENTED ALONGSIDE THE BIT SWITCH IN ROW ZERO FOR SOME OPTIONS                                 
   options[1][ 18] COUNTER OF INITIAL PASSES THROUGH MHA; BEGINS COUNT WITH VALUE -1                                                 
   options[1][0-9] WILL PERMANENTLY STORE 2-D WIDTH HISTORY AND IS THE ORIGINAL REASON OPTIONS WAS CODED AS LONG INT                 
                0 EQUALS ORIGINAL STRING         (o)                                                                                 
                1 EQUALS cleanseq PASS (1-D)     (c)                                                                                 
                2 EQUALS cinch_t 2-D PASS MAIN() (t)                                                                                 
                3 EQUALS cinch_l 2-D PASS        (l)                                                                                 
                4 EQUALS cinch_k 2-D PASS        (k)                                                                                 
                5 EQUALS nudgelize 2-D PASS      (n)  Some comments might refer to 'cyclelize', its original name and function.      
                6 EQUALS cinch_d  2-D PASS       (d)                                                                                 
                7 EQUALS relax_2D 2-D PASS       (r)                                                                                 
                8 EQUALS recovered 1-D from 2-D  (R)                                                                                 
               10 EQUALS passQ score / 1000      (A)                                                                                 

				 options[1][13] IS RESERVED FOR STORING SEQUENCE TYPE (DNA, RNA, PROTEIN, BABYLONIAN, etc.) 
				 options[1][17] RESERVED FOR 2D-ALIGNMENT HEIGHT 	
				 options[1][46] RESERVED FOR RECORDING LARGEST CINCH-T k-MER UNIT SIZE 
				 options[1][50] RESERVED FOR RECORDING CUMULATIVE BADSLIP TYPE PER RUN 
*/

char align2D[MAXROW][MAXROW] = {{0}};
char pathbox[MAXROW][MAXROW] = {{0}};
char consensus[MAXROW] = {0};
char file_name[255] = "internal_default";
char dev_notes[32] = "-";             /* STRING WRITTEN AS LAST FIELD IN OUTPUT FILE */
short unsigned int cinchled=0;			/* BIT FLAG FOR CINCH-L WRAPS */
char letr_unit[4] = {0};				/* UNIT STRING: "bp" FOR DNA, "nt" FOR RNA, 'aa' FOR PROTEINS, 'ch' FOR ALL OTHER; SET IN MAIN() */
FILE *fp_out;                           /* FILE FOR OUTPUT.LOG */

void 				clear_2D_ar(char wipe_align2D[][MAXROW]);
void 				clear_right(char swipe_align2D[][MAXROW]);
int 				col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown); 
unsigned int 		consensus_2D(int n_start, int n_width);
int 				count_wrap_blocks(int lcl_width, int lcl_opt_w);	/* lcl_width IS WIDTH OF 2-D MHA ARRAY */ 
int 				get_1Dz(int x, int y, int ignoreCheck);
void 				line_end(int type, int c, int lcl_width);
char 				mha_base62(int num);
void 				mha_head(int lcl_width);
void 				mha_UPPERback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW]);
void 				mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW]);
void 				mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW]);
void 				mha_writeconsensus(char align2D_one[][MAXROW], char consensus1D[MAXROW]);
void 				print1D(void);
short unsigned int	print_2Dseq(void);
void 				print_blockhead(int a, int b);	
short int 			pushdown(char pusharray[][MAXROW], int push_m, int push_n); 
int 				span_allrk(int point);
int 				span_rk(int point);
void 				warnhead(char l); 
short unsigned int 	recoverlen(void);
short unsigned int	cleanseq(char *s);
int 				get2Dtucknum(char arrayA[][MAXROW], char arrayB[][MAXROW]);
void 				mha_randomize1(char *input_seq);
void 				mha_randomize2(char *input_seq, int rsize);
void 				print_base62_table(void);
void 				shuffle(int *array, int n);
void 				usage(char *usage_version);	/* FOR PRINTING UNIFORM USAGE INFORMATION */
char 				*nmer_prefix(int i);		/* CONVERTS INTEGER TO N-MER PREFIX WRITTEN NAME */
void 				free_2D(int **p2D, int lenseq);

/*******************************************************************************************************/
 #include "microhom-devl.h"	/* maximal header: program development code, testing, and evaluation       */
/*******************************************************************************************************/

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
	int lenseq = options[1][1];

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
int lenseq = options[1][1];

	for (m=0; m <= lenseq; m++) {
		for (n=0; n <= lenseq; n++)
			wipe_align2D[m][n] = '\0';
	}
}
/*****************************************************************************************/

/*****************************************************************************************/
void clear_right(char swipe_align2D[][MAXROW])
{
int m=0, n=0;
int lenseq = options[1][1];	
int height = options[1][17];
char cropt_R_rght = options[1][27];
char letr;

	/* CLEAR TO THE RIGHT OF ROW TERMINATORS */
	for (m=0; m < height; m++) {
		n = 0;
		while ( (letr=swipe_align2D[m][n]) != '/' && letr != cropt_R_rght && letr != '>') {
			n++;
		}
		for (n = n+1; n <= lenseq; n++)
			swipe_align2D[m][n] = '\0';
		if (letr == '>') {
			for (m=m+1; m < height; m++) {
				for (n=0; n <=lenseq; n++) {
					swipe_align2D[m][n] = '\0';
				}
			}
			break;
		}
	}
}
/*****************************************************************************************/

/*****************************************************************************************/
unsigned int consensus_2D(int n_start, int n_width)
{
int badsites=0, m=0, n=0, n_end, x=1;
int con_width = options[1][32];
short unsigned int nuctype = options[1][13];		/* FOR SEQ TYPE, DNA=1, RNA=2, OTHER (NON-NA)=0 */
short unsigned int nuctransit = 0;					/* BIT FLAG FOR HANDLING NUCLEAR TRANSITIONS */ 
short unsigned int plustransit=0;					/* BIT FLAG ADDENDUM FOR COUNTING BADSITES AT COL */
short unsigned int checktransit=0;					/* BIT FLAG FOR CHECKING GOOD TRANSITION MARK */
char blnk = options[1][11]; 
char letr=blnk, ltr2=blnk, conletr=blnk;
int con_maxrows=26;
int consensus_ar[26][MAXROW] = {{0}};	 	/* COL n=0 FOR BIT FLAG */
                                       		/* ROW m=0 FOR COUNTER */
											/* ROW m=1 FOR CONSENSUS */
											/* ROWS m>1 FOR VARIANTS STORAGE */
	n_end = n_start + n_width;

	if (nuctype == 1)	/* IF DNA */
		nuctransit = 1;

	/* FILL CONSENSUS ARRAY WITH FIRST LETTER IN EACH COL */
	for (n = 0; n <= con_width; n++) {
		for (m = 0; m < MAXROW; m++) {
			if ( (isalpha(letr=align2D[m][n]))) {
				consensus_ar[0][n+1]++;
				if ( letr <= 'Z' || (nuctype && letr=='n') )  
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
			else conletr = blnk;
		}

		for (m = 0; align2D[m][0] != '\0'; m++) {
			if (isalpha(letr=align2D[m][n])) {
				ltr2=toupper(letr);
				if (nuctransit && ((conletr=='R' && (ltr2=='G' || ltr2=='A')) || (conletr=='Y' && (ltr2=='C' || ltr2=='T'))) ) {
					if (conletr=='R') {
						if (checktransit==18) {
							if (ltr2=='A')
								checktransit = checktransit -  1;
							else if (ltr2=='G')
								checktransit = checktransit - 17;
						}
						else if (checktransit==17 && ltr2=='G')
							checktransit = 0;
						else if (checktransit== 1 && ltr2=='A')
							checktransit = 0;
					}
					if (conletr=='Y') {
						if (checktransit==25) {
							if (ltr2=='C')
								checktransit = checktransit -  3;
							else if (ltr2=='T')
								checktransit = checktransit - 22;
						}
						else if (checktransit==22 && ltr2=='T')
							checktransit = 0;
						else if (checktransit== 3 && ltr2=='C')
							checktransit = 0;
					}
				}
				else if (letr!='n' && letr != consensus_ar[1][n+1]) {
					if (consensus_ar[0][n+1] == 1) {	/* IF THIS IS THE FIRST CONFLICT NOTED AT THIS POSITION */
						++badsites;						/* COUNT ADDITIONAL BAD SITE COLUMN */
						++consensus_ar[0][n+1];			/* INCREMENT COUNTER FOR DIFFERENT LETTERS AT COLUMN */
						consensus_ar[2][ 0 ] = 1;		/* TURN VARIANT ROW ON */
						consensus_ar[2][n+1] = letr;
					}
					else {		/* ELSE CHECK IF THIS IS A NEW LETTER AT THIS COLUMN */	
						for (x = 2; x < con_maxrows; x++) {
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
				else if (x >= con_maxrows)
					break;
			}
			else if (letr == '>') {
				while (col_isclear(align2D,n,m,-1)>-1)
					n++;
				consensus_ar[1][n+1] = letr;
				options[1][17] = m+1;		/* STORE HEIGHT IN HEIGHT SLOT */
			}
		} 
		if (checktransit) {		/* IF checktransit IS STILL POSITIVE THEN IT'S SUPPOSED TO NOT HAVE CHECKED OUT */
			consensus[n] = letr;
			plustransit = 0;
			if (options[1][18]) { /* IF PASS NUMBER */
				options[1][39] = 1; sprintf(dev_notes, "checktransit=%d at n=%d", checktransit, n);
				if (dev_print(LOGY,__LINE__)) {
					printf("checktransit=%d at n=%d.\n", checktransit, n);
				}
			}
		}
		else if (plustransit)
			++consensus_ar[0][n+1];
	}

	if (opt_K.bit) {							/* opt_K SHOW CONSENSUS ROW */	
		/* PRINT CONSENSUS ROWS */
		for (m = 1; consensus_ar[m][0] != '\0' && m < con_maxrows; m++) {
			line_end(BLOCKHEAD, 9, 9);
			printf(" ");								/* OPEN CONSENSUS PARENTHESES */
			for (n = n_start; n < n_end; n++) {
				if (isalpha(letr=consensus_ar[m][n+1]))
					printf("%c", letr); 
				else
					printf(".");
			}
			if (badsites == 0)
				printf(" = consensus\n");
			else
				printf("\n");
		}
		if (badsites > 0) {
			line_end(BLOCKHEAD, 9, 9);	
			printf(" ");
			for (n = n_start; n < n_end; n++) {
				letr = mha_base62(consensus_ar[0][n+1]);
				printf("%c", letr); 
			}
		    printf(" <= # (base 62)\n");
			warnhead('C');
			if (badsites == 1)
				printf("There is one site with mimatched letters.\n");
			else if (badsites > 1)
				printf("There are %d sites with mismatched letters.\n", badsites);
		}
	} 

	/* STORE CONSENSUS ROW */
	for (n = n_start; n <= n_end; n++) {
		consensus[n] = consensus_ar[1][n+1];
	}
	consensus[con_width] = '\0';

	if (badsites > 0) {
		return(badsites);				/* BAD CONSENSUS, REPORT IT 	*/
	}
	else								/* GOOD CONSENSUS, EARLY PASSES */
		return(0);
}
/*****************************************************************************************/

/*****************************************************************************************/
unsigned int connudge(char con_align2D[][MAXROW], int n_start, int n_width)
{
int badsites=0, m=0, n=0, n_end, x=1, nudge_row=0, nudge_col=0, nudge_span=0, frstletr;
int    lenseq = options[1][ 1];
int    height = options[1][17];
int con_width = options[1][32];
short unsigned int nuctype = options[1][13];	/* FOR SEQ TYPE, DNA=1, RNA=2, OTHER (NON-NA)=0 */
short unsigned int nuctransit = 0;					/* BIT FLAG FOR HANDLING NUCLEAR TRANSITIONS */ 
short unsigned int plustransit=0;					/* BIT FLAG ADDENDUM FOR COUNTING BADSITES AT COL */
short unsigned int checktransit=0;
char blnk = options[1][11]; 
char letr=blnk, conletr=blnk, chkletr=blnk, badletr=blnk;
int con_maxrows=26;
int consensus_ar[26][MAXROW] = {{0}};	 	/* COL n=0 FOR BIT FLAG */
                                       		/* ROW m=0 FOR COUNTER */
											/* ROW m=1 FOR CONSENSUS */
											/* ROWS m>1 FOR VARIANTS STORAGE */
	n_end = n_start + n_width;

	if (nuctype == 1)	/* IF DNA */
		nuctransit = 1;

	/* FILL CONSENSUS ARRAY WITH FIRST LETTER IN EACH COL */
	for (n = 0; n <= con_width; n++) {
		for (m = 0; m < lenseq; m++) {
			if ( (isalpha(letr=con_align2D[m][n]))) {
				consensus_ar[0][n+1]++;
				consensus_ar[1][n+1] = letr;
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
			if ((conletr=con_align2D[MAXROW][n])=='R' || conletr=='Y') {
				consensus_ar[1][n+1] = conletr;
				chkletr = blnk;
				plustransit = 1;
			}
			else {
				plustransit = 0;
				conletr = blnk;
			}
		}

		for (m = 1; con_align2D[m][0] != '\0'; m++) {
			if (isalpha(letr=con_align2D[m][n]) && isupper(letr)) {
				if (plustransit && ((conletr=='R' && (letr=='G' || letr=='A')) || (conletr=='Y' && (letr=='C' || letr=='T'))) ) {
					if (chkletr==blnk) {
						chkletr = letr;
						plustransit = 2;
					}
					else if (letr!=chkletr) {	
						plustransit = 3;	/* THREE LETTERS SEEN: R-A-G OR Y-C-T */
					}
				}
				else if (letr != consensus_ar[1][n+1]) {
					if (consensus_ar[0][n+1] == 1) {	/* IF THIS IS THE FIRST CONFLICT NOTED AT THIS POSITION */
						++badsites;						/* COUNT ADDITIONAL BAD SITE COLUMN */
						badletr=letr;					/* THIS IS THE LETTER THAT GOES AGAINST THE TRANSITION LETTER */
						++consensus_ar[0][n+1];			/* INCREMENT COUNTER FOR DIFFERENT LETTERS AT COLUMN */
						consensus_ar[2][ 0 ] = 1;		/* TURN VARIANT ROW ON */
						consensus_ar[2][n+1] = letr;
						if (nudge_row==0) {
							nudge_row=m;
							nudge_col=n;
							if (nuctransit && (conletr=='R' || conletr=='Y')) {
								if      (conletr=='R') {
									while (nudge_row <= height && con_align2D[nudge_row][nudge_col]!='A' && con_align2D[nudge_row][nudge_col]!='G')
										++nudge_row;
								}
								else if (conletr=='Y') {
									while (nudge_row <= height && con_align2D[nudge_row][nudge_col]!='C' && con_align2D[nudge_row][nudge_col]!='T')
										++nudge_row;
								}
								
								if (nudge_row > height)	/* REVERT IF MATCH NOT FOUND */
									nudge_row=m;
							}
						}
					}
					else {		/* ELSE CHECK IF THIS IS A NEW LETTER AT THIS COLUMN */	
						for (x = 2; x < con_maxrows; x++) {
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
				else if (x >= con_maxrows)
					break;
			}
			else if (letr == '>') {
				consensus_ar[1][n+1] = letr;
			}
		} /* END OF FOR m=1... */
		if (plustransit) {
			if (plustransit>1) {
				++consensus_ar[0][n+1];
				if (plustransit==2)
					consensus_ar[1][n+1] = chkletr;
			}
			else if (plustransit==1)
					consensus_ar[1][n+1] = badletr;
		}
	} /* END OF FOR n... */

	m=1;
	n=0;
	while ((consensus_ar[0][m]==1) || (consensus_ar[0][m]==2 && ((letr=consensus_ar[1][m])=='R' || letr=='Y')) ) {
		m++;
	}
	if (consensus_ar[0][m]==0) {
		return(0);
	}
	else {
		while (consensus_ar[0][m+n]!=1 && m+n <= con_width+1) {
			n++;
		}
		if (consensus_ar[0][m+n]==1) {
			if (m == nudge_col+1) {
				nudge_span = n;
			}
		}
	}

	/* CHECK FOR A NON-PRODUCTIVE NUDGE THAT WOULD CREATE AN EMPTY COLUMN */
	for (n = 0; !isalpha(con_align2D[nudge_row][n]) && n <= con_width; n++) {
		;
	}
	if (isalpha(con_align2D[nudge_row][n]))
		frstletr = n;
	else
		frstletr = -1;

	if (frstletr > -1) {
		while (con_align2D[nudge_row-1][frstletr  ]==con_align2D[nudge_row][frstletr] &&
			   con_align2D[nudge_row-1][frstletr-1]==blnk &&
			   con_align2D[nudge_row-1][frstletr+1]== '/') {
			--nudge_row;
		}
	}

	if (col_isclear(con_align2D, frstletr, nudge_row,-1)<0) {
		options[1][39]=4;
		return(0);
	}
	/* ****************************************************************** */

	/* NUDGE IT! */
	for (m=nudge_row; con_align2D[m][0] != '\0'; m++) {
		for (n = n_end+2; n > 0; n--) {
			con_align2D[m][n] = con_align2D[m][n-1];
		}
		con_align2D[m][0] = blnk;
	}

	/* ASSIGN CONSENSUS ROW LETTERS BUT VERIFY POST-NUDGE */
	for (n = n_end+2; n>=nudge_col+nudge_span; n--) {
		con_align2D[MAXROW][n] = consensus_ar[1][n];
	}
	con_align2D[MAXROW][n] = blnk;

	options[1][32]++;
	if (checktransit) {
		return(0);
	}
	else
		return(1);				/* RETURN SUCCESS */
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
				foam_ar[1][n] = '.';
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
	int m = 0;

	/* INVOKE updown AS -1 OR 1 TO CHECK BELOW OR ABOVE row, RESPECTIVELY */
    for (m = row + updown; m >= 0; m = m + updown) {
       	if (isalpha(check_array[m][at_n]))
    		return(m);
		if (check_array[m][0] == '\0') 
			break;
	}
    return(-1);		/* RETURN NEGATIVE VALUE TO INDICATE COLUMN IS CLEAR EITHER BELOW OR ABOVE */
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
	char zero_tick = (char) options[1][35];	/* opt_Z, Zero tick mark, default = 32 = ' ' */

	if (opt_B.val==2 && options[1][18]>1)	/* opt_B LEVELS FOR BLANKNESS IN FILLER & opt_I PASS NUM */
		zero_tick = '|';			/* 124 = '|'							*/

	if (type == 1)					/* FORMAT FOR LINE END. c IS CHAR. NUMBER */
		printf("%5d\n", c);

	if (opt_r.bit) {				/* FORMAT WITH LINE NUMBERING. DEFAULT */
		if (type == 0) {			/* FORMAT FOR LINE BEGINNING. c is LINE NUMBER (m+1) */
			if (c == 1)
				printf(" %4d. >", c);
			else 
				printf(" %4d. %c", c, zero_tick);
		}
		else if (type == 2) {			/* FORMAT FOR RULER. */
			printf("       %c%.*s\n", zero_tick, lcl_width, ruler + c);		/* DOUBLE USE OF c AS O-F-FSET FOR RULER */
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
				printf("   >");
			else           
				printf("   %c", zero_tick);
		}
		else if (type == 2) {			/* FORMAT FOR RULER OR SLIPRULER. */
			printf("   %c%.*s\n", zero_tick, lcl_width, ruler + c);		/* DOUBLE USE OF c AS O-F-F-SET FOR RULER */
		}
		else if (type == 3)			/* FORMAT FOR UNNUMBERED SLIP LOCATION LINES. */
			printf("    ");
		else if (type == 4)			/* FORMAT FOR PATHBOX HEADERS, STILL #'d */
			printf("      _|");
		else if (type == 5)			/* FORMAT FOR CONSENSUS BLOCKHEAD HEADERS, #'d */
			printf("   ");
		else if (type == 6) {		/* FORMAT FOR SLIPRULER. */
			ruler = rule2;
			printf("   %c%.*s\n", zero_tick, lcl_width, ruler + c);		/* DOUBLE USE OF c AS O-F-F-SET FOR RULER */
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
			 /*0123456789*/
char h1[]=	"\n\n_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//"
			    "_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//"
			    "_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//"
			    "_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//";
char h2[]=	"\n\n\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			    "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			    "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			    "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_";
char *h_rule = h1;						/* DEFAULT BANNER STYLE */
int min_len = 80;						/* MINIMUM LENGTH */
int med_len = 12*(lcl_width/10);		/* MEDIUM LENGTH, SCALING */
int max_len = options[1][58]+8;			/* MAXIMUM LENGTH */
int hr_len = min_len;					/* DEFAULT LENGTH OF HEADER BANNER */
int lcl_pass = options[1][18];			/* opt_I VALUE COUNTER FOR NUM OF PASSES */

	if (lcl_width > MAXROW && dev_print(LOGY, __LINE__)) {
		printf("Bad news bears: Unexpectedly, lcl_width > MAXROW. lcl_pass=%d, lcl_width=%d.\n", lcl_pass, lcl_width);
	}

	if (lcl_width+8 > hr_len && lcl_width+8 < (int) options[1][58])
		hr_len = med_len;
	else if (lcl_width+8 >= (int) options[1][58])
		hr_len = max_len;

	if (opt_X.bit)
		h_rule = h2;

	if (lcl_pass == 6) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: cinch-d [RESCUES TR's INTERRUPTED BY DE NOVO REPEATS OF REPEATS] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 5) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: nudgelize [RESCUES TR's OBSCURED BY INITIAL CYCLING FRAME] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 4 || lcl_pass == 6) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: cinch-k [INTRA-REPEAT *k-MERS > 0] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 3) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: cinch-l [*lONG HOMOPOLYMER TRACTS] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 2) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: cinch-t [LONGEST *tANDEM REPEATS, k>1] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 1) {
			printf("\n\n1-D sequence:\n");
	}
	else if (lcl_pass == 7) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: relax-2D [RELAXES HOMOPOLYMER RUNS THAT DID NOT AID cinch-d] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 0) {
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
	int seqlen = options[1][1];
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
	int seqlen = options[1][1];
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
void mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW])
{
char letr;
int lenseq       = options[1][ 1];
char monoL = options[1][26];		/* LHS character delimiter for homopolymer Run */
char monoR = options[1][27];		/* RHS character delimiter for homopolymer Run */
int m=0, n=0, widest_n=0;

	/* WRITE BACK TO align2D_prev */
	clear_right(lcl_align2D); 
	clear_2D_ar(align2D_prev);
	for (m = 0; lcl_align2D[m][0] != '\0' && m <= lenseq; m++) {
		for (n = 0; (letr=lcl_align2D[m][n]) != '/' && letr != '>' && letr != monoR; n++) {
			align2D_prev[m][n] = letr;
			if (letr == monoL && lcl_align2D[m][n+opt_M.val+1] == monoR ) 
				align2D_prev[m][n+opt_M.val+2] = '\0';
		}
		align2D_prev[m][n  ] = letr;	/* MHA-STANDARD TERMINATOR  */
		if (n > widest_n)
			widest_n = n;
		if (letr == '>') {				
			options[1][32] = widest_n;	/* ASSIGN 2-D WIDTH  */
			options[1][17] = m+1;		/* ASSIGN 2-D HEIGHT */
		}
	}
}

/*****************************************************************************************/
void mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW])
{
int i=0;

	/* COPY CONSENSUS ROW FROM ARRAY ONE TO ARRAY TWO */
	for (i = 0; i < options[1][32]; i++) {
		align2D_two[MAXROW][i] = align2D_one[MAXROW][i];
	}
}

/*****************************************************************************************/
void mha_writeconsensus(char align2D_one[][MAXROW], char consensus1D[MAXROW])
{
int i=0;

	/* COPY CONSENSUS ROW FROM ARRAY ONE TO 1D CONSENSUS ARRAY */
	for (i = 0; i < options[1][32]; i++) {
		consensus1D[i] = align2D_one[MAXROW][i];
	}
}

/*****************************************************************************************/
void mha_UPPERback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW])
{
char letr;
short unsigned 
     int nuctype = options[1][13];				/* EQUALS ONE IF DNA, TWO IF RNA */
char wopt_Q_left = options[1][26];				/* LHS character delimiter for homopolymer Run */
char wopt_R_rght = options[1][27];				/* RHS character delimiter for homopolymer Run */
int lenseq = options[1][1];
int i=0, m=0, n=0, widest_n=0;

	/* WRITE BACK TO align2D_prev */
	clear_right(lcl_align2D); 
	clear_2D_ar(align2D_prev);		/* DOES NOT CLEAR MAXROW ROW */
	for (m = 0; lcl_align2D[m][0] != '\0'; m++) {
		for (n = 0; (letr=lcl_align2D[m][n]) != '/' && letr != '>' && letr != wopt_R_rght; n++) {
			if (nuctype) { /* IF NON-ZERO LIKE DNA (1) OR RNA (2) */
				if (islower(letr) && letr != 'n')
					align2D_prev[m][n] = toupper(letr);
				else
					align2D_prev[m][n] = letr;
			}
			else if (islower(letr))
				align2D_prev[m][n] = toupper(letr);
			else
				align2D_prev[m][n] = letr;

			if (letr == wopt_Q_left && lcl_align2D[m][n+opt_M.val+1] == wopt_R_rght ) 
				align2D_prev[m][n+opt_M.val+2] = '\0';
		}
		align2D_prev[m][n  ] = letr;	/* MHA-STANDARD TERMINATOR  */
		if (n > widest_n)
			widest_n = n;
		if (letr == '>') {				/* ASSIGN CINCH-WIDTH TO CURRENT [32]	*/
			options[1][32] = widest_n;
		}
	}

	/* FLUSH POST-TERMINATOR ENDS OF align2D WITH NULLS */
	for (m = 0; lcl_align2D[m][0] != '\0'; m++) {
		n = 0;
		while ( (letr=lcl_align2D[m][n]) != '/' && letr != '>' && letr != wopt_R_rght) 
			n++;	
		for (i = 1; n+i <= lenseq; i++)
			align2D_prev[m][n+i] = '\0';	
	}
}

/*****************************************************************************************/
char *nmer_prefix(int p)
{
	static char *prefix[22] = { "X", 
		"Mono", "Di", "Tri", "Tetra", 
		"Penta", "Hexa", "Hepta", 
		"Octa", "Nona", "Deca",
		"11-", "12-", "13-", "14-", "15-", 
		"16-", "17-", "18-", "19-", "20-", "21-"
	};

	if (p == 0 || p > 21) 
		return prefix[0];

	else 
		return prefix[p];
}

/*****************************************************************************************/
short unsigned int print_2Dseq(void)
{
unsigned int foam_2D(int n_start, int n_width);
int all_clear;		/* COUNTER VARIABLE USED FOR CHECKING NEED TO PRINT BOTTOM ROWS */ 
int blocks2D=0, b=0, c=0, carry_over=0, d=0, fudge=0, g, h, i, j, j_start, j_end, m, m_start=0, n;
int mmsites=0, max_n=0;
char letr = 'B', next = 'B';			/* Begin the Beguine */
char blnk  = options[1][11];			/* opt_B blank character */
int cinchwidth = (int) options[1][32];		
int cip_linewidth = options[1][58];
char popt_Q_left = options[1][26]; 	/* LHS character delimiter for homopolymer Run */
char popt_R_rght = options[1][27];		/* RHS character delimiter for homopolymer Run */
int   lenseq = options[1][1];
int   height = options[1][17];
int head_start;							/* USE TO PASS RULER O-F-F-SET TO line_end() */
int scrimmageline;						/* USE TO INCREMENT AND TEST IF FILLER IS NEEDED, CAN BE OPTION TO DO SO */
char tick = ':'; 						/* OTHER POSSIBILITIES: |, ^ */
short unsigned int lcl_opt_F;

	if (opt_B.val > 1)
		tick = blnk;

	if (cinchwidth > MAXROW && dev_print(LOGY,__LINE__)) {
		printf("Bad news bears: Unexpectedly, cinchwidth > MAXROW. cinchwidth=%d", cinchwidth);
	}
	if (cinchwidth > MAXROW && dev_print(LOGY,__LINE__)) {
		printf("Bad news bears: Unexpectedly, cinchwidth > MAXROW. cinchwidth=%d.", cinchwidth);
	}
	mha_head(cinchwidth);

	blocks2D = count_wrap_blocks(cinchwidth, cip_linewidth);

	for (j = 0; j < blocks2D; j++) {
		if (blocks2D != 1)
			print_blockhead(j+1, blocks2D);

		/* FAST FORWARD TO INFORMATIVE ROW */
		m_start = 0;
		while (align2D[m_start][(j * options[1][58])] == '\0') {
			m_start++;
		}

		scrimmageline = 1;
				max_n = 0;

		for (m = m_start; align2D[m][0] != '\0'; ) {

			b = d = 0;		/* VAR b WILL COUNT BLANKS, VAR d WILL COUNT ALPHA CHAR ANEW FOR EACH ROW */

			j_start =   j   * options[1][58];
			j_end   = (j+1) * options[1][58];

			for (n = j_start; n < j_end && (letr=align2D[m][n])!='/' && letr!='>' && letr!=popt_R_rght; n++) {
				if (isalpha(letr)) {
					c++;
					d++;
					if (d+b > scrimmageline)
						scrimmageline++;
				}
				else if (letr == blnk || letr == popt_Q_left) {
					b++;
				}
			} /* END OF n SCAN LOOPS */
			if (b+d > max_n)
				max_n = b+d;	

			if (b == cip_linewidth) {
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

			for (n = j_start; n < j_end  && (letr=align2D[m][n])!='>' && letr!=popt_R_rght && letr!='\0'; n++) {
				if (d+b == 0) {		/* HANDLES BLOCK CONTINUATION LINES THAT HAVE ZERO CHARACTERS B/C OF ITS TUCK & LENGTH */
					break;
				}
				else if ((letr=align2D[m][n]) == blnk && (n+1) % 10 == 0 && opt_B.val < 2)		
					printf("%c", tick);		/* PRINT TICK MARKS AT 10 bp INTERVALS IF NOT BLANK SPACE */
				else if (options[1][11] == 32)		/* opt_B BLANK = SPACE ' ' */	
					printf("%c", letr);
				else					/* opt_B BLANK = FULLSTOP '.' */
					printf("%c", letr);

				/* TURN ON OPTION opt_F TO USE BLANK CHAR TO FILL < SCRIMMAGE IF EXCESSIVELY SHORT */
				if (!opt_F.bit && align2D[m][n] == '/' && scrimmageline-(b+d-2) > 10)
					lcl_opt_F = 1; 	/* TURN ON local opt_F IN THIS CASE */

				/* OPTION opt_F TO USE BLANK CHAR TO FILL < SCRIMMAGE */
				if ((opt_F.bit || lcl_opt_F) && align2D[m][n] == '/' && b+d-2 < scrimmageline) {
					for (i = 0; b+d+i < scrimmageline; i++) {
						if ( (b+d+2+i)%10 == 0)		/* WHY 2? +1 FOR STARTING AT 0, +1 FOR '/' CHAR */
							printf("%c", tick);		/* PRINT TICK MARKS AT 10 bp INTERVALS IF NOT BLANK SPACE */
						else
							printf("%c", blnk);		/* PRINT LINE-END FILLER CHARACTER */
					}
					lcl_opt_F = 0;
				}
			} /* END OF n PRINTING LOOPS */

			next = align2D[m][n];	/* TERMINAL CHARACTERS WILL ALSO PRINT AT END OF BLOCK IF...	*/ 
										/* ...THEY ARE FIRST CHARACTER OF NEXT BLOCK					*/

			if (next == '/' || next == popt_R_rght)		/* TO INDICATE ADJACENCY TO SLIP, WHICH ALSO WILL SHOW IN NEXT BLOCK */
				printf("%c", next);
			else if (next==blnk || (isalpha(next) && d+b != 0))	/* TRUE IF AT EDGE OF opt_w WINDOW */
				printf("=>");			/* TO INDICATE CONTINUATION TO NEXT BLOCK */
			else if (next == '\0' && n == j_start)		
				printf("%c", blnk);
			else if (next == '>') {		/* LAST CHARACTER */
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

		options[1][17] = m+1;	/* ASSIGN COUNTED HEIGHT (# OF ROWS) TO HEIGHT SLOT */

		/* PRINT RULER */
		head_start = (j * options[1][58]) % 10;
		if (j+1 < blocks2D) {
			line_end(RULER, head_start, cip_linewidth);
		}
		else {
			line_end(RULER, head_start, (cip_linewidth=max_n));
		}
		/* *********** */

		/* PRINT NUMBERS FOR CONSENSUS RULER */
		if (opt_B.val < 4) {
			line_end(SLIPS, head_start, cip_linewidth);
			fudge = carry_over % 10;
			printf("%.*s", 9-fudge, "          ");
			for (n = 1; n <= cip_linewidth+fudge; n++) {
				if ( (n+carry_over-fudge) % 10 == 0)
					printf("%-10d", n + carry_over-fudge);
			}
			carry_over = carry_over + n - fudge - 1;
			printf("\n");
		}
		/* *********************** */

		/* ADD TO COUNT OF MISMATCHED SITES */
		mmsites = mmsites + consensus_2D(j_start, cip_linewidth);

		if (opt_f.bit && options[1][18] > 6 && mmsites == 0) {	
			foam_2D(j_start, cip_linewidth);
		}
	} /* END OF FOR j PRINTING LOOP */

	if (c == lenseq && mmsites == 0) {
		options[0][10] = 1000;	
		printf("\n This %d %s sequence was auto-aligned correctly at this stage.\n", c, letr_unit);
		return(0);
	}
	else if (c == lenseq) {
		i=options[0][10] = round((1000*(cinchwidth-mmsites))/cinchwidth);	
		printf("\n This %d %s sequence was not auto-aligned correctly at this stage, but perhaps further cinching will rectify the 2-D alignment. \n", c, letr_unit);
		return(i);
	}
	else if (c < lenseq) {
		i=options[0][10] = round((1000*(cinchwidth-mmsites-(lenseq-c)))/cinchwidth);
		warnhead('-');
		printf(" 2-D auto-alignment is missing %d %s(s)!\n\n", lenseq-c, letr_unit); 

		/* ENSURE THERE ARE NO MISSING MHA ROW TERMINATORS */
		for (m=0; m<height; m++) {
			for (n=0; n<cinchwidth; n++) {
				while (align2D[m][n]==blnk) {
					n++;
				}
				if (align2D[m][n]==popt_Q_left) {
					n++;
				}
				while (isalpha(align2D[m][n])) {
					n++;
				}
				if (m<height-1 && (letr=align2D[m][n]) != '/' && letr!=popt_R_rght) {
					align2D[m][n  ] = '/';
					align2D[m][n+1] = '\0';
					if (dev_print(LOGY,__LINE__)) {
						printf("print_2Dseq() adding missing line terminator at row=%d, col=%d.", m,n);
					}
					break;
				}
				else if (m==height-1 && align2D[m][n] != '>') {
					align2D[m  ][n] = '>';
					align2D[m+1][0] = '\0';
					if (dev_print(LOGY,__LINE__)) {
						printf("print_2Dseq() adding missing final terminator at row=%d, col=%d.", m,n);
					}
					break;
				}
				else {
					break;
				}
			}
		}
		return(0);
	}
	else {
		i=options[0][10] = round((1000*(cinchwidth-mmsites-(c-lenseq)))/cinchwidth);	
		warnhead('+');
		printf(" 2-D auto-alignment contains an extra %d %s(s)!\n\n", c-lenseq, letr_unit);
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

/*****************************************************************************************/
void print_blockhead(int a, int b)              /**/
{                                               /**/
	if (a == 1)                                 /**/
	    printf("   Block %d of %d:\n", a, b);   /**/
	else                                        /**/
	    printf("\n   Block %d of %d:\n", a, b); /**/
}                                               /**/
/**************************************************/

/*****************************************************************************************/
short unsigned int recoverlen(void) 
{
int m=0, n=0;
int alpha_count=0;
char  blnk = options[1][11];
int  width = options[1][32];
int lenseq = options[1][ 1];
int height = options[1][17];
char monoL = options[1][26];
char monoR = options[1][27];
char letr;

	/* CHECK BASIC ASSUMPTIONS */
	if (lenseq>MAXROW || height>MAXROW || width>lenseq || height>lenseq) {
		return(0);
	}

	for (m=0; m<height; m++) {
		for (n=0; n<width; n++) {
			while (align2D[m][n]==blnk) {
				n++;
			}
			if (align2D[m][n]==monoL) {
				n++;
			}
			while(isalpha(align2D[m][n])) {
				n++;
				alpha_count++;
			}
			if (m<height-1 && ((letr=align2D[m][n]) != '/' && letr!=monoR)) {
				align2D[m][n  ] = '/';
				align2D[m][n+1] = '\0';
				if (dev_print(LOGY,__LINE__)) {
					printf("recoverlen() is adding a missing line terminator at m=%d, n=%d.", m,n);
				}
				break;
			}
			else if (m==height-1 && align2D[m][n] != '>') {
				align2D[m][n  ] = '>';
				align2D[m][n+1] = '\0';
				if (dev_print(LOGY,__LINE__)) {
					printf("recoverlen() is adding a missing final terminator at m=%d, n=%d.", m,n);
				}
				return(alpha_count);
			}
			else if (align2D[m][n] == '>') {
				return(alpha_count);
			}
			else {
				break;
			}
		} 
	}

	return(alpha_count);
}


/*****************************************************************************************/
void usage(char *usage_version)
{
	int i;

	printf("\nRunning maximal version %s, a program for micro-homology alignment (MHA).\n", usage_version);
	printf("\nUsage: ./maximal -[OPTIONS] sequence.txt\n");
	for (i=1; i<53; i++) {
		if (*stroptions[i]->description)
			printf("\n\t\t-%c\t%s",stroptions[i]->sym, stroptions[i]->description);
	}
	printf("\n\n");
	printf("\t\t\t* This option has a monotonically-increasing additive effect if repeated as an argument.\n");
	printf("\t\t\t  Multiple asterisks indicate the number of different option levels.\n");
	printf("\nExample usages: ./maximal -v         sequence_file.txt");
	printf("\n                ./maximal -KnO       sequence_file.txt");
	printf("\n                ./maximal -KnXXY 800 sequence_file.txt\n\n");
}

/*****************************************************************************************/
void warnhead(char l)
{ 
/*	printf("%c", 7); */	/* BELL CHARACTER */
	printf("\n * Notice (%c): ", l);
} 

/*************************************************************************************************/
short int pushdown(char pusharray[][MAXROW], int push_m, int push_n)
{
int i=0, j=0;
int h = options[1][17];		/* height slot */
char letr;
char blank = options[1][11];

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

/*** arrayA WILL USUALLY BE m2Dalig[] AND arrayB WILL BE align2D[] ********************/
int get2Dtucknum(char arrayA[][MAXROW], char arrayB[][MAXROW]) 
{
	int i=0, j=0, heightAB=0, heightA=0, heightB=0;
	int  width = options[1][32]; 
	char opt_R_rght = options[1][27];
	char letr; 
	int bottom[MAXROW] = {0};
	int    top[MAXROW] = {0};

	if (arrayA[MAXROW][width]!='\0') {
		for (j=width; (arrayA[MAXROW][j])!='\0' && j < MAXROW; j++) {
			;
		}
		width = j;
	}

	for (i=0; arrayA[i][0] != '\0' && i<MAXROW; i++) 
		;
	for (j=0; arrayB[j][0] != '\0' && j<MAXROW; j++) 
		;
	heightA = i;
	heightB = j;
	heightAB = i+j;

	/**** FILL BOTTOM BORDER ARRAY; BOTTOM BORDER IS FOR ARRAY ON TOP ****/
	bottom[width+1] = heightA;
	bottom[width  ] = heightA;
	j = width;
	while (isalpha(arrayA[heightA-1][j]) && j>=0) {
		bottom[j] = heightA;
		j--;
	}
	while (j>=0 && i>0) {
		if (isalpha(arrayA[i-1][j]) || isalpha(arrayA[i-1][j+1])) {
			bottom[j] = i;
			j--;
		}
		else {
			i--;
		}
	}

	/**** FILL TOP BORDER ARRAY; TOP BORDER IS FOR ARRAY ON BOTTOM  ****/
	top[0] = 0;
	for (j=0; j<=width; j++) {
		for (i=0; i<MAXROW; i++) {
			if (isalpha(letr=arrayB[i][j])) {
				top[j] = heightB - i;
				break;
			}
            else if (letr == '/' || letr == opt_R_rght) {
				top[  j] = heightB - i;
				top[++j] = heightB - i;
				break;
			}
			else if ((letr=arrayB[i+1][j-1]) == '/' || letr == opt_R_rght) {
				top[j] = heightB - (i+1);
				break;
			}
			else if (arrayB[i][j]=='>') {
				top[j] = heightB - i;
				j = width+1;
				break;
			}
		}
	}

	/**************************************************/
	if (opt_v.val>2) {	/* CODE DEVELOPMENT REPORTING */
		printf("\n ");	
		for (j=0; j <= width; j++)
			printf("%c", mha_base62(bottom[j]));	
			printf(" <-- bottom edge of top array\n ");
		for (j=0; j <= width; j++)
			printf("%c", mha_base62(top[j]));	
			printf(" <-- top edge of bottom array\n ");
	}
	/******************************************/

	bottom[0] = heightAB - bottom[0] - top[0];
	if (opt_v.val>2)
		printf("%c", mha_base62(bottom[0]));
	for (j=1; j <= width; j++) {
		if ((i=heightAB - bottom[j] - top[j]) < bottom[j-1]) {
			bottom[j] = i;
		}
		else
			bottom[j] = bottom[j-1];
		if (opt_v.val>2)
			printf("%c", mha_base62(bottom[j]));
	}
	if (opt_v.val>2)
		printf(" <-- MIN(heightAB - bottom edge - top edge)\n");
	i = heightA - bottom[j-1] + 1;
	if (opt_v.val>2)
		printf("\n 2Dtucknum = %d\n", i);
	if (i>0)
		return(i);
	else
		return(-1);
}

/*************************************************************************/
void print1D(void)
{
	int i, j, n;
	int lenseq = options[1][1];
	int blocks = count_wrap_blocks(lenseq, options[1][58]);
	int head_start = 0;			/* USE TO PASS RULER O-F-F-SET TO line_end() */
	char ch;

	mha_head(lenseq);

	for (j = 0; j < blocks; j++) {
		if (options[1][13] && opt_v.val) { /* IF DNA AND VERBOSITY */
			line_end(SLIPS, j+1, 0);	
   			for (n = j * options[1][58]; n < (j+1) * options[1][58] && tela[n].c!='>' && tela[n].c!='\0'; n++) {
				if ((ch=tela[n].t)=='R' || ch=='Y')
					printf("%c", ch);
				else
					printf("_");
			}
			printf("\n");
		}
		line_end(START, j+1, 0);	
   		for (n = j * options[1][58]; n < (j+1) * options[1][58] && tela[n].c!='>' && tela[n].c!='\0'; n++) {
			printf("%1c", tela[n].c);
		}
		if (tela[n].c == '>') {
			printf("%1c", tela[n].c);  /* PRINTS TERMINAL CHARACTER '>' */
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
			for (i = j * options[1][58]; i < (j+1) * options[1][58] && i<lenseq; i++) {
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
			for (i = j * options[1][58]; i < (j+1) * options[1][58] && i<lenseq; i++) {
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
			for (i = j * options[1][58]; i < (j+1) * options[1][58] && i<lenseq; i++) {
				printf("%c", tela[i].echoes);
			}
			if (j+1 == blocks)
				printf(" <==== (((( {  REVERB  } ))))\n");
			else
				printf("\n");
			/**********************************************/
			head_start = (j * options[1][58]) % 10;
			if (j+1 < blocks) {
				line_end(SLIPRULER, head_start, options[1][58]);
				printf("\n");
			}
			else 
				line_end(SLIPRULER, head_start, lenseq-options[1][58]*(blocks-1));
			/**********************************************/
		}   /* END OF opt_l PRINT MODULE */ 

	} /* END OF FOR j LOOP */
}

/*************************/
int span_allrk(int point)
{
	int product = tela[point].all_r * (tela[point].all_k);
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

#endif		/* !FILE_LOGY_SEEN */

/*************************************************************************************************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal */
/*************************************************************************************************************/
