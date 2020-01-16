/******************************************************************/
/* maximal's microhomology.h: header file since mha_v3.97.c       */
/******************************************************************/
/******************************************************************/
/******************************************************************/

#define MAXROW   2000       /* maximum input line size  */
#define WIDTH      72       /* BANDWIDTH: MAX WIDTH OF HEMIDIAGONAL OF PATHBOX; MAX TR UNIT SIZE */ 
#define CYCMAX     60       /* MAGIC NUMBER; SEARCH MAGIC TO FIND OTHER STOPGAPS */
#define FRAME_ROWS 25       /* NUMBER OF AVAILABLE ROWS FOR STORING OVERLAPPING REPEAT FRAMES; MULT. OF 4 - EXTRA */
#define START       0       /* FOR USE WITH line_end() */
#define END         1       /* FOR USE WITH line_end() */
#define RULER       2       /* FOR USE WITH line_end() */
#define SLIPS       3       /* FOR USE WITH line_end() */
#define PATHBOXHEAD 4       /* FOR USE WITH line_end() */
#define BLOCKHEAD   5       /* FOR USE WITH line_end() */
#define SLIPRULER   6       /* FOR USE WITH line_end() */
#define EXIT_GOOD	0		/* FOR STANDARD EXIT ERRORS */
#define EXIT_ERROR	1		/* FOR STANDARD EXIT ERRORS */
#define EXIT_EARLY	2		/* FOR STANDARD EXIT ERRORS */
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
} tela[MAXROW] = {0};
char align2D[MAXROW][MAXROW] = {{0}};
char pathbox[MAXROW][MAXROW] = {{0}};
char consensus[MAXROW] = {0};
char file_name[255];
FILE *fp_out;                           /* FILE FOR output.log */
char dev_notes[32] = "N/A";             /* STRING WRITTEN AS LAST FIELD IN OUTPUT FILE */

/***********************************************/
void signal_callback_handler(int signum) 
{
	printf("  )--- I caught signal %d before exiting (2=SIGINT, 11=SIGSEGV).\n\n",signum);
	fp_out = fopen("Surf_wavereport.mha", "a");		/* FOPEN RIGHT BEFORE WRITING TO MINIMIZE CHANCE OF CLOSING WITH OPEN FILES */
	fprintf(fp_out, "---->\tCanceled run for %s with signal=%d (2=SIGINT, 11=SIGSEGV). dev_notes: %s.\n", file_name, signum, dev_notes);
	fclose(fp_out);
	exit(signum);
}

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

long int options[4][62] = {
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,46, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0,10, 0, 0, 0,40,41, 0, 4, 0, 0, 0, 0, 0,32, 0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61},
/*s->cl>t->l->k->c->s->d->r        B  C  D     F     H  I     K  L  M     O  P  Q  R     T        W  X     Z  a--b--c---d---e---f---g---h   i       k   l   m   n   o   p       r   s   t   u   v   w   x       z*/
{48,49,50,51,52,53,54,55,56,57,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122}
}; /*                    |                                            |  |     |        |         |= Zero tick mark, default ' '= 32   */
                     /* 46 = "." (FULLSTOP) */					   /* (  ) */	   /* options[1][32] opt_W WILL STORE CURRENT CINCH-WIDTH */
                     /* 32 = " " (SPACE)    */ 					   /* L  R */	   /* options[1][0--9] WILL PERMANENTLY STORE 2-D WIDTH HISTORY */
/* options array rows */
/* options[0][n] = BIT MASK. 1 = OPTION CALLED AS ARG. W/ SOME EXCEPTIONS */
/* options[1][n] = OPTIONS VALUE */
/* options[2][n] = BASE 62 NUM VALUE, WILL USE AS opt_# INDEX */
/* options[3][n] = DECIMAL CHARACTER CODE OF LETTER SHOWN IN COMMENTS */
/* opt_M replaces long_homopolywrap = mwrap = 10 FOR WRAPPING LONG MONONUCLEOTIDE TRACTS */
								/* options[1][18] COUNTER OF INITIAL PASSES THROUGH MHA */
								/*  0 EQUALS ORIGINAL STRING			*/
								/*  1 EQUALS cleanseq PASS (1-D)    	*/
								/*  2 EQUALS cinch_t 2-D PASS MAIN()	*/
								/*  3 EQUALS cinch_l 2-D PASS			*/
								/*  4 EQUALS cinch_k 2-D PASS			*/
								/*  5 EQUALS nudgelize 2-D PASS			*/
								/*  6 EQUALS cinch_d  2-D PASS			*/
								/*  7 EQUALS relax_2D 2-D PASS			*/
								/*  8 EQUALS recovered 1-D from 2-D 	*/
								/* 10 EQUALS passQ score / 1000			*/
/* RESERVE options[1][26] (opt_Q) and options[1][27] (opt_R) for storing LEFT and RIGHT 'R'un delimiter characters */

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
short unsigned int	print_2Dseq(int print_lenseq2D);
void 				print_blockhead(int a, int b);	
short int 			pushdown(char pusharray[][MAXROW], int push_m, int push_n); 
int 				span_rk(int point);
void 				warnhead(char l); 
int               	cinch_k(void);  
int 				recover_1D(char recovered_1D[MAXROW]);
int 				recoverlen(void);
short unsigned int	cleanseq(char *s);
int 				get2Dtucknum(char arrayA[][MAXROW], char arrayB[][MAXROW]);
int					cinch_l(void);  
unsigned int       	nudgelize(void);
unsigned int       	cinch_d(short unsigned int cinch_d_opt);
void		 		relax_2D(void);
void 				mha_randomize1(char input_seq[MAXROW]);
void 				mha_randomize2(char input_seq[MAXROW], int rsize);
void 				print_base62_table(void);
void 				shuffle(int *array, int n);
void 				usage(char usage_version[], unsigned int FY_size);			/* FOR PRINTING UNIFORM USAGE INFORMATION */
char 				*nmer_prefix(int i);			/* CONVERTS INTEGER TO N-MER PREFIX WRITTEN NAME */
void 				free_2D(int **p2D, int lenseq);

