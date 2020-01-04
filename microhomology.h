/* microhomology.h: header file for mha_v#.##.c */

#define MAXROW   1600       /* maximum input line size  */
#define WIDTH      72       /* BANDWIDTH: MAX WIDTH OF HEMIDIAGONAL OF PATHBOX; MAX TR UNIT SIZE */ 
#define MATCH       8       /* MATCH SCORE */
#define CYCMAX     60       /* MAGIC NUMBER; SEARCH MAGIC TO FIND OTHER STOPGAPS */
#define FRAME_ROWS 29       /* NUMBER OF AVAILABLE ROWS FOR STORING OVERLAPPING REPEAT FRAMES; MULT. OF 4 - EXTRA */
#define PISO        2       /* FLOOR FOR TRANSITION MATCHING ABOVE THIS k-MER SIZE */
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
	int Dtr;		/* DIAGONAL TANDEM REPEAT (DTR) SCORE BY POSITION */
	char t;			/* IUPAC TRANSITIONS IN DNA USUALLY (RY) IN "IMPERFECT" TANDEM REPEATS */
	/*************************************************************************************************/
	int  cyc_F[FRAME_ROWS];	/* cycling frames; count-off column positions per unit; 32 - 3 = 29 */
							/* one row/frame; row 0 is row # locator; FRAME_ROWS IS BASED ON MEM AL. */
	/*************************************************************************************************/
	int cyc_Lf;		/* Left-side overlapping TR; 0=lenseq, which is also stored in options[1][1]     */
	int cyc_Rt;		/* Right-side overlapping TR */
	char cyc_o;		/* x => cinched; o => untaken cyclelizable option; !,** => CHECK_TELA VIOLATIONS */
	/*************************************************************************************************/
} tela[MAXROW] = {0};
char align2D[MAXROW+1][MAXROW] = {{0}};

/***********************************************/
void signal_callback_handler(int signum) 
{
	printf("  )--- I caught signal %d before exiting (2=SIGINT, 11=SIGSEGV).\n\n",signum);
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

int assign_tela(int eL, int eM, int eN, int mode, int pointA, int pointB);
int assign_transit(int n);
int check_tela(int eM, int eN, short unsigned int dim);
int cyclelize_tela(int cpos, int delta, int npos, long int options[][62]);
void clear_2D_ar(char wipe_align2D[][MAXROW]);
void clear_right(char swipe_align2D[][MAXROW], long int croptions[][62]);
int col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown); 
unsigned int consensus_2D(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width);
int count_wrap_blocks(int lcl_width, int lcl_opt_w);	/* lcl_width IS WIDTH OF 2-D MHA ARRAY */ 
int get_1Dz(int x, int y, int ignoreCheck);
void line_end(int type, int c, long int lend_options[][62], int lcl_width);
char mha_base62(int num);
void mha_head(int lcl_width, long int lcl_options[][62]);
void mha_UPPERback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
void mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
void mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
int push_tela(int n2, int n1);
void print1D(long int options[][62]);
short unsigned int print_2Dseq(int print_lenseq2D, long int poptions[][62]);
void print_blockhead(int a, int b);	
void print_tela(int a, int b);
short int pushdown(char pusharray[][MAXROW], int push_m, int push_n, long int push_options[0][62]);
int score_DTHR(int kmer, int squeeze);
int span_rk(int point);
int update_tela(void);
void warnhead(char l); 

short unsigned int 	user_query(unsigned int pass_num);

int               	cinch_k(long int koptions[][62]);  
int 				recover_1D(char recovered_1D[MAXROW], long int rec_options[][62]);
int 				recoverlen(long int rec_options[][62]);
short unsigned int 	cleanseq(char *s);

int 				get2Dtucknum(char arrayA[][MAXROW], char arrayB[][MAXROW], long int options[][62]);
unsigned int       	nudgelize(char cyc_align2D[][MAXROW], long int cyc_options[][62]);
unsigned int       	cinch_d(long int doptions[][62], short unsigned int cinch_d_opt);
short unsigned int 	cinch_l(long int loptions[][62]);  

void		 		relax_2D(long int roptions[0][62]);
void 				mha_randomize1(char input_seq[MAXROW]);
void 				mha_randomize2(char input_seq[MAXROW], int rsize);
void 				print_base62_table(long int boptions[][62]);
void 				shuffle(int *array, int n);
void 				usage(char usage_version[], unsigned int FY_size);			/* FOR PRINTING UNIFORM USAGE INFORMATION */
short int			tucksense(char tuckarray[][MAXROW], long int tuck_options[0][62]);
char 				*nmer_prefix(int i);			/* CONVERTS INTEGER TO N-MER PREFIX WRITTEN NAME */

