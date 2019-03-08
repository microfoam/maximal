/*******************************************************************************************************/
/***** The program "maximal" is a maximal homology alignment (MHA) program.                        *****/
/***** Designed and written by Dr. Albert J. Erives. 2017-2019. AGPL-3.0.                          *****/
/***** Code repository located at https://github.com/microfoam/maximal. Licensed under AGPL-3.0.   *****/
/***** This program renders a 1-D DNA sequence into a 2-D self-alignment to rescue micro-paralogy. *****/
/*******************************************************************************************************/

#include <ctype.h>  	/* isalpha() */
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 	/* system()  */
#include <string.h> 	/* strlen()  */
#include <time.h>		/* difftime() */

#define MAXROW	 1200		/* maximum input line size  */
#define WIDTH	  200		/* BANDWIDTH OF 2X MAX TR SIZE; 3x 100 bp MSR UNIT CEILING */ 
#define CYCMAX	   60		/* MAGIC NUMBER; SEARCH MAGIC TO FIND OTHER STOPGAPS */
#define START		0		/* FOR USE WITH line_end() */
#define END			1		/* FOR USE WITH line_end() */
#define RULER		2		/* FOR USE WITH line_end() */
#define SLIPS		3		/* FOR USE WITH line_end() */
#define PATHBOXHEAD	4		/* FOR USE WITH line_end() */
#define BLOCKHEAD	5		/* FOR USE WITH line_end() */
#define SLIPRULER	6		/* FOR USE WITH line_end() */
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })

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

int main(int argc, char *argv[])
{
	/* cinch_  t=tandem repeats (MAIN), l=long homopolymer runs, k=arbitrary intra-TR k-mers, d=de novo intra-TR repeats */
	short unsigned int cinch_l(char align2D_pass3[][MAXROW], long int loptions[][62]);  
	void               cinch_k(char align2D_pass4[][MAXROW], long int koptions[][62], int DTHR_lookie[WIDTH/2]);  
	unsigned int       cyclelize(char cyc_align2D[][MAXROW], char *cyc_Seq_name, long int cyc_options[][62]);
	short unsigned int cinch_s(char align2D_pass4[][MAXROW], long int soptions[][62]); 
	unsigned int       cinch_d(char align2D_pass7[][MAXROW], char *dptr_Seq_name, long int doptions[][62], short unsigned int cinch_d_opt);
	short int 			tucksense(char tuckarray[][MAXROW], long int tuck_options[0][62]);
	unsigned int 		consensus_2D(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width);
	short unsigned 		/* STRIPS NON-ALPHA AND CONVERTS TO UPPERCASE (EXCEPT n's IF DNA OR RNA */
			 int cleanseq(char *s, long int cloptions[][62], short unsigned int storebit);
	void 		 clear_right(char swipe_align2D[][MAXROW], long int croptions[][62]);
	int  		 count_wrap_blocks(int lcl_width, int lcl_opt_w);	/* lcl_width IS WIDTH OF 2-D MHA ARRAY */ 
	void 		 line_end(int type, int c, long int lend_options[][62], int lcl_width);
										/* PRINTS UNIFORM END OF LINE NUMBERS FOR 2-D ALIGNMENTS						*/
										/* line(START, m+1, options[][62], width) FOR LINE BEGINNINGS, W/ ROW NUMBER	*/
										/* line(END, c, options[][62], width) FOR LINE ENDS, W/ CHARACTER NUMBER		*/
										/* line(RULER, c, options[][62], width) FOR LINE ENDS, W/ CHARACTER NUMBER		*/
										/* line(PATHBOXHEAD, c, options[][62], width) FOR PATHBOX HEADER START			*/
										/* line(BLOCKHEAD, c, options[][62], width) FOR CONSENSUS LINE START			*/
	char 		 mha_base62(int num);	/* THIS IS A BASE 62 NUMBER CONVERTER, USES 0--9, A--Z, a--z		*/
										/*  ALLOWS REPRESENTATION OF TR # or k-MER SIZE USING SINGLE CHAR	*/
	void 		 mha_head(char lcl_seq_name, int lcl_width, long int lcl_options[][62]);
	void 		 mha_randomize1(char input_seq[MAXROW]);
	void 		 mha_randomize2(char input_seq[MAXROW], int rsize);
	void 	     mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
	void 		 mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
	char 		*nmer_prefix(int i);			/* CONVERTS INTEGER TO N-MER PREFIX WRITTEN NAME */
	short 
	unsigned int print_2Dseq(char align2D_print[][MAXROW], int print_lenseq2D, char *printptr_Seq_name, long int poptions[][62]);
	void 		 print_blockhead(int a, int b);	
	void 		 print_base62_table(long int boptions[][62]);
	void 		 recover_1D(char recovered_1D[MAXROW], char rec_align2D[][MAXROW], long int rec_options[][62]);
	void 		 relax_2D(char align2D_pass8[][MAXROW], long int roptions[0][62]);
	void 		 shuffle(int *array, int n);
	void 		 usage(char usage_version[], unsigned int FY_size);			/* FOR PRINTING UNIFORM USAGE INFORMATION */
	short 
	unsigned int user_query(unsigned int pass_num);
	void 		 warnhead(char l); 
	int get2Dtucknum(char arrayA[][MAXROW], char arrayB[][MAXROW], long int options[][62]);

	char version[] = "3.46";				/* current version number */
	char dev_notes[32] ={0};				/* STRING WRITTEN AS LAST FIELD IN OUTPUT FILE */
	int match      =  8;
	int mismatch   = -1;
	int transition =  match/2;			
	float fractransit = 0.08;				/* USED TO SET NUMBER OF ADDITIONAL ALLOWED TRANSITIONS FOR GIVEN k-MER */
											/* FOR 0.096, ALL k < 6 ALLOW ONLY 1 TRANSITION, 3 TRANSITIONS STARTING AT k > 15 */
											/* FOR 0.08, ALL k < 7 ALLOW ONLY 1 TRANSITION */
	int numtransit;
	int DTHR       = 101;					/* 101%. RESET IF opt_x, DIAGONAL THRESHOLD SCORES */ 
	int DTHR_lookup[1+WIDTH/2] = {0};		/* ARRAY OF PRE-COMPUTED k-DEPENDENT DTHR VALUES */
	short unsigned int continue_flag=0, imperfect_TR=0, Aimperfect_TR=0, nuctransit=0, seqtype=0;
	int alt_Did, alt_Dtr, alt_k, alt_m;		/* ALTERNATE k-MER COMPARISONS */
	unsigned int FY_size = 100;				/* SIZE OF FISHER-YATES RANDOMIZED STRING */
	char blank = '.';						/* DEFAULT BLANK CHARACTER FOR 2-D MHA. FULLSTOP = 46 */
											/* NEEDS TO BE SET IN options[][62] array */
	char transtype = 'X';					/* transtype WILL STORE EITHER 'R' OR 'Y' */
	char letr_unit[] = "alpha-characters";	/* DEFAULT STRING ALPHABET, LATER RESET TO "bp" FOR DNA AND "nt" FOR RNA */
	int cyc_runs=0, d_runs=0, r_runs=0, homopoly_flag=0, homopolyend_flag=0, i=0, j=0, k, m=0, mstop=0, n=0, oldbad=0; 
	int lenseq=MAXROW, r=0, scooch=0, TRcheck=0, z=0, sloc_z=0; 
	int h=0, l=0, overslip = 0; 
	float ratio1 = 1;						/* WIDTH CINCH RATIO (W.C.R.) post cinch-d, pre relax-2D 	*/
	float ratio2 = 1;						/* WIDTH CINCH RATIO (W.C.R.) post relax-2D 				*/

	FILE *fp_out;							/* FILE FOR output.log */
	FILE *fp_cons;							/* FILE FOR CONSENSUS STRING Surf_barrels.log */
  	FILE *fp_msa;							/* FILE FOR MHA MSA "TUBES.mha" */
	FILE *fp_tricksy;						/* DEV. FILE FOR IMPERFECT 2-D ALIGNED STRINGS waves/foam_and_chowder.mha */

	/* Variables for time calls */
	time_t lcl_time = time(NULL);

	long int cit_new2D_width = 0;
	long int a2D_n = 0;						/* NUMBER INDEX OF n FOR a2D_n */
	int delta_to_bad1Dn = 0;				/* TO KEEP TRACK OF TIME SINCE LAST SLIP */
	int bad_1Dn = 0;						/* POSITION OF BAD SLIP IN a2D_n COORDINATE */
	short unsigned int frst_badslip = 0;	/* BIT FLAG FOR FIRST BAD SLIP DECLARATION */
	short unsigned int msa = 0;				/* BIT FLAG FOR MSA CO-INUPUT */

	char Seq_a[] = "abcdefdefdefdefdefffffffghijkhhhijkhijhijklmnopnopxopxopxopopohyeaahhaahhopopopopxyzyzy"; 
	char Seq_b[] = "bwaahahatzatzxyxybwaahahatzatzxyxybwaahahatzatzxyxyyyopqqqqqqqqqqqqqqrstututuvwvwxyz"; 
	char Seq_c[] = "gatcgatcgatcgatcxoxohmmmatathmmmatathmmmatathmmmatathmmmatathmmmatathmmmatatopqrstututuvwvwxyzyyyzyzyz"; 
	char Seq_d[] = "databpttchmmmatathmmmatathmmmatathmmmmatathmomopakhmomopakhmomopakhmomopakabcdefghijklmnnnnnnnnnnnnnnnnnnnnnnopqrstuvwxyz"; 
	char Seq_e[] = "exxxamplexxxamplexxxammplhamhammmmmmmstringingyzyzzzzzzzzzzzzzzzzzzzzzzzzzzzzgattacagattacaaTGTTTTTGGGAGGGAGTGCATCATzzzzzzzzzzzzzzzzzzzzzzz"; 
	char Seq_f[] = "Fuuugh!Fuuggh!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxYZgatcgatcgatcgatcgatc"
				   "hmmmatatxyxyhmmmmatatxyxyhmmmatatxyxyhmmmatatxyxyopqrRRRRRRRRRstututu12345678./>|*"; 
	char Seq_g[] = "A-Goo Goo Goo, A-Ga Ga Ga, is all I wannnnnnt to say to youuuu!";
	char Seq_h[] = "TTGTGTTGCTGTGTGAGTGAGGTCGGATCGCTATGGGGGACCGGCTGACGAGCTANNNNNNRGATTACA";

	char Seq_i[MAXROW]={0}; 	/* INPUT SEQUENCE INSTEAD OF BUILT-IN EXAMPLES */
	char Seq_r[MAXROW]={0}; 	/* RANDOMIZED SEQUENCE */
	char SeqRY[MAXROW]={0};	/* TO STORE MICROPARALOGY TYPE */
	char *Seq = Seq_a;			/* POINTER TO EXAMPLE SEQUENCE, DEFAULT IS Seq_a, PRE-REFORMATTED */
	char Seq_name = 'a';		/* Example choice name */
	char Seq_head[100] = {0};	/* FASTA HEADER */
	char *ptr_Seq_name;

	ptr_Seq_name = &Seq_name;
	FILE *file_ptr;
	char file_name[255];
	char cycle[WIDTH];		/* THIS ARRAY HOLDS THE CYCLIC PATTERN OF TRs W/ >2 UNITS */
 	char m2Dalig[MAXROW+1][MAXROW] = {{0}};			

	/* IS THERE A FILE NAME ARGUMENT? */
	for (i = 1; i < argc; i++) {
		if (*argv[i] != '-') {	
			if (strcmp(argv[i],"TUBES.mha")) {		/* strcmp EVALUATES TO 0 ONLY IF STRINGS ARE THE SAME */
				/* USING j BELOW AS COUNTER TO NUMBER OF TIMES USER SPECFIES DIFFERENT SEQUENCES */

				if ( (file_ptr = fopen(argv[i], "r") ) == NULL) {
					printf("\n %2d. Error opening file '%s'.", ++j, argv[i]);
				}
				else {
					fseek(file_ptr, 0, SEEK_END);
					lenseq = ftell(file_ptr);
					if (lenseq > MAXROW) {		/* LAST ROW OF array2D WILL STORE CONSENSUS, SO NEED TO KEEP CLEAR */
						printf("\n %2d. Sequence (length %d) from file '%s' exceeds MAXROW limit (%d) by %d.", ++j, lenseq, argv[i], MAXROW, lenseq-MAXROW+1);
						printf("\n %2d. Exiting now. For help enter './maximal -H'. \n\n", ++j);
						exit(1);
					}
					else fseek(file_ptr, 0, SEEK_SET);
				 
					fscanf(file_ptr, "%[^!]c", Seq_i);
					fclose (file_ptr);
					strcpy(file_name, argv[i]);
					Seq = Seq_i;
					Seq_name = 'i';
					printf("\n %2d. Detecting sequence from file i: '%s'.", ++j, file_name);
	
					/* CHECK FOR FASTA HEADER AND SAVE IN Seq_head, THEN MASK IN Seq */
					if (Seq_i[0] == '>') {
						for (h = 0; Seq_i[h+1] != '\n' && Seq_i[h+1] != '\r' && h < 100; h++) {
							Seq_head[h] = Seq_i[h+1];
						}

						/* SCOOCH STRING INTO FASTA HEADER SPACE (ERASING IT) */
						i = (int) strlen(Seq_head) + 2;
						for (h = 0; Seq_i[h+i] != '\0'; h++)
							Seq_i[h] = Seq_i[h+i];
	
						Seq_i[h] = '\0';	
					}
					else if (Seq_name == 'i')
						strcpy(Seq_head, "input sequence");
				}
			}
			else {
				if ( (file_ptr = fopen(argv[i], "r") ) == NULL) {
					printf("\n*%2d. Error opening supporting file '%s'.", j, argv[i]);
				}
				else if (msa == 0){
					msa = 1;		/* USING THIS SLOT TO STORE BIT VALUE INDICATING TUBES.mha CO-INPUT */
					char ch=blank;

					printf("\n*%2d. Acknowledging requested use of supporting file 'TUBES.mha'.", j);

					for (m=0; m<MAXROW; m++) {
					    for (n=0; n<MAXROW; n++) {
					        fscanf(file_ptr, "%c", &ch);
					        m2Dalig[m][n] = ch;
							if (ch=='\r' && m2Dalig[m][n+1]=='\n') {		/* MS-DOS LINE ENDINGS */
								m2Dalig[m][n  ] = '\0';
								m2Dalig[m][n+1] = '\0';
								break;
							}
							else if (ch=='\r' || ch=='\n') {				/* macOS OR UNIX LINE ENDINGS */
								m2Dalig[m][n] = '\0';
								break;
							}
					    }
					}
					fclose (file_ptr);

					if (msa) {	
						printf("\n\nRead file 'TUBES.mha' into m2Dalig array as follows:\n");
						for (m=0; m2Dalig[m][0]!='\0'; m++)	
							printf("%s\n", m2Dalig[m]);
/*						for (n=0; m2Dalig[m-1][n]!='>'; n++)
							m2Dalig[m][n] = blank;
*/
					}
				}
			}

		}	/* END OF IF ARGV[I] != '-' */
	}	/* END OF FOR i = 1, i < argc, i++ */

	if (Seq_name != 'i') {	
		strcpy(file_name, "example sequence");
		strcpy(Seq_head, "example sequence");
	}

	char recovered[MAXROW] = {0};
	short unsigned int go_flag=0, cycle_flag=0;			/* USE THIS TYPE FOR TRUE BIT FLAG VARIABLES */
	unsigned int recovery_flag = 0;
	int cyclelen=0;
	int relax_length=0;			/* FOR USE WITH relax_2D CALL */
	int head_start = 0;			/* USE TO PASS RULER OFF-SET TO line_end() */
	int intraTR_reps_tot = 0; 	/* STORES INITIAL RETURN VALUE FROM cinch-d() */
	int intraTR_reps = 0;	 	/* STORES CURRENT RETURN VALUE FROM cinch-d() */
	int Did = 0;				/* Counter for identity (id) diagonal */
	int Dtr = 0;				/* Counter for tandem repeat (tr) diagonal */
	int Atr = 0;				/* Counter for additional repeats on the same diagonal */
	int row = 0;				/* Counter for row number in align2D box */
	int c = 0;					/* Counter of Seq alpha characters printed */
	int unislip = 0;			/* Counter of unique slips */
	int badslipspan = 0;			/* Counter of length of a conflicting prior slip starting from the second unit */ 
	int recslips= 0;			/* Counter of recent slips in region of first TR unit, derived from sliploc */

	int slips[WIDTH+1] = {0};	/* Array of counters for unique slips of WIDTH x	*/
								/*  1st element init. to 0, & by default rest to 0	*/
	unsigned int sliploc[MAXROW] = {0};	 		/* Array of counters for slips at all positions 			*/
	unsigned int sliploc_nmer[MAXROW] = {0};		/* Array of counters for slip type at all slip positions 	*/
	char sliploc_echoes[MAXROW] = {0};			/* Array of nucleotide reverb 								*/
	char pathbox[MAXROW+1][MAXROW] = {{0}};
	char align2D[MAXROW+1][MAXROW] = {{0}};			/* main() 2-D ARRAY MATRIX	*/
	char scratch[MAXROW+1][MAXROW];					/* STILL USED CODE DEVELOPMENT BUT CURRENT TESTS HARDLY REVERT */
	int passQ[10];									/* Array to store width cinch ratio through passes */	

	int o;						/* o IS CASE OPTION VARIABLE FOR SETTING options ARRAY (opt_ VARs) */ 

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
								/*  5 EQUALS cyclelize 2-D PASS			*/
								/*  6 EQUALS cinch_s  2-D PASS			*/
								/*  7 EQUALS cinch_d  2-D PASS			*/
								/*  8 EQUALS relax_2D 2-D PASS			*/
								/*  9 EQUALS recovered 1-D from 2-D 	*/
								/* 10 EQUALS passQ score / 1000			*/
/* RESERVE options[1][26] (opt_Q) and options[1][27] (opt_R) for storing LEFT and RIGHT 'R'un delimiter characters */
	char optR = options[1][27];
	options[1][58] = 80;	/* opt_w IS MAXIMUM WIDTH OF 1D/2D-PRINTED BLOCK, SCREEN WRAP LENGTH */
	int blocks;				/* Number of blocks for 1D output print */

	if (argc == 1) {
	 	system("clear"); 
		usage(version, FY_size);
		exit(1);
	}

	/**************************************/
	/* SET OPTIONS FROM ARGUMENTS  ********/
	while (--argc > 0 && (*++argv)[0] == '-') {
		while ((o = *++argv[0])) {
			switch (o) {
			case 'a':				/* OPTION TO USE EXAMPLE a */
					Seq = Seq_a;
					Seq_name = 'a';
					printf("\n %2d. You selected example sequence '%c'.", ++j, Seq_name);
					break;
			case 'b':				/* OPTION TO USE EXAMPLE b */
					Seq = Seq_b;
					Seq_name = 'b';
					printf("\n %2d. You selected example sequence '%c'.", ++j, Seq_name);
					break;
			case 'c':				/* OPTION TO USE EXAMPLE c */
					Seq = Seq_c;
					Seq_name = 'c';
					printf("\n %2d. You selected example sequence '%c'.", ++j, Seq_name);
					break;
			case 'd':				/* OPTION TO USE EXAMPLE d */
					Seq = Seq_d;
					Seq_name = 'd';
					printf("\n %2d. You selected example sequence '%c'.", ++j, Seq_name);
					break;
			case 'e':				/* OPTION TO USE EXAMPLE e */
					Seq = Seq_e;
					Seq_name = 'e';
					printf("\n %2d. You selected example sequence '%c'.", ++j, Seq_name);
					break;
			case 'f':				/* OPTION TO USE EXAMPLE f */
					Seq = Seq_f;
					Seq_name = 'f';
					printf("\n %2d. You selected example sequence '%c'.", ++j, Seq_name);
					break;
			case 'g':				/* OPTION TO USE EXAMPLE g */
					Seq = Seq_g;
					Seq_name = 'g';
					printf("\n %2d. You selected example sequence '%c'.", ++j, Seq_name);
					break;
			case 'h':				/* OPTION TO USE EXAMPLE h */
					Seq = Seq_h;
					Seq_name = 'h';
					printf("\n %2d. You selected example sequence '%c'.", ++j, Seq_name);
					break;
			case 'i':				/* OPTION TO USE Seq_i FROM INPUT FILE */
					if (strlen(Seq_i) == 0) {
						warnhead('i');
						printf ("No valid input file given.");
					}
					break;
			case 'k':						/* OPTION TO SHOW k-MER COUNTS */
					options[0][46] = 1;		/* opt_k ON   */
					break;
			case 'l':						/* OPTION TO SHOW SLIP LOCATIONS IN 1D SEQUENCE */
					options[0][47] = 1;		/* opt_l ON   */
					break;
			case 'n':						/* OPTION TO NOT DO RELAX-2D PASS */
					options[0][49] = 1;		/* opt_n ON   */
					/* options[1][48-49] RESERVED FOR STORING CONFLICT ROWS AND COLS */
					break;
			case 'o':						/* OPTION TO PRINT ORIGINAL STRING UNFORMATTED */
					options[0][50] = 1;		/* opt_o ON   */
					break;
			case 'p':						/* OPTION TO SHOW PARAMETERS */
					options[0][51] = 1;		/* opt_p ON   */
					break;
			case 'r':						/* OPTION TO SHOW ROW NUMBERING */
					options[0][53] = 1;		/* opt_r ON   */
					break;
			case 's':						/* OPTION TO SILENCE WRITE TO NORMAL OUTPUT FILE */
					options[0][54] = 1;		/* opt_s ON   */
					break;
			case 'u':						/* OPTION TO PRINT UNWRAPPED WHERE opt_w EQUALS lenseq */
					options[0][56] = 1;		/* opt_u ON 	*/
					++options[1][56];		/* ++opt_u VAL	*/
					break;
			case 'v':						/* OPTION FOR VERBOSITY */
					options[0][57] = 1;		/* opt_v ON 	*/
					++options[1][57];		/* opt_v INCREMENTED */

					/*      opt_F            opt_L    */
					options[0][15]  =options[0][21] = 1;
					++options[0][18];		/* opt_I */
					++options[0][27];		/* opt_R */
					options[0][47] =options[0][50] =options[0][51] =options[0][53] = 1;
					/*      opt_l           opt_o           opt_p           opt_r    */
					break;
			case 'x':						/* OPTION TO SQUEEZE DTHR VALUES BY 1 FOR k > 2 */
					options[0][59] = 1;		/* opt_x ON 		*/
					++options[1][59];		/* increment opt_x, NOT NECESSARILY IN USE IN CURRENT VERSION	*/
					break;
			case 'z':						/* OPTION FOR ZERO MISMATCH SCORE */
					options[0][61] = 1;		/* opt_z ON 	*/
					mismatch = mismatch - mismatch;
					break;
			case 'B':						/* USE SPACE FOR BLANK CHARACTER IN MHA's */
					++options[0][11];		/* opt_B ON & INCREMENTED */
					blank=options[1][11]=32;/* blank EQUALS ' ' (SPACE) */
					break;
			case 'C':						/* SHOW BASE 62 CODE */
					print_base62_table(options);
					exit(1);
					break;
			case 'D':
				   	options[0][13] = 1;		/* opt_D: SHOW DTHR_lookup values and early exit */
					/* optiosn[1][13] IS RESERVED FOR STORING SEQUENCE TYPE (DNA, RNA, PROTEIN, BABYLONIAN, etc.) */
					break;
			case 'F':						/* OPTION TO USE BLANK FILL CHAR W/ SCRIMMAGELINE */
					options[0][15] = 1;		/* opt_F ON: USE POST-TERMINATION FILLER	*/
					break;
			case 'H':						/* OPTION TO SHOW HELP */
					options[0][17] = 1;		/* opt_H ON 	*/
				/*  options[1][17] RESERVED FOR 2D-ALIGNMENT HEIGHT */	
					break;
			case 'I':						/* OPTION TO SHOW iNITIAL PASSES */
					++options[0][18];		/* opt_I = 1 is ON; opt_I > 1 is OFF (lets you toggle off after option v) */
				/*  options[1][18]            RESERVED FOR COUNTING iNITIAL PASSES. */
					break;
			case 'K':						/* OPTION TO SHOW CONSENSUS ROW */
					options[0][20] = 1;		/* opt_K ON */
					++options[1][20];		/* INCREMENT opt_K */
					break;
			case 'L':						/* OPTION TO SHOW LETTER #'s AT END OF LINES */
					options[0][21] = 1;		/* opt_L ON */
					break;
			case 'M':						/* OPTION TO DOUBLE LONG HOMOMONO WRAP */
					options[0][22] = 1;									/* opt_M ON		*/
					options[1][22] = options[1][22] + options[1][22];   /* DOUBLE opt_M  mwrap */
					break;
			case 'O':						/* OPTION TO OUTPUT CONSENSUS FILE */
					++options[0][24];		/* opt_O ON 	*/
					options[0][49] = 1;		/* opt_n NO RELAX 2-D */
					break;
			case 'P':						/* OPTION TO PRINT PATH BOX */
					options[0][25] = 1;		/* opt_P ON 	*/
					break;
			case 'R':						/* OPTION TO PRINT RECOVERED 1-D SEQUENCE FROM LAST 2-D */
					++options[0][27];		/* opt_R = 1 is ON; opt_R > 1 is OFF (lets you toggle off after option v) */
					break;
			case 'X':						/* OPTION TO SCRAMBLE SEQUENCE           */
					options[0][33] = 1;		/* opt_X ON 							 */
					++options[1][33];		/* INCREMENT opt_X to desired level 	 */
					break;					/*  X = rand() cheese, XX = FISHER-YATES */
			default:
					printf("maximal: Illegal option %c\n", o);
					argc = 0;
					usage(version, FY_size);
					exit(1);
					break;
			} /* END SWITCH o			*/
		} /* END WHILE o			*/
	} /* END WHILE argc argv	*/

	/*******************************************************/
	/* BEGIN OUTPUT OF MAXIMAL PROGRAM  ********************/
	if (options[0][17]) {		/* opt_H SHOW USAGE + EXAMPLES */
	 	system("clear"); 
		usage(version, FY_size);
		printf(  "Example a:\n\"%s\"", Seq_a);
		printf("\nExample b:\n\"%s\"", Seq_b);
		printf("\nExample c:\n\"%s\"", Seq_c);
		printf("\nExample d:\n\"%s\"", Seq_d);
		printf("\nExample e:\n\"%s\"", Seq_e);
		printf("\nExample f:\n\"%s\"", Seq_f);
		printf("\nExample g:\n\"%s\"", Seq_g);
		printf("\nExample h:\n\"%s\"\n\n", Seq_h);

		exit(1);
	}

	/* IF OPTION I TOGGLED BACK OFF FROM OPTION v (VERBOSE), then set opt_I back to zero for screen reporting purposes */
	if (options[0][18] > 1) {
		options[0][18] = 0;
	}

	/* OPTION TO APPEND TO GROWING 2-D MSA FILE */
	if (options[0][24] >= 2) {
		options[0][24] = 2;
		options[0][49] = 1;		/* opt_n NO RELAX 2-D */
	}

	/* IF OPTION R TOGGLED BACK OFF FROM OPTION v (VERBOSE), then set opt_R back to zero for screen reporting purposes */
	if (options[0][27] % 2 == 0) {		
		options[0][27] = 0;
	}

	if (j > 1) {
		printf("\n");
		warnhead('S');
		printf("Many sequences specfied. Using last sequence '%c'.\n", Seq_name);
	}
	else if (j == 0) {
		printf("\n");
		warnhead('S');
		printf("No sequences specfied. Using example sequence '%c'.\n", Seq_name);
	}
	printf("\n");
		
	mha_head(Seq_name, options[1][58], options);
	printf("maximal homology alignment (MHA) -%c", Seq_name);
	for (i = 10; i < 36; i++) {			/* UPPER-CASE LETTER OPTIONS 		*/
		if (options[0][i] > 0)
			printf("%c", (char) options[3][i]);
	}
	for (i = 45; i < 63; i++) {			/* LOWER-CASE LETTER OPTIONS > 'i' 	*/
		if (options[0][i] == 1)
			printf("%c", (char) options[3][i]);
	}
	if (options[0][11] > 1) {			/* opt_B */
		printf(" -B%ld", options[0][11]);
	}
	printf(" -M%ld", options[1][22]);	/* opt_M */
	if (options[0][33]) {				/* opt_X */
		if (options[1][33])
			printf(" -X [USE RANDOMIZED SEQUENCE]");
		else 
			printf(" -XX [USE FISHER-YATES RANDOMIZED SEQUENCE]");
	}
	if (options[1][56] > 1) {			/* opt_u */
		printf(" -u%ld", options[1][56]);
	}
	printf(" (version %s)\n", version);
	printf("%s", ctime(&lcl_time));

	if (options[0][51]) {		/* opt_p SHOW RUN PARAMETERS */
		printf(" Match: %d, Transition: %d, Mismatch: %d, Threshold: %d%%, Bandwidth: %d, Default block width: %ld",  
				 match,     transition,     mismatch,     DTHR,            WIDTH,         options[1][58]);
		if (options[1][56] == 1) {	/* opt_u */
			printf("*\n * Block width reset for one block of (unwrapped) sequence input.\n");
		}
		else if (options[1][56] > 1) {		/* opt_u */
			printf("*\n * Block width reset for %ld blocks of sequence input.\n", options[1][56]);
		}
		else
			printf("\n");
	}

	/**********************************************/
	/* INITIALIZATIONS/DECLARATIONS AFTER OPTIONS */
	lenseq = strlen(Seq);
	if (strlen(Seq) > MAXROW) {		/* LAST ROW OF array2D WILL STORE CONSENSUS, SO NEED TO KEEP CLEAR */
		warnhead('M');
		printf("Sequence %c (length %d) exceeds MAXROW size limit (%d) by %d.\n\n", Seq_name, lenseq, MAXROW, lenseq-MAXROW+1);
		exit(1);
	}

	++options[1][18];	/* INCREMENT opt_I COUNTER (pass_num) */
	if (Seq_name == 'i')
		lenseq = lenseq - 1;
	options[1][0] = options[1][32] = lenseq;	/* ASSIGN CINCH-WIDTH TO HISTORY [0--9] AND CURRENT [32]	*/
	if (options[1][56] >= 1) {					/* WILL CAUSE OUTPUT TO NOT BE WRAPPED (opt_u EQUALS 1),	*/ 
		options[1][58] = lenseq/options[1][56];	/*  OR WRAPPED INTO opt_u NUMBER OF BLOCKS.				    */
	}

	/* OPTION TO PRINT ORIGINAL STRING IN BLOCKS *******************************/
	if (options[0][50]) {	/* opt_o */
		blocks = count_wrap_blocks(lenseq, options[1][58]);

		mha_head(Seq_name, lenseq, options);
		if (strlen(Seq_head) > 0)
			printf(">%s\n", Seq_head);
		printf("\"");
		for (j = 0; j < blocks; j++) {
			for(n = j * options[1][58]; (n < (j+1) * options[1][58]) && (Seq[n] != '\0'); n++) {
				if (Seq[n] != 10 && Seq[n] != 13 && Seq[n] != EOF)
					printf("%c", Seq[n]);
				else
					printf("*");
			}
			if (n < lenseq)
				printf("\n ");
			else
				printf("\"\n");
		}
	}

	/***************************************************************************/
	seqtype = cleanseq(Seq, options, 1);	/* opt_D: use it to store type; 1 = DNA, 2 = RNA, 3 = PROTEIN, 0 = OTHER */

	if (seqtype == 1)
		strcpy(letr_unit, "bp");	/* BASE PAIRS */
	else if (seqtype == 2)
		strcpy(letr_unit, "nt");	/* NUCLEOTIDES */
	else if (seqtype == 3)
		strcpy(letr_unit, "aa");	/* AMINO ACIDS */
	else if (seqtype == 0)
		strcpy(letr_unit, "char");	/* OTHER */

	if (seqtype == 1) {			/* SET BIT VAR nuctransit IF DNA */
		nuctransit = 1;
	}

	if (options[0][33]) {		/* USE RANDOMIZED SEQUENCE */
		strcpy(Seq_r, Seq);
		srand(time(0));

		if (options[1][33] == 1)
			mha_randomize1(Seq_r);
		else 
			mha_randomize2(Seq_r, FY_size);

		printf("\nRandomized sequence: \"");
		for (i = 0; Seq_r[i] != '\0'; i++)
			printf("%c", Seq_r[i]);
		printf("\"\n");
		Seq = Seq_r;
	}

	++options[1][18];		/* INCREMENT ++pass_num */
	lenseq = strlen(Seq);
	options[1][1] = options[1][32] = lenseq;	/* ASSIGN CINCH-WIDTH TO HISTORY [0--9] AND CURRENT [32]	*/

	passQ[0] = 1000*options[1][0]/lenseq;
	passQ[1] = 1000;

	Seq[lenseq]   = '>';
	Seq[lenseq+1] = '\0';   

	for (i = 0; i <= lenseq; i++) {
		sliploc[i] = 0;		/* Initialize slip location array */
		sliploc_nmer[i] = 0;
		sliploc_echoes[i] = blank;
	}

	/************************************************************************************/
	/* INITIALIZE PATHBOX FOR SELF-MHA  *************************************************/

	homopoly_flag = 1;								/* FOR LONG HOMOPOLYMERIC RUN CASE **/
	homopolyend_flag = 0;

	for (n = 0; Seq[n] != '\0'; n++)				/* SET IDENTITY LINE ****************/
		pathbox[n][n] = match;

	for (n = 0; Seq[n] != '\0'; n++) {
		if  (n < WIDTH) {							/* SET VALUES FOR 1ST WIDTH COLS ****/
				for (m = 0; m < n; m++){
					if (seqtype && Seq[n] == 'n')			 
						pathbox[m][n] = mismatch;   /* SPECIAL TREATMENT FOR 'n' IN DNA**/
					else if (Seq[n] == Seq[m])
						pathbox[m][n] = match;		/* MATCH ****************************/
					else if (nuctransit) {			/* IF DNA AND CHECKING FOR TRANSITIONS */
						if      (Seq[n] == 'A' && Seq[m] == 'G')
							pathbox[m][n] = transition;  
						else if (Seq[n] == 'G' && Seq[m] == 'A')
							pathbox[m][n] = transition;   
						else if (Seq[n] == 'C' && Seq[m] == 'T')
							pathbox[m][n] = transition;   
						else if (Seq[n] == 'T' && Seq[m] == 'C')
							pathbox[m][n] = transition;   
						else
							pathbox[m][n] = mismatch;   /* MISMATCH IF NO TRANSITION ****/
					}	
					else 
						pathbox[m][n] = mismatch;   /* MISMATCH *************************/
				}
		}
		else {										/* SET VALUES FOR REST OF COLUMNS ***/
				for (m = n-WIDTH; m < n; m++){		/*  WITHIN BAND WIDTH ***************/
					if (seqtype && Seq[n] == 'n')			  
						pathbox[m][n] = mismatch;   /* DUE TO NUCLEOTIDE AMBIGUITY ******/
					else if (Seq[n] == Seq[m])
						pathbox[m][n] = match;		/* MATCH ****************************/
					else if (nuctransit) {			/* IF DNA AND CHECKING FOR TRANSITIONS */
						if      (Seq[n] == 'A' && Seq[m] == 'G')
							pathbox[m][n] = transition;  
						else if (Seq[n] == 'G' && Seq[m] == 'A')
							pathbox[m][n] = transition;   
						else if (Seq[n] == 'C' && Seq[m] == 'T')
							pathbox[m][n] = transition;   
						else if (Seq[n] == 'T' && Seq[m] == 'C')
							pathbox[m][n] = transition;   
						else
							pathbox[m][n] = mismatch;   /* MISMATCH IF NO TRANSITION ****/
					}	
					else 
						pathbox[m][n] = mismatch;   /* MISMATCH *************************/
				}
		}

		if (Seq[n] != Seq[n+1])	{					/* HERE CHECK FOR HOMOPOLYMERIC RUN */
			homopoly_flag = 1;						/*  OF LENGTH = (defined) WIDTH	 */
			if (homopolyend_flag == 1) {			/* TREAT LAST COL OF HOMOPOLY. RUN **/
				homopolyend_flag = 0;
				for (i = n - WIDTH + 1; i < n; i++)
					pathbox[i][n] = mismatch;
			}
		}
		else
			homopoly_flag++;

		if (homopoly_flag > WIDTH) {
			for (j = n-WIDTH+1; j < n+1 ; j++) {
				for (i = n-WIDTH+1; i < j; i++) {
					pathbox[i][j] = mismatch;
				}
			}
			homopolyend_flag = 1;
		} 
	} /**********************************************************************************/

	/**********************************************/
	/* PRINT VALUES OF PATH BOX IF OPTION SET *****/

	if (options[0][25] == 1) {	/* opt_P */
		blocks = count_wrap_blocks(lenseq, options[1][58]);

		printf("\nPATHBOX PASS (length = width = %d)\n\n", lenseq);
		for (j = 0; j < blocks; j++) {
			if (blocks != 1)
				print_blockhead(j+1, blocks);
			line_end(PATHBOXHEAD, 9, options, 9);	
			for(n = j * options[1][58]; (n < (j+1) * options[1][58]) && (Seq[n] != '\0') && Seq[n] != '>'; n++) 
				printf("%2c", Seq[n]);
			printf("\n");
			for(m = j * options[1][58]; (m < (j+1) * options[1][58]) && (Seq[m] != '\0') && Seq[m] != '>'; m++) {
				printf("%4d. %c ", m+1, Seq[m]);
					for (n = j * options[1][58]; (n < (j+1) * options[1][58]) && (Seq[n] != '\0') && Seq[n] != '>'; n++) {
						if (m > n)
							printf("%c%c", blank, blank);
						else
							printf("%2d", pathbox[m][n]);
				}
				printf("\n");
	   		 }
		}
	} /* END OF OPTION TO PRINT PATHBOX */

	/*********************************************/
	/*  USE PATHBOX TO BUILD FIRST 2-D ALIGNMENT */
	/*			  [cinch_t BEGINS]		         */

	a2D_n = row = 0;
	cit_new2D_width = lenseq;				/* START AT 1D LENGTH AND CONTRACT DURING CINCH-T */

	align2D[row][a2D_n++] = Seq[0];			/* FOR n=0, ENTER VALUE AT IDENTITY DIAGONAL, INCREMENT INDEX */

	/* IF NUCTRANSIT POPULATE DTHR MATRIX FOR ALL k <= WIDTH/2 */
	if (nuctransit) {
		DTHR_lookup[0] = match;		/* USE k=0 SLOT TO STORE MATCH VALUE */
		DTHR_lookup[1] = 100;
		DTHR_lookup[2] = 100;

		if (options[0][59]) {	/* IF opt_x EXTRA SQUEEZE: ALLOWS TRANSITION MATCHING FOR LOWER k AT CINCH_T */
			for (k=3; k <= WIDTH/2; k++) {
				numtransit = 1 + round(fractransit * k);
				DTHR_lookup[k] = 100*((k-numtransit)*match + numtransit*transition)/(k*match) - 1;
			}
		}
		else {					/* ELSE NOT */
			for (k=3; k <= WIDTH/2; k++) {
				numtransit = 1 + round(fractransit * k);
				DTHR_lookup[k] = 100*((k-numtransit)*match + numtransit*transition)/(k*match);
			}
		}
	}

	if (options[0][13]) {	/* opt_D: SHOW FIRST n DTHR VALUES */
		n = 30;	
		printf("\n Diagonal threshold (DTHR values) as a function of k (first %d values) and fractransit (%.2f):\n", n, fractransit);

		for (k = 1; k <= n; k++)
			printf("\n\tFor k = %2d and numtransit = %d, DTHR = %3d%%.", k, 1 + (int) round(fractransit*k), DTHR_lookup[k]);
		printf("\n\n Note: DTHR values are only populated if a sequence is specified.\n\n");
		exit(1);
	}

	for (n = 1; Seq[n] != '\0'; ) {
		for (m = 0; m < n; m++) {
			k = n-m;

			if (k > WIDTH/2) {
				m = n-WIDTH/2;
				k = n-m;
			}

			if (nuctransit) {
				DTHR = DTHR_lookup[k];
				imperfect_TR = 0;
			}

			homopoly_flag = 2;						/* HOMOPOLYMER RUN STATUS UNKNOWN */

			if (Seq[n] == '>')
				cit_new2D_width = a2D_n;

			if (k == 1) {	
				align2D[row][a2D_n++] = Seq[n++];		/* ENTER VALUE AT IDENTITY DIAGONAL, INCREMENT INDICES */
				break;									/* GO TO NEXT n */
			}
			else if (nuctransit == 0 && Seq[n] == Seq[m] && Seq[n] == Seq[n-1]) {
				align2D[row][a2D_n++] = Seq[n++];		/* ENTER VALUE AT IDENTITY DIAGONAL, INCREMENT INDICES */
                break;									/* GO TO NEXT n */
            } 

			if (pathbox[m][n] > mismatch && n+k <= lenseq) {
				for (i = m; i < n; i++) {
					Did = Did + pathbox[i][i];			/* COMPUTE SUM OF IDENTITY UNIT LINE */

					if (pathbox[i][i+k] == mismatch) {	/* STOP SHORT IF MISMATCH IS FOUND 		 		*/
						Dtr =  0;						/* B/C CURRENTLY ONLY CONSIDERING TRANSITIONS 	*/
						break;							
					}
					else
						Dtr = Dtr + pathbox[i][i+k];	/* COMPUTE SUM OF TANDEM REPEAT UNIT LINE */

					if (homopoly_flag && i >= m+1) {
						if ((nuctransit==0 || Seq[n-1]!=Seq[n]) && Seq[i] != Seq[i-1])
							homopoly_flag = 0;		/* COLLAPSE BIT FOR HOMOPOLYMER RUN: ABSENCE */
						else if (nuctransit && Seq[n-1]==Seq[n] && Seq[i+k] != Seq[i+k-1])
							homopoly_flag = 0;		/* COLLAPSE BIT FOR HOMOPOLYMER RUN: ABSENCE */
					}
				}

				if (homopoly_flag && i == n) {
					homopoly_flag = 1;					/* COLLAPSE BIT FOR HOMOPOLYMER RUN: DETECTED 	*/
														/* BIT IS THERE IF NEEDED BEYOND BREAK. 		*/
					align2D[row][a2D_n++] = Seq[n++];	/* ENTER VALUE AT IDENTITY DIAGONAL, INCREMENT INDICES */
					break;
				}

				if (nuctransit) { /* IF DNA AND IF CONSIDERING TRANSITIONS AS PARTIAL MATCHES */
					if (Dtr!=Did && 100*Dtr/Did > DTHR)	{	/* CHECK CASE OF IMPERFECT TR >= D THRESHOLD */
						imperfect_TR = 1;
						if (options[0][57]) {
							printf("\n DEV: Imperfect TR at n=%d.", n+1);
						}
						if (Seq[m] != Seq[n] && Seq[n] == Seq[n+k])	{ /* IMPERFECT TR AT n BUT PERHAPS PERFECT TR AT n+1 */
							imperfect_TR = 0;
						}
					}
				} /* END OF IF (CONSIDERING DNA TRANSITIONS) */

				/* ASCERTAIN ABSENCE OF BETTER DIAGONALS IN WINDOW OF TR, & SKIP TR IF SO. */
				/* CURRENTLY DOES NOT SEEM TO NEED Dtr == Did BUT LEAVING FOR CLARITY B/C UNSURE WHY YET */
				if (imperfect_TR || Dtr==Did) {
					for (l = 1; l < k; l++) {
						for (alt_k = 2*k; alt_k > k; alt_k--) {	/* MAGIC NUMBER 2*k FOR NOW TO SOLVE EXAMPLE; SEARCH MAGIC TO FIND OTHER STOPGAPS */
							alt_Dtr = 0;
							alt_Did = alt_k * match;
							alt_m = n+l - alt_k;

							for (j = 0; j < alt_k; j++) {
								if (pathbox[alt_m+j][(i=n+l+j)] > mismatch) {
									alt_Dtr = alt_Dtr + pathbox[alt_m+j][i];
								}
								else {
									alt_Dtr = 0;
									break;
								}
							} /* END OF FOR j LOOP TESTING ALT_TR AT n+l FOR ALT_k-MER */

							if (alt_Dtr && 100*alt_Dtr/alt_Did > 100*Dtr/Did) { /* TEST NEEDED FOR PERFECT OR IMPERFECT TR */
								if (options[1][57] > 1)
									printf("\n DEV: SKIPPING k-MER=%d at n=%d for alt_k-MER=%d at n=%d", k, n+1, alt_k, n+l+1);
								imperfect_TR = 0; 	/* BECAUSE ALTERNATIVE TR IS BETTER THAN CURRENT REGARDLESS OF PERFECTION  */
								Dtr = 0;			/* TO PREVENT ENTRY INTO TR LOOP */
								break;
							} 

						} /* END OF FOR alt_k LOOP */
					} /* END OF FOR l LOOP SCAN OF FUTURE SLIPS OF GREATER k WITHIN CURRENT k-MER TR WINDOW */
				} /* END OF IF IMPERFECT_TR */

				/* SKIP CINCH IF CAN AVOID CONFLICT VIA CYCLE FRAMING */
				if (nuctransit && Dtr==Did) {
					/* FIND IF THERE IS A PREVIOUS TR */
					for (i = n-1; i >= 0; i--) {
						if (sliploc[i])
							break;
					} 
					if ((i + (int) sliploc[i] * (j= (int) sliploc_nmer[i])) > m) {
						for (l = 0; n+l < lenseq && n+l < WIDTH; l++) {
							if (Seq[n+l] == Seq[n+k+l]) { 
								if (i%j == (n+l)%j) { 
									Dtr = imperfect_TR = 0;
									if (options[0][57]) {
										printf("\n DEV: Collapsing cinch to use later cycling frame at %d.", n+l+1);
									}
									break;
								}
							}
							else break;
						}
					}
				}

				/* SKIP CINCH IF IMPERFECT WHILE CONTAINING PERFECT TRs INSIDE */
				if (imperfect_TR) {
					for (l=k/2; l>1; l--) {
						for (i=0; i<l; i++) {
							if (Seq[n-l+i] != Seq[n+i])
								break;	
						}
						if (i==l) {
							imperfect_TR = 0;
							Dtr = 0;
							break;
						}
					}
				}

				if (Dtr==Did || imperfect_TR) {	/*  1st MEASUREMENT OF TANDEM REPEAT (TR) */
												/*  TR EQ. DIMER */
					unislip++;
					r = 1;
					sliploc[n]++;
					sliploc_nmer[n] = k = n-m;
					slips[k]++;
					TRcheck = 1;		/* FLAG TO KEEP CHECKING FOR ADD'L TRs */

					/* COUNT NUMBER of ADDITIONAL TANDEM REPEATS (TR) */
					while (TRcheck > 0) {
						if (m + (r+1)*k >= lenseq) {   	/* CHECK IF NEAR EDGE */
							Atr = Did = Dtr = 0;
							TRcheck = 0;					/* REINITIALIZE FLAG */
							break;
						}

						for (i = m; i < n; i++) {
							if ( (j=pathbox[i][(i + (r+1)*k)]) == mismatch) {	/* NEED TO BREAK IF PATHBOX POSITION IS MISMATCH =/= TRANSITION */
								Atr = 0;
								break;
							}
							else
								Atr = Atr + j;			/* j = pathbox[i][(i + (r+1)*k)] IS SET IN ABOVE IF TEST */
						}

						if (nuctransit) { /* IF DNA AND IF CONSIDERING TRANSITIONS AS PARTIAL MATCHES */
							if (Atr!=Did && 100*Atr/Did > DTHR)	{	/* CHECK CASE OF IMPERFECT TR >= D THRESHOLD */
								Aimperfect_TR = imperfect_TR = 1;
							}
							else
								Aimperfect_TR = 0;
						} /* END OF IF [CONSIDERING DNA TRANSITIONS] */

						if (Atr==Did || Atr==Dtr || Aimperfect_TR) {	/*  2nd MEASUREMENT OF A THIRD BLOCK */
												 						/*  ALSO ADD'L MEASUREMENTS > TRIMER */
							r++;
							sliploc[n]++;

							Atr = 0;				/* RESET TO COUNT POSSIBLE NEXT ONE */
						}
						else {
							Atr = Did = Dtr = 0;	/* RESET ALL DIAGONAL COUNTERS */
							TRcheck = 0;

							/* IF CYCLE REPEAT, STORE CYCLE RUN. CYCLIC REPEATS CAN BE REPEATS IN MORE THAN ONE FRAME. MUST BE >2k */
							i = 0;			/* CYCLE[] ARRAY INDEX */
							cycle_flag = 0;	
							for (j = -1; j < r; j++) {
								for (l = 0; l < k; l++) 
									cycle[i++] = Seq[(n + j*k + l)]; 	/* STORE WHOLE REPEAT */
							}
							for (l = 0; l < k-1; l++) {				/* STORE EXTENT OF PARTIAL REPEAT. CANNOT MATCH MORE THAN k */
								if (cycle[l] == Seq[(n + r*k + l)]) 	
									cycle[i++] = Seq[(n + r*k + l)];
								else
									break;
							}
							cycle[i] = '\0';
							cyclelen = i;

							if (cyclelen > 2*k) {
								cycle_flag = 1;		/* THIS BIT CAN BE USED TO ADD CYCLE NOTATION TO END OF LINE */
								if (options[0][57]) {
									printf("\n DEV: %2d-mer cycle sequence of length %2d starting at %4d: %s.", k, cyclelen, n-k+1, cycle);
								}
							}
							else {
								for (l = 0; l < WIDTH; l++)
									cycle[l] = '\0';
							}

							/* RECORD DNA "REVERB" IN SLIPLOC_ECHOES FOR ALL TR FRAMES */
							if (cyclelen == 2*k) 	h = 1;
							else 					h = k;

							for (j = 0; j < h; j++) {
								for (l = j; l+k <= i; l+=k) {
									if (sliploc_echoes[n-k+l] == blank) 		sliploc_echoes[n-k+l] = '(';
									else if (sliploc_echoes[n-k+l] == '(' ) 	sliploc_echoes[n-k+l] = '{';
									else if (sliploc_echoes[n-k+l] == '{' ) 	sliploc_echoes[n-k+l] = '[';
									else if (sliploc_echoes[n-k+l] == ')' ||
											 sliploc_echoes[n-k+l] == '}' ||
											 sliploc_echoes[n-k+l] == ']'   ) 	sliploc_echoes[n-k+l] = 'X';
	
									if (sliploc_echoes[n-1+l] == blank) 		sliploc_echoes[n-1+l] = ')';
									else if (sliploc_echoes[n-1+l] == ')' ) 	sliploc_echoes[n-1+l] = '}';
									else if (sliploc_echoes[n-1+l] == '}' ) 	sliploc_echoes[n-1+l] = ']';
									else if (sliploc_echoes[n-1+l] == '(' ||
											 sliploc_echoes[n-1+l] == '{' ||
											 sliploc_echoes[n-1+l] == '['   ) 	sliploc_echoes[n-1+l] = 'X';
								} /* END OF FOR l LOOP */
							} /* END OF FOR j LOOP */
						}

					} /* END OF else & TR_check = 0 PART */

					/* BELOW "SHADOW" IN COMMENTS REFERS TO REGION OF FIRST UNIT. I CALL IT THE UPSTREAM SHADOW B/C */
					/* OF THE WAY CINCH-T BEGINS MARKING A TR STARTING AT THE SECOND UNIT. 							*/
					/* THE OVERALL STRATEGY IS TO LOOP FROM n-1 DOWN TO m+1 AND COUNT VARIOUS USEFUL THINGS.		*/					
					/* COUNT MOST RECENT CONFLICTING SLIP LENGTH IN UPSTREAM SHADOW OF NEW k-MER 					*/
					/* THEN NEED TO COUNT OTHER BAD SLIPS IN SAME REGION IN CASE THESE WERE PREVIOUSLY TREATED		*/

					badslipspan = delta_to_bad1Dn = frst_badslip = bad_1Dn = overslip = recslips = scooch = oldbad = 0;
					mstop = m;
					for (i = n-1; i > mstop; i--) {							/* i WILL LOOP THROUGH TR SHADOW */
						if (sliploc[i]) {		/* DWELL ON POSITION W/ UNDAMPED SLIPS */
							/* BAD SLIP IS CLOSEST ONE TO n NOT HANDLED YET; BAD SLIP TYPE 1 */
							if (frst_badslip==0 && i - (int) sliploc_nmer[i] < m) {
								badslipspan = sliploc[i] * sliploc_nmer[i];	/* THIS IS THE SPAN LENGTH THAT WAS PREVIOUSLY TUCKED */
								delta_to_bad1Dn = n-i;						/* THIS IS DISTANCE	TO THE BAD SLIP */
								bad_1Dn = i;								/* STORE POSITION OF BAD SLIP */
								frst_badslip = 1;							/* TO MAKE SURE BADSLIPS ARE NOT COUNTED AGAIN LATER, */
								sliploc_echoes[i] = 'o';					/* MAKE A MARK INDICATING THIS ECHO'S BEEN DAMPENED   */
							}
                            /* THIS IS A BAD SLIP TYPE 2, BECAUSE OF NON-COINCIDENT TRANSITIONS  */
							else if (imperfect_TR && frst_badslip==0) {
	                            if (Seq[i]!=Seq[i+k] && Seq[i]==Seq[(i-(int) sliploc_nmer[i])]) {
	                                badslipspan = sliploc[i] * sliploc_nmer[i]; /* THIS IS THE SPAN LENGTH THAT WAS PREVIOUSLY TUCKED */
	                                delta_to_bad1Dn = n-i;                      /* THIS IS DISTANCE TO THE BAD SLIP */
	                                bad_1Dn = i;                                /* STORE POSITION OF BAD SLIP */
	                                frst_badslip = 1;                           /* TO MAKE SURE BADSLIPS ARE NOT COUNTED AGAIN LATER, */
	                                sliploc_echoes[i] = 'o';                    /* MAKE A MARK INDICATING THIS ECHO'S BEEN DAMPENED   */   
	                            }
								/* THIS IS BAD SLIP TYPE 2 EXCEPT MISMATCH IS OFFSET FWRD FROM sliploc[i] */
								else {
									for (l=1; l< (int) sliploc_nmer[i]; l++) {
										if ( Seq[i+l]!=Seq[i+l+k] && (Seq[i+l]==Seq[(i+l-(int) sliploc_nmer[i])]) )
											break;	
									}
									if (l< (int) sliploc_nmer[i]) {
										badslipspan = sliploc[i] * sliploc_nmer[i];	/* THIS IS THE SPAN LENGTH THAT WAS PREVIOUSLY TUCKED */
										delta_to_bad1Dn = n-i;						/* THIS IS DISTANCE	TO THE BAD SLIP */
										bad_1Dn = i;								/* STORE POSITION OF BAD SLIP */
										frst_badslip = 1;							/* TO MAKE SURE BADSLIPS ARE NOT COUNTED AGAIN LATER, */
										sliploc_echoes[i] = 'o';					/* MAKE A MARK INDICATING THIS ECHO'S BEEN DAMPENED   */
									}
								}
								/* THIS IS BAD SLIP TYPE 2 EXCEPT MISMATCH IS OFFSET BACK FROM sliploc[i] */
								if (frst_badslip==0) {
									for (l=1; l<= (int) sliploc_nmer[i]; l++) {
										if ( Seq[i-l]!=Seq[i-l+k] && (Seq[i-l]==Seq[(i-l+(int) sliploc_nmer[i])]) )
											break;	
									}
									if (l<= (int) sliploc_nmer[i]) {
										badslipspan = sliploc[i] * sliploc_nmer[i];	/* THIS IS THE SPAN LENGTH THAT WAS PREVIOUSLY TUCKED */
										delta_to_bad1Dn = n-i;						/* THIS IS DISTANCE	TO THE BAD SLIP */
										bad_1Dn = i;								/* STORE POSITION OF BAD SLIP */
										frst_badslip = 1;							/* TO MAKE SURE BADSLIPS ARE NOT COUNTED AGAIN LATER, */
										sliploc_echoes[i] = 'o';					/* MAKE A MARK INDICATING THIS ECHO'S BEEN DAMPENED   */
									}
								}
							}
							if (frst_badslip==0) {
								overslip = overslip + (sliploc[i] * sliploc_nmer[i]);	/* TOTAL SLIP LENGTH NOT UNDONE */
							}
							recslips = recslips + sliploc[i];				/* TOTAL SLIPS NOT UNDONE AS BADSLIPS */
						} /* END OF IF sliploc[i] */
					} /* END OF FOR i LOOP SWEEPING THROUGH POTENTIAL BADSLIP SHADOW */

					/* DEAL WITH OVERSLIPS NOT MARKED IN SHADOW */
					j=l=0;
					if (badslipspan) {
						for (i = bad_1Dn-1; i>=0; i--) {
							if (sliploc[i]) {
								j=(int) sliploc[i];
								break;
							}
						}
								if (j>1 && (i + (j-1) * (int) sliploc_nmer[i]) > m) {
									recslips = recslips + (j-1);
									a2D_n = a2D_n + (int) sliploc_nmer[i] * (j-1);
								}
								else if (j==1 && i>m && i + (int) sliploc_nmer[i] < bad_1Dn &&
														i - (int) sliploc_nmer[i] < m) {
									a2D_n = a2D_n - (int) sliploc_nmer[i];
								}
					}

					/* REWRITE UPSTREAM SHADOW TR AS SINGLE LINE IF ANY RECENT 'BAD' SLIPS FOUND THAT SPAN INTO NEW TR */
					if (badslipspan && bad_1Dn != m) {
						/* OVERWRITE THE OLD SLIP THAT IS IN THE WAY OF NEWER SLIP OF BIGGER k */
						a2D_n = a2D_n - delta_to_bad1Dn + overslip + sliploc_nmer[bad_1Dn]*(sliploc[bad_1Dn]-1);	/* REPOSITION TO OLD TR START */
						for (i = bad_1Dn-1; i > m; i--) {
							if (sliploc[i]) {
								a2D_n = a2D_n + sliploc_nmer[i]*(sliploc_nmer[i]-1);
							}
						}

						a2D_n = a2D_n + sliploc_nmer[bad_1Dn];                  /* REPOSITION TO AFTER OLD TR START */
						row = row - recslips;

						for (l = 0; l < delta_to_bad1Dn; l++) {					/* REFILL BADSLIPSPAN AS SINGLE ROW */
							align2D[row][a2D_n+l] = Seq[n-delta_to_bad1Dn + l];
						}
						a2D_n = a2D_n + delta_to_bad1Dn - overslip;	/* REPOSITION TO NEW STRETCHED POSITION */
					}	/* END OF IF BADSLIP */

                    if (sliploc[n-k] > 1) { 
                        overslip = overslip + (sliploc[n-k]-1) * sliploc_nmer[n-k];
                    }  

					sloc_z = sliploc[bad_1Dn];
                    for (l = 1; l < (z=sliploc_nmer[bad_1Dn]); l++) {
                        if (sliploc[bad_1Dn-z+l] > 0) { 
                            a2D_n = a2D_n - sliploc_nmer[bad_1Dn-z+l];
                        }    
                    }   

					for (i = 0; i < r; i++) {
						while (align2D[row][a2D_n + scooch] != '\0') {
                            scooch++;
                        }

						align2D[row  ][a2D_n+scooch  ] = '/'; 
						align2D[row  ][a2D_n+scooch+1] = '\0'; 	
						row++;										/* <==== ROW INCREMENTED HERE!		*/

                        if ( (z = k - overslip - (int) a2D_n) > 0) { 
							if (options[0][57]) { /* opt_v VERBOSITY */
    	                        printf("\n DEV: 2D-printing into upstream region by %d.", z);
							}
                            for (m = 0; m < row; m++) {
                                for (l = WIDTH; l >= 0; l--) 
                                    align2D[m][l+z] = align2D[m][l];    /* PUSH PREVIOUS ROWS TO THE RIGHT BY z    */   
                                for (l = 0; l < z; l++) {
                                    align2D[m][l] = blank;              /* FILL IN z COLUMNS WITH BLANKS           */   
                                }    
                            }    
                                a2D_n = a2D_n + z;                      /* ADJUST FUTURE a2D_n by + z              */   
                        }  
						else z=0;

						for (j = 0; j < (int) a2D_n; j++) 
							align2D[row][j] = blank;
						if (badslipspan && i==0) {			/* FLIP TO LOWERCASE IN FIRST UNIT */	
							if (islower(align2D[row-1][a2D_n - k + overslip - 1]))			/* IF OVERLAPPING OLD BADSLIP */
								align2D[row-1][a2D_n - k + overslip - 1] = Seq[m-1];		/* MAKE UPPERCASE */
							for (j = 0; j < k; j++) { 
								align2D[row-1][a2D_n - k + j + overslip] = Seq[m+j] + 32;
							}
 								align2D[row-1][a2D_n - k + j + overslip   ] = '/';
 								align2D[row-1][a2D_n - k + j + overslip +1] = '\0';
						}
						for (j = 0; j < k; j++) {
							align2D[row  ][a2D_n - k + j + overslip] = Seq[m+j+k*(i+1)];	
						}
					} /* END OF FOR i */
					if (badslipspan) {
						for (i=1; row+i < row+recslips; i++) {
							for (j=0; j < a2D_n + badslipspan; j++) {
								align2D[row+i][j] = '\0';
							}
						}
						sliploc[bad_1Dn] = 0;
					}

					if (imperfect_TR) {		/* FIND THE TRANSITIONS AND STORE TYPE (Y OR R) IN NEW CORRECTED POSITIONS */
						for (j=0; j < r; j++) {
							for (i=0; i < k; i++) {
								if (Seq[m+i] != Seq[n+i+j*k]) {
									if (Seq[m+i] == 'A' || Seq[m+i] == 'G') {
										align2D[MAXROW][a2D_n-k+i+overslip] = transtype = 'R';
										SeqRY[n+i+j*k] = transtype;
									}
									else if
									   (Seq[m+i] == 'C' || Seq[m+i] == 'T') {
										align2D[MAXROW][a2D_n-k+i+overslip] = transtype = 'Y';
										SeqRY[n+i+j*k] = transtype;
									}
									SeqRY[m+i] = transtype;
								}
							}
						} /* END OF FOR i SCAN FOR TRANSITIONS LOOP */
					} /* END OF IF IMPERFECT TANDEM REPEAT (TR) */

					n = n + r*k;
					a2D_n = a2D_n + overslip;
					r = 0;
					break;
				}
				else 				/* CASE WHEN Dtr != Did AND NOT IMPERFECT-TYPE TR */
					Dtr = Did = 0;
			}
		} /* END OF "for m" LOOP */
	} /* END OF "for n" LOOP => cinch_t PASS COMPLETE */

 	align2D[row][cit_new2D_width+1] = '\0';
	options[1][2] = options[1][32] = cit_new2D_width;	/* ASSIGN CINCH-WIDTH TO HISTORY [0--9] AND LCL CURRENT [32]	*/
	clear_right(align2D, options);

	/*************************************************************************/
	/* PRINT VALUES OF 1-D ALIGN BOX *****************************************/

	blocks = count_wrap_blocks(lenseq, options[1][58]);

	mha_head(Seq_name, lenseq, options);

	for (j = 0; j < blocks; j++) {
		line_end(START, j+1, options, 0);	
   		for (n = j * options[1][58]; (n < (j+1) * options[1][58]) && (Seq[n] != '>') && (Seq[n] != '\0'); n++) {
			printf("%1c", Seq[n]);
		}
		if (Seq[n] == '>') {
			printf("%1c", Seq[n]);  /* PRINTS TERMINAL CHARACTER '>' */
			if (options[0][21]) 
				line_end(END, n, options, 0);
			else
				printf("\n");
		}
		else {
			printf(" ");
			if (options[0][21]) 
				line_end(END, n, options, 0);
			else 
				printf("\n");
		}

		if (options[0][47]) {		   /* OPTION opt_l TO SHOW SLIP LOCATIONS */
			/**********************************************/
			line_end(SLIPS, 0, options, 0);
			for (i = j * options[1][58]; (i < (j+1) * options[1][58]) && (i < lenseq); i++) {
				if (sliploc[i] == 0)
					printf(".");
				else
					printf("%c", mha_base62(sliploc[i]));
			}
			if (j+1 == blocks)
				printf(" <==== # of TRs > 1 (base 62)\n");
			else
				printf("\n");
			/**********************************************/
			line_end(SLIPS, 0, options, 0);
			for (i = j * options[1][58]; (i < (j+1) * options[1][58]) && (i < lenseq); i++) {
				if (sliploc_nmer[i] == 0)
					printf(".");
				else {
					printf("%c", mha_base62(sliploc_nmer[i]));
				}
			}
			if (j+1 == blocks)
				printf(" <==== TR unit size (base 62)\n");
			else
				printf("\n");
			/**********************************************/
			line_end(SLIPS, 0, options, 0);
			for (i = j * options[1][58]; (i < (j+1) * options[1][58]) && (i < lenseq); i++) {
				printf("%c", sliploc_echoes[i]);
			}
			if (j+1 == blocks)
				printf(" <==== (((( {  REVERB  } ))))\n");
			else
				printf("\n");
			/**********************************************/
			head_start = (j * options[1][58]) % 10;
			if (j+1 < blocks) {
				line_end(SLIPRULER, head_start, options, options[1][58]);
				printf("\n");
			}
			else 
				line_end(SLIPRULER, head_start, options, lenseq-options[1][58]*(blocks-1));
			/**********************************************/
		}   /* END OF opt_l PRINT MODULE */ 

	} /* END OF FOR j LOOP */

	/********** 2. cinch_t MODULE: WRAPS LARGEST EXACT k-mers, IGNORES INTRA-TR TRs **********/
	i = ++options[1][18];	/* INCREMENT opt_I COUNTER (pass_num) */
	print_2Dseq(align2D, cit_new2D_width, ptr_Seq_name, options);
	passQ[i] = options[0][10];

	if (options[1][48]!=0 && options[1][49]!=0)
		tucksense(align2D, options);

	/********** 3. cinch_l MODULE: WRAPS HOMOPOLYMERIC RUNS IF >= 20 (2 * wrap VAR.) ********/
	i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */
	if (cinch_l(align2D, options))
		print_2Dseq(align2D, options[1][32], ptr_Seq_name, options);
	passQ[i] = options[0][10];

	/********* 4. cinch_k MODULE: HANDLES k-mers FROM SIZE WIDTH/2 DOWN TO k=1 ***********/
	i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */

	cinch_k(align2D, options, DTHR_lookup);
	cycle_flag = print_2Dseq(align2D, options[1][32], ptr_Seq_name, options);
	passQ[i] = options[0][10];

	if (cycle_flag) {
		while (tucksense(align2D, options))
			cycle_flag = print_2Dseq(align2D, options[1][32], ptr_Seq_name, options); 
	}

	/********* 5. cyclelize MODULE: HANDLES k=2 CYCLING, "FUDGE-CYCLELIZES"/PUSHES k>2 ***************/
	continue_flag = 1;

	if (continue_flag) {
		i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */
		if (cycle_flag || align2D[0][0] == blank) {
			cyc_runs++;
			mha_writeback(align2D, scratch, options);		/* 1. COPY align2D TO scratch ARRAY */
			mha_writecons(align2D, scratch, options);

			while ((go_flag=cyclelize(scratch, ptr_Seq_name, options)) != 0 && cyc_runs < CYCMAX) {	
				cyc_runs++;
			}
			if (go_flag == 0 || go_flag > cycle_flag) {
				mha_writeback(scratch, align2D, options);
				mha_writecons(scratch, align2D, options);
			}
			else {
				warnhead('R');
				printf("Reverting to post cinch-k due to micro-foam turbidity.\n");
				options[1][i] = options[1][32];
				cyc_runs = cyc_runs + CYCMAX*1000;		/* cyc_runs > CYCMAX IF PASS REVERTED */
			}
		}
		else {									/* 2. NO CYCLELIZE SO WRITE scratch BACK TO align2D */
			options[1][i] = options[1][i-1];
		}
		passQ[i] = options[0][10];
	}

	/********* 6. cinch_s MODULE: HANDLES TR k-mers IN A SINGLE LINE ***********/
	if (continue_flag) {
		i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */
		if (cinch_s(align2D, options)) {
			print_2Dseq(align2D, options[1][32], ptr_Seq_name, options); 
		}
		else {
			if (options[0][20]) {
				options[0][20] = 0;	/* TEMPORARY ASSIGNMENT TO PREVENT PRINTING OF CONSENSUS ROW */
				consensus_2D(align2D, options, 0, options[1][32]);
				options[0][20] = 1;	/* REASSIGN SETTING */
			}
			else
				consensus_2D(align2D, options, 0, options[1][32]);
		}
		passQ[i] = options[0][10];
	}

	/********* 7. cinch_d MODULE: HANDLES DE NOVO INTER-TR REPEATS *********************************/
	if (continue_flag) {
		i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */

		if (options[1][57]) {
			if (nuctransit) {
				mha_head(Seq_name, options[1][32], options);
				printf(" Pre-cinch_d report (p = perfect, i = imperfect tandem repeat): \n");
			}
			else {
				mha_head(Seq_name, options[1][32], options);
				printf(" Pre-cinch_d report: \n");
			}
		}
		intraTR_reps_tot = intraTR_reps = cinch_d(align2D, ptr_Seq_name, options, 0);
	
		if (intraTR_reps_tot == 0) {
			printf(" Nothing left for cinch-d to cinch! \n");
			print_2Dseq(align2D, options[1][32], ptr_Seq_name, options);
		}
	
		if (intraTR_reps_tot > 0) {
			while (intraTR_reps > 0) {
				intraTR_reps = cinch_d(align2D, ptr_Seq_name, options, 1);
				d_runs++;
			}
			options[0][7] = --d_runs;	/* STORE d_runs in PASS SLOT IN CASE FUTURE MODULES POST CINCH-D NEED IT */
							/* NEED TO COUNT LIKE THIS FOR FOLLOWING REASONS: 			*/
							/*  1. LAST RUN IS A CHECK RETURNING ZERO, SO RUNS NEED TO BE DECREMENTED BY 1. 	*/
							/*  2. APPARENT LAST EFFECTIVE CINCH-D RUN MAY REVEAL A NEW CINCH-D OPPORTUNITY. 	*/
							/* 	OTHERWISE COULD HAVE SET WHILE LOOP TO > 1 (DON'T DO THIS. )			*/
		}
		passQ[i] = options[0][10];
	}

	/********* 8. relax_2D MODULE: DE-CINCHES HOMOPOLYMER RUNS IF THEY DID NOT AID CINCH-D *******/
	if (continue_flag && options[0][49]<1) {		/* opt_n DO NOT DO RELAX-2D */
		i = ++options[1][18];		/* INCREMENT opt_I ++PASS NUM */
	
		do {
			relax_length = options[1][32] ;
			relax_2D(align2D, options);
			relax_length = options[1][32] - relax_length;
			r_runs++;
		}
		while (relax_length > 0);
	
		print_2Dseq(align2D, options[1][32], ptr_Seq_name, options);
		passQ[i] = options[0][10];
	}	

	/*************************************************************************/
	/* OPTION TO PRINT VALUES OF RECOVERED 1-D ALIGN BOX *********************/

	if (options[0][27]) { /* opt_R is 1 if ON */ 
		printf("\nChecking 1-D recovery from 2-D self-MHA:\n");
		recover_1D(recovered, align2D, options);
		blocks = count_wrap_blocks(lenseq, options[1][58]);

		k = lenseq;
		l = (int) strlen(recovered) - 1;
		m = max(k,l);
		z = 0;			/* USE TO COUNT RECOVERED LETTERS IDENTICAL TO 1D */

		if (k != l) {
			warnhead('d');
			printf("1-D sequence differs in length from 2-D recovered 1-D sequence by %d letters.\n", k-l); 
		}
	
		for (j = 0; j < blocks; j++) {
			printf("\n");
			line_end(START, j+1, options, 0);	
	   		for (n = j * options[1][58]; n < (j+1) * options[1][58] && Seq[n] != '>'; n++) {
				printf("%1c", Seq[n]);
			}
			if (Seq[n] == '>') {
				printf("%1c", Seq[n]);  /* PRINTS TERMINAL CHARACTER '>' */
				if (options[0][21]) 
					line_end(END, k, options, 0);
				else
					printf("  <== 1-D \n");
			}
			else {
				printf(" ");
				if (options[0][21]) 
					line_end(END, n, options, 0);
				else 
					printf("\n");
			}
	
			line_end(SLIPS, j+1, options, 0);	
	   		for (n = j * options[1][58]; n < (j+1) * options[1][58] && n < m; n++) {
				if (seqtype == 1 || seqtype == 2) {			/* IF DNA OR RNA */
					if (Seq[n] == 'n' || recovered[n] == 'n') {
						printf("?");
						z++;	/* Added in version 2.76 */
					}
					else if (Seq[n] == recovered[n]) {
						printf("|");
						z++;
					}
					else {
						printf("*");
						recovery_flag++;
					}
				}
				else {						/* ELSE IF NOT NA */
					if (Seq[n] == recovered[n]) {
						printf("|");
						z++;
					}
					else {
						printf("*");
						recovery_flag++;
					}
				}
			}
			printf("\n");
	
			line_end(START, j+1, options, 0);	
	   		for (n = j * options[1][58]; (n < (j+1) * options[1][58]) && recovered[n] != '>'; n++) {
				printf("%1c", recovered[n]);
			}
			if (recovered[n] == '>') {
				printf("%1c", recovered[n]);  /* PRINTS TERMINAL CHARACTER '>' */
				if (options[0][21]) {
					line_end(END, l, options, 0);
				}
				else
					printf("  <== 2-D\n");
			}
			else {
				printf(" ");
				if (options[0][21]) 
					line_end(END, n, options, 0);
				else 
					printf("\n");
			}
		} /* END OF FOR j LOOP */
		printf("\n");

		options[1][9] = z;			/* STORE NUMBER OF RECOVERED CHARACTERS */

		if (recovery_flag) {		/* LAST ROW OF array2D WILL STORE CONSENSUS, SO NEED TO KEEP CLEAR */
			warnhead('R');
			printf("Imperfect recovery of 1-D sequence from 2-D self-MHA.\n");
			passQ[9] = 1000*(lenseq-recovery_flag)/lenseq;
		}
		else {
			line_end(SLIPS, 0, options, 0);
			printf("Perfect recovery of 1-D sequence.\n");
			passQ[9] = 1000; 
		}
	} /* END of opt_R */

	/* PRINT OPTION FOR K-MER REPORT AFTER main() cinch_t **********************/
	if (options[0][46]) {	/* opt_k */

		int k_start=2;

		printf("\n Unique slips: %4d\n", unislip);

			if (WIDTH/2 <= 20) {
				for (i = k_start; i <= WIDTH/2; i++)
					printf(" %smers:%3d,", nmer_prefix(i), slips[i]);
				printf(" ...\n");
			}
			else {
				for (i = k_start; i <= 10; i++)
					printf(" %smers:%3d,", nmer_prefix(i), slips[i]);
				printf("\n");
				for (i = 11; i <= 20; i++)
					printf(" %smers:%2d,", nmer_prefix(i), slips[i]);
				printf(" ...\n");
			}
	}
	/***************************************************************************/

	c   = options[1][1];		/* REUSING c VAR FOR FORMATTED STRING LENGTH */
	if (options[1][33] > 1)		/* opt_XX FISHER-YATES RANDOMIZATIAN */
		row = FY_size;			/* USE IN PLACE OF ORIGINAL STRING LENGTH */
	else
		row = options[1][0];	/* RESUING row VAR FOR ORIGINAL STRING LENGTH */

	if (c != row && (options[0][57] || options[0][50]) ) {	/* opt_v VERBOSITY, OR		  */
		warnhead('F');										/* opt_o SHOW ORIGINAL STRING */
		printf("%d character(s) removed during formatting.\n", row-c); 
	}

	/***************************************************************************/
	if (seqtype == 1)		
		printf(   "\n   PASS      Width cinch history for %s %c (DNA)\n", Seq_head, Seq_name);
	else if (seqtype==2)	
		printf(   "\n   PASS      Width cinch history for %s %c (RNA)\n", Seq_head, Seq_name);
	else if (seqtype==3)	
		printf(   "\n   PASS      Width cinch history for %s %c (PROTEIN)\n", Seq_head, Seq_name);
	else if (seqtype==0)
		printf(   "\n   PASS      Width cinch history for %s %c (BABYLONIAN)\n", Seq_head, Seq_name);
	else
		printf(   "\n   PASS      Width cinch history:\n");

	for (i = 0; options[1][i] != '\0' && i < 9; i++) {
		printf("  %5d       => %4ld ", passQ[i], options[1][i]);
		switch (i) {
		case 0:
			if (options[0][33] == 0)
				printf("characters in original string\n");
			else if (options[1][33] == 1)
				printf("characters in original string => RANDOMIZED\n");
			else if (options[1][33] > 1)
				printf("characters in original string => FISHER-YATES RANDOMIZED TO FILL %d\n", FY_size);
			break;
		case 1:
			printf("%s post cleanseq  [pass #1]\n", letr_unit);
			break;
		case 2:
			printf("%s post cinch-t   [pass #2]\n", letr_unit);	/* STYLE: USE DASHED NAME FOR PRINTING, UNDERSCORED FOR CODING 	*/
			break;												/*  "cinch-x" VERSUS 'cinch_x' programming calls 				*/ 
		case 3:													/*  USEFUL FOR SEARCHING CODE.									*/
			printf("%s post cinch-l   [pass #3]\n", letr_unit);
			break;
		case 4:	
			if (options[0][4] == 0)
				printf("%s post cinch-k   [pass #4]\n", letr_unit);
			else if (options[0][4] == 1)
				printf("%s post cinch-k   [pass #4: one row added]\n", letr_unit);
			else if (options[0][4] > 1)
				printf("%s post cinch-k   [pass #4: %ld rows added]\n", letr_unit, options[0][4]);
			break;
		case 5:	
			k = options[0][5];
			if (k == 0)
				printf("%s post cyclelize [pass #5]\n", letr_unit);
			else if (cyc_runs == 1) {	
				if (k == 2)
					printf("%s post cyclelize [pass #5: %d run; cyclelized TR of type k=2]\n", letr_unit, cyc_runs);
				else if (k == 3)
					printf("%s post cyclelize [pass #5: %d run; nudge-cyclelized TR of type k>2]\n", letr_unit, cyc_runs);
				else if (k == 4)
					printf("%s post cyclelize [pass #5: %d run; tip-cyclelized TR]\n", letr_unit, cyc_runs);
			}
			else if (cyc_runs <= CYCMAX) {	/* IN WHICH CASE k WILL BE NON-ZERO */
				if (k == 2)
					printf("%s post cyclelize [pass #5: %d runs; last cyclelized TR was of type k=2]\n", letr_unit, cyc_runs);
				else if (k == 3)
					printf("%s post cyclelize [pass #5: %d runs; last TR was nudge-cyclelized (k>2)]\n", letr_unit, cyc_runs);
				else if (k == 4)
					printf("%s post cyclelize [pass #5: %d runs; last TR was tip-cyclelized]\n", letr_unit, cyc_runs);
			}
			else 							/* IN WHICH CASE cyc_runs >> CYCMAX */
				printf("%s post cyclelize [pass #5: reverted after %d runs due to gnarly micro-foam]\n", letr_unit, cyc_runs-CYCMAX*1000);
			break;
		case 6:	
			if (options[0][6] == 0)
				printf("%s post cinch-s   [pass #6]\n", letr_unit);
			else if (options[0][6] == 1)
				printf("%s post cinch-s   [pass #6: one row added]\n", letr_unit);
			else if (options[0][6] > 1)
				printf("%s post cinch-s   [pass #6: %ld rows added]\n", letr_unit, options[0][6]);
			break;
		case 7:	
			if (d_runs > 0)
				printf("%s post cinch-d   [pass #7: %d runs]\n", letr_unit, d_runs);
			else
				printf("%s post cinch-d   [pass #7]\n", letr_unit);
			break;
		case 8:	
			printf(    "%s post relax-2D  [pass #8: relaxed %d runs]\n", letr_unit, r_runs);
			break;
		}
	}
	if (options[0][27]) { /* opt_R */ 
		printf("  %5d       => %4ld ", passQ[9], options[1][9]);
			printf(    "%s recovered 1D   [final check pass]\n", letr_unit);
	}

	if (continue_flag) {
		printf("\n Width cinch ratio (post cinch-d):  %2.3f", ratio1=(float)options[1][7]/lenseq);
		if (options[0][49] < 1)
			printf("\n Width cinch ratio (post relax-2D): %.3f\n\n", ratio2=(float)options[1][32]/lenseq);
		else
			printf("\n\n");
	}
	else
		printf("\n Width cinch ratio (post cinch-k): %.3f\n\n", ratio1=ratio2=(float)options[1][4]/lenseq);

	if (options[0][24]) {								/* OPTION TO OUTPUT 2-D ALIGNMENT & CONSENSUS STRING TO FILE */
		align2D[MAXROW][options[1][32]] = '\0';		/* MAKE SURE CONSENSUS ROW IS TERMINATED AT CORRECT POSITION */
		fp_cons = fopen("Surf_barrels.log", "a");		/* FOPEN RIGHT BEFORE WRITING TO MINIMIZE CHANCE OF CLOSING W/ OPEN FILES */
		fprintf(fp_cons,">%s (%d > %d %s) x%d\n", 
			file_name, (int) options[1][1], (int) options[1][7], letr_unit, (int) options[1][59]);	

		for (m = 0; align2D[m][0] != '\0'; m++) {
			fprintf(fp_cons, " %s\n", align2D[m]);
		}

		fprintf(fp_cons, " %s\n\n", align2D[MAXROW]);
		fclose(fp_cons);
	}
	if (msa&&options[0][24]==2) {
		char ch=blank;
		n = get2Dtucknum(m2Dalig, align2D, options)+1;

				m2Dalig[  n][0] = '>';
		for (i=0; align2D[i][0]!='\0'; i++) {
			if (i>0)
				m2Dalig[i+n][0]= ' ';
			for (j=1; (ch=align2D[i][j-1])!='/' && ch!=optR && ch!='>'; j++) {
				m2Dalig[i+n][j] = ch;
			}
				m2Dalig[i+n][j] = ch;
		}
		fp_msa = fopen("TUBES.mha", "w");			/* FOPEN RIGHT BEFORE WRITING TO MINIMIZE CHANCE OF CLOSING W/ OPEN FILES */
			fprintf(fp_msa, "%s\n", m2Dalig[0]);
		for (m = 1; m2Dalig[m][1] != '\0'; m++) {
			fprintf(fp_msa, " %s\n", m2Dalig[m]);
		}

		fclose(fp_msa);
	}
	else if (options[0][24]==2) {					/* OPTION TO APPEND 2-D ALIGNMENT TO MSA FILE */
		fp_msa = fopen("TUBES.mha", "a");			/* FOPEN RIGHT BEFORE WRITING TO MINIMIZE CHANCE OF CLOSING W/ OPEN FILES */
			fprintf(fp_msa, ">%s\n", align2D[0]);
		for (m = 1; align2D[m][0] != '\0'; m++) {
			fprintf(fp_msa, " %s\n", align2D[m]);
		}

		fclose(fp_msa);
	}
	

	if (options[0][54] != 1) {		/* ONLY IF opt_s OPTION TO SILENCE OUTPUT IS NOT ON */
		if (options[1][39])			/* dev_notes always written here but can flag from within functions outside of main */
			strcpy(dev_notes,"connudge");

		fp_out = fopen("Surf_wavereport.mha", "a");		/* FOPEN RIGHT BEFORE WRITING TO MINIMIZE CHANCE OF CLOSING WITH OPEN FILES */
		fprintf(fp_out, "maximal v%s\t%.24s\t%4ld\t%.3f\t-x%d\tCYC:%2d (k=%ld)\tRND:-%.*s\t%c %s (%d %s) REC:%4d\t%d\t%s\n", 
				version, ctime(&lcl_time), options[0][10], ratio1, (int) options[1][59], cyc_runs, options[0][5], 
				(int) options[1][33], "XX", Seq_name, file_name, (int) options[1][1], letr_unit, passQ[9], (int) options[1][7], dev_notes);
		fclose(fp_out);

		/* IF IMPERFECT CONSENSUS OR IF CYCLELIZE REVERTED */
		if (options[0][10] != 1000 || cyc_runs > CYCMAX) {
			fp_tricksy = fopen("waves/foam_and_chowder.mha", "a");
			fprintf(fp_tricksy, "maximal v%s\t%.24s\t%4ld\t%.3f\t-x%d\tCYC:%2d (k=%ld)\tRND:-%.*s\t%c %s (%d %s) REC:%4d\t%s\n", 
					version, ctime(&lcl_time), options[0][10], ratio1, (int) options[1][59], cyc_runs, options[0][5], 
					(int) options[1][33], "XX", Seq_name, file_name, (int) options[1][1], letr_unit, passQ[9], dev_notes);
			for(n = 0; Seq[n] != '\0'; n++) {
				if (Seq[n] != 10 && Seq[n] != 13 && Seq[n] != EOF)
					fprintf(fp_tricksy, "%c", Seq[n]);
			}
			fprintf(fp_tricksy, "\n");
			fclose(fp_tricksy);
		}
	}

	exit(0);
} 
/************************** END OF MAIN() *********************************/


/*** FUNCTION 01 ************************************************************************************/
short unsigned int cinch_l(char align2D_pass3[][MAXROW], long int loptions[][62]) 
{
int cil_row=0, i=0, j=0, k=0, l=0, m=0, n=0, run=0, x=0;
char letr;
char lopt_Q_left = (char) loptions[1][26];		/* LHS character delimiter for homopolymer run */
char lopt_R_rght = (char) loptions[1][27];		/* RHS character delimiter for homopolymer run */
char blnk        = (char) loptions[1][11];		/* opt_B blank character		*/
int cil_mwrap    = (char) loptions[1][22];		/* opt_M long_homopolymer_run	*/
char cil_align2D[MAXROW+1][MAXROW];
int count_wrap_blocks(int lcl_width, int lcl_opt_w);
void line_end(int type, int c, long int lend_options[][62], int lcl_width);
void mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
void mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
short unsigned int print_2Dseq(char align2D_print[][MAXROW], int print_lenseq2D, char *printptr_Seq_name, long int poptions[][62]);
short unsigned int cinchled=0;	/* BIT FLAG TO SAY cinch_l DID SOMETHING */

	for (m = 0; align2D_pass3[m][0] != '\0'; m++) {
		for (n = 0; align2D_pass3[m][n] != '\0'; n++) {

			while (align2D_pass3[m][n] == blnk) {	/* MOVE WINDOW PAST BLANKS */
				cil_align2D[m+cil_row][n-x] = blnk;
				n++;			  
			}
			  /* letr ASSIGNED */
			if ((letr=align2D_pass3[m][n]) != '/' || letr != '\0') {	/* letr != '\0' USEFUL IF BLANKS FILLED TO SCRIMMAGE */ 
				cil_align2D[m+cil_row][n-x] = letr;	
				if (letr == align2D_pass3[m][n-1]) 
					run++;
				else
					run = 1;
			}
			else
				cil_align2D[m+cil_row][n-x] = letr;		/* WRITE TERMINAL MHA CHARACTERS '/' or '>' */
	
			if (run == 2*cil_mwrap) {			/* TRIGGER LENGTH MEASURE OF MONOMERIC RUN	*/
												/* AND WRITING OF (10)-MER BLOCKS			*/
				cinchled = 1;
				while (align2D_pass3[m][n+run-2*cil_mwrap] == letr) {
					run++;
				}
				
				cil_align2D[    m+cil_row][n-x-cil_mwrap+1  ] = lopt_R_rght;	/* WRITE SLIP AFTER FIRST 10 */
				for (l = 1; l < cil_mwrap; l++)								
					cil_align2D[m+cil_row][n-x-cil_mwrap+1+l] = '\0';
				cil_row++;														/* ADVANCE TO NEXT ROW, CONT.*/
	
				for (j = (run-cil_mwrap)/cil_mwrap; j > 0; j--) {
					for (i = 0; i < n-x-2*cil_mwrap; i++)
						cil_align2D[m+cil_row][i] = blnk;
					cil_align2D[m+cil_row][i] = lopt_Q_left;					/* MARK LEFT EDGE OF MWRAP RUN */
					for (k = 0; k < cil_mwrap; k++) {
						cil_align2D[m+cil_row][n-x-2*cil_mwrap+1+k] = letr; /* FILL W/ MONOMER LETTER  */
					}
					if (j > 1) {
						cil_align2D[m + cil_row  ][n-x-cil_mwrap+1] = lopt_R_rght;
                        cil_align2D[m + cil_row++][n-x-cil_mwrap+2] = '\0';             /* SHOULD NOT BE NECESSARY */
					}
					else if (j == 1) 
						break;			/* BREAK OUT OF FOR j LOOP */
				}

				n = n + cil_mwrap*(run/cil_mwrap) - 2*cil_mwrap - 1;		/* ADVANCE n TO JUST PAST LAST 10-MER BLOCK */
				x = x + cil_mwrap*(run/cil_mwrap) - cil_mwrap;				/* INCREMENT x TO REFLECT TUCK */
				
				run = 1;	/* RESET TO AVOID RETRIGGERING WRAPPING FOR REMAINDER OF RUN < mwrap */
			}	/* END OF MONOMERIC RUN BLOCK WRITES */ 
		}   /* END OF FOR n LOOPS */ 
	}   /* END OF FOR m LOOPS */

	loptions[1][3] = loptions[1][32];
	if (cinchled) {
		mha_writeback(cil_align2D, align2D_pass3, loptions);
		return(1);
	}
	else 
		return(0);
}


/*** FUNCTION 02 ************************************************************************************/
void cinch_k(char align2D_pass4[][MAXROW], long int koptions[][62], int DTHR_lookie[WIDTH/2]) 
{
int cik_row=0, i=0, k=0, l=0, m=0, n=0, scrimmage_line = -1, x=0, y=0, r=0; 
int first_mwrap_start=0, last_mwrap=0;
unsigned short int first_mwrap=0, keep_checking=1;
unsigned short int nuctype = koptions[1][13];		/* EQUALS ONE IF DNA STRING, TWO IF RNA, THREE IF PROTEIN */
unsigned short int nuctransit=0, check_imperf=0;	/* BIT FLAG FOR HANDLING NUCLEOTIDE TRANSITIONS SILENTLY (IGNORING) */
unsigned short int homopolyflag=0, imperfect_TR=0, piso_nuctransit=0;
int match 		= DTHR_lookie[0];
int sum4score;		/* SCORE VAR FOR IMPERFECT TR'S */
char letr, letr2, letr3;
char blnk        = (char) koptions[1][11];		/* opt_B blank character */
int  cik_mwrap   =        koptions[1][22];		/* opt_M long_homopolymer_run */
char kopt_Q_left = (char) koptions[1][26];		/* LHS character delimiter for homopolymer Run */
char kopt_R_rght = (char) koptions[1][27];		/* RHS character delimiter for homopolymer Run */
char cik_align2D[MAXROW+1][MAXROW] = {{0}};
int x_history[MAXROW] = {0};		/* STORE HISTORY OF x VARIABLE VIA POSITION n */
void               clear_2D_ar(char wipe_align2D[][MAXROW]);
unsigned int       consensus_2D(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width);
int                count_wrap_blocks(int lcl_width, int lcl_opt_w);
void               line_end(int type, int c, long int lend_options[][62], int lcl_width);
int 			   col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown); 
void 			   mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
void               mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
void               mha_UPPERback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
short unsigned int print_2Dseq(char align2D_print[][MAXROW], int print_lenseq2D, char *printptr_Seq_name, long int poptions[][62]);

	if (nuctype == 1) {		/* IF DNA */
		nuctransit = 1;
		piso_nuctransit = 2;		/* SET THE FLOOR >2 FOR k-mer LOOP */	
	}

	/* START AT BIGGEST k-MER POSSIBLE AT 2x */
	for (k = WIDTH/2; k > 0; k--) {

		cik_row = x = y = 0;	
		scrimmage_line = -1;

		for (m = 0; align2D_pass4[m][0] != '\0'; m++) {

			/* CHECK LINE AHEAD OF TIME FOR FIRST MONO-RUN TERMINATOR */
			for (i = 0; (letr=align2D_pass4[m][i]) != '\0'; i++) {	
				/* NOTE: APPARENT REDUNDANCY BELOW ALLOWS CHANGING L/R RUN DELIMITERS TO BE EQUAL */
				if (letr == kopt_R_rght && align2D_pass4[m][i-1] != blnk) {
					first_mwrap = 1;					/* TURN ON BIT FLAG FOR HANDLING FIRST MONO-RUN */
					last_mwrap = 1;						/* TURN ON BIT FLAG FOR HANDLING TO LAST MONO-RUN */
					first_mwrap_start = i-cik_mwrap;	/* SAVE POSITION OF START OF MWRAP */
					break;								/* BREAK OUT OF FOR i CHECK LOOP */
				}
			}
		
			for (n = 0; align2D_pass4[m][n] != '\0'; n++) { 
				keep_checking = 1;			/* THIS FLAG HANDLES THE CONTINUED NEED TO CHECK FOR INTRA-TR REPEATS   */
				imperfect_TR = 0;			/* THIS FLAG IS TURNED ON (SET TO ONE) WHEN TR W/ TRANSITION MISMATCHES IS FOUND */

				if (n == 0 && isalpha(align2D_pass4[m][0])) {
					x = 0;
					x_history[0] = x;
				}

				/* MOVE WINDOW PAST INITIAL BLANKS */
				while (align2D_pass4[m][n] == blnk) {
					cik_align2D[m+cik_row][n-x] = blnk;
					n++;
				}

				/* CHECK FOR NEED TO ADD ADDITIONAL BLANKS & ADJUST x */
				if (n <= scrimmage_line) {
					for (i = n-x; i < scrimmage_line; i++)
						cik_align2D[m+cik_row][i] = blnk;
					if (n > 1) {
						x = x_history[n-1];
					}
					else
						x = x_history[0];
				}

				/* CHECK FOR & DEAL WITH LONG HOMOPOLYMER RUN WRAPS (1ST ONE OR SUBSEQUENT ONES) */
				if (first_mwrap && n == first_mwrap_start) {		
					for (n = first_mwrap_start; (letr=align2D_pass4[m][n]) != '\0'; n++) {
						cik_align2D[m+cik_row][n-x] = letr;
						if (isalpha(letr)) {
							x_history[n] = x;					/* x_history WRITE-IN FOR MWRAP DONE ONCE */
						}
					}
					cik_align2D[m+cik_row][n-x] = '\0';	
					first_mwrap = 0;						/* RESET first_mwrap FLAG */
					break;									/* BREAK OUT OF FOR n LOOP */
				}
				else if (last_mwrap) {
					if (align2D_pass4[m][n] == kopt_Q_left && align2D_pass4[m][n+cik_mwrap+1] == kopt_R_rght) {		
						cik_align2D[m+cik_row][n-x] = kopt_Q_left;
						n++;
						while ( (letr=align2D_pass4[m][n]) != kopt_R_rght) {
							cik_align2D[m+cik_row][n-x] = letr;
							n++;
						}
						cik_align2D[m+cik_row][n-x] = letr;			/* WILL WRITE kopt_R_rght */
						for (i = cik_mwrap; i > 0; i--)				/* FAILSAFE:			*/
							cik_align2D[m+cik_row][n-x+i] = '\0';	/* OVERWRITE W/ NULLS	*/
						break;										/* BREAK OUT OF n LOOP	*/
					}
					else if (align2D_pass4[m][n] == kopt_Q_left && align2D_pass4[m][n+cik_mwrap+1] != kopt_R_rght) {
						last_mwrap = 0;
						for (i = 0; i <= cik_mwrap; i++) {
							cik_align2D[m+cik_row][n-x+i] = align2D_pass4[m][n+i];
						}
						letr = align2D_pass4[m][n+1];
						n = n + cik_mwrap;
						while (align2D_pass4[m][n] == letr) {
							cik_align2D[m+cik_row][n-x] = letr;
							x_history[n] = x;
							n++;
						}
					}
				}

				/* CHECK FOR & DEAL WITH LINE ENDS TOO SHORT TO HARBOR TR OF SIZE k */
				if ( isalpha(align2D_pass4[m][n+2*k-1]) == 0 ) {	/* TRUE IF WINDOW < 2x k-MER, WRITE REST OF LINE TO cik_align2D */
					for (i = n; align2D_pass4[m][i] != '\0'; i++) {
						cik_align2D[m+cik_row][i-x] = align2D_pass4[m][i];
						if (isalpha(align2D_pass4[m][i])) {
							x_history[i] = x;				/* x_history WRITE-IN FOR LINE ENDS TOO SHORT FOR TR */
						}
					}
					break;		/* BREAK OUT OF FOR n LOOP */
				}

				/* CHECK FOR TR OF SIZE k-MER */
				if (keep_checking) {
					if (k>1)
						homopolyflag = 1;	/* SET HOMOPOLYFLAG STATUS TO UNKNOWN/NEED TO CHECK */

					for (l = 0; l < k; l++) {
						/* BREAK EARLY IF WOULD-BE TRANSVERSION */
						if ( ((letr =align2D_pass4[m][n+l  ])=='A' || letr =='G') &&
							 ((letr2=align2D_pass4[m][n+l+k])=='C' || letr2=='T') ) {
							keep_checking = 0;
							break;
						}
						else if ( ((letr =align2D_pass4[m][n+l  ])=='C' || letr =='T') &&
							      ((letr2=align2D_pass4[m][n+l+k])=='A' || letr2=='G') ) {
							keep_checking = 0;
							break;
						}
						else {
							/* CHECK TO SEE IF THERE ARE n's */
							if (nuctype && (align2D_pass4[m][n+l]=='n' || align2D_pass4[m][n+l+k]=='n')) {
								keep_checking = 0;
								break;
							}
	
							/* CHECK TO SEE IF HOMOPOLYMER RUN CAN BE EXCLUDED */
							if (homopolyflag && l > 0) {	/* IF l > 0, THEN homopolyflag IS 1 */
								if ((align2D_pass4[m][n + l] != align2D_pass4[m][n+k+l  ]) ||
								    (align2D_pass4[m][n + l] != align2D_pass4[m][n + l-1]) ||
								    (align2D_pass4[m][n+k+l] != align2D_pass4[m][n+k+l-1]) ) {
									homopolyflag = 0;
								}
							}
	
							/* CHECK TO SEE IF THIS IS REGION OF AN IMPERFECT REPEAT W/ TRANSITIONS */
							if (nuctransit) {
								if (n == scrimmage_line || col_isclear(align2D_pass4,n,m,1) < 0) {
									y = 0;		/* RESET y VAR. B/C NO LONGER NEED TO ADJUST CONSENSUS COORDINATES */
								}

								if (     (letr=align2D_pass4[MAXROW][n-x+y+l  ]) == 'R' || letr == 'Y') {
									keep_checking = 0;
									if (k > piso_nuctransit) {
										check_imperf = 1;
									}
								}
								else if ((letr=align2D_pass4[MAXROW][n-x+y+l+k]) == 'R' || letr == 'Y') {
									keep_checking = 0;
									if (k > piso_nuctransit) {
										check_imperf = 1;
									}
								}
							}
	
							/* CHECK TO SEE IF THERE ARE MISMATCHES */
							if (keep_checking && align2D_pass4[m][n+l] != align2D_pass4[m][n+k+l]) {
								keep_checking = 0;
								if (nuctransit && k > piso_nuctransit)
									check_imperf = 1;	/* THUS WILL CHECK EVEN IF NOT ANNOTATED WITH R's & Y's BEFORE */
							}
						}
					} /* END OF FOR l SCAN LOOPS */

					if (homopolyflag && keep_checking) {		/* IF HOMOPOLYFLAG=1 WAS NOT SET TO ZERO */
						keep_checking = 0;
					}
					if (homopolyflag && check_imperf) {			/* TURN THIS OFF AS WELL */
						check_imperf = 0; 
					}
						homopolyflag = 0;		/* RESET */
				} 
				
				/* CHECK FOR TRANSITIONS IF DNA AND NEED REMAINS */
				if (check_imperf) { 
					/* LET MISMATCHES PASS THROUGH WITHOUT BREAKING IF ANNOTATED AS TRANSITIONS */
					sum4score = 0;
					for (l = 0; l < k; l++) {
						/* BREAK EARLY IF WOULD-BE TRANSVERSION */
						letr =align2D_pass4[m][n+l  ];
						letr2=align2D_pass4[m][n+l+k];
						if (isalpha(letr ) && !isupper(letr )) {
							letr = letr  - 32;	/* MAKE UPPER-CASE */
						}
						if (isalpha(letr2) && !isupper(letr2)) {
							letr2= letr2 - 32;	/* MAKE UPPER-CASE */
						}
						if ( (letr=='A' || letr=='G') && (letr2=='C' || letr2=='T') ) {
							break;
						}
						else if ( (letr=='C' || letr=='T') && (letr2=='A' || letr2=='G') ) {
							break;
						}

						if (align2D_pass4[m][n+l] != align2D_pass4[m][n+k+l]) { 
							if ( (letr =align2D_pass4[MAXROW][n-x+y + l]) != 'R' && letr  != 'Y' &&
								 (letr2=align2D_pass4[MAXROW][n-x+y+k+l]) != 'R' && letr2 != 'Y' ) {
								break;
							}

							/* MAKE SURE A MISMATCH IS NOT BEING GIVEN A PASS WHILE NOT CONFORMING TO TRANSITION TYPE */
							if ( (letr=align2D_pass4[MAXROW][n-x+y+l]) == 'R' && 
								( (letr2=align2D_pass4[m][n+l+k]) == 'C' || letr2 == 'T') ) {
								break;
							}
							else if (letr == 'Y' && 
								( (letr2=align2D_pass4[m][n+l+k]) == 'A' || letr2 == 'G')) {
								break;
							}

							if ( (letr=align2D_pass4[MAXROW][n-x+y+l+k]) == 'R' && 
								( (letr2=align2D_pass4[m][n+l]) == 'C' || letr2 == 'T') ) {
								break;
							}
							else if (letr == 'Y' && 
								( (letr2=align2D_pass4[m][n+l]) == 'A' || letr2 == 'G')) {
								break;
							}

							sum4score = sum4score + match;		/* JUSTIFICATION: MATCHES TRANSITION CALL */
						}
						/* IF LETTER AT n+l EQUAL TO LETTER AT n+l+k except 'n' */
						else if ((align2D_pass4[m][n+l] != 'n') && (align2D_pass4[m][n+l+k] != 'n'))	
							sum4score = sum4score + match;
					} /* END OF FOR l SCAN LOOPS */

					if (l == k && ((i=((100*sum4score)/(k*match))) >= DTHR_lookie[k]-1)) {	/*** INSPECT LATER: >= ***/
						imperfect_TR = 1;
					}
					else {
						imperfect_TR = 0;	/* IS REDUNDANT, BUT IS HERE AS A CODING FAIL-SAFE & FOR CLARITY */
					}						/* ALSO SET TO ZERO AT TOP OF n LOOP */

					check_imperf = 0;			/* RESET check_imperf HERE */
				} /* END OF IF check_imperf */

				/* BREAK CHECKING IF IT WILL PULL IN MISMATCHES IN REPEAT OR TO RIGHT OF REPEAT */
				if (nuctransit && (imperfect_TR || keep_checking) && n > scrimmage_line) { 
					for (l = 0; l < 2*k; l++) {
						letr2= align2D_pass4[m        ][n  +k+l];
						letr = align2D_pass4[MAXROW][n-x+k+l];
						if (imperfect_TR) {
							if (isalpha(align2D_pass4[m][n+l]) && (i=col_isclear(align2D_pass4,n+l,m,-1)) > -1 && 
								align2D_pass4[i][n+l] != letr2 && letr != 'R' && letr != 'Y') {
								keep_checking = imperfect_TR = 0;
								break;
							}
							if (letr == 'R' && (letr2 == 'C' || letr2 == 'T')) {
								imperfect_TR = 0;
								break;
							}
							else if (letr == 'Y' && (letr2 == 'A' || letr2 == 'G')) {
								imperfect_TR = 0;
								break;
							}
						}
					}
				}

				if (keep_checking && n > scrimmage_line) { 
					if ((l=col_isclear(cik_align2D,n-x+k,m,-1)) > -1 
						&& col_isclear(align2D_pass4,n+k,m,1)< 0) {

						/* CHECK IF WILL PULL IN ADJACENT MISMATCHES AFTER RUN OF REPEATS */
					    r = 1; 
						i = k;  /* VAR i SET TO k ONLY TO ENTER WHILE LOOP */
						while (i==k) {
							for (i = 0; i < k; i++) {
						    	if (align2D_pass4[m][n+i] != (letr2=align2D_pass4[m][n+(r+1)*k+i]) && isalpha(letr2))
									break;
					        }    
					        if (i == k)
								r++;        /* INCREMENT NUMBER OF REPEATS */
					    }    

						if (nuctransit) {
							if ((letr=cik_align2D[l][n-x+k+i]) != (letr2=align2D_pass4[m][n+(r+1)*k+i]) &&
								 letr!='R' && letr!='Y' && 
								 isalpha(letr) && isalpha(letr2) && (letr3=align2D_pass4[MAXROW][n-x+(r+1)*k+i])!='R' && letr3!='Y') {
						        keep_checking = 0; 
						    }    
						}
						else {
							if ((letr=cik_align2D[l][n-x+k+i]) != (letr2=align2D_pass4[m][n+(r+1)*k+i]) &&
								 isalpha(letr) && isalpha(letr2)) {
						        keep_checking = 0; 
						    }    
						}
			    	}
				}

				if (imperfect_TR && n > scrimmage_line && 
					isupper(align2D_pass4[m][n]) && align2D_pass4[m+1][n] == blnk) { 
					for (l = 0; l < k; l++) {
						/* CHECK MISMATCHES FROM PUSHING BOTTOM ROW TO LEFT OF REPEATS AFTER SLIP */
						if ((i=col_isclear(align2D_pass4,n+l,m,1)) > -1 &&
							(letr=align2D_pass4[i][n+l]) != align2D_pass4[MAXROW][n-x-k+l] && 
						 	                     letr != 32+align2D_pass4[MAXROW][n-x-k+l]) {
							imperfect_TR = 0;
							break;
						}
					}
				}

				if (keep_checking || imperfect_TR) {
					for (l = 0; l < k; l++) {
						cik_align2D[m+cik_row  ][n-x+l] = align2D_pass4[m][n+l  ];		
						cik_align2D[m+cik_row+1][n-x+l] = align2D_pass4[m][n+l+k];
						x_history[n+l] = x;					/* x_history WRITE-IN FOR NEW TR COLS */
					}
						cik_align2D[m+cik_row  ][n-x+k  ] = '/';
						cik_align2D[m+cik_row  ][n-x+k+1] = '\0';

 					for (i = 0; i < n-x; i++) {
						if (isalpha(cik_align2D[m+cik_row+1][i]) == 0)
							cik_align2D[m+cik_row+1][i] = blnk;
					}

					/* SCOOCH CONSENSUS ROW IF MINDING TRANSITIONS AND IF BOTTOM IS CLEAR (IS SAFE) */	
					if (nuctransit) {
						if (col_isclear(align2D_pass4,n,m,1) < 0) { 
							for (i = n-x; i < n-x+k; i++) {
								if ((letr=align2D_pass4[MAXROW][i]) != 'R' && letr != 'Y')
									align2D_pass4[MAXROW][i] = align2D_pass4[MAXROW][i+k];
							}
							for (i = n-x+k; i+k < koptions[1][32]; i++) {
								align2D_pass4[MAXROW][i] = align2D_pass4[MAXROW][i+k];
							}
							align2D_pass4[MAXROW][i] = '\0';
						}
						else if (n >= scrimmage_line) {
							y = y + k;	/* TO KEEP TRACK OF UNSHIFTED CONSENSUS ROW */
						}
					}
					x = x + k;			/* FUTURE SPACING TO BE SUBTRACTED B/C k-MER TUCKED UNDER 1st UNIT */
					scrimmage_line = n;
					x_history[n] = x;

					++cik_row;	 	
					n = n + k - 1;		/* ADVANCE ADJUSTMENT. NOTE UPCOMING n++ IN FOR n LOOP */

				}   /* END OF TR ASSIGN LOOPS */
				else {
					cik_align2D[m+cik_row][n-x] = align2D_pass4[m][n];
					x_history[n] = x;
				}

			}   /* END OF FOR n LOOPS */ 
		}   /* END OF FOR m LOOPS */

		if (cik_row > 0) {
			mha_writeback(cik_align2D, align2D_pass4, koptions); 
			if (koptions[0][57]) {			/* opt_v VERBOSITY */
				printf("\nk = %d:", k);
			}

			if (koptions[0][57] && k > 1)	/* opt_v VERBOSITY, k=1 WILL PRINT FROM MAIN */
				print_2Dseq(align2D_pass4, koptions[1][32], &letr, koptions);
		}

		koptions[0][4] = koptions[0][4] + cik_row;			/* STORE ROWS ADDED */

		if (k > 1)	/* NOT NEEDED AFTER k EQUALS ONE */
			clear_2D_ar(cik_align2D);	/* REMINDER: DOES NOT CLEAR MAXROW */

	} /* END OF FOR k LOOPS */ 

	mha_UPPERback(cik_align2D, align2D_pass4, koptions); /* THIS ALSO SAVES 2D-WIDTH in options[1][32] */

	i = koptions[1][18];
	koptions[1][i] = koptions[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS i WIDTH HISTORY */
}


/*** FUNCTION 03 **** CINCH sINGLE LINES ************************************************************/
short unsigned int cinch_s(char align2D_pass6[][MAXROW], long int soptions[][62]) 
{
int cis_row=0, i=0, k=0, l=0, m=0, n=0, next = -1, x=0; 
int first_mwrap_start=0, last_mwrap=0, nuctype;
unsigned short int first_mwrap=0, keep_checking=1, nuctransit=0;
char letr;
char blnk        = (char) soptions[1][11];		/* opt_B blank character */
int cis_mwrap    =        soptions[1][22];		/* opt_M long_homopolymer_run */
char sopt_Q_left = (char) soptions[1][26];		/* LHS character delimiter for homopolymer Run */
char sopt_R_rght = (char) soptions[1][27];		/* RHS character delimiter for homopolymer Run */
char               cis_align2D[MAXROW+1][MAXROW] = {{0}};
void               clear_2D_ar(char wipe_align2D[][MAXROW]);
int 			   col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown); 
int                count_wrap_blocks(int lcl_width, int lcl_opt_w);
void               line_end(int type, int c, long int lend_options[][62], int lcl_width);
void 			   mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
void               mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
void               mha_UPPERback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
int                next_transit(char lcl_align2D[][MAXROW], int last_transit);
short unsigned int print_2Dseq(char align2D_print[][MAXROW], int print_lenseq2D, char *printptr_Seq_name, long int poptions[][62]);

	nuctype = soptions[1][13];		/* EQUALS ONE IF DNA, TWO IF RNA */

	if (nuctype == 1)				/* IF DNA */
		nuctransit = 1;

	if (nuctransit) { 
		if ((next=next_transit(align2D_pass6, -1)) == -1) {
			nuctransit = 0;	/* THERE IS NO NEED TO KEEP LOOKING FOR TRANSITION MUTATIONS IN STRING */
		}
	} 

	/* START AT BIGGEST k-MER POSSIBLE AT 2x */
	for (k = soptions[1][32]/2; k > 0; k--) {
		cis_row = x = 0;

		clear_2D_ar(cis_align2D);

		for (m = 0; align2D_pass6[m][0] != '\0'; m++) {
			/* CHECK LINE AHEAD OF TIME FOR FIRST MONO-RUN TERMINATOR */
			for (i = 0; (letr=align2D_pass6[m][i]) != '\0'; i++) {	
				/* NOTE: APPARENT REDUNDANCY BELOW ALLOWS CHANGING L/R RUN DELIMITERS TO BE EQUAL */
				if (letr == sopt_R_rght && align2D_pass6[m][i-1] != blnk) {
					first_mwrap = 1;			/* TURN ON BIT FLAG FOR HANDLING FIRST MONO-RUN */
					last_mwrap = 1;				/* TURN ON BIT FLAG FOR HANDLING TO LAST MONO-RUN */
					first_mwrap_start = i-10;	/* SAVE POSITION OF START OF MWRAP */
					break;						/* BREAK OUT OF FOR i CHECK LOOP */
				}
			}
		
			for (n = 0; align2D_pass6[m][n] != '\0'; n++) { 

				/* MOVE WINDOW PAST INITIAL BLANKS */
				while (align2D_pass6[m][n] == blnk) {
					cis_align2D[m+cis_row][n-x] = blnk;
					n++;
				}
				keep_checking = 1;			/* THIS FLAG HANDLES THE POTENTIAL NEED TO CHECK FOR INTRA-TR REPEATS   */

				/* CHECK FOR & DEAL WITH LONG HOMOPOLYMER RUN WRAPS (1ST ONE OR SUBSEQUENT ONES) */
				if (first_mwrap && n == first_mwrap_start) {		
					for (n = first_mwrap_start; (letr=align2D_pass6[m][n]) != '\0'; n++) {
						cis_align2D[m+cis_row][n-x] = letr;
					}
					cis_align2D[m+cis_row][n-x] = '\0';	
					first_mwrap = 0;						/* RESET first_mwrap FLAG */
					break;									/* BREAK OUT OF FOR n LOOP */
				}
				else if (last_mwrap) {
					if (align2D_pass6[m][n] == sopt_Q_left && align2D_pass6[m][n+cis_mwrap+1] == sopt_R_rght) {		
						cis_align2D[m+cis_row][n-x] = sopt_Q_left;
						n++;
						while ( (letr=align2D_pass6[m][n]) != sopt_R_rght) {
							cis_align2D[m+cis_row][n-x] = letr;
							n++;
						}
						cis_align2D[m+cis_row][n-x] = letr;			/* WILL WRITE sopt_R_rght */
						for (i = cis_mwrap; i > 0; i--)				/* FAILSAFE:			*/
							cis_align2D[m+cis_row][n-x+i] = '\0';	/* OVERWRITE W/ NULLS	*/
						break;										/* BREAK OUT OF n LOOP	*/
					}
					else if (align2D_pass6[m][n] == sopt_Q_left && align2D_pass6[m][n+cis_mwrap+1] != sopt_R_rght) {
						last_mwrap = 0;
						for (i = 0; i <= cis_mwrap; i++) {
							cis_align2D[m+cis_row][n-x+i] = align2D_pass6[m][n+i];
						}
						letr = align2D_pass6[m][n+1];
						n = n + cis_mwrap;
						while (align2D_pass6[m][n] == letr) {
							cis_align2D[m+cis_row][n-x] = letr;
							n++;
						}
					}
				}

				/* CHECK FOR & DEAL WITH LINE ENDS TOO SHORT TO HARBOR TR OF SIZE k */
				if ( isalpha(align2D_pass6[m][n+2*k-1]) == 0 ) {	/* TRUE IF WINDOW < 2x k-MER, WRITE REST OF LINE TO cis_align2D */
					for (i = n; align2D_pass6[m][i] != '\0'; i++) {
						cis_align2D[m+cis_row][i-x] = align2D_pass6[m][i];
					}
					break;		/* BREAK OUT OF FOR n LOOP */
				}

				/* CHECK FOR TR OF SIZE k-MER ********************************************************************/
				if (keep_checking) {
					for (l = 0; l < k; l++) {
						if (align2D_pass6[m][n+l] != align2D_pass6[m][n+l+k] || align2D_pass6[m][n+l]=='n') {
							keep_checking = 0;  
							break;			  /* BREAK OUT OF FOR l LOOP */ 
						}
					}   /* END OF FOR l SCAN LOOPS */
				} /***********************************************************************************************/

				/* CHECK FOR POTENTIAL FOR RHS OVERSLIP INTO PRIOR TR */ 
				if (k > 1 && keep_checking && m > 0 && n > 0) {
					for (l = 1; l < k; l++) {
						if (align2D_pass6[m-1][n+k+l] == '/') {
							for (l = 0; l < k; l++) {
									if (align2D_pass6[m-1][n-1+l] != blnk) {
										keep_checking = 0;
										l = k; /* TO BREAK OUT OF OUTER FOR l LOOP */
										break; /* TO BREAK OUT OF INNER FOR l LOOP */
									}
							}
						}
					}
				}

				/* CHECK FOR EDGE HOMOPOLYMER EFFECTS: DON'T cinch_s UNTIL k=1 */
				if (k>1 && keep_checking) {
					if ((letr=align2D_pass6[ m ][n]) == align2D_pass6[m][n+1] && 
							  align2D_pass6[m+1][n] == blnk ) {
							cis_align2D[m+cis_row][n-x] = letr + 32;	/* MAKE LOWER-CASE */	
							keep_checking = 0;
						}
				}

				/* CHECK FOR CONFLICT ABOVE AND BELOW */
				for (i=0; i < m; i++) {
					if (isalpha(align2D_pass6[i][n+k])) {
						keep_checking = 0;
						break;
					}
				}
				for (i=m+1; align2D_pass6[i][0] != '\0'; i++) {
					if (isalpha(align2D_pass6[i][n+k-1])) {
						keep_checking = 0;
						break;
					}
				}
				
				/* BREAK CHECKING IF IT WILL PULL IN ADJACENT MISMATCH */
				if (keep_checking && isalpha(align2D_pass6[m][n+2*k])) {
					for (l = 0; l < k; l++) {

						if ((letr=align2D_pass6[MAXROW][n-x+k+l]) == 'R' || letr == 'Y') {	
							if (letr == 'R' && (align2D_pass6[m][n+2*k+l] != 'A' || align2D_pass6[m][n+2*k+l] != 'G') ) {
								keep_checking = 0;
							}
							else if (letr == 'Y' && (align2D_pass6[m][n+2*k+l] != 'C' || align2D_pass6[m][n+2*k+l] != 'T') ) {
								keep_checking = 0;
							}
						}
						else if (col_isclear(align2D_pass6,n+k+l,m,-1) > -1 && align2D_pass6[MAXROW][n-x+k+l] != align2D_pass6[m][n+2*k+l]) {
							keep_checking = 0;
						}
					}
				}
						
				if (keep_checking) {
					for (l = 0; l < k; l++) {
						cis_align2D[m+cis_row  ][n-x+l] = align2D_pass6[m][n+l  ];		
						cis_align2D[m+cis_row+1][n-x+l] = align2D_pass6[m][n+l+k];
					}
						cis_align2D[m+cis_row  ][n-x+k  ] = '/';
						cis_align2D[m+cis_row  ][n-x+k+1] = '\0';
					for (i = 0; i < n-x; i++) 
						if (isalpha(cis_align2D[m+cis_row+1][i]) == 0) {
							cis_align2D[m+cis_row+1][i] = blnk;
						}
						else
							cis_align2D[m+cis_row+1][i] = 'x';		/* FLAG FOR NON-PERMITTED CHARACTER */
					for (i = n+2*k; (letr=align2D_pass6[m][i]) != '\0'; i++)
						cis_align2D[m+cis_row+1][i-x-k] = letr;

					/* SCOOCH CONSENSUS ROW IF MINDING TRANSITIONS AND IF BOTTOM IS CLEAR (IS SAFE) */					
					if (nuctransit && col_isclear(align2D_pass6,n,m,1) < 0) { 
						for (i = n-x+k; i+k < soptions[1][32]; i++) {
							align2D_pass6[MAXROW][i] = align2D_pass6[MAXROW][i+k];
						}
						align2D_pass6[MAXROW][i] = '\0';
					}

					++cis_row;
					x = x + k;		
					n = n + 2*k - 1;		/* ADVANCE ADJUSTMENT. NOTE UPCOMING n++ IN FOR n LOOP */

				}   /* END OF TR ASSIGN LOOPS */
				else {
					cis_align2D[m+cis_row][n-x] = align2D_pass6[m][n];
				}

			}   /* END OF FOR n LOOPS */ 

		}   /* END OF FOR m LOOPS */
		
		if (cis_row > 0) {
			mha_writeback(cis_align2D, align2D_pass6, soptions);
			mha_writecons(align2D_pass6, cis_align2D, soptions);
		}
	
	soptions[0][6] = soptions[0][6] + cis_row;			/* STORE ROWS ADDED */

	} /* END OF FOR k LOOPS */ 

	i = soptions[1][18];
	soptions[1][i] = soptions[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS i WIDTH HISTORY */

	if (soptions[0][6]) {
		mha_UPPERback(cis_align2D, align2D_pass6, soptions); /* THIS ALSO SAVES 2D-WIDTH in options[1][32] */
		mha_writecons(align2D_pass6, cis_align2D, soptions);
		return(1);
	}
	else
		return(0);
}


/*** FUNCTION 04 **********************************************************/
unsigned int cinch_d(char align2D_pass7[][MAXROW], char *dptr_Seq_name, long int doptions[][62], short unsigned int cinch_d_opt)
{
int cid_mrow=0, cid_ncol=0, h=0, i=0, j=0, k=WIDTH, l=0, m=0, n=0, num=0, w=0, x=0;
int tot_repeats=0, uniq_TRs=0;
int cid_new2D_width = doptions[1][32]; 
unsigned short int nuctype=0, TR_check=0, first_write=1, mono_flag=1;		/* CHECK MONO IN ORDER TO KNOW TO SKIP IT */
unsigned short int nuctransit=0;						/* BIT FLAG FOR HANDLING NUCLEOTIDE TRANSITIONS SILENTLY (IGNORING) */
unsigned short int imperfect_TR=0;
int num_transits=0;
char letr = 'B';	/* Begin the Beguine */
char ltr2 = 'Z';	/* BZBZBZ				  */
char dopt_R_rght = (char) doptions[1][27];					/* RHS character delimiter for homopolymer Run */
char blnk        = (char) doptions[1][11];					/* opt_B blank character */
char cid_align2D[MAXROW+1][MAXROW];
unsigned int       consensus_2D(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width);
void               mha_head(char lcl_seq_name, int lcl_width, long int lcl_options[][62]);
void 			   mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
void               mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
short unsigned int print_2Dseq(char align2D_print[][MAXROW], int print_lenseq2D, char *printptr_Seq_name, long int poptions[][62]);
short int pushdown(char pusharray[][MAXROW], int push_m, int push_n, long int push_options[0][62]);
short int tucksense(char tuckarray[][MAXROW], long int tuck_options[0][62]);
	
	mha_writeback(align2D_pass7, cid_align2D, doptions); 
	nuctype = doptions[1][13];		/* EQUALS ONE IF DNA, TWO IF RNA */

	if (nuctype == 1) {		/* IF DNA */
		nuctransit = 1;
	}

	/* START AT BIGGEST k-MER POSSIBLE AT 2x */
	for (k = doptions[1][32]/2; k > 0; k--) {
		for (n=0; n < cid_new2D_width-2*k; n++) {	

			mono_flag = 1;			/* MONOMER RUN FLAG IS SET TO 0, WHEN NO LONGER POSSIBLE (ANY n != n+1) */
	
			if (TR_check == 0) 		/* RE-SET COUNTER FOR NUM */
				num = 0;

			if (nuctransit)			/* RE-SET COUNTER FOR NUMBER OF TRANSITIONS */
				num_transits = imperfect_TR = 0;

			for (l=0; l < k; l++) {
				if (nuctransit && k>2 && (letr=align2D_pass7[MAXROW][n+l]) != (ltr2=align2D_pass7[MAXROW][n+k+l])) {
					if (     (letr == 'A' || letr == 'G') && (ltr2 == 'C' || ltr2 == 'T'))
						break;
					else if ((letr == 'C' || letr == 'T') && (ltr2 == 'A' || ltr2 == 'G'))
						break;
					else if ((letr == 'R') && (ltr2 == 'C' || ltr2 == 'T' || ltr2 == 'Y'))
						break;
					else if ((letr == 'Y') && (ltr2 == 'A' || ltr2 == 'G' || ltr2 == 'R'))
						break;
					else if ((ltr2 == 'R') && (letr == 'C' || letr == 'T' || letr == 'Y'))
						break;
					else if ((ltr2 == 'Y') && (letr == 'A' || letr == 'G' || letr == 'R'))
						break;
					else 
						++num_transits;
				}
				else if (   (letr=align2D_pass7[MAXROW][n+l]) != (ltr2=align2D_pass7[MAXROW][n+k+l])) 
					break; 		/* BREAK OUT OF FOR l LOOP */

				if (letr == 'n' || ltr2 == 'n') {
					break;
				}

				if (l+1<k && letr!=align2D_pass7[MAXROW][n+l+1]) /* mono_flag IS FOR FAILING EARLY */
					mono_flag = 0;		/* CAN NO LONGER BE A HOMOPOLYMER RUN OF FOR THIS k-MER */

			} /* END OF FOR l LOOP CHECK OF k-MER TR (2x) */

			if (l==k && mono_flag==0) {
				++tot_repeats;
				num = 2;		/* THIS KEEPS COUNT OF HOW MANY REPEATS. WITH ONE RE-PEAT COUNTED, THERE ARE TWO */
				TR_check = 1;

				if (num_transits) {		/* ONLY POSSIBLE IF nuctransit MODE 2 */
					if (num_transits < 2)
						imperfect_TR = 1;
					else 
						TR_check = num_transits = imperfect_TR = 0;
				}

			}
			else 
				TR_check = 0;

			/* CHECK FOR COMPLEX cinch_d REPEATS THAT SHOULD NOT BE COUNTED/WRITTEN. */
			/*  THESE HAVE LETTERS IN NEXT ROW UNDERNEATH FIRST UNIT.                */
			if (TR_check) {
				m = 0;
				while (isalpha(align2D_pass7[m][n+k]) == 0) {
					m++;
				}
				letr=align2D_pass7[m][n+k];
				if (align2D_pass7[MAXROW][n]==letr && align2D_pass7[m][n+k+1]=='/' && 
					align2D_pass7[m+1][n]==blnk && isalpha(align2D_pass7[m+1][n+1])) {
						align2D_pass7[m  ][n+k  ]='/';
						align2D_pass7[m  ][n+k+1]='\0';
						align2D_pass7[m+1][n    ]=letr;
				}
				for (w=1; m+w < MAXROW && TR_check != 0; w++) {
					for (x=0; x < n+k; x++) {
						if (x==0 && align2D_pass7[m+w][0]=='\0') {
							w = MAXROW; 	/* TO BREAK FOR w LOOP */
							break;			/* TO BREAK FOR x LOOP */
						}
						else if (isalpha(align2D_pass7[m+w][x])) {
							--tot_repeats;	/* DECREMENT COUNTER  */
							TR_check = 0;	/* RESET BACK TO ZERO AND BREAK OUT OF FOR w LOOP */
							break; 			/* BREAK OUT OF FOR x LOOP     */
						}
					}
				}
			}
	
			while (TR_check) {
				++uniq_TRs;

				/* THIS PART JUST COUNTS REPEATS ADDITIONAL REPEATS, VAR num STARTS AT 2 */
				for (l=0; l < k; l++) {
					if (align2D_pass7[MAXROW][n+l] != align2D_pass7[MAXROW][n+num*k+l]) {
						TR_check = 0;	/* WILL BREAK OUT OF WHILE TR_check LOOP 		*/
						break; 			/* BREAK OUT OF FOR l LOOP */
					}
					if (l == k-1) {
						++num;		/* INCREMENT TR COUNT */
						l = -1;		/* RESET TO CHECK NEXT POSSIBLE UNIT REPEAT, l SET TO -1 B/C OF UPCOMING l++ */
					}
				} /* END OF FOR l LOOP */

				if (cinch_d_opt) {	/* cinch_d ENGINE **********************************************************************/
					if (first_write) {
						m = 0;
						while (isalpha(align2D_pass7[m][n+k]) == 0) {
							m++;
						}
						if (doptions[1][57])
							printf("\nWorking on next consensus TR: %dx %d-mer at consensus position %d; 2nd unit begins on row %d.", num, k, n+1, m+1);

						if (imperfect_TR == 1) {
							for (l=0; l < k; l++) {
								letr=align2D_pass7[MAXROW][n+l];
								if (letr != align2D_pass7[MAXROW][n+k+l]) {
									if (letr == 'A' || letr == 'G') 
										align2D_pass7[MAXROW][n+l] = 'R';
									else if (letr == 'C' || letr == 'T') 
										align2D_pass7[MAXROW][n+l] = 'Y';
								}
							} 
						}

						cid_align2D[m][n+k  ] = '/';
						cid_align2D[m][n+k+1] = '\0';
						first_write = 0;	/* TURN OFF NEED TO WRITE REMAINING PART OF 2-D ALIGNMENT */

						cid_ncol = k;
						cid_mrow = 1;

						/* DEAL WITH LOOSE SLIP CONNECTIONS PRODUCED BY FUDGE-CYCLELIZING */
						for (j = 0; j < n+k; j++) {
							if (cid_align2D[m][j] != blnk) {
								break;
							}
						}
						if (j == n+k) {
							if (cid_align2D[m-1][n] == '/')
								cid_mrow = -1;
							else 
								cid_mrow = 0;
						}

						for (i = m; i < MAXROW; i++) {
							for (j = n+k; j < doptions[1][32]+1; j++) {
								letr = cid_align2D[i+cid_mrow][j-cid_ncol] = align2D_pass7[i][j];
								if (letr=='/' || letr==dopt_R_rght) 
									cid_align2D[i+cid_mrow][j-cid_ncol+1] = '\0';
								else if (letr=='>') {
									for (h=0; h<n; h++)
										cid_align2D[i+cid_mrow][h] = blnk;
									cid_align2D[i+cid_mrow+1][0] = '\0';
									i = MAXROW;
									break;
								}
							}
						}

						if (nuctransit) {
							/* SCOOCH CONSENSUS ROW TO THE LEFT TO REFLECT CINCHED WIDTH */
							for (i = n+k; align2D_pass7[MAXROW][i+k] != '\0'; i++) {
								align2D_pass7[MAXROW][i] = align2D_pass7[MAXROW][i+k];
							}
							align2D_pass7[MAXROW][i] = '\0';
						} 

						if (letr == '>' && j-cid_ncol-1 < cid_new2D_width) {
							doptions[1][32] = j-cid_ncol-1;
							for (i = 0; (letr=align2D_pass7[MAXROW][i]) != '\0'; i++)
								cid_align2D[MAXROW][i] = letr;

							mha_writeback(cid_align2D, align2D_pass7, doptions);
							mha_writecons(align2D_pass7, cid_align2D, doptions);
						}
					} /* END OF IF first_write EQUALS ONE */
				} /*********************************************************************************************/
				else if (doptions[1][57]) {
					if (imperfect_TR)
						printf("  %4d. i-TR: %3dx %d-mer at consensus position %3d with %d transition(s).\n", uniq_TRs, num, k, n+1, num_transits);
					else if (nuctransit)
						printf("  %4d. p-TR: %3dx %d-mer at consensus position %3d.\n", uniq_TRs, num, k, n+1);
					else
						printf("  %4d. TR: %3dx %d-mer at consensus position %3d.\n", uniq_TRs, num, k, n+1);
				}

			} /* END OF WHILE TR_check EQUALS ONE LOOP */
		} /* END OF FOR n LOOP */
	} /* END OF FOR k LOOP */

	if (cinch_d_opt == 0 && tot_repeats == 0) {
		i = doptions[1][18];
		doptions[1][i] = doptions[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS WIDTH HISTORY */
		printf("\n");
	}
	else if (cinch_d_opt) {
		i = doptions[1][18];
		doptions[1][i] = doptions[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS WIDTH HISTORY */
		if (cid_new2D_width == doptions[1][32]) {
			print_2Dseq(align2D_pass7, cid_new2D_width, dptr_Seq_name, doptions);
			return(0);
		}
		else if (tot_repeats > 1 && doptions[1][57]) {
			cid_new2D_width = doptions[1][32];
			print_2Dseq(align2D_pass7, cid_new2D_width, dptr_Seq_name, doptions);
		}
		else if (tot_repeats > 1 && doptions[0][20]) {
			cid_new2D_width = doptions[1][32];
			doptions[0][20] = 0;	/* TEMPORARY ASSIGNMENT TO PREVENT PRINTING OF CONSENSUS ROW */
			consensus_2D(align2D_pass7, doptions, 0, doptions[1][32]);
			doptions[0][20] = 1;	/* REASSIGN SETTING */
		}
		else if (tot_repeats > 1) {
			cid_new2D_width = doptions[1][32];
			consensus_2D(align2D_pass7, doptions, 0, doptions[1][32]);
		}
		else { 
			cid_new2D_width = doptions[1][32];
			print_2Dseq(align2D_pass7, cid_new2D_width, dptr_Seq_name, doptions);
		}
	}
	return(tot_repeats);
}


/*** FUNCTION 05 ***********************************************************************/
short unsigned int cleanseq(char *s, long int cloptions[][62], short unsigned int storebit) 
{
int c=0, i=0, dna=0, rna=0, na=0, ens=0, prot=0, noprot=0, x=0, length=0;
int slop = 12;	/* PERCENT SLOP ALLOWED WHILE MATCHING SEQTYPE PROFILE; SLOP ACCOMODATES NON-DIAGNOSTIC ALLOWED LETTERS */
char letr;
short unsigned int stringtype = 0;	/* ALPHA TYPE IS DEFAULT UNTIL FOUND OTHERWISE				*/
short unsigned int nacheck = 1;		/* BIT FLAG. CHANGE TO ZERO WHEN NUCLEIC ACID IS NOT POSSIBLE 	*/

	/* MAKE STRING ALL UPPERCASE */
		while ((letr=s[c]) != '\0') {
			if (letr >= 'a' && letr <= 'z') 
				s[c] = letr - 32;
		c++;
		}

	/* REMOVE NON-LETTERS */
		c = 0; x = 0;
		while ((letr=s[c]) != '\0') {
			if (letr < 'A' || letr > 'Z') {
				x++;	/* EXTRA CHARACTERS REMOVED */
				c++;	/* POSITION */
			}
			else {
				if (letr == 'A' || letr == 'C' || letr == 'G') {
					++na;
				}
				else if (letr == 'T') {
					++dna;
				}
				else if (letr == 'U') {
					++rna;
				}
				else if (letr == 'E' || letr == 'F' || letr == 'I' || letr == 'L' || letr == 'P' || letr == 'Q') {
					++prot;
				}
				else if (letr == 'J' || letr == 'O' || letr == 'X' || letr == 'Z') {
					++noprot;
				}
				else if (letr == 'N') {
					++ens;
				}

				s[c-x] = letr;
				c++;
			}
		}

	length = c-x;

	if (prot+noprot > 0)
		nacheck = 0;

	/* ASSIGN STRING TYPE: NON-ZERO VALUES FOR NUCLEIC ACIDS DNA AND RNA */
	if (nacheck) {	/* SUMMARY: U's NECESSARY FOR CALLING RNA, BUT T's NOT NECESSARY FOR CALLING DNA (DEFAULT NA) 	*/
	 	if (rna > 0 && (100*(na+rna+ens)) >= length*(100-slop)) 	/* 0.90 EQUIVALENT TO 18/20 bp BEING A|C|G|U 	*/
			stringtype = 2;											/* STRING IS LIKELY RNA 						*/
	 	else if (      (100*(na+dna+ens)) >= length*(100-slop)) 	/* 0.90 EQUIVALENT TO 18/20 bp BEING A|C|G|T 	*/
			stringtype = 1;											/* STRING IS LIKELY DNA 						*/			
	}

	if (stringtype==0 && noprot+rna==0 && prot > 0 && (100*(prot+na+dna+ens))/length >= 50-slop)	/* 50 B/C COUNTING ~1/2 AMINO ACIDS */
		stringtype = 3;

	/* OVER-WRITE TERMINAL PART OF ORIGINAL STRING WITH NULLS */
	for (i = c-x; i <= length; i++)
		s[i] = '\0';

	/* IF DNA (STRICT, NOT IUPAC FULL), ALL N'S AND OTHERS TO LOWERCASE 'n' AS AMBIGUOUS DNA */
	if (stringtype == 1 || stringtype == 2) {	/* IF DNA OR RNA */
		c = 0;
		while (s[c] != '\0') {
			if (s[c] == 'N') 
				s[c] = 'n';
			else if (s[c]!='A' && s[c]!='G' && s[c]!='C' && s[c]!='T' && s[c]!='U')
				s[c] = 'n';
			c++;
		}
	}

	if (storebit)
		cloptions[1][13] = stringtype;

	return(stringtype); 

	/*	SIGNATURES OF DIFFERENT ALPHABETS: 
		1	A - C - - - G - - - - - - N - - - - - T - - - - - - = DNA 
		2	A - C - - - G - - - - - - N - - - - - - U - - - - - = RNA
		0	A B C D - - G H - - K - M N - - - R S T - V W - Y - = IUPAC DNA
		0	A B C D - - G H - - K - M N - - - R S - U V W - Y - = IUPAC RNA
		0	A - C D E F G H I - K L M N - P Q R S T - V W - Y - = PROTEIN (20 amino acids)
		0	A - C D E F G H I - K L M N O P Q R S T U V W - Y - = PROTEIN (22 a.a. w/ pyrolysine=O and selenocysteine=U)
		0	- - - - E F - - I - - L - - - P Q - - - - - - - - - = PROTEIN-ONLY (THESE 6 amino acids are diagnostic)
		0	- - - - - - - - - J - - - - O - - - - - - - - X - Z = NOT CANONICAL PROTEIN, NOT-IUPAC 
	 	0	A B C D E F G H I J K L M N O P Q R S T U V W X Y Z	= ALPHABET
	*/

}


/*** FUNCTION 06 *************************************************************/
void clear_2D_ar(char wipe_align2D[][MAXROW])
{
int m=0, n=0;

	/* CLEAR UP TO BUT NOT INCLUDING CONSENSUS ROW AT MAXROW -1 */
	for (m=0; m < MAXROW; m++) {
		for (n=0; n < MAXROW; n++)
			wipe_align2D[m][n] = '\0';
	}
}


/*** FUNCTION 07 *************************************************************/
void clear_right(char swipe_align2D[][MAXROW], long int croptions[][62])
{
int m=0, n=0;
int width = croptions[1][0];	/* opt_W */
int height = croptions[1][17];
char letr;
char cropt_R_rght = croptions[1][27];

	/* CLEAR TO THE RIGHT OF TERMINATORS */
	for (m=0; m < height; m++) {
		n = 0;
		while ( (letr=swipe_align2D[m][n]) != '/' && letr != cropt_R_rght && letr != '>') {
			n++;
		}
		for (n = n+1; n < width; n++)
			swipe_align2D[m][n] = '\0';
		if (letr == '>') {
			for (m=m+1; m < height; m++) {
				for (n=0; n < width; n++) {
					swipe_align2D[m][n] = '\0';
				}
			}
			break;
		}
	}
}


/*** FUNCTION 08a **********************************************************/
unsigned int consensus_2D(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width)
{
int badsites=0, m=0, n=0, n_end, x=1;
int magiclen=400; 	/* SEARCH MAGIC TO FIND OTHER STOPGAPS */
int con_width = con_options[1][32];
short unsigned int nuctype = con_options[1][13];	/* FOR SEQ TYPE, DNA=1, RNA=2, OTHER (NON-NA)=0 */
short unsigned int nuctransit = 0;					/* BIT FLAG FOR HANDLING NUCLEAR TRANSITIONS */ 
short unsigned int plustransit=0;					/* BIT FLAG ADDENDUM FOR COUNTING BADSITES AT COL */
char blnk = con_options[1][11];
char letr=blnk, conletr=blnk;
int con_maxrows=26;
int col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown); 
int consensus_ar[26][MAXROW] = {{0}};	 	/* COL n=0 FOR BIT FLAG */
                                       		/* ROW m=0 FOR COUNTER */
											/* ROW m=1 FOR CONSENSUS */
											/* ROWS m>1 FOR VARIANTS STORAGE */
void line_end(int type, int c, long int lend_options[][62], int lcl_width);
char mha_base62(int num);
void warnhead(char l);

	n_end = n_start + n_width;

	if (nuctype == 1)	/* IF DNA */
		nuctransit = 1;

	/* FILL CONSENSUS ARRAY WITH FIRST LETTER IN EACH COL */
	for (n = 0; n <= con_width; n++) {
		for (m = 0; m < MAXROW; m++) {
			if ( (isalpha(letr=con_align2D[m][n]))) {
				consensus_ar[0][n+1]++;
				if ( letr <= 'Z' || (nuctype && letr=='n') )  
					consensus_ar[1][n+1] = letr;
				else								/* ELSE MAKE UPPERCASE */
					consensus_ar[1][n+1] = letr-32;
				break; /* OUT OF FOR m LOOP */
			}
		}
	}
	consensus_ar[0][0] = 1;			/* COUNTER ROW ON */
	consensus_ar[0][n] = '\0';		/*  & TERMINATED. THOUGH ARRAY INITIALIZED ABOVE, COMPILER MAKES FASTER CODE WITH THIS! */
	consensus_ar[1][0] = 1;			/* CONSENSUS ROW ON */
	consensus_ar[1][n] = '\0';		/*  & TERMINATED. NOT REDUNDANT TO COMPILER VOODOO?  */

	for (n = n_start; n <= n_end; n++) {
		if (nuctransit) {
			plustransit = 0;		/* RE-INITIALIZE AT EACH COLUMN */
			if ((conletr=con_align2D[MAXROW][n])=='R' || conletr=='Y') {
				consensus_ar[1][n+1] = conletr;
				plustransit = 1;
			}
			else conletr = blnk;
		}

		for (m = 1; con_align2D[m][0] != '\0'; m++) {
			if (isalpha(letr=con_align2D[m][n]) && isupper(letr)) {
				if (nuctransit && ((conletr=='R' && (letr=='G' || letr=='A')) || (conletr=='Y' && (letr=='C' || letr=='T'))) ) {
				}
				else if (letr != consensus_ar[1][n+1]) {
					if (consensus_ar[0][n+1] == 1) {	/* IF THIS IS THE FIRST CONFLICT NOTED AT THIS POSITION */
						if (badsites==0) {				/* IF FIRST BADSITE COLUMN THEN NOTE COORDINATES */	
							con_options[1][48] = m;
							con_options[1][49] = n;
						}

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
				while (col_isclear(con_align2D,n,m,-1)>-1)
					n++;
				consensus_ar[1][n+1] = letr;
				con_options[1][17] = m+1;		/* STORE HEIGHT IN HEIGHT SLOT */
			}
		} /* END OF FOR m=1... */
		if (plustransit)
			++consensus_ar[0][n+1];
	}

	if (con_options[0][20]) {							/* opt_K SHOW CONSENSUS ROW */	
		/* PRINT CONSENSUS ROWS */
		for (m = 1; consensus_ar[m][0] != '\0' && m < con_maxrows; m++) {
			line_end(BLOCKHEAD, 9, con_options, 9);
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
			line_end(BLOCKHEAD, 9, con_options, 9);	
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
		con_align2D[MAXROW][n] = consensus_ar[1][n+1];
	}

	magiclen = 1000;		/* MAGIC LENGTH NUMBER: 100 bp per bad site; SEARCH MAGIC TO FIND OTHER STOPGAPS */
	if (con_options[1][18] == 7 && badsites>0 && con_width>magiclen && con_width/badsites < magiclen) {
		warnhead('C');
		printf("Too many non-consensus sites (%d) for cinch-d to continue?\n\n", badsites);
		return(badsites);	
	}
	else if (badsites > 0) {
		return(badsites);				/* BAD CONSENSUS, REPORT IT 	*/
	}
	else								/* GOOD CONSENSUS, EARLY PASSES */
		return(0);
}

/*** FUNCTION 08b ******consensus_2D NUDGE**********************************/
unsigned int connudge(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width)
{
int badsites=0, m=0, n=0, n_end, x=1, nudge_row=0, nudge_col=0, nudge_span=0, frstletr;
int    height = con_options[1][17];
int con_width = con_options[1][32];
short unsigned int nuctype = con_options[1][13];	/* FOR SEQ TYPE, DNA=1, RNA=2, OTHER (NON-NA)=0 */
short unsigned int nuctransit = 0;					/* BIT FLAG FOR HANDLING NUCLEAR TRANSITIONS */ 
short unsigned int nudge_flag = 0;					/* BIT FLAG FOR PROCEEDING TO NUDGE */
short unsigned int plustransit=0;					/* BIT FLAG ADDENDUM FOR COUNTING BADSITES AT COL */
short unsigned int checktransit=0, badcol=0;
char blnk = con_options[1][11];
char letr=blnk, conletr=blnk, chkletr=blnk, badletr=blnk;
int con_maxrows=26;
int consensus_ar[26][MAXROW] = {{0}};	 	/* COL n=0 FOR BIT FLAG */
                                       		/* ROW m=0 FOR COUNTER */
											/* ROW m=1 FOR CONSENSUS */
											/* ROWS m>1 FOR VARIANTS STORAGE */
char mha_base62(int num);
void warnhead(char l);
int col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown); 

	n_end = n_start + n_width;

	if (nuctype == 1)	/* IF DNA */
		nuctransit = 1;

	/* FILL CONSENSUS ARRAY WITH FIRST LETTER IN EACH COL */
	for (n = 0; n <= con_width; n++) {
		for (m = 0; m < MAXROW; m++) {
			if ( (isalpha(letr=con_align2D[m][n]))) {
				consensus_ar[0][n+1]++;
				consensus_ar[1][n+1] = letr;
				break; /* OUT OF FOR m LOOP */
			}
		}
	}
	consensus_ar[0][0] = 1;			/* COUNTER ROW ON */
 	consensus_ar[0][n] = '\0';		/*  & TERMINATED. THOUGH ARRAY INITIALIZED ABOVE, COMPILER MAKES FASTER CODE WITH THIS! */
	consensus_ar[1][0] = 1;			/* CONSENSUS ROW ON */
 	consensus_ar[1][n] = '\0';		/*  & TERMINATED. NOT REDUNDANT TO COMPILER VOODOO?  */

	for (n = n_start; n <= n_end; n++) {
		badcol = 0;
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
						badcol = 1;
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
				nudge_flag = 1;
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
		if (con_options[1][57]>1) 
			printf("\n DEV: connudge top is clear at nudge_col=%d from nudge_row=%d", nudge_col+1, nudge_row+1);
		con_options[1][39]=4;
		return(0);
	}
	/* ****************************************************************** */

	if (con_options[1][57]>1) {	
		printf("\n DEV connudge: nudge_row=%d, nudge_col=%d, nudge_span=%d", nudge_row+1, nudge_col+1, nudge_span);
	}

	/* NUDGE IT! */
	for (m=nudge_row; con_align2D[m][0] != '\0'; m++) {
		for (n = n_end+2; n > 0; n--) {
			con_align2D[m][n] = con_align2D[m][n-1];
		}
	}

	/* ASSIGN CONSENSUS ROW LETTERS BUT VERIFY POST-NUDGE */
	for (n = n_end+1; n>=nudge_col+nudge_span; n--) {
		con_align2D[MAXROW][n] = consensus_ar[1][n];
	}
	con_align2D[MAXROW][n] = blnk;

	con_options[1][32]++;
	if (checktransit) {
		return(0);
	}
	else
		return(1);				/* RETURN SUCCESS */
}


/*** FUNCTION 09 *** CURRENTLY ONLY CYCLELIZES FOR k=2, AND FUDGES IT FOR k > 2 (PUSHES RIGHT) ******************/
unsigned int cyclelize(char cyc_align2D[][MAXROW], char *cyc_Seq_name, long int cyc_options[][62])
{
int cyc_col=0, cyc_row=0, a, b, i, j, k=2, kmer=0, m=0, n=0, threshold, r, s, x=0, y=0;	/* conflict_position */
int cyc_len   = cyc_options[1][1];
int cyc_width = cyc_options[1][32];					/* THIS IS opt_W SLOT TO STORE CURRENT 2-D WIDTH */
short unsigned int edge0=0, mono=0, mslip=0;
unsigned short int nuctype = cyc_options[1][13];	/* EQUALS ONE IF DNA STRING, TWO IF RNA, THREE IF PROTEIN */
unsigned short int nuctransit=0, dud_nudge=0;		/* BIT FLAG FOR HANDLING NUCLEOTIDE TRANSITIONS SILENTLY (IGNORING) */
unsigned short int tipcyc_flag=0;					/* BIT FLAG FOR TIP CYCLING OPPORTUNITY */
char blnk = cyc_options[1][11];
char letr, conletr, topletr;
char cyc_ar[MAXROW+1][MAXROW] = {{0}};
void clear_right(char swipe_align2D[][MAXROW], long int croptions[][62]);
unsigned int consensus_2D(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width);
unsigned int connudge(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width);
void mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
void mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);
short unsigned int print_2Dseq(char align2D_print[][MAXROW], int print_lenseq2D, char *printptr_Seq_name, long int poptions[][62]);
int col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown); 

	if (nuctype == 1)	/* IF DNA */
		nuctransit = 1;

	mha_writeback(cyc_align2D, cyc_ar, cyc_options);
	mha_writecons(cyc_align2D, cyc_ar, cyc_options);
 
	/* FLAG SPECIAL CASE OF CYCLING NEED AT n=0 COLUMN */
	if (cyc_options[0][5] == 0 && cyc_ar[0][0] == blnk)
		edge0 = 1;

	for (n = 0; n <= cyc_width; n++) {

		conletr = cyc_align2D[MAXROW][n];
		for (m = 1; cyc_align2D[m][0] != '\0' && m <= cyc_len; m++) {
			if (isalpha(letr=cyc_align2D[m][n])) {
				if (nuctransit && (col_isclear(cyc_align2D,n,m,-1) == -1)) {
					topletr = cyc_align2D[m][n];
				}
				else if (nuctransit && ( (conletr=='R' && (letr=='A'||letr=='G') && (topletr=='A'||topletr=='G')) ||
								         (conletr=='Y' && (letr=='C'||letr=='T') && (topletr=='C'||topletr=='T'))  )) {
					;	/* NOTHING: GO TO NEXT m */
				}
				else if (letr != cyc_ar[MAXROW][n] || edge0) {
					if (edge0) {
						while (isalpha(cyc_align2D[m][0]) == 0 && m <= cyc_len) {
							m++;
						}
						cyc_row = m;	/* THIS IS ROW COORDINATE OF NON-CONSENSUS */
						cyc_col = 0;	/* THIS IS COLUMN COORDINATE OF NON-CONSENSUS */
						while (cyc_align2D[0][n] == blnk) {	/* SCOOCH RIGHT */
								n++;
						}
						m = 0;			/* NEED TO RESET THIS TO FIRST ROW */
					}
					else {
						cyc_row = m;	/* THIS IS ROW COORDINATE OF NON-CONSENSUS */
						cyc_col = n;	/* THIS IS COLUMN COORDINATE OF NON-CONSENSUS */

						m=0;			/* RE-FIND FIRST LETTER OF COLUMN */
						while (isalpha(cyc_align2D[m][cyc_col]) == 0 && m < MAXROW) {
							m++;
						}
	
						/* SOMETIMES NEED TO SKIP ROWS AT BLEEDING EDGE OF SLIPS */	
						if (cyc_align2D[m][cyc_col+1] == '/') {
							while (cyc_align2D[m][cyc_col+1] == '/') {
								m++;
							}
							mono = 1;
						}
						else {
							while (cyc_align2D[m+1][n] != blnk) {	/* MOVE DOWN ****/
								m++;
							}
	
							while (letr != '0') {
								letr = cyc_align2D[m][n];			/* STORE THIS LETTER TO CHECK FOR INTERVENING REPEATS */
								for (i = m+1; i < cyc_row; i++) {
									if (cyc_align2D[i][n] == letr) {
										m = i;	/* SKIP TO THIS ROW */
										letr = cyc_align2D[m][n];	 
										break;	/* BREAK OUT OF FOR i LOOP & RECHECK W/IN WHILE LOOP */
									}
								}
								letr = '0';		/* JUST USING AS FLAG FOR WHILE LOOP */
							}

							while (cyc_align2D[m+1][n] == blnk) {	/* SCOOCH RIGHT */
								n++;
							}
						}
					} /* END OF ELSE (IF NOT edge0) */

					/* TIP-CYCLELIZE AS SOON AS DETECTED */
					tipcyc_flag = 0;
					for (i = m; i < cyc_row; i++) {
						if (cyc_align2D[i][cyc_col] != '/')
							break;
					}
					if (i == cyc_row && cyc_align2D[m-1][cyc_col+1] == '/') {
						j = 0;
						while (cyc_align2D[m][j] == blnk) {
							j++;
						}

						if (cyc_align2D[m][j] == cyc_align2D[m-1][cyc_col]) {
							tipcyc_flag = 1;		
							kmer = 4;		/* NOT ACTUAL kmer, JUST USING VAR TO CODE CYC TYPE */
							for (a = 0; a < m-1; a++) {
								for (b=0; (letr=cyc_align2D[a][b]) != '\0'; b++) {
									cyc_ar[a][b] = letr;
								}
								cyc_ar[m-2][b] = '\0';
							}
							for (b=0; b < cyc_col; b++) {
								cyc_ar[m-1][b  ] = cyc_align2D[m-1][b];
							}
								cyc_ar[m-1][b  ] = '/';
								cyc_ar[m-1][b+1] = '\0';
							for (b = 0; b < j; b++) {
								cyc_ar[m  ][b  ] = blnk;
							}
								cyc_ar[m  ][j  ] = cyc_align2D[m-1][cyc_col];
								cyc_ar[m  ][j+1] = '/';
								cyc_ar[m  ][j+2] = '\0';
							for (a = m; cyc_align2D[a][0] != '\0' && a < MAXROW; a++) {
								for (b=0; (letr=cyc_align2D[a][b]) != '\0'; b++) {
									cyc_ar[a+1][b] = letr;
								}
							}
							cyc_ar[MAXROW][cyc_col] = cyc_align2D[cyc_row][cyc_col];
							if (cyc_options[0][57])
								printf("\n DEV: TIP CYCLING OPPORTUNITY FOR cyc_col = %d; j=%d.", cyc_col+1, j+1);
						}
					}

					if (tipcyc_flag == 0) {
						/* NEXT: FIND WHAT k-MER TYPE THIS IS (CODED ONLY FOR k=2, k>2 REQUIRES CONSIDERATION OF INTRA-TR SLIPS) */
						k = 2;
						threshold = 0;
						for (i = m; i < cyc_row; i++) {
							if (i == m) {
								s = r = 1;			/* LOOK TO NEXT ROW */
							}
							else if (i == m+1) {
								r = -1;				/* LOOK TO NEXT ROW */
								s =  0;
							}
							else {
								s = r = -1;			/* LOOK TO PREV ROW */
							}

							if (threshold < 3*k && cyc_align2D[i][n] == blnk) {
								break;
							}

							for (j = n; j < n+k; j++) {
								/* BREAK OUT OF FOR j LOOP IF ROW NOT A k=2 TR ROW */
								if (j==n && (cyc_align2D[ i ][ n ] == blnk ||
											 cyc_align2D[i+s][n-1] != blnk || 
											 cyc_align2D[i+r][n+2] != '/' ) ) {
									break;
								}

								if (isalpha(letr=cyc_align2D[ i ][ j ]) && 
										         cyc_align2D[i+r][ j ]  == letr ) {
									threshold++;
								}
								else
									threshold = 0;
							} /* END OF FOR j COL LOOP */
						} /* END OF FOR i ROW LOOP */

						if (threshold >= 3*k) {
							kmer = k;
							if (cyc_options[0][57])		/* opt_v VERBOSITY */
								printf("\nDetected cyclelizable k-mer=%d TR (threshold=%d), cyc_row=%d, cyc_col=%d.", kmer, threshold, cyc_row+1, cyc_col+1);
						}

                    	if (kmer != 2)      	/* HACK UNTIL cyclelize is generalized for arbitrary k */ 
							kmer = 3;       	/* IF k=9, THEN BELOW WILL FUDGE CYCLELIZE BY PUSHING RIGHT. HACK WORKS FOR ALL k */ 
					}

					cyc_options[0][5] = kmer;	/* USING THE 0 ROW ABOVE PASS WIDTH ROW TO STORE cyclelize kmer VAR. */

					/* CYCLELIZE k=2 */	
					if (kmer == 2) {
						if (nuctransit) {
							for (j=0; j < kmer; j++) {	
								if ((letr=cyc_align2D[MAXROW][n+j]) == 'R' || letr == 'Y') {
									m = cyc_len+1;		/* BREAKS OUT OF FOR m LOOP */
									n = n+j;			/* ADVANCE n */
									cyc_options[1][5] = cyc_options[1][32] = cyc_width;	/* ASSIGN [32] CURRENT WIDTH and PASS x WIDTH HISTORY */
									break;				/* BREAKS OUT OF FOR j LOOP */
								}
							}
						}
						if (m != cyc_len+1) {
							/* CYCLELIZE FOR k=2 */
							for (i = (threshold/2 - 1); i > 0; i--) {
								cyc_ar[m][n+2] = cyc_align2D[m+1][n];
												      cyc_ar[m+1][n] = blnk;	/* STRICT CONSERVATION OF LETRZ */
								cyc_ar[m][n+3] = '/';
								cyc_ar[m][n+4] = '\0';
								m++;
							}

							/* WRITE LAST LINE OF k=2 TR */
							x=1;
							while ( (letr=cyc_align2D[m][n+x]) != '/') {
								cyc_ar[m-1][n + 2 + x] = letr;
								x++;
							}
							cyc_ar[m-1][n+2+x] = '/';
		
							/* SOMETIMES LAST LINE IS A NEW k=2 TR => NEW SLIP */
							if (cyc_ar[m-1][n+1] == cyc_ar[m-1][n+3] &&
							    cyc_ar[m-1][n+2] == cyc_ar[m-1][n+4]   ) {
								/* WRITE BOTTOM ROW UNDERNEATH FIRST IN ORDER TO NOT LOSE SEQUENCE */
								for (x = 0; (letr=cyc_ar[m-1][n+3+x]) != '/'; x++) {
									cyc_ar[m][n+1+x] = letr;
								}
								cyc_ar[m-1][n+3  ] = '/';
								cyc_ar[m-1][n+4  ] = '\0';
								cyc_ar[m  ][n    ] = blnk;
								cyc_ar[m  ][n+1+x] = '/';
								cyc_ar[m  ][n+2+x] = '\0';
								mslip = 1;
								y = 2;
							}
							m++;
		
							while (m < cyc_row) {
								for (n = 0; (letr=cyc_ar[m][n]) != '\0'; n++) 
									cyc_ar[m-1+mslip][n+2-y] = letr;
								m++;
							}
	
							for (n = 0; cyc_align2D[cyc_row][n] != blnk; n++)
								cyc_ar[cyc_row][n] = blnk;
	
							while (cyc_align2D[m][0] != '\0') {
								for (n = 0; (letr=cyc_align2D[m][n]) != '\0'; n++) 
									cyc_ar[m-1+mslip][n+2-mono] = letr;
								m++;
							}
							if (nuctransit) {	/* IF DOING CONSENSUS MATCHING PUSH CONSENSUS ROW TO RIGHT */
								for (n = cyc_width+2-mono; n > cyc_col+2-mono; n--)
									cyc_ar[MAXROW][n] = cyc_ar[MAXROW][n-2+mono];
								cyc_ar[MAXROW][cyc_col] = blnk;
							}
	
							cyc_width = cyc_width+2-mono;		
							cyc_options[1][5] = cyc_options[1][32] = cyc_width;	/* ASSIGN [32] CURRENT WIDTH and PASS x WIDTH HISTORY */
		
							for (n=0; n < cyc_width; n++)
								cyc_ar[m-1+mslip][n] = '\0';
							printf("\nCyclelizing (%c%c) to (%c%c) to resolve conflict at m=%d, n=%d.", 
										cyc_align2D[cyc_row][cyc_col+1], cyc_align2D[cyc_row][cyc_col  ],
										cyc_align2D[cyc_row][cyc_col  ], cyc_align2D[cyc_row][cyc_col+1], cyc_row+1, cyc_col+1);
						}
						n = cyc_width+1;	/* BREAKS OUT OF n LOOP */
						break;				/* BREAKS OUT OF m LOOP */
					} /* END OF IF kmer == 2 */

					/* TIP-CYCLELIZE */
					else if (tipcyc_flag) {
						/* WILL TIP-CYCLELIZE ABOVE AT FLAG CALL */
						cyc_options[1][5] = cyc_options[1][32] = cyc_width;	/* ASSIGN [32] CURRENT WIDTH and PASS x WIDTH HISTORY */
						tipcyc_flag = kmer = 0;
						n = cyc_width+1;	/* BREAKS OUT OF n LOOP */
						break;				/* BREAKS OUT OF m LOOP */
					}

					/* NUDGE-CYCLELIZE: */
					else {	
						if (connudge(cyc_ar, cyc_options, 0, cyc_width) == 0) {
							if (cyc_options[0][57])
								printf("\n DEV: dud_nudge");
							dud_nudge = 1;
							i = cyc_options[1][18];
							cyc_options[1][i] = cyc_width = cyc_options[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS x WIDTH HISTORY */
							n = cyc_width+1; 		/* BREAK OUT OF FOR n LOOP AFTER BREAKING OUT OF FOR m LOOP */
							break; 					/* BREAK OUT OF FOR m LOOP */
						}

						clear_right(cyc_ar, cyc_options);

						i = cyc_options[1][18];
						cyc_options[1][i] = cyc_width = cyc_options[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS x WIDTH HISTORY */
						n = cyc_width+1; 		/* BREAK OUT OF FOR n LOOP AFTER BREAKING OUT OF FOR m LOOP */
						break; 					/* BREAK OUT OF FOR m LOOP */
					} /* END OF ELSE (IF kmer != 2) */

				} /* END OF IF LETTER != CONSENSUS */
			} /* END OF IF ISALPHA */
		} /* END OF FOR m LOOP */
	} /* END OF FOR n LOOP */

	/* PUSH LEFT IF EMPTY: SHOULD MOVE TO ITS OWN fUNCTION IF NEEDED ELSEWHERE */
	if (edge0 && !dud_nudge) {
		i = 0; /* COUNTER FOR AMOUNT TO PUSH LEFT */
		for (n=0; n < cyc_width; n++) {
			for (m=0; cyc_ar[m][n] != '\0' && m < MAXROW; m++) {
				if (cyc_ar[m][n] != blnk) {
					n = MAXROW; 	/* TO BREAK OUT OF FOR n LOOP */
					break;			/* TO BREAK OUT OF FOR m LOOP */
				}
			}
			if (n < MAXROW)
				i++;
		} 
		for (m = 0; cyc_ar[m][0] != '\0' && m < MAXROW; m++) {
			for (n = 0; n < cyc_width+i; n++) {
				cyc_ar[m][n] = cyc_ar[m][n+i];
			}
		}
	} /*************** END OF PUSH LEFT GIVEN edge0 **************************/

	if (!dud_nudge)	{
		mha_writeback(cyc_ar, cyc_align2D, cyc_options);
		mha_writecons(cyc_ar, cyc_align2D, cyc_options);
		return(print_2Dseq(cyc_align2D, cyc_width, cyc_Seq_name, cyc_options));
	}
	else {
		return (0);
	}
}

/*** FUNCTION 10 **********************************************************/
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

/*** FUNCTION 11 **********************************************************/
unsigned int foam_2D(char foam_align2D[][MAXROW], long int foam_options[][62], int n_start, int n_width)
{
int m=0, n=0, n_end; 
char letr;
int foam_ar[3][MAXROW] = {{0}};	 	/* ROW m=0 FOR COLUMN LETR COUNTER 	*/
										/* ROW m=1 FOR CONSENSUS 			*/
										/* ROW m=2 FIRST LETTER ROW			*/
void line_end(int type, int c, long int lend_options[][62], int lcl_width);
char mha_base62(int num);

	n_end = n_start + n_width;

	for (n = n_start; n <= n_end; n++) {
		for (m = 0; m < MAXROW; m++) {
			if (isalpha(letr=foam_align2D[m][n])) {
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
		for (m = foam_ar[2][n]+1; foam_align2D[m][0] != '\0'; m++) {
			if ( isalpha(foam_align2D[m][n]) ) {
				++foam_ar[0][n];		/* INCREMENT COUNTER */
				foam_ar[1][n] = '.';
			}
		}
	}

	if (foam_options[1][20] > 1) {	/* opt_K CONSENSUS ROW */	
		line_end(BLOCKHEAD, 9, foam_options, 9);
		printf(" ");
		for (n = n_start; n < n_end; n++) {
			printf("%c", foam_ar[1][n]); 
		}
		printf(" = foam-free\n");		
	} 

	return(0);
}


/*** FUNCTION 12 ***** RETURN -1 IF ROW DOES NOT HARBOR NEXT FOAM-FREE POSITION *********************/
int next_foamfree(char check_array[][MAXROW], int row, int at_n) 
{
int col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown);

	int n = at_n;

	while ( isalpha(check_array[row][n]) ) {
		if (col_isclear(check_array,n,row,-1) == -1 && col_isclear(check_array,n,row,1) == -1)
			return(n);
		else
			n++;
	}	

	return(-1);

}		


/*** FUNCTION 14 **** RETURN -1 IF ALL CLEAR ABOVE OR BELOW ROW, OTHERWISE RETURN CLOSEST ROW WITH ALPHA **********************/
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
    return(-1);		/* RETURN SIGN INDICATED COLUMN IS CLEAR EITHER BELOW OR ABOVE */
}

/*** FUNCTION 14 **** RETURN -1 IF ALL CLEAR ABOVE AND BELOW ROW  **********************/
int col_totsclear(char charray[][MAXROW], unsigned int at_n, int row) 
{
int col_isclear(char check_array[][MAXROW], unsigned int at_n, int row, short int updown);

	if (col_isclear(charray, at_n, row, -1) + col_isclear(charray, at_n, row, 1) == -2)
		return(1);
	else
		return(0);
}


/*** FUNCTION 15 **** FOR UNIFORM FORMATTING OF 2-D LINE STARTS AND ENDS ****/
void line_end(int type, int c, long int lend_options[][62], int lcl_width)
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
char zero_tick = (char) lend_options[1][35];	/* opt_Z, Zero tick mark, default = 32 = ' ' */
int lcl_opt_r  = lend_options[0][53];

if (lend_options[0][11]==2 && lend_options[1][18]>1)	/* opt_B LEVELS FOR BLANKNESS IN FILLER & opt_I PASS NUM */
	zero_tick = '|';				/* 124 = '|'							*/

	if (type == 1)						/* FORMAT FOR LINE END. c IS CHAR. NUMBER */
		printf("%5d\n", c);

	if (lcl_opt_r == 1) {			/* FORMAT WITH LINE NUMBERING. DEFAULT */
		if (type == 0) {			/* FORMAT FOR LINE BEGINNING. c is LINE NUMBER (m+1) */
			if (c == 1)
				printf(" %4d. >", c);
			else 
				printf(" %4d. %c", c, zero_tick);
		}
		else if (type == 2) {			/* FORMAT FOR RULER. */
			printf("       %c%.*s\n", zero_tick, lcl_width, ruler + c);		/* DOUBLE USE OF c AS OFF-SET FOR RULER */
		}
		else if (type == 3)			/* FORMAT FOR NUMBERED SLIP LOCATION LINES. */
			printf("        ");
		else if (type == 4)			/* FORMAT FOR PATHBOX HEADERS, #'d */
			printf("      _|");
		else if (type == 5)			/* FORMAT FOR CONSENSUS BLOCKHEAD HEADERS, #'d */
			printf("       ");
		else if (type == 6) {		/* FORMAT FOR SLIP RULER. */
			ruler = rule1;
			printf("        %.*s\n", lcl_width, ruler + c);		/* DOUBLE USE OF c AS OFF-SET FOR RULER */
		}
	}
	else if (lcl_opt_r == 0) {	  	/* FORMAT WITHOUT LINE NUMBERS. r = REMOVED */
		if (type == 0) {		   	/* FORMAT FOR LINE BEGINNING. c is LINE NUMBER (m+1) */
			if (c == 1)
				printf("   >");
			else           
				printf("   %c", zero_tick);
		}
		else if (type == 2) {			/* FORMAT FOR RULER OR SLIPRULER. */
			printf("   %c%.*s\n", zero_tick, lcl_width, ruler + c);		/* DOUBLE USE OF c AS OFF-SET FOR RULER */
		}
		else if (type == 3)			/* FORMAT FOR UNNUMBERED SLIP LOCATION LINES. */
			printf("    ");
		else if (type == 4)			/* FORMAT FOR PATHBOX HEADERS, STILL #'d */
			printf("      _|");
		else if (type == 5)			/* FORMAT FOR CONSENSUS BLOCKHEAD HEADERS, #'d */
			printf("   ");
		else if (type == 6) {		/* FORMAT FOR SLIPRULER. */
			ruler = rule2;
			printf("   %c%.*s\n", zero_tick, lcl_width, ruler + c);		/* DOUBLE USE OF c AS OFF-SET FOR RULER */
		}
	}
}


/*** FUNCTION 16 ******** TO REPRESENT NUMBERS WITH SINGLE CHAR. ******/
char mha_base62(int num)
{
char c;
	if (num < 10) 
		c = '0' + num; 
	else if (num >=10 && num < 36) 
		c = 'A'+ num-10;
	else if (num >=36 && num < 58) 
		c = 'A' + num-4;
	else 
		c = '!';
	return(c);
}


/*** FUNCTION 17 ***************************************************************************************/
void mha_head(char lcl_seq_name, int lcl_width, long int lcl_options[][62])
{
			 /*0123456789*/
char h1[]=	"\nMHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//_"
			  "MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//_"
			  "MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//_"
			  "MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//_";
char h2[]=	"\n\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			  "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			  "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			  "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_";
char *h_rule = h1;						/* DEFAULT BANNER STYLE */
int min_len = 80;						/* MINIMUM LENGTH */
int med_len = 12*(lcl_width/10);		/* MEDIUM LENGTH, SCALING */
int max_len = lcl_options[1][58]+8;		/* MAXIMUM LENGTH */
int hr_len = min_len;					/* DEFAULT LENGTH OF HEADER BANNER */
int lcl_pass = lcl_options[1][18];		/* opt_I VALUE COUNTER FOR NUM OF PASSES */

	if (lcl_width+8 > hr_len && lcl_width+8 < lcl_options[1][58])
		hr_len = med_len;
	else if (lcl_width+8 >= lcl_options[1][58])
		hr_len = max_len;

	if (lcl_options[0][33])	/* opt_X == 1 */
		h_rule = h2;

	if (lcl_pass == 8) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: relax-2D [RELAXES HOMOPOLYMER RUNS THAT DID NOT AID cinch-d] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 7) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: cinch-d [RESCUES TR's INTERRUPTED BY DE NOVO REPEATS OF REPEATS] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 6) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: cinch_s  [CINCHES NON-COMPLEX k-mer TR's IN A SINGLE LINE] (width = %d)\n\n", lcl_pass, lcl_width);
	}
	else if (lcl_pass == 5) {
		printf("%.*s\n", hr_len, h_rule);
		printf("2-D PASS #%d: cyclelize [RESCUES TR's OBSCURED BY INITIAL CYCLING FRAME] (width = %d)\n\n", lcl_pass, lcl_width);
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
	else if (lcl_pass == 0) {
		printf("\nOriginal string (sequence %c, length = %d):\n", lcl_seq_name, lcl_width);
	}
	else 
		printf("%.*s\n", max_len, h_rule);
}


/*** FUNCTION 18 **************************************/
void mha_randomize1(char input_seq[MAXROW]) 
{
    int i, rand_num;
	int rand_seqlen = strlen(input_seq);
	char random_seq[MAXROW];

    for (i = 0 ; i < rand_seqlen; i++)
    {
        rand_num = rand() % rand_seqlen;
		random_seq[i] = input_seq[rand_num];
    }

	for (i = 0; i < rand_seqlen; i++) 
		input_seq[i] = random_seq[i];
	input_seq[i] = '\0';
}


/*** FUNCTION 19 **************************************/
void mha_randomize2(char input_seq[MAXROW], int rsize) 
{
int i, j, tmp;
int rnumbers[MAXROW];
int seqlen = strlen(input_seq);
char random_seq[MAXROW];

	for (i = 0; i < rsize; i++)
		rnumbers[i]= i;
  
	for (i = rsize - 1; i > 0; i--) {
    	j = random_i(i + 1);
        tmp = rnumbers[j];
        rnumbers[j] = rnumbers[i];
        rnumbers[i] = tmp;
	}

	for (i = 0; i < rsize; i++) {
		random_seq[i] = input_seq[(rnumbers[i] % seqlen)];
	}		

	for (i = 0; i < rsize; i++) 
		input_seq[i] = random_seq[i];
	input_seq[i] = '\0';	/* B/C SOMETIMES input_seq > size */ 
}


/*** FUNCTION 20 ******************************************************************************/
void mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62])
{
char letr;
int wb_mwrap     = woptions[1][22];		/* opt_M long_homopolymer_run */
char wopt_Q_left = woptions[1][26];		/* LHS character delimiter for homopolymer Run */
char wopt_R_rght = woptions[1][27];		/* RHS character delimiter for homopolymer Run */
int m=0, n=0, widest_n=0;
void clear_2D_ar(char wipe_align2D[][MAXROW]);
void clear_right(char swipe_align2D[][MAXROW], long int croptions[][62]);

	/* WRITE BACK TO align2D_prev */
	clear_right(lcl_align2D, woptions); 
	clear_2D_ar(align2D_prev);
	for (m = 0; lcl_align2D[m][0] != '\0' && m < MAXROW; m++) {
		for (n = 0; (letr=lcl_align2D[m][n]) != '/' && letr != '>' && letr != wopt_R_rght; n++) {
			align2D_prev[m][n] = letr;
			if (letr == wopt_Q_left && lcl_align2D[m][n+wb_mwrap+1] == wopt_R_rght ) 
				align2D_prev[m][n+wb_mwrap+2] = '\0';
		}
		align2D_prev[m][n  ] = letr;	/* MHA-STANDARD TERMINATOR  */
		if (n > widest_n)
			widest_n = n;
		if (letr == '>') {				/* ASSIGN CINCH-WIDTH TO CURRENT [32]	*/
			woptions[1][32] = widest_n;
		}
	}
}

/*** FUNCTION 21 ******************************************************************************/
void mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62])
{
int i=0;

	/* COPY CONSENSUS ROW FROM ARRAY ONE TO ARRAY TWO */
	for (i = 0; i < wroptions[1][32]; i++) {
		align2D_two[MAXROW][i] = align2D_one[MAXROW][i];
	}
}


/*** FUNCTION 22 ******* LIKE mha_writeback EXCEPT CHECKS TO MAKE UPPERCASE *******************/
void mha_UPPERback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62])
{
char letr;
short unsigned 
     int nuctype = woptions[1][13]; 			/* EQUALS ONE IF DNA, TWO IF RNA */
int wb_mwrap     = woptions[1][22];				/* opt_M long_homopolymer_run */
char wopt_Q_left = woptions[1][26];				/* LHS character delimiter for homopolymer Run */
char wopt_R_rght = woptions[1][27];				/* RHS character delimiter for homopolymer Run */
int i=0, m=0, n=0, widest_n=0;
void clear_2D_ar(char wipe_align2D[][MAXROW]);
void clear_right(char swipe_align2D[][MAXROW], long int croptions[][62]);

	/* WRITE BACK TO align2D_prev */
	clear_right(lcl_align2D, woptions); 
	clear_2D_ar(align2D_prev);		/* DOES NOT CLEAR MAXROW ROW */
	for (m = 0; lcl_align2D[m][0] != '\0'; m++) {
		for (n = 0; (letr=lcl_align2D[m][n]) != '/' && letr != '>' && letr != wopt_R_rght; n++) {
			if (nuctype) { /* IF NON-ZERO LIKE DNA (1) OR RNA (2) */
				if (isalpha(letr) && letr >= 'a' && letr != 'n')
					align2D_prev[m][n] = letr - 32;		/* MAKE UPPERCASE IF LOWER & NOT 'n' */
				else									/* ELSE KEEP SAME CASE */
					align2D_prev[m][n] = letr;
			}
			else if (isalpha(letr) && letr >= 'a')
				align2D_prev[m][n] = letr - 32;		/* MAKE UPPERCASE IF LOWER & NOT DNA */
			else									/* ELSE KEEP SAME CASE */
				align2D_prev[m][n] = letr;

			if (letr == wopt_Q_left && lcl_align2D[m][n+wb_mwrap+1] == wopt_R_rght ) 
				align2D_prev[m][n+wb_mwrap+2] = '\0';
		}
		align2D_prev[m][n  ] = letr;	/* MHA-STANDARD TERMINATOR  */
		if (n > widest_n)
			widest_n = n;
		if (letr == '>') {				/* ASSIGN CINCH-WIDTH TO CURRENT [32]	*/
			woptions[1][32] = widest_n;
		}
	}

	/* FLUSH POST-TERMINATOR ENDS OF align2D WITH NULLS */
	for (m = 0; lcl_align2D[m][0] != '\0'; m++) {
		n = 0;
		while ( (letr=lcl_align2D[m][n]) != '/' && letr != '>' && letr != wopt_R_rght) 
			n++;	
		for (i = 1; n+i <= woptions[1][1]; i++)
			align2D_prev[m][n+i] = '\0';	
	}
}


/*** FUNCTION 23 **************************************************************/
/* WHEN CALLING THIS fUNCTION FOR THE FIRST TRANSIT USE last_transit = -1 */
int next_transit(char lcl_align2D[][MAXROW], int last_transit) 
{
int i;
char letr;
	
	for (i = last_transit + 1; (letr=lcl_align2D[MAXROW][i]) != '>' && letr != '\0'; i++) {
		if (letr == 'R' || letr == 'Y')
			return(i);		
	}
	return(last_transit);	/* WILL RETURN -1 IF CALLED FOR THE FIRST TIME WITH -1 AND NO TRANSITIONS FOUND */
							/* WILL RETURN last_transit IF NO OTHERS FOUND SINCE LAST TRANSTION */
}


/*** FUNCTION 24 ********************************/
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


/*** FUNCTION 25 ********************************************************************************************/
short unsigned int print_2Dseq(char align2D_print[][MAXROW], int print_lenseq2D, char *printptr_Seq_name, long int poptions[][62])
{
unsigned int consensus_2D(char con_align2D[][MAXROW], long int con_options[][62], int n_start, int n_width);
int count_wrap_blocks(int lcl_width, int lcl_opt_w);
unsigned int foam_2D(char foam_align2D[][MAXROW], long int foam_options[][62], int n_start, int n_width);
void line_end(int type, int c, long int lend_options[][62], int lcl_width);
void mha_head(char lcl_seq_name, int lcl_width, long int lcl_options[][62]);
char mha_base62(int num);
void print_blockhead(int a, int b);
void warnhead(char l);
int all_clear;		/* COUNTER VARIABLE USED FOR CHECKING NEED TO PRINT BOTTOM ROWS */ 
int blocks2D, b=0, c=0, carry_over=0, d=0, fudge=0, g, h, i, j, j_start, j_end, m, m_start=0, n, mmsites=0, max_n=0;
char letr = 'B', next = 'B';			/* Begin the Beguine */
char blnk  = poptions[1][11];			/* opt_B blank character */
int cinchwidth = poptions[1][32];		
int cip_linewidth = poptions[1][58];
char popt_Q_left = poptions[1][26]; 	/* LHS character delimiter for homopolymer Run */
char popt_R_rght = poptions[1][27];		/* RHS character delimiter for homopolymer Run */
int head_start;							/* USE TO PASS RULER OFF-SET TO line_end() */
int scrimmageline;						/* USE TO INCREMENT AND TEST IF FILLER IS NEEDED, CAN BE OPTION TO DO SO */
int blnk_lvl = poptions[0][11];			/* USE TO STORE DEGREE OF BLANKNESS */
int length = poptions[1][1];
char tick = ':'; 						/* OTHER POSSIBILITIES: |, ^ */
short unsigned int lcl_opt_F;

	if (poptions[0][11] > 1)
		tick = blnk;

	mha_head(*printptr_Seq_name, print_lenseq2D, poptions);

	blocks2D = count_wrap_blocks(print_lenseq2D, poptions[1][58]);

	for (j = 0; j < blocks2D; j++) {
		if (blocks2D != 1)
			print_blockhead(j+1, blocks2D);

		/* FAST FORWARD TO INFORMATIVE ROW */
		m_start = 0;
		while (align2D_print[m_start][(j * poptions[1][58])] == '\0') {
			m_start++;
		}

		scrimmageline = 1;
				max_n = 0;

		for (m = m_start; align2D_print[m][0] != '\0'; ) {

			b = d = 0;		/* VAR b WILL COUNT BLANKS, VAR d WILL COUNT ALPHA CHAR ANEW FOR EACH ROW */

			j_start =   j   * poptions[1][58];
			j_end   = (j+1) * poptions[1][58];

			for (n = j_start; n < j_end && (letr=align2D_print[m][n])!='/' && letr!='>' && letr!=popt_R_rght; n++) {
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
				for (g = m+1; align2D_print[g][0] != '\0' && g < MAXROW; g++) {
					for (h = 0; h < j_end; h++) {
						if (isalpha(align2D_print[g][h]))
							all_clear = 0;
					}
				}
				if (all_clear) 
					break;
			}

			line_end(START, m+1, poptions, 0);

			for (n = j_start; n < j_end  && (letr=align2D_print[m][n])!='>' && letr!=popt_R_rght && letr!='\0'; n++) {
				if (d+b == 0) {		/* HANDLES BLOCK CONTINUATION LINES THAT HAVE ZERO CHARACTERS B/C OF ITS TUCK & LENGTH */
					break;
				}
				else if ((letr=align2D_print[m][n]) == blnk && (n+1) % 10 == 0 && blnk_lvl < 2)		
					printf("%c", tick);		/* PRINT TICK MARKS AT 10 bp INTERVALS IF NOT BLANK SPACE */
				else if (poptions[1][11] == 32)		/* opt_B BLANK = SPACE ' ' */	
					printf("%c", letr);
				else					/* opt_B BLANK = FULLSTOP '.' */
					printf("%c", letr);

				/* TURN ON OPTION opt_F TO USE BLANK CHAR TO FILL < SCRIMMAGE IF EXCESSIVELY SHORT */
				if (poptions[0][15] != 1 && align2D_print[m][n] == '/' && scrimmageline-(b+d-2) > 10)
					lcl_opt_F = 1; 	/* TURN ON local opt_F IN THIS CASE */

				/* OPTION opt_F TO USE BLANK CHAR TO FILL < SCRIMMAGE */
				if ((poptions[0][15] || lcl_opt_F) && align2D_print[m][n] == '/' && b+d-2 < scrimmageline) {
					for (i = 0; b+d+i < scrimmageline; i++) {
						if ( (b+d+2+i)%10 == 0)		/* WHY 2? +1 FOR STARTING AT 0, +1 FOR '/' CHAR */
							printf("%c", tick);		/* PRINT TICK MARKS AT 10 bp INTERVALS IF NOT BLANK SPACE */
						else
							printf("%c", blnk);		/* PRINT LINE-END FILLER CHARACTER */
					}
					lcl_opt_F = 0;
				}
			} /* END OF n PRINTING LOOPS */

			next = align2D_print[m][n];	/* TERMINAL CHARACTERS WILL ALSO PRINT AT END OF BLOCK IF...	*/ 
										/* ...THEY ARE FIRST CHARACTER OF NEXT BLOCK					*/

			if (next == '/' || next == popt_R_rght)		/* TO INDICATE ADJACENCY TO SLIP, WHICH ALSO WILL SHOW IN NEXT BLOCK */
				printf("%c", next);
			else if (next==blnk || (isalpha(next) && d+b != 0))	/* TRUE IF AT EDGE OF opt_w WINDOW */
				printf("=>");			/* TO INDICATE CONTINUATION TO NEXT BLOCK */
			else if (next == '\0' && n == j_start)		
				printf("%c", blnk);
			else if (next == '>') {		/* LAST CHARACTER */
				printf("%c", next);
				if (poptions[0][21])	/* opt_L = LINE END NUMBERING (ON) */
					line_end(END, c, poptions, 0);
				else
					printf("\n");
				break;	/* BREAK OUT OF m=m_start FOR LOOP */
			}

			if (poptions[0][21])	/* opt_L = LINE END NUMBERING ON */
				line_end(END, c, poptions, 0);
			else
				printf("\n");
			m++;

		} /* END OF m=m_start FOR LOOP */

		poptions[1][17] = m+1;	/* ASSIGN COUNTED HEIGHT (# OF ROWS) TO HEIGHT SLOT */

		/* PRINT RULER */
		if (poptions[0][11] < 5) {
			head_start = (j * poptions[1][58]) % 10;
			if (j+1 < blocks2D) {
				line_end(RULER, head_start, poptions, cip_linewidth);
			}
			else {
				line_end(RULER, head_start, poptions, (cip_linewidth=max_n));
			}
		}
		/* *********** */

		/* PRINT NUMBERS FOR CONSENSUS RULER */
		if (poptions[0][11] < 4) {
			line_end(SLIPS, head_start, poptions, cip_linewidth);
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
		mmsites = mmsites + consensus_2D(align2D_print, poptions, j_start, cip_linewidth);

		if (poptions[1][20] > 1 && poptions[1][18] > 6 && mmsites == 0) {
			foam_2D(align2D_print, poptions, j_start, cip_linewidth);
		}
	} /* END OF FOR j PRINTING LOOP */

	if (c == poptions[1][1] && mmsites == 0) {
		poptions[0][10] = 1000;	
		printf("\n Successfully 2-D self-aligned all %d letters from formatted sequence. \n", c);
		return(0);
	}
	else if (c == length) {
		i=poptions[0][10] = round((1000*(cinchwidth-mmsites))/cinchwidth);	
		printf("\n 2-D self-aligned all %d letters from formatted sequence but consensus indicates more cinching required. \n", c);
		return(i);
	}
	else if (c < length) {
		i=poptions[0][10] = round((1000*(cinchwidth-mmsites-(length-c)))/cinchwidth);
		warnhead('-');
		printf("2-D printing is missing %d letter(s)!\n", length-c); 
		return(0);
	}
	else {
		i=poptions[0][10] = round((1000*(cinchwidth-mmsites-(c-length)))/cinchwidth);	
		warnhead('+');
		printf("2-D printing left %d extra letter(s)!\n", c-length);
		return(0);
	}
}


/*** FUNCTION 26 *************************************/
void print_base62_table(long int boptions[][62])
{
unsigned int i=0;
void mha_head(char lcl_seq_name, int lcl_width, long int lcl_options[][62]);

	mha_head(0, 80, boptions);
	printf("\nMHA base 62 single-letter representation:\n\n Base 10: ");
	for (i = 0; i < 31; i++)
		printf("%2d ", i);
	printf("\n Base 62: ");
	for (i = 0; i < 31; i++)
		printf("%2c ", (char) boptions[3][i]);
	printf("\n\n Base 10: ");
	for (i = 31; i < 62; i++)
		printf("%2d ", i);
	printf(">61\n Base 62: ");
	for (i = 31; i < 62; i++)
		printf("%2c ", (char) boptions[3][i]);
	printf("  !\n\n");
}


/*** FUNCTION 27 **********************************/
void print_blockhead(int a, int b)              /**/
{                                               /**/
	if (a == 1)                                 /**/
	    printf("   Block %d of %d:\n", a, b);   /**/
	else                                        /**/
	    printf("\n   Block %d of %d:\n", a, b); /**/
}                                               /**/
/**************************************************/


/*** FUNCTION 28 *******************************************************/
void recover_1D(char recovered_1D[MAXROW], char rec_align2D[][MAXROW], long int rec_options[][62]) 
{
int m=0, n=0, x=0;
int rec_lenseq = rec_options[1][32];	/* opt_W CURRENT 2-D WIDTH */
char letr;
char rec_R_rght = rec_options[1][27];

	for (m=0; (rec_align2D[m][0] != '\0') || (m < MAXROW); ) {
		for (n=0; (letr=rec_align2D[m][n]) != '/' && letr!=rec_R_rght && letr!='>' && n < rec_lenseq; n++) {
			if (isalpha(letr)) {
				recovered_1D[x] = letr;
				x++;
			}
		} /* END OF FOR n LOOP */
		if (letr=='/' || letr==rec_R_rght)
			m++;
		else if (letr=='>') {
			recovered_1D[ x ] = letr;
			break; /* BREAK OUT OF FOR m LOOP */
		}
	} /* END OF FOR m LOOP */

}


/*** FUNCTION 29 *************************************/
void relax_2D(char align2D_pass8[][MAXROW], long int roptions[0][62])
{
int height=0, i, j, m, n, rlx_col=0, v=0, w=0, z=0;
int width = roptions[1][32]; /* CURRENT WIDTH */
char blnk = roptions[1][11];
unsigned short int nuctype = roptions[1][13], nuctransit=0;
char rlx_opt_R_rght = roptions[1][27];
char letr = 'B'; 	/* Begin the Beguine */
char rlx_align2D[MAXROW+1][MAXROW];
void mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
void mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);

	if (nuctype == 1)		/* IF DNA */
		nuctransit = 1;

	mha_writeback(align2D_pass8, rlx_align2D, roptions); 
	mha_writecons(align2D_pass8, rlx_align2D, roptions); 
 
	while (align2D_pass8[height][0] != '\0') {
		height++;
	}
	roptions[1][17] = height;	/* 9 PLUS 8 (H) EQUALS 17 (base 62 H) */

	for (n=0; n < width; n++) {
		v = w = z = m = 0;

		while (align2D_pass8[m][n] == '\0') {
			m++;
		}

		while (isalpha(letr=align2D_pass8[m][n]) && isalpha(align2D_pass8[m][n+1])) {
			n++;
		}

		if (align2D_pass8[m][n+1] == '/') {		/* EDGE DETECTED */
			while (isalpha(align2D_pass8[m+v][n-1]) && 
				   isalpha(align2D_pass8[m+v][ n ]) &&
						   align2D_pass8[m+v][n+1] == '/'   ) {
				v++;
			}
			while (align2D_pass8[m+v+w][n-1] == blnk && 
				   align2D_pass8[m+v+w][ n ] == letr && 
				   align2D_pass8[m+v+w][n+1] == '/'     ) {
				w++;
			}
			if (align2D_pass8[m+v+w][n-1] == blnk && align2D_pass8[m+v+w][n] == letr) 
				w++;		/* TOTAL w EQUALS LENGTH OF HOMOPOLYMER TO RELAX */

			for (i = m+v+w; i < height; i++) {
				for (j = 0; j <= n; j++) {
					if (isalpha(align2D_pass8[i][j]))
						w = 0;
				}
			}

			if (w > 0) {
				while (align2D_pass8[m+v+w+z][n] == blnk) {	
					z++;
				}
				while (m+v+w+z < height) {
					if (isalpha(align2D_pass8[m+w+z][n]))
						break; 	/* BREAK OUT OF WHILE m+z < height LOOP */
					else
						z++;	
				}
				if (m+v+w+z == height) {		/* WRITE TO LOCAL 2D ARRAY */
					for (j = 0; j < height; j++) 
						rlx_align2D[m+v - rlx_col + j][n + rlx_col] = blnk;	
					for (i = 0; i < w; i++) {
						rlx_align2D[m+v-1 - rlx_col][n + i + rlx_col + 1] = letr;
						for (j = 0; j < height; j++) 
							rlx_align2D[m+v - rlx_col + j][n + i + rlx_col + 1] = blnk;
					}
					rlx_col = rlx_col + w;

					/* REWRITE REST OF ARRAY TO NEW COORDINATES */
					for (i = m+v+w-1; i < height; i++) {
						for (j = n+1; (letr=align2D_pass8[i][j]) != '\0'; j++) 
							rlx_align2D[i-rlx_col][j+rlx_col] = letr;
					}

					if (nuctransit) {
						j = roptions[1][32] + rlx_col + w;
						rlx_align2D[MAXROW][j] = '\0';
						while (j > n+rlx_col-w) {
							rlx_align2D[MAXROW][j] = rlx_align2D[MAXROW][j-w];
							j--;
						}
						for (j = n+rlx_col-w+1; j <= n+rlx_col; j++) 
							rlx_align2D[MAXROW][j] = '\0';
					}

				}
			} /* END OF IF w > 1 */
		}
		else if (align2D_pass8[m][n+1] == '>' || align2D_pass8[m][n+1] == rlx_opt_R_rght) {
			rlx_align2D[m - rlx_col][n + rlx_col] = letr;
			rlx_align2D[m - rlx_col][n + rlx_col + 1] = align2D_pass8[m][n+1];
			rlx_align2D[m - rlx_col][n + rlx_col + 2] = '\0';

			if (align2D_pass8[m][n+1] == '>') {
				roptions[1][32] = n + rlx_col + 1;
			    align2D_pass8[m][n+2] = '\0';
			}

		}

	} /* END OF FOR n LOOP */

	mha_writeback(rlx_align2D, align2D_pass8, roptions);
	mha_writecons(rlx_align2D, align2D_pass8, roptions); 

	i = roptions[1][18];
	roptions[1][i] = roptions[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS [9] WIDTH HISTORY */
}


/*** FUNCTION 30 **************************************/
void usage(char usage_version[], unsigned int FY_size)
{
	printf("This is maximal version %s, a program for maximal homology alignment (MHA) and surfing the foamy waves!\n", usage_version);
	printf("Usage: ./maximal -[OPTIONS] sequence.txt (OPTIONAL FASTA HEADER AND NON-ALPHA SEQUENCE CHARACTERS ARE IGNORED.)\n"
							"\t\t -(a-h)   USE ONE OF 8 BUILT-IN EXAMPLE SEQUENCES INSTEAD OF SEQUENCE FILE.\n"
							"\t\t -i       USE INPUT SEQUENCE FROM FILE, NAMED AS COMMAND LINE ARGUMENT. '-i' IS OPTIONAL.\n"
							"\t\t -k       SHOW k-MER COUNT.\n"
							"\t\t -l       SHOW CINCH-T SLIP LOCATIONS, NUMBER (BASE 62), AND K-MER SIZES (BASE 62).\n"
							"\t\t -n       DO NOT DO RELAX-2D PASS.\n"
							"\t\t -o       SHOW ORIGINAL INPUT STRING.\n"
							"\t\t -p       SHOW PARAMETERS.\n"
							"\t\t -r       SHOW ROW NUMBER.\n"
							"\t\t -s       SILENCE WRITING TO NORMAL OUTPUT FILE (Surf_wavereport.mha).\n"
							"\t\t -u(u...) DO NOT WRAP OUTPUT (WRAP IN ONE BLOCK); EACH '-u' WRAPS OUTPUT INTO ANOTHER BLOCK.\n"
							"\t\t -v(v)    VERBOSE MODE: \"maximal -FLIRlopr\" + VERBOSITY. TOGGLE OFF I & R WITH -vI & -vR.\n"
							"\t\t -x       EXTRA SQUEEZE: REDUCE CINCH-T THRESHOLDS FOR TRANSITION MATCHING BY ONE.\n"
							"\t\t -z       ENSURE MISMATCH SCORE IS ZERO; ALTERS PATHBOX DISPLAY.\n"
							"\t\t -B       USE EMPTY SPACE FOR BLANK CHARACTER INSTEAD OF PERIOD.\n"
							"\t\t -B(BBBB) DO NOT SHOW TICK MARKS, ZERO TICKLINE, RULER #'s, & RULER, RESPECTIVELY.\n"
							"\t\t -C       SHOW BASE 62 SINGLE-LETTER CODE (USED FOR REPEAT #'s).\n"
							"\t\t -D       SHOW DIAGNONAL THRESHOLD VALUES FOR TRANSITIONS HANDLING.\n"
							"\t\t -F       DO NOT FILL SHORT LINE ENDS TO SCRIMMAGE LINE.\n"
							"\t\t -H       SHOW HELP, INCLUDING BUILT-IN EXAMPLE SEQUENCES.\n"
							"\t\t -I       SHOW INITIAL PASSES.\n"
							"\t\t -K       SHOW CONSENSUS ROWS, -KK FOR FOAM-FREE SEGMENTS.\n"
							"\t\t -L       SHOW POSITION AT END OF EACH ROW.\n"
							"\t\t -M       DOUBLE THE DEFAULT LONG HOMOMONOMER TRACT WRAP LENGTH.\n"
							"\t\t -O       SAVE (APPEND) FORMATTED 2-D ALIGNMENT TO CONSENSUS FILE (Surf_barrels.log).\n"
							"\t\t -OO      SAVE (APPEND) RAW 2-D ALIGNMENT TO SPECIAL 2-D INPUT FILE (TUBES.mha).\n"
							"\t\t -P       SHOW PATHBOX.\n"
							"\t\t -R       RECOVER AND CHECK 1-D SEQUENCE FROM 2-D SELF-MHA.\n"
							"\t\t -X       RUN ON SCRAMBLED SEQUENCE OF SAME LENGTH.\n"
							"\t\t -XX      RUN ON FISHER-YATES SHUFFLED SEQUENCE OF LENGTH %d.\n", FY_size);
	printf("Example usage: ./maximal -KnO sequence-file.txt\n\n");
}


/*** FUNCTION 31 *** RETURNS 1 FOR EARLY EXIT OPTIONED BY USER REQUEST ***********/
short unsigned int user_query(unsigned int pass_num)
{
char key_input[100];
char letr;
int i=0;

	printf("\n MAXIMAL: May I continue to the next pass (PASS #%d)?\tYOU: ", pass_num+1);
	scanf("%s", key_input);

	while (key_input[i] == ' ') {
		i++;
	}
	if((letr=key_input[i]) == 'N' || letr == 'n' || letr == 'Q' || letr == 'q' || letr == 'I' || letr == 'i' || letr == 'F' || letr == 'f') {
		printf(" MAXIMAL: Okay. I'll finish up nicely right away.\n");
		return(0);
	}
	else if (letr=='Y'||letr=='y'||
			 letr=='S'||letr=='s'||
			 letr=='C'||letr=='C'||
			 letr=='X'||letr=='x'||
			 letr=='K'||letr=='k'||
			 letr=='H'||letr=='h'||
			 letr=='O'||letr=='o'||
			 letr=='D'||letr=='d'||
			 letr=='A'||letr=='a') {
		/* YES, SURE, CONTINUE, SHI DE, XING, KE YI, HAO DE, MEI CUO, HAI, OUI, DA, ORALE, OK, K, SI, ALLRIGHT */ 
		printf(" MAXIMAL: Great! I will continue on to the next passes.\n");
		return(1);
	}
	else {
		printf(" MAXIMAL: I do not understand you, so I will continue to fulfill my purpose in life.\n");
		return(1);
	}
}


/*** FUNCTION 32 *************************************/
void warnhead(char l)
{ 
/*	printf("%c", 7); */	/* BELL CHARACTER */
	printf("\n * Notice (%c): ", l);
} 

/*** FUNCTION 33 *************************************/
short int pushdown(char pusharray[][MAXROW], int push_m, int push_n, long int push_options[0][62])
{
int i=0, j=0;
int h = push_options[1][17];		/* height slot */
char letr;
char blank = push_options[1][11];

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

/**************************************************************************************/
int get2Dtucknum(char arrayA[][MAXROW], char arrayB[][MAXROW], long int options[][62]) 
{
int i=0, j=0, height=0;
int width = options[1][32];
char letr, opt_R_rght = options[1][27];
int bottom[MAXROW] = {0};
int    top[MAXROW] = {0};
char mha_base62(int num);	

	for (i=0; arrayA[i][0] != '\0'; i++) {
		;
	}
	height = i;

	/**** FILL BOTTOM BORDER ARRAY; BOTTOM BORDER IS FOR ARRAY ON TOP ****/
	j = width;
	while (j>0) {
		if (i>0 && (isalpha(arrayA[i-1][j]) || isalpha(arrayA[i-1][j+1]))) {
			bottom[j-1] = i;
			j--;
		}
		else {
			i--;
		}
	}

	/**** FILL TOP BORDER ARRAY; TOP BORDER IS FOR ARRAY ON BOTTOM  ****/
	top[0] = 0;
	for (j=0; j<width; j++) {
		for (i=0; i<MAXROW; i++) {
			if (isalpha(letr=arrayB[i][j])) {
				top[j] = i;
				break;
			}
            else if (letr == '/' || letr == opt_R_rght) {
				top[  j] = i;
				top[++j] = i;
				break;
			}
			else if ((letr=arrayB[i+1][j-1]) == '/' || letr == opt_R_rght) {
				top[j] = i+1;
				break;
			}
			else if (arrayB[i][j]=='>') {
				i = MAXROW;
				break;
			}
		}
	}

	short unsigned int devReport=0;	/* FLIP TO REPORT */
	/**************************************************/
	if (devReport) {	/* CODE DEVELOPMENT REPORTING */
		printf("\n");	
		for (j=0; j < width; j++)
			printf("%c", mha_base62(bottom[j]));	
			printf(" < bottom\n");
		for (j=0; j < width; j++)
			printf("%c", mha_base62(top[j]));	
			printf(" < top\n");
	}
	/******************************************/

	bottom[0] = top[0] + height - bottom[0] - 1;
	if (devReport)
		printf("%c", mha_base62(bottom[0]));
	for (j=1; j < width; j++) {
		if ((i=top[j] + height - bottom[j] - 1) < bottom[j-1]) {
			bottom[j] = i;
		}
		else
			bottom[j] = bottom[j-1];
		if (devReport)
			printf("%c", mha_base62(bottom[j]));
	}
	if (devReport)
		printf(" < height - bottom + top - 1\n");
	i=bottom[j-1];
	return(i);
}


/*** FUNCTION 34 **********************************************************/ 
short int tucksense(char tuckarray[][MAXROW], long int tuck_options[0][62])
{
int i=0, j=0, j_stop=0;
char letr;
int bad_m = tuck_options[1][48];	/* THESE ARE FILLED BY consensus FUNCTION */
int bad_n = tuck_options[1][49];	/* THESE ARE FILLED BY consensus FUNCTION */
char badletr = tuckarray[bad_m][bad_n];
char blank = tuck_options[1][11];
short int tuckcall=0;
char lcl_align2D[MAXROW+1][MAXROW] = {{0}};
unsigned short int nuctype = tuck_options[1][13], nuctransit=0;
short int pushdown(char pusharray[][MAXROW], int push_m, int push_n, long int push_options[0][62]);
short unsigned int print_2Dseq(char align2D_print[][MAXROW], int print_lenseq2D, char *printptr_Seq_name, long int poptions[][62]);
void mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW], long int wroptions[][62]);
void mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW], long int woptions[][62]);

	if (nuctype == 1)		/* IF DNA */
		nuctransit = 1;

	if (tuck_options[1][57] > 1)
		printf("\n DEV TUCK: First bad site (%c) at row=%d, column=%d.", badletr, bad_m+1, bad_n+1);

	mha_writeback(tuckarray, lcl_align2D, tuck_options);
	mha_writecons(tuckarray, lcl_align2D, tuck_options);

	if (lcl_align2D[bad_m][bad_n-1] == blank) {
		for (j = bad_n; lcl_align2D[bad_m-1][j] != '/' && j < MAXROW; j++) {
			;
		}	
		j_stop = j-1;
		if (lcl_align2D[bad_m-1][j_stop] == badletr &&
			isalpha(lcl_align2D[bad_m][j_stop]) && 
			lcl_align2D[bad_m][j_stop] == badletr) {
				tuckcall = -1;									/* CAN TUCK LETTER ABOVE */
		} 
		else if (nuctransit && ((letr=lcl_align2D[MAXROW][j_stop])=='R' || letr=='Y')) {
			if      (letr=='R' && (badletr=='A' || badletr=='G'))
				tuckcall = -1;
			else if (letr=='Y' && (badletr=='C' || badletr=='T'))
				tuckcall = -1;
		}
	}
	if (tuckcall != -1) {
		for (i=0; i<bad_m; i++) {
			if (isalpha(lcl_align2D[i][bad_n]) && lcl_align2D[i][bad_n+1]=='/') {
				badletr = lcl_align2D[i][bad_n];
				bad_m = i;
				tuckcall = 1;
				break;
			}
		}
		if (tuckcall && lcl_align2D[bad_m+1][bad_n-1] != blank) {
			for (j=bad_n-2; j >= 0; j--) {
				if (lcl_align2D[bad_m+1][j] == blank)
					break;
			}
			j_stop = j+1;
			if (lcl_align2D[bad_m+1][j_stop]==badletr && lcl_align2D[bad_m+1][j]==blank) {
				tuckcall = 1;									/* CAN TUCK LETTER BELOW */ 
			}
			else if (nuctransit && ((letr=lcl_align2D[MAXROW][j_stop])=='R' || letr=='Y')) {
				if      (letr=='R' && (badletr=='A' || badletr=='G')) {
					tuckcall = 1;
				}
				else if (letr=='Y' && (badletr=='C' || badletr=='T')) {
					tuckcall = 1;
				}
				else
					tuckcall = 0;
			}
			else {
				return(0);
			}
		}
		else {
			tuckcall = 0;
			return(0);
		}
	}

	if (tuckcall==-1) {
		printf("\n Winner! Winner! Chicken dinner! Tuck up!");
		if (pushdown(lcl_align2D, bad_m, bad_n+1, tuck_options)) {
			mha_writeback(lcl_align2D, tuckarray, tuck_options);
			mha_writecons(lcl_align2D, tuckarray, tuck_options);
			for (j=0; j < j_stop; j++) {
				tuckarray[bad_m][j  ] = blank;
			}
				tuckarray[bad_m][j  ] = badletr;
				tuckarray[bad_m][j+1] = '/';
				tuckarray[bad_m][j+2] = '\0';
			return(1);			/* RETURN THAT TUCKSENSE SHOULD BE INVOKED AGAIN */
		}
		else 
			return(0);
	}
	else if (tuckcall==1) {
		printf("\n Winner! Winner! Chicken dinner! Tuck down!");
		if (pushdown(lcl_align2D, bad_m+1, 0, tuck_options)) {
			mha_writeback(lcl_align2D, tuckarray, tuck_options);
			mha_writecons(lcl_align2D, tuckarray, tuck_options);
				tuckarray[bad_m  ][bad_n   ] = '/';
				tuckarray[bad_m  ][bad_n+1 ] = '\0';
				tuckarray[bad_m+1][j_stop+1] = '/';
				tuckarray[bad_m+1][j_stop+2] = '\0';
			return(1);			/* RETURN THAT TUCKSENSE SHOULD BE INVOKED AGAIN */
		}
		else
			return(0);
	}
	else
		return(0);			/* RETURN NOTHING TUCKABLE SENSED */
}

/**************************************************************** END OF MHA fUNCTIONS **********************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal			*/

