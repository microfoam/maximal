/*******************************************************************************************************/
/***** The program "maximal" is a micro homology alignment (MHA) program.                          *****/
/***** Designed and written by Dr. Albert J. Erives. 2017-2020. AGPL-3.0.                          *****/
/***** Code repository located at https://github.com/microfoam/maximal. Licensed under AGPL-3.0.   *****/
/***** This program renders a 1-D DNA sequence into a 2-D self-alignment to rescue micro-paralogy. *****/
/*******************************************************************************************************/

#include <stdio.h>			/********************/
#include <math.h>			/* mathematical!    */
#include <ctype.h>  		/* isalpha() 		*/
#include <stdlib.h> 		/* system()  		*/
#include <string.h> 		/* strlen()  		*/
#include <time.h>			/* difftime() 		*/
#include <signal.h>			/* signal() 		*/
#include <unistd.h>			/* signal()-related */
/*******************************************************************************************************/
#include "microhomology.h"	/* maximal header: main header file                                        */
#include "microhom-devl.h"	/* maximal header: program development code, testing, and evaluation       */
#include "microhom-skor.h"	/* maximal header: alignment scoring definitions and functions             */
#include "microhom-tela.h"	/* maximal header: verb_tela() functions and coord struct-related          */
#include "microhom-cinc.h"	/* maximal header: original cinch modules except for cinch-t (in main)     */
/*******************************************************************************************************/

int main(int argc, char *argv[])
{
	int match      = MATCH;
	int transition = TRANSITION;			
	int mismatch   = MISMATCH;
	char version[] = "4.27";				/* current development version number */

	int c=0, f=0, i=0; 
	short unsigned int continue_flag=0;

	int j=0, h=0, k=0; 
	short unsigned int imperfect_TR=0; 

	int l=0, m=0, n=0; 
	short unsigned int Aimperfect_TR=0; 

	int o=0, p=0, q=0;
	short unsigned int nuctransit=0; 

	int reps=0, r=0, z=0;	
	short unsigned int seqtype=0; 

	int DTHR = 100;
	int m2Da_height = 1;	
	int tuck;
	short unsigned int conflict_flag=0;

	int number=0;
	int sumspan=0;
	int homopoly_flag=0;
	short unsigned int msa = 0;	

	int homopolyend_flag=0, overslip=0, TRcheck = 0;
	unsigned int FY_size = 100;				/* DEFAULT SIZE OF FISHER-YATES RANDOMIZED STRING */

	int a2D_n = 0;							/* NUMBER INDEX OF n FOR a2D_n */
	int lenseq = 0;
	int citwidth = 0;
	char blank = '.';						/* DEFAULT BLANK CHARACTER FOR 2-D MHA. FULLSTOP = 46 */
											/* NEEDS TO BE SET IN options[][62] array */
	int badslip_type = 0;
	int scooch = 0;
	int alt_k = 0;
	char ch = blank;

	float ratio1 = 1;						/* WIDTH CINCH RATIO (W.C.R.) post cinch-d, pre relax-2D 	*/
	float ratio2 = 1;						/* WIDTH CINCH RATIO (W.C.R.) post relax-2D 				*/

	int passQ[16] = {0};        			/* PASS QUALITY */
	int passR[16] = {0};        			/* PASS RUNS */
 	char m2Dalig[MAXROW+1][MAXROW] = {{0}};			
	char letr_unit[8] = {0};				/* UNIT STRING: "bp" FOR DNA, "nt" FOR RNA, 'aa' FOR PROTEINS, 'ch' FOR ALL OTHER */
	char cycle[WIDTH+1];		/* THIS ARRAY HOLDS THE CYCLIC PATTERN OF TRs W/ >2 UNITS */
	char numstring[8] = {0};
	char Seq_head[100] = {0};	/* FASTA HEADER */
	char Seq_i[MAXROW] = "TGTGTGAGTGAnnnnnnTGTGTGAGTGAGnnnnnTGTGTGAGTGAGTGAnnTGTGTGAGTGAGTGAGT"; 	/* INPUT SEQUENCE W/ DEFAULT */
	char Seq_r[MAXROW] = {0}; 	/* RANDOMIZED SEQUENCE */
	char *Seq = Seq_i;			/* POINTER TO INPUT SEQUENCE */

	FILE *file_ptr;
	FILE *fp_cons;							/* FILE FOR CONSENSUS STRING Surf_barrels.log */
  	FILE *fp_msa;							/* FILE FOR MHA MSA "TUBES.mha" */
	FILE *fp_tricksy;						/* DEV. FILE FOR IMPERFECT 2-D ALIGNED STRINGS waves/foam_and_chowder.mha */

	time_t lcl_time = time(NULL);			/* START TIME */
	char time0[26];							/* START TIME STRING */
	strcpy(time0,ctime(&lcl_time));			/* TEXT-READABLE START TIME */

	signal(SIGINT, signal_callback_handler);	/*  2 */
	signal(SIGBUS, signal_callback_handler);	/* 10 */
	signal(SIGSEGV, signal_callback_handler);	/* 11 */

	/* IS THERE A FILE NAME ARGUMENT? */
	for (i = 1; i < argc; i++) {
		if (*argv[i] != '-') {
			if (*argv[i] > '0' && *argv[i] <= '9') {
				strcpy(numstring, argv[i]);
				l = strlen(numstring);
				for (r=0; r<l; r++) {
					if (!isdigit(numstring[r])) {
						warnhead('d');
						printf("Command line arguments starting with numbers need to be strings with digits-only.\n\n");
						usage(version, FY_size);
						exit(EXIT_ERROR);
					}
					else
						number = number * 10 + (numstring[r]-'0');
				}
				printf("\n %2d. Reading an argument as the number %d. ", j, number);
			}
			else if (strcmp(argv[i],"TUBES.mha")) {		/* strcmp EVALUATES TO 0 ONLY IF STRINGS ARE THE SAME */
				/* USING j BELOW AS COUNTER TO NUMBER OF TIMES USER SPECFIES DIFFERENT SEQUENCES */

				if ( (file_ptr = fopen(argv[i], "r") ) == NULL) {
					printf("\n %2d. Error opening file '%s'. Exiting now.\n\n", ++j, argv[i]);
					exit(EXIT_ERROR);
				}
				else {
					fseek(file_ptr, 0, SEEK_END);
					lenseq = ftell(file_ptr);
					if (lenseq > MAXROW) {		/* LAST ROW OF array2D WILL STORE CONSENSUS, SO NEED TO KEEP CLEAR */
						printf("\n %2d. Sequence (length %d) from file '%s' > MAXROW limit (%d) by %d.", ++j, lenseq, argv[i], MAXROW, lenseq-MAXROW+1);
						printf("\n %2d. Exiting now. For help enter './maximal -h'. \n\n", ++j);
						exit(EXIT_ERROR);
					}
					else fseek(file_ptr, 0, SEEK_SET);
				 
					fscanf(file_ptr, "%[^!]c", Seq_i);
					fclose (file_ptr);
					strcpy(file_name, argv[i]);
					printf("\n %2d. Detecting sequence from file: '%s'.", ++j, file_name);
	
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
					else 
						strcpy(Seq_head, "input sequence");
				}
			}
			else {
				if ( (file_ptr = fopen(argv[i], "r") ) == NULL) {
					printf("\n*%2d. Error opening supporting file '%s'.", j, argv[i]);
				}
				else if (msa == 0){
					msa = 1;		/* USING THIS SLOT TO STORE BIT VALUE INDICATING TUBES.mha CO-INPUT */

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
						for (m=0; m2Dalig[m][0]!='\0'; m++)	{
							printf("%s\n", m2Dalig[m]);
						}
						m2Da_height = m;
					}
				}
			}
		}	/* END OF IF ARGV[I] != '-' */
	}	/* END OF FOR i = 1, i < argc, i++ */

	char recovered[MAXROW] = {0};
	short unsigned int go_flag=0, cycle_flag=0;			/* USE THIS TYPE FOR TRUE BIT FLAG VARIABLES */
	unsigned int recovery_flag = 0;
	int relax_length=0;			/* FOR USE WITH relax_2D CALL */
	int intraTR_reps_tot = 0; 	/* STORES INITIAL RETURN VALUE FROM cinch-d() */
	int intraTR_reps = 0;	 	/* STORES CURRENT RETURN VALUE FROM cinch-d() */
	int Did = 0;				/* Counter for identity (id) diagonal */
	int Dtr = 0;				/* Counter for tandem repeat (tr) diagonal */
	int Atr = 0;				/* Counter for additional repeats on the same diagonal */
	int row = 0;				/* Counter for row number in align2D box */
	int recslips= 0;			/* Counter of recent slips in region of first TR unit, derived from tela[].r */
	int slips[WIDTH+1] = {0};	/* Array of counters for unique slips of WIDTH x	*/
	int opt;					/* opt IS CASE OPTION VARIABLE FOR SETTING options ARRAY (opt_ VARs) */ 
    char optR = options[1][27];
	options[1][58] = 80;    /* opt_w IS MAXIMUM WIDTH OF 1D/2D-PRINTED BLOCK, SCREEN WRAP LENGTH */
	int blocks;				/* Number of blocks for 1D output print */

	if (argc == 1) {
	 	system("clear"); 
		usage(version, FY_size);
		exit(EXIT_EARLY);
	}

	/**************************************/
	/* SET OPTIONS FROM ARGUMENTS  ********/
	while (--argc > 0 && (*++argv)[0] == '-') {
		while ((opt = *++argv[0])) {
			switch (opt) {
			case 'c':						/* SHOW BASE 62 CODE */
					print_base62_table();
					exit(EXIT_EARLY);
					break;
			case 'd':						/* OPTION TO SKIP CINCH-D CINCHING */
					options[0][39] = 1;		/* opt_d ON		*/
					break;
			case 'f':						/* OPTION TO SHOW FOAM-FREE SEGMENTS BELOW CONSENSUS ROW*/
					options[0][48] = 1;		/* opt_m ON 	*/
					options[0][41] = 1;		/* opt_f ON		*/
					break;
			case 'g':
					options[0][42] = 1;
					--options[1][48];		/* opt_g GEL-UP, COUNTERACT DEFAULT MELTAGE */
					break;
			case 'h':						/* OPTION TO SHOW HELP */
					options[0][43] = 1;		/* opt_h ON 	*/
					break;
			case 'k':						/* OPTION TO SHOW k-MER COUNTS */
					options[0][46] = 1;		/* opt_k ON   */
				/*  options[1][46] RESERVED FOR RECORDING LARGEST CINCH-T k-MER UNIT SIZE */
					break;
			case 'l':						/* OPTION TO SHOW SLIP LOCATIONS IN 1D SEQUENCE */
					options[0][47] = 1;		/* opt_l ON   */
					break;
			case 'm':						/* OPTION TO SPLIT, OPEN, AND MELT */
					options[0][48] = 1;		/* opt_m ON 	*/
					++options[1][48];		/* ++opt_m VAL	*/
					break;
			case 'n':						/* OPTION TO NOT DO RELAX-2D PASS */
					options[0][49] = 1;		/* opt_n ON   */
					break;
			case 'o':						/* OPTION TO PRINT ORIGINAL STRING UNFORMATTED */
					options[0][50] = 1;		/* opt_o ON   */
				/*  options[1][50] RESERVED FOR RECORDING CUMULATIVE BADSLIP TYPE PER RUN */
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
			case 't':						/* OPTION TO SKIP CINCH-T */
					options[0][55] = 1;		/* opt_s ON   */
					break;
			case 'u':						/* OPTION TO PRINT UNWRAPPED WHERE opt_w EQUALS lenseq */
					options[0][56] = 1;		/* opt_u ON 	*/
					++options[1][56];		/* ++opt_u VAL	*/
					break;
			case 'v':						/* OPTION FOR VERBOSITY */
					options[0][57] = 1;		/* opt_v ON 	*/
					++options[1][57];		/* opt_v INCREMENTED: 1-2 FOR USERS, 3-4 FOR DEVELOPERS */
											/*   1=EXTRA INFO; 2=BUFFER; 3=DEV-ACTIVE; 4=DEV-LEGACY */
					++options[0][18];		/* opt_I */
					++options[0][27];		/* opt_R */
					options[0][15]=options[0][21] = 1; /* opt_F, opt_L */
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
			case 'C':						/* OPTION TO USE REVERSE COMPLEMENT */
					options[0][12] = 1;		/* opt_C ON		*/
					break;
			case 'D':
				   	options[0][13] = 1;		/* opt_D: TURN DEV_PROMPTS ON */
					/* options[1][13] IS RESERVED FOR STORING SEQUENCE TYPE (DNA, RNA, PROTEIN, BABYLONIAN, etc.) */
					break;
			case 'F':						/* OPTION TO USE BLANK FILL CHAR W/ SCRIMMAGELINE */
					options[0][15] = 1;		/* opt_F ON: USE POST-TERMINATION FILLER	*/
					break;
			case 'H':						/* OPTION TO SHOW HELP */
					options[0][43] = 1;		/* opt_h ON 	*/
					options[0][17] = 1;		/* opt_H ON 	*/
				/*  options[1][17] = 1 RESERVED FOR 2D-ALIGNMENT HEIGHT */	
					break;
			case 'I':						/* OPTION TO SHOW iNITIAL PASSES */
					++options[0][18];		/* opt_I = 1 is ON; opt_I > 1 is O-F-F (lets you toggle off after option v) */
				/*  options[1][18]            RESERVED FOR COUNTING iNITIAL PASSES. */
					break;
			case 'K':						/* OPTION TO SHOW CONSENSUS ROW */
					options[0][20] = 1;		/* opt_K ON */
					break;
			case 'L':						/* OPTION TO SHOW POSITIONS AT END OF LINES */
					options[0][21] = 1;		/* opt_L ON */
					break;
			case 'M':						/* OPTION TO DOUBLE LONG HOMOMONO WRAP */
					options[0][22] = 1;									/* opt_M ON		*/
					options[1][22] = options[1][22] + options[1][22];   /* DOUBLE opt_M  mwrap */
					break;
			case 'O':						/* OPTION TO OUTPUT CONSENSUS FILE */
					++options[0][24];		/* opt_O ON 	*/
					break;
			case 'P':						/* OPTION TO PRINT PATH BOX */
					options[0][25] = 1;		/* opt_P ON 	*/
					break;
			case 'R':						/* OPTION TO PRINT RECOVERED 1-D SEQUENCE FROM LAST 2-D */
					++options[0][27];		/* opt_R = 1 is ON; opt_R > 1 is O-F-F (lets you toggle off after option v) */
					break;
			case 'T':
				   	options[0][29] = 1;		/* opt_T: SHOW DTHR_lookup values (later below). Will exit early. */
					break;
			case 'X':						/* OPTION TO SCRAMBLE SEQUENCE           */
					options[0][33] = 1;		/* opt_X ON 							 */
					++options[1][33];		/* INCREMENT opt_X to desired level 	 */
					break;					/*  X = rand() cheese, XX = FISHER-YATES */
			case 'Y':						/* OPTION TO SPECIFY FY_size	*/
					options[0][34] = 1;
					if (number) {
						if (number < MAXROW) {
							FY_size = number;
						}
						else {
							warnhead('Y');
							printf("Option 'Y' to specify FY_size, but it must be less than %d.", MAXROW);
						}
					}
					else {
						warnhead('Y');
						printf("Option 'Y' to specify FY_size, but none specified; using default FY_size = %d.", FY_size);
					}
					break;
			default:
					printf("maximal: Illegal option %c\n", opt);
					argc = 0;
					usage(version, FY_size);
					exit(1);
					break;
			} /* END SWITCH opt			*/
		} /* END WHILE opt			*/
	} /* END WHILE argc argv	*/

	/* IF CERTAIN OPTIONS ARE ON, SKIP RELAX-2D */
	if (options[0][39] || options[0][24]) {
		options[0][49] = 1;		/* opt_n NO RELAX 2-D */
	}

	/*******************************************************/
	/* BEGIN OUTPUT OF MAXIMAL PROGRAM  ********************/
	if (options[0][17] || options[0][43]) {		/* opt_H SHOW USAGE */
	 	system("clear"); 
		usage(version, FY_size);
		exit(2);
	}

	/* IF OPTION I TOGGLED BACK O-F-F FROM OPTION v (VERBOSE), then set opt_I back to zero for screen reporting purposes */
	if (options[0][18] > 1) {
		options[0][18] = 0;
	}

	/* OPTION TO APPEND TO GROWING 2-D MSA FILE */
	if (options[0][24] >= 2) {
		options[0][24] = 2;
		options[0][49] = 1;		/* opt_n NO RELAX 2-D */
	}

	/* IF OPTION R TOGGLED BACK O-F-F FROM OPTION v (VERBOSE), then set opt_R back to zero for screen reporting purposes */
	if (options[0][27] % 2 == 0) {		
		options[0][27] = 0;
	}

	if (j > 1) {
		printf("\n");
		warnhead('S');
		printf("Many sequences specfied. Using last sequence.\n");
	}
	else if (j == 0) {
		printf("\n");
		warnhead('S');
		printf("No sequences specfied. Using example sequence.\n");
	}
	printf("\n");
		
	mha_head(options[1][58]);
	printf("micro homology alignment (MHA) -");
	for (i = 10; i < 36; i++) {			/* UPPER-CASE LETTER OPTIONS */
		if (options[0][i] > 0)
			printf("%c", (char) options[3][i]);
	}
	for (i = 36; i < 62; i++) {			/* LOWER-CASE LETTER OPTIONS */
		if (options[0][i])
			printf("%c", (char) options[3][i]);
	}
	if (options[0][11] > 1) {			/* opt_B */
		printf(" -B%ld", options[0][11]);
	}
	printf(" -M%ld", options[1][22]);	/* opt_M */
	if (options[0][33]) {				/* opt_X */
		if (options[1][33]) {
			printf(" -X [USE RANDOMIZED SEQUENCE]");
		}
		else {
			printf(" -XX [USE FISHER-YATES RANDOMIZED SEQUENCE]");
		}
	}
	if (options[1][56] > 1) {			/* opt_u */
		printf(" -u%ld", options[1][56]);
	}
	printf(" (version %s)\n", version);
	printf("%s", time0);

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
	lenseq = strlen(Seq) - 1;
	if (lenseq > MAXROW) {		/* LAST ROW OF array2D WILL STORE CONSENSUS, SO NEED TO KEEP CLEAR */
		warnhead('M');
		printf("Sequence (length %d) exceeds MAXROW size limit (%d) by %d.\n\n", lenseq, MAXROW, lenseq-MAXROW+1);
		exit(1);
	}

	++options[1][18];	/* INCREMENT opt_I COUNTER (pass_num) */

	options[1][0] = options[1][32] = lenseq;	/* ASSIGN CINCH-WIDTH TO HISTORY [0--9] AND CURRENT [32]	*/

	if ((i=(int) options[1][56]) >= 1) {					/* WILL CAUSE OUTPUT TO NOT BE WRAPPED (opt_u EQUALS 1),	*/ 
		options[1][58] = (lenseq/i);	/*  OR WRAPPED INTO opt_u NUMBER OF BLOCKS.				    */
	}

	/* OPTION TO PRINT ORIGINAL STRING IN BLOCKS *******************************/
	if (options[0][50]) {	/* opt_o */
		blocks = count_wrap_blocks(lenseq, options[1][58]);

		mha_head(lenseq);
		if (strlen(Seq_head) > 0)
			printf(">%s\n", Seq_head);
		printf("\"");
		for (j = 0; j < blocks; j++) {
			for(n = j * ((int) options[1][58]); (n < (j+1) * ((int) options[1][58])) && (Seq[n] != '\0'); n++) {
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
	options[1][13] = seqtype = cleanseq(Seq);	/* opt_D: STORES SEQTYPE: 1=DNA, 2=RNA, 3=PROTEIN, 0=OTHER */
	lenseq = strlen(Seq);
	options[1][1] = options[1][32] = lenseq;	/* ASSIGN CINCH-WIDTH TO HISTORY [0--9] AND CURRENT [32]	*/

	if (options[0][29]) {	/* opt_D: SHOW DTHR VALUES */
		show_DTHR_table();
	}
	else if (seqtype == 1) {
		nuctransit = 1;
		strcpy(letr_unit, "bp");	/* BASE PAIRS */
	}
	else if (seqtype == 2)
		strcpy(letr_unit, "nt");	/* NUCLEOTIDES */
	else if (seqtype == 3)
		strcpy(letr_unit, "aa");	/* AMINO ACIDS */
	else if (seqtype == 0)
		strcpy(letr_unit, "ch");	/* OTHER */

	if (nuctransit && options[0][12]) {		/* opt_c USE REVERSE COMPLEMENT */
		for (i=1; i<=lenseq; i++) {
			if ( (ch=Seq[lenseq-i]) == 'A')
				Seq_r[i-1] = 'T';
			else if (ch == 'T')
				Seq_r[i-1] = 'A';
			else if (ch == 'C')
				Seq_r[i-1] = 'G';
			else if (ch == 'G')
				Seq_r[i-1] = 'C';
			else
				Seq_r[i-1] = 'n';
		}
		Seq_r[lenseq] = '>';
		strcpy(Seq, Seq_r);
		if (options[1][57]) {			/* opt_v VERBOSITY */
			printf("\nReverse Complement: \"");
			for (i = 0; Seq_r[i] != '\0'; i++)
				printf("%c", Seq_r[i]);
			printf("\"\n");
		}
		Seq = Seq_r;
		if (options[0][33])
			options[0][33] = options[1][33] = 0;		/* USE RANDOMIZED SEQUENCE DISABLED IF R.C. REQUESTED */
	}
	else if (options[0][33]) {		/* USE RANDOMIZED SEQUENCE */
		strcpy(Seq_r, Seq);
		srand(time(0));

		if (options[1][33] == 1) {
			mha_randomize1(Seq_r);
		}
		else {
			mha_randomize2(Seq_r, FY_size);
			lenseq = options[1][1] = FY_size;
		}

		printf("\nRandomized sequence: \"");
		for (i = 0; Seq_r[i] != '\0'; i++) {
			printf("%c", Seq_r[i]);
		}
		printf("\"\n");
		Seq = Seq_r;
	}

	++options[1][18];		/* INCREMENT ++pass_num */
	passQ[0] = passQ[1] = 1000;

	Seq[lenseq] = tela[lenseq].t = tela[lenseq].c = '>';

	for (i = 0; i <= lenseq; i++) {
		tela[i].x = tela[i].X = i;
		tela[i].t = tela[i].c = Seq[i];
		tela[i].cyc_o = tela[i].echoes = blank;

		if (nuctransit) {
			if 		(tela[i].c=='A' || tela[i].c=='G')
				tela[i].e = 'R';
			else if (tela[i].c=='C' || tela[i].c=='T')
				tela[i].e = 'Y';
			else
				tela[i].e = tela[i].c;
		}
		else
			tela[i].e = tela[i].c;
	}

	/************************************************************************************/
	/* INITIALIZE PATHBOX FOR SELF-MHA  *************************************************/

	homopoly_flag = 1;								/* FOR LONG HOMOPOLYMERIC RUN CASE **/
	homopolyend_flag = 0;

	for (n = 0; tela[n].c != '\0'; n++)				/* SET IDENTITY LINE ****************/
		pathbox[n][n] = MATCH;

	for (n = 0; tela[n].c != '\0'; n++) {
		if  (n <= WIDTH) {							/* SET VALUES FOR 1ST WIDTH COLS ****/
				for (m = 0; m < n; m++){
					if (seqtype && tela[n].c == 'n')			 
						pathbox[m][n] = mismatch;   /* SPECIAL TREATMENT FOR 'n' IN DNA**/
					else if (tela[n].c == tela[m].c)
						pathbox[m][n] = MATCH;		/* MATCH ****************************/
					else if (nuctransit) {			/* IF DNA AND CHECKING FOR TRANSITIONS */
						if      (tela[n].c == 'A' && tela[m].c == 'G')
							pathbox[m][n] = transition;  
						else if (tela[n].c == 'G' && tela[m].c == 'A')
							pathbox[m][n] = transition;   
						else if (tela[n].c == 'C' && tela[m].c == 'T')
							pathbox[m][n] = transition;   
						else if (tela[n].c == 'T' && tela[m].c == 'C')
							pathbox[m][n] = transition;   
						else
							pathbox[m][n] = mismatch;   /* MISMATCH IF NO TRANSITION ****/
					}	
					else 
						pathbox[m][n] = mismatch;   /* MISMATCH *************************/
				}
		}
		else {										/* SET VALUES FOR REST OF COLUMNS ***/
				for (m = n-WIDTH; m < n+1; m++){	/*  WITHIN BAND WIDTH ***************/
					if (seqtype && tela[n].c == 'n')			  
						pathbox[m][n] = mismatch;   /* DUE TO NUCLEOTIDE AMBIGUITY ******/
					else if (tela[n].c == tela[m].c)
						pathbox[m][n] = MATCH;		/* MATCH ****************************/
					else if (nuctransit) {			/* IF DNA AND CHECKING FOR TRANSITIONS */
						if      (tela[n].c == 'A' && tela[m].c == 'G')
							pathbox[m][n] = transition;  
						else if (tela[n].c == 'G' && tela[m].c == 'A')
							pathbox[m][n] = transition;   
						else if (tela[n].c == 'C' && tela[m].c == 'T')
							pathbox[m][n] = transition;   
						else if (tela[n].c == 'T' && tela[m].c == 'C')
							pathbox[m][n] = transition;   
						else
							pathbox[m][n] = mismatch;   /* MISMATCH IF NO TRANSITION ****/
					}	
					else 
						pathbox[m][n] = mismatch;   /* MISMATCH *************************/
				}
		}

		if (tela[n].c != tela[n+1].c)	{					/* HERE CHECK FOR HOMOPOLYMERIC RUN */
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
	if (options[0][25]) {	/* opt_P */
		blocks = count_wrap_blocks(lenseq, options[1][58]);

		printf("\nPATHBOX FILL-IN PASS (length = width = %d)\n\n", lenseq);
		for (j = 0; j < blocks; j++) {
			if (blocks != 1)
				print_blockhead(j+1, blocks);
			line_end(PATHBOXHEAD, 9, 9);	
			for(n = j * options[1][58]; (n < (j+1) * options[1][58]) && (tela[n].c != '\0') && tela[n].c != '>'; n++) 
				printf("%2c", tela[n].c);
			printf("\n");
			for(m = j * options[1][58]; (m < (j+1) * options[1][58]) && (tela[m].c != '\0') && tela[m].c != '>'; m++) {
				printf("%4d. %c ", m+1, tela[m].c);
					for (n = j * options[1][58]; (n < (j+1) * options[1][58]) && (tela[n].c != '\0') && tela[n].c != '>'; n++) {
						if (m > n) {
							if (pathbox[m][n])
								printf("%2d", pathbox[m][n]);
							else 
								printf("%2c", blank);
						}
						else if (n-m <= WIDTH)
							printf("%2d", pathbox[m][n]);
						else 
							printf("%2c", blank);
				}
				printf("\n");
	   		 }
		}
	} /* END OF OPTION TO PRINT PATHBOX */

	/*********************************************************/
	/*        USE PATHBOX TO BUILD FIRST 2-D ALIGNMENT       */
	/*        	          [cinch_t BEGINS]                   */
	citwidth = lenseq;					/* START AT 1D LENGTH AND CONTRACT DURING CINCH-T ("cit")     */
	a2D_n = row = 0;					
	align2D[row][a2D_n++] = tela[0].c;	/* FOR n=0, ENTER VALUE AT IDENTITY DIAGONAL, INCREMENT INDEX */

	mark_tela();			/* WILL MARK ALL TRs WITHOUT CINCHING AND RECORD IN tela[].all_k, all_r, all_S, all_L/R */
/*	dev_prompt(MAIN,__LINE__,file_name);
*/
	for (n = 1; n<=lenseq; ) {
		/* FOR COLUMN n LOOP 1/3 */
		if (!options[0][55]) {		/* SKIP TO NEXT MARKED TR */	
			while (!(tela[n].all_S) && n!=lenseq) {
				assign_tela(n++, row, a2D_n++, ONE);	/* MODES ZERO O-F-F, NON-ZERO ASSIGN  */
			}
		}
		else {						/* ELSE opt_t: SKIP CINCH-T */ 
			strcpy(align2D[0],Seq);
			break;
		}

		/* FOR COLUMN n LOOP 2/3: SKIP PRESENT TR IF CONFLICT AND CAN CYCLE WITH SAME SCORE */
		if (tela[n].all_L && tela[n].all_S == tela[n+1].all_S && !tela[n+1].all_L && tela[(tela[n].all_L)].cyc_o == 'x') 
			assign_tela(n++, row, a2D_n++, ONE);	/* MODES ZERO O-F-F, NON-ZERO ASSIGN  */

		/* FOR COLUMN n LOOP 3/3 */
		for (m = 0; m < n; m++) {
			/* FOR ROW m LOOP 1/6: UPDATE VAR CITWIDTH AT END */
			if (tela[n].c == '>') {
				citwidth = a2D_n;
				align2D[row+1][0] = '\0';
				options[1][17] = row+1;
			}
			else if (tela[n].cyc_o == 'o') {	/* IF THIS POSITION HAD A BIGGER k-MER SQUASHED HERE (E.G., FOR LATER CYCLING) */
				assign_tela(n++, row, a2D_n++, ONE);	/* MODES ZERO O-F-F, NON-ZERO ASSIGN  */
				break;	/* GO TO NEXT n */
			}

			/* FOR ROW m LOOP 2/6: SLIDE DOWN TO ROW WITHIN POPULATED HEMIDIAGONAL */
			if (n-m > WIDTH+1) 
				m = n-WIDTH;

			/* FOR ROW m LOOP 3/6: SET K-MER SIZE AND DTHR SCORE THRESHOLD */
			k = n-m;
			if (nuctransit) {
				DTHR = score_DTHR(k);
				imperfect_TR = 0;
			}

			/* FOR ROW m LOOP 4/6: SKIP k=ONE */
			if (k == 1) {	
				assign_tela(n++, row, a2D_n++, ONE);		/* MODES ZERO O-F-F, NON-ZERO ASSIGN  */
				break;	/* GO TO NEXT n */
			}

			/* FOR ROW m LOOP 5/6: SET HOMOPOLYMER RUN STATUS UNKNOWN; USED TO RULE OUT k>1 MONONUCLEOTIDE "REPEATS" */
			homopoly_flag = 2;
			if (tela[n].c != tela[n-1].c)
				homopoly_flag = 0;

			/* FOR ROW m LOOP 6/6: START COUNTING SCORE IF PATHBOX POSITION HAS VALUE > MISMATCH */
			if (pathbox[m][n] > mismatch && n+k <= lenseq) {
				Dtr = 0;

				/* IF SUMMING PATHBOX DIAGONAL 1/: COMPUTE SCORES OF IDENTITY LINE AND 1st REPEAT DIAGONAL*/
				Did = k*MATCH;
				for (i = m; i < n; i++) {
					if (pathbox[i][i+k] == mismatch) {	/* STOP SHORT IF MISMATCH IS FOUND 		 		*/
						Dtr =  0;						/* B/C CURRENTLY ONLY CONSIDERING TRANSITIONS 	*/
						break;							
					}
					else
						Dtr = Dtr + pathbox[i][i+k];	/* COMPUTE SUM OF TANDEM REPEAT UNIT LINE */

					/* SET HOMOPOLYMERIC RUN BIT TO FALSE IF NOT A POSSIBILITY */
					if (homopoly_flag && i > m && tela[i].c != tela[i-1].c)
						homopoly_flag = 0;
				}

				/* IF SUMMING PATHBOX DIAGONAL 2/: SET HOMOPOLYMERIC RUN BIT TO TRUE IF DETECTED 	*/
				if (homopoly_flag && i == n) {
					homopoly_flag = 1;				/* BIT IS THERE IF NEEDED BEYOND BREAK. 		*/
					assign_tela(n++, row, a2D_n++, ONE);	/* MODES ZERO O-F-F, NON-ZERO ASSIGN  */
					Dtr = 0;
					break;							/* GO TO NEXT n */
				}

				/* IF SUMMING PATHBOX DIAGONAL 3/: IF CONSIDERING NUCL. TRANSITIONS AS PARTIAL MATCHES */
				if (nuctransit && Dtr && Dtr!=Did) { 
					if (k>PISO && 100*Dtr/Did > DTHR)	{	
						imperfect_TR = 1;		/* CALLING TR W/ TRANSITIONS FOR n BLOCK VS m BLOCK */
						tela[n].Dtr = Dtr;
					}
					else {
						Dtr = 0;
					}
				} 

				if (Dtr && dev_print(MAIN,__LINE__)) {
					printf("cinch_t evaluating k=%d-mer all_r=%d at n=%d, imperfect_TR=%d.", k, tela[n].all_r, n, imperfect_TR);
				}

				/* IF SUMMING PATHBOX DIAGONAL 4/: FIND AND STORE POSITION OF LEFT-MOST OVERLAPPING TRs */
				/* SKIP CINCH IF CAN AVOID CONFLICT WITH A LATER CYCLE & CYCLING PREVIOUS TR NOT AN OPTION */
				if (Dtr && (Dtr==Did || imperfect_TR)) {
					for (i = n-1; i > 1; i--) {
						if (tela[i].r && i + tela[i].k * (tela[i].r - 1) > m) {
							tela[n].X = i;		/* UPDATE LEFT-MOST OVERLAPPING TR */

							if (dev_print(MAIN,__LINE__)) {
								printf("Setting tela[%d].X to i=%d. m=%d, k=%d, Dtr=%d, imperfect_TR=%d.", n, i, m, k, Dtr, imperfect_TR);
							}
						}
					} 
					if (!imperfect_TR && (l=tela[n].X) > m && l + span_rk(l) <= n)
						;
					else if (tela[n].X != n) {
						q = tela[n].X;
						if (tela[q+1].cyc_o != 'o') {
							alt_k = tela[q].k;
							j = q + alt_k*(tela[q].r - 1);
							for (l = n; l+k <= lenseq && l+k <= m+WIDTH; l++) {
								if (tela[l].c == tela[l+k].c) { 
									if (j <= l+1) {
										pull_tela(n);
										Dtr = imperfect_TR = 0;
										assign_tela(n++, row, a2D_n++, ONE);	/* MODES ZERO O-F-F, NON-ZERO ASSIGN  */
										tela[n].cyc_o = 'o';
										if (dev_print(MAIN,__LINE__)) {
											printf("Skipping cinch here to use later cycling frame at %d with print_tela:.", l+1);
											print_tela(prtela_A, prtela_B);
										}
										break;
									}
								}
								else break;
							}
						}
					}
				}

				/* IF SUMMING PATHBOX DIAGONAL 5/: SKIP CINCH IF IMPERFECT WHILE CONTAINING PERFECT TRs INSIDE */
				if (imperfect_TR) {
					for (l=k/2; l>1; l--) {
						for (i=0; i<l; i++) {
							if (tela[n-l+i].c != tela[n+i].c)
								break;	
						}
						if (i==l) {
							Dtr = imperfect_TR = 0;
							break;
						}
					}
				}

				/* IF SUMMING PATHBOX DIAGONAL 6/: START COUNTING REPEATS */
				if (Dtr && (Dtr==Did || imperfect_TR)) {
					/* COUNT NUMBER OF REPEATS ALBERT-STYLE */
					TRcheck = 1;
					reps = 1;
					while (TRcheck) {
						Atr = Aimperfect_TR = 0;
						if (m + (reps+1)*k >= lenseq) { 
							Atr = 0;
							break;
						}

						if (nuctransit) {
							Atr = score_kmer(n+k*reps,k,TWO);
							if (Atr!=Did && (100*Atr)/Did > DTHR) {
								Aimperfect_TR = 1;
							}
							else
								Aimperfect_TR = 0;
						}
						else
							Atr = score_kmer(n+k*reps,k,ONE);

						if (Atr==Did || Aimperfect_TR) {
							reps++;
						}
						else {		/* ELSE FINAL NUMBER OF REPEATS (REPS) IS NOW KNOWN *****************/
							if (reps > 1)
								tela[n].cyc_l = k;		/* STORE # OF FRAMES CAN CYCLE THROUGH: AN ENTIRE UNIT-LENGTH */
							break;
						}
					}

					badslip_type = 0;
					/* SKIP CINCH IF TR PRIOR TO m SPANS INTO PRESENT TR */
					if (tela[n].X != n) {
						for (l=n-1; l>=tela[n].X; l--) {
							if ((q=tela[l].k) && (l + (p=span_rk(l)) - q) > m) { 
								if (p > reps*k) {
									badslip_type = 1;					/* FROM SEQUENCE IN TYPES: (1) -3-5-10-30-50-100-300-500 */
									options[1][50] += badslip_type;
									options[1][39] = 1; sprintf(dev_notes, "bslip sum %d", (int) options[1][50]);
									if (dev_print(MAIN,__LINE__)) {
									  	printf("badslip type %d at n=%d for k=%d with TR at l=%d.", badslip_type, n, k, l);
									}
									assign_tela(n++, row, a2D_n++, ONE);	/* MODES ZERO O-F-F, NON-ZERO ASSIGN  */
									Dtr = imperfect_TR = 0; /* SKIP PRESENT SLIP */
									break;
								}
							}    
						}
					}    
				}
				
				if (nuctransit && Dtr && tela[n].all_S != tela[n].all_r *k*match) {
					imperfect_TR++;
				}

				j = n + tela[n].all_k * tela[n].all_r;
				if (tela[n].all_L && check_tela(tela[n].all_L, j, ONE)!=3) {
					Dtr = imperfect_TR = 0;
				}
				else if (check_tela(m, j, ONE)!=3) {
					Dtr = imperfect_TR = 0;
				}

				/* IF SUMMING PATHBOX DIAGONAL 7/: COMMITTING TO CINCH AFTER AXIOMATIC TEST BY PUSH_TELA  ************************/
				if (Dtr) {
					if (imperfect_TR) {
						 o = push_tela(n,m, ONE);			/* WITH MODE ONE WILL ASSIGN TRANSITIONS WITHIN PUSH_TELA */
					}
					else {
						 o = push_tela(n,m, THREE);
					}
					if (o) {
						Dtr = imperfect_TR = 0;
						pull_tela(n);
						/* OFF flatline_after_TR(n);*/	/* WILL FLATLINE AT PREVIOUS TR */
						if (dev_print(MAIN,__LINE__)) {
							printf("push_tela violations=%d (+1 CONT, +2 EQUIV). Skipping k=%d-mer at n=%d.", o, k, n);
							print_tela(prtela_A, prtela_B);
						}
					}
				}
				/* IF SUMMING PATHBOX DIAGONAL 8/:  1st MEASUREMENT OF TANDEM REPEAT (TR) */ 
				if (Dtr==Did || imperfect_TR) {	
					tela[n].Dtr = Dtr;		/* SAVE Dtr SCORE */

					passR[2]++;

					for (i=0; i<k; i++)
						pathbox[n+i][m+i] = 114; 	/* "r" */

					r = 1;
					tela[n].r = reps;
					tela[n].k = k = n-m;
					slips[k]++;
					TRcheck = 1;

					while (TRcheck) {
						Atr = 0;

						if (r<reps) {
							if (nuctransit) { 
								for (i=m; i<n; i++) 
									Atr = Atr + pathbox[i][(i+(r+1)*k)];
								if (k>PISO && Atr!=Did && (100*Atr)/Did > DTHR) 
									imperfect_TR = 1;
							} 
						}
						else 
							Atr = TRcheck = 0;
						
						if (r<reps) {
							z=r*k;

							push_tela(n+z,m+z, THREE);

							for (i=0; i<k; i++) {
								pathbox[n+r*k+i][m+i] = 82; 	/* "R" */
							}

							r++;
							Atr = 0;
						}
						else {
							Atr = Did = Dtr = TRcheck = sumspan = conflict_flag = 0;

							if (imperfect_TR) {
								assign_transit(n,THREE); 	/* O-F-F; ONE=ALL_K/R; TWO=CYC_K/R; THREE=K/R */
							}

							/* IF CYCLE REPEAT, STORE CYCLE RUN. CYCLIC REPEATS CAN BE REPEATS IN MORE THAN ONE FRAME. MUST BE >2k */
							i = 0;			/* CYCLE[] ARRAY INDEX */
							cycle_flag = 0;	
							for (j = -1; j < r; j++) {				/* r = reps BECAUSE THIS IS IN ELSE EXIT LOOP */
								for (l = 0; l < k; l++) 
									cycle[i++] = tela[(n + j*k + l)].c; 	/* STORE WHOLE REPEAT */
							}
							for (l = 0; l < k-1; l++) {				/* STORE EXTENT OF PARTIAL REPEAT. CANNOT MATCH MORE THAN k */
								if (cycle[l] == tela[(o = n + r*k + l)].c) {
									cycle[i++] = tela[o].c;
								}
								else {
									if (reps == 1) 
										tela[n].cyc_l = l+1;		/* STORE LENGTH OF PARTIAL REPEAT */
									break;
								}
							}
							cycle[i] = '\0';
							tela[n].o = i;	/* STORE CYCLE LENGTH */

							if (!badslip_type) {
								/* NUMBER POSITIONS OF COLUMNS IN FRAME */
								for (l = 0; l < tela[n].cyc_l; l++) {
									tela[n+l].cyc_k = k;

									if (l == 0) {
										f = 1;	/* ROW NUMBER IN FRAMES ARRAY; OTHERWISE KEEP INCREMENTING */
										while (tela[m+l].cyc_F[f] && f < FRAME_ROWS)		/* FIND FIRST AVAILABLE ROW */
											f++;
										if (f==1) {
											tela[n].cyc_o = 'x';		/* NO CONFLICT SO WILL BE TAKING THIS FRAME */
										}
										else {
											conflict_flag = 1;
											tela[n].cyc_o = 'o';
											sumspan = tela[n].cyc_l;
											j=1;
											while (n-j >= 0 && tela[n-j].cyc_l == 0)
												++j;
											if (tela[n-j].cyc_l > tela[n].cyc_l) 
												sumspan = -tela[(z=n-j)].cyc_l;		/* POS. z IS WHERE TO START STORING PRODUCTS & SUMS OF PRODUCTS */
											else 
												z = n;								/* POS. z IS WHERE TO START STORING PRODUCTS & SUMS OF PRODUCTS */
										}
										tela[n].cyc_F[0] = f;	/* USE 0 ROW TO STORT LOCATION OF INDEXED UNIT TRs */	
									}
									else
										tela[n+l].cyc_F[0] = ++f;	/* USE ROW 0 TO STORE ROW # OF FRAME */

									for (j = 0; j < tela[n+l].all_r; j++) {
										if (j==0) {		/* WRITE FOR UNIT REPEAT STARTING AT m ONETIME */
											for (o = 0; o < k; o++) 
												tela[m + l + o].cyc_F[f] = o+1;
											if (!conflict_flag && l>0)
												tela[n+l].cyc_o = 'o';
										}
										for (o = 0; o < k; o++) {
											if (tela[n+j*k+l+o].c==tela[n-k+l+o].c) {
												tela[n+j*k+l+o].cyc_F[f] = o+1;
											}
											else if (imperfect_TR && tela[n+j*k+l+o].t==tela[n-k+l+o].t) {
												tela[n+j*k+l+o].cyc_F[f] = o+1;
											}
											else {	/* ELSE ERASE LAST PARTIAL UNIT */
												while (o >= 0) {
													tela[n + j*k + l + o].cyc_F[f] = 0;
													o--;
												}
												break;
											}
										}
									}
									tela[n+l].cyc_P = tela[n+l].cyc_k * tela[n+l].all_r;
									tela[n+l].cyc_r = tela[n+l].all_r;
								} 

								/* SUM UP COMPATIBLE TR PRODUCTS IN WINDOW OF LENGTH SUMSPAN BEGINNING AT POSITION z OF TR B */
								if (sumspan > 0 && tela[n].X != n) {	/* SUMSPAN IS LENGTH OF WINDOWS FOR WHICH SUMS OF PRODUCTS ARE RECORDED */
									for (j = 0; j < sumspan; j++) {
										for (f = tela[z+j].cyc_F[0] - 1 - j; f > 0; f--) {
											l=z+j-k;
											if (tela[l].cyc_F[f] == 1 && tela[l-1].cyc_F[f] != 0) {

												while (tela[l].cyc_o != 'x' && tela[l].cyc_o != 'o' && l>0) 
													l--;

												while (tela[l].cyc_F[f] != 1 && l>0) 
													l--;

												if (dev_print(MAIN,__LINE__)) {
													printf("n=%d, k=%d, l=%d, z=%d, j=%d of sumspan=%d.", n, k, l, z, j, sumspan);
												}

												tela[z+j].cyc_S = tela[z+j].cyc_P + tela[l].cyc_P;	
												tela[l  ].cyc_Rt = z+j;
												tela[z+j].cyc_Lf = l;
												break;
											}
											else if (tela[l].cyc_F[f] == 1 && tela[l-1].cyc_F[f] == 0) {
												while ((tela[l].cyc_F[f] != 1 || !tela[l].k) && l<n)
													l++;

												if (dev_print(MAIN,__LINE__)) {
													printf("n=%d, k=%d, l=%d, z=%d, j=%d of sumspan=%d.", n, k, l, z, j, sumspan);
												}

												tela[z+j].cyc_S = tela[z+j].cyc_P + tela[l].cyc_P;	
												tela[l  ].cyc_Rt = z+j;
												tela[z+j].cyc_Lf = l;
												break;
											}
											else if (tela[l].cyc_F[f] == 0) {
												tela[z+j].cyc_S = tela[z+j].cyc_P;
												break;
											}		
										}
									}
									/* FIND BEST CINCH SET */
									c = tela[(l=z)].cyc_S;	/* RUNNING BEST SCORE FOR CINCH SETS AT POSTION l */
									for (j = 1; j < sumspan; j++) {
										if (tela[z+j].cyc_S > c) {
											c = tela[ (l=z+j) ].cyc_S;
										}
									}
									if (l == z && tela[(j=tela[l].cyc_Lf)].cyc_o == 'o') {	
										i = j-1;	/* SAVE VAR j, CYCLING POSITION; VAR i TO COUNT DOWN TO POSITION THAT NEEDS TO BE CYCLED AWAY FROM */
										while (tela[i].cyc_o != 'x' && tela[i].cyc_o != blank) 
											i--;
										if (tela[i].cyc_o == 'x') {
											if (cyclelize_tela(i, j-i, n)) {
												badslip_type = 30;					/* FROM SEQUENCE IN TYPES: 1-3-5-10- (30) -50-100-300-500 */
												options[1][50] += badslip_type;
												options[1][39] = 1; sprintf(dev_notes, "bslip sum %d", (int) options[1][50]);
												if (dev_print(MAIN,__LINE__)) {
													printf("badslip type %d at n=%d for k=%d with TR at l=%d, delta=%d.", badslip_type, n, k, l, j-i);
												}

												a2D_n = tela[n].x+k; 
												row   = tela[n].y-1;
											}
										}
									}
								}
								else if (sumspan < 0 && tela[n].X != n) {
									for (j = 0; j < -sumspan; j++) {
										l = z+j + tela[z+j].cyc_k * (tela[z+j].cyc_r - 1);	/* VAR l IS POSITION 1 OF LAST UNIT OF REPEAT A */
										for (i = tela[n].cyc_F[0]; i < FRAME_ROWS; i++) {
											if (tela[l].cyc_F[i] == 1 && z+j != l+k) {
												tela[z+j].cyc_S = tela[z+j].cyc_P + tela[z+j+k].cyc_P;	
												tela[z+j].cyc_Rt = l+k;
												tela[l+k].cyc_Lf = z+j;
												tela[l+k].cyc_o = 'x';
												break;
											}
											else if (tela[l].cyc_F[i] == 0) { 				/* IF ZERO, THIS IS INCOMPATIBLE WITH ANY CYCLE OF B */
												tela[z+j].cyc_S = tela[z+j].cyc_P + 0;		/* PLUS ZERO IS FOR CODING CLARITY */
												break;
											}
											else if (i == FRAME_ROWS-1) {
												tela[z+j].cyc_S = tela[z+j].cyc_P + 0;		/* PLUS ZERO IS FOR CODING CLARITY */
												break;
											}
										}
									}
									/* FIND BEST CINCH SET */
									c = tela[z].cyc_S;				/* RUNNING BEST SCORE FOR CINCH SETS AT POSTION l */
									int cycto = z;					/* CYCLING POSITION, FOR CLARITY. WILL NOT STAY AT z */
									for (j = 1; j < -sumspan; j++) {
										if (tela[z+j].cyc_S > c) {
											c = tela[(cycto=z+j)].cyc_S;
										}
									}
									if (cycto != z) {
										if (cyclelize_tela(z, cycto-z, n)) {	/* REMINDER: cyclelize_tela(int cpos, int delta, int npos) */
											badslip_type = 50;					/* FROM SEQUENCE IN TYPES: 1-3-5-10-30- (50) -100-300-500 */
											options[1][50] += badslip_type;
											options[1][39] = 1; sprintf(dev_notes, "bslip sum %d", (int) options[1][50]);
											if (dev_print(MAIN,__LINE__)) {
												printf("badslip type %d at n=%d for k=%d with TR at cycto=%d, z=%d.", badslip_type, n, k, cycto, z);
											}

											a2D_n = tela[n].x+k; 
											row   = tela[n].y-1; 
										}
									}
								}
							}

							if (tela[n].o > 2*k) {
								cycle_flag = 1;		/* THIS BIT CAN BE USED TO ADD CYCLE NOTATION TO END OF LINE */
								if (dev_print(MAIN,__LINE__)) {
									printf("%d-mer cycle sequence of length %2d starting at %4d (n=%d): %s.", k, tela[n].o, n-k, n, cycle);
								}
							}
							else {
								for (l = 0; l <= WIDTH; l++)
									cycle[l] = '\0';
							}

							/* RECORD DNA "REVERB" IN SLIPLOC_ECHOES FOR ALL TR FRAMES */
							if (options[0][47] && tela[n].o) {    /********** OPTION TO SHOW SLIP LOCATIONS */
								if (tela[n].o > 2*k)
									h = k;
								else
									h = 1;

								for (j = 0; j < h; j++) {
									for (l = j; l+k <= tela[n].o; l+=k) {
										if 		(tela[n-k+l].echoes == blank)	tela[n-k+l].echoes = '(';
										else if (tela[n-k+l].echoes == '('  )	tela[n-k+l].echoes = '{';
										else if (tela[n-k+l].echoes == '{'  )	tela[n-k+l].echoes = '[';
										else if (tela[n-k+l].echoes == ')'||
												 tela[n-k+l].echoes == '}'||
												 tela[n-k+l].echoes == ']'  )	tela[n-k+l].echoes = 'X';
		
										if 		(tela[n-1+l].echoes == blank)	tela[n-1+l].echoes = ')';
										else if (tela[n-1+l].echoes == ')'  )	tela[n-1+l].echoes = '}';
										else if (tela[n-1+l].echoes == '}'  )	tela[n-1+l].echoes = ']';
										else if (tela[n-1+l].echoes == '('||
												 tela[n-1+l].echoes == '{'||
												 tela[n-1+l].echoes == '['  )	tela[n-1+l].echoes = 'X';
									}
								}
							}/**********************************************************************************/

						}
					} /* END OF else & TR_check = 0 PART */

					/* THE TR "SHADOW" IN COMMENTS REFERS TO REGION OF FIRST UNIT, SO-CALLED BECAUSE           	*/
					/* OF THE WAY CINCH-T BEGINS MARKING A TR STARTING AT THE SECOND UNIT. 						*/
					/* THE OVERALL STRATEGY IS TO LOOP FROM n-1 DOWN TO m+1 AND COUNT VARIOUS USEFUL THINGS.	*/					
					/* COUNT MOST RECENT CONFLICTING SLIP LENGTH IN UPSTREAM SHADOW OF NEW k-MER 				*/

					overslip = recslips = scooch = 0;

					if (!badslip_type) {
                        for (i = n-1; i > m; i--) {                            /* i WILL LOOP THROUGH TR SHADOW */
                            if (tela[i].r) {
								recslips = recslips + tela[i].r;
								overslip = overslip + span_rk(i);
							}
						}
					}
					
					if (tela[m].r) {
						recslips = recslips + tela[m].r;
						if (tela[m].r > 1) { 
							overslip = overslip + (tela[m].r-1) * tela[m].k;
                    	}  
					}

					for (i = 0; i < r; i++) {
						if ((ch=align2D[row][a2D_n]) != '\0') {
							scooch = overslip;
						}

						align2D[row  ][a2D_n+scooch  ] = '/'; 
						align2D[row  ][a2D_n+scooch+1] = '\0'; 	
						row++;										/* <==== ROW INCREMENTED HERE!		*/

						for (j = 0; j < (int) a2D_n; j++) {
							align2D[row][j] = blank;
						}
						for (j = 0; j < k; j++) {
							o = n + k*i + j;
							q = a2D_n - k + overslip + j;
							assign_tela(o, row, q, ONE);	/* MODES ZERO O-F-F, NON-ZERO ASSIGN  */
						}
					} /* END OF FOR i LOOP OF r */

					if (dev_print(MAIN,__LINE__)) {
						p = (int) a2D_n; q = lenseq;
						printf("\n               | via 2-D check_tela() = %d (+1 CONTINUITY, +2 EQUIVALENCE).\n", check_tela(0,p, TWO));
						printf("\n               | via 1-D check_tela() = %d (+1 CONTINUITY, +2 EQUIVALENCE).\n", check_tela(0,q, ONE));
						printf("\n               | print_tela for k=%d and r=%d at n=%d.", k, r, n);
						print_tela(0,56);
/*						print_tela(prtela_A, prtela_B);
*/
					}
					
					n = n + reps*k;
					a2D_n = a2D_n + overslip;
					r = 0;
					break;
				}
				else 				/* CASE WHEN Dtr != Did AND NOT IMPERFECT-TYPE TR */
					Dtr = Did = 0;
			}
		} /* END OF "for m" LOOP */
	} /* END OF "for n" LOOP => cinch_t PASS COMPLETE */

	/**********************************************/
	/* PRINT VALUES OF PATH BOX IF OPTION SET *****/

	if (options[0][25] == 1) {	/* opt_P */
		blocks = count_wrap_blocks(lenseq, options[1][58]);

		printf("\n\nPATHBOX CINCH PASS (length = width = %d)\n\n", lenseq);
		for (j = 0; j < blocks; j++) {
			if (blocks != 1)
				print_blockhead(j+1, blocks);
			line_end(PATHBOXHEAD, 9, 9);	
			for(n = j * options[1][58]; (n < (j+1) * options[1][58]) && (tela[n].c != '\0') && tela[n].c != '>'; n++) 
				printf("%2c", tela[n].c);
			printf("\n");
			for(m = j * options[1][58]; (m < (j+1) * options[1][58]) && (tela[m].c != '\0') && tela[m].c != '>'; m++) {
				printf("%4d. %c ", m+1, tela[m].c);
					for (n = j * options[1][58]; (n < (j+1) * options[1][58]) && (tela[n].c != '\0') && tela[n].c != '>'; n++) {
						if (m > n) {
							if (pathbox[m][n])
								printf("%2c", (char) pathbox[m][n]);
							else 
								printf("%2c", blank);
						}
						else if (n-m <= WIDTH)
							printf("%2d", pathbox[m][n]);
						else 
							printf("%2c", blank);
				}
				printf("\n");
	   		 }
		}
	} /* END OF OPTION TO PRINT PATHBOX */

 	align2D[row][citwidth+1] = '\0';
	options[1][2] = options[1][32] = citwidth;	/* ASSIGN CINCH-WIDTH TO HISTORY [0--9] AND LCL CURRENT [32]	*/
	clear_right(align2D);
	print1D();

	/********** 2. cinch_t MODULE: WRAPS LARGEST EXACT k-mers, IGNORES INTRA-TR TRs **********/
	i = ++options[1][18];	/* INCREMENT opt_I COUNTER (pass_num) */
	print_2Dseq(citwidth);
	passQ[i] = options[0][10];
/*	dev_prompt(MAIN,__LINE__,file_name);
*/
	if (recoverlen()==lenseq)
		update_tela();

	if (OFF && dev_print(MAIN,__LINE__)) {
		n = check_tela(0,lenseq, ONE);
		if (n!=3) {
			printf("\n  check_tela via 1-D coords, axioms = %2d (<3).", n);
		}
		n = check_tela(0,citwidth, TWO);
		if (n!=3) {
			printf("\n  check_tela via 2-D coords, axioms = %2d (<3).", n);
		}
		print_tela(prtela_A, prtela_B);
	}

	/********** 3. cinch_l MODULE: WRAPS HOMOPOLYMERIC RUNS IF >= 20 (2 * wrap VAR.) ********/
	i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */
	passR[i] = (int) cinch_l();
	if (cinchled) {
		print_2Dseq(options[1][32]);
	}
	passQ[i] = options[0][10];

	/********* 4. cinch_k MODULE: HANDLES k-mers FROM SIZE WIDTH DOWN TO k=1 ***********/
	i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */

	passR[i] = cinch_k();
	cycle_flag = print_2Dseq(options[1][32]);
	passQ[i] = options[0][10];
	dev_prompt(MAIN,__LINE__,file_name);

	/********* 5. nudgelize MODULE: "NUDGES" CONFLICT BY PUSHING COLS TO RIGHT ***************/
	continue_flag = 1;

	if (continue_flag) {
		i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */
		if (cycle_flag || align2D[0][0] == blank) {
			passR[5]++;

			while ((go_flag=nudgelize())!=0 && passR[5] < CYCMAX) {	
				passR[5]++;
			}
		}
		else {	
			options[1][i] = options[1][i-1];
		}
		passQ[i] = options[0][10];
	}

	/********* 6. cinch_d MODULE: HANDLES DE NOVO INTER-TR REPEATS *********************************/
	if (continue_flag) {
		i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */

		if (dev_print(MAIN,__LINE__)) {
			if (nuctransit) {
				printf("Pre-cinch_d report (p = perfect, i = imperfect tandem repeat):");
			}
			else {
				printf("Pre-cinch_d report:");
			}
		}
		intraTR_reps_tot = intraTR_reps = cinch_d(0);
	
		if (intraTR_reps_tot == 0) {
			printf(" Nothing left for cinch-d to cinch! \n");
			print_2Dseq(options[1][32]);
		}
	
		if (intraTR_reps_tot > 0) {
			while (intraTR_reps > 0) {
				intraTR_reps = cinch_d(1);
				passR[6]++;
			}
			options[0][6] = --passR[6];	/* STORE d runs in PASS SLOT IN CASE FUTURE MODULES POST CINCH-D NEED IT */
							/* NEED TO COUNT LIKE THIS FOR FOLLOWING REASONS: 			*/
							/*  1. LAST RUN IS A CHECK RETURNING ZERO, SO RUNS NEED TO BE DECREMENTED BY 1. 	*/
							/*  2. APPARENT LAST EFFECTIVE CINCH-D RUN MAY REVEAL A NEW CINCH-D OPPORTUNITY. 	*/
							/* 	OTHERWISE COULD HAVE SET WHILE LOOP TO > 1 (DON'T DO THIS. )			*/
		}
		passQ[i] = options[0][10];
	}

	/********* 7. relax_2D MODULE: DE-CINCHES HOMOPOLYMER RUNS IF THEY DID NOT AID CINCH-D *******/
	if (continue_flag && options[0][49]<1) {		/* opt_n DO NOT DO RELAX-2D */
		i = ++options[1][18];		/* INCREMENT opt_I ++PASS NUM */
	
		do {
			relax_length = options[1][32] ;
			relax_2D();
			relax_length = options[1][32] - relax_length;
			passR[8]++;
		}
		while (relax_length > 0);
	
		print_2Dseq(options[1][32]);
		passQ[i] = options[0][10];
	}	

	/********** 8. RECOVER_1D*************************************************/
	/* OPTION TO PRINT VALUES OF RECOVERED 1-D ALIGN BOX *********************/

	if (options[0][27]) { /* opt_R is 1 if ON */ 
		printf("\nChecking 1-D recovery from 2-D self-MHA:\n");
		recover_1D(recovered);

		r = (int) strlen(recovered) - 1;
		m = max(lenseq,r);
		blocks = count_wrap_blocks(m, options[1][58]);
		z = 0;			/* USE TO COUNT RECOVERED LETTERS IDENTICAL TO 1D */
		if (lenseq != r) {
			l = r - lenseq;
			for (i=0; i < m; i++) {
				if (tela[i].c != recovered[i]) {
					k = i;
					break;
				}
			}
			if (l>0) {
				for (i=m+l; i >= k; i--)
					tela[i].c = tela[i-l].c;	/* (WELL...) USE Seq[i] TO STORE GAPS FOR COMPARISON; NEVER tela[i].c */
				for (i=k; i < k+l; i++)
					tela[i].c = '-';
			}
			else {
				for (i=m-l; i >= k+l; i--)
					recovered[i-l] = recovered[i];
				for (i=k; i < k-l; i++)
					recovered[i] = '-';
			}
		}
		else
			l = 0;

		if (l) {
			warnhead('d');
			printf("1-D sequence differs in length from 2-D recovered 1-D sequence by %d letter(s).\n", l); 
		}

		q = 0;	/* ADJUSTS LINE NUMBERS IF SEQS DIFFERENT LENGTHS */
		for (j = 0; j < blocks; j++) {
			printf("\n");
			line_end(START, j+1, 0);	
	   		for (n = j * options[1][58]; n < (j+1) * options[1][58] && tela[n].c != '>'; n++) {
				printf("%1c", tela[n].c);
			}
			if (tela[n].c == '>') {
				printf("%1c", tela[n].c);  /* PRINTS TERMINAL CHARACTER '>' */
				if (options[0][21]) 
					line_end(END, lenseq, 0);
				else
					printf("  <== 1-D \n");
			}
			else {
				printf(" ");
				if (options[0][21]) {
					if (l>0 && n>k)
						q = l;	
					line_end(END, n-q, 0);
				}
				else 
					printf("\n");
			}
	
			line_end(SLIPS, j+1, 0);	
	   		for (n = j * options[1][58]; n < (j+1) * options[1][58] && n < m; n++) {
				if (seqtype == 1 || seqtype == 2) {			/* IF DNA OR RNA */
					if (tela[n].c == 'n' || recovered[n] == 'n') {
						printf("?");
						z++;	
					}
					else if (tela[n].c == recovered[n]) {
						printf("|");
						z++;
					}
					else {
						printf("*");
						recovery_flag++;
					}
				}
				else {						/* ELSE IF NOT NA */
					if (tela[n].c == recovered[n]) {
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
	
			line_end(START, j+1, 0);	
	   		for (n = j * options[1][58]; (n < (j+1) * options[1][58]) && recovered[n] != '>'; n++) {
				printf("%1c", recovered[n]);
			}
			if (recovered[n] == '>') {
				printf("%1c", recovered[n]);  /* PRINTS TERMINAL CHARACTER '>' */
				if (options[0][21]) {
					line_end(END, r, 0);
				}
				else
					printf("  <== 2-D\n");
			}
			else {
				printf(" ");
				if (options[0][21]) 
					line_end(END, n, 0);
				else 
					printf("\n");
			}
		} /* END OF FOR j LOOP */
		printf("\n");

		options[1][8] = z;			/* STORE NUMBER OF RECOVERED CHARACTERS */

		if (recovery_flag) {		/* LAST ROW OF array2D WILL STORE CONSENSUS, SO NEED TO KEEP CLEAR */
			warnhead('R');
			printf("Imperfect recovery of 1-D sequence from 2-D self-MHA.\n");
			passQ[8] = 1000*(lenseq-recovery_flag)/lenseq;
		}
		else {
			line_end(SLIPS, 0, 0);
			printf("Perfect recovery of 1-D sequence.\n");
			passQ[8] = 1000; 
		}
	} /* END of opt_R */

	/* PRINT OPTION FOR K-MER REPORT AFTER cinch_t **********************/
	if (options[0][46]) {	/* opt_k */

		int k_start=2;

		printf("\n Unique slips: %4d\n", passR[2]);

		for (i = k_start; i <= 10; i++)
			printf(" %smers:%3d \n", nmer_prefix(i), slips[i]);
		for (i = 11; i <= WIDTH; i++)
			printf(" %d-mers:%2d\n", i, slips[i]);
	}
	/***************************************************************************/
	c   = options[1][1];		/* REUSING c VAR FOR FORMATTED STRING LENGTH */
	if (options[1][33] > 1)		/* opt_XX FISHER-YATES RANDOMIZATIAN */
		row = FY_size;			/* USE IN PLACE OF ORIGINAL STRING LENGTH */
	else
		row = options[1][0];	/* RESUING row VAR FOR ORIGINAL STRING LENGTH */
	/***************************************************************************/

	if (seqtype == 1)		
		printf(   "\n   PASS      Width cinch history for %s (DNA)\n", Seq_head);
	else if (seqtype==2)	
		printf(   "\n   PASS      Width cinch history for %s (RNA)\n", Seq_head);
	else if (seqtype==3)	
		printf(   "\n   PASS      Width cinch history for %s (PROTEIN)\n", Seq_head);
	else if (seqtype==0)
		printf(   "\n   PASS      Width cinch history for %s (BABYLONIAN)\n", Seq_head);
	else
		printf(   "\n   PASS      Width cinch history:\n");

	for (i = 0; options[1][i] != '\0' && i < 8; i++) {
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
			if (passR[2]>1)
				printf("%s post cinch-t   [pass #2: %d cinches]\n", letr_unit, passR[2]);	/* STYLE: USE DASHED NAME FOR PRINTING, UNDERSCORED FOR CODING 	*/
			else if (passR[2])
				printf("%s post cinch-t   [pass #2: %d cinch]\n", letr_unit, passR[2]);	/* STYLE: USE DASHED NAME FOR PRINTING, UNDERSCORED FOR CODING 	*/
			else if (options[0][55]) 
				printf("%s post cinch-t   [pass #2: SKIPPED BY REQUEST]\n", letr_unit);	/* STYLE: USE DASHED NAME FOR PRINTING, UNDERSCORED FOR CODING 	*/
			else 
				printf("No effective cinch-t cinches taken.\n");
			break;												/*  "cinch-x" VERSUS 'cinch_x' programming calls 				*/ 
		case 3:													/*  USEFUL FOR SEARCHING CODE.									*/
			printf("%s post cinch-l   [pass #3: %d run(s)]\n", letr_unit, passR[3]);
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
				printf("%s post nudgelize [pass #5]\n", letr_unit);
			else if (passR[5] == 1) {	
				if (k == 3)
					printf("%s post nudgelize [pass #5: %d run; nudge-cyclelized TR]\n", letr_unit, passR[5]);
				else if (k == 4)
					printf("%s post nudgelize [pass #5: %d run; tip-cyclelized TR]\n", letr_unit, passR[5]);
			}
			else if (passR[5] <= CYCMAX) {	/* IN WHICH CASE k WILL BE NON-ZERO */
				if (k == 3)
					printf("%s post nudgelize [pass #5: %d runs; last TR was nudge-cyclelized]\n", letr_unit, passR[5]);
				else if (k == 4)
					printf("%s post nudgelize [pass #5: %d runs; last TR was tip-cyclelized]\n", letr_unit, passR[5]);
			}
			else 							/* IN WHICH CASE cyc runs >> CYCMAX */
				printf("%s post nudgelize [pass #5: reverted after %d runs due to gnarly micro-foam]\n", letr_unit, passR[5]-CYCMAX*1000);
			break;
		case 6:	
			if (passR[6] > 0)
				printf("%s post cinch-d   [pass #6: %d runs]\n", letr_unit, passR[6]);
			else if (options[0][39])
				printf("%s post cinch-d   [pass #6: SKIPPED BY REQUEST]\n", letr_unit);
			else
				printf("%s post cinch-d   [pass #6]\n", letr_unit);
			break;
		case 7:	
			printf(    "%s post relax-2D  [pass #7: relaxed %d runs]\n", letr_unit, passR[7]);
			break;
		}
	}
	if (options[0][27]) { /* opt_R */ 
		printf("  %5d       => %4ld ", passQ[8], options[1][8]);
			printf(    "%s recovered 1D   [final check pass]\n", letr_unit);
	}

	if (continue_flag) {
		printf("\n Width cinch ratio (post cinch-d):  %2.3f", ratio1=(float)options[1][6]/lenseq);
		if (options[0][49] < 1)
			printf("\n Width cinch ratio (post relax-2D): %.3f\n\n", ratio2=(float)options[1][32]/lenseq);
		else
			printf("\n\n");
	}
	else
		printf("\n Width cinch ratio (post cinch-k): %.3f\n\n", ratio1=ratio2=(float)options[1][4]/lenseq);

	if (options[0][24]) {							/* OPTION TO OUTPUT 2-D ALIGNMENT & CONSENSUS STRING TO FILE */
		fp_cons = fopen("Surf_barrels.log", "a");
		fprintf(fp_cons,">%s (%d > %d %s) x%d\n", 
			file_name, (int) options[1][1], (int) options[1][6], letr_unit, (int) options[1][59]);	

		for (m = 0; align2D[m][0] != '\0'; m++) {
			fprintf(fp_cons, " %s\n", align2D[m]);
		}

		fprintf(fp_cons, " %s\n\n", consensus);
		fclose(fp_cons);
	}
	if (msa && options[0][24]==2) {
		tuck = get2Dtucknum(m2Dalig, align2D);
		if (tuck<0) {
			tuck = m2Da_height;	/* GET2DTUCKNUM IS BUST AND MUST BE DEBUGGED, USE ZERO TUCK */
			warnhead('t');
			printf("There was a problem computing 2Dtucknum, so it was set to height = %d.\n\n", tuck);
		}
		/* SUPERIMPOSE SECOND 2D ARRAY OVER FIRST FOR COMPARISON */
		m2Dalig[tuck][0] = '>';
		for (i=0; align2D[i][0]!='\0'; i++) {
			if (i>0)
				m2Dalig[tuck+i][0]= ' ';
			for (j=1; (ch=align2D[i][j-1])!='/' && ch!=optR && ch!='>'; j++) {
				m2Dalig[tuck+i][j] = ch;
			}
				m2Dalig[tuck+i][j] = ch;
		}
		fp_msa = fopen("TUBES.mha", "w");
			fprintf(fp_msa, "%s\n", m2Dalig[0]);
		for (m = 1; m2Dalig[m][1] != '\0'; m++) {
			fprintf(fp_msa, "%s\n", m2Dalig[m]);
		}

		fclose(fp_msa);
	}
	else if (options[0][24]==2) {
		fp_msa = fopen("TUBES.mha", "a");
			fprintf(fp_msa, ">%s\n", align2D[0]);
		for (m = 1; align2D[m][0] != '\0'; m++) {
			fprintf(fp_msa, " %s\n", align2D[m]);
		}

		fclose(fp_msa);
	}
	
	if (options[0][54] != 1) {			/* ONLY IF opt_s OPTION TO SILENCE OUTPUT IS NOT ON */
		fp_out = fopen("Surf_wavereport.mha", "a");
		fprintf(fp_out, "v%s\t%.20s\t x%ld\t%4ld\t%.3f\tCYC:%3d (t=%ld)\tRND:%.*s\t%38s -%ld (%4ld %s) REC:%4d\t%3ld\t%s\n", 
				version, time0+4, options[1][59], options[0][10], ratio1, passR[5], options[0][5], (int) options[1][33], "XX", 
				file_name, options[0][29], options[1][1], letr_unit, passQ[8], options[1][6], dev_notes);
		fclose(fp_out);

		/* IF IMPERFECT CONSENSUS OR IF CYCLELIZE REVERTED */
		if (options[0][10] != 1000 || passR[5] > CYCMAX) {
			fp_tricksy = fopen("waves/foam_and_chowder.mha", "a");
			fprintf(fp_tricksy, "v%s\t%.20s\t x%ld\t%4ld\t%.3f\tCYC:%2d (t=%ld)\tRND:-%.*s\t%s -%ld (%ld %s) REC:%4d\t%s\n", 
					version, time0+4, options[1][59], options[0][10], ratio1, passR[5], options[0][5], (int) options[1][33], "XX", 
					file_name, options[0][29], options[1][1], letr_unit, passQ[8], dev_notes);
			for(n = 0; n<lenseq; n++) {
				fprintf(fp_tricksy, "%c", tela[n].c);
			}
			fprintf(fp_tricksy, "\n");
			fclose(fp_tricksy);
		}
	}
	/* dev_prompt(MAIN,__LINE__,file_name); */

	exit(EXIT_GOOD);	/* Exit main(). */
} 


/*************************************************************************************************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal */
/*************************************************************************************************************/
