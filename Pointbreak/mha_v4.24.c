/*******************************************************************************************************/
/***** The program "maximal" is a micro homology alignment (MHA) program.                          *****/
/***** Designed and written by Dr. Albert J. Erives. 2017-2020. AGPL-3.0.                          *****/
/***** Code repository located at https://github.com/microfoam/maximal. Licensed under AGPL-3.0.   *****/
/***** This program renders a 1-D DNA sequence into a 2-D self-alignment to rescue micro-paralogy. *****/
/*******************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <ctype.h>  		/* isalpha() 		*/
#include <stdlib.h> 		/* system()  		*/
#include <string.h> 		/* strlen()  		*/
#include <time.h>			/* difftime() 		*/
#include <signal.h>			/* signal() 		*/
#include <unistd.h>			/* signal()-related */
#include "microhomology.h"	/* maximal header: main header file                                  */
#include "microhom-devl.h"	/* maximal header: program development code, testing, and evaluation */
#include "microhom-skor.h"	/* maximal header: alignment scoring definitions and functions       */
#include "microhom-tela.h"	/* maximal header: verb_tela() functions and coord struct-related    */


int main(int argc, char *argv[])
{
	int match      = MATCH;
	int transition = TRANSITION;			
	int mismatch   = MISMATCH;
	char version[] = "4.24";				/* current version number */

	int c=0, f=0, i=0; 
	short unsigned int continue_flag=0;

	int j=0, h=0, k=0; 
	short unsigned int imperfect_TR=0; 

	int l=0, m=0, n=0; 
	short unsigned int Aimperfect_TR=0; 

	int o=0, p=0, q=0;						/* DEV: RESERVE opq FOR chk_tela */
	short unsigned int nuctransit=0; 

	int reps=0, r=0, z=0;	
	short unsigned int seqtype=0; 

	int DTHR = 100;							/* 100%. RESET IF opt_x, k-DEPENDENT DIAGONAL THRESHOLDS < SCORES */ 
	int m2Da_height = 1;	
	int tuck;
	short unsigned int conflict_flag=0;

	int number=0;
	int sumspan=0;
	int homopoly_flag=0;
	short unsigned int msa = 0;				/* BIT FLAG FOR MSA CO-INUPUT */

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

 	char m2Dalig[MAXROW+1][MAXROW] = {{0}};			
	char letr_unit[8] = {0};				/* UNIT STRING: "bp" FOR DNA, "nt" FOR RNA, 'aa' FOR PROTEINS, 'ch' FOR ALL OTHER */
	char cycle[WIDTH+1];		/* THIS ARRAY HOLDS THE CYCLIC PATTERN OF TRs W/ >2 UNITS */
	char numstring[8] = {0};
	char Seq_head[100] = {0};	/* FASTA HEADER */
	char Seq_i[MAXROW] = {0}; 	/* INPUT SEQUENCE */
	char Seq_r[MAXROW] = {0}; 	/* RANDOMIZED SEQUENCE */
	char *Seq = Seq_i;			/* POINTER TO INPUT SEQUENCE */
	int passQ[16] = {0};        			/* PASS QUALITY */
	int passR[16] = {0};        			/* PASS RUNS */

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
			case 'c':						/* OPTION TO USE REVERSE COMPLEMENT */
					options[0][38] = 1;		/* opt_c ON		*/
					break;
			case 'f':						/* OPTION TO SHOW FOAM-FREE SEGMENTS BELOW CONSENSUS ROW*/
					options[0][41] = 1;		/* opt_f ON		*/
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
			case 'C':						/* SHOW BASE 62 CODE */
					print_base62_table();
					exit(EXIT_EARLY);
					break;
			case 'D':
				   	options[0][13] = 1;		/* opt_D: SHOW DTHR_lookup values (later below). Will exit early. */
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
					options[0][49] = 1;		/* opt_n NO RELAX 2-D */
					break;
			case 'P':						/* OPTION TO PRINT PATH BOX */
					options[0][25] = 1;		/* opt_P ON 	*/
					break;
			case 'R':						/* OPTION TO PRINT RECOVERED 1-D SEQUENCE FROM LAST 2-D */
					++options[0][27];		/* opt_R = 1 is ON; opt_R > 1 is O-F-F (lets you toggle off after option v) */
					break;
			case 'X':						/* OPTION TO SCRAMBLE SEQUENCE           */
					options[0][33] = 1;		/* opt_X ON 							 */
					++options[1][33];		/* INCREMENT opt_X to desired level 	 */
					break;					/*  X = rand() cheese, XX = FISHER-YATES */
			case 'Y':						/* OPTION TO SPECIFY FY_size	*/
					options[0][34] = 1;
					if (number) {
						if (number < MAXROW)
							FY_size = number;
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
		if (options[1][33])
			printf(" -X [USE RANDOMIZED SEQUENCE]");
		else 
			printf(" -XX [USE FISHER-YATES RANDOMIZED SEQUENCE]");
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

	if (options[0][13]) {	/* opt_D: SHOW DTHR VALUES */
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

	if (nuctransit && options[0][38]) {		/* opt_c USE REVERSE COMPLEMENT */
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
				assign_tela(n++, row, a2D_n++, TWO);	/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */
			}
		}
		else {						/* ELSE opt_t: SKIP CINCH-T */ 
			strcpy(align2D[0],Seq);
			break;
		}

		/* FOR COLUMN n LOOP 2/3: SKIP PRESENT TR IF CONFLICT AND CAN CYCLE WITH SAME SCORE */
		if (tela[n].all_L && tela[n].all_S == tela[n+1].all_S && !tela[n+1].all_L && tela[(tela[n].all_L)].cyc_o == 'x') 
			assign_tela(n++, row, a2D_n++, TWO);	/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */

		/* FOR COLUMN n LOOP 3/3 */
		for (m = 0; m < n; m++) {
			/* FOR ROW m LOOP 1/6: UPDATE VAR CITWIDTH AT END */
			if (tela[n].c == '>') {
				citwidth = a2D_n;
				align2D[row+1][0] = '\0';
				options[1][17] = row+1;
			}
			else if (tela[n].cyc_o == 'o') {	/* IF THIS POSITION HAD A BIGGER k-MER SQUASHED HERE (E.G., FOR LATER CYCLING) */
				assign_tela(n++, row, a2D_n++, TWO);	/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */
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
				assign_tela(n++, row, a2D_n++, TWO);	/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */
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
					assign_tela(n++, row, a2D_n++, TWO);	/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */
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
										assign_tela(n++, row, a2D_n++, TWO);	/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */
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
									assign_tela(n++, row, a2D_n++, TWO);	/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */
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
						assign_tela(n++, row, a2D_n++, ONE);	/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */
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
											for (p=n; p<lenseq; p++) 
												tela[p].y = row;
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
							assign_tela(o, row, q, TWO);		/* MODES O-F-F, ONE--FIVE: ONE=FLATLINE AT N, TWO=ASSIGN  */
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
	dev_prompt(MAIN,__LINE__,file_name);

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
	if (passR[i])
		print_2Dseq(options[1][32]);
	passQ[i] = options[0][10];

	/********* 4. cinch_k MODULE: HANDLES k-mers FROM SIZE WIDTH DOWN TO k=1 ***********/
	i = ++options[1][18];	/* INCREMENT opt_I ++PASS NUM */

	passR[i] = cinch_k();
	cycle_flag = print_2Dseq(options[1][32]);
	passQ[i] = options[0][10];

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
					printf("%s post nudgelize [pass #5: %d run; nudge-cyclelized TR of type k>2]\n", letr_unit, passR[5]);
				else if (k == 4)
					printf("%s post nudgelize [pass #5: %d run; tip-cyclelized TR]\n", letr_unit, passR[5]);
			}
			else if (passR[5] <= CYCMAX) {	/* IN WHICH CASE k WILL BE NON-ZERO */
				if (k == 3)
					printf("%s post nudgelize [pass #5: %d runs; last TR was nudge-cyclelized (k>2)]\n", letr_unit, passR[5]);
				else if (k == 4)
					printf("%s post nudgelize [pass #5: %d runs; last TR was tip-cyclelized]\n", letr_unit, passR[5]);
			}
			else 							/* IN WHICH CASE cyc runs >> CYCMAX */
				printf("%s post nudgelize [pass #5: reverted after %d runs due to gnarly micro-foam]\n", letr_unit, passR[5]-CYCMAX*1000);
			break;
		case 6:	
			if (passR[6] > 0)
				printf("%s post cinch-d   [pass #6: %d runs]\n", letr_unit, passR[6]);
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
		fp_cons = fopen("Surf_barrels.log", "a");	/* FOPEN RIGHT BEFORE WRITING TO MINIMIZE CHANCE OF CLOSING W/ OPEN FILES */
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
		fp_out = fopen("Surf_wavereport.mha", "a");		/* FOPEN RIGHT BEFORE WRITING TO MINIMIZE CHANCE OF CLOSING WITH OPEN FILES */
		fprintf(fp_out, "v%s\t%.20s\t x%ld\t%4ld\t%.3f\tCYC:%3d (k=%ld)\tRND:%.*s\t%38s -%ld (%4ld %s) REC:%4d\t%3ld\t%s\n", 
				version, time0+4, options[1][59], options[0][10], ratio1, passR[5], options[0][5], (int) options[1][33], "XX", 
				file_name, options[0][38], options[1][1], letr_unit, passQ[8], options[1][6], dev_notes);
		fclose(fp_out);

		/* IF IMPERFECT CONSENSUS OR IF CYCLELIZE REVERTED */
		if (options[0][10] != 1000 || passR[5] > CYCMAX) {
			fp_tricksy = fopen("waves/foam_and_chowder.mha", "a");
			fprintf(fp_tricksy, "v%s\t%.20s\t x%ld\t%4ld\t%.3f\tCYC:%2d (k=%ld)\tRND:-%.*s\t%s -%ld (%ld %s) REC:%4d\t%s\n", 
					version, time0+4, options[1][59], options[0][10], ratio1, passR[5], options[0][5], (int) options[1][33], "XX", 
					file_name, options[0][38], options[1][1], letr_unit, passQ[8], dev_notes);
			for(n = 0; tela[n].c != '\0'; n++) {
				if (tela[n].c != 10 && tela[n].c != 13 && tela[n].c != EOF)
					fprintf(fp_tricksy, "%c", tela[n].c);
			}
			fprintf(fp_tricksy, "\n");
			fclose(fp_tricksy);
		}
	}
/*	dev_prompt(MAIN,__LINE__,file_name);
*/	exit(EXIT_GOOD);
} 
/******************************************* END OF MAIN() **********************************************************************/
/********************************************************************************************************************************/

/*****************************************/
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
/*		warnhead('z');
		printf("get_1Dz count_check = %d\n", count_check);
*/
		return(-z);		/* KEEP NEGATIVE TO PRESERVE INFO/ AVOID RETURNED VALUE LETTING A FOR LOOP GO */
	}
}

/*** FUNCTION 01 ************************************************************************************/
int cinch_l(void) 
{
int cil_row=0, i=0, j=0, k=0, l=0, m=0, n=0, run=0, x=0;
char letr;
char lopt_Q_left = (char) options[1][26];		/* LHS character delimiter for homopolymer run */
char lopt_R_rght = (char) options[1][27];		/* RHS character delimiter for homopolymer run */
char blnk        = (char) options[1][11];		/* opt_B blank character		*/
int cil_mwrap    = (char) options[1][22];		/* opt_M long_homopolymer_run	*/
int cinchled=0;					/* BIT FLAG TO SAY cinch_l DID SOMETHING */
char lclalign2D[MAXROW][MAXROW] = {{0}};

	for (m = 0; align2D[m][0] != '\0'; m++) {
		for (n = 0; align2D[m][n] != '\0'; n++) {
			while (align2D[m][n] == blnk) {	/* MOVE WINDOW PAST BLANKS */
				lclalign2D[m+cil_row][n-x] = blnk;
				n++;			  
			}
			/* letr ASSIGNED */
			if ((letr=align2D[m][n]) != '/' || letr != '\0') {	/* letr != '\0' USEFUL IF BLANKS FILLED TO SCRIMMAGE */ 
				lclalign2D[m+cil_row][n-x] = letr;	
				if (letr == align2D[m][n-1]) 
					run++;
				else
					run = 1;
			}
			else
				lclalign2D[m+cil_row][n-x] = letr;		/* WRITE TERMINAL MHA CHARACTERS '/' or '>' */
	
			if (run == 2*cil_mwrap) {			/* TRIGGER LENGTH MEASURE OF MONOMERIC RUN	*/
												/* AND WRITING OF (10)-MER BLOCKS			*/
				++cinchled;
				while (align2D[m][n+run-2*cil_mwrap] == letr) {
					run++;
				}
				
				lclalign2D[    m+cil_row][n-x-cil_mwrap+1  ] = lopt_R_rght;	/* WRITE SLIP AFTER FIRST 10 */
				for (l = 1; l < cil_mwrap; l++)								
					lclalign2D[m+cil_row][n-x-cil_mwrap+1+l] = '\0';
				cil_row++;														/* ADVANCE TO NEXT ROW, CONT.*/
	
				for (j = (run-cil_mwrap)/cil_mwrap; j > 0; j--) {
					for (i = 0; i < n-x-2*cil_mwrap; i++)
						lclalign2D[m+cil_row][i] = blnk;
					lclalign2D[m+cil_row][i] = lopt_Q_left;					/* MARK LEFT EDGE OF MWRAP RUN */
					for (k = 0; k < cil_mwrap; k++) {
						lclalign2D[m+cil_row][n-x-2*cil_mwrap+1+k] = letr; /* FILL W/ MONOMER LETTER  */
					}
					if (j > 1) {
						lclalign2D[m + cil_row  ][n-x-cil_mwrap+1] = lopt_R_rght;
                        lclalign2D[m + cil_row++][n-x-cil_mwrap+2] = '\0';             /* SHOULD NOT BE NECESSARY */
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
	options[1][3] = options[1][32];

	if (cinchled) {
		mha_writeback(lclalign2D, align2D); 
	}
	return(cinchled);
}

/*** FUNCTION 02 ************************************************************************************/
int cinch_k(void) 
{
int cik_row=0, i=0, k=0, l=0, m=0, n=0, scrimmage_line = -1, x=0, y=0, r=0, q=0; 
int first_mwrap_start=0, last_mwrap=0;
unsigned short int first_mwrap=0, keep_checking=1;
unsigned short int nuctype = options[1][13];		/* EQUALS ONE IF DNA STRING, TWO IF RNA, THREE IF PROTEIN */
unsigned short int nuctransit=0, check_imperf=0;	/* BIT FLAG FOR HANDLING NUCLEOTIDE TRANSITIONS SILENTLY (IGNORING) */
unsigned short int homopolyflag=0, imperfect_TR=0;
int sum4score;		/* SCORE VAR FOR IMPERFECT TR'S */
char letr, letr2, letr3;
char blnk        = (char) options[1][11];		/* opt_B blank character */
int  cik_mwrap   =        options[1][22];		/* opt_M long_homopolymer_run */
char kopt_Q_left = (char) options[1][26];		/* LHS character delimiter for homopolymer Run */
char kopt_R_rght = (char) options[1][27];		/* RHS character delimiter for homopolymer Run */
int max_k = (int) options[1][46]/2;				/* MAX k-SIZE FROM mark_tela() */
int lenseq = options[1][1];
int *x_history = NULL;
int symbol_count = 0;
char cik_align2D[MAXROW][MAXROW] = {{0}};
/* char **cik_align2D = NULL; */

	x_history = (int *)calloc(lenseq, sizeof(int));
/*	cik_align2D = (char **)calloc(lenseq, sizeof(int));
	for (i=0; i<lenseq; i++) {
		cik_align2D[i] = (char *)calloc(lenseq, sizeof(int));
	}
*/
	if (nuctype == 1) {		/* IF DNA */
		nuctransit = 1;
	}	
	if (dev_print(MAIN,__LINE__)) {
		printf("Post-cinch-t max_k = %d.", max_k);
	}
	if (!max_k)
		max_k = WIDTH;

	/* START AT BIGGEST k-MER POSSIBLE AT 2x */
	for (k = max_k; k > 0; k--) {

		cik_row = x = symbol_count = 0;	
		scrimmage_line = -1;

		for (m = 0; align2D[m][0] != '\0'; m++) {

			y = 0;

			/* CHECK LINE AHEAD OF TIME FOR FIRST MONO-RUN TERMINATOR */
			for (i = 0; (letr=align2D[m][i]) != '\0'; i++) {	
				/* NOTE: APPARENT REDUNDANCY BELOW ALLOWS CHANGING L/R RUN DELIMITERS TO BE EQUAL */
				if (letr == kopt_R_rght && align2D[m][i-1] != blnk) {
					first_mwrap = 1;					/* TURN ON BIT FLAG FOR HANDLING FIRST MONO-RUN */
					last_mwrap = 1;						/* TURN ON BIT FLAG FOR HANDLING TO LAST MONO-RUN */
					first_mwrap_start = i-cik_mwrap;	/* SAVE POSITION OF START OF MWRAP */
					break;								/* BREAK OUT OF FOR i CHECK LOOP */
				}
			}
		
			for (n = 0; align2D[m][n] != '\0'; n++) { 
				keep_checking = 1;			/* THIS FLAG HANDLES THE CONTINUED NEED TO CHECK FOR INTRA-TR REPEATS   */
				imperfect_TR = 0;			/* THIS FLAG IS TURNED ON (SET TO ONE) WHEN TR W/ TRANSITION MISMATCHES IS FOUND */

				if (n == 0 && isalpha(align2D[m][0])) {
					x = 0;
					x_history[0] = x;
				}

				/* MOVE WINDOW PAST INITIAL BLANKS */
				while (align2D[m][n] == blnk) {
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
					symbol_count += cik_mwrap;
					for (n = first_mwrap_start; (letr=align2D[m][n]) != '\0'; n++) {
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
					if (align2D[m][n] == kopt_Q_left && align2D[m][n+cik_mwrap+1] == kopt_R_rght) {		
						symbol_count += cik_mwrap;
						cik_align2D[m+cik_row][n-x] = kopt_Q_left;
						n++;
						while ( (letr=align2D[m][n]) != kopt_R_rght) {
							cik_align2D[m+cik_row][n-x] = letr;
							n++;
						}
						cik_align2D[m+cik_row][n-x] = letr;			/* WILL WRITE kopt_R_rght */
						for (i = cik_mwrap; i > 0; i--)				/* FAILSAFE:			*/
							cik_align2D[m+cik_row][n-x+i] = '\0';	/* OVERWRITE W/ NULLS	*/
						break;										/* BREAK OUT OF n LOOP	*/
					}
					else if (align2D[m][n] == kopt_Q_left && align2D[m][n+cik_mwrap+1] != kopt_R_rght) {
						symbol_count += cik_mwrap;
						last_mwrap = 0;
						for (i = 0; i <= cik_mwrap; i++) {
							cik_align2D[m+cik_row][n-x+i] = align2D[m][n+i];
						}
						letr = align2D[m][n+1];
						n = n + cik_mwrap;
						while (align2D[m][n] == letr) {
							symbol_count++;
							cik_align2D[m+cik_row][n-x] = letr;
							x_history[n] = x;
							n++;
						}
						symbol_count--;		/* EACH BLOCK OF MWRAPS OVERCOUNTS BY ONE AFTER THE LAST WHILE LOOP */
					}
				}

				/* CHECK FOR & DEAL WITH LINE ENDS TOO SHORT TO HARBOR TR OF SIZE k */
				if (!isalpha(align2D[m][n+2*k-1])) {	/* TRUE IF WINDOW < 2x k-MER, WRITE REST OF LINE TO cik_align2D */
					for (i = n; align2D[m][i] != '\0'; i++) {
						cik_align2D[m+cik_row][i-x] = align2D[m][i];
						if (isalpha(align2D[m][i])) {
							symbol_count++;
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
						if ( ((letr =align2D[m][n+l  ])=='A' || letr =='G') &&
							 ((letr2=align2D[m][n+l+k])=='C' || letr2=='T') ) {
							keep_checking = 0;
							break;
						}
						else if ( ((letr =align2D[m][n+l  ])=='C' || letr =='T') &&
							      ((letr2=align2D[m][n+l+k])=='A' || letr2=='G') ) {
							keep_checking = 0;
							break;
						}
						/* FOR AESTHETICS AND MORE? DON'T SCOOCH MONOS INTO TRANSITS */
						else if (k==1 && ((consensus[n-x+y  ]=='R' || consensus[n-x+y  ]=='Y') ||
										  (consensus[n-x+y+1]=='R' || consensus[n-x+y+1]=='Y'))) {
							keep_checking = 0;
							break;
						}
						else {
							/* CHECK TO SEE IF THERE ARE n's */
							if (nuctype && (align2D[m][n+l]=='n' || align2D[m][n+l+k]=='n')) {
								keep_checking = 0;
								break;
							}
	
							/* CHECK TO SEE IF HOMOPOLYMER RUN CAN BE EXCLUDED */
							if (homopolyflag && l > 0) {	/* IF l > 0, THEN homopolyflag IS 1 */
								if ((align2D[m][n + l] != align2D[m][n+k+l  ]) ||
								    (align2D[m][n + l] != align2D[m][n + l-1]) ||
								    (align2D[m][n+k+l] != align2D[m][n+k+l-1]) ) {
									homopolyflag = 0;
								}
							}
	
							if (nuctransit && keep_checking) {
								if (n == scrimmage_line || col_isclear(align2D,n,m,1) < 0) {
									y = 0;		/* RESET y VAR. B/C NO LONGER NEED TO ADJUST CONSENSUS COORDINATES */
								}
								if ( (letr2=consensus[n-x+y+l]) == 'R' || letr2 == 'Y') {
									keep_checking = 0;
								}
								else if (k>0 && ((letr3=consensus[(q=n-x+y+l+k)]) == 'R' || letr3 == 'Y')) {
									if (k==1 && col_isclear(align2D,q,m,1)>0)
										;
									else
										keep_checking = 0;
								}
							}
							if (keep_checking && align2D[m][n+l] != align2D[m][n+k+l]) {
								keep_checking = 0;
								if (nuctransit && k >= PISO) {
									check_imperf = 1;
								}
							}
						}
					} /* END OF FOR l SCAN LOOPS */

					if (homopolyflag && keep_checking) {		/* IF HOMOPOLYFLAG=1 WAS NOT SET TO ZERO */
						keep_checking = 0;
					}
					if (homopolyflag && check_imperf) {			/* TURN THIS O-F-F AS WELL */
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
						letr =align2D[m][n+l  ];
						letr2=align2D[m][n+l+k];
						if (islower(letr)) {
							letr = toupper(letr);
						}
						if (islower(letr2)) {
							letr2= toupper(letr2);
						}
						if ( (letr=='A' || letr=='G') && (letr2=='C' || letr2=='T') ) {
							break;
						}
						else if ( (letr=='C' || letr=='T') && (letr2=='A' || letr2=='G') ) {
							break;
						}

						if (align2D[m][n+l] != align2D[m][n+k+l]) { 
							if ( (letr =consensus[n-x+y + l]) != 'R' && letr  != 'Y' &&
								 (letr2=consensus[n-x+y+k+l]) != 'R' && letr2 != 'Y' ) {
								break;
							}

							/* MAKE SURE A MISMATCH IS NOT BEING GIVEN A PASS WHILE NOT CONFORMING TO TRANSITION TYPE */
							if ( (letr=consensus[n-x+y+l]) == 'R' && 
								( (letr2=align2D[m][n+l+k]) == 'C' || letr2 == 'T') ) {
								break;
							}
							else if (letr == 'Y' && 
								( (letr2=align2D[m][n+l+k]) == 'A' || letr2 == 'G')) {
								break;
							}

							if ( (letr=consensus[n-x+y+l+k]) == 'R' && 
								( (letr2=align2D[m][n+l]) == 'C' || letr2 == 'T') ) {
								break;
							}
							else if (letr == 'Y' && 
								( (letr2=align2D[m][n+l]) == 'A' || letr2 == 'G')) {
								break;
							}

							sum4score += MATCH;		/* JUSTIFICATION: MATCHES TRANSITION CALL */
						}
						/* IF LETTER AT n+l EQUAL TO LETTER AT n+l+k except 'n' */
						else if ((align2D[m][n+l] != 'n') && (align2D[m][n+l+k] != 'n'))	
							sum4score += MATCH;
					} /* END OF FOR l SCAN LOOPS */

					if (l == k && 100*sum4score/(k*MATCH) > score_DTHR(k)) {
						imperfect_TR = 1;
					}
					check_imperf = 0;			/* RESET check_imperf HERE */
				} /* END OF IF check_imperf */

				/* BREAK CHECKING IF IT WILL PULL IN MISMATCHES IN REPEAT OR TO RIGHT OF REPEAT */
				if (nuctransit && (imperfect_TR || keep_checking) && n > scrimmage_line) { 
					for (l = 0; l < 2*k; l++) {
						letr2= align2D[m        ][n  +k+l];
						letr = consensus[n-x+k+l];
						if (imperfect_TR) {
							if (isalpha(align2D[m][n+l]) && (i=col_isclear(align2D,n+l,m,-1)) > -1 && 
								align2D[i][n+l] != letr2 && letr != 'R' && letr != 'Y') {
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

				if (keep_checking && tela[symbol_count+k].stat == '-') {
					keep_checking = 0;
				}

				if (keep_checking && n > scrimmage_line) { 
					if ((l=col_isclear(cik_align2D,n-x+k,m,-1)) > -1 
						&& col_isclear(align2D,n+k,m,1)< 0) {

						/* CHECK IF WILL PULL IN ADJACENT MISMATCHES AFTER RUN OF REPEATS */
					    r = 1; 
						i = k;  /* VAR i SET TO k ONLY TO ENTER WHILE LOOP */
						while (i==k) {
							for (i = 0; i < k; i++) {
						    	if (align2D[m][n+i] != (letr2=align2D[m][n+(r+1)*k+i]) && isalpha(letr2))
									break;
					        }    
					        if (i == k)
								r++;        /* INCREMENT NUMBER OF REPEATS */
					    }    

						if (nuctransit) {
							if ((letr=cik_align2D[l][n-x+k+i]) != (letr2=align2D[m][n+(r+1)*k+i]) &&
								 letr!='R' && letr!='Y' && 
								 isalpha(letr) && isalpha(letr2) && (letr3=consensus[n-x+(r+1)*k+i])!='R' && letr3!='Y') {
						        keep_checking = 0; 
						    }    
						}
						else {
							if ((letr=cik_align2D[l][n-x+k+i]) != (letr2=align2D[m][n+(r+1)*k+i]) &&
								 isalpha(letr) && isalpha(letr2)) {
						        keep_checking = 0; 
						    }    
						}
			    	}
				}

				if (imperfect_TR && n > scrimmage_line && 
					isupper(align2D[m][n]) && align2D[m+1][n] == blnk) { 
					for (l = 0; l < k; l++) {
						/* CHECK MISMATCHES FROM PUSHING BOTTOM ROW TO LEFT OF REPEATS AFTER SLIP */
						if ((i=col_isclear(align2D,n+l,m,1)) > -1 &&
							(letr=align2D[i][n+l]) != consensus[n-x-k+l] && 
							 letr != tolower(consensus[n-x-k+l])) {
							imperfect_TR = 0;
							break;
						}
					}
				}

				/*** LAST CHECK TO MAKE SURE NO BAD SLIPS CREATED OUT OF PREVIOUS SLIPS */
				if (keep_checking || imperfect_TR) {
					for (i=n+1; i<lenseq; i++) {
						if (tela[i].k && tela[i].x > n && tela[i].x < n+k && tela[i].y > m+cik_row) {
							keep_checking = imperfect_TR = 0;
							break;
						}
					}
				}

				if (keep_checking || imperfect_TR) {
					if (dev_print(MAIN,__LINE__)) {
						printf("cinch-k taking k-mer=%2d at symbol_count=%3d (lenseq = %3d).", k, symbol_count, lenseq);
					}

					for (l = 0; l < k; l++) {
						cik_align2D[m+cik_row  ][n-x+l] = align2D[m][n+l  ];	
						cik_align2D[m+cik_row+1][n-x+l] = align2D[m][n+l+k];
						x_history[n+l] = x;					/* x_history WRITE-IN FOR NEW TR COLS */
					}
					symbol_count += k;
					cik_align2D[m+cik_row  ][n-x+k  ] = '/';
					cik_align2D[m+cik_row  ][n-x+k+1] = '\0';

 					for (i = 0; i < n-x; i++) {
						if (!isalpha(cik_align2D[m+cik_row+1][i]))
							cik_align2D[m+cik_row+1][i] = blnk;
					}

					/* SCOOCH CONSENSUS ROW IF MINDING TRANSITIONS AND IF BOTTOM IS CLEAR (IS SAFE) */	
					if (nuctransit) {
						if (col_isclear(align2D,n,m,1) < 0) { 
							for (i = n-x; i < n-x+k; i++) {
								if ((letr=consensus[i]) != 'R' && letr != 'Y')
									consensus[i] = consensus[i+k];
							}
							for (i = n-x+k; i+k < lenseq; i++) {
								consensus[i] = consensus[i+k];
							}
							consensus[i] = '\0';
						}
						else if (n >= scrimmage_line) {
							y = y + k;	/* TO KEEP TRACK OF UNSHIFTED CONSENSUS ROW */
						}
					}
					x = x + k;			/* FUTURE SPACING TO BE SUBTRACTED B/C k-MER TUCKED UNDER 1st UNIT */
					x_history[n] = x;
					scrimmage_line = n;
					n = n + k - 1;		/* ADVANCE ADJUSTMENT. NOTE UPCOMING n++ IN FOR n LOOP */
					++cik_row;

				}   /* END OF TR ASSIGN LOOPS */
				else {
					letr = cik_align2D[m+cik_row][n-x] = align2D[m][n];
					x_history[n] = x;
					if (isalpha(letr)) {
						symbol_count++;
					}
				}

			}   /* END OF FOR n LOOPS */ 
		}   /* END OF FOR m LOOPS */

		if (cik_row > 0) {
			mha_writeback(cik_align2D, align2D); 
			printf("\n k = %d:", k);

			if (k > 1)	/* k=1 WILL PRINT FROM MAIN */
				print_2Dseq(options[1][32]);
		}

		options[0][4] = options[0][4] + cik_row;			/* STORE ROWS ADDED */
		if (dev_print(MAIN,__LINE__)) {
			printf("Post-cinch-k for k-mers=%2d: symbol_count=%3d (lenseq = %3d).", k, symbol_count, lenseq);
		}

		if (k > 1)	/* NOT NEEDED AFTER k EQUALS ONE */
			clear_2D_ar(cik_align2D);

	} /* END OF FOR k LOOPS */ 

	mha_UPPERback(cik_align2D, align2D); /* THIS ALSO SAVES 2D-WIDTH in options[1][32] */

	options[1][4] = options[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS i WIDTH HISTORY */
	free(x_history);
/*	free_2D(cik_align2D, lenseq); */

	return (options[0][4]);
}


/*** FUNCTION 04 **********************************************************/
unsigned int cinch_d(short unsigned int cinch_d_opt)
{
int cid_mrow=0, cid_ncol=0, h=0, i=0, j=0, k=WIDTH, l=0, m=0, n=0, num=0, w=0, x=0, tot_repeats=0, uniq_TRs=0, num_transits=0;
int cidwidth = options[1][32]; 
int height = options[1][17];		/* height slot */
unsigned short int nuctype=0, TR_check=0, first_write=1, mono_flag=1;		/* CHECK MONO IN ORDER TO KNOW TO SKIP IT */
unsigned short int nuctransit=0;						/* BIT FLAG FOR HANDLING NUCLEOTIDE TRANSITIONS SILENTLY (IGNORING) */
unsigned short int imperfect_TR=0;
char letr = 'B', ltr2 = 'Z';
char dopt_R_rght = (char) options[1][27];
char blnk        = (char) options[1][11];
char cid_align2D[MAXROW][MAXROW];

	nuctype = options[1][13];		/* EQUALS ONE IF DNA, TWO IF RNA */
	if (nuctype == 1) {		/* IF DNA */
		nuctransit = 1;
	}

	/* START AT BIGGEST k-MER POSSIBLE AT 2x */
	for (k = options[1][32]/2; k > 0; k--) {
		for (n=0; n <= cidwidth-2*k; n++) {	

			mono_flag = 1;			/* MONOMER RUN FLAG IS SET TO 0, WHEN NO LONGER POSSIBLE (ANY n != n+1) */
	
			if (TR_check == 0) 		/* RE-SET COUNTER FOR NUM */
				num = 0;

			if (nuctransit)			/* RE-SET COUNTER FOR NUMBER OF TRANSITIONS */
				num_transits = imperfect_TR = 0;

			for (l=0; l < k; l++) {
				if (isalpha(letr=consensus[n  +l]) && 
					isalpha(ltr2=consensus[n+k+l])) {
					if (l+1<k && letr!=consensus[n+l+1]) 
						mono_flag = 0;		/* CAN NO LONGER BE A HOMOPOLYMER RUN OF FOR THIS k-MER */
					if (letr == 'n' || ltr2 == 'n') 
						break;
					else if (nuctransit && k>2 && letr!=ltr2) {
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
					else if (letr!=ltr2)
						break; 	
				}
				else if (letr=='>' || ltr2=='>')
					break;
				else if (!cinch_d_opt && !isalpha(letr) && isalpha(consensus[n+l+1])) {
					for (i=0; align2D[i][0]!='\0'; i++) {
						for (j=n+l; (letr=align2D[i][j+1])!='\0'; j++) 
							align2D[i][j] = letr;
					}
					for (j=n+l; (letr=consensus[j+1])!='\0'; j++) 
						consensus[j] = letr;
					options[1][32]--;
					break; 
				}
				else if (letr!=ltr2)
					break; 	
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

			/* TESTING IN v3.64: CHECK FOR SUPER-SHORT k=1 HOMODINUCLEOTIDES AND SKIP */
			if (options[1][59]==5 && TR_check && k==1) {		
				if (isalpha(letr=align2D[m][n]) && letr!= align2D[m][n-1] && letr!=align2D[m][n+2]) {
					TR_check = 0;
				}
			}

			/* CHECK FOR COMPLEX cinch_d REPEATS THAT SHOULD NOT BE COUNTED/WRITTEN. */
			/*  THESE HAVE LETTERS IN NEXT ROW UNDERNEATH FIRST UNIT.                */
			if (TR_check) {
				m = 0;
				while (!isalpha(align2D[m][n+k])) {
					m++;
				}
				for (w=1; m+w < height && TR_check != 0; w++) {
					for (x=0; x < n+k; x++) {
						if (x==0 && align2D[m+w][0]=='\0') {
							w = height; 	/* TO BREAK FOR w LOOP */
							break;			/* TO BREAK FOR x LOOP */
						}
						else if (isalpha(align2D[m+w][x])) {
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
					if (consensus[n+l] != consensus[n+num*k+l]) {
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
						while (isalpha(align2D[m][n+k]) == 0) {
							m++;
						}
						if (dev_print(MAIN,__LINE__)) {
							printf("Working on %d-mer consensus TR (%dx) at position %d; 2nd unit begins on row %d.", k, num, n+1, m+1);
						}

						if (imperfect_TR == 1) {
							for (l=0; l < k; l++) {
								letr=consensus[n+l];
								if (letr != consensus[n+k+l]) {
									if (letr == 'A' || letr == 'G') 
										consensus[n+l] = 'R';
									else if (letr == 'C' || letr == 'T') 
										consensus[n+l] = 'Y';
								}
							} 
						}

						mha_writeback(align2D, cid_align2D); 
						cid_align2D[m][n+k  ] = '/';
						cid_align2D[m][n+k+1] = '\0';
						first_write = 0;	/* TURN O-F-F NEED TO WRITE REMAINING PART OF 2-D ALIGNMENT */

						cid_ncol = k;
						cid_mrow = 1;

						/* DEAL WITH LOOSE SLIP CONNECTIONS PRODUCED BY NUDGELIZING */
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
							for (j = n+k; j < options[1][32]+1; j++) {
								letr = cid_align2D[i+cid_mrow][j-cid_ncol] = align2D[i][j];
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
							for (i = n+k; consensus[i+k] != '\0'; i++) {
								consensus[i] = consensus[i+k];
							}
							consensus[i] = '\0';
						} 

						if (letr == '>' && j-cid_ncol-1 < cidwidth) {
							options[1][32] = j-cid_ncol-1;
							mha_writeback(cid_align2D, align2D);
						}
					} /* END OF IF first_write EQUALS ONE */
				} /*********************************************************************************************/
				else if (dev_print(MAIN,__LINE__)) {
					if (imperfect_TR) 
						printf("%4d. i-TR: %3dx %d-mer at consensus position %3d with %d transition(s).", uniq_TRs, num, k, n+1, num_transits);
					else if (nuctransit) 
						printf("%4d. p-TR: %3dx %d-mer at consensus position %3d.", uniq_TRs, num, k, n+1);
					else 
						printf("%4d. TR: %3dx %d-mer at consensus position %3d.", uniq_TRs, num, k, n+1);
				}

			} /* END OF WHILE TR_check EQUALS ONE LOOP */
		} /* END OF FOR n LOOP */
	} /* END OF FOR k LOOP */

	if (cinch_d_opt == 0 && tot_repeats == 0) {
		i = options[1][18];
		options[1][i] = options[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS WIDTH HISTORY */
		printf("\n");
	}
	else if (cinch_d_opt) {
		i = options[1][18];
		options[1][i] = options[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS WIDTH HISTORY */
		if (cidwidth == options[1][32]) {
			print_2Dseq(cidwidth);
			return(0);
		}
		else if (tot_repeats > 1 && options[0][20]) {
			cidwidth = options[1][32];
			options[0][20] = 0;	/* TEMPORARY ASSIGNMENT TO PREVENT PRINTING OF CONSENSUS ROW */
			consensus_2D(0, options[1][32]);
			options[0][20] = 1;	/* REASSIGN SETTING */
		}
		else if (tot_repeats > 1) {
			cidwidth = options[1][32];
			consensus_2D(0, options[1][32]);
		}
		else { 
			cidwidth = options[1][32];
			print_2Dseq(cidwidth);
		}
	}
	return(tot_repeats);
}

/*** FUNCTION 05 ***********************************************************************/
short unsigned int cleanseq(char *s) 
{
int i=0, x=0, length=0, orig_length=0;
int dna=0, rna=0, na=0, eNs=0, prot=0, noprot=0;	/* COUNTERS FOR SIGNATURE CUES OF STRING TYPE */ 
int slop = 12;										/* PERCENT MAGIC SLOP TO ACCOMODATE SOME NON-DIAGNOSTIC SYMBOLS */
short unsigned int stringtype = 0;					/* ALPHA TYPE IS DEFAULT UNTIL FOUND OTHERWISE				*/
short unsigned int nacheck = 1;						/* BIT FLAG. CHANGE TO ZERO WHEN NUCLEIC ACID IS NOT POSSIBLE 	*/
char letr;

	/* MAKE STRING ALL UPPERCASE */
	while ((letr=s[i]) != '\0') {
		if (islower(letr)) 
			s[i] = toupper(letr);
		i++;
	}

	orig_length = i;
	i = 0;

	/* REMOVE NON-LETTERS */
	while ((letr=s[i]) != '\0') {
		if (letr < 'A' || letr > 'Z') {
			x++;	/* EXTRA CHARACTERS REMOVED */
			i++;	/* POSITION */
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
				++eNs;
			}

			s[i-x] = letr;
			i++;
		}
	}

	length = i-x;

	if (prot+noprot > 0)
		nacheck = 0;

	/* ASSIGN STRING TYPE: NON-ZERO VALUES FOR NUCLEIC ACIDS DNA AND RNA */
	if (nacheck) {	/* SUMMARY: U's NECESSARY FOR CALLING RNA, BUT T's NOT NECESSARY FOR CALLING DNA (DEFAULT NA) 	*/
	 	if (rna > 0 && (100*(na+rna+eNs)) >= length*(100-slop)) 	/* 0.90 EQUIVALENT TO 18/20 bp BEING A|C|G|U 	*/
			stringtype = 2;											/* STRING IS LIKELY RNA 						*/
	 	else if (      (100*(na+dna+eNs)) >= length*(100-slop)) 	/* 0.90 EQUIVALENT TO 18/20 bp BEING A|C|G|T 	*/
			stringtype = 1;											/* STRING IS LIKELY DNA 						*/			
	}

	if (stringtype==0 && noprot+rna==0 && prot > 0 && (100*(prot+na+dna+eNs))/length >= 50-slop)	/* 50 B/C COUNTING ~1/2 AMINO ACIDS */
		stringtype = 3;

	/* OVER-WRITE TERMINAL PART OF ORIGINAL STRING WITH NULLS */
	for (i = length; i <= orig_length; i++)
		s[i] = '\0';

	/* IF DNA (STRICT, NOT IUPAC FULL), ALL N'S AND OTHERS TO LOWERCASE 'n' AS AMBIGUOUS DNA */
	if (stringtype == 1 || stringtype == 2) {	/* IF DNA OR RNA */
		i = 0;
		while (s[i] != '\0') {
			if (s[i] == 'N') 
				s[i] = 'n';
			else if (s[i]!='A' && s[i]!='G' && s[i]!='C' && s[i]!='T' && s[i]!='U')
				s[i] = 'n';
			i++;
		}
	}

	return(stringtype);		/* RETURNS 0, 1, 2, OR 3 */

	/*	SIGNATURES OF DIFFERENT ALPHABETS: 
	 	1	A - C - - - G - - - - - - N - - - - - T - - - - - - = DNA 
		2	A - C - - - G - - - - - - N - - - - - - U - - - - - = RNA
		0	A B C D - - G H - - K - M N - - - R S T - V W - Y - = IUPAC DNA
		0	A B C D - - G H - - K - M N - - - R S - U V W - Y - = IUPAC RNA
		3	A - C D E F G H I - K L M N - P Q R S T - V W - Y - = PROTEIN (20 amino acids)
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
int lenseq = options[1][1];

	for (m=0; m <= lenseq; m++) {
		for (n=0; n <= lenseq; n++)
			wipe_align2D[m][n] = '\0';
	}
}

/*** FUNCTION 07 *************************************************************/
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

/*** FUNCTION 08 **********************************************************/
unsigned int consensus_2D(int n_start, int n_width)
{
int badsites=0, m=0, n=0, n_end, x=1;
int con_width = options[1][32];
short unsigned int nuctype = options[1][13];		/* FOR SEQ TYPE, DNA=1, RNA=2, OTHER (NON-NA)=0 */
short unsigned int nuctransit = 0;					/* BIT FLAG FOR HANDLING NUCLEAR TRANSITIONS */ 
short unsigned int plustransit=0;					/* BIT FLAG ADDENDUM FOR COUNTING BADSITES AT COL */
short unsigned int checktransit=0;					/* BIT FLAG FOR CHECKING GOOD TRANSITION MARK */
char blnk = options[1][11], letr=blnk, ltr2=blnk, conletr=blnk;
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
				if (dev_print(MAIN,__LINE__)) {
					printf("checktransit=%d at n=%d.\n", checktransit, n);
				}
			}
		}
		else if (plustransit)
			++consensus_ar[0][n+1];
	}

	if (options[0][20]) {							/* opt_K SHOW CONSENSUS ROW */	
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

	if (badsites > 0) {
		return(badsites);				/* BAD CONSENSUS, REPORT IT 	*/
	}
	else								/* GOOD CONSENSUS, EARLY PASSES */
		return(0);
}

/*** FUNCTION 09 ******consensus_2D NUDGE**********************************/
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
char blnk = options[1][11], letr=blnk, conletr=blnk, chkletr=blnk, badletr=blnk;
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

/*** FUNCTION 10 *** PUSHES CONFLICT RIGHT; CONTAINS LEGACY CYCLELIZE CODE THAT CAN BE DELETED WITH CARE  ******************/
unsigned int nudgelize(void)
{
int cyc_col=0, cyc_row=0, a, b, i, j, kmer=0, m=0, n=0;
int lenseq   = options[1][1];
int cyc_width = options[1][32];						/* THIS IS opt_W SLOT TO STORE CURRENT 2-D WIDTH */
short unsigned int edge0=0;
unsigned short int nuctype = options[1][13];		/* EQUALS ONE IF DNA STRING, TWO IF RNA, THREE IF PROTEIN */
unsigned short int nuctransit=0, dud_nudge=0;		/* BIT FLAG FOR HANDLING NUCLEOTIDE TRANSITIONS SILENTLY (IGNORING) */
unsigned short int tipcyc_flag=0;					/* BIT FLAG FOR TIP CYCLING OPPORTUNITY */
char blnk = options[1][11], letr, conletr, topletr;
char cyc_ar[MAXROW+1][MAXROW] = {{0}};
unsigned int connudge(char con_align2D[][MAXROW], int n_start, int n_width);

	if (nuctype == 1)	/* IF DNA */
		nuctransit = 1;

	mha_writeback(align2D, cyc_ar);
	for (j=0; j<=lenseq; j++)
		cyc_ar[MAXROW][j] = consensus[j];
 
	/* FLAG SPECIAL CASE OF CYCLING NEED AT n=0 COLUMN */
	if (options[0][5] == 0 && cyc_ar[0][0] == blnk)
		edge0 = 1;

	for (n = 0; n <= cyc_width; n++) {

		conletr = align2D[MAXROW][n];
		for (m = 1; align2D[m][0] != '\0' && m <= lenseq; m++) {
			if (isalpha(letr=align2D[m][n])) {
				if (nuctransit && (col_isclear(align2D,n,m,-1) == -1)) {
					topletr = align2D[m][n];
				}
				else if (nuctransit && ( (conletr=='R' && (letr=='A'||letr=='G') && (topletr=='A'||topletr=='G')) ||
								         (conletr=='Y' && (letr=='C'||letr=='T') && (topletr=='C'||topletr=='T'))  )) {
					;	/* NOTHING: GO TO NEXT m */
				}
				else if (letr != cyc_ar[MAXROW][n] || edge0) {
					if (edge0) {
						while (isalpha(align2D[m][0]) == 0 && m <= lenseq) {
							m++;
						}
						cyc_row = m;	/* THIS IS ROW COORDINATE OF NON-CONSENSUS */
						cyc_col = 0;	/* THIS IS COLUMN COORDINATE OF NON-CONSENSUS */
						while (align2D[0][n] == blnk) {	/* SCOOCH RIGHT */
								n++;
						}
						m = 0;			/* NEED TO RESET THIS TO FIRST ROW */
					}
					else {
						cyc_row = m;	/* THIS IS ROW COORDINATE OF NON-CONSENSUS */
						cyc_col = n;	/* THIS IS COLUMN COORDINATE OF NON-CONSENSUS */

						m=0;			/* RE-FIND FIRST LETTER OF COLUMN */
						while (isalpha(align2D[m][cyc_col]) == 0 && m < MAXROW) {
							m++;
						}
	
						/* SOMETIMES NEED TO SKIP ROWS AT BLEEDING EDGE OF SLIPS */	
						if (align2D[m][cyc_col+1] == '/') {
							while (align2D[m][cyc_col+1] == '/') {
								m++;
							}
						}
						else {
							while (align2D[m+1][n] != blnk && align2D[m+1][n] != '\0') {	/* MOVE DOWN ****/
								m++;
							}
	
							while (letr != '0') {
								letr = align2D[m][n];			/* STORE THIS LETTER TO CHECK FOR INTERVENING REPEATS */
								for (i = m+1; i < cyc_row; i++) {
									if (align2D[i][n] == letr) {
										m = i;	/* SKIP TO THIS ROW */
										letr = align2D[m][n];	 
										break;	/* BREAK OUT OF FOR i LOOP & RECHECK W/IN WHILE LOOP */
									}
								}
								letr = '0';		/* JUST USING AS FLAG FOR WHILE LOOP */
							}

							while (align2D[m+1][n] == blnk) {	/* SCOOCH RIGHT */
								n++;
							}
						}
					} /* END OF ELSE (IF NOT edge0) */

					/* TIP-CYCLELIZE AS SOON AS DETECTED */
					tipcyc_flag = 0;
					for (i = m; i < cyc_row; i++) {
						if (align2D[i][cyc_col] != '/')
							break;
					}
					if (i == cyc_row && align2D[m-1][cyc_col+1] == '/') {
						j = 0;
						while (align2D[m][j] == blnk) {
							j++;
						}

						if (align2D[m][j] == align2D[m-1][cyc_col]) {
							tipcyc_flag = 1;		
							kmer = 4;		/* NOT ACTUAL kmer, JUST USING VAR TO CODE CYC TYPE */
							for (a = 0; a < m-1; a++) {
								for (b=0; (letr=align2D[a][b]) != '\0'; b++) {
									cyc_ar[a][b] = letr;
								}
								cyc_ar[m-2][b] = '\0';
							}
							for (b=0; b < cyc_col; b++) {
								cyc_ar[m-1][b  ] = align2D[m-1][b];
							}
								cyc_ar[m-1][b  ] = '/';
								cyc_ar[m-1][b+1] = '\0';
							for (b = 0; b < j; b++) {
								cyc_ar[m  ][b  ] = blnk;
							}
								cyc_ar[m  ][j  ] = align2D[m-1][cyc_col];
								cyc_ar[m  ][j+1] = '/';
								cyc_ar[m  ][j+2] = '\0';
							for (a = m; align2D[a][0] != '\0' && a < MAXROW; a++) {
								for (b=0; (letr=align2D[a][b]) != '\0'; b++) {
									cyc_ar[a+1][b] = letr;
								}
							}
							cyc_ar[MAXROW][cyc_col] = align2D[cyc_row][cyc_col];
							if (dev_print(MAIN,__LINE__)) {
								printf("TIP CYCLING OPPORTUNITY FOR cyc_col = %d; j=%d.", cyc_col+1, j+1);
							}
						}
					}

                    if (kmer != 4)      	/* LEGACY NAMING OF VARIABLE FROM WHEN NUDGELIZE USED TO BE CYCLELIZE; NOW SEPARATE FUNC. */
						kmer = 3;       	/* IF k=9, THEN BELOW WILL FUDGE CYCLELIZE BY PUSHING RIGHT. HACK WORKS FOR ALL k */ 

					if (options[0][5]!=3)		/* GIVING PRECEDENCE TO THE MEMORY OF HAVING NUDGE-CYCLED */
						options[0][5] = kmer;	/* USING THE 0 ROW ABOVE PASS WIDTH ROW TO STORE cyclelize kmer VAR. */

					/* TIP-CYCLELIZE */
					if (tipcyc_flag) {
						/* WILL TIP-CYCLELIZE ABOVE AT FLAG CALL */
						options[1][5] = options[1][32] = cyc_width;	/* ASSIGN [32] CURRENT WIDTH and PASS x WIDTH HISTORY */
						tipcyc_flag = kmer = 0;
						n = cyc_width+1;	/* BREAKS OUT OF n LOOP */
						break;				/* BREAKS OUT OF m LOOP */
					}

					/* NUDGE-CYCLELIZE: */
					if (1) {	
						if (connudge(cyc_ar, 0, cyc_width) == 0) {
							if (dev_print(MAIN,__LINE__)) {
								printf("dud_nudge");
							}
							dud_nudge = 1;
							i = options[1][18];
							options[1][i] = cyc_width = options[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS x WIDTH HISTORY */
							n = cyc_width+1; 		/* BREAK OUT OF FOR n LOOP AFTER BREAKING OUT OF FOR m LOOP */
							break; 					/* BREAK OUT OF FOR m LOOP */
						}

						clear_right(cyc_ar);

						i = options[1][18];
						options[1][i] = cyc_width = options[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS x WIDTH HISTORY */
						n = cyc_width+1; 		/* BREAK OUT OF FOR n LOOP AFTER BREAKING OUT OF FOR m LOOP */
						break; 					/* BREAK OUT OF FOR m LOOP */
					} /* END OF ELSE (IF kmer != 2) */

				} /* END OF IF LETTER != CONSENSUS */
			} /* END OF IF ISALPHA */
		} /* END OF FOR m LOOP */
	} /* END OF FOR n LOOP */

	/* PUSH LEFT IF EMPTY: SHOULD MOVE TO ITS OWN FUNCTION IF NEEDED ELSEWHERE */
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
		mha_writeback(cyc_ar, align2D);
		mha_writeconsensus(cyc_ar, consensus);
		return(print_2Dseq(cyc_width));
	}
	else {
		return (0);
	}
}

/*** FUNCTION 11 **********************************************************/
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

/*** FUNCTION 12 **********************************************************/
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

	if (options[0][41]) {	/* opt_f FOAM-FREE CONSENSUS ROW */	
		line_end(BLOCKHEAD, 9, 9);
		printf(" ");
		for (n = n_start; n < n_end; n++) {
			printf("%c", foam_ar[1][n]); 
		}
		printf(" = foam-free\n");		
	} 

	return(0);
}

/*** FUNCTION 13 ***** RETURN -1 IF ROW DOES NOT HARBOR NEXT FOAM-FREE POSITION *********************/
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
    return(-1);		/* RETURN NEGATIVE VALUE TO INDICATE COLUMN IS CLEAR EITHER BELOW OR ABOVE */
}


/*** FUNCTION 15 **** FOR UNIFORM FORMATTING OF 2-D LINE STARTS AND ENDS ****/
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
	int lcl_opt_r  = options[0][53];

	if (options[0][11]==2 && options[1][18]>1)	/* opt_B LEVELS FOR BLANKNESS IN FILLER & opt_I PASS NUM */
		zero_tick = '|';			/* 124 = '|'							*/

	if (type == 1)					/* FORMAT FOR LINE END. c IS CHAR. NUMBER */
		printf("%5d\n", c);

	if (lcl_opt_r == 1) {			/* FORMAT WITH LINE NUMBERING. DEFAULT */
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
	else if (lcl_opt_r == 0) {	  	/* FORMAT WITHOUT LINE NUMBERS. r = REMOVED */
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

/*****************************************************************************************/
void mha_head(int lcl_width)
{
			 /*0123456789*/
char h1[]=	"\n_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//"
			  "_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//"
			  "_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//"
			  "_MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//__MHA___//";
char h2[]=	"\n\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			  "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			  "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_"
			  "\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_\\____//^_";
char *h_rule = h1;						/* DEFAULT BANNER STYLE */
int min_len = 80;						/* MINIMUM LENGTH */
int med_len = 12*(lcl_width/10);		/* MEDIUM LENGTH, SCALING */
int max_len = options[1][58]+8;			/* MAXIMUM LENGTH */
int hr_len = min_len;					/* DEFAULT LENGTH OF HEADER BANNER */
int lcl_pass = options[1][18];			/* opt_I VALUE COUNTER FOR NUM OF PASSES */

	if (lcl_width > MAXROW && dev_print(MAIN, __LINE__)) {
		printf("Bad news bears: Unexpectedly, lcl_width > MAXROW. lcl_pass=%d, lcl_width=%d.\n", lcl_pass, lcl_width);
	}

	if (lcl_width+8 > hr_len && lcl_width+8 < (int) options[1][58])
		hr_len = med_len;
	else if (lcl_width+8 >= (int) options[1][58])
		hr_len = max_len;

	if (options[0][33])	/* opt_X == 1 */
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
void mha_writeback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW])
{
char letr;
int lenseq       = options[1][ 1];
int wb_mwrap     = options[1][22];		/* opt_M long_homopolymer_run */
char wopt_Q_left = options[1][26];		/* LHS character delimiter for homopolymer Run */
char wopt_R_rght = options[1][27];		/* RHS character delimiter for homopolymer Run */
int m=0, n=0, widest_n=0;

	/* WRITE BACK TO align2D_prev */
	clear_right(lcl_align2D); 
	clear_2D_ar(align2D_prev);
	for (m = 0; lcl_align2D[m][0] != '\0' && m <= lenseq; m++) {
		for (n = 0; (letr=lcl_align2D[m][n]) != '/' && letr != '>' && letr != wopt_R_rght; n++) {
			align2D_prev[m][n] = letr;
			if (letr == wopt_Q_left && lcl_align2D[m][n+wb_mwrap+1] == wopt_R_rght ) 
				align2D_prev[m][n+wb_mwrap+2] = '\0';
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

/*** FUNCTION 21a ******************************************************************************/
void mha_writecons(char align2D_one[][MAXROW], char align2D_two[][MAXROW])
{
int i=0;

	/* COPY CONSENSUS ROW FROM ARRAY ONE TO ARRAY TWO */
	for (i = 0; i < options[1][32]; i++) {
		align2D_two[MAXROW][i] = align2D_one[MAXROW][i];
	}
}

/*** FUNCTION 21b ******************************************************************************/
void mha_writeconsensus(char align2D_one[][MAXROW], char consensus1D[MAXROW])
{
int i=0;

	/* COPY CONSENSUS ROW FROM ARRAY ONE TO 1D CONSENSUS ARRAY */
	for (i = 0; i < options[1][32]; i++) {
		consensus1D[i] = align2D_one[MAXROW][i];
	}
}

/*** FUNCTION 22 ******* LIKE mha_writeback EXCEPT CHECKS TO MAKE UPPERCASE *******************/
void mha_UPPERback(char lcl_align2D[][MAXROW], char align2D_prev[][MAXROW])
{
char letr;
short unsigned 
     int nuctype = options[1][13];				/* EQUALS ONE IF DNA, TWO IF RNA */
int wb_mwrap     = options[1][22];				/* opt_M long_homopolymer_run */
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

			if (letr == wopt_Q_left && lcl_align2D[m][n+wb_mwrap+1] == wopt_R_rght ) 
				align2D_prev[m][n+wb_mwrap+2] = '\0';
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

/*** FUNCTION 23 ********************************/
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

/*** FUNCTION 24-LITE ********************************************************************************************/
short unsigned int print_2Dseq(int print_lenseq2D)
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
int head_start;							/* USE TO PASS RULER O-F-F-SET TO line_end() */
int scrimmageline;						/* USE TO INCREMENT AND TEST IF FILLER IS NEEDED, CAN BE OPTION TO DO SO */
int blnk_lvl = options[0][11];			/* USE TO STORE DEGREE OF BLANKNESS */
int   length = options[1][1];
char tick = ':'; 						/* OTHER POSSIBILITIES: |, ^ */
short unsigned int lcl_opt_F;

	if (options[0][11] > 1)
		tick = blnk;

	if (cinchwidth > MAXROW && dev_print(MAIN,__LINE__)) {
		printf("Bad news bears: Unexpectedly, cinchwidth > MAXROW. cinchwidth=%d", cinchwidth);
	}
	if (print_lenseq2D > MAXROW && dev_print(MAIN,__LINE__)) {
		printf("Bad news bears: Unexpectedly, print_lenseq2D > MAXROW. print_lenseq2D=%d.", print_lenseq2D);
	}
	mha_head(print_lenseq2D);

	blocks2D = count_wrap_blocks(print_lenseq2D, cip_linewidth);

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
				else if ((letr=align2D[m][n]) == blnk && (n+1) % 10 == 0 && blnk_lvl < 2)		
					printf("%c", tick);		/* PRINT TICK MARKS AT 10 bp INTERVALS IF NOT BLANK SPACE */
				else if (options[1][11] == 32)		/* opt_B BLANK = SPACE ' ' */	
					printf("%c", letr);
				else					/* opt_B BLANK = FULLSTOP '.' */
					printf("%c", letr);

				/* TURN ON OPTION opt_F TO USE BLANK CHAR TO FILL < SCRIMMAGE IF EXCESSIVELY SHORT */
				if (options[0][15] != 1 && align2D[m][n] == '/' && scrimmageline-(b+d-2) > 10)
					lcl_opt_F = 1; 	/* TURN ON local opt_F IN THIS CASE */

				/* OPTION opt_F TO USE BLANK CHAR TO FILL < SCRIMMAGE */
				if ((options[0][15] || lcl_opt_F) && align2D[m][n] == '/' && b+d-2 < scrimmageline) {
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
				if (options[0][21])	/* opt_L = LINE END NUMBERING (ON) */
					line_end(END, c, 0);
				else
					printf("\n");
				break;	/* BREAK OUT OF m=m_start FOR LOOP */
			}

			if (options[0][21])	/* opt_L = LINE END NUMBERING ON */
				line_end(END, c, 0);
			else
				printf("\n");
			m++;

		} /* END OF m=m_start FOR LOOP */

		options[1][17] = m+1;	/* ASSIGN COUNTED HEIGHT (# OF ROWS) TO HEIGHT SLOT */

		/* PRINT RULER */
		if (options[0][11] < 5) {
			head_start = (j * options[1][58]) % 10;
			if (j+1 < blocks2D) {
				line_end(RULER, head_start, cip_linewidth);
			}
			else {
				line_end(RULER, head_start, (cip_linewidth=max_n));
			}
		}
		/* *********** */

		/* PRINT NUMBERS FOR CONSENSUS RULER */
		if (options[0][11] < 4) {
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

		if (options[0][41] && options[1][18] > 6 && mmsites == 0) {	
			foam_2D(j_start, cip_linewidth);
		}
	} /* END OF FOR j PRINTING LOOP */

	if (c == length && mmsites == 0) {
		options[0][10] = 1000;	
		printf("\n Successfully 2-D self-aligned all %d letters from formatted sequence. \n", c);
		return(0);
	}
	else if (c == length) {
		i=options[0][10] = round((1000*(cinchwidth-mmsites))/cinchwidth);	
		printf("\n 2-D self-aligned all %d letters from formatted sequence, but consensus indicates more cinching required. \n", c);
		return(i);
	}
	else if (c < length) {
		i=options[0][10] = round((1000*(cinchwidth-mmsites-(length-c)))/cinchwidth);
		warnhead('-');
		printf(" 2-D printing is missing %d letter(s)!\n\n", length-c); 
		return(0);
	}
	else {
		i=options[0][10] = round((1000*(cinchwidth-mmsites-(c-length)))/cinchwidth);	
		warnhead('+');
		printf(" 2-D printing left %d extra letter(s)!\n\n", c-length);
		return(0);
	}
}

/*** FUNCTION 25 *************************************/
void print_base62_table(void)
{
unsigned int i=0;

	mha_head(80);
	printf("\nMHA base 62 single-letter representation:\n\n Base 10: ");
	for (i = 0; i < 31; i++)
		printf("%2d ", i);
	printf("\n Base 62: ");
	for (i = 0; i < 31; i++)
		printf("%2c ", (char) options[3][i]);
	printf("\n\n Base 10: ");
	for (i = 31; i < 62; i++)
		printf("%2d ", i);
	printf(">61\n Base 62: ");
	for (i = 31; i < 62; i++)
		printf("%2c ", (char) options[3][i]);
	printf("  !\n\n");
}

/*** FUNCTION 26 **********************************/
void print_blockhead(int a, int b)              /**/
{                                               /**/
	if (a == 1)                                 /**/
	    printf("   Block %d of %d:\n", a, b);   /**/
	else                                        /**/
	    printf("\n   Block %d of %d:\n", a, b); /**/
}                                               /**/
/**************************************************/

/*** FUNCTION 27 *******************************************************/
int recover_1D(char recovered_1D[MAXROW]) 
{
int m=0, n=0, x=0;
int lenseq = options[1][1];
char rec_R_rght = options[1][27], letr;

	for (m=0; (align2D[m][0] != '\0') || (m < MAXROW); ) {
		for (n=0; (letr=align2D[m][n]) != '/' && letr!=rec_R_rght && letr!='>' && n < lenseq; n++) {
			if (isalpha(letr)) {
				recovered_1D[x] = letr;
				x++;
			}
		} /* END OF FOR n LOOP */
		if (letr=='/' || letr==rec_R_rght)
			m++;
		else if (letr=='>') {
			recovered_1D[x] = letr;
			return(x);	/* RETURNS LENGTH OF RECOVERED_1D[] IN ADDITION TO POPULATING IT W/ SEQUENCE */
		}
	} /* END OF FOR m LOOP */
	return(0);

}

/*** FUNCTION 27b *******************************************************/
int recoverlen(void) 
{
int m=0, n=0, x=0;
int rec_lenseq = options[1][32];
char rec_R_rght = options[1][27], letr;

	for (m=0; (align2D[m][0] != '\0') || (m < MAXROW); ) {
		for (n=0; (letr=align2D[m][n]) != '/' && letr!=rec_R_rght && letr!='>' && n < rec_lenseq; n++) {
			if (isalpha(letr)) {
				x++;
			}
		} /* END OF FOR n LOOP */
		if (letr=='/' || letr==rec_R_rght)
			m++;
		else if (letr=='>') {
			return(x);
		}
	} /* END OF FOR m LOOP */
	return(0);

}

/*** FUNCTION 28 *************************************/
void relax_2D(void)
{
int height=0, i, j, m, n, rlx_col=0, v=0, w=0, z=0;
int width = options[1][32]; 
char blnk = options[1][11], letr;
unsigned short int nuctype = options[1][13], nuctransit=0;
char rlx_opt_R_rght = options[1][27];
char rlx_align2D[MAXROW][MAXROW];

	if (nuctype == 1)		/* IF DNA */
		nuctransit = 1;

	mha_writeback(align2D, rlx_align2D); 
 
	while (align2D[height][0] != '\0') {
		height++;
	}
	options[1][17] = height;

	for (n=0; n < width; n++) {
		v = w = z = m = 0;

		while (align2D[m][n] == '\0') {
			m++;
		}

		while (isalpha(letr=align2D[m][n]) && isalpha(align2D[m][n+1])) {
			n++;
		}

		if (align2D[m][n+1] == '/') {		/* EDGE DETECTED */
			while (isalpha(align2D[m+v][n-1]) && 
				   isalpha(align2D[m+v][ n ]) &&
						   align2D[m+v][n+1] == '/'   ) {
				v++;
			}
			while (align2D[m+v+w][n-1] == blnk && 
				   align2D[m+v+w][ n ] == letr && 
				   align2D[m+v+w][n+1] == '/'     ) {
				w++;
			}
			if (align2D[m+v+w][n-1] == blnk && align2D[m+v+w][n] == letr) 
				w++;		/* TOTAL w EQUALS LENGTH OF HOMOPOLYMER TO RELAX */

			for (i = m+v+w; i < height; i++) {
				for (j = 0; j <= n; j++) {
					if (isalpha(align2D[i][j]))
						w = 0;
				}
			}

			if (w > 0) {
				while (align2D[m+v+w+z][n] == blnk) {	
					z++;
				}
				while (m+v+w+z < height) {
					if (isalpha(align2D[m+w+z][n]))
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
						for (j = n+1; (letr=align2D[i][j]) != '\0'; j++) 
							rlx_align2D[i-rlx_col][j+rlx_col] = letr;
					}

					if (nuctransit) {
						j = options[1][32] + rlx_col + w;
						consensus[j] = '\0';
						while (j > n+rlx_col-w) {
							consensus[j] = consensus[j-w];
							j--;
						}
						for (j = n+rlx_col-w+1; j <= n+rlx_col; j++) 
							consensus[j] = '\0';
					}

				}
			} /* END OF IF w > 1 */
		}
		else if (align2D[m][n+1] == '>' || align2D[m][n+1] == rlx_opt_R_rght) {
			rlx_align2D[m - rlx_col][n + rlx_col] = letr;
			rlx_align2D[m - rlx_col][n + rlx_col + 1] = align2D[m][n+1];
			rlx_align2D[m - rlx_col][n + rlx_col + 2] = '\0';

			if (align2D[m][n+1] == '>') {
				options[1][32] = n + rlx_col + 1;
			    align2D[m][n+2] = '\0';
			}

		}

	} /* END OF FOR n LOOP */

	mha_writeback(rlx_align2D, align2D);

	i = options[1][18];
	options[1][i] = options[1][32];	/* ASSIGN [32] CURRENT WIDTH and PASS [9] WIDTH HISTORY */
}

/*** FUNCTION 29 **************************************/
void usage(char usage_version[], unsigned int FY_size)
{
	printf("\nRunning maximal version %s, a program for micro-homology alignment (MHA).\n", usage_version);
	printf("\nUsage: ./maximal -[OPTIONS] sequence.txt (OPTIONAL FASTA HEADER AND NON-ALPHA SEQUENCE CHARACTERS ARE IGNORED.)\n"
							"\t\t -c       USE REVERSE COMPLEMENT.\n"
							"\t\t -f       SHOW FOAM-FREE SEGMENTS (REQUIRES ALLOWING DEFAULT RELAX OPTION).\n"
							"\t\t -h       SHOW HELP, USAGE.\n"
							"\t\t -k       SHOW k-MER COUNT.\n"
							"\t\t -l       SHOW CINCH-T SLIP LOCATIONS, NUMBER (BASE 62), AND K-MER SIZES (BASE 62).\n"
							"\t\t -n       DO NOT DO RELAX-2D PASS.\n"
							"\t\t -o       SHOW ORIGINAL INPUT STRING.\n"
							"\t\t -p       SHOW PARAMETERS.\n"
							"\t\t -r       SHOW ROW NUMBER.\n"
							"\t\t -s       SILENCE WRITING TO NORMAL OUTPUT FILE (Surf_wavereport.mha).\n"
							"\t\t -t       SKIP CINCH-T MODULE, PROCEED TO REMAINING MODULES.\n"
							"\t\t -u(u...) DO NOT WRAP OUTPUT (WRAP IN ONE BLOCK); EACH '-u' WRAPS OUTPUT INTO ANOTHER BLOCK.\n"
							"\t\t -v(v)    VERBOSE MODE: \"maximal -FLIRlopr\" + VERBOSITY. TOGGLE off I & R WITH -vI & -vR.\n"
							"\t\t -x       EXTRA SQUEEZE: REDUCE CINCH-T THRESHOLDS FOR TRANSITION MATCHING BY ONE.\n"
							"\t\t -z       ENSURE MISMATCH SCORE IS ZERO; ALTERS PATHBOX DISPLAY.\n"
							"\t\t -B       USE EMPTY SPACE FOR BLANK CHARACTER INSTEAD OF PERIOD.\n"
							"\t\t -B(BBBB) DO NOT SHOW TICK MARKS, ZERO TICKLINE, RULER #'s, & RULER, RESPECTIVELY.\n"
							"\t\t -C       SHOW BASE 62 SINGLE-LETTER CODE (USED FOR REPEAT #'s).\n"
							"\t\t -D       SHOW DIAGNONAL THRESHOLD VALUES FOR TRANSITIONS HANDLING.\n"
							"\t\t -F       DO NOT FILL SHORT LINE ENDS TO SCRIMMAGE LINE.\n"
							"\t\t -H       SHOW HELP, USAGE.\n"
							"\t\t -I       SHOW INITIAL PASSES.\n"
							"\t\t -K       SHOW CONSENSUS ROW.\n"
							"\t\t -L       SHOW POSITION AT END OF EACH ROW.\n"
							"\t\t -M       DOUBLE THE DEFAULT LONG HOMOMONOMER TRACT WRAP LENGTH.\n"
							"\t\t -O       SAVE (APPEND) FORMATTED 2-D ALIGNMENT TO CONSENSUS FILE (Surf_barrels.log).\n"
							"\t\t -OO      SAVE (APPEND) RAW 2-D ALIGNMENT TO SPECIAL 2-D INPUT FILE (TUBES.mha).\n"
							"\t\t -P       SHOW PATHBOX.\n"
							"\t\t -R       RECOVER AND CHECK 1-D SEQUENCE FROM 2-D SELF-MHA.\n"
							"\t\t -X       RUN ON SCRAMBLED SEQUENCE OF SAME LENGTH.\n"
							"\t\t -XX      RUN ON FISHER-YATES SHUFFLED SEQUENCE OF LENGTH %d.\n"
							"\t\t -Y       USE NUMBER ARGUMENT AS FY_SIZE INSTEAD OF DEFAULT (%d).\n\n", FY_size, FY_size);
	printf("Example usages: ./maximal -a\n");
	printf("                ./maximal -KnO sequence-file.txt\n");
	printf("                ./maximal -KnXXY 1200 sequence-file.txt\n\n");
}

/*** FUNCTION 31 *************************************/
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
	short unsigned int devReport=options[1][57];

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
	if (devReport>2) {	/* CODE DEVELOPMENT REPORTING */
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
	if (devReport>2)
		printf("%c", mha_base62(bottom[0]));
	for (j=1; j <= width; j++) {
		if ((i=heightAB - bottom[j] - top[j]) < bottom[j-1]) {
			bottom[j] = i;
		}
		else
			bottom[j] = bottom[j-1];
		if (devReport)
			printf("%c", mha_base62(bottom[j]));
	}
	if (devReport>2)
		printf(" <-- MIN(heightAB - bottom edge - top edge)\n");
	i = heightA - bottom[j-1] + 1;
	if (devReport>2)
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
		if (options[1][13] && options[1][57]) { /* IF DNA AND VERBOSITY */
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
			if (options[0][21]) 
				line_end(END, n, 0);
			else
				printf("\n");
		}
		else {
			printf(" ");
			if (options[0][21]) 
				line_end(END, n, 0);
			else 
				printf("\n");
		}

		if (options[0][47]) {		   /* OPTION opt_l TO SHOW SLIP LOCATIONS */
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

/*********************/
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


/**************************************************************** END OF MHA fUNCTIONS **********************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal			*/
