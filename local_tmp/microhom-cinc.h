/******************************************************************/
/* microhom-cinc.h header file, since mha_v4.26.c                 */
/* This file has MHA the original cinch modules called in         */
/* sequence as listed below.                                      */
/******************************************************************/

#ifndef FILE_CINC_SEEN
#define FILE_CINC_SEEN

short unsigned int 	cleanseq(char *s);
/*					cinch_t = main() */
int 				cinch_l(void); 
int 				cinch_k(short unsigned int mode);
unsigned int 		nudgelize(void);
unsigned int 		cinch_d(short unsigned int cinch_d_opt);
void 				relax_2D(void);
int 				recover_1D(char *recovered_1D);
void 				clear_pathbox(void);

/******************************************************************/
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
				s[i] = ambig.sym;
			else if (s[i]!='A' && s[i]!='G' && s[i]!='C' && s[i]!='T' && s[i]!='U')
				s[i] = ambig.sym;
			i++;
		}
	}
	else if (stringtype == 3 && opt_x.bit) {
		i = 0;
		while ((letr=s[i]) != '\0') {
			if (letr=='L' || letr=='V') 					/* WORKHORSE HYDROPHOBICITY,  PT. 1 */
				s[i] = 'O';
			else if (letr=='S' || letr=='T') 				/* -OH 		*/
				s[i] = 'B';
			else if (letr=='E' || letr=='D' || letr=='Q') 	/* POLAR -, OR CARBOXYMIDE STRUCTURE-FRIENDLY (Q) */
				s[i] = 'U';
			else if (letr=='R' || letr=='K') 				/* POLAR +	*/
				s[i] = 'Z';
			else if (letr=='I' || letr=='M') 				/* WORKHORSE HYDROPHOBICITY, PT. 2  */
				s[i] = 'O';
			else if (letr=='G' || letr=='A') 				/* SMALL 	*/
				s[i] = 'J';
			else if (letr=='F' || letr=='Y' || letr=='W') 	/* HEAVY HYDROPHOBIC MOTHERS */
				s[i] = 'X';
			i++;
		}
		print_protein_waxes();
	}

	return(stringtype);		/* RETURNS 0, 1, 2, OR 3 */

	/*	PROFILES OF DIFFERENT ALPHABETS: 
	 	1	A - C - - - G - - - - - - N - - - - - T - - - - - - = DNA 
		2	A - C - - - G - - - - - - N - - - - - - U - - - - - = RNA
		0	A B C D - - G H - - K - M N - - - R S T - V W - Y - = IUPAC DNA
		0	A B C D - - G H - - K - M N - - - R S - U V W - Y - = IUPAC RNA
		3	A - C D E F G H I - K L M N - P Q R S T - V W - Y - = PROTEIN (20 amino acids)
		3x	J - C U U X J H O - Z O O N - P U Z B B - O X - X - = DOC OZ'S BLOSUM-90 RAREFIED WAXES
		0	A - C D E F G H I - K L M N O P Q R S T U V W - Y - = PROTEIN (22 a.a. w/ pyrolysine=O and selenocysteine=U)
		0	- - - - E F - - I - - L - - - P Q - - - - - - - - - = PROTEIN-ONLY (THESE 6 amino acids are diagnostic)
		0	- - - - - - - - - J - - - - O - - - - - - - - X - Z = NOT CANONICAL PROTEIN, NOT-IUPAC 
	 	0	A B C D E F G H I J K L M N O P Q R S T U V W X Y Z	= ALPHABET
	*/
}
/******************************************************************/

void print_protein_waxes(void)
{
		printf("\nProtein sequence waxed into BLOSUM-90 equivalence classes:\n");
		printf(" O = [LVIM]\t Workhorse hydrophobicity\n");
		printf(" X = [FYW] \t Heavy hydrophobic mothers\n");
		printf(" U = [EDQ] \t Polar, negative or structure-friendly carboxymide (Q)\n");
		printf(" Z = [RK]  \t Polar, positive\n");
		printf(" B = [ST]  \t Small -OH\n");
		printf(" J = [GA]  \t Small\n");
}

/** INVOKE PRIOR TO RE-USING PATHBOX **/
void clear_pathbox(void)
{
	int i=0, j=0;

	for (i=0; i<=Clean.pass_W; i++)
		for (j=0; j<=Clean.pass_W; j++)
			pathbox[i][j] = '\0';

	for (j=0; j<=Clean.pass_W; j++)
		pathbox[MAXROW][j] = '\0';
}

/******************************************************************/
int cinch_l(void) 
{
int cil_row=0, i=0, j=0, k=0, l=0, m=0, n=0, run=0, x=0;
char letr;
char blnk = Fill->sym;		/* opt_B blank character		*/

	clear_pathbox();
	for (m = 0; align2D[m][0] != '\0'; m++) {
		for (n = 0; align2D[m][n] != '\0'; n++) {
			while (align2D[m][n] == blnk) {	/* MOVE WINDOW PAST BLANKS */
				pathbox[m+cil_row][n-x] = blnk;
				n++;			  
			}
			/* letr ASSIGNED */
			if ((letr=align2D[m][n]) != slip.sym || letr != '\0') {	/* letr != '\0' USEFUL IF BLANKS FILLED TO SCRIMMAGE */ 
				pathbox[m+cil_row][n-x] = letr;	
				if (n>1 && letr==align2D[m][n-1])
					run++;
				else
					run = 1;
			}
			else
				pathbox[m+cil_row][n-x] = letr;		/* WRITE LINE TERMINAL MHA CHARACTERS */
	
			if (run == 2*opt_M.val) {			/* TRIGGER LENGTH MEASURE OF MONOMERIC RUN	*/
												/* AND WRITING OF (10)-MER BLOCKS			*/
				++cinchled;
				while (align2D[m][n+run-2*opt_M.val] == letr) {
					run++;
				}
				
				pathbox[    m+cil_row][n-x-opt_M.val+1  ] = monoR.sym;	/* WRITE SLIP AFTER FIRST 10 */
				for (l = 1; l < opt_M.val; l++)								
					pathbox[m+cil_row][n-x-opt_M.val+1+l] = '\0';
				cil_row++;														/* ADVANCE TO NEXT ROW, CONT.*/
	
				for (j = (run-opt_M.val)/opt_M.val; j > 0; j--) {
					for (i = 0; i < n-x-2*opt_M.val; i++)
						pathbox[m+cil_row][i] = blnk;
					pathbox[m+cil_row][i] = monoL.sym;					/* MARK LEFT EDGE OF MWRAP RUN */
					for (k = 0; k < opt_M.val; k++) {
						pathbox[m+cil_row][n-x-2*opt_M.val+1+k] = letr; /* FILL W/ MONOMER LETTER  */
					}
					if (j > 1) {
						pathbox[m + cil_row  ][n-x-opt_M.val+1] = monoR.sym;
                        pathbox[m + cil_row++][n-x-opt_M.val+2] = '\0';             /* SHOULD NOT BE NECESSARY */
					}
					else if (j == 1) 
						break;			/* BREAK OUT OF FOR j LOOP */
				}

				for (i = n-opt_M.val; i <= Current.pass_W; i++)
					consensus[i] = consensus[i+opt_M.val];

				n = n + opt_M.val*(run/opt_M.val) - 2*opt_M.val - 1;		/* ADVANCE n TO JUST PAST LAST 10-MER BLOCK */
				x = x + opt_M.val*(run/opt_M.val) - opt_M.val;				/* INCREMENT x TO REFLECT TUCK */
				
				run = 1;	/* RESET TO AVOID RETRIGGERING WRAPPING FOR REMAINDER OF RUN < mwrap */
			}	/* END OF MONOMERIC RUN BLOCK WRITES */ 
		}   /* END OF FOR n LOOPS */ 
	}   /* END OF FOR m LOOPS */
	Cinch_L.pass_W = Current.pass_W;

	if (cinchled) {
		mha_writeback(pathbox, align2D); 
	}
	return(cinchled);
}
/******************************************************************/


/******************************************************************/
int cinch_k(short unsigned int mode) 
{

	if (!mode) {
		Cinch_K.pass_W = Current.pass_W;	/* ASSIGN CURRENT WIDTH and PASS i WIDTH HISTORY 			*/
		return(Cinch_K.pass_V);				/* RETURNS 0; Cinch_K.pass_V holds cumulative cinch-k runs. */
	}

	int cik_row=0, i=0, k=0, l=0, m=0, n=0, scrimmage_line = -1, x=0, y=0;
	int first_mwrap_start=0, last_mwrap=0;
	unsigned short int first_mwrap=0, keep_checking=1;
	unsigned short int nuctype = Clean.pass_V;			/* EQUALS ONE IF DNA STRING, TWO IF RNA, THREE IF PROTEIN */
	unsigned short int nuctransit=0, check_imperf=0;	/* BIT FLAG FOR HANDLING NUCLEOTIDE TRANSITIONS SILENTLY (IGNORING) */
	unsigned short int homopolyflag=0, imperfect_TR=0;
	char letr, letr2, letr3;
	char blnk = Fill->sym;				/* opt_B fill character */
	int max_k = WIDTH/2;
	int lenseq = Clean.pass_W;
	int tela_m = 0, tela_n = 0;
	int *x_history = NULL;
	char unshifted[MAXROW] = {0};

	for (i=0; i<=Current.pass_W; i++)
		unshifted[i] = consensus[i];

	x_history = (int *)calloc(lenseq, sizeof(int));
	clear_pathbox();

	if (nuctype == 1) {		/* IF DNA */
		nuctransit = 1;
	}	
	if (dev_print(CINCH,__LINE__)) {
		printf("Post cinch-t max_k = %d. Cinch_T.pass_Q = %d", max_k, Cinch_T.pass_Q);
	}
	if (Cinch_T.pass_V<1 || mode==1)
		max_k = 1;

	if (dev_print(CINCH,__LINE__)) {
		printf("Using max_k = %d.\n", max_k);
	}

	/* START AT BIGGEST k-MER POSSIBLE AT 2x */
	for (k = max_k; k > 0; k--) {
		cik_row = x = tela_m = 0;
		scrimmage_line = -1;

		for (m = 0; align2D[m][0] != '\0'; m++) {
			y = 0;

			if (cinchled) {
				/* CHECK LINE AHEAD OF TIME FOR FIRST MONO-RUN TERMINATOR */
				for (i = 0; (letr=align2D[m][i]) != '\0'; i++) {	
					/* NOTE: APPARENT REDUNDANCY BELOW ALLOWS CHANGING L/R RUN DELIMITERS TO BE EQUAL */
					if (letr == monoR.sym && align2D[m][i-1] != blnk) {
						first_mwrap = 1;					/* TURN ON BIT FLAG FOR HANDLING FIRST MONO-RUN */
						last_mwrap = 1;						/* TURN ON BIT FLAG FOR HANDLING TO LAST MONO-RUN */
						first_mwrap_start = i-opt_M.val;	/* SAVE POSITION OF START OF MWRAP */
						break;								/* BREAK OUT OF FOR i CHECK LOOP */
					}
				}
			}
		
			for (n = 0; align2D[m][n] != '\0'; n++) { 
				keep_checking = 1;			/* THIS FLAG HANDLES THE CONTINUED NEED TO CHECK FOR INTRA-TR REPEATS   */
				imperfect_TR = 0;			/* THIS FLAG IS TURNED ON (SET TO ONE) WHEN TR W/ TRANSITION MISMATCHES IS FOUND */
				tela_n = tela_m + k;	/* SHORTER AND AVOIDS HAVING TO DO THIS SUM MANY TIMES */

				if (!n && isalpha(align2D[m][0])) {
					x = 0;												/* x-VAR INITIALIZED, TOP OF FOR n LOOP */
					x_history[0] = x;
				}

				/* MOVE WINDOW PAST INITIAL BLANKS */
				while (align2D[m][n] == blnk) {
					pathbox[m+cik_row][n-x] = blnk;
					n++;
				}

				/* CHECK FOR NEED TO ADD ADDITIONAL BLANKS & ADJUST x */
				if (n <= scrimmage_line) {								/* x-VAR PROXIMAL B.S. FIX */
					for (i = n-x; i < scrimmage_line; i++)
						pathbox[m+cik_row][i] = blnk;
					if (n > 1) {
						x = x_history[n-1];
					}
					else
						x = x_history[0];
				}

				if (cinchled) {
					/* CHECK FOR & DEAL WITH LONG HOMOPOLYMER RUN WRAPS (1ST ONE OR SUBSEQUENT ONES) */
					if (first_mwrap && n == first_mwrap_start) {		
						tela_m += opt_M.val;
						for (n = first_mwrap_start; (letr=align2D[m][n]) != '\0'; n++) {
							pathbox[m+cik_row][n-x] = letr;
							if (isalpha(letr)) {
								x_history[n] = x;					/* x_history WRITE-IN FOR MWRAP DONE ONCE */
							}
						}
						pathbox[m+cik_row][n-x] = '\0';	
						first_mwrap = 0;						/* RESET first_mwrap FLAG */
						break;									/* BREAK OUT OF FOR n LOOP */
					}
					else if (last_mwrap) {
						if (align2D[m][n] == monoL.sym && align2D[m][n+opt_M.val+1] == monoR.sym) {		
							tela_m += opt_M.val;
							pathbox[m+cik_row][n-x] = monoL.sym;
							n++;
							while ( (letr=align2D[m][n]) != monoR.sym) {
								pathbox[m+cik_row][n-x] = letr;
								n++;
							}
							pathbox[m+cik_row][n-x] = letr;			/* WILL WRITE monoR.sym */
							for (i = opt_M.val; i > 0; i--)				/* FAILSAFE:			*/
								pathbox[m+cik_row][n-x+i] = '\0';	/* OVERWRITE W/ NULLS	*/
							break;										/* BREAK OUT OF n LOOP	*/
						}
						else if (align2D[m][n] == monoL.sym && align2D[m][n+opt_M.val+1] != monoR.sym) {
							tela_m += opt_M.val;
							last_mwrap = 0;
							for (i = 0; i <= opt_M.val; i++) {
								pathbox[m+cik_row][n-x+i] = align2D[m][n+i];
							}
							letr = align2D[m][n+1];
							n = n + opt_M.val;
							while (align2D[m][n] == letr) {
								tela_m++;
								pathbox[m+cik_row][n-x] = letr;
								x_history[n] = x;
								n++;
							}
							tela_m--;		/* EACH BLOCK OF MWRAPS OVERCOUNTS BY ONE AFTER THE LAST WHILE LOOP */
						}
					}
				} /* END OF A CINCHLD BLOCK */

				/* CHECK FOR & DEAL WITH LINE ENDS TOO SHORT TO HARBOR TR OF SIZE k */
				if (align2D[m][n+2*k] == Term->sym && col_isclear(align2D, n,m,-1)>=0)
					keep_checking = 0;
				else if (k==1 && align2D[m+1][n]==monoL.sym)
					keep_checking = 0;

				if (!isalpha(align2D[m][n+2*k-1])) {	/* TRUE IF WINDOW < 2x k-MER, WRITE REST OF LINE TO pathbox */
					for (i = n; (letr=align2D[m][i]) != '\0'; i++) {
						pathbox[m+cik_row][i-x] = align2D[m][i];
						if (isalpha(letr)) {
							tela_m++;
							x_history[i] = x;				/* x_history WRITE-IN FOR LINE ENDS TOO SHORT FOR TR */
						}
					}
					break;		/* BREAK OUT OF FOR n LOOP */
				}

				if (nuctransit && keep_checking && k<=opt_b.val) {
					/* SKIP IF CANNOT BE CINCHED B/C OF SUB-THRESHOLD FRACTALS AT ANY PARALOGOUS POSITION */
					for (i=0; i<k; i++) {
						if (tela[tela_m+i].t!=tela[tela_m+i].c || tela[tela_n+i].t!=tela[tela_n+i].c) {
							keep_checking = 0;
							break;
						}
					}

					for (l=n; l<n+2*k; l++) {
						if ((letr=unshifted[l])=='R' || letr=='Y')
							keep_checking = 0;
					}
				}

				/* CHECK FOR TR OF SIZE k-MER */
				if (keep_checking) {
					if (k>1)
						homopolyflag = 1;	/* SET HOMOPOLYFLAG STATUS TO UNKNOWN/NEED TO CHECK */
					else if (nuctransit) {
						if (tela[tela_m].t!=tela[tela_m].c || tela[tela_n].t!=tela[tela_n].c) {
							keep_checking = 0;
						}
					}

					for (i = 0; i < k; i++) {
						/* BREAK EARLY IF WOULD-BE TRANSVERSION */
						if ( ((letr =align2D[m][n+i  ])=='A' || letr =='G') &&
							 ((letr2=align2D[m][n+i+k])=='C' || letr2=='T') ) {
							keep_checking = 0;
							break;
						}
						else if ( ((letr =align2D[m][n+i  ])=='C' || letr =='T') &&
							      ((letr2=align2D[m][n+i+k])=='A' || letr2=='G') ) {
							keep_checking = 0;
							break;
						}
						else {
							/* CHECK TO SEE IF THERE ARE n's */
							if (nuctype && (align2D[m][n+i]==ambig.sym || align2D[m][n+i+k]==ambig.sym)) {
								keep_checking = 0;
								break;
							}
	
							/* CHECK TO SEE IF HOMOPOLYMER RUN CAN BE EXCLUDED */
							if (homopolyflag && i > 0) {	/* IF i > 0, THEN homopolyflag IS 1 */
								if ((align2D[m][n + i] != align2D[m][n+k+i  ]) ||
								    (align2D[m][n + i] != align2D[m][n + i-1]) ||
								    (align2D[m][n+k+i] != align2D[m][n+k+i-1]) ) {
									homopolyflag = 0;
								}
							}
							if (nuctransit && keep_checking) {
								if ((letr=tela[tela_m+i].t)=='R' || letr=='Y')
									keep_checking = 0;
								else if (!y && ((letr3=consensus[n-x+k+i])=='R' || letr3=='Y'))
									keep_checking = 0;
							}

							if (keep_checking && align2D[m][n+i] != align2D[m][n+k+i]) {
								keep_checking = 0;
								if (nuctransit && k >= opt_b.val) {
									check_imperf = 1;
								}
							}
						}
					}

					if (homopolyflag) 		/* IF HOMOPOLYFLAG=1 WAS NOT SET TO ZERO */
						keep_checking = check_imperf = 0;

					homopolyflag = 0;		/* RESET */
				} 

				/* 9/9/2020 v4.35 mucho chowder but dramatic drop in WCR avg's */
				if (nuctransit && keep_checking) {
					if (k==1 && (col_isclear(align2D,n  ,m,1)>0 || col_isclear(align2D,n  ,m,-1)>=0) &&
						        (col_isclear(align2D,n+1,m,1)>0 || col_isclear(align2D,n+1,m,-1)>=0)) {
						keep_checking = 0;
					}
					else if (k>2 && count_unique_chars(align2D[m]+n, 2*k)==2) {
						keep_checking = 0;
					}
				}

				/* THIS BLOCK SPOTS NON-FRACTAL TRs AT NEXUS OF TWO OVERLAPPING AND/OR ABUTTING TRs AND SKIPS THEM */
				/* 9/9/2020 one extra chowder bit at b=3 and some other wiggle when put in OFF */
				if ((keep_checking || check_imperf) && col_isclear(align2D,n      ,m, 1)<0 
													&& col_isclear(align2D,n+2*k-1,m,-1)<0) {

					int case_X = 1;
					if (align2D[m][n]!=align2D[m-1][n]) 
						case_X=0;
					if (case_X) {
						i=n+2*k-2;
						if (align2D[m][i]!=align2D[m+1][i]) 
							case_X=0;
					}
					if (case_X) {
						keep_checking = check_imperf = 0;
					}
				}

				/* 9/09/2020 v4.35 expansion of WCR avg's for b=3 benchmarks w/o chowder when OFF */
				/* 9/17/2020 v4.35 same about expansion but not squirt_17 is dependent on this, plus this is main mechanism to squash things in cinch_k */
				if ((keep_checking || check_imperf) && k>1 && tela[tela_n].echoes==cyc_skip.sym && tela[tela_n].statf!=st_fract.sym) {
					keep_checking = check_imperf = 0;
				}

				/* CHECK TO MAKE SURE NO BAD SLIPS CREATED OUT OF PREVIOUS SLIPS */
				/* 9/9/2020 v4.35, chowder in more than one benchmark test script when this block is OFF */
				if ((keep_checking || imperfect_TR) && k>1 && isalpha(letr=align2D[m+1][n+k-1]) && !isalpha(align2D[m+1][n])) {
					keep_checking = imperfect_TR = 0;
				}

				/* CHECK TO MAKE SURE NO BAD SLIPS CREATED OUT OF PREVIOUS SLIPS */
				/* 9/9/2020 mucho chowder in cleanup_set-all when this block is OFF; did not test rest */
				if ((keep_checking || imperfect_TR) && k>1 && isalpha(letr=align2D[m-1][n+k])) {
					
					int toprow = m-1;
					int j = k;	/* DO THIS TO ENTER WHILE LOOP FIRST TIME */
					while (j==k && toprow>=0) {
	
						for (j=0; j<k; j++) {
							if (align2D[toprow][n+k+j]!=align2D[m][n+j])
								break;
						}
						if (j==k) {
							if (isalpha(letr=align2D[toprow-1][n+k]))
								toprow--;
							else
								break;
						}
					}
					if (j<k) {
						int checkcol = n+2*k;
						for (i=m; i<=lenseq; i++) {
							if (isalpha(letr2=align2D[i][checkcol]) && letr2 != letr) {
								keep_checking = imperfect_TR = 0;
								break;
							}
							else if (align2D[i][0]=='\0') {
								keep_checking = imperfect_TR = 0;
								break;
							}
						}
					}
				}

				/* Handles block of cinching fractal TRs in the first row if they overlay cryptic overlapping TRs in lower rows; churly11 is index case */
				/* 9/9/2020 v4.35, churly11-13 are only strings for which this block matters; 9/28/2020 churly16 too  */
				if (keep_checking||imperfect_TR) {
					i = m+1;
					int j = n+k-1;
					while (align2D[i][n]!=blnk || !isalpha(align2D[i][j]) || col_isclear(align2D,j,i,1)) {
						if (col_isclear(align2D,j,i+1,1)<0) {
							i++;
							break;
						}
						else
							i++;
					}
					if (align2D[i][0] && align2D[i][n]==blnk && isalpha(align2D[i][j])) {
						keep_checking = imperfect_TR = 0;
					}
				}

				/**************************************************************************************************/
				if (keep_checking || imperfect_TR) {
					if (k>0 && dev_print(CINCH,__LINE__)) {
						printf("cinch-k taking k-mer=%2d at tela_m=%3d; m=%2d, n=%2d, x=%d, y=%d.", k, tela_m, m,n,x,y);
					}

					push_gPnt_kmer(tela_n,k,1);

					for (i=0; i<k; i++) {
						pathbox[m+cik_row  ][n-x+i] = align2D[m][n+i  ];
						pathbox[m+cik_row+1][n-x+i] = align2D[m][n+i+k];
						x_history[n+i] = x;					/* x_history WRITE-IN FOR NEW TR COLS */
					}

					/* UPDATE TELA FOR CINCH-K CINCH; SHOULD BE OWN FUNCTION? */
					if (ON) {
						for (i=tela_n; i<=lenseq; i++) 
							++tela[i].y;

						int row = m+cik_row+1;
						for (i=tela_n; i<= lenseq && tela[i].y==row; i++) 
							tela[i].x -= k;

						if (col_isclear(align2D,n+k,m,1)<0) {	
							for ( ; i<=lenseq; i++) 
								tela[i].x -= k;
						}
					}

					pathbox[m+cik_row  ][n-x+k  ] = slip.sym;
					pathbox[m+cik_row  ][n-x+k+1] = '\0';

 					for (i = 0; i < n-x; i++) {
						if (!isalpha(pathbox[m+cik_row+1][i]))
							pathbox[m+cik_row+1][i] = blnk;
					}

					/* SCOOCH CONSENSUS ROW IF MINDING TRANSITIONS AND IF BOTTOM IS CLEAR (IS SAFE) */	
					if (nuctransit) {
						if (col_isclear(align2D,n,m,1) < 0) { 
							for (i = n-x; i < n-x+k; i++) {
								if ((letr=consensus[i]) != 'R' && letr != 'Y')
									consensus[i] = consensus[i+k];
							}
							for (i = n-x+k; i <= Current.pass_W; i++)
								consensus[i] = consensus[i+k];
						}
						else if (n >= scrimmage_line) 
							y += k;		/* TO KEEP TRACK OF UNSHIFTED CONSENSUS ROW */
					}

					x += k;				/* x-VAR UPDATED. FUTURE SPACING TO BE SUBTRACTED B/C k-MER TUCKED UNDER 1st UNIT */
					for (i=n; i<=n+k; i++)
						x_history[i] = x;

					if (tela[tela_n].statf==st_fract.sym) {		/* CHURLY-25 ONLY? 7/18/2020 v3.44 series; but also speed bump */
						for (i=tela_n; i<tela_n+k; i++) {
							tela[i].echoes = cyc_skip.sym;
						}
					}

					scrimmage_line = n;
					tela_m = tela_n;		/* SAME AS "tela_m += k" */
					n += k-1;				/* -1 BECAUSE UPCOMING n++ IN FOR n LOOP */
					cik_row++;

				}   /* END OF TR ASSIGN LOOPS */
				else {
					letr = pathbox[m+cik_row][n-x] = align2D[m][n];
					x_history[n] = x;
					if (isalpha(letr)) {
						tela_m++;
					}
				}

			}   /* END OF FOR n LOOPS */ 
		}   /* END OF FOR m LOOPS */

		if (cik_row > 0) {
			printf("\n Next: cinch-k for k = %d...\n", k);
			mha_writeback(pathbox, align2D); 
			print_2Dseq();
		}

		Cinch_K.pass_V += cik_row;			/* STORE ROWS ADDED */
		if (cik_row && dev_print(CINCH,__LINE__)) {
			printf("Post cinch-k k=%d loop: tela_m=%3d (lenseq = %3d).\n", k, tela_m, lenseq);
		}

		if (k > 1)	/* NOT NEEDED AFTER k EQUALS ONE */
			clear_2D_ar(pathbox);

	} /* END OF FOR k LOOPS */ 

	mha_UPPERback(pathbox, align2D); /* THIS ALSO SAVES 2D-WIDTH in Current.pass_W */

	Cinch_K.pass_W = Current.pass_W;	/* ASSIGN CURRENT WIDTH and PASS i WIDTH HISTORY */
	free(x_history);

	return (Cinch_K.pass_V);
}
/******************************************************************/


/******************************************************************/
unsigned int nudgelize(void)
{
int cyc_col=0, cyc_row=0, i, j, m=0, n=0;
int lenseq = Clean.pass_W;
int cyc_width = Current.pass_W;						/* THIS IS opt_W SLOT TO STORE CURRENT 2-D WIDTH */
short unsigned int edge0=0;
unsigned short int dud_nudge=0;
char blnk = Fill->sym, letr;
unsigned int connudge(char con_align2D[][MAXROW], int n_start, int n_width);

	clear_pathbox();
	mha_writeback(align2D, pathbox);
	for (j=0; j<=lenseq; j++)
		pathbox[MAXROW][j] = consensus[j];
 
	/* FLAG SPECIAL CASE OF CYCLING NEED AT n=0 COLUMN */
	if (!Nudge.pass_V && pathbox[0][0] == blnk)
		edge0 = 1;

	for (n = 0; n <= cyc_width; n++) {
		for (m = 1; align2D[m][0] != '\0' && m <= lenseq; m++) {
			if (isalpha(letr=align2D[m][n])) {
				if (letr != pathbox[MAXROW][n] || edge0) {
					if (edge0) {
						while (isalpha(align2D[m][0])==0 && m<=lenseq)
							m++;
						while (align2D[0][n] == blnk)
							n++;

						cyc_row = m;	/* THIS IS ROW COORDINATE OF NON-CONSENSUS */
						cyc_col = 0;	/* THIS IS COLUMN COORDINATE OF NON-CONSENSUS */
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
						if (align2D[m][cyc_col+1] == slip.sym) {
							while (align2D[m][cyc_col+1] == slip.sym) {
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

					Nudge.pass_V = 3;	/* WILL CHAGE EVENTUALLY; IT'S LIKE THIS FOR DEPRECATED REASONS */

					/* NUDGE-CYCLELIZE: */
					if (connudge(pathbox, 0, cyc_width) == 0) {
						if (dev_print(CINCH,__LINE__)) {
							printf("dud_nudge");
						}
						dud_nudge = 1;
						i = Current.pass_V;
						Cinches[i]->pass_W = cyc_width = Current.pass_W;	/* ASSIGN CURRENT WIDTH and PASS x WIDTH HISTORY */
						n = cyc_width+1; 		/* BREAK OUT OF FOR n LOOP AFTER BREAKING OUT OF FOR m LOOP */
						break; 					/* BREAK OUT OF FOR m LOOP */
					}

					clear_right(pathbox);

					i = Current.pass_V;
					Cinches[i]->pass_W = cyc_width = Current.pass_W;	/* ASSIGN CURRENT WIDTH and PASS x WIDTH HISTORY */
					n = cyc_width+1; 		/* BREAK OUT OF FOR n LOOP AFTER BREAKING OUT OF FOR m LOOP */
					break; 					/* BREAK OUT OF FOR m LOOP */

				} /* END OF IF LETTER != CONSENSUS */
			} /* END OF IF ISALPHA */
		} /* END OF FOR m LOOP */
	} /* END OF FOR n LOOP */

	/* PUSH LEFT IF EMPTY: SHOULD MOVE TO ITS OWN FUNCTION IF NEEDED ELSEWHERE */
	if (edge0 && !dud_nudge) {
		i = 0; /* COUNTER FOR AMOUNT TO PUSH LEFT */
		for (n=0; n < cyc_width; n++) {
			for (m=0; pathbox[m][n] != '\0' && m < MAXROW; m++) {
				if (pathbox[m][n] != blnk) {
					n = MAXROW; 	/* TO BREAK OUT OF FOR n LOOP */
					break;			/* TO BREAK OUT OF FOR m LOOP */
				}
			}
			if (n < MAXROW)
				i++;
		} 
		for (m = 0; pathbox[m][0] != '\0' && m < MAXROW; m++) {
			for (n = 0; n < cyc_width+i; n++) {
				pathbox[m][n] = pathbox[m][n+i];
			}
		}
	} 

	if (!dud_nudge)	{
		mha_writeback(pathbox, align2D);
		mha_writeconsensus(pathbox, consensus);
		return(print_2Dseq());
	}
	else {
		return (0);
	}
}
/******************************************************************/


/******************************************************************/
unsigned int cinch_d(short unsigned int cinch_d_opt)
{
int delta_mrow=0, delta_ncol=0, h=0, i=0, j=0, k=WIDTH, l=0, m=0, n=0, num=0, w=0, x=0, tot_repeats=0, uniq_TRs=0, num_transits=0;
int cidwidth = Current.pass_W; 
int height = Current.pass_H;		/* height slot */
int translimit = 0;
unsigned short int nuctype=0, TR_check=0, first_write=1, mono_flag=1;		/* CHECK MONO IN ORDER TO KNOW TO SKIP IT */
unsigned short int nuctransit=0;						/* BIT FLAG FOR HANDLING NUCLEOTIDE TRANSITIONS SILENTLY (IGNORING) */
unsigned short int imperfect_TR=0;
char letr, ltr2;
char blnk = Fill->sym;
int lenseq = Clean.pass_W;

	clear_pathbox();
	nuctype = Clean.pass_V;		/* EQUALS ONE IF DNA, TWO IF RNA */
	if (nuctype == 1) {			/* IF DNA */
		nuctransit = 1;
	}

	/* START AT BIGGEST k-MER POSSIBLE AT 2x */
	for (k = Current.pass_W/2; k > 0; k--) {
		if (nuctransit) {
			if (k > opt_b.val) {
				if (opt_m.bit || opt_g.bit) {					/* opt_m OR opt_g ELECTED MAGIC MELTAGE OR GELLING */
					translimit = opt_m.val + TEMP;
				}
				else {
					translimit = allowed_transits(k);
				}
			}
			else {
				translimit = 0;
			}
		}

		for (n=0; n <= cidwidth-2*k; n++) {	
			mono_flag = 1;			/* MONOMER RUN FLAG IS SET TO 0, WHEN NO LONGER POSSIBLE (ANY n != n+1) */
	
			if (TR_check == 0) 		/* RE-SET COUNTER FOR NUM */
				num = 0;

			if (nuctransit)	{		/* RE-SET COUNTER FOR NUMBER OF TRANSITIONS */
				num_transits = imperfect_TR = 0;
			}

			for (l=0; l < k; l++) {
				letr=consensus[n  +l];
				ltr2=consensus[n+k+l];
				if (isalpha(letr=consensus[n+l]) && isalpha(ltr2=consensus[n+k+l])) {
					if (l+1<k && letr!=consensus[n+l+1])
						mono_flag = 0;		/* CAN NO LONGER BE A HOMOPOLYMER RUN OF FOR THIS k-MER */

					if (letr==ambig.sym || ltr2==ambig.sym)
						break;
					else if (num_transits > translimit)
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
						else {
							imperfect_TR = 1;
							++num_transits;
						}
					}
					else if (letr!=ltr2)
						break; 	
				}
				else if (letr==Term->sym || ltr2==Term->sym) {
					break;
				}
				else if (!cinch_d_opt && !isalpha(letr) && isalpha(consensus[n+l+1])) {
					for (i=0; align2D[i][0]!='\0'; i++) {
						for (j=n+l; (letr=align2D[i][j+1])!='\0'; j++) {
							align2D[i][j] = letr;
						}
					}
					for (j=n+l; (letr=consensus[j+1])!='\0'; j++) {
						consensus[j] = letr;
					}
					--Current.pass_W;
					break; 
				}
				else if (letr!=ltr2) {
					break; 	
				}
			} /* END OF FOR l LOOP CHECK OF k-MER TR (2x) */

			if (l==k && mono_flag==0) {
				if (nuctransit) {
					if (opt_m.bit && num_transits > translimit) { 	/* opt_m (OR opt_g) ELECTED MAGIC MELTAGE OR GEL */
						if (dev_print(CINCH,__LINE__)) {
							printf("At n=%4d: Skipping k=%d, num_transits=%d, and translimit=%d as set by opt_m/g and opt_b.val=%d.", 
										n+1, k, num_transits, translimit, opt_b.val);
						}
					}
					else if (num_transits <= translimit) {
						TR_check = 1;
						++tot_repeats;
						num = 2;		/* THIS KEEPS COUNT OF HOW MANY REPEATS. WITH ONE RE-PEAT COUNTED, THERE ARE TWO */
					}
				}
				else {
					TR_check = 1;
					++tot_repeats;
					num = 2;		/* THIS KEEPS COUNT OF HOW MANY REPEATS. WITH ONE RE-PEAT COUNTED, THERE ARE TWO */
				}
			}
			else {
				TR_check = 0;
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

				if (cinch_d_opt && !opt_d.bit) {	/* CINCH-D ENGINE IF NOT opt_d (SKIP-CINCH-D CINCHING) */
printf("\n 1. cinch_d. Am here. Linux. n=%2d\n", n);
					if (first_write) {
printf("\n 2. cinch_d. Am here. Linux. n=%2d\n", n);
						++Cinch_D.pass_R;

						m = 0;
						while (!(isalpha(align2D[m][n+k]))) {
							m++;
						}

						/* SCAN TELA STRUCT FOR COORDINATES */
						for (l=0; l< lenseq; l++) {
							if (tela[l].y==m && tela[l].x==n) {
								break;
							}
						}
						if (l!=lenseq && tela[l+k].echoes==cyc_skip.sym)
							break;

						if (dev_print(CINCH,__LINE__)) {
							printf("%4d. Working on %2d-mer consensus TR (%dx) at position %4d, row %4d.", 
									Cinch_D.pass_R, k, num, n, m);
						}

						if (imperfect_TR) {
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

						mha_writeback(align2D, pathbox); 
						pathbox[m][n+k  ] = slip.sym;

						first_write = 0;	/* TURN O-F-F NEED TO WRITE REMAINING PART OF 2-D ALIGNMENT */

						delta_ncol = k;
						delta_mrow = 1;

						/* DEAL WITH LOOSE SLIP CONNECTIONS PRODUCED BY NUDGELIZING */
						for (j = 0; j < n+k; j++) {
							if (pathbox[m][j] != blnk) {
								break;
							}
						}
						if (j == n+k) {
							if (pathbox[m-1][n] == slip.sym)
								delta_mrow = -1;
							else 
								delta_mrow = 0;
						}

						for (i = m; i < MAXROW; i++) {
							for (j = n+k; j <= Current.pass_W+1; j++) {
								letr = pathbox[i+delta_mrow][j-delta_ncol] = align2D[i][j];
								if (letr==slip.sym || letr==monoR.sym) 
									pathbox[i+delta_mrow][j-delta_ncol+1] = '\0';
								else if (letr==Term->sym) {
									for (h=0; h<n; h++)
										pathbox[i+delta_mrow][h] = blnk;
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

						if (letr == Term->sym && j-delta_ncol-1 <= cidwidth) {
							Current.pass_W = j-delta_ncol-1;
							mha_writeback(pathbox, align2D);
						}

					} /* END OF IF first_write EQUALS ONE */
				} /*********************************************************************************************/
				else if (dev_print(CINCH,__LINE__)) {		/* ELSE IF PRE-CINCH-D AND DEV_PRINT OPTION */
					if (imperfect_TR) 
						printf("%4d. i-TR: %3dx %2d-mer at consensus position %3d with %d transition(s).", uniq_TRs, num, k, n, num_transits);
					else if (nuctransit) 
						printf("%4d. p-TR: %3dx %2d-mer at consensus position %3d.", uniq_TRs, num, k, n);
					else 
						printf("%4d. TR: %3dx %2d-mer at consensus position %3d.", uniq_TRs, num, k, n);
				}
			} /* END OF WHILE TR_check */

		} /* END OF FOR n LOOP */
	} /* END OF FOR k LOOP */

	i = Current.pass_V;
	if (!cinch_d_opt) {
		if (!tot_repeats) {
			Cinches[i]->pass_W = Current.pass_W;	/* ASSIGN CURRENT WIDTH and PASS WIDTH HISTORY */
		}
	}
	else {
		Cinches[i]->pass_W = Current.pass_W;	/* ASSIGN CURRENT WIDTH and PASS WIDTH HISTORY */
		if (!delta_ncol && !opt_v.bit) {
			opt_K.bit = 1;						/* EVEN IF NOT OPTIONED, GOOD TO SHOW FOR LAST RUN */
			print_2Dseq();
			return(0);
		}
		else if (tot_repeats > 1 && opt_K.bit && !opt_v.bit) {
			cidwidth = Current.pass_W;
			opt_K.bit = 0;					/* TMP ASSIGNMENT TO PREVENT PRINTING OF CONSENSUS ROW */
			consensus_2D(0, Current.pass_W);
			opt_K.bit = 1;					/* REASSIGN SETTING */
		}
		else if (tot_repeats&& opt_v.bit) {
			cidwidth = Current.pass_W;
			print_2Dseq();
		}
		else { 
			cidwidth = Current.pass_W;
			opt_K.bit = 0;					/* TMP ASSIGNMENT TO PREVENT PRINTING OF CONSENSUS ROW */
			consensus_2D(0, Current.pass_W);
			opt_K.bit = 1;					/* REASSIGN SETTING */
		}
	}
	return(tot_repeats);
}
/******************************************************************/


/******************************************************************/
void relax_2D(void)
{
int height=0, i, j, m, n, rlx_col=0, v=0, w=0, z=0;
int width = Current.pass_W; 
char blnk = Fill->sym;
char letr;
unsigned short int nuctype = Clean.pass_V;
unsigned short int nuctransit=0;

	if (nuctype == 1)		/* IF DNA */
		nuctransit = 1;

	clear_pathbox();
	mha_writeback(align2D, pathbox); 
 
	while (align2D[height][0] != '\0') {
		height++;
	}
	Current.pass_H = height;

	for (n=0; n < width; n++) {
		v = w = z = m = 0;

		while (align2D[m][n] == '\0') {
			m++;
		}

		while (isalpha(letr=align2D[m][n]) && isalpha(align2D[m][n+1])) {
			n++;
		}

		if (align2D[m][n+1] == slip.sym) {		/* EDGE DETECTED */
			while (isalpha(align2D[m+v][n-1]) && 
				   isalpha(align2D[m+v][ n ]) &&
						   align2D[m+v][n+1] == slip.sym   ) {
				v++;
			}
			while (align2D[m+v+w][n-1] == blnk && 
				   align2D[m+v+w][ n ] == letr && 
				   align2D[m+v+w][n+1] == slip.sym     ) {
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
						pathbox[m+v - rlx_col + j][n + rlx_col] = blnk;	
					for (i = 0; i < w; i++) {
						pathbox[m+v-1 - rlx_col][n + i + rlx_col + 1] = letr;
						for (j = 0; j < height; j++) 
							pathbox[m+v - rlx_col + j][n + i + rlx_col + 1] = blnk;
					}
					rlx_col = rlx_col + w;

					/* REWRITE REST OF ARRAY TO NEW COORDINATES */
					for (i = m+v+w-1; i < height; i++) {
						for (j = n+1; (letr=align2D[i][j]) != '\0'; j++) 
							pathbox[i-rlx_col][j+rlx_col] = letr;
					}

					if (nuctransit) {
						j = Current.pass_W + rlx_col + w;
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
		else if (align2D[m][n+1] == Term->sym || align2D[m][n+1] == monoR.sym) {
			pathbox[m - rlx_col][n + rlx_col] = letr;
			pathbox[m - rlx_col][n + rlx_col + 1] = align2D[m][n+1];
			pathbox[m - rlx_col][n + rlx_col + 2] = '\0';

			if (align2D[m][n+1] == Term->sym) {
				for (i=0; i<Current.pass_W; i++)
				    pathbox[m-rlx_col+1][i] = '\0';
				Current.pass_W = n + rlx_col + 1;
			}

		}

	} /* END OF FOR n LOOP */

	mha_writeback(pathbox, align2D);

	i = Current.pass_V;
	Cinches[i]->pass_W = Current.pass_W;	/* ASSIGN CURRENT WIDTH and PASS [9] WIDTH HISTORY */
}
/******************************************************************/


/******************************************************************/
int recover_1D(char *recovered_1D) 
{
int m=0, n=0, x=0;
int 	height 		= Current.pass_H;	/* height slot 			*/
int 	lenseq 		= Clean.pass_W;		/* length slot 			*/
char	letr;

	for (m=0; m<height; ) {
		for (n=0; (letr=align2D[m][n]) != slip.sym && letr!=monoR.sym && letr!=Term->sym && n < lenseq; n++) {
			if (isalpha(letr)) {
				recovered_1D[x] = letr;
				x++;
			}
		}
		if (letr==slip.sym || letr==monoR.sym)
			m++;
		else if (letr==Term->sym) {
			recovered_1D[x] = letr;
			return(x);	/* RETURN LENGTH OF RECOVERED_1D[] AFTER POPULATING IT W/ SEQUENCE */
		}
	}

	return(0);			/* RETURN FAIL */

}
/******************************************************************/

#endif		/* !FILE_CINC_SEEN */

/*************************************************************************************************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal */
/*************************************************************************************************************/