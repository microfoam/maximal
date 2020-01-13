/******************************************************************/
/* microhom-tela.h header file, since mha_v4.11.c                 */
/* This file has MHA functions either in the form of verb_tela(), */
/* or verb_object(), where object is related to tela struct.      */
/******************************************************************/

int  assign_tela(int eL, int eM, int eN, int mode);
void assign_transit(int n, int kr_src);
int  check_tela(int eM, int eN, short unsigned int mode_dim);
void clearall_tela(int n, int span, int keep_score, int mode);
int  cyclelize_tela(int cpos, int delta, int npos);
void mark_tela(void);
void print_tela(int a, int b);
void pull_tela(int n);
int  push_tela(int n2, int n1, short unsigned int axioms);
int  settle_tiescores(int n, int span, int max_score, int iteration);
int  update_tela(void);

int assign_tela(int eL, int eM, int eN, int mode)
{
	if (!mode)
		return(0);
	else {
		int i=0, j=0, l=0, conflict_flag=0;
		int lenseq = options[1][1];
		int start=0, end=lenseq;

		start = 0;
		end   = lenseq;
	
		if (mode==2) {					/* MODE 2: ASSIGN SAME COLUMN TO IDENTICAL LETTERS */
			align2D[eM][eN] = tela[eL].c;
			tela[eL].y = eM;
			tela[eL].x = eN;
		}
		else if (mode==1) {				/* MODE 1: FLAT-LINE TELA STARTING AT POINT eL */
			if (eL > 0 && eL <= lenseq) {
				i = tela[eL-1].y;
				j = tela[eL-1].x + 1; 
			}
			else if (eL == 0) {
				i = tela[eL].y;
				j = 0;
			}
			else
				conflict_flag++;
			
			if (!conflict_flag) {		
				if (dev_print(TELA,__LINE__)) {
					printf("assign_tela() flat-lining tela starting at n=%d.", eL);
				}
				for (l = eL; l <= lenseq; l++) { 
					tela[l].y = i;
					tela[l].x = j++;
				}
			}
		}
		else if (dev_print(TELA,__LINE__)) {
			printf("Undefined mode invoked.");
			return(0);
		}
	
		if (conflict_flag)
			return(0);
		else
			return(1);
	}
}


/********************************************************************/
/* ASSIGNS AND PROPAGATES TRANSITIONS TO ALL POSITIONS AT EACH UNIT */
/* MUST SAVE k AND r AT SOURCE BEFORE CALLING THIS FUNCTION AT n.   */
void assign_transit(int n, int kr_src)
{
	if (!kr_src) {					/* kr_src = ZERO MODE */
		return;						/* THIS IS A DEVELOPMENT FEATURE: EASY TO TURN O-F-F */
	}

	int i=0, j=0;
	int k, r;

	if      (kr_src == 1) {		/* kr_src = ONE */
		k = tela[n].all_k;
		r = tela[n].all_r;
	}
	else if (kr_src == 2) {		/* kr_src = TWO */
		k = tela[n].cyc_k;
		r = tela[n].cyc_r;
	}
	else {
		k = tela[n].k;			/* kr_src = THREE (OR > TWO) */
		r = tela[n].r;
	}

	int m = n - k;

	/* ASSIGN TRANSITIONS TO .t IF IMPERFECT_TR, ALL IN REFERENCE TO FIRST UNIT STARTING AT m */
	for (i=0; i<r; i++) {
		for (j=0; j<k; j++) {	/* 1ST TIME TO NOTE THE TRANSITION POSITIONS BY DIFFERENCES */
			if (tela[n+i*k+j].c != tela[m+j].c && tela[n+i*k+j].e == tela[m+j].e) {
				tela[n+i*k+j].t = tela[m+j].t = tela[m+j].e;
			}
		}
	}
	for (i=0; i<r; i++) {
		for (j=0; j<k; j++) {	/* 2ND TIME TO PROPAGATE TRANSITIONS TO ALL PARALAGOUS POSITIONS */
			if (tela[n+i*k+j].t != tela[m+j].t) {
				tela[n+i*k+j].t = tela[m+j].t;
			}
		}
	}
	/* WRITE TRANSITIONS TO CONSENSUS ROW OF ALIGN2D */
	for (i=0; i<k; i++) {
		if (consensus[tela[m+i].x] != 'R' && consensus[tela[m+i].x] != 'Y')
			consensus[tela[m+i].x] = tela[m+i].t;
	}
}


/****** TELA: A FABRIC, UNDER AXIOMATIC LAWS **********************/
int check_tela(int eM, int eN, short unsigned int mode_dim) 
{
	if (!mode_dim)	/* check_tela called in O-F-F mode, but will return 1+2 = success */
		return(3);
	else {
		int i=0, j=0, lineM=0, lineN=0, axioms=0, badflag=0;
		int lenseq = options[1][1];
	
		if (eM>=eN) {
			axioms-=5;
			if (dev_print(TELA,__LINE__)) {
				printf("Need to call check_tela explicitly with %d-D positions eN > eM.", mode_dim);
			}
		}
	
		if (mode_dim==1) {
			/* GET 2-D COORDINATES FROM 1-D COORDINATES */
			lineM = tela[eM].x;	
			lineN = tela[eN].x;	
		}
		else if (mode_dim==2) {
			/* SAVE 2-D COORDINATES */
			lineM = eM;	
			lineN = eN;	
			/* TRANSLATE 2-D COORDINATES INTO 1-D COORDINATES */
			while (tela[i].y != lineM && i<lenseq)
				i++;
			eM = i;
	
			j=lenseq;
			while (tela[j].x != lineN && j>0)
				j--;
			eN = j;
		}
		else {
			axioms-=7;
			if (dev_print(TELA,__LINE__)) {
				printf("Need to call check_tela explicitly with dimension dim=1 or dim=2.");
			}
		}
	
		if (axioms==0) {
			/* AXIOM ONE: CONTINUITY */
			for (i=eM; i<eN; i++) {
				if (     tela[i+1].x == tela[i].x + 1 &&
					     tela[i+1].y == tela[i].y) {
					;
				}
				else if (tela[i+1].y == tela[i].y + 1 &&
						 tela[i+1].x <=  tela[i].x) {
					;
				}
				else {
					tela[i].cyc_o = '!';						/* MARK EDGE OF DISCONTINUITY */
					break;
				}
			}
			if (i==eN)
				axioms = 1;
			else if (dev_print(TELA,__LINE__)) {
				printf("check_tela(mode_dim=%d): Problem of continuity at 1-D positions %d --> %d (columns %d and %d)", 
									mode_dim, i, i+1, tela[i].x, tela[i+1].x);
			}
	
			/* AXIOM TWO: EQUIVALENCE */
			for (i=eM; i<eN && !badflag; i++) {
				for (j=i+1; j<eN; j++) {
					if (tela[j].x == tela[i].x && 
						tela[j].e != tela[i].e) {
						tela[j].cyc_o = tela[i].cyc_o = '*';	/* MARK PAIR OF NON-EQUIVALENT SITES SHARING SAME COLUMN */
						badflag++;
						break;		/* TO BREAK FOR j LOOP */
					}	
				}
				if (badflag)
					break;
			}
			if (!badflag) 
				axioms+=2;
			else if (dev_print(TELA,__LINE__)) {
				printf("check_tela(mode_dim=%d): Problem of equivalence at 1-D positions %d and %d (both in column %d)", 
									mode_dim, i, j, tela[i].x);
			}
		}
		return(axioms);	/* 0 IF BOTH FAIL; +1 IF ONLY ONE PASSES; +2 IF ONLY TWO PASSES; +3 IF BOTH PASS */ 
	}
}


/* FUNCTION TO CLEAR _ALL ELEMENTS IN TELA WITHIN A NON-CONFLICTED CYCLING ISLAND W/ A SINGLE MAX-SCORE */
void clearall_tela(int n, int span, int keep_score, int mode)
{
	int i;

	if (!mode)				/* MODE 0=O-F-F, 1=CLEAR S ONLY; 2=CLEAR ALL */
		return;

	for (i=n; i< n+span; i++) {
		if (tela[i].all_Z != keep_score) {
			tela[i].all_S = 0; 
			if (mode>1) {
				tela[i].all_k = tela[i].all_r = tela[i].all_Z = 0; 
				tela[i].all_L = tela[i].all_R = 0;
			}
		}
	}
}


/**** FUNCTION TO CYCLELIZE REPEAT AT POSITION BY DELTA *****/
int cyclelize_tela(int cpos, int delta, int npos)
{
	int    k = tela[cpos].k;
	int reps = tela[cpos].r;
	int i, j, m, n, r, z=0;
	int lenseq = options[1][1];
	char c;
	char blnk = (char) options[1][11];		/* opt_B blank character */

	if      ( cpos>lenseq ||  cpos<0 || tela[cpos].cyc_l < 2)
		return(0);
	else if (delta>lenseq || delta<0 || delta > tela[cpos].cyc_l)
		return(0);

	if (k && reps && tela[cpos].cyc_o == 'x') {
		pull_tela(cpos);
		pull_tela(npos);
		z = cpos;
		for (r=0; r<reps; r++) {
			for (j=0; j<delta; j++) {
				z++;					/* VAR z is 1D cycling position */
				c = tela[(i=cpos+r*k)].c;
				m = tela[i].y;
				n = tela[i].x;

				align2D[m][n] = blnk;
				m = m-1;
				n = n+k;
				tela[i].y = m;
				tela[i].x = n;
				align2D[m][n] = c;
			}
			if (r < reps - 1) {
				align2D[m][n+1] = '/';
				align2D[m][n+2] = '\0';
			}
			else if (r == reps-1) {
/*				n = n - delta;
*/				for (j = z+1; j < npos; j++) {
					align2D[m][++n] = tela[j].c;
					tela[j].y = m;
					tela[j].x = n;
				}
				m++;
				tela[npos].y = m = tela[(npos - tela[npos].k)].y + 1;
				tela[npos].x = n = tela[(npos - tela[npos].k)].x;
			}
		}	
		for (j=0; j<lenseq; j++)
			align2D[m][j] = '\0';
		
		for (j=npos+1; j<=lenseq; j++) {
			tela[j].x = ++n; 
			tela[j].y = m; 
		}

		tela[cpos+delta].cyc_o = 'x';
		return (1);		/* RETURN SUCCESS, BUT EVENTUALLY ADD A CHECK_TELA CALL IN HERE */
	}
	else
		return (0);
}


/**************** FUNCTION TO MARK ALL POSSIBLE k-MERs BEFORE LEGACY CINCH-T PASS ******************************/
void mark_tela(void) 
{
	int i, j, m, n, k, reps, span; 
	int threshold=0, max_score=0, max_count=0, min_k;
	int lenseq = options[1][1];
	unsigned short int nuctype = options[1][13], nuctransit=0, TRcheck=0, imperfect_TR=0, Aimperfect_TR=0, gapcheck=0;
	int homopoly_flag=0, Did=0, Dtr=0, Atr=0;
	int mismatch   = -1;	/* MOVE ME TO HEADER FILE */
	unsigned short int checkconflict=0;

	if (nuctype == 1)		/* IF DNA */
		nuctransit = 1;

	for (n = 1; n<=lenseq; n++) {
		for (m = 0; m < n; m++) {
			/* FOR ROW m LOOP 1/5: SLIDE DOWN TO ROW WITHIN POPULATED HEMIDIAGONAL */
			if (n-m > WIDTH+1) 
				m = n-WIDTH;

			/* FOR ROW m LOOP 2/5: SET K-MER SIZE AND DTHR SCORE THRESHOLD */
			k = n-m;
			if (nuctransit) {
				threshold = score_DTHR(k);
			}

			/* FOR ROW m LOOP 3/5: SKIP k=ONE */
			if (k == 1) {
				break;	/* GO TO NEXT n */
			}

			/* FOR ROW m LOOP 4/5: SET HOMOPOLYMER RUN STATUS UNKNOWN; USED TO RULE OUT k>1 MONONUCLEOTIDE "REPEATS" */
			homopoly_flag = 2;
			if (tela[n].c != tela[n-1].c)
				homopoly_flag = 0;

			/* FOR ROW m LOOP 5/5: START COUNTING SCORE IF PATHBOX POSITION HAS VALUE > MISMATCH */
			if (pathbox[m][n] > mismatch && n+k <= lenseq) {
				Dtr = 0;

				/* IF SUMMING PATHBOX DIAGONAL 1/4: COMPUTE SCORES OF IDENTITY LINE AND REPEAT DIAGONAL*/
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

				/* IF SUMMING PATHBOX DIAGONAL 2/4: SET HOMOPOLYMERIC RUN BIT TO TRUE IF DETECTED 	*/
				if (homopoly_flag && i == n) {
					homopoly_flag = 1;				/* BIT IS THERE IF NEEDED BEYOND BREAK. 		*/
					Dtr = 0;
					break;							/* GO TO NEXT n */
				}

				/* IF SUMMING PATHBOX DIAGONAL 3/4: IF CONSIDERING NUCL. TRANSITIONS AS PARTIAL MATCHES */
				if (nuctransit && Dtr && Dtr!=Did) { 
					if (k>PISO && 100*Dtr/Did > threshold)	{	
						imperfect_TR = 1;		/* CALLING TR W/ TRANSITIONS FOR n BLOCK VS m BLOCK */
					}
					else 
						Dtr = 0;
				} 

				/* IF SUMMING PATHBOX DIAGONAL 4/4: START COUNTING REPEATS */
				if (Dtr && (Dtr==Did || imperfect_TR)) {
					/* COUNT NUMBER OF REPEATS ALBERT-STYLE */
					TRcheck = 1;
					reps = 1;
					tela[n].all_k = k;
					tela[n].all_r = reps;
					tela[n].all_S = Dtr;	/* SAVE INITIAL UNIT SCORE */
					while (TRcheck) {
						Atr = 0;
						if (m + (reps+1)*k >= lenseq) { 
							Atr = 0;
							break;
						}
						else {
							for (i = m; i < n  ; i++) {		/* COMPARE TO FIRST UNIT */
								if ( (j=pathbox[i][(i + (reps+1)*k)]) == mismatch) {
									Atr = 0;
									TRcheck = 0;
									break;
								}
								else
									Atr = Atr + j;
							}
						}
						if (nuctransit) { 
							if (Atr!=Did && (100*Atr)/Did > threshold)
								Aimperfect_TR = 1;
							else
								Aimperfect_TR = 0;
						} 
						if (Atr==Did || Aimperfect_TR) {
							reps++;
							tela[n].all_S += Atr;
							Atr = 0;
						}
						else {		/* ELSE FINAL NUMBER OF REPEATS (REPS) IS NOW KNOWN *****************/
							tela[n].all_r = reps;
							break;
						}
					}
					break;		/* OTHERWISE MAY OVERWRITE TR WITH ONE OF SMALLER K */
				}
			}
		} /* END OF FOR m */
	} /* END OF FOR n */

	/* FILL IN CYCLING GAPS CAUSED BY BELOW THRESHOLD FRAMES: THIS SUPPRESSES INTRA-TR CONFLICT REPORTING */
	for (n=0; n<=lenseq; n++) {
		if (tela[n].all_k && !tela[n+1].all_k) {
			k = tela[n].all_k;
			gapcheck = 0;
			for (i=n+2; i <= n + k * tela[n].all_r; i++) {
				if (tela[i].all_k && tela[i].all_k != k) {
					gapcheck = 0;
					break;
				}
				else if (tela[i].all_k == k) {
					gapcheck = 1;
					break;
				}
			}
			if (gapcheck) {
				for (j=i-1; j>n; j--) {
					tela[j].all_k = k;
					tela[j].all_r = 0;
				}
			}
		}
	}

	/* NOW MARK ALL CONFLICTING TRs */
	for (n=lenseq; n>0; n--) {
		if (tela[n].all_k) {
			k = tela[n].all_k;
			m = n - k;
			j = n - 1;
			while (tela[j].all_k)		/* SKIP CYCLE COLUMNS OF SAME K-MER */
				j--;
			for (i=j; i>0; i--) {
				if (tela[i].all_k && (i + tela[i].all_k * (tela[i].all_r-1)) > m) {
					tela[n].all_L = i;		/* UPDATE LEFT-MOST OVERLAPPING & CONFLICTING TR */
					tela[i].all_R = n;		/* UPDATE RIGHT-MOST OVERLAPPING & CONFLICTING TR */
				}
			}
		}
	}

	/* IDENTIFY CYCLING ISLANDS WITHOUT CONFLICT AND DETERMINE TIE-BREAKER SCENARIOS */
	/* IDENTIFY CYCLING ISLANDS WITH SOME TYPES OF CONFLICT AND RESOLVE PICK */
	for (n=0; n<=lenseq; n++) {
		if (tela[n].all_S) {
			checkconflict = span = 1;
			max_count = 0;
			min_k = k = tela[n].all_k;
			max_score = tela[n].all_S;
			if (tela[n].all_L || tela[n].all_R) {
				checkconflict = 0;
			}
			for (i=n+1; tela[i].all_k; i++) {
				span++;
				if (checkconflict) {
					if (tela[i].all_L || tela[i].all_R) {
						checkconflict = 0;
					}
					else {
						if (tela[i].all_S > max_score) 
							max_score = tela[i].all_S;
						if (tela[i].all_k < min_k) 
							min_k= tela[i].all_k;
					}
				}
			}
			if (!checkconflict) {
				/* CONFLICT SCENARIO ONE */
				if (span==1 && !(tela[n].all_L) && tela[n].all_R) {
					j = tela[n].all_R;
					int k2 = tela[j].all_k;
					for (i=j+1; tela[i].all_k == k2; i++) {
						if (!tela[i].all_L && !tela[i].all_R && tela[i].all_S > tela[j].all_S) {
							clearall_tela(j, i-j, tela[i].all_S, TWO);		/* O-F-F, ONE, OR TWO */
							tela[i].all_Z = tela[i].all_S;
							tela[n].all_Z = tela[n].all_S;
							tela[n].all_R = 0;
						}
					}
				}
				/* CONFLICT SCENARIO TWO */
				else if (span>1 && tela[n].all_L && !(tela[n].all_R) && tela[n+1].all_k < tela[n].all_k &&
							!(tela[n+1].all_L) && !(tela[n+1].all_R) && tela[n].all_k % tela[n+1].all_k==0) {
					clearall_tela(n, 2, tela[n+1].all_S, TWO);		/* O-F-F, ONE, OR TWO */
					/* POSSIBLE THIS CASE COULD BE GENERALIZED...FOR A RAINY DAY */
				}
			}
			else if (span==1) {
				tela[n].all_Z = tela[n].all_S;
				if (tela[n].all_S != tela[n].r*k*MATCH)
					assign_transit(n,OFF);	/* O-F-F; ONE=ALL_K/R; TWO=CYC_K/R; THREE=K/R */
			}
			else {
				for (i=n; i<n+span; i++) {
					if (tela[i].all_S == max_score && tela[i].all_k == min_k) {
						max_count++;
						tela[i].all_Z = tela[i].all_S;
					}
				}
				if (max_count > 1) {
					j=1;	/* VAR j WILL COUNT ITERATIONS REQUIRED TO BREAK TIES */
					while (max_count > 1)
						max_count = settle_tiescores(n, span, max_score, j++);
				}
				else
					clearall_tela(n, span, max_score, ONE);			/* O-F-F, ONE, OR TWO */
			}
			n = n + span - 1;
		}
	}
}

/**** FUNCTION TO CONSIDER GHOST FLANKING REPEAT UNITS AS TIE-BREAKERS */
/**** -GHOST --> SO-NAMED B/C THESE ARE VIRTUAL SUB-THRESHOLD TR UNITS */
/**** RETURNS THE NUMBER OF TIES (MAX_COUNT)                           */
int settle_tiescores(int n, int span, int max_score, int iteration)
{
int i=0, j=0, k=0, up, dn, m, o;
int match = MATCH;
int transition = TRANSITION;
int lenseq = options[1][1];
int max_count=0, ratchet=0;

	for (i=n; i<n+span; i++) {
		if (tela[i].all_S == max_score) {
			k = tela[i].all_k;
			m = i-k;					/* DEFINES START OF FIRST REPEAT UNIT */
			o = i+k*(tela[i].all_r-1);	/* DEFINES START OF LAST REPEAT UNIT */
			up = m - k*iteration;		/* DEFINES THE GHOST FLANKING UNIT STARTING AT m-1 */
			dn = o + k*iteration;		/* DEFINES THE GHOST FLANKING UNIT STARTING AFTER REPEATS */
			for (j=0; j<k; j++) {
				if (up+j >= 0 && tela[up+j].e == tela[m+j].e) {
					if (tela[up+j].c == tela[m+j].c) {
						tela[i].all_Z += match;
						ratchet++;
					}
					else {
						tela[i].all_Z += transition;
						ratchet++;
					}
				}
				if (dn+j<=lenseq && tela[dn+j].e == tela[o+j].e) {
					if (tela[dn+j].c == tela[o+j].c) {
						tela[i].all_Z += match;
						ratchet++;
					}
					else {
						tela[i].all_Z += transition;
						ratchet++;
					}
				}
			}
		}
	}
	if (!ratchet)		/* IF NO RATCHETING, THEN GHOST ITERATIONS ARE OUT-OF-BOUNDS AND UNEFFECTIVE TIE-BREAKERS */
		return(0);
	max_score = tela[n].all_Z;
	for (i=n+1; i<n+span; i++) {
		if (tela[i].all_Z > max_score)
			max_score = tela[i].all_Z;
	}
	for (i=n; i< n+span; i++) {
		if (tela[i].all_Z == max_score)
			max_count++;
	}
	if (max_count == 1) 				/* THERE IS ONE OPTIMAL CINCH LOCATION AND NO CONFLICT SO ERASE ALL OTHERS */
		clearall_tela(n, span, max_score, TWO);		/* O-F-F, ONE, OR TWO */

	return(max_count);
}

/********************* PRINT TELA FROM a TO b *********************/
void print_tela(int a, int b)
{
int i=0, f=0;
int width = 50;	    /* 60 x 3 = 180 COLS, PRACTICAL DISPLAY WIDTH */
int lenseq = options[1][1];

	if (a<0 || a>lenseq) 
		a = 0;
	if (b>lenseq)
		b = lenseq;

	if (width > lenseq) {
		a = 0;
		b = lenseq;
	}
	else if (b>a && b-a < width) {
		if (a + width < lenseq)
			b = a + width;
		else
			b = lenseq;
	}
	else if (b <= 0 && a >=0) {
		if (a + width > lenseq)
			b = lenseq;
		else
			b = width;
	}
	/************* BEGIN PRINTING LINES *******************/
	printf("\n e:");
	for (i=a; i<=b; i++) {
		printf("%3c", tela[i].e);
	}
	printf("\n t:");
	for (i=a; i<=b; i++) {
		if (tela[i].c != tela[i].t)
			printf("%3c", tela[i].t);
		else
			printf("  .");
	}
	printf("\n c:");
	for (i=a; i<=b; i++)
		printf("%3c", tela[i].c);

	printf("\n n:");
	for (i=a; i<=b; i++)
		printf("%3d", i);

	printf("\nLf:");
	for (i=a; i<=b; i++) {
		if (tela[i].cyc_Lf)
			printf("%3d", tela[i].cyc_Lf);
		else
			printf("  <");
	}

	printf("\nRt:");
	for (i=a; i<=b; i++) {
		if (tela[i].cyc_Rt)
			printf("%3d", tela[i].cyc_Rt);
		else
			printf("  >");
	}
	printf("\n X:");
	for (i=a; i<=b; i++) {
		if (tela[i].X != i)
			printf("%3d", tela[i].X);
		else
			printf("  .");
	}
	printf("\n y:");
	for (i=a; i<=b; i++)
		printf("%3d", tela[i].y);

	printf("\n x:");
	for (i=a; i<=b; i++)
		printf("%3d", tela[i].x);

	printf("\nak:");
	for (i=a; i<=b; i++) {
		if (tela[i].all_k)
			printf("%3d", tela[i].all_k);
		else
			printf("  .");
	}
	printf("\nar:");
	for (i=a; i<=b; i++) {
		if (tela[i].all_k)					/* THIS IS MEANT TO BE ALL_K */
			printf("%3d", tela[i].all_r);
		else
			printf("  .");
	}
	printf("\naS:");
	for (i=a; i<=b; i++) {
		if (tela[i].all_S)
			printf("%3d", tela[i].all_S);
		else
			printf(" __");
	}
	printf("\naZ:");
	for (i=a; i<=b; i++) {
		if (tela[i].all_Z)
			printf("%3d", tela[i].all_Z);
		else
			printf(" __");
	}
	printf("\naR:");
	for (i=a; i<=b; i++) {
		if (tela[i].all_R)
			printf("%3d", tela[i].all_R);
		else
			printf("  .");
	}
	printf("\naL:");
	for (i=a; i<=b; i++) {
		if (tela[i].all_L)
			printf("%3d", tela[i].all_L);
		else
			printf(" __");
	}

	printf("\n k:");
	for (i=a; i<=b; i++) {
		if (tela[i].k)
			printf("%3d", tela[i].k);
		else
			printf("  .");
	}
	printf("\n r:");
	for (i=a; i<=b; i++) {
		if (tela[i].r)
			printf("%3d", tela[i].r);
		else
			printf(" __");
	}

	printf("\nDt:");
	for (i=a; i<=b; i++) {
		if (tela[i].Dtr)
			printf("%3d", tela[i].Dtr);
		else
			printf("  .");
	}

	printf("\n o:");
	for (i=a; i<=b; i++) {
		if (tela[i].o)
			printf("%3d", tela[i].o);
		else
			printf("  .");
	}
	printf("\n E:");
	for (i=a; i<=b; i++)
		printf("%3c", tela[i].echoes);

	printf("\ncO:");
	for (i=a; i<=b; i++)
		printf("%3c", tela[i].cyc_o);

	printf("\ncL:");
	for (i=a; i<=b; i++) {
		if (tela[i].cyc_l) 
			printf("%3d", tela[i].cyc_l);
		else
			printf("  .");
	}

	printf("\ncK:");
	for (i=a; i<=b; i++) {
		if (tela[i].cyc_k)
			printf("%3d", tela[i].cyc_k);
		else
			printf("  .");
	}

	printf("\ncR:");
	for (i=a; i<=b; i++) {
		if (tela[i].cyc_r)
			printf("%3d", tela[i].cyc_r);
		else
			printf("  .");
	}

	printf("\ncP:");
	for (i=a; i<=b; i++) {
		if (tela[i].cyc_P)
			printf("%3d", tela[i].cyc_P);
		else
			printf("  .");
	}

	printf("\ncS:");
	for (i=a; i<=b; i++) {
		if (tela[i].cyc_S)
			printf("%3d", tela[i].cyc_S);
		else
			printf("  .");
	}

	printf("\n   ");
	for (i=a; i<=b; i++)
		printf(" --");

	/* PRINT TOP FRAME ROWS */
	for (f=1; f<9; f++) {
		printf("\nf%d:", f);
		for (i=a; i<=b; i++) {
			if (tela[i].cyc_F[f])
				printf("%3d", tela[i].cyc_F[f]);
			else
				printf("  .");
		}
	}

	printf("\n");
}


/* ERASE REPEAT AT n, INCLUDING ANY TRANSITIONS, LEAVES k MARK ****/
void pull_tela(int n)
{
	int i=0;
	int k = tela[n].k;
	int r = tela[n].r;
	int m = n-k;

	for (i=m; i<(n+k*r); i++) {
		consensus[(tela[i].x)] = tela[i].t = tela[i].c;
	}

/* THIS CAUSED ME LOTS OF PAIN AND I LEAVE IT HERE TO REMIND ME. -AJE */
/*	for (i=0; i<k*r; i++) {				
		tela[n+i].y = tela[m    ].y;
		tela[n+i].x = tela[m+k+i].x;
	}
*/
	tela[n].r = 0;
	tela[n].echoes = tela[n].cyc_o = 'o';
}


/**** FUNCTION TO TOKENIZE A 2-D ROW IF THERE ARE NO VIOLATIONS BASED ON TESTS OF AXIOMS 1 OR 2, NONE (0), OR BOTH (3) *****/
int push_tela(int n2, int n1, short unsigned int axioms) 
{
	int coord_A_x=0, coord_B_x=0;
	int coord_A_y=0, coord_B_y=0;
	int i = 0;
	int k = n2 - n1;
	int lenseq=options[1][1];
	int violation=0;

	/* CHECK VALIDITY OF THE INPUT */
	if (k<=0 || lenseq==0 || axioms<0 || axioms>3)
		return(10);

	if (axioms==1 || axioms==3) {
		/* CHECK PRINCIPLE OF CONTINUITY */
		for (i=0; i<k; i++) {
			coord_A_x = tela[n1+i    ].x;
			coord_A_y = tela[n1+i    ].y;

			coord_B_x = tela[n1+i + 1].x;
			coord_B_y = tela[n1+i + 1].y;

			if      (coord_B_y == coord_A_y && coord_B_x == coord_A_x + 1) 
				;
			else if (coord_B_y == coord_A_y + 1 && coord_B_x <= coord_A_x) 
				;	
			else {
				violation = 1;
				if (dev_print(TELA,__LINE__)) {
					printf("push_tela() viol-1 for n2=%d k=%d (i=%d): A_x=%d, A_y=%d, B_x=%d, B_y=%d.", 
												n2, k, i, coord_A_x, coord_A_y, coord_B_x, coord_B_y); 
				}
				break;
			}
		}
	}

	if (axioms==2 || axioms==3) {
		/* CHECK PRINCIPLE OF EQUIVALENCE */
		for (i=0; i<k; i++) {
			if      (tela[n1+i].c == tela[n2+i].c)
				;
			else if (tela[n1+i].t == tela[n2+i].t) 
				;
			else if (OFF && tela[n1+i].e == tela[n2+i].e) { 
				;
			}
			else {
				violation += 2;
				if (dev_print(TELA,__LINE__)) {
					printf("push_tela() viol-2 at %d and %d for k=%d.", n1+i, n2+i, k);
				}
				break;	
			}
		}
	}
	else if (!violation) {		/* ELSE ASSIGN EQUIVALENCE HERE PROVIDED AXIOM OF CONTINUITY CHECKED OUT (1) OR WAS NOT CHECKED (0) */
		for (i=0; i<k; i++) {
			if (tela[n1+i].c != tela[n2+i].c &&			/* ASSIGN SAME TRANSITS */
				tela[n1+i].e == tela[n2+i].e) {
					tela[n1+i].t = tela[n1+i].e;
					tela[n2+i].t = tela[n2+i].e;
			}
			else if (tela[n1+i].c == tela[n2+i].c &&	/* PROPAGATE TRANSIT COLUMNS THROUGH THE UNITS */
				tela[n1+i].t != tela[n2+i].t) {
					tela[n1+i].t = tela[n1+i].e;
					tela[n2+i].t = tela[n2+i].e;
			}
		}
	}

	/* UPDATE TELA STRUCTURE IF AND ONLY IF THERE ARE NO VIOLATIONS */
	if (!violation) {
		for (i=0; i<k; i++) { 
			tela[n2 + i].x = tela[n1+i].x;
			tela[n2 + i].y++;
		}
		for (i=n2+k; i<=lenseq; i++) {
			tela[i].x -= k;
			tela[i].y++;
		}
	}
	return(violation);	/* RETURNS 1 IF CONTINUITY FAILS, 3 IF EQUIVALENCE FAILS, 4 IF BOTH FAIL */
}


/* CALL ONLY IF align2D is 100% ***********************************/
int update_tela(void)
{
	int c=0, i=0, j=0;
	int lenseq=options[1][1];
	char letr;

	for (i=0; align2D[i][0]!='\0'; i++) {
		for (j=0; align2D[i][j]!='\0'; j++) {
			while(!isalpha(align2D[i][j]))
				j++;
			while(isalpha(letr=align2D[i][j])) {
				if (letr!=tela[c].c && letr!=tolower(tela[c].c)) {
					return(c);		/* 1-D COORDINATE OF DISCREPANCY */
				}
				else {
					c++; 
					j++;
				}
			}
			break;
		}
	}

	if (c==lenseq) {
		for (i=0; align2D[i][0]!='\0'; i++) {
			for (j=0; align2D[i][j]!='\0'; j++) {
				while(!isalpha(align2D[i][j]))
					j++;
				while(isalpha(letr=align2D[i][j])) {
					if (letr==tela[c].c) {
						tela[c].y = i;
						tela[c].x = j;
						c++; 
						j++;
					}
					else
						break;
				}
				if (letr=='>') {
					tela[c].y = i;
					tela[c].x = j;
				}
				else
					break;
			}
		}
	}
	return(c);	/* SHOULD BE lenseq IF SUCCESSFUL */
}


