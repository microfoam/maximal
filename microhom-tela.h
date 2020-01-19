/******************************************************************/
/* microhom-tela.h header file, since mha_v4.11.c                 */
/* This file has MHA functions either in the form of verb_tela(), */
/* or verb_object(), where object is related to tela struct.      */
/******************************************************************/

#ifndef FILE_TELA_SEEN
#define FILE_TELA_SEEN

int  		assign_tela(int pos, int eM, int eN, int mode);
void 		assign_transit(int n, int kr_src);
short int 	check_solo(int pos);
int  		check_tela(int eM, int eN, short unsigned int mode_dim);
void 		clearall_tela(int n, int span, int keep_score, int mode);
int  		cyclelize_tela(int cpos, int delta, int npos);
void 		flatline_after_TR(int pos);
void 		mark_tela(void);
void 		print_tela(int a, int b);
void 		pull_tela(int n);
int  		push_tela(int n2, int n1, short unsigned int axioms);
int 		score_kmer(int n, int k, short unsigned int mode);
int  		settle_tiescores(int n, int span, int max_score, int iteration);
int			update_tela(void);

int assign_tela(int pos, int eM, int eN, int mode)
{
	if (!mode)
		return(0);
	else {
	
		/* MODE >=1: ASSIGN COORDINATES */
		align2D[eM][eN] = tela[pos].c;
		tela[pos].y = eM;
		tela[pos].x = eN;
	
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
			if (dev_print(TELA,__LINE__)) {
				printf("Need to call check_tela explicitly with %d-D positions eN > eM.", mode_dim);
			}
			exit(1);
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
			if (dev_print(TELA,__LINE__)) {
				printf("Need to call check_tela() explicitly with dimension dim=1 or dim=2.\n");
			}
			exit(1);
		}
	
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
				tela[i].cyc_o = '>';						/* MARK EDGE OF DISCONTINUITY */
				if (dev_count < dev_limit && dev_print(TELA,__LINE__)) {
					printf("check_tela() marking edge of discontinuity at i=%d.", i);
					dev_count++;
				}
				break;
			}
		}
		if (i==eN)
			axioms = 1;
		else if (dev_count < dev_limit && dev_print(TELA,__LINE__)) {
			printf("check_tela(mode_dim=%d): Problem of continuity at 1-D positions %d --> %d (columns %d and %d)", 
								mode_dim, i, i+1, tela[i].x, tela[i+1].x);
			dev_count++;
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
		else if (dev_count < dev_limit && dev_print(TELA,__LINE__)) {
			printf("check_tela(mode_dim=%d): Problem of equivalence at 1-D positions %d and %d (both in column %d)", 
								mode_dim, i, j, tela[i].x);
			dev_count++;
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

		flatline_after_TR(npos);	

		tela[cpos+delta].cyc_o = 'x';
		return (1);		/* RETURN SUCCESS, BUT EVENTUALLY ADD A CHECK_TELA CALL IN HERE */
	}
	else
		return (0);
}

/******************************/
void flatline_after_TR(int pos)
{
	int lenseq = options[1][1];
	int i, start, x, y;

	if (!tela[pos].all_k) {
		for (i=pos-1; i>0; i--) {
			if (tela[i].all_k) {
				pos = i;
				break;
			}
		}
		if (i==0) {
			return;
		}
	}

	x = tela[pos].x + tela[pos].all_k;						/* START OF FIRST COLUMN AFTER REPEAT */
	y = tela[pos].y + tela[pos].all_r - 1;					/* FINISHING ROW */

	start = pos + tela[pos].all_k * tela[pos].all_r;		/* 1D-POINT TO START */

	for (i=start; i<=lenseq; i++) {
		tela[i].x = x++; 
		tela[i].y = y; 
	}
	if (dev_print(TELA,__LINE__)) {
		printf("Function cyclelize_tela() finished by flat-lining tela at position=%d with coordinates (x,y) = (%d,%d).", i,x,y);
	}
}


/**************** FUNCTION TO MARK ALL POSSIBLE k-MERs BEFORE LEGACY CINCH-T PASS ******************************/
void mark_tela(void) 
{
	int i, j, l, m, n, k, reps, span, min_k; 
	int threshold=0, max_score=0, max_count=0, projection=0;
	int lenseq = options[1][1];
	unsigned short int nuctype = options[1][13], nuctransit=0, TRcheck=0, imperfect_TR=0, Aimperfect_TR=0, gapcheck=0;
	int homopoly_flag=0, Did=0, Dtr=0, Atr=0;
	int mismatch   = -1;	/* MOVE ME TO HEADER FILE */
	unsigned short int checkconflict=0;

	if (nuctype == 1)		/* IF DNA */
		nuctransit = 1;

	if (dev_print(TELA,__LINE__)) {
		printf("Function mark_tela() engaged.");
	}

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
					if (dev_print(TELA,__LINE__)) {
						printf("         mark_tela() counting repeats at n=%d for k-mer=%d.", n,k);
					}
					/* COUNT NUMBER OF REPEATS ALBERT-STYLE */
					TRcheck = 1;
					tela[n].all_r = reps = 1;
					tela[n].all_k = k;
					tela[n].all_S = Dtr;	/* SAVE INITIAL UNIT SCORE */
					while (TRcheck) {
						Atr = 0;
						if (m + (reps+1)*k >= lenseq) { 
							Atr = 0;
							break;
						}

						if (nuctransit) {
							Atr = score_kmer(n+k*reps,k,TWO);
							if (Atr!=Did && (100*Atr)/Did > threshold) {
								Aimperfect_TR = 1;
							}
							else
								Aimperfect_TR = 0;
						}
						else
							Atr = score_kmer(n+k*reps,k,ONE);

						if (Atr==Did || Aimperfect_TR) {
							reps++;
							tela[n].all_S += Atr;
						}
						else {		/* ELSE FINAL NUMBER OF REPEATS (REPS) IS NOW KNOWN *****************/
							tela[n].all_r = reps;
							if (n+k*reps > projection) {
								projection = n + k*reps;
							}
							else {
								if (tela[n-1].all_k == k) {
									tela[n].stat = 'c';		/* c FOR TRIVIAL-CASE OF CYCLING FRAME TYPE REPEAT */
								}
								else {
									tela[n].stat = 'f';		/* FRACTAL REPEAT = EMBEDDED IN ANOTHER REPEAT */
								}
							}
							TRcheck = 0;
							break;
						}
					}
					if (dev_print(TELA,__LINE__)) {
						printf("                     repeats=%d at n=%d for k-mer=%d.", reps,n,k);
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
				if (dev_print(TELA,__LINE__)) {
					printf("         mark_tela() filling in gap between %d and %d, inclusive of these points, for k-mer=%d.", i-1,n-1,k);
				}
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
			while (tela[j].all_k) {		/* SKIP CYCLE COLUMNS OF SAME K-MER */
				j--;
			}
			for (i=j; i>0; i--) {
				int recslips = 0;	/* COUNTS RECENT FRACTAL SLIPS IN UPSTREAM TR SHADOW */
				if (tela[i].all_k && (i + tela[i].all_k * (tela[i].all_r-1)) > m-recslips) {
					/* CASE OF NON-CONFLICTING FRACTAL REPEATS */
					if (i>=m && i<n && tela[i].all_S == tela[i+k].all_S) {
						tela[i].stat = tela[i+k].stat = tela[n].stat = 'f';
						tela[n].all_L = i;		/* UPDATE LEFT-MOST OVERLAPPING & CONFLICTING TR */
						tela[i].all_R = n;		/* UPDATE RIGHT-MOST OVERLAPPING & CONFLICTING TR */
						recslips += span_allrk(i); 
					}
					else {
						tela[n].all_L = i;		/* UPDATE LEFT-MOST OVERLAPPING & CONFLICTING TR */
						tela[i].all_R = n;		/* UPDATE RIGHT-MOST OVERLAPPING & CONFLICTING TR */
						if (dev_print(TELA,__LINE__)) {
							printf("         mark_tela() marking conflict between i=%d and n=%d.", i,n);
						}
					}
				}
			}
		}
	}
	/* NOW MARK ALL CONFLICTING TRs AT m's W.R.T. n's */
	for (n=lenseq; n>0; n--) {
		if (tela[n].all_k && tela[(m=(n-(k=tela[n].all_k)))].all_k) {
			short unsigned int check_shadow = 2;
			int recslips = 0;
			j = n - 1;
			while (tela[j].all_k) {		/* SKIP CYCLE COLUMNS OF SAME K-MER */
				j--;
			}
			for (i=j; i>=m; i--) {
				if (tela[i].all_k) {
					check_shadow = 1;
					recslips += span_allrk(i);
				}
			}
			if (recslips == 0) {
				check_shadow = 0;
			}
			else if ( m + tela[m].all_k * tela[m].all_r >= m - recslips && tela[m].stat != 'f' && 
						tela[n-tela[m].all_k].all_k != tela[n].all_k && check_solo(m) && check_solo(n) ) {
				tela[m].stat = '_';
				clearall_tela(m, 1, -1, TWO);		/* O-F-F, ONE, OR TWO */
			}
		}	
	}

	/* IDENTIFY CYCLING MARKS FOR PERFECT REPEATS LAID OVER INTERNAL REPEATS WITH r>1 */
	/* CAN CLEAN UP HOW THESE LOOK AND PERHAPS ACT ON THEM PRIOR TO CINCH-K */
	for (n = 1; n<=lenseq; n++) {
		if (tela[n].all_r > 1 && (!nuctransit || tela[n].all_S == tela[n].all_k * tela[n].all_r)) {
			k    = tela[n].all_k;
			reps = tela[n].all_r;
			m = n-k;

			/* CHECK FIRST AND LAST REPEATS FOR CONFLICTS */
			for (j=0; j<k; j++) {
				if (tela[m+j].all_L || tela[m+reps*k+j].all_R)
					break;
				else if (tela[m+j].all_R && tela[m+j].all_R != n)
					break;
			}
			if (j==k) {		/* IF THERE IS NO CONFLICT ON EDGE REPEATS, THEN CONTINUE */
				for (j=1; j<k; j++) {
					tela[n+j].all_k = tela[m+j].all_k;
					tela[n+j].all_r = tela[m+j].all_r;
					tela[n+j].all_S = tela[m+j].all_S;
					tela[n+j].all_Z = tela[m+j].all_Z;
					tela[n+j].all_L = tela[m+j].all_L;
				}
				for (j=0; j<k; j++) {
					for (i=1; i<reps; i++) {
						tela[n+i*k+j].all_k = tela[m+j].all_k;
						tela[n+i*k+j].all_r = tela[m+j].all_r;
						tela[n+i*k+j].all_S = tela[m+j].all_S;
						tela[n+i*k+j].all_Z = tela[m+j].all_Z;
						tela[n+i*k+j].all_L = tela[m+j].all_L;
					}
				}
				n = n + k*reps - 1; 	/* NEED TO ADVANCE BEYOND CONFLICT-FREE LAST UNIT */
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
					if (OFF && tela[n].stat == '!') {	/* SO MARKED IN THE CONFLICT TR LOOP */
						clearall_tela(n, 1, -1, TWO);		/* O-F-F, ONE, OR TWO */
						if (dev_print(TELA,__LINE__)) {
							printf("mark_tela() calling clearall_tela at n=%d.", n);
						}
					}
					else {
						j = tela[n].all_R;
						int k2 = tela[j].all_k;
						if (dev_print(TELA,__LINE__)) {
							printf("mark_tela() at n=%d, span=%d with no left-conflict, but right_conflict=%d.", n, span, tela[n].all_R);
						}
						for (i=j+1; tela[i].all_k == k2; i++) {
							if (!tela[i].all_L && !tela[i].all_R && tela[i].all_S > tela[j].all_S) {
								clearall_tela(j, i-j, tela[i].all_S, TWO);		/* O-F-F, ONE, OR TWO */
								tela[i].all_Z = tela[i].all_S;
								tela[n].all_Z = tela[n].all_S;
								tela[n].all_R = 0;
							}
						}
					}
				}
				/* CONFLICT SCENARIO TWO */
				else if (span>1 && tela[n].all_L && !(tela[n].all_R) && tela[n+1].all_k < tela[n].all_k &&
							!(tela[n+1].all_L) && !(tela[n+1].all_R) && tela[n].all_k % tela[n+1].all_k==0) {
					clearall_tela(n, 1, tela[n+1].all_S, TWO);		/* O-F-F, ONE, OR TWO */
					/* POSSIBLE THIS CASE COULD BE GENERALIZED...FOR A RAINY DAY */
					if (dev_print(TELA,__LINE__)) {
						printf("mark_tela() at n=%d, span=%d with left-conflict=%d, and no right_conflict, " 
								"and n+1 has smaller k that is a multiple of k-size at n with no conflicts.", n, span, tela[n].all_L);
					}
				}
				/* CONFLICT SCENARIO THREE: FIRST TR CAN BE CYCLED BUT NOT THE SECOND, WHICH HAS HIGHER SCORE ANYWAYS */
				else if (span>1 && tela[(l=tela[n].all_R)].all_S > tela[n].all_S && !tela[l+1].all_S) {
					k = tela[n].all_k;
					for (j=n+1; j<n+span; j++) {
						if (tela[j].all_k == k && !tela[j].all_R && !tela[j].all_L) {
							if (dev_print(TELA,__LINE__)) {
								printf("mark_tela() at n=%d, scenario three", n);
							}
							for (i=j-1; i>=n; i--) {
								clearall_tela(i, j-n, -1, TWO);		/* O-F-F, ONE, OR TWO */
							}
							for (i=j+1; i<n+span; i++) {
								if (tela[i].all_k == k && tela[i].all_R == l) {
									clearall_tela(i,1,-1, TWO);		/* O-F-F, ONE, OR TWO */
								}
							}
						}
					}
				}
			}
			else if (span==1) {
				tela[n].all_Z = tela[n].all_S;	
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
				else {
					clearall_tela(n, span, max_score, ONE);			/* O-F-F, ONE, OR TWO */
					if (dev_print(TELA,__LINE__)) {
						printf("         mark_tela() engaging clearall_tela(ONE) at n=%d.", n);
					}
				}
			}
			n = n + span - 1;
		}
	}

	/* CHECK TO SEE IF THERE ARE CYCLING DIFFERENCES BETWEEN INTERNAL REPEATS */
	for (n=2; n+1<lenseq; n++) {
		if (!tela[n-1].all_k && tela[n].all_k && !tela[n+1].all_k) {
			int p,q;
			k = tela[n].all_k;
			m = n - k;
			for (i=0; i<tela[n].all_r; i++) {
				for (j = 1; j < tela[n].all_k-1; j++) {
					if (i==0 && tela[(p=m+j)].all_k != tela[(q=n+j)].all_k 
								&& tela[p].stat == 'f' && tela[n].stat != 'f') {
						clearall_tela(n,1,-1, TWO);		/* O-F-F, ONE, OR TWO */
						tela[n].stat = '-';
						if (dev_print(TELA,__LINE__)) {
							printf("mark_tela() at n=%d evaluating cycling differences "
									"but calling clearall only at n to save fractal TR at p=%d.", n, p);
						}
					}
					else if (tela[(p=m+i*k+j)].all_k != tela[(q=n+i*k+j)].all_k) { 
						if (dev_print(TELA,__LINE__)) {
							printf("mark_tela() at n=%d evaluating cycling differences "
									"and calling clearall at p=%d and q=%d.", n, p, q);
						}
						if (tela[p].all_k)
							clearall_tela(p,1,-1, TWO);		/* O-F-F, ONE, OR TWO */
						if (tela[q].all_k)
							clearall_tela(q,1,-1, TWO);		/* O-F-F, ONE, OR TWO */
					}
				}
			}
		}
	}

	/* FIND AND RECORD MAX_K SIZE */
	int max_k = 0;
	for (n=0; n<lenseq; n++) {
		if (tela[n].all_k > max_k) {
			max_k = tela[n].all_k;
		}
	}
	options[1][46] = max_k;  

	if (dev_print(TELA,__LINE__)) {
		printf("                     Finishing. print_tela() follows.");
		print_tela( 0,  56);
	}
}

short int check_solo(int pos)
{
	if (!tela[pos].all_k) {
		return(-1);					/* NEGATIVE SIGNALS NO REPEAT AT GIVEN POSITION */
	}
	else if ( tela[pos-1].all_k || 
			  tela[pos+1].all_k ) {
		return(0);					/* ZERO SIGNALS NOT SOLO REPEAT POSITION */
	}
	else {
		return(1);					/* POSITIVE SIGNALS IS SOLO */
	}
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
	printf("\n t:");
	for (i=a; i<=b; i++) {
		if (tela[i].c != tela[i].t)
			printf("__%c", tela[i].t);
		else
			printf("___");
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

	printf("\nst:");
	for (i=a; i<=b; i++) {
		if (tela[i].stat)
			printf("%3c", tela[i].stat);
		else
			printf("  .");
	}
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
	for (f=1; f<=10; f++) {
		printf("\nf%c:", mha_base62(f));
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

	/* IF IMPERFECT REPEAT, ERASE. THIS IS MOSTLY RIGHT, BUT EVENTUALLY NEED TO CHECK OVERLAPPING TRANSITIONS */
	if (k*r > tela[n].all_S) {	
		for (i=m; i<(n+k*r); i++) {
			consensus[(tela[i].x)] = tela[i].t = tela[i].c;
		}
	}

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
			if      (tela[n1+i].c == tela[n2+i].c) {
				;
			}
			else if (OFF && tela[n1+i].e == tela[n2+i].e) {
				;
			}
			else if (tela[n1+i].t == tela[n2+i].t) {	/* THIS IS REDUNDANT IF ].e IS EVALUATED BUT HERE FOR TESTING PURPOSES */
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

	if (!violation) {		/* ELSE ASSIGN EQUIVALENCE HERE PROVIDED AXIOM OF CONTINUITY CHECKED OUT (1) OR WAS NOT CHECKED (0) */
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

/*********************************************/
int score_kmer(int n, int k, short unsigned int mode)
{
	if (!mode) {
		return(0);		/* O-F-F MODE */
	}
	else {
		int i;
		int score = 0;
		int match = MATCH;
		int transition = TRANSITION;

		if (mode == 2) {
			for (i=0; i<k; i++) {
				if      (tela[n+i].c == tela[n-k+i].c)
					score += match;
				else if (tela[n+i].e == tela[n-k+i].e)
					score += transition;
				else {
					score = 0;
					break;
				}
			}
			return(score);
		}
		else {
			for (i=0; i<k; i++) {
				if      (tela[n+i].c == tela[n-k+i].c)
					score += match;
				else {
					score = 0;
					break;
				}
			}
			return(score);
		}
		return(score);
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

	if (dev_print(TELA,__LINE__)) {
		printf("settle_tiescores() engaged at n=%d, max_score=%d, iteration=%d.", n, max_score, iteration);
	}

	for (i=n; i<n+span; i++) {
		if (tela[i].all_S == max_score) {
			k = tela[i].all_k;
			m = i-k;					/* DEFINES START OF FIRST REPEAT UNIT */
			o = i+k*(tela[i].all_r-1);	/* DEFINES START OF LAST REPEAT UNIT */
			up = m - k*iteration;		/* DEFINES THE GHOST FLANKING UNIT STARTING AT m-1 */
			dn = o + k*iteration;		/* DEFINES THE GHOST FLANKING UNIT STARTING AFTER REPEATS */
			for (j=0; j<k; j++) {
				if (up+j >= 0 && tela[up+j].e == tela[m-k*(iteration-1)+j].e) {
					if (tela[up+j].c == tela[m-k*(iteration-1)+j].c) {
						tela[i].all_Z += match;
						ratchet++;
					}
					else {
						tela[i].all_Z += transition;
						ratchet++;
					}
				}
				if (dn+j<=lenseq && tela[dn+j].e == tela[o+k*(iteration-1)+j].e) {
					if (tela[dn+j].c == tela[o+k*(iteration-1)+j].c) {
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

#endif		/*!FILE_TELA_SEEN */

/*************************************************************************************************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal */
/*************************************************************************************************************/
