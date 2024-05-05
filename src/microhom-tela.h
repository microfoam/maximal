/******************************************************************/
/* microhom-tela.h header file                                    */
/* This file has MHA functions either in the form of verb_tela(), */
/* or verb_object(), where object is related to tela struct.      */
/******************************************************************/

#ifndef FILE_TELA_SEEN
#define FILE_TELA_SEEN

short unsigned int  assign_tela(int pos, int eM, int eN, int mode);
void 				assign_transit(int n, int kr_src);
short int 			check_solo(int pos);
int  				check_tela(int eM, int eN, short unsigned int mode_dim);
void 				clearall_tela(int n, int span, int keep_score, int mode);
void				clearmem(void);
int  				cyclelize_tela(int cpos, int delta, int npos);
void				push_tela_or(int n);
void 				flatline_after_TR(int pos);
void				mark_tela(void);
void				push_mem(int pos, int row);
int 				push_gPnt(short unsigned int ymode, int pos, int prev_par);
int 				push_gPnt_kmer(int pos, int kmer, int reps);
void 				print_tela(int a, int b);
void 				pull_tela(int n);
int  				push_tela(int n2, int n1, short unsigned int axioms);
int 				score_kmer(int n, int k, short unsigned int mode);
int  				settle_tiescores(int n, int span, int max_score, int iteration);
int					update_tela(void);


/******************************************/
void clearmem(void) {
int i=0, j=0;
	for (i=0; i<MEMROWS; i++) {
		for (j=0; j<=Clean.pass_W; j++)
			tela[j].mem[i] = '\0';
	}
}


/******************************* MODE ZERO (+x) OR ONE (+y) *****************************************/
int push_gPnt(short unsigned int ymode, int pos, int prev_par)
{
	if (!ymode && pos == prev_par) {	
		/* No paralogy, increment in the x-direction */
		tela[pos].gPnt.rel_xy = 0;
		tela[pos].gPnt.prevPar = tela[pos].gPnt.topPar = pos;
	}
	else if (prev_par < pos) {
		if (ymode) {
			/* Paralogy, increment in the y-direction */
			tela[pos].gPnt.rel_xy = 1;
		}
		else {
			/* Paralogy, increment in the x-direction */
			tela[pos].gPnt.rel_xy = 0;
		}
		tela[pos].gPnt.topPar = tela[prev_par].gPnt.topPar;
		tela[pos].gPnt.prevPar = prev_par;
	}

	return((tela[pos].gPnt.topPar));
}


/***********************************************/
int push_gPnt_kmer(int pos, int kmer, int reps)
{
	int top_left = push_gPnt(YDIR, pos, pos-kmer);
	int top_right = 0, r=0, i=0;

	for (r=0; r<reps; r++) {
		if ( (top_right=push_gPnt(YDIR, pos+r*kmer, pos-kmer)) < 0 )
			break;
		
		for (i=1; i<kmer; i++) {
			if ( (top_right=push_gPnt(XDIR, pos+r*kmer+i, pos-kmer+i)) < 0 )
				break;
		}
	}
	
	if (top_left<0 || top_right<0)
		return(-1);
	else
		return(top_left);
}


/****************************************************************/
short unsigned int assign_tela(int pos, int eM, int eN, int mode)
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
	if (!kr_src) {				/* kr_src = ZERO MODE; DEV. FEATURE */
		return;
	}

	int i=0, j=0;
	int k, r;

	if      (kr_src == 1) {		/* kr_src = ONE */
		k = tela[n].ok;
		r = tela[n].or;
	}
	else if (kr_src == 2) {		/* kr_src = TWO */
		k = tela[n].cyc_k;
		r = tela[n].cyc_r;
	}
	else {
		k = tela[n].k;			/* kr_src = 3 (OR > 2) */
		r = tela[n].r;
	}

	int m = n - k;
	int m_pos=0, n_pos=0;

	/* ASSIGN TRANSITIONS TO .t IF IMPERFECT_TR, ALL IN REFERENCE TO FIRST UNIT STARTING AT m */
	for (i=0; i<r; i++) {
		for (j=0; j<k; j++) {	/* 1ST TIME TO NOTE THE TRANSITION POSITIONS BY DIFFERENCES */
			m_pos = m     + j;
			n_pos = n+i*k + j;
			if (tela[n_pos].c != tela[m_pos].c && tela[n_pos].e == tela[m_pos].e) {
				if (tela[m_pos].t != tela[m_pos].c)
					tela[n_pos].t = tela[n_pos].e;
				else
					tela[m_pos].t = tela[m_pos].e;
			}
		}
	}
	for (i=0; i<r; i++) {
		for (j=0; j<k; j++) {	/* 2ND TIME TO PROPAGATE TRANSITIONS TO ALL PARALAGOUS POSITIONS */
			m_pos = m     + j;
			n_pos = n+i*k + j;
			if (tela[n_pos].t != tela[m_pos].t && tela[n_pos].e == tela[m_pos].e) {
				tela[n_pos].t = tela[m_pos].t = tela[n_pos].e;
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
	if (mode_dim > 2)
		mode_dim = 2;	/* Unspecified non-zero modes --> reset to 2 */

	if (!mode_dim)		/* Check_tela called in O-F-F mode; will return 1+2 = success */
		return(3);
	else {
		int i=0, j=0, axioms=0, badflag=0;
	
		if (eM>=eN) {
			if (dev_print(TWO,__LINE__)) {
				printf("Need to call Check_tela explicitly with %d-D positions eN > eM.", mode_dim);
			}
			return(0);
		}
	
		if (mode_dim==2) {					/* Check_tela called with 2D row and column coordinates */ 
			/* SAVE 2-D COORDINATES */
			int lenseq = Clean.pass_W;

			/* TRANSLATE 2-D COORDINATES INTO 1-D COORDINATES */
			while (tela[i].y != eM && i<lenseq)
				i++;
			eM = i;
	
			j=lenseq;
			while (tela[j].x != eN && j)
				j--;
			eN = j;
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
			else
				break;
		}
		if (i==eN)
			axioms = 1;
		else if (dev_print(TWO,__LINE__)) {
			printf("Check_tela(mode_dim=%d): Problem of continuity at 1D positions %d-->%d (columns %d and %d)",mode_dim,i,i+1,tela[i].x,tela[i+1].x);
		}

		/* AXIOM TWO: EQUIVALENCE */
		for (i=eM; i<eN && !badflag; i++) {
			for (j=i+1; j<eN; j++) {
				if (tela[j].x == tela[i].x && 
					tela[j].e != tela[i].e) {
					badflag++;
					break;		/* TO BREAK FOR j LOOP */
				}	
			}
			if (badflag)
				break;
		}
		if (!badflag) 
			axioms+=2;
		else if (dev_print(TWO,__LINE__)) {
			printf("Check_tela(mode_dim=%d): Problem of equivalence at 1-D positions %d and %d (both in column %d)", 
								mode_dim, i, j, tela[i].x);
		}

		return(axioms);	/* 0 IF BOTH FAIL; +1 IF ONLY ONE PASSES; +2 IF ONLY TWO PASSES; +3 IF BOTH PASS */ 
	}
}


/* FUNCTION TO CLEAR_ALL ELEMENTS IN TELA WITHIN A NON-CONFLICTED CYCLING ISLAND W/ A SINGLE MAX-SCORE.        */
/* The use of clearall-tela has evolved in practice. As written it scans a window starting at n of length span */
/* and then clears elements at tela[n] if tela[n].all_Z did not match 'keep_score'. In practice, it was        */
/* very easy to start specifying spans of 1 (meaning only scan at n), and to set the keep_score to -1, which   */
/* is not a value that is ever assigned. Almost always it is now run as clearall_tela (n, 1, -1, TWO).         */
void clearall_tela(int n, int span, int keep_score, int mode)
{
	int i;
	if (!mode)		/* MODE 0 = O-F-F; 1=CLEAR S ONLY; 2=CLEAR ALL; 3=CLEAR ALL & SHIFT k-mers ok<- k2 ->k1->k0 */
		return;

	for (i=n; i< n+span; i++) {
		if (tela[i].all_Z != keep_score) {
			tela[i].all_S = 0;										/* MODES > 0 */					
			if (mode>1) {
				tela[i].k0 = tela[i].k1;							/* MODES > 1 */
				tela[i].ok = tela[i].or = tela[i].all_Z = 0;		/* MODES > 1 */  
				tela[i].all_L = tela[i].all_R = 0;					/* MODES > 1 */ 

				if (mode>2 && tela[i].k1 && !tela[i].k0) {
					tela[i].ok = tela[i].k2;						/* MODES > 2 */ 
					tela[i].k0 = tela[i].k1;						/* MODES > 2 */ 
					tela[i].k1 = tela[i].k2;						/* MODES > 2 */ 
					tela[i].k2 = 0;									/* MODES > 2 */ 
					push_tela_or(i);
				}
			}
		}
	}
}


/* FUNCTION TO CYCLELIZE UPSTREAM TR MARKED AT CPOS BY DELTA SO DOESN'T CONFLICT WITH TR MARKED AT NPOS */
int cyclelize_tela(int cpos, int delta, int npos)
{
	int lenseq = Clean.pass_W;

	if      ( cpos>lenseq ||  cpos<0 || tela[cpos].cyc_l < 2)
		return(0);
	else if (delta>lenseq || delta<1 || delta > tela[cpos].cyc_l)
		return(0);

	int    k = tela[cpos].k;
	int reps = tela[cpos].r;
	int i, j, m, n, r;
	char blnk = Fill->sym;		/* opt_B blank character */
	char c;
	short int nuctype = Clean.pass_V;
	short int nuctransit = 0;

	if (nuctype==1)
		nuctransit++;

	for (i=cpos; i<npos; i++) {
		tela[i].gPnt.rel_xy = XDIR;									/* CLEARS rel_xy --> 0; YDIR = 1 */
		tela[i].gPnt.prevPar = tela[i].gPnt.topPar = i;				/* RESETS prev_Par and topPar to default no paralogy */
	}
	push_gPnt_kmer(cpos+delta, k, tela[cpos+delta].or);

	if (k && reps) {
		int cycle_end = npos;		/* TMP SAFE-INITIALIZATION */
		for (r=0; r<reps; r++) {
			for (j=0; j<delta; j++) {
				c = tela[(i=cpos+r*k+j)].c;
				m = tela[i].y;
				n = tela[i].x;

				align2D[m][n] = blnk;
				m -= 1;
				n += k;
				tela[i].y = m;
				tela[i].x = n;
				align2D[m][n] = c;
			}
			if (r < reps - 1) {
				align2D[m][n+1] = slip.sym;
				align2D[m][n+2] = '\0';
			}
			else if (r == reps-1) {
				cycle_end = cpos+delta+k*tela[(cpos+delta)].or;
				for (j = cycle_end; j < npos; j++) {
					align2D[m][++n] = tela[j].c;
					tela[j].y = m+1;
				}
			}
		}
		m = tela[npos].y = tela[(npos - tela[npos].k)].y + 1;
		    tela[npos].x = tela[(npos - tela[npos].k)].x;

		int mpos = npos - tela[npos].k;
		int x_start = tela[mpos].x;
		int y_start = tela[mpos].y;
		for (i=1; i< tela[npos].k; i++) {
			tela[mpos+i].x = tela[npos+i].x = x_start + i;
			tela[mpos+i].y = y_start;
			tela[npos+i].y = y_start + 1;
		}

		for (j=0; j<lenseq; j++)
			align2D[m][j] = '\0';

		flatline_after_TR(npos);	
		tela[cpos].cyc_o = cyc_skip.sym;
		tela[cpos+delta].cyc_o = cyc_take.sym;

		/* DEAL WITH TRANSITION MARKS IF ANY */
		if (nuctransit) {
			for (j=npos-tela[npos].k; j<npos; j++) {
				if (tela[j].c != tela[j].t)
					consensus[(tela[j].x)] = tela[j].t;
				else
					consensus[(tela[j].x)] = '\0';
			}
		}

		return (1);		/* RETURN SUCCESS, BUT EVENTUALLY ADD A CHECK_TELA CALL IN HERE */
	}
	else
		return (0);
}

/******************************/
void flatline_after_TR(int pos)
{
	int lenseq = Clean.pass_W;
	int i, start, x, y;

	if (!tela[pos].ok) {
		for (i=pos-1; i>0; i--) {
			if (tela[i].ok) {
				pos = i;
				break;
			}
		}
		if (i==0) {
			return;
		}
	}

	x = tela[pos].x + tela[pos].ok;						/* START OF FIRST COLUMN AFTER REPEAT */
	y = tela[pos].y + tela[pos].or - 1;					/* FINISHING ROW */

	start = pos + tela[pos].ok * tela[pos].or;		/* 1D-POINT TO START */

	for (i=start; i<=lenseq; i++) {
		tela[i].x = x++; 
		tela[i].y = y; 
		tela[i].gPnt.rel_xy = 0;
		tela[i].gPnt.prevPar = tela[i].gPnt.topPar = i;
	}
	if (dev_print(TWO,__LINE__)) {
		printf("Function cyclelize_tela() finished by flat-lining tela at position=%d with coordinates (x,y) = (%d,%d).", i,x,y);
	}
}

/* Get repeats at n, Albert-style (units start at m, annotated at n, beginning of second unit); runs on operational k (ok) */
/* Prototype will focus only on perfect repeats. */
void push_tela_or(int n)
{
	int i=0, reps=0;
	int k=tela[n].ok;
	short unsigned int rcheck = 1;
	int m=n-k;
	int lenseq = Clean.pass_W;

	if (m<0 || n+k>lenseq) {
		tela[n].or = 0;
		return;
	}

	while (rcheck && n + k*(reps+1) < lenseq) {
		for (i=0; i<k; i++) {
			if (tela[m+reps*k + i].c != tela[n+reps*k + i].c) {
				rcheck = 0;
				break;	
			}
		}
		if (i==k)
			reps++;
	}
	tela[n].or = reps;
	tela[n].all_S = reps * k * MATCH;
	return;
}


/* Return k-mer repeat size smaller than k if it exists, otherwise return 0 */
int next_k(int n, int k1, short unsigned int seqtype) 
{
	int i=0, m=0, k=0, runs;

	for (k=k1-1; k>1; k--) {
		m = n-k;

		/* CHECK MONO-CHARACTER RUNS */
		runs = 0;
		for (i=0; i<k; i++) {
			if (tela[m+i].c!=tela[n+i].c)
				break;
			if (i && tela[n+i].c == tela[n+i-1].c) {
				if (!runs)
					runs=2;
				else
					runs++;
			}
		}
		if (runs==k)
			return(0);

		if (seqtype==1) {
			int transits=0;
			int maxtransits = allowed_transits(k);
			short int pur, pyr;
			pur = pyr = 0;

			for (i=0; i<k; i++) {
				if (tela[m+i].c=='n'||tela[n+i].c=='n')
					break;
				if (tela[m+i].c != tela[n+i].c) {
					if (opt_x.bit && k>opt_b.val && tela[m+i].e == tela[n+i].e) {
						if (++transits > maxtransits) 
							break;
					}
					else if (k>opt_b.val && tela[m+i].e == tela[n+i].e) {
						if (++transits >= maxtransits) 
							break;
					}
					else
						break;
				}

				if (tela[m+i].c == ambig.sym || tela[n+i].c == ambig.sym)
					break;

				if      (tela[m+i].c=='G')
					pur++;
				else if (tela[m+i].c=='C')
					pyr++;
				else if (tela[m+i].c=='A')
					pur++;
				else if (tela[m+i].c=='T')
					pyr++;

				if      (tela[n+i].c=='G')
					pur++;
				else if (tela[n+i].c=='C')
					pyr++;
				else if (tela[n+i].c=='A')
					pur++;
				else if (tela[n+i].c=='T')
					pyr++;
			}
			if (i==k) {
				if (transits && ((pur && !pyr) || (pyr && !pur)) && !next_k(n,k,seqtype)) {
					tela[n].k0 = k;		/* SAVE k-MER FOR CINCH-K; STORING FOR DE-BUGGING ONLY, NOT NEEDED */
					return(0);
				}
				else if (transits)
					return(-k);			/* RETURN NEGATIVE VALUE TO INDICATE IMPERFECT TR */
				else 
					return(k);
			}
		}
		else if (seqtype==3 && !opt_x.bit) {	/* PROTEIN */
			for (i=0; i<k; i++) {
				if (tela[m+i].e != tela[n+i].e)
					break;
			}
			if (i==k && k<k1)
				return(k);
		}
		else {									/* ELSE IF RNA (2) OR NON-BIOLOGICAL (0), OR opt_x */
			for (i=0; i<k; i++) {
				if (seqtype==2 && (tela[m+i].c=='n'||tela[n+i].c=='n'))
					break;
				if (tela[m+i].c != tela[n+i].c) 
					break;
			}
			if (i==k && k<k1)
				return(k);
		}
	}
	return(0);
}


/**************** FUNCTION TO CHECK IF POSITION i IS FRACTAL TO PARENT AT n ************************************/
/**************** RETURN(0) IF NOT; RETURN(candk) IF YES.                   ************************************/
int isfractal(int i, int n, int k) {

		int candk=tela[i].ok;					/* candidate k: k-mer size of candidate fractal */
		short int mfract=0;						/* switches --- perfect=1; imperfect=-1; nfract not necessary */

		if (!candk) {
			if (tela[i].impk<0) {				/* candidate k is an imperfect TR */
				candk = -tela[i].impk;
				mfract = -1;
			}
			else
				return(0);
		}

		if (candk>=k)
			return(0);
		else if (i-candk<n-k)
			return(0);
		else if (i+span_ork(i)>n)
			return(0);
		else if (candk==tela[i+k].k1 || candk==tela[i+k].k2)	/* Note: tela[i+k].ok might not be filled yet. */
			return(candk);
		else if (candk == -tela[i+k].impk) {
			if (mfract>0)
				return (candk);
			else {								/* ELSE both are imperfect so check to see they are concordantly so */
				int transits=0, j=0;
				for (j=i-candk; j<i+candk; j++) {
					if (tela[j].c!=tela[j+k].c && tela[j].e==tela[j+k].e)
						++transits;
					if (transits>2)				/* APPROX. MOSTLY CORRECT. TAGGED: <MAGIC> */
						return(0);
				}
				return(candk);
			}
		}
		else
			return(0);
}

/****** MARK FRACTALS AND PARENT AT POSITIONS mf (m-unit fractal), mf+k, AND p, RESPECTIVELY */
void makefract(int p, int k, int mf)
{
	push_mem(p,0); push_mem(mf,0); push_mem(mf+k,0);
	push_mem(p,5); push_mem(mf,6); push_mem(mf+k,6);
	tela[p].stat = st_parent.sym;
	tela[mf].statf  = tela[mf+k].statf  = st_fract.sym;
}

/******** BOTTOM-UP k LOOP ****************/
int get_unitk(int n)
{
	int i,k;

	for (k=1; k<WIDTH/2; k++) {
		if (n-k<0 || n+k>Clean.pass_W)
			return(0);
		for (i=n; i<n+k; i++) {
			if (tela[i].c != tela[i-k].c)
				break;
		}
		if (i==n+k)
			return(k);
	}
	return(0);
}

/**************** FUNCTION TO MARK ALL POSSIBLE k-MERs BEFORE LEGACY CINCH-T PASS ******************************/
void mark_tela(void) 
{
	int i, j, l, m, n, k, p, reps, span, smallest_k; 
	int threshold=0, max_score=0, max_count=0;
	int projection=0, projector=0, proj_k=0, fract_k=0;
	unsigned short int skip_break=0;
	int lenseq = Clean.pass_W;
	unsigned short int nuctype = Clean.pass_V;
	unsigned short int nuctransit=0, TRcheck=0, imperfect_TR=0, Aimperfect_TR=0, gapcheck=0;
	int homopoly_flag=0, Did=0, Dtr=0, Atr=0;
	unsigned short int checkconflict=0;
	int prev_k;
	int k1=0, k2=0, k_tmp=0;
	short unsigned int min_k=2; 		/* MINIMUM LIMIT k-MER SIZE MARKED; ALWAYS USE THIS INSTEAD OF VALUE */
	int first_imp=0, last_imp=0;

	if (nuctype == 1) {		/* IF DNA */
		nuctransit = 1;
	}

	/* Annotations of the largest (k1) and second largest (k2) k-mers at position n are permanent. 		*/
	/* Annotation at all_k (now ok) is eraseable and is the operational k-mer used at that position.	*/
	for (n=1; n<lenseq; n++) {
		k2 = k1 = k_tmp=0;
		for (k=WIDTH; k>=min_k; k--) {
			k = k_tmp = next_k(n,k,nuctype);
			if (k_tmp<0) {
				k = abs(k_tmp);
				k1 = next_k(n,k,nuctype);
				if (!k1 || (k1>0 && 3*k1<k)) {
					tela[n].impk = k_tmp;
					if (k1)
						k1 = 0;
					else
						break;
				}
				else if (k1>0)
					k = tela[n].k1 = k1;
			}
			else if (!k_tmp)
				break;

			if  (k_tmp>=min_k || k1>=min_k) {

				/* CHECK FOR HIGHER ORDER ARTIFACTS. EG., FOR [(GA)2]2 SKIP k=4 AND GO TO k=2 */
				int unit_k;
				if ((unit_k=get_unitk(n)) && unit_k<k_tmp && !(k_tmp%unit_k)) {
					m = n-k_tmp;
					int loops = k_tmp/unit_k;
					for (l=1; l<loops; l++) {
						if (tela[m+l*unit_k].k1 != unit_k)
							break;
					}
					if (l==loops) {
						tela[n].k1 = unit_k;
						break;
					}
				}

				if ( (k_tmp && !checkfractals_in_imperfect(k_tmp,n)) ||
					 (k1    && !checkfractals_in_imperfect(k1   ,n)) ) {
					if (k1>=min_k)
						tela[n].k1 = k_tmp = k1;
					else
						tela[n].k1 = k1 = k_tmp;

					if ((k_tmp=next_k(n,k1,0))>=min_k) {
						if (k_tmp==tela[n-1].k1 && k1==2*k_tmp && !(k1%k_tmp)) {
							for (i=0; i<k_tmp; i++) {
								if (tela[n+i].c!=tela[n+i+k_tmp].c)
									break;
							}
							if (i==k1-k_tmp)
								tela[n].k1 = k1 = k_tmp;	/* BECAUSE LARGER k1 IS A HIGH-REPEAT NUMBER ARTIFACT */
							else {
								tela[n].k2 = k2 = k_tmp;
								break;						/* BREAK BECAUSE CURRENTLY ONLY STORING TOP TWO k-MERS AT n */
							}
						}
						else {
							/* CLEANS UP TELA ANNOTATION IN RARE CASES; NOT MEASUREABLY HELPFUL ATM v4.34 7/25/2020 */
							if (k_tmp==tela[n-k1].k1 && tela[n-1].k2==k_tmp && tela[n-1-k1].k1!=k_tmp)
								tela[n-1].k2 = '\0';

							tela[n].k2 = k2 = k_tmp;
							break;						/* BREAK BECAUSE CURRENTLY ONLY STORING TOP TWO k-MERS AT n */
						}
					}
					else
						break;
				}
			}
		}
	}

	if (nuctransit) {
		for (n = 1; n<lenseq; n++) {
			if (tela[n].impk) {

				/* CANCEL IMPERFECT CYCLING ISLANDS CONTAINING A PERFECT COLUMN OF SIZE > k/2 */
				if (tela[n].k1*2 > -tela[n].impk) {
					k = tela[n].impk;
					
					i = 0;	/* CANCEL DOWNSTREAM PART OF ISLAND */
					while (tela[n-i].impk == k)
						tela[n - i++].impk = '\0';
					
					i = 1;	/* CANCEL UPSTREAM PART OF ISLAND */
					while (tela[n+i].impk == k && i<lenseq)
						tela[n + i++].impk = '\0';
				}
				/* CANCEL IMPERFECT IF ADJACENT COLUMN IS PERFECT */
				else if (tela[n-1].k1 || tela[n+1].k1)
					tela[n].impk = '\0';
			}
		}
	}

	tela[0].isl = 0;	/* WILL STORE NUMBER OF INDEPENDENT ISLANDS OF OVERLAPPING TR's. */

	for (n = 1; n<=lenseq; n++) {
		for (m = 0; m < n; m++) {
			/* SLIDE DOWN TO ROW WITHIN POPULATED HEMIDIAGONAL */
			if (n-m > WIDTH+1) 
				m = n-WIDTH;
			
			/* SET K-MER SIZE AND DTHR SCORE THRESHOLD */
			k = n-m;

			if (tela[n-1].ok && tela[n-1].k1 && tela[n-1].k1 == tela[n].k1 && k>tela[n].k1 && !(k%tela[n-1].ok)) {
				k = tela[n].k1;
				m = n - k;
			}

			/* SKIP k = O N E */
			if (k==1 || k==tela[n].k0)
				break;	/* GO TO NEXT n */
			else if (nuctransit) 
				threshold = score_DTHR(k);

			/* SET HOMOPOLYMER RUN STATUS UNKNOWN; USED TO RULE OUT k>1 MONONUCLEOTIDE "REPEATS" */
			homopoly_flag = 2;
			if (tela[n].c != tela[n-1].c)
				homopoly_flag = 0;

			/* START COUNTING SCORE (EQUIV. TO: IF PATHBOX POSITION HAS VALUE > MISMATCH) */
			if (n+k <= lenseq) {
				Dtr = imperfect_TR = 0;		/* INITIALIZATION */

				/* IF SUMMING PATHBOX DIAGONAL 1/4: COMPUTE SCORES OF IDENTITY LINE AND REPEAT DIAGONAL*/
				Did = k*MATCH;
				for (j = 0; j < k; j++) {
					if (nuctransit && (tela[m+j].c==ambig.sym || tela[n+j].c==ambig.sym)) {
						Dtr =  0;
						break;							
					}
					else if (tela[m+j].c == tela[n+j].c) 
						Dtr += MATCH;	
					else if (nuctransit && tela[m+j].e == tela[n+j].e) 
						Dtr += TRANSITION;
					else {
						Dtr = 0;
						break;
					}

					/* SET HOMOPOLYMERIC RUN BIT TO FALSE IF NOT A POSSIBILITY */
					if (homopoly_flag && j > 0 && tela[m+j].c != tela[m+j-1].c)
						homopoly_flag = 0;
				}

				/* IF SUMMING PATHBOX DIAGONAL 2/4: SET HOMOPOLYMERIC RUN BIT TO TRUE IF DETECTED 	*/
				if (j==k && homopoly_flag) {
					homopoly_flag = 1;				/* BIT IS THERE IF NEEDED BEYOND BREAK. 		*/
					Dtr = 0;
					if (tela[n].k1==k)
						tela[n].k1 = '\0';
					break;							/* GO TO NEXT n b/c will be true for all remaining rows */
				}

				/* IF SUMMING PATHBOX DIAGONAL 3/4: IF CONSIDERING NUCL. TRANSITIONS AS PARTIAL MATCHES */
				if (nuctransit && Dtr && Dtr!=Did) {
					if (k>opt_b.val && 100*Dtr/Did > threshold && tela[n].impk==-k && !tela[n].k1) {
						if (tela[n+1].impk==-k) {
							imperfect_TR = 1;
							i = 0;
							while (tela[n+i].impk==-k) {	/* CHECK IMPERFECT ISLAND FOR PERFECT PEAK */
								if (tela[n+1 + i++].k1==k) {
									int perfect_col = n+i;
									imperfect_TR = 0;
									Dtr = 0;
									while (tela[n+i-1].impk==-k && n+i>0)
										tela[n-1 + i--].impk = '\0';
									while ((tela[++perfect_col].k1==k || tela[perfect_col].impk==-k) && perfect_col<lenseq) {
										if (tela[perfect_col].impk==-k)
											tela[perfect_col].impk = '\0';
									}
									break;
								}
							}
						}

						if (Dtr) {
							for (i=m+1; i<n; i++) {
								if (tela[i].k1 && i-tela[i].k1<m)	/* GIVEN PERFECT NON-FRACTAL k-MER IN m-UNIT, SKIP IMPERFECT */
									break;
							}
							if (i==n)
								imperfect_TR = 1;
							else
								Dtr = 0;
						}
					}
					else 
						Dtr = 0;
				}

				/** MOD TESTS. Original example for if part: seq-146-v344_33-snippet.txt    **/
				/**            Original example for if else part: seq-15-cycle4-snippet.txt **/
				if (Dtr && (prev_k=tela[n-1].ok) && prev_k!=k) {
					short unsigned int k2_check=1;
					for (i=n; i<n+k; i++) {
						if (tela[i].k2) {
							k2_check=0;
							break;
						}
					}
					if (!imperfect_TR && k>prev_k && k%prev_k && !tela[n-1].k1) {
						clearall_tela(n-1,1,-1,TWO);
						push_mem(n-1,4);
					}
					else if (k2_check && k>prev_k && k%prev_k) {
						push_mem(n-1,1);
						tela[n-1].stat = st_Fract.sym;		/* st_Fract = orPHan Fractal */
						push_mem(n,0);						/* ROW ZERO IS FOR ALL MARKS, NOT JUST THOSE SLATED FOR CLEARALL */
						push_mem(n,1);
						if (k<9 && k<3*prev_k) {			/* HEURISTIC: IF SMALL ks ARE CLOSE IN SIZE BETTER TO SKIP LARGER k-MER. SEE ABBA-ZABBA */
							n += k-1;
							Dtr = 0;
						}
					}
					else if (prev_k>k && prev_k % k && tela[n-prev_k].ok != k) {
						if (tela[n-prev_k].k2==k) {
							push_mem(n  ,0);		/* ROW ZERO IS FOR ALL MARKS, NOT JUST THOSE SLATED FOR CLEARALL */
							push_mem(n  ,2);
							tela[n-prev_k].statf = st_fract.sym;
							tela[n       ].statf = st_fract.sym;
						}
						else {						/* CAN BE TURNED O F F WITH WIGGLE IN WCR (WIDER) BUT NOT EXTRA CHOWDER. v4.33, 7/12/2020 */
							push_mem(n  ,0);		/* ROW ZERO IS FOR ALL MARKS, NOT JUST THOSE SLATED FOR CLEARALL */
							push_mem(n  ,2);
							tela[n  ].stat = st_Fract.sym;
							if (tela[n-1-prev_k].stat == st_cycle.sym) {
								push_mem(n-1,2);
								tela[n-1].stat = st_Fract.sym;
							}
						}
						n += prev_k-2;
						Dtr = 0; 
					}
				}

				/* CHECK TO SEE IF THERE ARE CYCLING ISLANDS THAT ARE PARTIALLY FRACTAL TO A PARENT & CANCEL NON-FRACTAL PART IF SIMPLE. squid~34 */
				if (Dtr) {
					for (i=n+k-1; i>=n+min_k; i--) {
						if (tela[i].k1) {
							if ((fract_k=tela[i].k1)<k && i-fract_k>=n && tela[i-k].ok==fract_k && i-k+tela[i-k].k1<=n) {
								makefract(n,k,i-k);
								if (tela[i-k-1].ok==fract_k) {
									j = 1;
									while (OFF && tela[i-k - j].ok == fract_k) {
										push_mem(i-k - j,12);
										clearall_tela(i-k - j++,1,-1,TWO);
									}
									j = 1;
									while (tela[n+j].k1==k && i-fract_k<n+j) {
										push_mem(n + j,13);
										tela[n + j++].k1 = '\0';
									}
								}
							}
						}
					}
				}

				/* CHECK TO SEE IF THERE ARE FRACTAL REPEATS WITH BELOW THRESHOLD DOPPELGANGERS. EXAMPLE: GTGT IN ONE UNIT, GCGT IN THE ADJACENT UNIT */
				if (Dtr && imperfect_TR) {
					for (i=m+min_k; i<n; i++) {
						if (isfractal(i,n,k)) {
							makefract(n,k,i);
						}
                        else if (OFF && tela[i].statf!=st_fract.sym && (fract_k=tela[i].k1) && fract_k<=opt_b.val && !tela[i+k].k1 && i
                                    && tela[n].stat!=st_parent.sym && tela[n+1].k1!=k) {
                            push_mem(n,0);
                            push_mem(n,7);
                            tela[n].stat = st_parent.sym;
                            push_mem(i,8);
                            tela[i].stat   = st_Fract.sym;
                            tela[i].echoes = cyc_skip.sym;
                        }    
						else if ((fract_k=tela[i+k].k1) && fract_k<=opt_b.val && !tela[i].k1 && i+k-fract_k>=n && i+k+tela[i+k].k1<=n+fract_k 
								&& tela[n].stat!=st_parent.sym && tela[n+1].k1!=k) {
							push_mem(n,0);
							push_mem(n,7);
							tela[n].stat = st_parent.sym;
							push_mem(i+k,8);
							tela[i+k].stat   = st_Fract.sym;
							tela[i+k].echoes = cyc_skip.sym;
						}
					}
				}

				/* IF SUMMING PATHBOX DIAGONAL 4/4: START COUNTING REPEATS */
				if (Dtr && (Dtr==Did || imperfect_TR)) {
					/* COUNT NUMBER OF REPEATS ALBERT-STYLE */
					TRcheck = 1;
					tela[n].ok = k;
					tela[n].or = reps = 1;
					tela[n].all_S = Dtr;	/* SAVE INITIAL UNIT SCORE */
					first_imp = last_imp = 0;

					if (imperfect_TR)
						first_imp = last_imp = reps;

					while (TRcheck) {
						Atr = 0;
						if (m + (reps+1)*k > lenseq) { 
							Atr = 0;
							tela[n].or = reps;
							push_mem(n,0);		/* ROW ZERO IS FOR ALL MARKS, NOT JUST THOSE SLATED FOR CLEARALL */
							break;
						}

						if (nuctransit) {
							Atr = score_kmer(n+k*reps,k,TWO);
							if (Atr!=Did && (100*Atr)/Did>threshold)
								Aimperfect_TR = 1;
							else
								Aimperfect_TR = 0;
						}
						else
							Atr = score_kmer(n+k*reps,k,ONE);

						if (Atr==Did) {
							reps++;
							tela[n].all_S += Atr;
						}
						else if (Aimperfect_TR) {
							int test_k = 1;
							if (reps<2) {
								for (i=n+reps*k; i<n+reps*k+WIDTH && i<lenseq; i++) {
									if (tela[i].k1 || tela[i].impk) {
										if (!tela[i].k1)
											test_k = - tela[i].impk;
										else
											test_k = tela[i].k1;

										if (test_k>k && i-test_k < n+reps*k) {
											test_k = -1;
											break;
										}
									}
								}
							}
							if (test_k>0) {
								reps++;
								tela[n].all_S += Atr;
								if (!first_imp)
									first_imp = last_imp = reps;
								else
									last_imp=reps;
							}
							else
								Aimperfect_TR = TRcheck = 0;
						}
						else {		/* ELSE FINAL NUMBER OF REPEATS (REPS) IS NOW KNOWN *****************/
							if (tela[n].all_S) {

								/* CANCEL SINGLE IMPERFECT REP AT BEGINING OR END OF MULTIPLE PERFECT REPEATS */
								if (nuctransit && reps>2 && first_imp && first_imp==last_imp) {
									if (first_imp==reps) {
										reps--; 
									}
									else if (first_imp==1) {
										tela[n].or = tela[n].ok = tela[n].all_S = 0;
										break;
									}
								}

								tela[n].or = reps;
								push_mem(n,0);		/* ROW ZERO IS FOR ALL MARKS, NOT JUST THOSE SLATED FOR CLEARALL */

								/* BLOCK TO NUMBER INDEPENDENT ISLANDS. ISLAND = SPAN OF OVERLAPPING TRs */
								if (!tela[m].isl) {
									for (i=m+1; i<n; i++) {
										if (tela[i].isl)
											break;
									}
									if (i==n)
										++tela[0].isl;
								}
								else
									tela[0].isl = tela[m].isl;
								for (i=m; i<n+span_ork(n); i++)
									tela[i].isl = tela[0].isl;
	
								if (n+k*reps > projection) {
									/* BEFORE ADVANCING PROJECTION, CHECK TO SEE IF THIS TR CALL IS COVERING A FRACTAL REPEAT OF SMALLER K.  */	
									/* RECALL THAT A FRACTAL TR (k>1) CAN ONLY EXIST STARTING AT min_k+1 COLUMN OF EACH UNIT OF PARENT TR. */
									/* RECALL THAT A CYCLING TR IS CALLED AND MARKED ONLY WHEN THE SECOND CYCLING POSITION IS CALLED.        */
									if (k==tela[n-1].ok && tela[n-1].stat==st_cycle.sym && tela[n-proj_k].ok && tela[n-proj_k].ok<proj_k) {
										for (i=n; i<n+span_ork(n); i++) {
											if (tela[n].c != tela[n-proj_k].c)
												break;
										}
										if (i==n+span_ork(n))
											skip_break = 1;
									}
									if (!skip_break) {
										projector = n;
										projection = n + k*reps;
										proj_k = k;
										if (k == tela[n-1].ok && tela[n].stat!=st_parent.sym) {
											tela[n-1].stat = st_cycle.sym; 		/* c FOR TRIVIAL-CASE OF CYCLING FRAME TYPE REPEAT */
											tela[n  ].stat = st_cycle.sym;
										}
									}
								}
								else {
									if (tela[n-1].ok && (k == tela[n-1].ok || tela[n-1].ok % k == 0)) {
										if (tela[n-1].stat!=st_parent.sym)
											tela[n-1].stat = st_cycle.sym; 		/* c FOR TRIVIAL-CASE OF CYCLING FRAME TYPE REPEAT */
										if (tela[n  ].stat!=st_parent.sym)
											tela[n  ].stat = st_cycle.sym;
									}
									else if ((tela[n-proj_k].ok==k || tela[n-proj_k].k1==k) && n-k>=projector) {
										makefract(projector,proj_k,n-proj_k);
										i = n-proj_k-1;
										while (tela[i].ok == k && tela[i].stat == st_cycle.sym && i>=projector-proj_k && tela[i].or<2) {
											clearall_tela(i,1,-1,TWO);
											push_mem(i, 4);
											i--;
										}
									}
									else {
										tela[n].stat = st_Fract.sym;
										push_mem(i, 4);
									}
								}
							}
							TRcheck = 0;
						}
					}

					/* v4.30: MARK FRACTAL TR'S FOR CINCH-T TO SKIP, AND LEAVE FOR CINCH-K */
					if (n>=2*min_k) {
						for (i=m+1; i<n; i++) {	
							if ((fract_k=tela[i].ok) && tela[n].all_S > tela[i].all_S && i+span_ork(i)<=n && tela[i].ok!=k)  {
										if (i>=n+fract_k && tela[i].all_S==tela[i+k].all_S) {
											makefract(n,k,i);
										}
										else {
											tela[i].stat = st_Fract.sym;
											push_mem(i, 4);			/* MARKING IN CLEARALL ROW BUT NOT CLEARING */
										}
							}
						} /* END OF FOR LOOP THROUGH SHADOW STARTING AT m+1 */
					}
					if (!skip_break)
						break;				/* OTHERWISE MAY OVERWRITE TR WITH ONE OF SMALLER K */
					else
						skip_break = 0;		/* TO RE-INITIALIZE STATE */
				}
			}
		} /* END OF FOR m */
	} /* END OF FOR n */

	for (n=1; n<lenseq; n++) {
		if (tela[n].ok) {
			k = tela[n].ok;

			/* FILL IN CYCLING GAPS CAUSED BY BELOW THRESHOLD FRAMES: THIS SUPPRESSES INTRA-TR CONFLICT REPORTING */
			if (!tela[n+1].ok) {
				gapcheck = 0;
				for (i=n+2; i <= n + k * tela[n].or; i++) {
					if (tela[i].ok && tela[i].ok != k) {
						gapcheck = 0;
						break;
					}
					else if (tela[i].ok == k) {
						gapcheck = 1;
						break;
					}
				}
				if (gapcheck) {
					if (dev_print(TWO,__LINE__)) {
						printf("         Mark_tela filling in gap between %d and %d, inclusive of these points, for k-mer=%d.", i-1,n-1,k);
					}
					for (j=i-1; j>n; j--) {
						tela[j].ok = k;
						tela[j].or = 0;
					}
				}
			}
		}
	}

	for (n=2; n<lenseq; n++) {
		if (tela[n].ok) {
			/* CANCEL MARKS FOR k-MERS W/ TRANSITIONS WHEN EMBEDDED IN A CYCLING ISLAND WITH AT LEAST ONE COLUMN W/ PERFECT REPEATS */
			if (tela[n].k1 && tela[n].k1==tela[n].ok && (tela[n].k1 == -tela[n-1].impk || tela[n].k1 == -tela[n+1].impk)) {
				int lkmer = tela[n].ok;
				int lreps = tela[n].or;
				for (i=1; ;i++) {
					if (tela[n-i].impk == -lkmer && tela[n-i].or==lreps && tela[n-i].k1<lkmer) {
						tela[n-i].cyc_o = cyc_skip.sym;
						push_mem(n-i, 9);
					}
					else
						break;
				}
				for (i=1; ;i++) {
					if (tela[n+i].impk == -lkmer && tela[n+i].or==lreps && tela[n+i].ok==lkmer) {
						clearall_tela(n+i, 1, -1, TWO);
						push_mem(n+i, 9);
					}
					else if (tela[n+i].impk != -lkmer)
						break;
				}
			}

			/* CANCEL MARKS AT EDGES OF CYCLING ISLANDS THAT HAVE MORE TRANSITIONS THAN ADJACENT COLUMN */
			if ((tela[n-1].ok || tela[n+1].ok) && (!tela[n-1].ok || !tela[n+1].ok) && 
				tela[n].ok   == tela[n-1].ok    + tela[n+1].ok &&
				tela[n].or   == tela[n-1].or    + tela[n+1].or &&
				tela[n].all_S < tela[n-1].all_S + tela[n+1].all_S) {
					clearall_tela(n, 1, -1, TWO);		/* O-F-F, ONE, OR TWO */
					push_mem(n, 9);
			}
			/* CANCEL MARKS AT UPSTREAM EDGE OF CYCLING ISLAND IF THEY HAVE A TRANSITION OVERLAPPING AN UPSTREAM TR AND CAN CYCLE */
			else if (!tela[n-1].ok && tela[n].stat==st_cycle.sym && tela[(m=n-tela[n].ok)].ok && tela[m].ok<tela[n].ok && 
						tela[m].stat!=st_cycle.sym && tela[m].all_S==MATCH*span_ork(m) && tela[m].c != tela[n].c) {
				clearall_tela(n, 1, -1, TWO);		/* O-F-F, ONE, OR TWO */
				push_mem(n, 9);
			}
		}
	}

	/* MARK LOW-COMPLEXITY REGIONS WITH TWO OR MORE k-mers MARKED AT COLUMN. */
	/* HERE, LOW-COMPLEXITY MEANS REPEATS ARE BASED ON TWO SYMBOLS OR LESS.  */
	for (n=2; n<lenseq; n++) {
		if (tela[n].k1 && (tela[n].k2||(tela[n+1].k0 && tela[n+1].k0!=tela[n].k1))) {
			int symb1=0, symb2=0;		/* WILL USE INTS AND RELY ON CHAR VALUES */
			k = tela[n].k1;
			m = n-k;
			for (i=m; i<n+k; i++) {
				if (!symb1)
					symb1 = tela[i].c;
				else if (!symb2)
					symb2 = tela[i].c;
				else if (tela[i].c != symb1 && tela[i].c != symb2) {
					symb1 = symb2 = 0;
					break;
				}
			}
			if (i==n+k) {
				int max_S = tela[n].all_S;
				for (i=m; i<n+k; i++) {
					tela[i].statl = st_lowcm.sym;		/* MARK AS LOW-COMPLEXITY */
					if (tela[i].all_S > max_S)
						max_S = tela[i].all_S;
				}
				/* EXTEND MARKS UPSTREAM AND DOWNSTREAM TO MARK EXTENT OF LOW-COMPLEXITY "ISLAND" */
				/* THEN POKE A HOLE AT MAX_S */
				for (i=m-1; i>=0; i--) {
					if (tela[i].c == symb1 || tela[i].c == symb2) {
						tela[i].statl = st_lowcm.sym;
						if (tela[i].all_S > max_S)
							max_S = tela[i].all_S;
					}
					else
						break;
				}
				for (j=n+k; j<lenseq; j++) {
					if (tela[j].c == symb1 || tela[j].c == symb2) {
						tela[j].statl = st_lowcm.sym;
						if (tela[j].all_S > max_S)
							max_S = tela[j].all_S;
					}
					else
						break;
				}
				for (l=i+1; l<j; l++) {
					if (tela[l].all_S == max_S) {
						if (tela[l].stat != st_parent.sym)
							tela[l].statl = '\0';
						else
							tela[l].cyc_o = cyc_skip.sym;
					}
				}			
			}
		}
	}

	/* FRACTAL SPLITTING */
	for (n=lenseq; n>1; n--) {
		if (tela[n].ok && tela[n-1].ok<tela[n].ok) {
			k = tela[n].ok;
			m = n - k;
			if (k>4 && ((tela[m].or && tela[m].ok<k) || (tela[m-1].or && tela[m-1].ok<k) || tela[m].k0)) {
				int splitcol = m;
				if (tela[m].or>1 && tela[m].all_S<=tela[n].all_S && tela[m-1].ok!=tela[m].ok && splitcol+span_ork(splitcol)<=n) {
					int k_at_m = tela[splitcol].ok;
					int less_r = tela[splitcol].or - 1;
					int jump = 0;
					tela[splitcol].or = 1;
					i = 1;
					while (tela[(jump=(splitcol-i*k_at_m))].ok==k_at_m && tela[jump].stat==st_cycle.sym) {
						tela[jump].or -= less_r;
						++i;
					}
					jump = splitcol - (--i)*k_at_m;
					if (tela[--jump].ok == k_at_m)
						tela[jump].or -= less_r;
				}
				if (tela[m].ok && (tela[m-1].ok==tela[m].ok || (!tela[m].ok && tela[(splitcol=m-1)].ok))) {
					while (splitcol && tela[splitcol].ok && tela[splitcol-1].ok <= tela[splitcol].ok) {
						while (splitcol + span_ork(splitcol) - m > tela[splitcol].ok && tela[splitcol].or>1)
							tela[splitcol].or--;
						splitcol--;
					}
				}
				else if (tela[--splitcol].ok) {
					while (splitcol && tela[splitcol-1].ok==tela[m-1].ok)
						splitcol--;
					if (tela[splitcol].all_S <= tela[n].all_S) {
						while (splitcol + span_ork(splitcol) - m > tela[splitcol].ok && tela[splitcol].or>1)
							tela[splitcol].or--;
					}
				}
				else if (tela[m].k0 && tela[m].k0<k) {
					while (splitcol && (tela[splitcol-1].ok==tela[m].k0 || tela[splitcol-1].k0==tela[m].k0))
						splitcol--;
					if (tela[splitcol].all_S <= tela[n].all_S) {
						while (splitcol + span_ork(splitcol) - m > tela[splitcol].ok && tela[splitcol].or>1)
							tela[splitcol].or--;
					}
				}
			}
		}
	}

	/* NOW MARK ALL CONFLICTING TRs */
	for (n=lenseq; n>0; n--) {
		if (tela[n].ok) {
			k = tela[n].ok;
			m = n - k;
			j = n - 1;

			while (tela[j].ok==k) {		/* SKIP CYCLE COLUMNS OF SAME K-MER */
				j--;
			}
			for (i=j; i>0; i--) {
				if (tela[i].ok) {
					if (i + tela[i].ok*(tela[i].or-1) > m || (tela[n].statf!=st_fract.sym && i+tela[i].ok>n)) {
						tela[n].all_L = i;				/* UPDATE LEFT-MOST OVERLAPPING & CONFLICTING TR */
						if (!tela[i].all_R)
							tela[i].all_R = n;			/* UPDATE RIGHT-MOST OVERLAPPING & CONFLICTING TR */
						/* CASE OF NON-CONFLICTING FRACTAL REPEATS */
						if (tela[i].all_S<tela[n].all_S && (!tela[i].k1 || tela[i].k2)) {
							clearall_tela(i, 1, -1, TWO);		/* O-F-F, ONE, OR TWO */
							push_mem(i, 13);
						}
					}
				}
			}
			/* CANDIDATE FOR DELETION? */
			if (tela[m].or && tela[m].all_S < tela[n].all_S && m+tela[m].or*(tela[m].ok)>n) {
				clearall_tela(m, 1, -1, TWO);		/* O-F-F, ONE, OR TWO */
				push_mem(m, 9);
			}
		}
	}

	for (n=1; n<lenseq; n++) {
		if (tela[n].ok) {
			k = tela[n].ok;

			/* 9/30/2020 v4.35 WHEN DELETED CAUSES ONLY MINOR WCR EXPANSIONS IN SOME TESTS */
			/* 9/ 1/2020 v4.25 FIXES ONLY churly-14. */
			if (!tela[n-1].ok) {
				m = n-k;
				int prev_k;
				for (i=m+1; i<n-1; i++) {
					if (tela[i].ok && tela[i].statf!=st_fract.sym && (prev_k=tela[i].ok)<k && i-prev_k<m && 
							!tela[i].all_L && tela[i].all_S<tela[n].all_S && i+span_ork(i)<n) {
						clearall_tela(i, 1, -1, TWO);
						push_mem(i, 9);
						tela[i].echoes = cyc_skip.sym;
					}
				}
			}

			/* CANCEL MARK IF OVERLAPS PARENT W/ FRACTALS AND PARENT K > k */
			if (!tela[n-1].ok && tela[n].stat==st_cycle.sym) {
				for (i=n-1; i>n-k; i--) {
					if (tela[i].statf==st_fract.sym) {
						if (tela[i].ok*2>=k) {
							clearall_tela(n, 1, -1, TWO);		/* O-F-F, ONE, OR TWO */
							push_mem(n, 9);
							break;
						}
						for (j=i-1; j>0; j--) {
							if (tela[j].stat==st_parent.sym && tela[j].ok > k) {
								clearall_tela(n, 1, -1, TWO);		/* O-F-F, ONE, OR TWO */
								push_mem(n, 9);
								break;
							}
						}
					}
				}
			}
		}
	}

	/* IDENTIFY CYCLING MARKS FOR PERFECT REPEATS LAID OVER INTERNAL REPEATS WITH r>1 */
	/* CAN CLEAN UP HOW THESE LOOK AND PERHAPS ACT ON THEM PRIOR TO CINCH-K */
	for (n = 1; n<=lenseq; n++) {
		if (tela[n].or > 1 && (tela[n].all_S == tela[n].ok * tela[n].or)) {
			k    = tela[n].ok;
			reps = tela[n].or;
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
					tela[n+j].ok = tela[m+j].ok;
					tela[n+j].or = tela[m+j].or;
					tela[n+j].all_S = tela[m+j].all_S;
					tela[n+j].all_Z = tela[m+j].all_Z;
//					tela[n+j].all_L = tela[m+j].all_L;			/* 2024.05.05 APPEARS UNNECESSARY */
				}
				for (j=0; j<k; j++) {
					for (i=1; i<reps; i++) {
						tela[n+i*k+j].ok = tela[m+j].ok;
						tela[n+i*k+j].or = tela[m+j].or;
						tela[n+i*k+j].all_S = tela[m+j].all_S;
						tela[n+i*k+j].all_Z = tela[m+j].all_Z;
//						tela[n+i*k+j].all_L = tela[m+j].all_L;	/* 2024.05.05 APPEARS UNNECESSARY */
					}
				}
				n = n + k*reps - 1; 	/* NEED TO ADVANCE BEYOND CONFLICT-FREE LAST UNIT */
			}
		}
	}

	/* IDENTIFY IRREDEEMABLE CONFLICT ISLANDS IN THE SHADOW OF DOWNSTREAM CYCLING ISLANDS */
	for (n=0; n<=lenseq; n++) {
		if (tela[n].all_L && (k=tela[n].ok) && !tela[n-1].ok) {
			m = n-k;
			prev_k=0;	/* INITIALIZE TO ZERO STATE TO USE AS TEST FOR SETTING ONLY ONCE */

			for (i=m+1; i<n; i++) {
				if (!prev_k && tela[i].ok && tela[i].ok<k && i-tela[i].ok<m && tela[i].statf!=st_fract.sym) {
					while (tela[i-1].ok==tela[i].ok) {
						i--;
						if (tela[i].all_S>tela[n].all_S) {
							prev_k = -1;	/* USING AS FLAG TO BREAK OUT OF BOTH WHILE LOOP AND FOR LOOP THROUGH SHADOW */
							if (tela[n].k2) {
								clearall_tela(n, 1, -1, ONE);	/* O-F-F, ONE, OR TWO */
								push_mem(n, 15);				/* TMP NUMERIC ASSIGNMENT */
							}
							break; 
						}
					}
					if (prev_k<0)
						break;				/* BREAK OUT OF FOR LOOP THROUGH SHADOW */
					else {
						prev_k = tela[i].ok;			/* IDEA IS THAT PREV_K GETS ASSIGNED ONLY ONCE */
						if (tela[i].all_R==n) {
							clearall_tela(i, 1, -1, ONE);	/* O-F-F, ONE, OR TWO */
							push_mem(i, 15);				/* TMP NUMERIC ASSIGNMENT */
						}
					}
				}
				else if (prev_k && tela[i].ok==prev_k && i-tela[i].ok<m && tela[i].statf!=st_fract.sym && tela[i].all_R==n) {
					clearall_tela(i, 1, -1, ONE);	/* O-F-F, ONE, OR TWO */
					push_mem(i, 16);				/* TMP NUMERIC ASSIGNMENT */
				}
				else if (prev_k && tela[i].ok != prev_k)
					break;
				else if (tela[i].impk && -tela[i].impk<k && i+tela[i].impk<m) {
					int impspan = 1;
					while (tela[i-impspan].impk==tela[i].impk)
						impspan++;
					if (tela[i-impspan].all_S<tela[n].all_S) {
						for (j=i-impspan+1; j<m; j++) {
							clearall_tela(j, 1, -1, TWO);	/* O-F-F, ONE, OR TWO */
							push_mem(j, 16);				/* TMP NUMERIC ASSIGNMENT */
						}
					}
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
			smallest_k = k = tela[n].ok;
			max_score = tela[n].all_S;
			if (tela[n].all_L || tela[n].all_R) {
				checkconflict = 0;
			}
			for (i=n+1; tela[i].ok; i++) {
				span++;
				if (checkconflict) {
					if (tela[i].all_L || tela[i].all_R) {
						checkconflict = 0;
					}
					else {
						if (tela[i].all_S > max_score) 
							max_score = tela[i].all_S;
						if (tela[i].ok < smallest_k) 
							smallest_k = tela[i].ok;
					}
				}
			}
			if (!checkconflict) {						/* MEANING CHECKCONFLICT FLAG WAS NEVER TURNED O-F-F, I.E., THERE IS CONFLICT */
				/* CONFLICT SCENARIO 2 */
				if (span>1 && tela[n].stat!=st_parent.sym && tela[n].all_L && !(tela[n].all_R) && tela[n+1].ok < tela[n].ok &&
							!(tela[n+1].all_L) && !(tela[n+1].all_R) && tela[n].ok % tela[n+1].ok==0) {
					clearall_tela(n, 1, tela[n+1].all_S, TWO);		/* O-F-F, ONE, OR TWO */
					push_mem(n, 10);
					/* POSSIBLE THIS CASE COULD BE GENERALIZED...FOR A RAINY DAY */
					if (dev_print(TWO,__LINE__)) {
						printf("Mark_tela at n=%d, span=%d with left-conflict=%d, and no right_conflict, " 
								"and n+1 has smaller k that is a multiple of k-size at n with no conflicts.", n, span, tela[n].all_L);
					}
				}
				/* CONFLICT SCENARIO 3: FIRST TR CAN BE CYCLED BUT NOT THE SECOND, WHICH HAS HIGHER SCORE ANYWAYS */
				else if (span>1 && tela[(l=tela[n].all_R)].all_S > tela[n].all_S && !tela[l+1].all_S) {
					k = tela[n].ok;
					for (j=n+1; j<n+span; j++) {
						if (tela[j].ok == k && !tela[j].all_R && !tela[j].all_L) {
							if (dev_print(TWO,__LINE__)) {
								printf("Mark_tela at n=%d, scenario three", n);
							}
							for (i=j-1; i>=n; i--) {
								clearall_tela(i, j-n, -1, TWO);		/* O-F-F, ONE, OR TWO */
								for (p=i; p<=i+j-n; p++)
									push_mem(p, 10);
							}
							for (i=j+1; i<n+span; i++) {
								if (tela[i].ok == k && tela[i].all_R == l) {
									clearall_tela(i,1,-1, TWO);		/* O-F-F, ONE, OR TWO */
									push_mem(i, 11);
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
					if (tela[i].all_S == max_score && tela[i].ok == smallest_k) {
						max_count++;
						tela[i].all_Z = tela[i].all_S;
					}
				}
				if (max_count > 1) {
					j=1;	/* VAR j WILL COUNT ITERATIONS REQUIRED TO BREAK TIES */
					while (max_count > 1)
						max_count = settle_tiescores(n, span, max_score, j++);
				}
			}
			n = n + span - 1;
		}
	}

	/* CHECK TO SEE IF THERE ARE CYCLING DIFFERENCES BETWEEN INTERNAL? REPEATS */
	for (n=2; n+1<lenseq; n++) {
		if (tela[n].ok && !tela[n-1].ok && !tela[n+1].ok && tela[n].stat!=st_parent.sym) {
			k = tela[n].ok;
			m = n - k;

			if (tela[n].or==1 && tela[n].c==tela[n-1].c && tela[n].statf!=st_fract.sym) { 
				if (k==3 && !tela[m+1].ok && !tela[n+2].ok) {
				}
				else {
					tela[n].echoes = cyc_skip.sym;
					push_mem(n, 14);
				}
			}

			int p,q;
			for (i=0; i<tela[n].or; i++) {
				for (j=min_k; j<tela[n].ok-1; j++) {
					if (tela[(p=m+i*k+j)].ok != tela[(q=n+i*k+j)].ok && tela[n].statl != st_lowcm.sym && !tela[p].k2) {
						if (tela[p].ok && tela[p].statf!=st_fract.sym) {
							clearall_tela(p,1,-1, TWO);
							push_mem(p, 15);
						}
						if (tela[q].ok) {
							clearall_tela(q,1,-1, TWO);
							push_mem(q, 13);
						}
					}
				}
			}
		}
	}

	for (n=2; n<lenseq; n++) {
		if ((k=tela[n].k1) && tela[n].stat==st_parent.sym && tela[n].statf==st_fract.sym) {
			tela[n].stat = st_Fract.sym;
			i = 1;
			while (tela[n+i].ok==k)
				tela[n + i++].stat = st_Fract.sym;
		}
	}

	/* CLEAR OUT FRACTALS OF m UNIT OF PARENT k-MERS */
	for (n=lenseq; n>0; n--) {
		if ((k=tela[n].ok) && !tela[(m=n-k)].ok) {
			for (i=n-1; i>n-k; i--) {
				if (tela[i].ok && tela[i].ok<k && i-tela[i].ok>m && (tela[i].all_S<tela[n].all_S || tela[i].impk)) {
					clearall_tela(i,1,-1, TWO);
					push_mem(i, 17);
				}
			}
		}
	}

	/* FIND AND RECORD MAX_K SIZE */
	int max_k = 0;
	for (n=0; n<lenseq; n++) {
		if (tela[n].ok > max_k) {
			max_k = tela[n].ok;
		}
	}
	Cinch_T.pass_V = max_k;  

	if (opt_D.val==1 || dev_print(TWO,__LINE__)) {
		printf("\n\nPre-marking of tela completed as shown below. (To shift output window, change O-F-F-SET definition in 'src/microhom-devl.h'.)\n");
		print_tela(prtela_A, prtela_B);
	}
}


/*****************************************/
/* ORIGINAL IDEA FOR THIS FUNCTION WAS TO CALL TRs BY UNIQUE NUMBERED CASES AND THEN LATER CLEAR THEM (ERASE) IN AN EASILY MODIFIABLE ORDER */
/* CURRENTLY ONLY BEING USED FOR TRACKING HOW A TR WAS MARKED */
void push_mem(int pos, int row)
{
	if (row == 0) {
		tela[0].mem[0]++;			/* BIT SLOT TO INDICATE THIS ROW WAS WRITTEN TO */
		tela[pos].mem[0] = 18;		/* 18 for r-epeat, and convenient if non-zero so can find */
	}
	else if (row < MEMROWS) {
		tela[0].mem[row]++;			/* BIT SLOT TO INDICATE THIS ROW WAS WRITTEN TO */
		tela[pos].mem[row] = row;
	}
}

/*****************************/
short int check_solo(int pos)
{
	if (!tela[pos].ok) {
		return(-1);					/* NEGATIVE SIGNALS NO REPEAT AT GIVEN POSITION */
	}
	else if ( tela[pos-1].ok || 
			  tela[pos+1].ok ) {
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
int width = 56;	    /* 60 x 3 = 180 COLS, PRACTICAL DISPLAY WIDTH DEPENDS ON SCREEN SIZE */
int lenseq = Clean.pass_W;

	a += OFFSET;
	b += OFFSET;
	if (b>lenseq && lenseq > b-a) {
		a = lenseq-width;
		b = lenseq;
	}
	else if (width > lenseq) {
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

	if (b > lenseq)
		b = lenseq+1;

	/************* BEGIN PRINTING LINES *******************/
	if (Clean.pass_Q==1000) {	/* post-Mark_tela() */
		printf("\nxy:");
		for (i=a; i<b; i++)
			printf("%3d", tela[i].gPnt.rel_xy);
		printf("\ntP:");
		for (i=a; i<b; i++)
			printf("%3d", tela[i].gPnt.topPar);
		printf("\npP:");
		for (i=a; i<b; i++)
			printf("%3d", tela[i].gPnt.prevPar);
	
		printf("\n t:");
		for (i=a; i<b; i++) {
			if (tela[i].c != tela[i].t)
				printf("__%c", tela[i].t);
			else
				printf("___");
		}
		if (Clean.pass_V==3) { 	/* IF PROTEIN */
			printf("\n e:");
			for (i=a; i<b; i++)
				printf("  %c", tela[i].e);
		}
	}
	printf("\n c:");
	for (i=a; i<b; i++)
		printf("%3c", tela[i].c);

	printf("\n n:");
	if (b<100) {
		for (i=a; i<b; i++)
			printf("%3d", i);
	}
	else {
		for (i=a; i<b; i++) {
			if (i<b-1) {
				if (i%10)
					printf("   ");
				else
					printf("  |");
			}
			else
				printf("  %d", i);
		}
	}

	if (Clean.pass_Q==1000) {
		printf("\nLf:");
		for (i=a; i<b; i++) {
			if (tela[i].cyc_Lf)
				printf("%3d", tela[i].cyc_Lf);
			else
				printf("  .");
		}
	
		printf("\nRt:");
		for (i=a; i<b; i++) {
			if (tela[i].cyc_Rt)
				printf("%3d", tela[i].cyc_Rt);
			else
				printf("  .");
		}
		printf("\n X:");
		for (i=a; i<b; i++) {
			if (tela[i].X != i)
				printf("%3d", tela[i].X);
			else
				printf("  .");
		}
		printf("\n y:");
		for (i=a; i<b; i++)
			printf("%3d", tela[i].y);
	
		printf("\n x:");
		for (i=a; i<b; i++)
			printf("%3d", tela[i].x);
	}

	printf("\nst:");
	for (i=a; i<b; i++) {
		if (tela[i].stat)
			printf("%3c", tela[i].stat);
		else
			printf("  .");
	}
	printf("\nfr:");
	for (i=a; i<b; i++) {
		if (tela[i].statf)
			printf("%3c", tela[i].statf);
		else
			printf("  .");
	}
	printf("\nlc:");
	for (i=a; i<b; i++) {
		if (tela[i].statl)
			printf("%3c", tela[i].statl);
		else
			printf("  .");
	}
	printf("\nk0:");
	for (i=a; i<b; i++) {
		if (tela[i].k0)
			printf("%3d", tela[i].k0);
		else
			printf("  .");
	}
	printf("\nk1:");
	for (i=a; i<b; i++) {
		if (tela[i].k1)
			printf("%3d", tela[i].k1);
		else
			printf("  .");
	}
	printf("\nk2:");
	for (i=a; i<b; i++) {
		if (tela[i].k2)
			printf("%3d", tela[i].k2);
		else
			printf("  .");
	}
	printf("\nik:");
	for (i=a; i<b; i++) {
		if (tela[i].impk<0)
			printf("%3d", tela[i].impk);
		else
			printf("  .");
	}
	printf("\nok:");
	for (i=a; i<b; i++) {
		if (tela[i].ok)
			printf("%3d", tela[i].ok);
		else
			printf("  .");
	}
	printf("\nor:");
	for (i=a; i<b; i++) {
		if (tela[i].ok)					/* THIS IS MEANT TO BE ALL_K */
			printf("%3d", tela[i].or);
		else
			printf("  .");
	}
	printf("\naS:");
	for (i=a; i<b; i++) {
		if (tela[i].all_S)
			printf("%3d", tela[i].all_S);
		else
			printf(" __");
	}
	printf("\naZ:");
	for (i=a; i<b; i++) {
		if (tela[i].all_Z)
			printf("%3d", tela[i].all_Z);
		else
			printf(" __");
	}
	printf("\naR:");
	for (i=a; i<b; i++) {
		if (tela[i].all_R)
			printf("%3d", tela[i].all_R);
		else
			printf("  .");
	}
	printf("\naL:");
	for (i=a; i<b; i++) {
		if (tela[i].all_L)
			printf("%3d", tela[i].all_L);
		else
			printf(" __");
	}

	if (Clean.pass_Q==1000) {
		printf("\n k:");
		for (i=a; i<b; i++) {
			if (tela[i].k)
				printf("%3d", tela[i].k);
			else
				printf("  .");
		}
		printf("\n r:");
		for (i=a; i<b; i++) {
			if (tela[i].r)
				printf("%3d", tela[i].r);
			else
				printf(" __");
		}
	
		printf("\nDt:");
		for (i=a; i<b; i++) {
			if (tela[i].Dtr)
				printf("%3d", tela[i].Dtr);
			else
				printf("  .");
		}
	
		printf("\n o:");
		for (i=a; i<b; i++) {
			if (tela[i].o)
				printf("%3d", tela[i].o);
			else
				printf("  .");
		}
	}
	printf("\n E:");
	for (i=a; i<b; i++)
		printf("%3c", tela[i].echoes);

	printf("\nisl");
		for (i=a; i<b; i++)
			if (tela[i].isl) {
				if (tela[i].isl<10)
					printf("  %d", tela[i].isl);
				else
					printf("  %c", mha_base62(tela[i].isl));
			}
			else
				printf(" ~~");
	
	if (Clean.pass_Q==1000) {
		printf("\ncO:");
		for (i=a; i<b; i++)
			printf("%3c", tela[i].cyc_o);
	
		printf("\ncL:");
		for (i=a; i<b; i++) {
			if (tela[i].cyc_l) 
				printf("%3d", tela[i].cyc_l);
			else
				printf("  .");
		}
	
		printf("\ncK:");
		for (i=a; i<b; i++) {
			if (tela[i].cyc_k)
				printf("%3d", tela[i].cyc_k);
			else
				printf("  .");
		}
		printf("\ncR:");
		for (i=a; i<b; i++) {
			if (tela[i].cyc_r)
				printf("%3d", tela[i].cyc_r);
			else
				printf("  .");
		}

		printf("\ncP:");
		for (i=a; i<b; i++) {
			if (tela[i].cyc_P)
				printf("%3d", tela[i].cyc_P);
			else
				printf("  .");
		}
	
		printf("\ncS:");
		for (i=a; i<b; i++) {
			if (tela[i].cyc_S)
				printf("%3d", tela[i].cyc_S);
			else
				printf("  .");
		}
	}
	printf("\n   ");
	for (i=a; i<b; i++)
		printf(" --");

	/* PRINT TOP MEM ROWS: Mark_tela mem OR frame rows in cinch-t */
	for (f=0; f<MEMROWS; f++) {
		if (tela[0].mem[f]) {
			printf("\nm%c:", mha_base62(f));
			for (i=a; i<b; i++) {
				if (tela[i].mem[f])
					printf("%3d", tela[i].mem[f]);
				else
					printf("  .");
			}
		}
	}

	printf("\n");
}


/* ERASE REPEAT AT n, INCLUDING ANY TRANSITIONS, LEAVES k MARK ****/
void pull_tela(int n)
{
	int i=0;
	int k = tela[n].ok;
	int r = tela[n].or;
	int m = n-k;

	/* IF IMPERFECT REPEAT, ERASE. THIS IS MOSTLY RIGHT, BUT EVENTUALLY NEED TO CHECK OVERLAPPING TRANSITIONS */
	if (k*r*MATCH > tela[n].all_S) {	
		for (i=m; i<(n+k*r); i++) {
			consensus[(tela[i].x)] = tela[i].t = tela[i].c;
		}
	}

	tela[n].r = 0;
}


/**** FUNCTION TO TOKENIZE A 2-D ROW IF THERE ARE NO VIOLATIONS BASED ON TESTS OF AXIOMS 1 OR 2, NONE (0), OR BOTH (3) *****/
int push_tela(int n2, int n1, short unsigned int axioms) 
{
	int coord_A_x=0, coord_B_x=0;
	int coord_A_y=0, coord_B_y=0;
	int i = 0;
	int k = n2 - n1;
	int lenseq=Clean.pass_W;
	int violation=0;

	if (axioms==1 || axioms==3) {
		/* CHECK PRINCIPLE OF CONTINUITY */
		for (i=0; i<k; i++) {
			coord_A_x = tela[n1+i    ].x;
			coord_A_y = tela[n1+i    ].y;

			coord_B_x = tela[n1+i + 1].x;
			coord_B_y = tela[n1+i + 1].y;

			if      (coord_B_y==coord_A_y     && coord_B_x==coord_A_x + 1)
				;
			else if (coord_B_y==coord_A_y + 1 && coord_B_x<=coord_A_x    )
				;	
			else {
				violation = 1;
				sprintf(dev_notes, "viol1 tela-%d", __LINE__);
				if (dev_print(TWO,__LINE__)) {
					printf("push_tela() viol-1 for n2=%d, n1=%d, k=%d (i=%d): A_x=%d, A_y=%d, B_x=%d, B_y=%d. axioms=%d",
												n2, n1, k, i, coord_A_x, coord_A_y, coord_B_x, coord_B_y, axioms);
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
			else if (tela[n1+i].e == tela[n2+i].e) {	/* THIS IS REDUNDANT IF ].e IS EVALUATED BUT HERE FOR TESTING PURPOSES */
				;
			}
			else {
				violation += 2;
				sprintf(dev_notes, "viol2 tela-%d", __LINE__);
				if (dev_print(TWO,__LINE__)) {
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
		/* UPDATE TELA STRUCTURE IF AND ONLY IF THERE ARE NO VIOLATIONS */
			tela[n2].x = tela[n2-1].x-k+1;
			tela[n2].y = tela[n2-1].y +1;
		for (i=n2+1; i<=lenseq; i++) {
			tela[i].x -= k;
			tela[i].y = tela[n2].y;
		}
	}

	return(violation);	/* RETURNS 1 IF CONTINUITY FAILS, 3 IF EQUIVALENCE FAILS, 4 IF BOTH FAIL */
}

/****************************************************/
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

		if (mode == 2) { /* SCORE TRANSITIONS */
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
	int lenseq = Clean.pass_W;
	int max_count=0, ratchet=0, ratchet_up=0, ratchet_dn=0;
	unsigned short int nuctype = Clean.pass_V;
	unsigned short int nuctransit = 0;
	int up_score=0, dn_score=0;

	if (nuctype == 1) {		/* IF DNA */
		nuctransit = 1;
	}

	for (i=n; i<n+span; i++) {
		if (tela[i].all_S == max_score) {
			ratchet_up = ratchet_dn = up_score = dn_score = 0;
			k = tela[i].ok;
			m = i-k;						/* DEFINES START OF FIRST REPEAT UNIT */
			o = i+k*(tela[i].or-1);			/* DEFINES START OF LAST REPEAT UNIT */
			up = m-1 - k*(iteration-1);		/* DEFINES THE GHOST FLANKING UNIT STARTING AT m-1, ADJUSTED FOR MOVING LEFTWARDS */
			dn = o + k*iteration;			/* DEFINES THE GHOST FLANKING UNIT STARTING AFTER REPEATS */
			for (j=0; j<k; j++) {
				if (ratchet_up >= 0) {
					if (up-j<0 || (nuctransit && tela[up-j].c == ambig.sym))
						ratchet_up = -1;
					else if (tela[up-j].c == tela[up-j+k].c) {
						up_score += match;
						ratchet_up++;
						ratchet++;
					}
					else if (tela[up-j].e == tela[up-j+k].e){
						up_score += transition;
						ratchet_up++;
						ratchet++;
					}
					else {
						ratchet_up++;
						ratchet++;
					}
				}
				if (ratchet_dn >= 0) {
					if (dn+j > lenseq || (nuctransit && tela[dn+j].c == ambig.sym))
						ratchet_dn = -1;
					else if (tela[dn+j].c == tela[dn+j-k].c) {
						dn_score += match;
						ratchet_dn++;
						ratchet++;
					}
					else if (tela[dn+j].e == tela[dn+j-k].e) {
						dn_score += transition;
						ratchet_dn++;
						ratchet++;
					}
					else {
						ratchet_dn++;
						ratchet++;
					}
				}
			}
			tela[i].all_Z += up_score;
			tela[i].all_Z += dn_score;
		}
	}
	if (!ratchet)		/* IF NO RATCHETING, THEN GHOST ITERATIONS ARE OUT-OF-BOUNDS AND INEFFECTIVE TIE-BREAKERS */
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
	if (max_count == 1) {				/* THERE IS ONE OPTIMAL CINCH LOCATION AND NO CONFLICT SO ERASE ALL OTHERS */
		clearall_tela(n, span, max_score, ONE);		/* O-F-F, ONE, OR TWO */
		for (i=n; i<n+span; i++) {
			if (tela[i].all_Z != max_score)
				push_mem(i, 17);
		}
	}

	return(max_count);
}


/* CALL ONLY IF align2D is 100% ***********************************/
int update_tela(void)
{
	int c=0, i=0, j=0;
	int lenseq=Clean.pass_W;
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
				if (letr==Term->sym) {
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
