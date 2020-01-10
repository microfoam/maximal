/******************************************************************/
/* microhom-tela.h header file, since mha_v4.11.c                 */
/* This file has MHA functions either in the form of verb_tela(), */
/* or verb_object(), where object is related to tela struct.      */
/******************************************************************/


int assign_tela(int eL, int eM, int eN, int mode, int pointA, int pointB)
{
	int i=0, j=0, l=0, conflict_flag=0;
	int lenseq = options[1][1];
	int start=0, end=lenseq;

	if (pointB > pointA) {	/* ASSIGN LINE VALUES IF CHECKING LIMITED SPAN */
		start = pointA;
		end   = pointB;
	}
	else {
		start = 0;
		end   = lenseq;
	}

	if (!mode) {					/* MODE 0: ASSIGN SAME COLUMN TO IDENTICAL LETTERS */
		align2D[eM][eN] = tela[eL].c;
		tela[eL].y = eM;
		tela[eL].x = eN;
	}
	else if (mode==1) {				/* MODE 1: CHECK AND ASSIGN SAME COLUMN TO IDENTICAL LETTERS */
		for (i=start; i<end; i++) {
			if (tela[i].x == eN && tela[i].c != tela[eN].c) {
				conflict_flag = 1;
				break;
			}
		}
		if (!conflict_flag) {
			align2D[eM][eN] = tela[eL].c;
			tela[eL].y = eM;
			tela[eL].x = eN;
		}
	}
	else if (mode==2) {				/* MODE 2: ASSIGN SAME COLUMN TO EQUIVALENT LETTERS */
		align2D[eM][eN] = tela[eL].c;
		tela[eL].y = eM;
		tela[eL].x = eN;
		for (i=start; i<=end; i++) {
			if (tela[i].x == eN && tela[i].c != tela[eL].c && tela[i].e == tela[eL].e) {
				j = tela[i].x;
				consensus[j] = tela[i].t = tela[i].e;
			}
		}	
	}
	else if (mode==3) {				/* MODE 3: CHECK AND ASSIGN SAME COLUMN TO EQUIVALENT LETTERS */
		for (i=start; i<end; i++) {
			if (tela[i].x == eN && tela[i].c != tela[eL].c && tela[i].e != tela[eL].e) {
				conflict_flag = 1;
				break;
			}
		}

		if (!conflict_flag) {
			align2D[eM][eN] = tela[eL].c;
			tela[eL].y = eM;
			tela[eL].x = eN;
			for (i=start; i<=end; i++) {
				if (tela[i].x == eN && tela[i].c != tela[eL].c && tela[i].e == tela[eL].e) {
					j = tela[i].x;
					consensus[j] = tela[i].t = tela[i].e;
				}
			}	
		}
	}
	else if (mode==4) {			/* MODE 4: FLAT-LINE TELA STARTING AT POINT eL */
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
			for (l = eL; l <= lenseq; l++) { 
				tela[l].y = i;
				tela[l].x = j++;
			}
		}
	}
	else {
		warnhead('m');
		printf("\ntela-DEV-0093: undefined mode invoked");
	}

	if (conflict_flag)
		return(0);
	else
		return(1);
}


/****** TELA: A FABRIC, UNDER AXIOMATIC LAWS **********************/
int check_tela(int eM, int eN, short unsigned int dim) 
{
int i=0, j=0, lineM=0, lineN=0, princeps=0, badflag=0;
int lenseq = options[1][1];

	if (eM>=eN) {
		printf("\ntela-DEV-0110: Need to call check_tela explicitly with %d-D positions eN > eM.", dim);
		princeps-=5;
	}

	if (dim==1) {
		/* GET 2-D COORDINATES FROM 1-D COORDINATES */
		lineM = tela[eM].x;	
		lineN = tela[eN].x;	
	}
	else if (dim==2) {
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
		printf("\ntela-DEV-0134: Need to call check_tela explicitly with dimension dim=1 or dim=2.");
		princeps-=7;
	}

	if (princeps==0) {
		/* ANGEL ONE: THE PRINCE OF CONTINUITY */
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
			princeps = 1;
		else
			printf("\ntela-DEV-0157: check_tela(): Problem of continuity at 1-D positions %d --> %d (columns %d and %d)", i, i+1, tela[i].x, tela[i+1].x);

		/* ANGEL TWO: THE PRINCE OF EQUIVALENCE */
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
			princeps+=2;
		else
			printf("\ntela-DEV-0175: check_tela(): Problem of equivalence at 1-D positions %d and %d (both in column %d)", i, j, tela[i].x);
	}

	return(princeps);	/* 0 IF BOTH FAIL; +1 IF ONLY ONE PASSES; +2 IF ONLY TWO PASSES; +3 IF BOTH PASS */ 
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


/* PROPAGATE NEW TELA COORDINATES TO THE RIGHT FROM GIVEN START ***/
void propagate_tela(int start, int row, int a2D_n) 
{
	int i=0;
	int lenseq = options[1][1];

	for (i=start; i<=lenseq; i++) {
		tela[i].y = row;
		tela[i].x = a2D_n++;
	}
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
				printf("\ntela-DEV-0521: push_tela() viol-1 for n2=%d k=%d (i=%d): A_x=%d, A_y=%d, B_x=%d, B_y=%d.", 
														n2, k, i, coord_A_x, coord_A_y, coord_B_x, coord_B_y); 
				break;
			}
		}
	}

	if (axioms==2 || axioms==3) {
		/* CHECK PRINCIPLE OF EQUIVALENCE */
		for (i=0; i<k; i++) {
			if (0 && tela[n1+i].t != tela[n2+i].t && tela[n1+i].c == tela[n2+i].c) {
				tela[n1+i].t = tela[n1+i].e;
				tela[n2+i].t = tela[n2+i].e;
			}
			else if (tela[n1+i].c == tela[n2+i].c)
				;
			else if (tela[n1+i].t == tela[n2+i].t) 
				;
			else if (0 && tela[n1+i].e == tela[n2+i].e) { 
				tela[n1+i].t = tela[n1+i].e;
				tela[n2+i].t = tela[n2+i].e;
			}
			else {
				violation += 2;
				printf("\ntela-DEV-0543: push_tela() viol-2 at %d and %d for k=%d.", n1+i, n2+i, k);
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
		for (i=n2+k ; i<=lenseq; i++) {
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


