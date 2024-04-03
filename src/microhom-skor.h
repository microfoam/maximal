/******************************************************************/
/* microhom-skor.h header file, since mha_v4.12.c                 */
/* This file has MHA functions and definions related to scoring.  */
/******************************************************************/

#ifndef FILE_SKOR_SEEN
#define FILE_SKOR_SEEN

#define MATCH       8		/* MATCH SCORE */
#define MISMATCH   -1		/* MISMATCH SCORE */
#define TRANSITION  4		/* TRANSITION = HALF MATCH SCORE */

int 	allowed_transits(int k);
int 	score_DTHR(int kmer);
int 	score_transits(int k, int numtransits);
void 	show_DTHR_table(void);

int thr_table[WIDTH+1] = {0};


/*******************************/
int allowed_transits(int k)
{
	int numtransit = 0;	

	if (ON) {
		numtransit = (int) round(0.32*log2( (float) k ));		/* EQUIVALENT TO 0.64 * LOG_BASE_4 (k) */
		/* FUNCTION ONE: LOGARITHMIC: SLOW, ORGANIC INCREASE IN ALLOWED TRANSITIONS; MANIFESTS TRANSITIONS-FLOOR NATURALLY */
	} 
	else {
		float fractransit = 0.08;	/* SETS NUMBER OF ADDITIONAL ALLOWED TRANSITIONS FOR GIVEN k-MER */
		numtransit = 1 + round(fractransit * k);    /* FRACTRANSIT 0.08 ALLOWS k<7 KMERS ONLY ONE TRANSITION */
		/* FUNCTION TWO: FLAT FRACTION OF TRANSITIONS: VERY HAIRY, BUT GREAT FOR DEBUGGING CRAZY KNOTS       */
		/* ---THIS IS A LEGACY DEVELOPMENT-PHASE FUNCTION SAVED FOR DIDACTIC & HISTORICAL PURPOSES.          */
	} 
	return(numtransit);
}
/*******************************/


/************************/
int score_DTHR(int kmer)
{	/* RETURNS SCORE THRESHOLDS FOR K-MER FOR ALLOWED TRANSITIONS) */
	int k=0, numtransit = 0, threshold = 0;

	if (!thr_table[0]) {
		for (k=0; k<=opt_b.val; k++)
			thr_table[k] = 100;
	
	    for (k = opt_b.val+1; k <= WIDTH; k++) {
			numtransit = allowed_transits(k);
	        thr_table[k] = 100*((k-numtransit)*MATCH + numtransit*TRANSITION)/(k*MATCH) - opt_x.bit;
	    }
	}    
	threshold = thr_table[kmer];	
	return(threshold);   
}
/************************/


/********************************************/
int score_transits(int k, int numtransits)
{
	int score;

	score = 100*((k-numtransits)*MATCH + numtransits*TRANSITION)/(k*MATCH);
	return(score);
}
/********************************************/


/**************************/
void show_DTHR_table(void)
{	
	int k, maxtransits, max_score, table_score;
	int match = MATCH;
	int transition = TRANSITION;
	short unsigned int seqtype = Clean.pass_V;

	printf("\n Diagonal thresholds as a function of k (must exceed threshold):\n");

	printf("\n k-mer\t Max. trans.\t Threshold\t Score with maximum transitions"); 
	for (k = 1; k <= WIDTH; k++) {
		if (seqtype==1 && k>opt_b.val) {
			maxtransits = allowed_transits(k);
			table_score = score_DTHR(k);
			max_score = (int) 100*((k-maxtransits)*match + maxtransits*transition)/(k*match);
		}
		else {
			maxtransits = 0;
			table_score = max_score = 100;
		}
		printf("\n %3d\t %d\t %12d\t\t %9d", k, maxtransits, table_score, max_score);
		if (max_score>table_score) {
			printf("\t* Above threshold");
		}
		else {
			printf("\t (< threshold)");
		}
	}
	printf("\n\n");
	exit(EXIT_EARLY);
}
/**************************/


#endif		/* !FILE_SKOR_SEEN */

/*************************************************************************************************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal */
/*************************************************************************************************************/
