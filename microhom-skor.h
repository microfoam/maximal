/******************************************************************/
/* microhom-skor.h header file, since mha_v4.12.c                 */
/* This file has MHA functions and definions related to scoring.  */
/*                                                                */
/******************************************************************/

#define MATCH       8       /* MATCH SCORE */
#define MISMATCH   -1		/* MISMATCH SCORE */
#define TRANSITION	4		/* TRANSITION = HALF MATCH SCORE */
#define PISO        5       /* FLOOR FOR TRANSITION MATCHING ABOVE THIS k-MER SIZE */

int 	score_DTHR(int kmer);
void 	show_DTHR_table(void);

/***********************/
int score_DTHR(int kmer)
{	/* RETURNS SCORE THRESHOLDS FOR K-MER FOR ALLOWED TRANSITIONS) */
	int k=0, numtransit = 0;
	float fractransit = 0.08;	/* SETS NUMBER OF ADDITIONAL ALLOWED TRANSITIONS FOR GIVEN k-MER */
	int transition = TRANSITION; 
	static int thr_table[WIDTH+1];
	int squeeze = (int) options[1][59];

	for (k=0; k<=PISO; k++)
		thr_table[k] = 100;

    for (k = PISO+1; k <= WIDTH; k++) {
        /* FUNCTION ONE: LOGARITHMIC: SLOW, ORGANIC INCREASE IN ALLOWED TRANSITIONS; ALSO MANIFESTS PISO (FLOOR) NATURALLY */
        if (ONE) {
			numtransit = (int) round(0.32*log2( (float) k ));		/* EQUIVALENT TO 0.64 * LOG_BASE_4 (k) */
       	} 
 
        else {
            numtransit = 1 + round(fractransit * k);    /* FRACTRANSIT 0.08 ALLOWS k<7 KMERS ONLY ONE TRANSITION */
        	/* FUNCTION TWO: FLAT FRACTION OF TRANSITIONS: VERY HAIRY, BUT GREAT FOR DEBUGGING CRAZY KNOTS       */
        	/* ---THIS IS A LEGACY DEVELOPMENT-PHASE FUNCTION SAVED FOR DIDACTIC & HISTORICAL PURPOSES.          */
        }    

        thr_table[k] = 100*((k-numtransit)*MATCH + numtransit*transition)/(k*MATCH) - squeeze;
    }    

	return(thr_table[kmer]);   
}
/***********************/


/*************************/
void show_DTHR_table(void)
{	
	int k, maxtransits, max_score, table_score;
	int match = MATCH;
	int transition = TRANSITION;
	short unsigned int seqtype = options[1][13];

	printf("\n Diagonal thresholds as a function of k (must exceed threshold):\n");

	printf("\n k-mer\t Max. trans.\t Threshold\t Score with maximum transitions"); 
	for (k = 1; k <= WIDTH; k++) {
		if (seqtype==1 && k>PISO) {
			maxtransits = (int) round(0.32*log2((float) k));
			table_score = score_DTHR(k);
			max_score = (int) 100*((k-maxtransits)*match + maxtransits*transition)/(k*match);
		}
		else {
			maxtransits = 0;
			table_score = max_score = 100;
		}
		printf("\n %3d\t %d\t %12d\t\t %d", k, maxtransits, table_score, max_score);
		if (max_score>table_score)
			printf("\t* Above threshold");
		else
			printf("\t (< threshold)");
	}
	printf("\n\n Note: DTHR values are only populated if a sequence is specified.\n\n");
	exit(EXIT_EARLY);
}
/*************************/

