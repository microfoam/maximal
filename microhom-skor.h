/******************************************************************/
/* microhom-skor.h header file, since mha_v4.12.c                 */
/* This file has MHA functions and definions related to scoring.  */
/*                                                                */
/******************************************************************/

#define MATCH       8       /* MATCH SCORE */
#define TRANSITION	4		/* TRANSITION = HALF MATCH SCORE */

int score_DTHR(int kmer);

/**** FUNCTION TO RETURN SCORE THRESHOLDS FOR K-MER (NUMBER OF ALLOWED TRANSITIONS) *****/
int score_DTHR(int kmer)
{
	int k=0, numtransit = 0;
	int transition = TRANSITION; 
	float fractransit = 0.08;		/* USED TO SET NUMBER OF ADDITIONAL ALLOWED TRANSITIONS FOR GIVEN k-MER FOR FUNC. TWO */
	static int thr_table[WIDTH+1];
	int squeeze = options[1][59];

	for (k=0; k<=PISO; k++)
		thr_table[k] = 100;

    for (k = PISO+1; k <= WIDTH; k++) {
        /* FUNCTION ONE: LOGARITHMIC: SLOW, ORGANIC INCREASE IN ALLOWED TRANSITIONS; ALSO MANIFESTS PISO (FLOOR) NATURALLY */
        if (1) {
			numtransit = round(0.32*log2( (float) k ));		/* EQUIVALENT TO 0.64 * LOG_BASE_4 (k) */
        	}                                               /* ONE TRANSITION ALLOWED STARTING AT k=4 OR k=3 WITH -x SQUEEZE */

        /* FUNCTION TWO: FLAT FRACTION OF TRANSITIONS: VERY HAIRY, BUT GREAT FOR DEBUGGING CRAZY KNOTS (EVO. GEN. PROGRAMMING) */
        /* ---THIS IS A LEGACY DEVELOPMENT-PHASE FUNCTION SAVED FOR DIDACTIC & HISTORICAL PURPOSES.                            */
        else {
            numtransit = 1 + round(fractransit * k);    /* FRACTRANSIT 0.08 ALLOWS k<7 KMERS ONLY ONE TRANSITION */
        }    

        thr_table[k] = 100*((k-numtransit)*MATCH + numtransit*transition)/(k*MATCH) - squeeze;
    }    

	return(thr_table[kmer]);   
}

