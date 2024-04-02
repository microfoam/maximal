/******************************************************************/
/* microhom-devl.h header file, since mha_v4.23.c                 */
/* This file has MHA functions related to program development,    */
/* testing, error emission, and evaluation.                       */
/******************************************************************/

#ifndef FILE_DEVL_SEEN
#define FILE_DEVL_SEEN

#define DEVX		1		/* TO COORDINATE TESTING OF NEW CODE DISTRIBUTED THROUGHOUT FILES */
#define OFFSET		0		/* FOR SETTING START POSITION OF PRINT TELA */
#define OFF			0		/* FOR CLARITY OF MODE ARGUMENTS */
#define SKIP		0		/* FOR CLARITY OF MODE ARGUMENTS; DEV-USE */
#define ON 			1		/* FOR CLARITY OF MODE ARGUMENTS */
#define ONE			1		/* USE ONLY FOR FUNCTIONS THAT CAN BE TURNED OFF WITH MODE OFF */
#define TWO			2		/* USE ONLY FOR FUNCTIONS THAT CAN BE TURNED OFF WITH MODE OFF */
#define THREE		3		/* USE ONLY FOR FUNCTIONS THAT CAN BE TURNED OFF WITH MODE OFF */
#define MAIN		1		/* FOR USE WITH dev_print() */
#define CINCH		3		/* FOR USE WITH dev_print() */
#define LOGY		4		/* FOR USE WITH dev_print() */
#define EXIT_GOOD	0		/* FOR STANDARD EXIT ERRORS */
#define EXIT_ERROR	1		/* FOR STANDARD EXIT ERRORS */
#define EXIT_EARLY	2		/* FOR STANDARD EXIT ERRORS */

void 				dev_linehead(int mode, int line_no);
short unsigned int	dev_print(short unsigned int mode, int line_no);
void 				dev_prompt(short unsigned int mode, int line_no, char *filename);
void 				signal_callback_handler(int signum);

int prtela_A =   0;		/* DEVELOPMENT VARIABLE TO GLOBALLY SET PRINT TELA START AND END POINT FOR DIAGNOSTICS */
int prtela_B =  47;		/* SET TO SOMETHING APPROPRIATE FOR GIVEN COMPUTER SCREEN */

/*****************************************************************/
short unsigned int dev_print(short unsigned int mode, int line_no)
{
	if (!mode) {
		/***  CODED O-F-F  ***/	
		return(0);
	}
	/***  DEV-LEVEL VERBOSITY  ***/	
	else if (opt_v.val > 2) {		/* opt_v: reserve 1-2 for user level verbosity */
		dev_linehead(mode, line_no);		
		return(1);
	}
	else {
		/***  DEV-VERBOSITY NOT OPTIONED  ***/	
		return(0);
	}
}
/*****************************************************************/

/****************************************************/
void dev_prompt(short unsigned int mode, int line_no, char *filename)
{
	if (!mode) {
		/***  CODED O-F-F  ***/	
		return;
	}
	else if (opt_D.bit) {	/* IF opt_D DEV_PROMPTS ON */
		if (opt_v.val>1)
			dev_linehead(mode, line_no);		
		printf("\n./maximal: Hit <enter> to continue run of %s (%c), or enter 'q' to quit this run.\n", filename, Strand->sym);
		char key=getchar();
		if (key=='q' || key=='Q')
			exit(17);
		return;
	}
}
/****************************************************/

/***************************************/
void dev_linehead(int mode, int line_no) 
{
	char section_label[8] = "unknown";
	if (mode == 1)
		strcpy(section_label, "main");
	else if (mode == 2)
		strcpy(section_label, "tela");
	else if (mode == 3)
		strcpy(section_label, "cinc");
	else if (mode == 4)
		strcpy(section_label, "logy");

	printf("\nDEV-%s-%4d: ", section_label, line_no);

}
/***************************************/

/***************************************/
void signal_callback_handler(int signum) 
{
	printf("  )--- Signal %d caught (2=SIGINT, 8=SIGFPE, 10=SIGBUS, 11=SIGSEGV).\n\n",signum);
	fp_out = fopen("Surf_wavereport.mha", "a");
	fprintf(fp_out, "---->\tCanceled run for %s with signal=%d (2=SIGINT, 8=SIGFPE, 10=SIGBUS, 11=SIGSEGV). Current.pass_V=%d, %s.\n", file_name, signum, Current.pass_V, dev_notes);
	fclose(fp_out);
	exit(signum);
}
/***************************************/

#endif		/* !FILE_DEVL_SEEN */

/*************************************************************************************************************/
/* WRITTEN BY DR. ALBERT J. ERIVES, AGPL-3.0 license. Code available at https://github.com/microfoam/maximal */
/*************************************************************************************************************/
