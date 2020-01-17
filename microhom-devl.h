/******************************************************************/
/* microhom-devl.h header file, since mha_v4.23.c                 */
/* This file has MHA functions related to program development,    */
/* testing, and evaluation.                                       */
/******************************************************************/

#define DEVBIT1		0		/* TURN ON dev_prompt()'s IF ON */
#define OFF			0		/* FOR CLARITY OF MODE ARGUMENTS */
#define ON 			1		/* FOR CLARITY OF MODE ARGUMENTS */
#define ONE			1		/* USE ONLY FOR FUNCTIONS THAT CAN BE TURNED OFF WITH MODE OFF */
#define TWO			2		/* USE ONLY FOR FUNCTIONS THAT CAN BE TURNED OFF WITH MODE OFF */
#define THREE		3		/* USE ONLY FOR FUNCTIONS THAT CAN BE TURNED OFF WITH MODE OFF */
#define FOUR		4		/* USE ONLY FOR FUNCTIONS THAT CAN BE TURNED OFF WITH MODE OFF */
#define FIVE		5		/* USE ONLY FOR FUNCTIONS THAT CAN BE TURNED OFF WITH MODE OFF */
#define MAIN		1		/* FOR USE WITH dev_print() */
#define TELA		2		/* FOR USE WITH dev_print() */

void 				dev_linehead(int mode, int line_no);
short unsigned int	dev_print(short unsigned int mode, int line_no);
void 				dev_prompt(short unsigned int mode, int line_no, char *filename);
void 				signal_callback_handler(int signum);

int prtela_A =   0;		/* DEVELOPMENT VARIABLE TO GLOBALLY SET PRINT TELA START AND END POINT FOR DIAGNOSTICS */
int prtela_B =  56;		/* SET TO SOMETHING APPROPRIATE FOR GIVEN COMPUTER SCREEN */	

/*****************************************************************/
short unsigned int dev_print(short unsigned int mode, int line_no)
{
	if (!mode) {
		/***  CODED O-F-F  ***/	
		return(0);
	}
	/***  DEV-LEVEL VERBOSITY  ***/	
	else if (options[1][57] > 2) {		/* opt_v: reserve 1-2 for user level verbosity */
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
	else if (DEVBIT1) {
		dev_linehead(mode, line_no);		
		printf("Press <enter> to continue with %s. ", filename);
		getchar();
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

	printf("\nDEV-%s-%4d: ", section_label, line_no);

}
/***************************************/


/***************************************/
void signal_callback_handler(int signum) 
{
	printf("  )--- I caught signal %d before exiting (2=SIGINT, 10=SIGBUS, 11=SIGSEGV).\n\n",signum);
	fp_out = fopen("Surf_wavereport.mha", "a");		/* FOPEN RIGHT BEFORE WRITING TO MINIMIZE CHANCE OF CLOSING WITH OPEN FILES */
	fprintf(fp_out, "---->\tCanceled run for %s with signal=%d (2=SIGINT, 10=SIGBUS, 11=SIGSEGV). dev_notes: %s.\n", file_name, signum, dev_notes);
	fclose(fp_out);
	exit(signum);
}
/***************************************/


