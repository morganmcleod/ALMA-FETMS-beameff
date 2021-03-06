#include <stdio.h>
#include <string.h>
#include <math.h>
#include "iniparser.h"
#include "utilities.h"
extern int DEBUGGING;

char warning[WARNING_LENGTH];  /* this is global, so rename if conflict */ 
char printmessage[MESSAGE_LENGTH];  //* this is global, so rename if conflict */ 
char stdOutDirectory[MESSAGE_LENGTH];

void warn(void) {
  static FILE *warnfile = NULL;
  char fn[MESSAGE_LENGTH];
  if (warnfile == NULL) {
      strcpy(fn, stdOutDirectory);
      strcat(fn, "/stderr.txt");
      warnfile = fopen(fn, "a");
      printf("warnfile=%s\n", fn);
  }
  printf(warning); /* to the stdout, i.e. to the screen */
  if (warnfile) {
      fprintf(warnfile, warning);  /* and now to the file */
      fflush(warnfile);  /* make sure the file gets updated */
  }
} 

void print_stdout(void) {
  static FILE *stdout_file = NULL;
  char fn[MESSAGE_LENGTH];
  if (stdout_file == NULL) {
      strcpy(fn, stdOutDirectory);
      strcat(fn, "/stdoutput.txt");
      stdout_file = fopen(fn, "a");
      printf("stdout_fileme=%s\n", fn);
  }
  printf(printmessage); /* to the stdout, i.e. to the screen */
  if (stdout_file) {
      fprintf(stdout_file, printmessage);  /* and now to the file */
      fflush(stdout_file);  /* make sure the file gets updated */
  }
} 


#ifdef LINUX
/* the following 2 functions are from Kernighan and Ritchie */

/* kreverse:  reverse string s in place */
void kreverse(char s[])
{
    int i, j;
    char c;

    for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

/* itoa:  convert n to characters in s */
char *itoa(int n, char s[], int unusedLength)
{
    int i, sign;

    if ((sign = n) < 0)  /* record sign */
        n = -n;          /* make n positive */
    i = 0;
    do {       /* generate digits in reverse order */
        s[i++] = n % 10 + '0';   /* get next digit */
    } while ((n /= 10) > 0);     /* delete it */
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    kreverse(s);
    return(s);
} 
#endif

void GetScanData(dictionary *scan_file_dict, const char *sectionname, SCANDATA *ResultScanData) {
       /*
       This function creates a SCANDATA structure. It reads the inputfile,
       looks the sectionname provided to the function, and assigns values
       to the SCANDATA structure from the keys.
       Section name is "scan_x", where x is an integer
       Return vlaue is a SCANDATA struct.
       */
       char section_key[200];

       //band
       strcpy(section_key,sectionname);
       strcat(section_key,":band");
       ResultScanData->band = iniparser_getint (scan_file_dict, section_key, 0);
       //scanset
       strcpy(section_key,sectionname);
       strcat(section_key,":scanset");
       ResultScanData->scanset = iniparser_getint (scan_file_dict, section_key, 0);
       //type
       strcpy(section_key,sectionname);
       strcat(section_key,":type");
       strcpy(ResultScanData->type,iniparser_getstring (scan_file_dict, section_key, "null"));
       //pol
       strcpy(section_key,sectionname);
       strcat(section_key,":pol");
       ResultScanData->pol = iniparser_getint (scan_file_dict, section_key, 0);
       //tilt
       strcpy(section_key,sectionname);
       strcat(section_key,":tilt");
       ResultScanData->tilt = iniparser_getint (scan_file_dict, section_key, 0);
       //f
       strcpy(section_key,sectionname);
       strcat(section_key,":f");
       ResultScanData->f = iniparser_getint (scan_file_dict, section_key, 0);
       //sb
       strcpy(section_key,sectionname);
       strcat(section_key,":sb");
       ResultScanData->sb = iniparser_getint (scan_file_dict, section_key, 0);
       //probe zdistance
       strcpy(section_key,sectionname);
       strcat(section_key,":zdistance");
       ResultScanData->zdistance = iniparser_getint (scan_file_dict, section_key, 260);
       if (ResultScanData->zdistance == 0)
           ResultScanData->zdistance = 260;
       //ifatten
       strcpy(section_key,sectionname);
       strcat(section_key,":ifatten");
       ResultScanData->ifatten = iniparser_getint (scan_file_dict, section_key, 0);
       //nf
       strcpy(section_key,sectionname);
       strcat(section_key,":nf");
       strcpy(ResultScanData->nf,iniparser_getstring (scan_file_dict, section_key, "null"));
       //nf_startrow
       strcpy(section_key,sectionname);
       strcat(section_key,":nf_startrow");
       ResultScanData->nf_startrow = iniparser_getint (scan_file_dict, section_key, 0);
       //ff
       strcpy(section_key,sectionname);
       strcat(section_key,":ff");
       strcpy(ResultScanData->ff,iniparser_getstring (scan_file_dict, section_key, "null"));
       //ff_startrow
       strcpy(section_key,sectionname);
       strcat(section_key,":ff_startrow");
       ResultScanData->ff_startrow = iniparser_getint (scan_file_dict, section_key, 0);
       
       //nf2
       strcpy(section_key,sectionname);
       strcat(section_key,":nf2");
       strcpy(ResultScanData->nf2,iniparser_getstring (scan_file_dict, section_key, "null"));
       //nf2_startrow
       strcpy(section_key,sectionname);
       strcat(section_key,":nf2_startrow");
       ResultScanData->nf2_startrow = iniparser_getint (scan_file_dict, section_key, 0);
       //ff2
       strcpy(section_key,sectionname);
       strcat(section_key,":ff2");
       strcpy(ResultScanData->ff2,iniparser_getstring (scan_file_dict, section_key, "null"));
       //ff2_startrow
       strcpy(section_key,sectionname);
       strcat(section_key,":ff2_startrow");
       ResultScanData->ff2_startrow = iniparser_getint (scan_file_dict, section_key, 0);
       

	  strcpy(section_key,sectionname);
	  strcat(section_key,":scanset_id");
	  strcpy(ResultScanData->scanset_id,iniparser_getstring (scan_file_dict, section_key, "null"));

	  strcpy(section_key,sectionname);
	  strcat(section_key,":scan_id");
	  strcpy(ResultScanData->scan_id,iniparser_getstring (scan_file_dict, section_key, "null"));

	  strcpy(section_key,sectionname);
	  strcat(section_key,":ts");
	  strcpy(ResultScanData->ts,iniparser_getstring (scan_file_dict, section_key, "null"));

	  strcpy(section_key,sectionname);
	  strcat(section_key,":fecfg");
	  strcpy(ResultScanData->fecfg,iniparser_getstring (scan_file_dict, section_key, "null"));
       
       //notes
       strcpy(section_key,sectionname);
       strcat(section_key,":notes");
       strcpy(ResultScanData->notes,iniparser_getstring (scan_file_dict, section_key, "-1"));

       //is it a 4545 scan?
       strcpy(section_key,sectionname);
       strcat(section_key,":4545_scan");
       strcpy(ResultScanData->is4545_scan,iniparser_getstring (scan_file_dict, section_key, "FALSE"));
       
       /*
       printf("test1\n");
       getchar();
       printf("%s\n",iniparser_getstring (scan_file_dict, "scan_1:datetime", "-1"));
       getchar();
    */
    
    
       //datetime
       strcpy(section_key,sectionname);
       strcat(section_key,":datetime");
       strcpy(ResultScanData->datetime,iniparser_getstring (scan_file_dict, section_key, "-1"));
   

       strcpy(section_key,sectionname);
       strcat(section_key,":ff_xcenter");
       ResultScanData->ff_xcenter = iniparser_getdouble (scan_file_dict, section_key, 0);
       
       strcpy(section_key,sectionname);
       strcat(section_key,":ff_ycenter");
       ResultScanData->ff_ycenter = iniparser_getdouble (scan_file_dict, section_key, 0);

       strcpy(ResultScanData->sectionname,sectionname);
       return;
}

int GetNumberOfScans(dictionary *scan_file_dict) {
    // Number of scans is assumed to be equal to the number of sections -1 for the [settings] section:
    int result = iniparser_getnsec(scan_file_dict)-1;
    return result;   
}

int RemoveKeys(dictionary *scan_file_dict) {
    int i;
    const char *sectionname;
    char section_key[50];
    for (i=0;i<GetNumberOfScans(scan_file_dict)+1;i++) {

        sectionname = iniparser_getsecname(scan_file_dict,i);

        //do for all sections except "settings"
        if(strcmp(sectionname,"settings")) {
            sprintf(section_key,"%s:type",sectionname);
    
            //Remove copol keys from crosspol sections
            if(!strcmp(iniparser_getstring(scan_file_dict, section_key,"null"),"xpol")){

                sprintf(section_key,"%s:eta_spillover",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:eta_taper",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:eta_illumination",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:ff_xcenter",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:ff_ycenter",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:nf_xcenter",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:nf_ycenter",sectionname);
                iniparser_unset(scan_file_dict, section_key);

                sprintf(section_key,"%s:az_nominal",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:el_nominal",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:delta_x",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:delta_y",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:delta_z",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:eta_phase",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:ampfit_amp",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:ampfit_width_deg",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:ampfit_u_off_deg",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:ampfit_v_off_deg",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:ampfit_d_0_90",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:ampfit_d_45_135",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:plot_copol_nfamp",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:plot_copol_nfphase",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:plot_copol_ffamp",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:plot_copol_ffphase",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:squint",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:squint_arcseconds",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:eta_total_nofocus",sectionname);
                iniparser_unset(scan_file_dict, section_key);
            }                                                   

            //Remove crosspol keys from copol sections
            if(!strcmp(iniparser_getstring(scan_file_dict, section_key,"null"),"copol")) {

                //Remove squint key for pol 0, copol scan
                sprintf(section_key,"%s:pol",sectionname);

                if(iniparser_getint(scan_file_dict, section_key,0) == 0) {
                    sprintf(section_key,"%s:squint",sectionname);
                    iniparser_unset(scan_file_dict, section_key);  
                    sprintf(section_key,"%s:squint_arcseconds",sectionname);
                    iniparser_unset(scan_file_dict, section_key);  
                }
                sprintf(section_key,"%s:max_dbdifference",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:eta_spill_co_cross",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:eta_pol_on_secondary",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:eta_pol_spill",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:eta_total_nofocus",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:plot_xpol_nfamp",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:plot_xpol_nfphase",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:plot_xpol_ffamp",sectionname);
                iniparser_unset(scan_file_dict, section_key);
                sprintf(section_key,"%s:plot_xpol_ffphase",sectionname);
                iniparser_unset(scan_file_dict, section_key);
            }
        }
    }
    return 1;
}


int GetNumberOfScanSetsForBand(dictionary *scan_file_dict, int band){
    int result=0;
    int tempmin=0;
    int num_sets=0;
    int i;
    char ibuf[5];
    char sectionname[10];
    char section_key[30];
    
    for (i=0;i<(GetNumberOfScans(scan_file_dict)+1);i++){
        strcpy(sectionname,"scan_");
        strcat(sectionname,itoa(i, ibuf, 10));  
        strcpy(section_key,sectionname);
        strcat(section_key,":band");

        if (iniparser_getint (scan_file_dict, section_key, 0)==band){
        //only look at scanset key if band is right value                     
            strcpy(sectionname,"scan_");
            strcat(sectionname,itoa(i, ibuf, 10));  
            strcpy(section_key,sectionname);
            strcat(section_key,":scanset");    
            num_sets=iniparser_getint (scan_file_dict, section_key, 0);
            if (num_sets > tempmin){
               tempmin=num_sets;
               result++;
            }
        }
    }
    return result;   
}


int GetNumberOfScanSets(dictionary *scan_file_dict){
    int result=0;
    int num_sets=0;
    int scancount=0;
    int i;
    char ibuf[5];
    char sectionname[10];
    char section_key[30];

    for (i=0;i<(GetNumberOfScans(scan_file_dict)+1);i++){    
            strcpy(sectionname,"scan_");
            strcat(sectionname,itoa(i, ibuf, 10));  
            strcpy(section_key,sectionname);
            if (DEBUGGING) {
              printf("Parsing section_key = %s\n",section_key);
            }
            strcat(section_key,":scanset");    
            num_sets=iniparser_getint (scan_file_dict, section_key, 0);
            if (DEBUGGING) {
              printf("Parsing scanset = %s\n",section_key);
            }
            if (num_sets > result){
               result=num_sets;
               scancount++;
            }
    }
    return scancount;   
}

int GetScanSetNumberArray(dictionary *scan_file_dict, int scansetarray[], int scansetarray_size) {
    int scancount=0;        // Accumulate size of array of valid [scan_n] sections seen.
    int scansetcount=0;     // To return total number of scansets found in [scan_n] sections.
    int scanset;            // Scanset value found in section.
    int i;
    char ibuf[10];
    char sectionname[10];   // "scan_n"
    char section_key[30];   // "scan_n:scanset"

    for (i = 0; i < iniparser_getnsec(scan_file_dict); i++) {
        // Build section name:
        strcpy(sectionname, "scan_");
        strcat(sectionname, itoa(i, ibuf, 10));
        
        // Build key name within section:
        strcpy(section_key, sectionname);
        strcat(section_key, ":scanset");

        // Look up the key:
        scanset = iniparser_getint(scan_file_dict, section_key, -1);
        if (scanset != -1) {
            // Found it.  Append to the array of found scansets
            scansetarray[scancount] = scanset;
            scancount++;
        }
    }

//    printf("GetUniqueArrayInt...\n");

    // Reduce the array to unique values for scanset.
    scansetcount = GetUniqueArrayInt(scansetarray, scancount);
    return scansetcount;   
}



int tokenizeDelimiter(char *input, char *tokenArray[MAX_TOKENS], char *delimiter) {
#define DEBUG 0
  int i;
  int non_blanks = 0;
  int tokens = 0;
  int unprintable;
  int length;

  if (input == NULL) {
    return(0);
  }
  length = strlen(input);
  if (DEBUG) {
    fprintf(stderr,"input=%s, length=%d\n",input,length);
  }
  if (DEBUG) {
    fprintf(stderr,"done.\n");
  }
  if (length > 0) {
    if (length < 5) {
      if (DEBUG) {
        fprintf(stderr,"about to call 'for' loop\n");
      }
      unprintable = 0;
      for (i=0; i<length; i++) {
        if ((input[i] >= ' ' && input[i] <= 'z') || input[i]==0xa) {
        } else {
          if (DEBUG) {
            fprintf(stderr,"Found unprintable character at %d\n",i);
          }
          unprintable = 1;
        }
      }
      if (unprintable == 1) {
        memset(input, '\0', strlen(input));
        length = 0;
      }
    }
    for (i=0; i<length; i++) {
      if (input[i] > ' ' && input[i] <= 'z') {
        non_blanks = 1;
        break;
      }
    }
    if (non_blanks == 0) {
      if (DEBUG) {
        fprintf(stderr,"non_blanks = 0, length=%d\n",length);
      }
      return(0);
    } 
    if (DEBUG) {
      fprintf(stderr,"about to run strtok(%s,%s)\n",input,delimiter);
    }
    tokenArray[tokens++] = strtok(input,delimiter);
    if (DEBUG) {
      fprintf(stderr,"tokenArray[0] = %s\n",tokenArray[0]);
    }
    while ((tokenArray[tokens] = strtok(NULL,delimiter)) != NULL) {
      tokens++;
    }
  }
  return(tokens);
}

int PickNominalAngles(int almaBand, float *xtarget, float *ytarget, PointingOptionType pointing) {
    switch (pointing) {
        case ACTUAL:
            *xtarget = 0;
            *ytarget = 0;
            break;

        case NOMINAL:
            switch (almaBand) {
                case 1:  *xtarget = -1.7553; *ytarget = -1.7553; break;
                case 2:  *xtarget = -1.7553; *ytarget = +1.7553; break;
                case 3:  *xtarget = +0.3109; *ytarget = +1.7345; break;
                case 4:  *xtarget = +0.3109; *ytarget = -1.7345; break;
                case 5:  *xtarget = +1.6867; *ytarget = +1.6867; break;
                case 6:  *xtarget = +1.6867; *ytarget = -1.6867; break;
                case 7:  *xtarget = +0.9740; *ytarget = +0.0000; break;
                case 8:  *xtarget = +0.0000; *ytarget = +0.9740; break;
                case 9:  *xtarget = +0.0000; *ytarget = -0.9740; break;
                case 10: *xtarget = -0.9740; *ytarget = +0.0000; break;
                default:
                    printf("Illegal band number = %d\n", almaBand);
                    return(-1);
            }
            break;

        case ACA_7METER:
            switch (almaBand) {
                case 1:  *xtarget = -2.943499; *ytarget = -2.943499; break;
                case 2:  *xtarget = -2.898850; *ytarget = +2.898850; break;
                case 3:  *xtarget = +0.521949; *ytarget = +2.918507; break;
                case 4:  *xtarget = +0.549947; *ytarget = -3.120381; break;
                case 5:  *xtarget = +1.874227; *ytarget = +1.874227; break;
                case 6:  *xtarget = +1.990572; *ytarget = -1.990572; break;
                case 7:  *xtarget = +0.764417; *ytarget = +0.000000; break;
                case 8:  *xtarget = +0.000000; *ytarget = +0.757841; break;
                case 9:  *xtarget = +0.000000; *ytarget = -0.735064; break;
                case 10: *xtarget = -0.735064; *ytarget = +0.000000; break;
                default:
                    printf("Illegal band number = %d\n", almaBand);
                    return(-1);
            }
            break;

        case BAND1_TEST:
            switch (almaBand) {
                case 1:  *xtarget = 0; *ytarget = +2.48; break;

                default:
                    printf("Illegal band number for pointingOption BAND1_TEST = %d\n", almaBand);
                    return(-1);
            }
            break;

        default:
            printf("Illegal value for pointingOption = %d\n", (int) pointing);
            return(-1);
    }
    return(0);
}


int GetUniqueArrayInt(int invals[], int arrsize) {
    // Input is an array of integers
    // Duplicate values are deleted from the array
    // Resulting array consists of all uniqe values, followed
    // by "-1" in the remaining array slots.
    // Return value is number of unique elements in array.
    char tempstr1[900];
    char tempstr2[900];
    char tempval[20];
    int i, ucount = 0;
    char delim[] = "_";
    char *result = NULL;

    // Initialize the two string buffers to empty:
    tempstr1[0] = tempstr2[0] = '\0';

    // Print all the values into the two temporary strings, separated by the delimiter:
    for (i = 0; i < arrsize; i++) {
        sprintf(tempval,"%d%s", invals[i], delim);

        // Copy all values into tempstr1
        strcat(tempstr1, tempval);

        // Copy only unique values into tempstr2
        if (strstr(tempstr2, tempval) == NULL) {
            strcat(tempstr2, tempval);
            // accumulate a count of unique values seen:
            ucount++;
        }
    }

//    printf("ucount=%d\n"
//           "tempstr1=%s\n"
//           "tempstr2=%s\n", ucount, tempstr1, tempstr2);

    //Replace input array with array containing only unique values, filled with -1 for the rest:
    result = strtok(tempstr2, delim);
    for (i = 0; i < arrsize; i++) {
        if (i < ucount) {
            invals[i] = atoi(result);
            result = strtok(NULL, delim);
        } else {
            invals[i] = -1;
        }
    }
    return ucount;
}

int ReplaceDelimiter(char input[400], const char *olddelim, const char *newdelim) {
    char resultstr[400];
    char *temptoken;

    printf("input= %s\n",input);
    
    resultstr[0] = '\0';
    temptoken = strtok(input,olddelim);
    sprintf(resultstr,"%s%s%s",resultstr,temptoken,newdelim);
    while(temptoken) {
        temptoken = strtok(NULL,olddelim);
        if (temptoken){             
           sprintf(resultstr,"%s%s%s",resultstr,temptoken,newdelim);
        }
    }
    strncpy(input,resultstr,strlen(resultstr)-1);
    
    
    
    return 1;
}
