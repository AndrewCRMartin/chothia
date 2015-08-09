/*************************************************************************

   Program:    Chothia
   File:       chothia.c
   
   Version:    V1.2
   Date:       30.11.95
   Function:   Assign canonical classes and display reasons for 
               mismatches.
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 387 7050 X 3284
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  17.05.95 Original
   V1.1  18.08.95 No longer fails if residues missing. Just reports it.
                  Fixed bug in FindRes(). When asking for xxB was finding
                  xx rather than xxA
   V1.2  30.11.95 Fixed problem in reading blank lines in Chothia data 
                  file if they had spaces rather than being truly blank

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>


#include "bioplib/macros.h"
#include "bioplib/general.h"


/************************************************************************/
/* Defines and macros
*/
#define ENV_KABATDIR "KABATDIR"  /* Environment variable for Kabat      */
                                 /* directory                           */
#define MAXCHOTHRES  80          /* Max number of key residues per class*/
#define MAXBUFF      160         /* General buffer size                 */
#define MAXSEQ       300         /* Max length of light + heavy chains  */
#define NCDR         5           /* Number of CDRs to process           */

#define TERMALPHA(x) do {  int _termalpha_j;                  \
                        for(_termalpha_j=0;                   \
                            (x)[_termalpha_j];                \
                            _termalpha_j++)                   \
                        {  if(isalpha((x)[_termalpha_j]))     \
                           {  (x)[_termalpha_j] = '\0';       \
                              break;                          \
                     }  }  }  while(0)
      
typedef struct _chothia
{
   struct _chothia *next;
   int             length;
   char            LoopID[8],
                   class[8],
                   source[MAXBUFF],
                   resnum[MAXCHOTHRES][8],
                   restype[MAXCHOTHRES][24];
}  CHOTHIA;

typedef struct
{
   char            resnum[8],
                   seq;
}  SEQUENCE;

typedef struct
{
   char name[8],
        start[8],
        stop[8];
}  LOOP;


/************************************************************************/
/* Globals
*/
CHOTHIA   *gChothia = NULL;        /* Linked list of Chothia data       */
BOOL      gCanonChothNum = FALSE;  /* Data file uses Chothia numbering  */

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ReadChothiaData(char *filename);
int ReadInputData(FILE *in, SEQUENCE *Sequence);
void ReportCanonicals(FILE *out, SEQUENCE *Sequence, int NRes, 
                      BOOL verbose);
int FindRes(SEQUENCE *Sequence, int NRes, char *res);
void ReportACanonical(FILE *out, char *LoopName, int LoopLen, 
                      SEQUENCE *Sequence, int NRes, BOOL verbose);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *ChothiaFile, BOOL *verbose);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for assigning canonicals

   16.05.95 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   char     InFile[MAXBUFF],
            OutFile[MAXBUFF],
            ChothiaFile[MAXBUFF];
   FILE     *in  = stdin,
            *out = stdout;
   SEQUENCE Sequence[MAXSEQ];
   int      NRes;
   BOOL     verbose;

   strcpy(ChothiaFile,"chothia.dat");

   if(ParseCmdLine(argc, argv, InFile, OutFile, ChothiaFile, &verbose))
   {
      if(OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if(ReadChothiaData(ChothiaFile))
         {
            if((NRes = ReadInputData(in, Sequence)) != 0)
            {
               ReportCanonicals(out, Sequence, NRes, verbose);
            }
            else
            {
               fprintf(stderr,"Chothia: Error in input data\n");
               Usage();
               return(1);
            }
         }
         else
         {
            fprintf(stderr,"Chothia: Unable to read Chothia datafile\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"Chothia: Unable to open i/o files\n");
         return(1);
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>BOOL ReadChothiaData(char *filename)
   ------------------------------------
   Input:   char *filename      The Chothia data filename
   Returns: BOOL                Success?
   Globals: CHOTHIA *gChothia   Linked list of Chothia data
            BOOL gCanonChothNum Chothia (rather than Kabat) numbering 
                                used in file

   Reads a Chothia canonical definition file. This file has the format:
   LOOP loopid class length
  [SOURCE ............................ ]
   resid types
   resid types
   ...

   16.05.95 Original based on ReadChothiaData() from KabatMan
   30.11.95 Remove leading spaces from strings read from file
   07.05.96 Handles the CHOTHIANUMBERING keyword
*/
BOOL ReadChothiaData(char *filename)
{
   FILE    *fp;
   char    buffer[MAXBUFF],
           word[40],
           *chp,
           *buffp;
   CHOTHIA *p = NULL;
   int     count = 0;
   BOOL    NoEnv;
   
   /* Open the data file                                                */
   if((fp=OpenFile(filename,ENV_KABATDIR,"r",&NoEnv))==NULL)
   {
      return(FALSE);
   }

   gCanonChothNum = FALSE;

   while(fgets(buffer,160,fp))
   {
      TERMINATE(buffer);
      buffp = buffer;
      while(isspace(*buffp))
         buffp++;
      
      if(strlen(buffp) && buffp[0] != '!' && buffp[0] != '#')
      {
         /* Handle the SOURCE keyword                                   */
         if(!upstrncmp(buffp,"SOURCE",6))
         {
            if(p!=NULL)
            {
               /* Strip out the SOURCE keyword                          */
               chp = GetWord(buffp,word);
               /* Store the text                                        */
               strcpy(p->source, chp);
            }
         }
         else if(!upstrncmp(buffp,"CHOTHIANUM",10))
         {
            gCanonChothNum = TRUE;
         }
         else if(!upstrncmp(buffp,"LOOP",4))   /* Start of entry       */
         {
            /* Terminate the previous list of resnums                   */
            if(p!=NULL)
               strcpy(p->resnum[count],"-1");
            
            /* Allocate space in linked list                            */
            if(gChothia == NULL)
            {
               INIT(gChothia,CHOTHIA);
               p = gChothia;
            }
            else
            {
               ALLOCNEXT(p,CHOTHIA);
            }
            if(p==NULL) return(FALSE);
            
            /* Strip out the word LOOP                                  */
            chp = GetWord(buffp,word);
            /* Get the loop id                                          */
            chp = GetWord(chp,p->LoopID);
            /* Get the class name                                       */
            chp = GetWord(chp,p->class);
            /* Get the loop length                                      */
            chp = GetWord(chp,word);
            sscanf(word,"%d",&(p->length));
            
            /* Set the resnum counter to zero                           */
            count = 0;
         }
         else
         {
            /* Not the start of an entry, so must be a resid/type pair  */
            if(p!=NULL)
            {
               chp = GetWord(buffp,p->resnum[count]);
               chp = GetWord(chp,p->restype[count]);
               if(++count > MAXCHOTHRES)
               {
                  fprintf(stderr,"Error: (Reading Chothia file) Too many \
Chothia key residues.\n");
                  return(FALSE);
               }
            }
         }
      }
   }
   
   /* Terminate the previous list of resnums                            */
   if(p!=NULL)
      strcpy(p->resnum[count],"-1");
   
   return(TRUE);
}


/************************************************************************/
/*>int ReadInputData(FILE *in, SEQUENCE *Sequence)
   -----------------------------------------------
   Input:   FILE     *in          Input data file pointer
   Output:  SEQUENCE *Sequence    Sequence array
   Returns: int                   Length of sequence

   16.05.95 Original    By: ACRM
*/
int ReadInputData(FILE *in, SEQUENCE *Sequence)
{
   char buffer[MAXBUFF],
        word[16],
        *chp;
   int  count = 0;
   
   while(fgets(buffer, MAXBUFF, in))
   {
      if((buffer[0] == 'L' || buffer[0] == 'H') &&
         isdigit(buffer[1]))
      {
         chp = GetWord(buffer, word);
         strcpy(Sequence[count].resnum, word);
         
         chp = GetWord(chp, word);
         if(strlen(word) == 0)
            return(0);
         
         Sequence[count].seq = word[0];

         count++;
      }
   }
   
   return(count);
}

      
/************************************************************************/
/*>void ReportCanonicals(FILE *out, SEQUENCE *Sequence, int NRes, 
                         BOOL verbose)
   ---------------------------------------------------------------
   Input:   FILE     *out          Output file pointer
            SEQUENCE *Sequence     Sequence array
            int      NRes          Length of sequence
            BOOL     verbose       Flag to display reasons

   Reports the canonical classes for all 6 loops. Calls ReportACanonical()
   to do the work.

   16.05.95 Original    By: ACRM
   18.08.95 No longer fails if unable to find a residue; just reports the
            fact. Now returns void.
*/
void ReportCanonicals(FILE *out, SEQUENCE *Sequence, int NRes, 
                      BOOL verbose)
{
   int         loop,
               len,
               start,
               stop;
   static LOOP LoopDef[] = 
   {  {  "L1", "L24", "L34"  },
      {  "L2", "L50", "L56"  },
      {  "L3", "L89", "L97"  },
      {  "H1", "H26", "H35B" },
      {  "H2", "H50", "H58"  },
      {  "H3", "H95", "H102" }
   }  ;
   
   for(loop=0; loop<NCDR; loop++)
   {
      if((start = FindRes(Sequence, NRes, LoopDef[loop].start)) == (-1))
      {
         fprintf(stderr,"Chothia: Unable to find residue %s in input\n",
                 LoopDef[loop].start);
         fprintf(out,"CDR %s  Missing Residues\n",LoopDef[loop].name);
         continue;
      }

      if((stop = FindRes(Sequence, NRes, LoopDef[loop].stop)) == (-1))
      {
         fprintf(stderr,"Chothia: Unable to find residue %s in input\n",
                 LoopDef[loop].stop);
         fprintf(out,"CDR %s  Missing Residues\n",LoopDef[loop].name);
         continue;
      }
         
      len   = 1 + stop - start;
      
      ReportACanonical(out, LoopDef[loop].name, len, Sequence, NRes,
                       verbose);
   }
}


/************************************************************************/
/*>int FindRes(SEQUENCE *Sequence, int NRes, char *res)
   ----------------------------------------------------
   Input:   SEQUENCE   *Sequence      The sequence array
            int        NRes           Length of sequence array
            char       *res           The res number to find
   Returns: int                       Offset into Sequence array
                                      -1 if not found

   Finds a residue label in the sequence array. If no exact match is found
   it runs through again checking without any insert codes. This assumes
   that residues which we search for can be substituted by an earlier 
   residue in the sequence.

   16.05.95 Original    By: ACRM
   17.05.95 Added no-insert checking
   18.08.95 Fixed failure when 35B specified and had 35A in sequence
            (Used to find 35 instead of 35A).
*/
int FindRes(SEQUENCE *Sequence, int NRes, char *InRes)
{
   int  i,
        InsPos = (-1);
   char res[16];

   strcpy(res,InRes);
   
   /* Check full residue names                                          */
   for(i=0; i<NRes; i++)
   {
      if(!strcmp(Sequence[i].resnum, res))
         return(i);
   }

   /* Exact match failed, see if there is an insert code. Counts from 1
      not zero as 0 is the chain name
   */
   for(i=1; i<strlen(res); i++)
   {
      if(isalpha(res[i]))
      {
         InsPos = i;
         break;
      }
   }

   /* If there wasn't an insert code, then we've definitely failed      */
   if(InsPos < 0)
      return(-1);

   /* Step the insert code down                                         */
   res[InsPos]--;
   
   while(res[InsPos] >= 'A')
   {
      for(i=0; i<NRes; i++)
      {
         if(!strcmp(Sequence[i].resnum, res))
            return(i);
      }
   
      /* Step the insert code down                                      */
      res[InsPos]--;
   }

   /* Finally try without the insert code                               */
   res[InsPos] = '\0';
   for(i=0; i<NRes; i++)
   {
      char buff[16];

      strcpy(buff, Sequence[i].resnum);
      TERMALPHA(buff+1);

      if(!strcmp(buff, res))
         return(i);
   }
   
   /* Tried everything!                                                 */
   return(-1);
}


/************************************************************************/
/*>void ReportACanonical(FILE *out, char *LoopName, int LoopLen, 
                         SEQUENCE *Sequence, int NRes, BOOL verbose)
   --------------------------------------------------------------
   Input:   FILE     *out          Output file pointer
            char     *LoopName     Name of a loop (e.g. L1)
            int      LoopLen       Length of the loop
            SEQUENCE *Sequence     Sequence array
            int      NRes          Length of sequence
            BOOL     verbose       Flag to display reasons
   Returns: BOOL                   Success?

   Reports the canonical class for an individual loop

   16.05.95 Original    By: ACRM
   17.05.95 Only prints source data if verbose
*/
void ReportACanonical(FILE *out, char *LoopName, int LoopLen, 
                      SEQUENCE *Sequence, int NRes, BOOL verbose)
{
   CHOTHIA *p,
           *best = NULL;
   int     i,
           res,
           NMismatch,
           MinMismatch = 10000;
   BOOL    OK = FALSE;   /* Assume this is not a canonical              */
   
   /* Run through the linked list of Canonical definitions              */
   for(p=gChothia; p!=NULL; NEXT(p))
   {
      /* If the Loop name and length match                              */
      if(!strcmp(p->LoopID, LoopName) && (p->length == LoopLen))
      {
         OK = TRUE;   /* Now assume it is a canonical                   */
         NMismatch = 0;
         
         /* Check each residue specified by this canonical definition   */
         for(i=0; strcmp(p->resnum[i], "-1"); i++)
         {
            res = FindRes(Sequence, NRes, p->resnum[i]);

            /* As soon as we find a disallowed residue type, set the OK
               flag to false and increment the mismatch counter
            */
            if(!strchr(p->restype[i], Sequence[res].seq))
            {
               OK = FALSE;
               NMismatch++;
            }
         }

         if(OK)   /* If everything OK, break out of the linked list     */
         {
            break;
         }
         else     /* See if this has the fewest mismatches              */
         {
            if(NMismatch < MinMismatch)
            {
               MinMismatch = NMismatch;
               best = p;
            }
         }
         
      }
   }

   if(OK)
   {
      fprintf(out,"CDR %s  Class %-3s", LoopName, p->class);
      if(verbose && strlen(p->source))
         fprintf(out," %s", p->source);
      fprintf(out,"\n");
   }
   else
   {
      fprintf(out,"CDR %s  Class ?  \n", LoopName);
   
      if(verbose)
      {
         if(best==NULL)
         {
            fprintf(out, "! Length mismatch\n");
         }
         else
         {
            fprintf(out, "! Similar to class %s, but:\n", best->class);

            /* Display each mismatch for this canonical definition      */
            for(i=0; strcmp(best->resnum[i], "-1"); i++)
            {
               res = FindRes(Sequence, NRes, best->resnum[i]);
               if(!strchr(best->restype[i], Sequence[res].seq))
               {
                  fprintf(out, "!    %s = %c (allows: %s)\n", 
                          best->resnum[i],
                          Sequence[res].seq,
                          best->restype[i]);
               }
            }
         }
      }
   }   
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   16.05.95 Original    By: ACRM
   18.08.95 V1.1   
   30.11.95 V1.2
*/
void Usage(void)
{
   fprintf(stderr,"\nChothia V1.2 (c) 1995, Dr. Andrew C.R. Martin, \
UCL\n\n");

   fprintf(stderr,"Usage: chothia [-c filename] [-v] [input.seq \
[output.dat]]\n");
   fprintf(stderr,"               -c Specify Chothia datafile (Default: \
chothia.dat)\n");
   fprintf(stderr,"               -v Verbose; give explanations when \
no canonical found\n");
   fprintf(stderr,"       I/O is through stdin/stdout if files are not \
specified.\n\n");

   fprintf(stderr,"Chothia is a program to assign canonical classes to \
an antibody sequence.\n");
   fprintf(stderr,"Input to the program is a listing of Kabat residue \
numbers and the\n");
   fprintf(stderr,"one-letter code amino acid name for the residue at \
each position. Such\n");
   fprintf(stderr,"a file may be generated from a PIR file using the \
program KabatSeq.\n\n");

   fprintf(stderr,"The program will look for the datafile first in the \
current directory\n");
   fprintf(stderr,"and then in the directory specified by the %s \
environment variable.\n", ENV_KABATDIR);
   fprintf(stderr,"This data file is also used by the KabatMan database \
software.\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *ChothiaFile, BOOL *verbose)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            char   *ChothiaFile Chothia data file
            BOOL   *verbose     Flag to show details of mismatches
   Returns: BOOL                Success?

   Parse the command line
   
   16.05.95 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *ChothiaFile, BOOL *verbose)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *verbose = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
            argc--;
            argv++;
            strcpy(ChothiaFile, argv[0]);
            break;
         case 'v':
            *verbose = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


