/*************************************************************************

   Program:    Chothia
   File:       chothia.c
   
   Version:    V2.3
   Date:       12.10.21
   Function:   Assign canonical classes and display reasons for 
               mismatches.
   
   Copyright:  (c) Prof. Andrew C. R. Martin, UCL 1995-2021
   Author:     Prof. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
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

   Must be linked with KabCho.c from KabatMan


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
   V1.3  08.05.96 Now handles data files and sequence data with Chothia 
                  numbering as well as Kabat numbering.
   V1.4  09.05.96 Mismatches weren't being reported correctly if the
                  datafile and sequence file didn't use the same numbering
   V1.5  30.05.96 Reports numbering scheme when printing mismatches and
                  correctly handles deleted residues
   V1.6  19.12.08 Modified to use new GetWord() and various bug-avoiding
                  changes.
   V1.7  13.02.14 All errors and warnings now printed in consistent 
                  format.
                  Now allows 1-letter or 3-letter code sequence files
   V2.0  14.02.11 Modified to allow new PRIORITY and SUBORDINATE keywords
                  in data file to deal with breakdown in classical 
                  canonicals
   V2.1  09.08.15 Added -L and -H switches to do single chains
   V2.2  14.12.16 File reading terminates at CR as well as LF to deal 
                  with windows files.
                  Changed to new blXXX() Bioplib functions
   V2.3  12.10.21 MAXBUFF bumped to 240 and MAXSEQ to 3000 (inherited 
                  from abYsis version)

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "bioplib/seq.h"

/************************************************************************/
/* Defines and macros
*/
#define ENV_KABATDIR "KABATDIR"  /* Environment variable for Kabat      */
                                 /* directory                           */
#define MAXCHOTHRES  80          /* Max number of key residues per class*/
#define MAXBUFF      240         /* General buffer size                 */
#define MAXSEQ       3000        /* Max length of light + heavy chains  */
#define MAXEXPSEQ    300         /* Expected max light + heavy          */
#define NCDR         5           /* Number of CDRs to process           */
#define MAXWORD      40          /* Max length of an extracted word     */
#define SMALLWORD    16          /* Length of small extracted word      */

/* Terminates a string at the first alphabetic character                */
#define TERMALPHA(x) do {  int _termalpha_j;                  \
                        for(_termalpha_j=0;                   \
                            (x)[_termalpha_j];                \
                            _termalpha_j++)                   \
                        {  if(isalpha((x)[_termalpha_j]))     \
                           {  (x)[_termalpha_j] = '\0';       \
                              break;                          \
                     }  }  }  while(0)
      
/* Linked list to store information on canonical class definitions      */
typedef struct _chothia
{
   struct _chothia *next,                           /* Linked list      */
                   *priority_over,                  /* Priority over which
                                                       other classes when
                                                       key residues 
                                                       clash            */
                   *subordinate_to;                 /* Suborinate to which
                                                       other classes when
                                                       key residues
                                                       clash            */
   int             length;                          /* Loop length      */
   char            LoopID[SMALLWORD],               /* CDR (L1, L2, etc */
                   class[SMALLWORD],                /* Class name       */
                   source[MAXBUFF],                 /* Info on class
                                                       maybe including PDB
                                                       code in []       */
                   resnum[MAXCHOTHRES][SMALLWORD],  /* Key positions    */
                   restype[MAXCHOTHRES][MAXWORD],   /* Allowed residues */
                   subordinate[SMALLWORD],          /* Class name to
                                                       which this class is
                                                       subordinate
                                                       (0 or 1)         */
                   priority[SMALLWORD];             /* Class name over
                                                       which this class
                                                       takes priority
                                                       (0 or 1)         */
   int             npriority,                       /* Number over which
                                                       this class takes
                                                       priority (0 or 1)*/
                   nsubordinate;                    /* Number of classes
                                                       to which this is
                                                       subordinate 
                                                       (0 or 1)         */
   
}  CHOTHIA;

/* Input sequence data (array) - residue number label and amino acid    */
typedef struct
{
   char            resnum[SMALLWORD],
                   seq;
}  SEQUENCE;

/* Definitions of CDR loop boundaries (array)                           */
typedef struct
{
   char name[SMALLWORD],
        start[SMALLWORD],
        stop[SMALLWORD];
}  LOOP;


/************************************************************************/
/* Globals
*/
CHOTHIA   *gChothia = NULL;         /* Linked list of Chothia data      */
BOOL      gCanonChothNum = FALSE,   /* Data file uses Chothia numbering?*/
          gChothiaNumbered = FALSE; /* Sequence data uses Chothia 
                                       numbering?                       */

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ReadChothiaData(char *filename);
int  ReadInputData(FILE *in, SEQUENCE *Sequence);
void ReportCanonicals(FILE *out, SEQUENCE *Sequence, int NRes, 
                      BOOL verbose, char chain);
int  FindRes(SEQUENCE *Sequence, int NRes, char *res);
void ReportACanonical(FILE *out, char *LoopName, int LoopLen, 
                      SEQUENCE *Sequence, int NRes, BOOL verbose,
                      char *cdr, int cdrlen);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *ChothiaFile, BOOL *verbose, char *chain);
char *KabCho(char *cdr, int length, char *kabspec);
char *ChoKab(char *cdr, int length, char *kabspec);
int TestThisCanonical(CHOTHIA *p, char *LoopName, int LoopLen,
                      SEQUENCE *Sequence, int NRes, 
                      char *cdr1, int cdr1len);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for assigning canonicals

   16.05.95 Original    By: ACRM
   19.12.08 Changed strcpy() to strncpy()
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
   char     chain = ' ';

   strncpy(ChothiaFile,"chothia.dat", MAXBUFF);

   if(ParseCmdLine(argc, argv, InFile, OutFile, ChothiaFile, &verbose,
                   &chain))
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         if(ReadChothiaData(ChothiaFile))
         {
            if((NRes = ReadInputData(in, Sequence)) != 0)
            {
               ReportCanonicals(out, Sequence, NRes, verbose, chain);
            }
            else
            {
               fprintf(stderr,"Error (chothia): Error in input data\n");
               return(1);
            }
         }
         else
         {
            fprintf(stderr,"Error (chothia): Unable to read Chothia \
datafile\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"Error (chothia): Unable to open i/o files\n");
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
   19.12.08 Uses MAXWORD and updated for new GetWord()
            Fixed explicit 160 in fgets to MAXBUFF
            Changed strcpy() to strncpy()
   14.02.11 Added PRIORITY and SUBORDINATE keywords
   14.12.16 Changed to blGetWord()
*/
BOOL ReadChothiaData(char *filename)
{
   FILE    *fp;
   char    buffer[MAXBUFF],
           word[MAXWORD],
           *chp,
           *buffp;
   CHOTHIA *p = NULL;
   int     count = 0;
   BOOL    NoEnv;
/*           GotSubPri = FALSE; */
   
   /* Open the data file                                                */
   if((fp=blOpenFile(filename,ENV_KABATDIR,"r",&NoEnv))==NULL)
   {
      return(FALSE);
   }

   gCanonChothNum = FALSE;

   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);
      buffp = buffer;
      while(isspace(*buffp))
         buffp++;
      
      if(strlen(buffp) && buffp[0] != '!' && buffp[0] != '#')
      {
         /* Handle the SOURCE keyword                                   */
         if(!blUpstrncmp(buffp,"SOURCE",6))
         {
            if(p!=NULL)
            {
               /* Strip out the SOURCE keyword                          */
               chp = blGetWord(buffp,word,MAXWORD);
               /* Store the text                                        */
               strncpy(p->source, chp, MAXWORD);
            }
         }
         else if(!blUpstrncmp(buffp,"PRIORITY",8))
         {
/*            GotSubPri = TRUE; */
            /* Strip out the PRIORITY keyword                           */
            chp = blGetWord(buffp,word,MAXWORD);
            /* And grab the class name over which this takes priority   */
            chp = blGetWord(chp,p->priority,SMALLWORD);
            p->npriority = 1;
         }
         else if(!blUpstrncmp(buffp,"SUBORDINATE",11))
         {
/*            GotSubPri = TRUE; */
            /* Strip out the SUBORDINATE keyword                        */
            chp = blGetWord(buffp,word,MAXWORD);
            /* And grab the class name to which this is subordinate     */
            chp = blGetWord(chp,p->subordinate,SMALLWORD);
            p->nsubordinate = 1;
         }
         else if(!blUpstrncmp(buffp,"CHOTHIANUM",10))
         {
            gCanonChothNum = TRUE;
         }
         else if(!blUpstrncmp(buffp,"LOOP",4))    /* Start of entry     */
         {
            /* Terminate the previous list of resnums                   */
            if(p!=NULL)
               strncpy(p->resnum[count],"-1",SMALLWORD);
            
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
            
            /* 14.02.11 Initialize the PRIORITY and SUBORDINATE fields  */
            p->priority_over = p->subordinate_to = NULL;
            p->npriority     = p->nsubordinate   = 0;

            /* Strip out the word LOOP                                  */
            chp = blGetWord(buffp,word,MAXWORD);
            /* Get the loop id                                          */
            chp = blGetWord(chp,p->LoopID,SMALLWORD);
            /* Get the class name                                       */
            chp = blGetWord(chp,p->class,SMALLWORD);
            /* Get the loop length                                      */
            chp = blGetWord(chp,word,MAXWORD);
            sscanf(word,"%d",&(p->length));
            
            /* Set the resnum counter to zero                           */
            count = 0;
         }
         else
         {
            /* Not the start of an entry, so must be a resid/type pair  */
            if(p!=NULL)
            {
               chp = blGetWord(buffp,p->resnum[count],SMALLWORD);
               chp = blGetWord(chp,p->restype[count],MAXWORD);
               if(++count > MAXCHOTHRES)
               {
                  fprintf(stderr,"Error: (chothia) Too many key \
residues when reading Chothia file.\n");
                  return(FALSE);
               }
            }
         }
      }
   }
   
   /* Terminate the previous list of resnums                            */
   if(p!=NULL)
      strncpy(p->resnum[count],"-1",SMALLWORD);

   /* 14.02.11 If we have any PRIORITY/SUBORDINATEs then set the 
      information for the pointers rather than simple text labels
   */
   for(p=gChothia; p!=NULL; NEXT(p))
   {
      /* See if this takes priority over anything else                  */
      if(p->npriority)
      {
         int     nmatch = 0;
         CHOTHIA *q = NULL;
         for(q=gChothia; q!=NULL; NEXT(q))
         {
            if(!strcmp(p->priority, q->class))
            {
               p->priority_over = q;
               nmatch++;
            }
         }
         if(nmatch > 1)
         {
            fprintf(stderr,"Chothia: Error 1, Loop %s takes priority \
over %s, but %s matches %d classes\n", 
                    p->class, p->priority, p->priority, nmatch);
            return(FALSE);
         }
         if(nmatch < 1)
         {
            fprintf(stderr,"Chothia: Error 2, Loop %s takes priority \
over %s, but %s not found as a valid canonical name\n", 
                    p->class, p->priority, p->priority);
            return(FALSE);
         }
         if(p->length != p->priority_over->length)
         {
            fprintf(stderr,"Chothia: Error 5, Loop %s takes priority \
over %s, but lengths do not match\n", 
                    p->class, p->priority);
            return(FALSE);
         }
         
         
      }
      
      /* See if this is subordinate to anything else                    */
      if(p->nsubordinate)
      {
         int     nmatch = 0;
         CHOTHIA *q = NULL;
         for(q=gChothia; q!=NULL; NEXT(q))
         {
            if(!strcmp(p->subordinate, q->class))
            {
               p->subordinate_to = q;
               nmatch++;
            }
         }
         if(nmatch > 1)
         {
            fprintf(stderr,"Chothia: Error 3, Loop %s is subordinate \
to %s, but %s matches %d classes\n", 
                    p->class, p->subordinate, p->subordinate, nmatch);
            return(FALSE);
         }
         if(nmatch < 1)
         {
            fprintf(stderr,"Chothia: Error 4, Loop %s is subordinate \
to %s, but %s not found as a valid canonical name\n", 
                    p->class, p->subordinate, p->subordinate);
            return(FALSE);
         }
         if(p->length != p->subordinate_to->length)
         {
            fprintf(stderr,"Chothia: Error 6, Loop %s is subordinate \
to %s, but lengths do not match\n", 
                    p->class, p->priority);
            return(FALSE);
         }
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>int ReadInputData(FILE *in, SEQUENCE *Sequence)
   -----------------------------------------------
   Input:   FILE     *in          Input data file pointer
   Output:  SEQUENCE *Sequence    Sequence array
   Returns: int                   Length of sequence

   16.05.95 Original    By: ACRM
   19.12.08 Changed strcpy() to strncpy()
            Changed word[16] to word[MAXWORD]
            New GetWord()
            Checks number of residues in file
            Added check on residue names of '-'
   14.12.16 Changed to blGetWord()
*/
int ReadInputData(FILE *in, SEQUENCE *Sequence)
{
   char buffer[MAXBUFF],
        word[MAXWORD],
        *chp;
   int  count = 0;
   
   while(fgets(buffer, MAXBUFF, in))
   {
      TERMINATE(buffer);  /* 13.02.14 Added this                        */
      TERMINATECR(buffer);/* 14.12.16 Added this                        */
      
      if((buffer[0] == 'L' || buffer[0] == 'H') &&
         isdigit(buffer[1]))
      {
         chp = blGetWord(buffer, word, MAXWORD);
         strncpy(Sequence[count].resnum, word, SMALLWORD);
         
         chp = blGetWord(chp, word, MAXWORD);
         if(strlen(word) == 0)
            return(0);
         
         if(word[0] != '-')
         {
            if(strlen(word) == 3)
            {
               Sequence[count].seq = blThrone(word);
            }
            else if(strlen(word) == 1)
            {
               Sequence[count].seq = word[0];
            }
            else
            {
               fprintf(stderr,"Warning (chothia): illegal residue name: \
%s\n", word);
               fprintf(stderr,"                   residue ignored.\n");
            }
            
            
            if(++count >= MAXSEQ)
            {
               fprintf(stderr,"Error (chothia): Too many residues in \
sequence file\n");
               return(0);
            }
         }
      }
   }

   if(count > MAXEXPSEQ)
   {
      fprintf(stderr,"Warning (chothia): %d residues in input file. \
Expect <%d. Maybe two antibodies?\n", count, MAXEXPSEQ);
   }
   
   return(count);
}

      
/************************************************************************/
/*>void ReportCanonicals(FILE *out, SEQUENCE *Sequence, int NRes, 
                         BOOL verbose, char chain)
   ---------------------------------------------------------------
   Input:   FILE     *out          Output file pointer
            SEQUENCE *Sequence     Sequence array
            int      NRes          Length of sequence
            BOOL     verbose       Flag to display reasons
            char     chain         Chain to handle (both if eq '')

   Reports the canonical classes for all 6 loops. Calls ReportACanonical()
   to do the work.

   16.05.95 Original    By: ACRM
   18.08.95 No longer fails if unable to find a residue; just reports the
            fact. Now returns void.
   08.05.96 Modified to determine length of CDR1 in each chain and pass
            it to the ReportACanonical() routine
   19.12.08 Changed strcpy() to strncpy()
   09.08.15 Added chain handling
*/
void ReportCanonicals(FILE *out, SEQUENCE *Sequence, int NRes, 
                      BOOL verbose, char chain)
{
   int         loop,
               len,
               cdr1len,
               start,
               stop,
               firstLoop,
               lastLoop;
   char        cdr1[SMALLWORD];
   static LOOP LoopDef[] = 
   {  {  "L1", "L24", "L34"  },
      {  "L2", "L50", "L56"  },
      {  "L3", "L89", "L97"  },
      {  "H1", "H26", "H35B" },
      {  "H2", "H50", "H58"  },
      {  "H3", "H95", "H102" }
   }  ;
   
   /* Default to all CDRs                                               */
   firstLoop = 0;
   lastLoop  = NCDR;
   
   /* Update it we have specified to do only one chain                  */
   if(chain == 'L')
   {
      firstLoop = 0;
      lastLoop  = 3;
   }
   else if(chain == 'H')
   {
      firstLoop = 3;
      lastLoop  = NCDR;
   }

   for(loop=firstLoop; loop<lastLoop; loop++)
   {
      if((start = FindRes(Sequence, NRes, LoopDef[loop].start)) == (-1))
      {
         fprintf(stderr,"Warning (chothia): Unable to find residue %s \
in input\n", LoopDef[loop].start);
         fprintf(out,"CDR %s  Missing Residues\n",LoopDef[loop].name);
         continue;
      }

      if((stop = FindRes(Sequence, NRes, LoopDef[loop].stop)) == (-1))
      {
         fprintf(stderr,"Warning (chothia): Unable to find residue %s \
in input\n", LoopDef[loop].stop);
         fprintf(out,"CDR %s  Missing Residues\n",LoopDef[loop].name);
         continue;
      }
         
      len   = 1 + stop - start;

      if(loop==0 || loop==3) 
      {
         strncpy(cdr1, LoopDef[loop].name, SMALLWORD);
         cdr1len = len;
      }
      
      ReportACanonical(out, LoopDef[loop].name, len, Sequence, NRes,
                       verbose, cdr1, cdr1len);
   }
}


/************************************************************************/
/*>int FindRes(SEQUENCE *Sequence, int NRes, char *InRes)
   ------------------------------------------------------
   Input:   SEQUENCE   *Sequence      The sequence array
            int        NRes           Length of sequence array
            char       *InRes         The res number to find
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
   30.05.96 Checks for InRes being --- and returns -1
   19.12.08 Changed strcpy() to strncpy()
            Changed res[16] and buff[16] to use MAXWORD
*/
int FindRes(SEQUENCE *Sequence, int NRes, char *InRes)
{
   int  i,
        InsPos = (-1);
   char res[MAXWORD];

   if(!strncmp(InRes,"---",3))
      return(-1);

   strncpy(res,InRes,MAXWORD);
   
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
      char buff[MAXWORD];

      strncpy(buff, Sequence[i].resnum, MAXWORD);
      TERMALPHA(buff+1);

      if(!strcmp(buff, res))
         return(i);
   }
   
   /* Tried everything!                                                 */
   return(-1);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   16.05.95 Original    By: ACRM
   18.08.95 V1.1   
   30.11.95 V1.2
   08.05.96 V1.3 Added -n
   09.05.96 V1.4  
   30.05.96 V1.5
   19.12.08 V1.6
   13.02.14 V1.7
   09.08.15 V2.0
   09.08.15 V2.1 Added -L and -H
   14.12.16 V2.2 
   12.10.21 V2.3
*/
void Usage(void)
{
   fprintf(stderr,"\nChothia V2.3 (c) 1995-2021, Prof. Andrew C.R. \
Martin, UCL\n\n");

   fprintf(stderr,"Usage: chothia [-c filename] [-L|-H] [-v] [-n] \
[input.seq [output.dat]]\n");
   fprintf(stderr,"               -c Specify Chothia datafile (Default: \
chothia.dat)\n");
   fprintf(stderr,"               -L Input only contains light chain\n");
   fprintf(stderr,"               -H Input only contains heavy chain\n");
   fprintf(stderr,"               -v Verbose; give explanations when \
no canonical found\n");
   fprintf(stderr,"               -n The sequence file has Chothia \
(rather than Kabat) numbering\n");
   fprintf(stderr,"       I/O is through stdin/stdout if files are not \
specified.\n\n");

   fprintf(stderr,"Chothia is a program to assign canonical classes to \
an antibody sequence.\n");
   fprintf(stderr,"Input to the program is a listing of Kabat residue \
numbers and the\n");
   fprintf(stderr,"1-letter or 3-letter code name for the residue at \
each position. Such\n");
   fprintf(stderr,"a file may be generated from a PIR file using the \
program KabatSeq.\n");
   fprintf(stderr,"The numbering in this file is normally Kabat \
numbering; if the -n switch is\n");
   fprintf(stderr,"specified on the command line, the file must have \
Chothia numbering.\n\n");

   fprintf(stderr,"The program will look for the datafile first in the \
current directory\n");
   fprintf(stderr,"and then in the directory specified by the %s \
environment variable.\n", ENV_KABATDIR);
   fprintf(stderr,"This data file is also used by the KabatMan database \
software.\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *ChothiaFile, BOOL *verbose, char *chain)
   ---------------------------------------------------------------------
   Input:   int  argc             Argument count
            char **argv           Argument array
   Output:  char *infile          Input file (or blank string)
            char *outfile         Output file (or blank string)
            char *ChothiaFile     Chothia data file
            BOOL *verbose         Flag to show details of mismatches
            char *chain           Chain to  handle (default both)
   Returns: BOOL                  Success?
   Globals: BOOL gChothiaNumbered The sequence data is Chothia numbered

   Parse the command line
   
   16.05.95 Original    By: ACRM
   08.05.96 Added -n
   19.12.08 Changed strcpy() to strncpy()
   09.08.15 Added -l and -h for chain specification
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *ChothiaFile, BOOL *verbose, char *chain)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *verbose = FALSE;
   *chain   = ' ';

   gChothiaNumbered = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
            argc--;
            argv++;
            strncpy(ChothiaFile, argv[0], MAXBUFF);
            break;
         case 'v':
            *verbose = TRUE;
            break;
         case 'n':
            gChothiaNumbered = TRUE;
            break;
         case 'L':
            if(*chain != ' ')
               return(FALSE);
            *chain = 'L';
            break;
         case 'H':
            if(*chain != ' ')
               return(FALSE);
            *chain = 'H';
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
         strncpy(infile, argv[0], MAXBUFF);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strncpy(outfile, argv[0], MAXBUFF);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}

/************************************************************************/
/*>int TestThisCanonical(CHOTHIA *p, char *LoopName, int LoopLen,
                         SEQUENCE *Sequence, char *cdr1, int cdr1len)
   ------------------------------------------------------------------
   16.02.11 Extracted from ReportACanonical()
*/
int TestThisCanonical(CHOTHIA *p, char *LoopName, int LoopLen,
                      SEQUENCE *Sequence, int NRes, 
                      char *cdr1, int cdr1len)
{
   int  NMismatch = 10000, /* Return this if loop length/name wrong     */
        res,
        i;
   
   /* If the Loop name and length match                                 */
   if(!strcmp(p->LoopID, LoopName) && (p->length == LoopLen))
   {
      NMismatch = 0;  /* Assume we are OK                               */
         
      /* Check each residue specified by this canonical definition      */
      for(i=0; strcmp(p->resnum[i], "-1"); i++)
      {
         if(gCanonChothNum == gChothiaNumbered)
         {
            /* Both the datafile and the sequence data use the same
               numbering scheme (Kabat or Chothia)
            */
            res = FindRes(Sequence, NRes, p->resnum[i]);
         }
         else if(gCanonChothNum)
         {
            /* Datafile uses Chothia numbering while the sequence data
               uses Kabat numbering
            */
            res = FindRes(Sequence, NRes, 
                          ChoKab(cdr1, cdr1len, p->resnum[i]));
         }
         else
         {
            /* Datafile uses Kabat numbering while the sequence data
               uses Chothia numbering
            */
            res = FindRes(Sequence, NRes, 
                          KabCho(cdr1, cdr1len, p->resnum[i]));
         }

         /* This is a disallowed residue type, so increment the mismatch 
            counter
            30.05.96 Added check on -1
         */
         if((res==(-1)) || (!strchr(p->restype[i], Sequence[res].seq)))
         {
            NMismatch++;
         }
      }
   }

   return(NMismatch);
}

/************************************************************************/
/*>void ReportACanonical(FILE *out, char *LoopName, int LoopLen, 
                         SEQUENCE *Sequence, int NRes, BOOL verbose,
                         char *cdr1, int cdr1len)
   -----------------------------------------------------------------
   Input:   FILE     *out          Output file pointer
            char     *LoopName     Name of a loop (e.g. L1)
            int      LoopLen       Length of the loop
            SEQUENCE *Sequence     Sequence array
            int      NRes          Length of sequence
            BOOL     verbose       Flag to display reasons
            char     *cdr1         Name of CDR1 (L1 or H1)
            char     *cdr1len      Length of CDR1
   Returns: BOOL                   Success?

   Reports the canonical class for an individual loop

   16.05.95 Original    By: ACRM
   17.05.95 Only prints source data if verbose
   08.05.96 Converts between Chothia and Kabat numbering if required
   09.05.96 Fixed bug in reporting mismatches with numbering conversion
   30.05.96 Mismatches report numbering scheme in data file
            Added check for -1 return from FindRes(); reports deleted
            residues
   16.02.11 Moved actual canonical finding code out into 
            TestThisCanonical()
   17.02.11 Re-written to deal with priority chains
*/
void ReportACanonical(FILE *out, char *LoopName, int LoopLen, 
                      SEQUENCE *Sequence, int NRes, BOOL verbose,
                      char *cdr1, int cdr1len)
{
   CHOTHIA *p,
           *theMatch,
           *best = NULL;
   int     i,
           res,
           NMismatch,
           MinMismatch = 10000;
   
   /* Run through the linked list of Canonical definitions              */
   for(p=gChothia; p!=NULL; NEXT(p))
   {
      /* If this is subordinate to something else                       */
      if(p->nsubordinate)
      {
         /* If it hasn't also got things it has priority over (i.e. it's
            not in the middle of a priority chain
         */
         if(p->npriority == 0)
         {
            CHOTHIA *q = p;
         
            /* Walk to the highest priority class                       */
            while(q->subordinate_to != NULL)
            {
               NEXT(q);
            }
            /* q now points to the highest priority class, walk back to 
               the lowest priority class, testing for a perfect match
            */
            do {
               NMismatch = TestThisCanonical(q, LoopName, LoopLen, 
                                             Sequence, NRes, 
                                             cdr1, cdr1len);
               if(NMismatch == 0)
               {
                  break;
               }
               q=q->priority_over;
            }  while(q != NULL);
            
            /* If we found a match then use that, otherwise, use p      
               In other words we only accept mismatches against the 
               lowest priority class.
             */
            theMatch = (q==NULL)?p:q;
         }
      }
      else
      {
         theMatch = p;
         NMismatch = TestThisCanonical(p, LoopName, LoopLen, Sequence, 
                                       NRes, cdr1, cdr1len);
      }

      if(NMismatch == 0)  /* We've found the canonical                  */
      {
         break;
      }
      else
      {
         if(NMismatch < MinMismatch)
         {
            MinMismatch = NMismatch;
            best = theMatch;
         }
      }
   }

   if(NMismatch == 0)
   {
      fprintf(out,"CDR %s  Class %-3s", LoopName, theMatch->class);
      if(verbose && strlen(theMatch->source))
         fprintf(out," %s", theMatch->source);
      fprintf(out,"\n");
   }
   else
   {
      fprintf(out,"CDR %s  Class ?  \n", LoopName); 
   
      if(verbose)
      {
         if(best==NULL)
         {
            fprintf(out, "! No canonical of the same loop length\n");
         }
         else
         {
            fprintf(out, "! Similar to class %s, but:\n", best->class);

            /* Display each mismatch for this canonical definition      */
            for(i=0; strcmp(best->resnum[i], "-1"); i++)
            {
               if(gCanonChothNum == gChothiaNumbered)
               {
                  /* Both the datafile and the sequence data use the same
                     numbering scheme (Kabat or Chothia)
                  */
                  res = FindRes(Sequence, NRes, best->resnum[i]);
               }
               else if(gCanonChothNum)
               {
                  /* Datafile uses Chothia numbering while the sequence 
                     data uses Kabat numbering
                  */
                  res = FindRes(Sequence, NRes, 
                                ChoKab(cdr1, cdr1len, best->resnum[i]));
               }
               else
               {
                  /* Datafile uses Kabat numbering while the sequence data
                     uses Chothia numbering
                  */
                  res = FindRes(Sequence, NRes, 
                                KabCho(cdr1, cdr1len, best->resnum[i]));
               }

               /* 30.05.96 Added check on -1                            */
               if(res==(-1))
               {
                  fprintf(out, "!    %s (%s Numbering) is deleted.\n", 
                          best->resnum[i],
                          (gCanonChothNum?"Chothia":"Kabat"));
               }
               else if(!strchr(best->restype[i], Sequence[res].seq))
               {
                  fprintf(out, "!    %s (%s Numbering) = %c \
(allows: %s)\n", 
                          best->resnum[i],
                          (gCanonChothNum?"Chothia":"Kabat"),
                          Sequence[res].seq,
                          best->restype[i]);
               }
            }
         }
      }
   }   
}


