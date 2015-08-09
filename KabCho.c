/*************************************************************************

   Program:    
   File:       KabCho.c
   
   Version:    V1.2
   Date:       30.05.96
   Function:   Convert Kabat antibody numbering to Chothia numbering
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 1996
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
   EMail:      martin@biochem.ucl.ac.uk
               andrew@stagleys.demon.co.uk
               
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
   V1.0  07.05.96 Original
   V1.1  08.05.96 Corrected 10-residue L1s
   V1.2  30.05.96 Changed variable name in ChoKab()

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>

/************************************************************************/
/* Defines and macros
*/
#define MAXLEN_L1 17
#define MAXLEN_H1 12

/************************************************************************/
/* Globals
*/

/* The following are conversion tables between Chothia (in row 0)
   and Kabat numbering. The row number corresponds to the loop
   length.

   For short loops where residues other than the insertion code
   residues are deleted, we currently assume that the output ID is the
   same as the input; this can be changed by modifying the appropriate
   rows of the tables below.

   Updated for 10-residue CDR-L1

   Note that the second array dimensions must be altered if the loop
   definitions are ever changed.
*/
static char *sL1Table[MAXLEN_L1+1][18] = 
{{"L24","L25","L26","L27","L28", "L29", "L30", "L30A","L30B","L30C",
     "L30D","L30E","L30F","L31","L32","L33","L34",NULL}, /* Chothia     */
 {"L24","L25","L26","L27","L28", "L29", "----","----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length =  1 */
 {"L24","L25","L26","L27","L28", "L29", "----","----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length =  2 */
 {"L24","L25","L26","L27","L28", "L29", "----","----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length =  3 */
 {"L24","L25","L26","L27","L28", "L29", "----","----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length =  4 */
 {"L24","L25","L26","L27","L28", "L29", "----","----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length =  5 */
 {"L24","L25","L26","L27","L28", "L29", "----","----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length =  6 */
 {"L24","L25","L26","L27","L28", "L29", "----","----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length =  7 */
 {"L24","L25","L26","L27","L28", "L29", "----","----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length =  8 */
 {"L24","L25","L26","L27","L29", "---", "----","----","----","----",
     "----","----","----","---","L32","L33","L34",NULL}, /* Length =  9 */
 {"L24","L25","L26","L27","L29", "L30", "L31", "----","----","----",
     "----","----","----","---","L32","L33","L34",NULL}, /* Length = 10 */
 {"L24","L25","L26","L27","L28", "L29", "L30", "----","----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length = 11 */
 {"L24","L25","L26","L27","L27A","L28", "L29", "L30", "----","----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length = 12 */
 {"L24","L25","L26","L27","L27A","L27B","L28", "L29", "L30", "----",
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length = 13 */
 {"L24","L25","L26","L27","L27A","L27B","L27C","L28", "L29", "L30", 
     "----","----","----","L31","L32","L33","L34",NULL}, /* Length = 14 */
 {"L24","L25","L26","L27","L27A","L27B","L27C","L27D","L28", "L29", 
     "L30", "----","----","L31","L32","L33","L34",NULL}, /* Length = 15 */
 {"L24","L25","L26","L27","L27A","L27B","L27C","L27D","L27E","L28" ,
     "L29", "L30", "----","L31","L32","L33","L34",NULL}, /* Length = 16 */
 {"L24","L25","L26","L27","L27A","L27B","L27C","L27D","L27E","L27F",
     "L28", "L29", "L30", "L31","L32","L33","L34",NULL}  /* Length = 17 */
};
   
static char *sH1Table[MAXLEN_H1+1][13] = 
{{"H26","H27","H28","H29","H30","H31","H31A","H31B","H32","H33","H34",
     "H35",NULL},                                        /* Chothia     */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  1 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  2 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  3 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  4 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  5 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  6 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  7 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  8 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length =  9 */
 {"H26","H27","H28","H29","H30","H31","----","----","H32","H33","H34",
     "H35",NULL},                                        /* Length = 10 */
 {"H26","H27","H28","H29","H30","H31","H32", "----","H33","H34","H35",
     "H35A",NULL},                                       /* Length = 11 */
 {"H26","H27","H28","H29","H30","H31","H32", "H33", "H34","H35","H35A",
     "H35B",NULL}                                        /* Length = 12 */
};

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>char *KabCho(char *cdr, int length, char *kabspec)
   --------------------------------------------------
   Input:   char   *cdr      The CDR-ID (L1, L2, L3, H1, H2 or H3)
            int    length    The CDR length (AbM length i.e. normal length
                             for CDR-H1 is 10 residues)
            char   *kabspec  The Kabat residue spec (e.g. L27A)
   Returns: char   *         The equivalent Chothia residue spec

   Given a CDR-ID (L1, L2, L3, H1, H2 or H3), the length of the CDR and
   a Kabat residue spec. (e.g. L27A), returns the equivalent Chothia
   residue spec. Note that the supplied loop length must be AbM style
   loop definition.

   07.05.96 Original   By: ACRM
*/
char *KabCho(char *cdr, int length, char *kabspec)
{
   int i;

   if(!strncmp(cdr,"L1",2))
   {
      if(length <= MAXLEN_L1)
      {
         for(i=0; sL1Table[length][i] != NULL; i++)
         {
            if(!strcmp(kabspec,sL1Table[length][i]))
            {
               return(sL1Table[0][i]);
            }
         }
      }
   }
   else if(!strncmp(cdr,"H1",2))
   {
      if(length <= MAXLEN_H1)
      {
         for(i=0; sH1Table[length][i] != NULL; i++)
         {
            if(!strcmp(kabspec,sH1Table[length][i]))
            {
               return(sH1Table[0][i]);
            }
         }
      }
   }

   /* If we get here, it's one of the other loops where the specs are the
      same, or this residue wasn't found in the tables, or the loop was
      too long. Just return the specification string as entered.
   */
   return(kabspec);
}


/************************************************************************/
/*>char *ChoKab(char *cdr, int length, char *chospec)
   --------------------------------------------------
   Input:   char   *cdr      The CDR-ID (L1, L2, L3, H1, H2 or H3)
            int    length    The CDR length (AbM length i.e. normal length
                             for CDR-H1 is 10 residues)
            char   *chospec  The Chothia residue spec (e.g. L30A)
   Returns: char   *         The equivalent Kabat residue spec

   Given a CDR-ID (L1, L2, L3, H1, H2 or H3), the length of the CDR and
   a Chothia residue spec. (e.g. L30A), returns the equivalent Kabat
   residue spec. Note that the supplied loop length must be AbM style
   loop definition.

   07.05.96 Original   By: ACRM
   30.05.96 Changed variable name kabspec to chospec; it was very
            confusing before...
*/
char *ChoKab(char *cdr, int length, char *chospec)
{
   int i;

   if(!strncmp(cdr,"L1",2))
   {
      if(length <= MAXLEN_L1)
      {
         for(i=0; sL1Table[0][i] != NULL; i++)
         {
            if(!strcmp(chospec,sL1Table[0][i]))
            {
               return(sL1Table[length][i]);
            }
         }
      }
   }
   else if(!strncmp(cdr,"H1",2))
   {
      if(length <= MAXLEN_H1)
      {
         for(i=0; sH1Table[0][i] != NULL; i++)
         {
            if(!strcmp(chospec,sH1Table[0][i]))
            {
               return(sH1Table[length][i]);
            }
         }
      }
   }

   /* If we get here, it's one of the other loops where the specs are the
      same, or this residue wasn't found in the tables, or the loop was
      too long. Just return the specification string as entered.
   */
   return(chospec);
}
