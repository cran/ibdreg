/* $Author: sinnwell $ */
/* $Date: 2006/11/28 21:20:19 $ */
/* $Header: /people/biostat3/sinnwell/Projects/IBDReg/Make/RCS/sim.mark.prop.c,v 1.3 2006/11/28 21:20:19 sinnwell Exp $ */
/* $Locker:  $ */
/*
 * $Log: sim.mark.prop.c,v $
 * Revision 1.3  2006/11/28 21:20:19  sinnwell
 * define long to unsigned int
 *
 * Revision 1.2  2005/10/04 18:42:48  folie
 * Edited comments
 * 
 * Revision 1.1  2005/09/28 18:37:57  folie
 * Initial revision
 * * 
 */
#include "mcmemory.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* JPS, building the package for R 2.15.2 11/1/2012:
   note that printped .c should be copied into here if we want to be able 
   to print a whole pedigree. It is left out of the R package
   because of the printf statements, which could be changed to REprintf, but
   that requires including R.h, which causes many problems with the free and
   alloc functions that Dan Folie defines here instead of there
*/

/* This file contains the C routines called by the sim.mark.prop() S-PLUS     */
/* function                                                                   */

#define long unsigned

#define MALE           1
#define FEMALE         0
static long REACHED  = 1;
static long COMPUTED = 1;
static long FLAG     = 0;
struct pedigreeData *pedData = NULL;

/* Represents a person in a generated pedigree                                */

static struct person
{ struct marriage  *parents;              /* Parents' marriage node           */
  struct marriage  *nuclearFamily;        /* Person's first marriage node     */
  long              sex;		  /* M or F depending on sex          */
  long              id;			  /* Person ID                        */
  long              chrom1;               /* Markers for 1st chromosome       */
  long              chrom2;               /* Markers for 2nd chromosome       */
  long              traverseStatus;       /* Has this node been reached       */
  long              computeStatus;        /* Has this node had data generated */
};

/* Represents a marriage in a generated pedigree                              */

static struct marriage
{ struct person    *husband;	          /* Husband in this marriage         */
  struct person    *wife;                 /* Wife in this marriage	      */
  struct children  *child;                /* Children from this marriage      */
  struct marriage  *nextHusbandMarriage;  /* Husband's other marriage(s)      */
  struct marriage  *nextWifeMarriage;     /* Wife's other marriage(s)	      */
};

/* Represents a linked list of children in a marriage                         */

static struct children
{ struct person    *firstChild;		  /* First child		      */
  struct children  *nextChild;		  /* Next child			      */
};

/* Represents a pedigree                                                      */

static struct pedigree
{ struct person    *proband;              /* Proband of the pedigree          */
  long              n;                    /* Number of people in the pedigree */
  struct person   **personNodes;          /* People in the pedigree           */
};

/* Represents a conglomeration of all the data required/generated for a       */
/* pedigree                                                                   */

static struct pedigreeData
{ long     *father;                  /* Father ID's                           */
  long     *mother;		     /* Mother ID's                           */
  long     *person;                  /* Person ID's                           */
  long     *sex;                     /* Sex status                            */
  long      n;                       /* Length of the above five vectors      */
  long      nMark;                   /* Number of marker positions            */
  long      XProbandIndex;           /* Proband index in X covariate matrix   */
  long    **CM;                      /* n-by-2*nMark matrix of marker data    */
  long     *cm;                      /* The vector form of CM passed from S   */
  struct pedigree *ped;              /* Points to current pedigree structure  */
};

static struct pedigree *generate_pedigree(struct pedigreeData *pedData);
static void generate_markers_autosomal(struct person *subject,
                                       struct pedigreeData *pedData);
static void generate_markers_xlinked(struct person *subject,
                                     struct pedigreeData *pedData);
static void traverse(struct pedigreeData *pedData, void (*func)());
static void traverse_engine(struct person *subject, 
                            struct pedigreeData *pedData,
                            void(*func)());
static long runifAS183_seed(long iseed1, long iseed2, long iseed3);
static double runifAS183();
static void qs(double array_sort[], long left, long right);
static long rpoisson(double mu);

/* SIM_MARK_PROP_GEN_PED

                                   FUNCTION ARGUMENTS                         

sim_mark_prop_gen_ped() takes the following arguments:

-father:          An array of long integers such that fatherID[i] denotes the 
                  person that is the father of the person given by personID[i].  
                  An entry of zero denotes that the person identifed by 
                  personID[i] has no father in this pedigree.

-mother:          An array of long integers such that motherID[i] denotes the 
                  person that is the mother of the person given by personID[i].  
                  An entry of zero denotes that the person identifed by 
                  personID[i] has no mother in this pedigree.

-person:          An array of long such that for each integer i between 1 and n,
                  inclusive, we have personID[i] as the logical identifier for 
                  person i in the pedigree.

-sex:             An array of long such that sex[i] denotes the sex of the 
                  person identified by personID[i].  A value of 1 denotes 
                  male, and a value of zero denotes female.

-n:               A long integer pointer such that *n gives the length of each 
                  of the preceding n-by-1 arrays. 

-nMark:           A long pointer such that *nMark is the number of markers
                  per chromosome fragment and 2 * *nMark is the number of 
                  markers per pair of chromosome fragments (i.e. per person).

-cm:              A long pointer to a vector representing the column-major
                  form of the n-by-2*nMark data matrix containing simulated 
                  marker segregation data for the generated pedigree.

-iseed1:          A long integer pointer to the first of three seeds used by
                  the Wichman and Hill uniform random number generator.

-iseed2:          A long integer pointer to the second of three seeds used by
                  the Wichman and Hill uniform random number generator.

-iseed3:          A long integer pointer to the third of three seeds used by
                  the Wichman and Hill uniform random number generator.



                                   SIDE EFFECTS

The global struct pedigreeData *pedData data structure will be constructed.

                                   RETURN VALUE

sim_mark_prop_gen_ped() returns no value.

                                   SYNOPSIS

sim_mark_prop_gen_ped() generates the pedigree specified by its five 
arguments and instantiates non-founders with there marker data (generated in 
S-Plus).  The random seeds for the runifAS183() random number generator are
also set.

*/

void sim_mark_prop_gen_ped(long   *person,
                           long   *father,    
                           long   *mother,
                           long   *sex,
                           long   *n,
                           long   *nMark,
                           long   *cm,
                           long   *iseed1,
                           long   *iseed2,
                           long   *iseed3)
{
 struct person *subject;
 long i, j, row, col;

 /* Set the seeds for the uniform random number generator:                   */  
 runifAS183_seed(*iseed1, *iseed2, *iseed3);

 /* Store input parameters into a single struct which can then be passed to  */
 /* all called functions - simplifies the calling conventions.               */
 pedData = (struct pedigreeData *)Salloc(1,struct pedigreeData);
 pedData->person = person; 
 pedData->father = father;
 pedData->mother = mother;
 pedData->sex    = sex;
 pedData->n      = *n;
 pedData->nMark  = *nMark;
 pedData->cm     = cm;

 /* Generate the pedigree and instantiate founder markers: */
 generate_pedigree(pedData);

 return;
}

/* SIMULATE_MARKER_PROPAGATION

                                   FUNCTION ARGUMENTS                         

simulate_marker_propagation() takes the following arguments:


-numIter:         A long integer pointer such that *numIter gives the number of
                  times that markers are propagated throughout the pedigree.

-ret:             A long pointer to a vector representing the column-major
                  form of the n-by-2*numIter data matrix returned to S-PLUS 
                  containing simulated marker segregation data for the 
                  generated pedigree.

-proband:         A long pointer listing the proband to start traversing from.
                  The starting proband should not affect the inheritance
                  performed (although this is hard to tell due to the use of 
                  a random number generator in the inheritance computations). 
                  This variable is primarily used for testing/debugging 
                  purposes.

-xlinked:         A long pointer giving a logical value (zero/nonzero) that
                  choses whether marker inheritance will be performed as for
                  X-linked markers or for autosomal markers.

                                   SIDE EFFECTS

The variable ret will be modified by reference.  Also, the static global 
variable FLAG will be toggled to the value (FLAG+1) mod 2.

                                   RETURN VALUE

simulate_marker_propagation() returns no value.

                                   SYNOPSIS

simulate_marker_propagation() proceeds to traverse the pedigree and 
propagates markers from all pedigree founders to non-founders.  After each
traversal, the generated markers are stored in the corresponding entries 
(columns) of ret.

*/

void simulate_marker_propagation(long *numIter,
                                 long *ret,
                                 long *proband,
                                 long *xlinked)
{
 struct person *subject;
 long i, j, n, nrow, col;

 n = *numIter;
 nrow = pedData->n;

 if(*xlinked)
   {
    for(i = 0; i < n; i++)
       {
        /* Segregate the markers from founders to non-founders (X-linked):  */
        if((0 < *proband) && (*proband <= pedData->n) )
            pedData->ped->proband = pedData->ped->personNodes[*proband];
         traverse(pedData, generate_markers_xlinked);

         /*** Copy data back to S-Plus ***/
         for(j = 1; j <= nrow; j++)
           {
            subject                 = pedData->ped->personNodes[j];
            ret[2*i*nrow + j-1]     = subject->chrom1;
            ret[(2*i+1)*nrow + j-1] = subject->chrom2;
          } /*** end for(j)       ***/
        }  /**** end for(i)       ***/
   }      /***** end if(*xlinked) ***/
 else    /****** (!(*xlinked))    ***/
   {
    for(i = 0; i < n; i++)
       {
        /* Segregate the markers from founders to non-founders (autosomal): */
        if((0 < *proband) && (*proband <= pedData->n) )
            pedData->ped->proband = pedData->ped->personNodes[*proband];
         traverse(pedData, generate_markers_autosomal);

         /*** Copy data back to S-Plus ***/
         for(j = 1; j <= nrow; j++)
           {
            subject                 = pedData->ped->personNodes[j];
            ret[2*i*nrow + j-1]     = subject->chrom1;
            ret[(2*i+1)*nrow + j-1] = subject->chrom2;
          } /*** end for(j)            ***/
        }  /**** end for(i)            ***/
   }      /***** end else(!(*xlinked)) ***/
 return;
}

/* SIM_MARK_PROP_FREE_MEM

                                   FUNCTION ARGUMENTS                         

sim_mark_prop_free_mem() takes no arguments.

                                   SIDE EFFECTS

The variable ret will be modified by reference.  Also, the static global 
variable FLAG will be toggled to the value (FLAG+1) mod 2.

                                   RETURN VALUE

sim_mark_prop_free_mem() returns no value.

                                   SYNOPSIS

sim_mark_prop_free_mem() frees dynamically allocated memory


*/

void sim_mark_prop_free_mem()
{
 mc_free_all();
 return;
}

/* GENERATE_PEDIGREE

                                   FUNCTION ARGUMENTS                         

generate_pedigree() takes a pointer to a struct of type pedigreeData called 
"pedData" and the following fields of pedData are used:

pedData->father:   An array of long integers such that fatherID[i] denotes 
                   the person that is the father of the person given by 
                   personID[i].  An entry of zero denotes that the person 
                   identifed by personID[i] has no father in this pedigree.

pedData->mother:   An array of long integers such that motherID[i] denotes 
                   the person that is the mother of the person given by 
                   personID[i].  An entry of zero denotes that the person 
                   identifed by personID[i] has no mother in this pedigree.

pedData->person:   An array of long such that for each integer i between 1 
                   and n, inclusive, we have personID[i] as the logical 
                   identifier for person i in the pedigree.

pedData->sex:      An array of long such that sex[i] denotes the sex of the 
                   person identified by personID[i].  A value of 1 denotes 
                   male, and a value of zero denotes female.

pedData->n:        A long integer pointer such that *n gives the length of 
                   each of the preceding five arrays. 

pedData->nMark:    A long pointer such that *nMark is the number of markers
                   per chromosome fragment and 2 * *nMark is the number of 
                   markers per pair of chromosome fragments (i.e. per person).

pedData->cm:       A double pointer to a vector representing the column-major
                   form of the n-by-2*nMark data matrix containing simulated 
                   marker segregation data for the generated pedigree.

                                   SIDE EFFECTS

The pedData->CM, pedData->XProbandIndex, and pedData->ped fields will be set 
to the the marker matrix, the index of the row in CM corresponding to the male
founder in the pair of founders in the top-most pedigree generation, and 
pedigree that are each constructed by generate_pedigree().  None the the seven
struct fields listed above will be modifed.

Note:  This pedigree generator is for use with both cyclic and acyclic pedigrees
with one or more pairs of ancestoral founders.

                                   RETURN VALUE

generate_pedigree() returns a pointer to a struct of type pedigree.
*/

static struct pedigree *generate_pedigree(struct pedigreeData *pedData)
{
 long             *father   = pedData->father;
 long             *mother   = pedData->mother; 
 long             *person   = pedData->person;
 long             *sex      = pedData->sex;
 long              n        = pedData->n;
 long              nMark    = pedData->nMark;
 long             *cm       = pedData->cm;
 long              i, j, action;
 long            **CM;
 struct person   **personNodes;
 struct marriage  *marriageNode;
 struct children  *childNode;
 struct person    *currentPerson;
 struct person    *currentFather;
 struct person    *currentMother;
 struct marriage  *currentSpouse;
 struct marriage  *currentMarriage;
 struct marriage  *currentFatherMarriage;
 struct marriage  *currentMotherMarriage;
 struct children  *currentChild;
 struct pedigree  *Pedigree;

/* Set up the n-by-2*nMark chromosome marker matrix CM    */
 CM = (long **) Salloc(n+1, long *);
 for(i = 0; i <= n; i++)
    CM[i] = (long *) Salloc(2*nMark, long);

/* Fill in the data from S-Plus for the founders          */
 for(i = 1; i <= n; i++)
   for(j = 0; j < 2*nMark; j++)
      CM[i][j] = cm[n*j + i - 1];
    
 pedData->CM = CM;

/* Set up person nodes and fill in basic information      */
 personNodes = (struct person **) Salloc(n+1, struct person *);
 personNodes[0] = NULL;
 for(i=1; i <= n; i++)
   {
    personNodes[i] = (struct person *) Salloc(1L, struct person);
    currentPerson                 = personNodes[i];
    currentPerson->id             = person[i-1];
    currentPerson->sex            = sex[i-1];
    currentPerson->traverseStatus = REACHED;
    currentPerson->computeStatus  = COMPUTED;
    currentPerson->parents        = NULL;
    currentPerson->nuclearFamily  = NULL;
    currentPerson->chrom1         = *(CM[i]);
    currentPerson->chrom2         = *(CM[i] + nMark);
   }


/**********************Link person nodes into a pedigree***********************/

 for(i = 1; i <= n; i++)
   {
    currentFather = personNodes[father[i-1]];
    currentMother = personNodes[mother[i-1]];
    currentPerson = personNodes[i];

/* The variable "action" may take on the values 0,1,2,...,6.  For a given     */
/* person, these values are interpreted as follows:                           */
/*                                                                            */
/* 0: The current person has no parents in this pedigree.                     */
/* 1: This is the first marriage and child for both parents.                  */
/* 2: This is the first marriage for the father, the mother has been          */
/*    previously married, and this is the first child in this marriage.       */
/* 3: This is the first marriage for the mother, the father has been          */
/*    previously married, and this is the first child in this marriage.       */
/* 4: The father and mother are already married to each other, each parent's  */
/*    first marriage points to this marriage, and this is the kth child in    */
/*    this marriage, with k >= 2.                                             */
/*                                                                            */
/* For value 5, the father and mother have each been married before and at    */
/* least one of them, possibly both, has multiple spouses.                    */
/*                                                                            */
/* 5: At least one of the mother and father has multiple marriages.  This     */
/*    may be the first child in this marriage, if the marriage doesn't        */
/*    already exist, or the kth child in an existing marriage, k >= 2.        */



  if(   (currentFather != NULL) &&
        (currentMother != NULL)   )
           if(   (currentFather->nuclearFamily == NULL) &&
                 (currentMother->nuclearFamily == NULL)   )
            action = 1;
           else if(   (currentFather->nuclearFamily == NULL) &&
                      (currentMother->nuclearFamily != NULL)   ) 
                action = 2;
           else if(   (currentFather->nuclearFamily != NULL) &&
                      (currentMother->nuclearFamily == NULL)   ) 
                action = 3;
/**************  We now know that the following is true:         **************/
/**************  (   (currentFather->nuclearFamily != NULL) &&   **************/
/**************      (currentMother->nuclearFamily != NULL)   )  **************/
           else if(currentFather->nuclearFamily == currentMother->nuclearFamily)
	        action = 4;
/*************  We now know that at least one of the father and mother    *****/
/*************  has more than one marriage:                               *****/
           else 
	        action = 5;
  else
        action = 0;

 switch(action)
   {
    case 0:
      break;

    case 1:
      {
       marriageNode = (struct marriage *)Salloc(1L, struct marriage);
       marriageNode->husband = currentFather;
       marriageNode->wife = currentMother;
       marriageNode->nextHusbandMarriage = NULL;
       marriageNode->nextWifeMarriage = NULL;

       childNode = (struct children *)Salloc(1L, struct children);
       childNode->firstChild = currentPerson;
       childNode->nextChild = NULL;
       marriageNode->child = childNode;
             
       currentFather->nuclearFamily = marriageNode;
       currentMother->nuclearFamily = marriageNode;
       currentPerson->parents = marriageNode;
      }  
       break;

    case 2: 
      {
       marriageNode = (struct marriage *)Salloc(1L, struct marriage);
       marriageNode->husband = currentFather;
       marriageNode->wife = currentMother;
       marriageNode->nextHusbandMarriage = NULL;
       marriageNode->nextWifeMarriage = NULL;
 
       childNode = (struct children *)Salloc(1L, struct children);
       childNode->firstChild = currentPerson;
       childNode->nextChild = NULL;
       marriageNode->child = childNode;
       
       currentFather->nuclearFamily = marriageNode;
       /*** Put this marriage at the end of the mother's marriage list ***/
       currentSpouse = currentMother->nuclearFamily;
       while(currentSpouse->nextWifeMarriage != NULL)
         { 
	   currentSpouse = currentSpouse->nextWifeMarriage; 
         }
       currentSpouse->nextWifeMarriage = marriageNode;
       currentPerson->parents = marriageNode;
      }
       break;
 
    case 3:
      {
       marriageNode = (struct marriage *)Salloc(1L, struct marriage);
       marriageNode->husband = currentFather;
       marriageNode->wife = currentMother;
       marriageNode->nextHusbandMarriage = NULL;
       marriageNode->nextWifeMarriage = NULL;
 
       childNode = (struct children *)Salloc(1L, struct children);
       childNode->firstChild = currentPerson;
       childNode->nextChild = NULL;
       marriageNode->child = childNode;
         
       /*** Put this marriage at the end of the father's marriage list ***/
       currentSpouse = currentFather->nuclearFamily;
       while(currentSpouse->nextHusbandMarriage != NULL)
         { 
	   currentSpouse = currentSpouse->nextHusbandMarriage; 
         }
       currentSpouse->nextHusbandMarriage = marriageNode;
       currentMother->nuclearFamily = marriageNode;
       currentPerson->parents = marriageNode;
      }
       break;

    case 4:
      {
       childNode = (struct children *)Salloc(1L, struct children);
       childNode->firstChild = currentPerson;
       childNode->nextChild = NULL;
       currentChild = currentFather->nuclearFamily->child;
       while(currentChild->nextChild != NULL)
         { 
          currentChild = currentChild->nextChild; 
         }
       currentChild->nextChild = childNode;   
       currentPerson->parents = currentFather->nuclearFamily;
      }
       break;

    case 5:
      { 
       currentMarriage = NULL;
       currentFatherMarriage = currentFather->nuclearFamily;
       while(currentFatherMarriage != NULL && currentMarriage == NULL)
         {
          currentMotherMarriage = currentMother->nuclearFamily;
          while(currentMotherMarriage != NULL && currentMarriage == NULL)
            {
	      if(currentFatherMarriage == currentMotherMarriage)
                 currentMarriage = currentFatherMarriage;
              currentMotherMarriage = currentMotherMarriage->nextWifeMarriage;
            }
          currentFatherMarriage = currentFatherMarriage->nextHusbandMarriage;
         }

       if(currentMarriage != NULL)
         {
	  /*** Marriage already exists - add the child: ***/
          childNode = (struct children *)Salloc(1L, struct children);
          childNode->firstChild = currentPerson;
          childNode->nextChild = NULL;
          currentChild = currentMarriage->child;
          while(currentChild->nextChild != NULL)
            { 
             currentChild = currentChild->nextChild; 
            }
          currentChild->nextChild = childNode;   
          currentPerson->parents = currentMarriage;
         }
       else
         {
	  /*** Marriage doesn't exist yet- create it: ***/
          marriageNode = (struct marriage *)Salloc(1L, struct marriage);
          marriageNode->husband = currentFather;
          marriageNode->wife = currentMother;
          marriageNode->nextHusbandMarriage = NULL;
          marriageNode->nextWifeMarriage = NULL;

	  /*** Add the child: ***/
          childNode = (struct children *)Salloc(1L, struct children);
          childNode->firstChild = currentPerson;
          childNode->nextChild = NULL;
          marriageNode->child = childNode;

	  /*** Put this marriage at the end of the father's marriage list ***/
          currentFatherMarriage = currentFather->nuclearFamily;
          while(currentFatherMarriage->nextHusbandMarriage != NULL)
	    currentFatherMarriage = currentFatherMarriage->nextHusbandMarriage;
          currentFatherMarriage->nextHusbandMarriage = marriageNode;

	  /*** Put this marriage at the end of the mother's marriage list ***/
          currentMotherMarriage = currentMother->nuclearFamily;
          while(currentMotherMarriage->nextWifeMarriage != NULL)
	    currentMotherMarriage = currentMotherMarriage->nextWifeMarriage;   
          currentMotherMarriage->nextWifeMarriage = marriageNode;

          currentPerson->parents = marriageNode;
         }
       break;
      }

   default:
     {}
    }
  }

/***********************Set up the pedigree return value***********************/

 Pedigree = (struct pedigree *)Salloc(1L, struct pedigree);
 Pedigree->n = n;
 Pedigree->personNodes = personNodes;

/*******Set the proband for the male founder at the top of the pedigree********/
 for(i = 1; i <=n; i++)
   {
    currentPerson = personNodes[i];
    if(currentPerson->parents != NULL)
      if(currentPerson->parents->husband->parents == NULL)
        if(currentPerson->parents->wife->parents == NULL)
          {
	    Pedigree->proband = currentPerson->parents->husband;
            pedData->XProbandIndex = Pedigree->proband->id;
            break;
          }
   }

/**********************Set pedData->ped to the above pedigree******************/

 pedData->ped = Pedigree;
 return Pedigree;
}

/* GENERATE_MARKERS_AUTOSOMAL

                                   FUNCTION ARGUMENTS                         

generate_markers_autosomal() takes a struct person pointer subject that points
to the person whose marker data will be populated and a pointer to a struct 
pedigreeData data struct.

                                   SIDE EFFECTS

The "chrom1" and "chrom2" fields of non-founder subjects will be overwritten
with an inherited marker

                                   SYNOPSIS

This function is passed to the pedigree traversal algorithm and propogates 
autosomal markers throughout a pedigree.

                                   RETURN VALUE

generate_markers_autosomal() returns no value.

*/

static void generate_markers_autosomal(struct person *subject, 
                                       struct pedigreeData *pedData)
{
  struct person *father, *mother;
  double u;

/****** Compute chromosome markers inherited from the child's parents ********/

  /*** Check if this is a founder: ***/
  if(subject->parents == NULL)
    return;

  /*** Simulate marker inheritance from the father: ***/
  father = subject->parents->husband;
  u = runifAS183();
  if(u < 0.5)
    subject->chrom1 = father->chrom1;
  else /*** (u >= 0.5) ***/
    subject->chrom1 = father->chrom2;

  /*** Simulate marker inheritance from the mother: ***/
  mother = subject->parents->wife;
  u = runifAS183();
  if(u < 0.5)
    subject->chrom2 = mother->chrom1;
  else /*** (u >= 0.5) ***/
    subject->chrom2 = mother->chrom2;

  return;
}

/* GENERATE_MARKERS_XLINKED

                                   FUNCTION ARGUMENTS                         

generate_markers_xlinked() takes a struct person pointer subject that points
to the person whose marker data will be populated and a pointer to a struct 
pedigreeData data struct.

                                   SIDE EFFECTS

The "chrom1" and "chrom2" fields of non-founder subjects will be overwritten
with an inherited marker

                                   SYNOPSIS

This function is passed to the pedigree traversal algorithm and propogates 
X-linked markers throughout a pedigree.

                                   RETURN VALUE

generate_markers_xlinked() returns no value.

*/

static void generate_markers_xlinked(struct person *subject, 
                                     struct pedigreeData *pedData)
{
  struct person *father, *mother;
  double u;

/****** Compute chromosome markers inherited from the child's parents ********/

  /*** Check if this is a founder: ***/
  if(subject->parents == NULL)
    return;

  if(subject->sex == MALE)
    {
     /*** Simulate X-linked marker inheritance from the mother, storing the ***/ 
     /*** inherited marker from the mother in both marker locations to keep ***/
     /*** the male homozygous at the X-linked locus:                        ***/
     mother = subject->parents->wife;
     u = runifAS183();
     if(u < 0.5)
       {
        subject->chrom1 = mother->chrom1;
        subject->chrom2 = mother->chrom1;
       }   /*** end if(u < 0.5) ***/
     else /**** (u >= 0.5)      ***/
       {
        subject->chrom1 = mother->chrom2;
        subject->chrom2 = mother->chrom2;
       } /*** end else (u >= 0.5)          ***/
    }   /**** end if(subject->sex == MALE) ***/
  else /***** (subject ->sex == FEMALE)    ***/
    {
     /*** Simulate marker inheritance from the father: ***/
     /*** Father has only one X-linked locus:          ***/
     father = subject->parents->husband;
     subject->chrom1 = father->chrom2;

     /*** Simulate marker inheritance from the mother: ***/
     mother = subject->parents->wife;
     u = runifAS183();
     if(u < 0.5)
       subject->chrom2 = mother->chrom1;
     else /*** (u >= 0.5) ***/
       subject->chrom2 = mother->chrom2;
   } /*** end else(subject->sex == FEMALE) ***/

  return;
}




/* TRAVERSE

                                   FUNCTION ARGUMENTS                         

traverse() takes the following arguments:

   -pedData:  A struct pedigreeData pointer to the pedigree structure being 
              processed.

   -func:     A function pointer to a function that returns type void.  The 
              function that is assigned to this argument should accept a struct 
              person pointer and a struct pedigreeData pointer, respectively, 
              such that the struct pedigreeData pointer points to a pedigree 
              that contains the person being pointed to by the struct person 
              pointer.

                                   SIDE EFFECTS

The traverse() function will update the values of the static global variables 
"REACHED", "COMPUTED", and "FLAG".  This is done to allow the a pedigree to be 
traversed more than one time.  This traversal is actually performed by 
traverse_engine().

Note:  This pedigree traversal algorithm is designed to work with cyclic as 
       well as acyclic pedigrees.

                                   RETURN VALUE

traverse() returns no value.

*/

static void traverse(struct pedigreeData *pedData,
                     void (*func)())
{
 if(FLAG == 0)
   {
    REACHED++;  COMPUTED++;
    traverse_engine(pedData->ped->proband,pedData,func);
    REACHED--;  COMPUTED--;
    FLAG = 1;
      }
 else
   {
    REACHED--;  COMPUTED--;
    traverse_engine(pedData->ped->proband,pedData,func);
    REACHED++;  COMPUTED++;
    FLAG = 0;
   }
 return;
}

/* TRAVERSE_ENGINE

                                   FUNCTION ARGUMENTS                         

traverse_engine() takes the following arguments:

   -proband:  A struct person pointer to the proband in a pedigree structure.

   -pedData:  A struct pedigreeData pointer to the pedigree structure being 
              processed.

   -func:     A function pointer to a function that returns type void.  The 
              function that is assigned to this argument should accept a struct 
              person pointer and a struct pedigreeData pointer, respectively, 
              such that the struct pedigreeData pointer points to a pedigree 
              that contains the person being pointed to by the struct person 
              pointer.

	                           SIDE EFFECTS

traverse_engine() updates the value of the struct person data member 
person.traverseStatus to the current global value of "REACHED" when a person
is reached, and the value of person.computeStatus will be set to the current
global calue of "COMPUTED" when the person has marker data computed.  Any 
side effects of functions assigned the the "func" function pointer argument 
above will also be side effects of the traverse_engine() function.

Note:  This pedigree traversal algorithm is designed to work with cyclic as 
       well as acyclic pedigrees.

                                   RETURN VALUE

traverse_engine() returns no value.

                                   SYNOPSIS

traverse_engine() is a recursive function that traverses through a pedigree and 
performs arbitrary operations upon the members of the pedigree, with these 
actions being specified by the function assigned to the "func" parameter.  

*/

static void traverse_engine(struct person *subject, 
                            struct pedigreeData *pedData,
                            void(*func)())
{
 struct children *currentChild;
 struct marriage *currentSpouse;

 /*** Mark this person so that recursive traversal calls will  not ***/
 /*** re-branch to this pedigree node:                             ***/

 subject->traverseStatus = REACHED;

 /******** Parents: ********/

 /*** Traverse to a parents only if they exist and neither one of them has ***/
 /*** been reached by a previous traverse_engine() call.                   ***/

 if( subject->parents                         != NULL       && 
    (subject->parents->husband->computeStatus != COMPUTED ||
     subject->parents->wife->computeStatus    != COMPUTED   )  )
   { 
    if(subject->parents->husband->computeStatus != COMPUTED)
      {
       traverse_engine(subject->parents->husband, pedData, func);
       return;
      }
    else
      {
       traverse_engine(subject->parents->wife, pedData, func);
       return;
      }
   }

 /*** segregate this subject's markers:  ***/ 

 if(subject->computeStatus != COMPUTED)
   {
    (*func)(subject, pedData);
    subject->computeStatus = COMPUTED;
   }

 /*** Spouse: ***/

 if(subject->sex == MALE)
   { currentSpouse = subject->nuclearFamily;
     while(currentSpouse != NULL)
       { 
	if(currentSpouse->wife->computeStatus != COMPUTED)
	  {
	   traverse_engine(currentSpouse->wife, pedData, func);
	  }
        currentSpouse = currentSpouse->nextHusbandMarriage;
       }
   } 

 if(subject->sex == FEMALE)
   { currentSpouse = subject->nuclearFamily;
     while(currentSpouse != NULL)
       { 
	if(currentSpouse->husband->computeStatus != COMPUTED)
          {
	   traverse_engine(currentSpouse->husband, pedData, func);
	  }
        currentSpouse = currentSpouse->nextWifeMarriage;
       }
   }

 /*** Children: ***/

 if(subject->sex == MALE)
   { currentSpouse = subject->nuclearFamily;
     while(currentSpouse != NULL)
       { 
         if(currentSpouse->wife->computeStatus == COMPUTED)
           {
	    currentChild = currentSpouse->child;
            while(currentChild != NULL)
              { 
               if(currentChild->firstChild->computeStatus != COMPUTED)
                  traverse_engine(currentChild->firstChild, pedData, func);
               currentChild = currentChild->nextChild;
              }
	   }
         currentSpouse = currentSpouse->nextHusbandMarriage;
       }
   }

 if(subject->sex == FEMALE)
   { currentSpouse = subject->nuclearFamily;
     while(currentSpouse != NULL)
       { 
         if(currentSpouse->husband->computeStatus == COMPUTED)
           {
	    currentChild = currentSpouse->child;
            while(currentChild != NULL)
              { 
               if(currentChild->firstChild->computeStatus != COMPUTED)
                  traverse_engine(currentChild->firstChild, pedData, func);
               currentChild = currentChild->nextChild;
              }
	   }
         currentSpouse = currentSpouse->nextWifeMarriage;
       }
   }

return;
}




/****************************************************************************/
/***               quicksort algorithm for large arrays                   ***/

static void qs(double array_sort[], long left, long right)
{
  long i,j;
  double mid, temp_sort;

  i = left; j = right;
  mid = array_sort[(left + right) / 2];

  do {
    while(array_sort[i] < mid && i < right) i++;
    while(mid < array_sort[j] && left < j) j--;
    
    if(i<=j) {
      temp_sort = array_sort[i];		/* swap elements	       */
      array_sort[i] = array_sort[j];
      array_sort[j] = temp_sort;
      
      i++; j--;
    }
  } while(i <= j);

  if(left < j) 
    qs(array_sort, left, j);
  if(i < right) 
    qs(array_sort, i, right);

  return; 
}

/*
     Algorithm AS 183 Appl. Statist. (1982) vol.31, no.2

     Returns a pseudo-random number uniformly (rectangularly) distributed
     between 0 and 1.   The cycle length is 6.95E+12 (See page 123
     of Applied Statistics (1984) vol.33), not as claimed in the
     original article.

     ix, iy and iz should be set to integer values between 1 and
     30000 before the first entry. To do this, 
     first call ranAS183_seed(iseed1,iseed2,iseed3), where iseed#
     are 3 long int seeds between 1 and 30000. The 3  seeds are
     saved, but ix,iy,iz can change.

*/

static long ix, iy, iz;

static long runifAS183_seed(long iseed1, long iseed2, long iseed3)
{
  long error;

  error=1;
  if( (iseed1 >=1 && iseed1 <=30000) && (iseed2 >=1 && iseed2 <=30000) && 
      (iseed3 >=1 && iseed3 <=30000)) error=0;
  if(error) return (error);
  ix = iseed1;
  iy = iseed2;
  iz = iseed3;
  return (error);
}

static double runifAS183()
{
   double u;

   ix = (171*ix) % 30269;
   iy = (172*iy) % 30307;
   iz = (170*iz) % 30323;
   u  = (double)ix/30269.0 + (double)iy/30307.0 + (double)iz/30323.0;
   return ( u - (int) u );
}


/****************************************************************************/

/* Poisson random variates (Ripley, Alg 3.3 ,page 55 )
   Prior to calling rpoisson, the function runifAS183_seed must be 
   called with three input starting seed values.

   Requires functions:  runifAS183 and runifAS183_seed

*/


static long rpoisson(double mu){
   
  long n;
  double c, p;

  c=exp(-mu); 
  p=1.0;
  n=0;
  do {
    p *=  runifAS183();
    n++;
  } while (p >= c);
  return n-1;
}

