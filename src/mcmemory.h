#include <stdlib.h>

#define Salloc(l,s) mc_calloc((l),sizeof(s))
#define Calloc(l,s) mc_calloc((l),(s))

/*
 #define Free(p) mc_free(p)  
rename mc_free to just Free so this line not needed, so this not needed  JPS 9/2021
*/

struct mcMemory
{ void *firstMemory;
  struct mcMemory  *nextMemory;
};

static struct mcMemory* memoryStart = NULL;
static struct mcMemory* memoryEnd = NULL;

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* The function mc_calloc allocates dynmamic memory consisting of num*size  */
/* bytes of memory initialized to zero.  Be careful to check for NULL       */
/* pointers.                                                                */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void *mc_calloc(size_t num, size_t size)
{ 
  void *memoryAllocation;
  struct mcMemory *temp;
  memoryAllocation = calloc(num, size);
  if(memoryAllocation == NULL)
    return memoryAllocation;
  if(memoryStart == NULL)
    {
     /* This is the first memory allocation for this program. */
     memoryStart = (struct mcMemory *)malloc(sizeof(struct mcMemory));
     if(memoryStart == NULL)
       { 
        free(memoryAllocation);
        return NULL;
       }
     memoryStart->firstMemory = memoryAllocation;
     memoryStart->nextMemory = NULL;
     memoryEnd = memoryStart;
     return memoryAllocation;
    }
  else
    {
     /* Memory has been previously allocated. */
     temp = (struct mcMemory *)malloc(sizeof(struct mcMemory));
     if(temp == NULL)
       { 
        free(memoryAllocation);
        return NULL;
       }
     temp->firstMemory = memoryAllocation;
     temp->nextMemory = NULL;
     memoryEnd->nextMemory = temp;
     memoryEnd = temp;
     return memoryAllocation;
    }
}

 


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/ 
/*                                                                          */
/* The function mc_free_all frees all dynamically allocated memory - use    */
/* immediately before program termination.  If no memory has been           */
/* allocated, no action will be taken.                                      */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void mc_free_all(void) // JPS added void 11/2022
{ 
  struct mcMemory *currentMemory;
  struct mcMemory *temp;
  if(memoryStart == NULL)
    return;

  /* Free all the allocated memory and the memory structs, which were also */
  /* dynamically allocated:                                                */
  
  currentMemory = memoryStart;
  while(currentMemory != NULL)
    {
      free(currentMemory->firstMemory);
      temp = currentMemory;
      currentMemory = currentMemory->nextMemory;
      free(temp);
    }
  memoryStart = NULL;
  memoryEnd = NULL;
  return;
}
