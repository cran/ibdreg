#include <stdlib.h>

#define Salloc(l,s) mc_calloc((l),sizeof(s))
#define S_alloc(l,s) mc_calloc((l),(s))
#define Calloc(l,s) mc_calloc((l),(s))
#define Free(p) mc_free(p)

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
/* The function mc_malloc allocates dynmamic memory consisting of num*size  */
/* bytes of memory.  Be careful to check for NULL pointers                  */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void *mc_malloc(size_t size)
{ 
  void *memoryAllocation;
  struct mcMemory *temp;
  memoryAllocation = malloc(size);
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
/* The function mc_free frees a single dynamically allocated unit of        */
/* memory.  If the memory pointed to does not exist, no action is taken.    */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void mc_free(void *memory2free)
{
  struct mcMemory *previousMemory = NULL;
  struct mcMemory *currentMemory, *oldMemoryEnd;
  long error = 0;

  if(memoryStart == NULL)
    return;

  if(memoryStart == memoryEnd)
    {
     if(memory2free == memoryStart->firstMemory)
       {
	 free(memoryStart->firstMemory);
         free(memoryStart);
         memoryStart = NULL;
         memoryEnd   = NULL;
         return;
       }
     else
       {
	 return;
       }
    }

  /* Use psuedo-data: */
  memoryEnd->nextMemory = (struct mcMemory *)malloc(sizeof(struct mcMemory));

  /* Cannot use psuedo-data */
  if(memoryEnd->nextMemory == NULL)
    {
     currentMemory = memoryStart;
     while( (currentMemory->firstMemory != memory2free) && 
            (currentMemory->nextMemory  != NULL)          )
       {
        previousMemory = currentMemory;
        currentMemory = currentMemory->nextMemory;
       }
     if(currentMemory->firstMemory != memory2free)
       {
        /* Memory Allocation not found - we're at the end of the list */
        return;
       }
     else
       {
        /* Perform the requested memory deallocation: */
        free(currentMemory->firstMemory);
        if(previousMemory != NULL)
          previousMemory->nextMemory = currentMemory->nextMemory;
        else /* (currentMemory == memoryStart) */
          memoryStart = memoryStart->nextMemory;
        if(currentMemory == memoryEnd)
          memoryEnd = previousMemory;
        free(currentMemory);
        return;
       } 
    }

  /* Can use psuedo-data */
  oldMemoryEnd = memoryEnd;
  memoryEnd = memoryEnd->nextMemory;
  memoryEnd->firstMemory = memory2free;
  memoryEnd->nextMemory  = NULL;

  currentMemory = memoryStart;
  while(currentMemory->firstMemory != memory2free)
    {
     previousMemory = currentMemory;
     currentMemory = currentMemory->nextMemory;
    }
  if(currentMemory == memoryEnd)
    {
      /* Memory Allocation not found */
      free(currentMemory);
      memoryEnd = oldMemoryEnd;
      memoryEnd->nextMemory = NULL;
      return;
    }
  else
    {
      /* Found the memory allocation - free psuedo-data first: */
      free(memoryEnd);
      memoryEnd = oldMemoryEnd;
      memoryEnd->nextMemory = NULL;
      /* Now perform the requested memory deallocation: */
      free(currentMemory->firstMemory);
      previousMemory->nextMemory = currentMemory->nextMemory;
      if(currentMemory == memoryStart)
        memoryStart = memoryStart->nextMemory;
      if(currentMemory == memoryEnd)
        memoryEnd = previousMemory;
     free(currentMemory);
     return;
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

static void mc_free_all()
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

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/*                   Local Memory Management functions                      */
/*                                                                          */
/* The following functions allow memory management to be handled at a local */
/* level, thus allowing functions to clean up memory that is no longer      */
/* needed without freeing memory containing data that is still required.    */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void *mc_calloc_local(size_t num, size_t size, 
                             struct mcMemory **memoryStart, 
                             struct mcMemory **memoryEnd)
{ 
  void *memoryAllocation;
  struct mcMemory *temp;
  memoryAllocation = calloc(num, size);
  if(memoryAllocation == NULL)
    return memoryAllocation;
  if(*memoryStart == NULL)
    {
     /* This is the first memory allocation for this program. */
     *memoryStart = (struct mcMemory *)malloc(sizeof(struct mcMemory));
     if(*memoryStart == NULL)
       { 
        free(memoryAllocation);
        return NULL;
       }
     (*memoryStart)->firstMemory = memoryAllocation;
     (*memoryStart)->nextMemory = NULL;
     *memoryEnd = *memoryStart;
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
     (*memoryEnd)->nextMemory = temp;
     *memoryEnd = temp;
     return memoryAllocation;
    }
}

static void *mc_malloc_local(size_t size, 
                             struct mcMemory **memoryStart, 
                             struct mcMemory **memoryEnd)
{ 
  void *memoryAllocation;
  struct mcMemory *temp;
  memoryAllocation = malloc(size);
  if(memoryAllocation == NULL)
    return memoryAllocation;
  if(*memoryStart == NULL)
    {
     /* This is the first memory allocation for this program. */
     *memoryStart = (struct mcMemory *)malloc(sizeof(struct mcMemory));
     if(*memoryStart == NULL)
       { 
        free(memoryAllocation);
        return NULL;
       }
     (*memoryStart)->firstMemory = memoryAllocation;
     (*memoryStart)->nextMemory = NULL;
     *memoryEnd = *memoryStart;
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
     (*memoryEnd)->nextMemory = temp;
     *memoryEnd = temp;
     return memoryAllocation;
    }
}

static void mc_free_local(void *memory2free,
                          struct mcMemory **memoryStart, 
                          struct mcMemory **memoryEnd)
{
  struct mcMemory *previousMemory = NULL;
  struct mcMemory *currentMemory, *oldMemoryEnd;
  long error = 0;

  if(*memoryStart == NULL)
    return;

  if(*memoryStart == *memoryEnd)
    {
     if(memory2free == (*memoryStart)->firstMemory)
       {
	 free((*memoryStart)->firstMemory);
         free(*memoryStart);
         *memoryStart = NULL;
         *memoryEnd   = NULL;
         return;
       }
     else
       {
	 return;
       }
    }

  /* Use psuedo-data: */
  (*memoryEnd)->nextMemory = (struct mcMemory *)malloc(sizeof(struct mcMemory));

  /* Cannot use psuedo-data */
  if((*memoryEnd)->nextMemory == NULL)
    {
     currentMemory = *memoryStart;
     while( (currentMemory->firstMemory != memory2free) && 
            (currentMemory->nextMemory  != NULL)          )
       {
        previousMemory = currentMemory;
        currentMemory = currentMemory->nextMemory;
       }
     if(currentMemory->firstMemory != memory2free)
       {
        /* Memory Allocation not found - we're at the end of the list */
        return;
       }
     else
       {
        /* Perform the requested memory deallocation: */
        free(currentMemory->firstMemory);
        if(previousMemory != NULL)
          previousMemory->nextMemory = currentMemory->nextMemory;
        else /* (currentMemory == memoryStart) */
          *memoryStart = (*memoryStart)->nextMemory;
        if(currentMemory == *memoryEnd)
          *memoryEnd = previousMemory;
        free(currentMemory);
        return;
       } 
    }

  /* Can use psuedo-data */
  oldMemoryEnd = *memoryEnd;
  *memoryEnd = (*memoryEnd)->nextMemory;
  (*memoryEnd)->firstMemory = memory2free;
  (*memoryEnd)->nextMemory  = NULL;

  currentMemory = *memoryStart;
  while(currentMemory->firstMemory != memory2free)
    {
     previousMemory = currentMemory;
     currentMemory = currentMemory->nextMemory;
    }
  if(currentMemory == *memoryEnd)
    {
      /* Memory Allocation not found */
      free(currentMemory);
      *memoryEnd = oldMemoryEnd;
      (*memoryEnd)->nextMemory = NULL;
      return;
    }
  else
    {
      /* Found the memory allocation - free psuedo-data first: */
      free(*memoryEnd);
      *memoryEnd = oldMemoryEnd;
      (*memoryEnd)->nextMemory = NULL;
      /* Now perform the requested memory deallocation: */
      free(currentMemory->firstMemory);
      previousMemory->nextMemory = currentMemory->nextMemory;
      if(currentMemory == *memoryStart)
        *memoryStart = (*memoryStart)->nextMemory;
      if(currentMemory == *memoryEnd)
        *memoryEnd = previousMemory;
     free(currentMemory);
     return;
    }
}

static void mc_free_all_local(struct mcMemory **memoryStart, 
                              struct mcMemory **memoryEnd)
{ 
  struct mcMemory *currentMemory;
  struct mcMemory *temp;
  if(*memoryStart == NULL)
    return;

  /* Free all the allocated memory and the memory structs, which were also */
  /* dynamically allocated:                                                */
  
  currentMemory = *memoryStart;
  while(currentMemory != NULL)
    {
      free(currentMemory->firstMemory);
      temp = currentMemory;
      currentMemory = currentMemory->nextMemory;
      free(temp);
    }
  *memoryStart = NULL;
  *memoryEnd = NULL;
  return;
}

