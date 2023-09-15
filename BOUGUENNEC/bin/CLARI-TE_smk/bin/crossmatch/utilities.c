/*****************************************************************************
#   Copyright (C) 1993-2000 by Phil Green.                          
#   All rights reserved.                           
#                                                                           
#   This software is part of a beta-test version of the swat/cross_match/phrap 
#   package.  It should not be redistributed or
#   used for any commercial purpose, including commercially funded
#   sequencing, without written permission from the author and the
#   University of Washington.
#   
#   This software is provided ``AS IS'' and any express or implied
#   warranties, including, but not limited to, the implied warranties of
#   merchantability and fitness for a particular purpose, are disclaimed.
#   In no event shall the authors or the University of Washington be
#   liable for any direct, indirect, incidental, special, exemplary, or
#   consequential damages (including, but not limited to, procurement of
#   substitute goods or services; loss of use, data, or profits; or
#   business interruption) however caused and on any theory of liability,
#   whether in contract, strict liability, or tort (including negligence
#   or otherwise) arising in any way out of the use of this software, even
#   if advised of the possibility of such damage.
#                                                                         
#*****************************************************************************/

#if vms
#include stdio
#else
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#endif

char *plural(int_var)
     int int_var;
{
  return int_var == 1 ? "" : "s";
}

/* N.B. following needs to be tested for non-unix systems */
fatalError(message)
     char *message;
{
  fprintf(stdout, "\nFATAL ERROR: %s\n", message);
  fprintf(stderr, "\nFATAL ERROR: %s\n", message);
  exit(1);
}

notify(message)
     char *message;
{
  fprintf(stderr, message);
  if (!strcmp(message, " Done\n")) check_alloc_list();
  fflush(stderr);
  fflush(stdout);
}

FILE *fopenRead(fileName)
  char *fileName;
{
  FILE *fp;
  FILE *fopen();

  if(!(fp = fopen(fileName,"r"))) {
    fprintf(stderr, "\n ERROR: FILE %s NOT FOUND \n",fileName);
    exit(1);
  }
  else return fp;
}

FILE *fopenWrite(fileName)
  char *fileName;
{
  FILE *fp;
  FILE *fopen();

  if(!(fp = fopen(fileName,"w"))) {
    fprintf(stderr, "\n ERROR: FILE %s NOT AVAILABLE FOR WRITING \n",fileName);
    exit(1);
  }
  else return fp;
}

FILE *fopenAppend(fileName)
  char *fileName;
{
  FILE *fp;
  FILE *fopen();

  if(!(fp = fopen(fileName,"a"))) {
    fprintf(stderr, "\n ERROR: FILE %s NOT AVAILABLE FOR WRITING \n",fileName);
    exit(1);
  }
  else return fp;
}

/* N.B. The following is modified from Kernighan & Ritchie, The C Programming
Language, 1st edition */

typedef size_t ALLOC; /* typedef unsigned int ALLOC; */
typedef short SHORT;

typedef double ALIGN;

union header {
	struct {
		union header *ptr;
		ALLOC size;
	} s;
	ALIGN  x;
};

typedef union header HEADER;

static HEADER base;
static HEADER *allocp = 0;	
static unsigned long int total_allocated;

char *our_alloc(nbytes)
        ALLOC nbytes;
{
	HEADER *morecore();
	HEADER *p;
	HEADER *q;
	ALLOC nunits;
 
	nunits = 1 +(nbytes + sizeof(HEADER) - 1)/ sizeof(HEADER);
	if (!(q = allocp)) {
		base.s.ptr = allocp = q = &base;
		base.s.size = 0;
	}
	for (p = q->s.ptr;;q = p, p = p->s.ptr) {
 		if (p->s.size >= nunits){
			if(p->s.size == nunits) q->s.ptr = p->s.ptr;
			else {
				p->s.size -= nunits;
				p += p->s.size;
				p->s.size = nunits;
			}
			allocp = q;
			return ((char *)(p + 1));
		}
		if (p == allocp)
			if (!(p = morecore(nunits))) return 0;
	}
}

check_alloc_list()
{
	long int num_blocks, num_ends, num_units;
	HEADER *p, *q;

	p = allocp;
	num_blocks = num_ends = num_units = 0;
	do {
	  q = p->s.ptr;
/* Now that ALLOC is unsigned, following is vacuous
	  if (p->s.size < 0)
	    fatalError(" our_alloc linked list corruption (check_alloc_list1)");
*/
	  num_units += p->s.size;
	  if (q > p && p + p->s.size > q)
	    fatalError(" our_alloc linked list corruption (check_alloc_list2)");
	  if (q < p) num_ends++;
	  p = q;
	  num_blocks++;
	} while (p != allocp);
	if (num_ends > 1 || num_ends == 0 && num_blocks != 1) {
	  fprintf(stderr, "\nnum_ends: %ld, num_blocks: %ld\n", num_ends, num_blocks);
	  fatalError(" our_alloc linked list corruption (check_alloc_list3)");
	}
	/*
	if (num_ends != 1)
	    fatalError(" our_alloc linked list corruption (check_alloc_list3)");
	*/
	fprintf(stderr, "Total space allocated: %.3f Mbytes; currently free: %.3f Mbytes in %ld blocks\n",
		total_allocated / 1000000.0, num_units * sizeof(HEADER) / 1000000.0, num_blocks);
}

find_frags()
{
	SHORT num_frags;
	HEADER *p;

	for (p = allocp->s.ptr, num_frags = 1; p != allocp;
            p = p->s.ptr, num_frags++);
        fprintf(stderr, "\n number of allocation fragments = %d", num_frags);
}

#define  NALLOC  63000

extern void *malloc();
HEADER *morecore(nu)
        ALLOC nu;
{
	char *cp;
	HEADER *up;
	ALLOC rnu;
        SHORT our_free();
   
	rnu = NALLOC * ((nu + NALLOC - 1) / NALLOC);
	cp = malloc(rnu * sizeof(HEADER));
	if (!cp) {
	  fprintf(stderr, "\n\n%.3f Mbytes requested but unavailable", 
		  rnu * sizeof(HEADER) / 1000000.0);
	  check_alloc_list();
	  fatalError("REQUESTED MEMORY UNAVAILABLE");
	}
	total_allocated += (unsigned long)(rnu * sizeof(HEADER));
        fprintf(stderr, "%.3f Mbytes allocated -- total %.3f Mbytes\n", 
		rnu * sizeof(HEADER) / 1000000.0, total_allocated / 1000000.0);
	up = (HEADER *)cp;
	up->s.size = rnu;
	our_free((char *)(up + 1));
	return allocp;
}

SHORT our_free(ap)
     char *ap;
{
	HEADER *p;
	HEADER *q;
	
	if (!ap) return 0; 

	p = (HEADER *)ap - 1;
	for (q = allocp; !(p > q && p < q->s.ptr); q = q->s.ptr)
		if (q>= q->s.ptr && (p > q || p < q->s.ptr)) break;

	if ( p + p->s.size == q->s.ptr){
		p->s.size += q->s.ptr->s.size;
		p->s.ptr = q->s.ptr->s.ptr;
	} 
	else if (p + p->s.size > q->s.ptr && p <= q->s.ptr) {
	  fatalError(" our_alloc linked list corruption (our_free1)");
	}
	else p->s.ptr = q->s.ptr;
	if (q + q->s.size == p){
		q->s.size += p->s.size;
		q->s.ptr = p->s.ptr;
	} 
	else if (q + q->s.size > p && p >= q) {
	  fatalError(" our_alloc linked list corruption (our_free2)");
	}
	else q->s.ptr = p;
	allocp = q;
        return 0;
}
 
