/*****************************************************************************
#   Copyright (C) 1994-1999 by Phil Green.                          
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

#include "swat.h"

/* make free_segment empty, and restore original form of alloc_segment, for
   faster (but memory-hungry) performance */

/* insert segment into a linked list, and revise it to ensure non-overlap cond''n */

Segment *insert_segment(head, start, end)
     Segment *head;
     int start, end;
{
  Segment *prev_segment, *segment, *segment2;
  Segment *alloc_segment();
  
  if (!head) { /* construct head if it doesn''t already exist */
    head = alloc_segment(); /* (Segment *)our_alloc(sizeof(Segment)); */
    head->start = start;
    head->end = end;
    head->next = 0;
    return head;
  }
  for (prev_segment = 0, segment = head; segment && start - 1 > segment->end;
       prev_segment = segment, segment = segment->next);
  if (!segment || end + 1 < segment->start) {
    segment2 = alloc_segment(); /* (Segment *)our_alloc(sizeof(Segment)); */
    segment2->start = start;
    segment2->end = end;
    segment2->next = segment;
    if (prev_segment) prev_segment->next = segment2;
    else return segment2;
  }
  else {
    if (start < segment->start) segment->start = start;
    if (end > segment->end) {
      segment->end = end;
      for (segment2 = segment->next; segment2 && end + 1 >= segment2->start;
	   segment2 = segment2->next) {
	if (end < segment2->end) segment->end = segment2->end;
	segment->next = segment2->next;
/*	free_segment(segment2);  */
      }
    }
  }
  return head;
}

/* as above, but allow overlaps -- remove segment only if completely
   contained in another one. In overlapping list, starts are strictly increasing;
   ends are strictly increasing; but start can be less than end of preceding
   segment */
Segment *weak_insert_segment(head, start, end)
     Segment *head;
     int start, end;
{
  Segment *prev_segment, *segment, *segment2;
  Segment *alloc_segment();

  if (end < start) fatalError("weak inserted segment");
  if (!head) { /* construct head if it doesn''t already exist */
    head = alloc_segment(); /* (Segment *)our_alloc(sizeof(Segment)); */
    head->start = start;
    head->end = end;
    head->next = 0;
    return head;
  }
  for (prev_segment = 0, segment = head; 
       segment && start > segment->start;
       prev_segment = segment, segment = segment->next);
  if (prev_segment && prev_segment->end >= end
      || segment && segment->start == start && segment->end >= end) 
    return head; /* new segment is contained in an old one */
  if (!segment || segment->end > end) { /* does not contain an old segment */
    segment2 = alloc_segment(); /* (Segment *)our_alloc(sizeof(Segment)); */
    segment2->start = start;
    segment2->end = end;
    segment2->next = segment;
    if (prev_segment) prev_segment->next = segment2;
    else return segment2;
  }
  else { /* contains an old segment */
    segment->start = start;
    segment->end = end;
    for (segment2 = segment->next; segment2 && end >= segment2->end;
	 segment2 = segment2->next) {
      segment->next = segment2->next;
    }
  }
  return head;
}

print_segments(head, index, complement) 
     Segment *head;
     int index;
     int complement;
{
  Segment *segment;

  if (!head) return;
  printf("Contig %2d   %c ", index, complement ? 'C' : ' ');
  for (segment = head; segment; segment = segment->next) {
    printf("%6d %6d   ", 
	   complement ? -segment->start : segment->start, 
	   complement ? -segment->end : segment->end);
  } 
}

old_print_segments(head, index1, index2, offset, rel_orient) 
     Segment *head;
     int index1, index2, offset;
     char rel_orient;
{
  Segment *segment;

  if (!head) return;
  printf("\nContig %2d     ", index1);
  for (segment = head; segment; segment = segment->next) 
    printf("%6d %6d   ", segment->start, segment->end);
  printf("matches \nContig %2d   %c ", index2, rel_orient ? 'C' : ' ');
  for (segment = head; segment; segment = segment->next) {
      printf("%6d %6d   ", segment->start - offset, segment->end - offset);
  } /* N.B. segments won't be in correct order! */
  printf("\n\n");
}

/* finds maximum gap size for a linked list of segments */	
int find_max_gap(head, start, end)
     Segment *head;
     int start, end;
{
  int max_gap, last_end;
  Segment *segment;

  max_gap = 0;
  last_end = start;
  for (segment = head; segment && segment->start <= end; segment = segment->next) {
    if (segment->start - last_end > max_gap) max_gap = segment->start - last_end;
    if (last_end < segment->end + 1) last_end = segment->end + 1;
  }
  if (end + 1 - last_end > max_gap) max_gap = end + 1 - last_end;
  return max_gap;
}

/* finds maximum gap size (in 2d, with respect to first) for two linked list of segments */	
int find_max_gap_list(head1, head2)
     Segment *head1, *head2;
{
  int max_gap, gap;
  Segment *segment;

  gap = max_gap = 0;
  for (segment = head2; segment; segment = segment->next) {
    gap = find_max_gap(head1, segment->start, segment->end);
    if (gap > max_gap) max_gap = gap;
  }
  return max_gap;
}

/* finds parts of segments in 2d list that don't occur in first, for two linked list of segments */	
Segment *find_gap_segment_list(head1, head2)
     Segment *head1, *head2;
{
  int last_end;
  Segment *segment, *segment2, *gap_segs;
  Segment *insert_segment();

  gap_segs = 0;
  for (segment2 = head2; segment2; segment2 = segment2->next) {
    last_end = segment2->start;
    for (segment = head1; segment && segment->start <= segment2->end; segment = segment->next) {
      if (segment->start > last_end) 
	gap_segs = insert_segment(gap_segs, last_end, segment->start - 1);
      if (last_end < segment->end + 1) last_end = segment->end + 1;
    }
    if (segment2->end >= last_end)
      gap_segs = insert_segment(gap_segs, last_end, segment2->end);
  }
  return gap_segs;
}

check_segments(head)
     Segment *head;
{
  Segment *segment;
  int prev_end;

  for (segment = head; segment; segment = segment->next) {
    if (segment->start > segment->end)
      fatalError("Segment corruption1");
    if (segment != head && segment->start < prev_end)
      fatalError("Segment corruption2");
    prev_end = segment->end;
  }
}

print_segment_gaps(fp, head, start, end)
     FILE *fp;
     Segment *head;
     int start, end;
{
  int last_end, flag, s_flag;
  Segment *segment;

  last_end = start;
  flag = 0;
  s_flag = 1;
  for (segment = head; segment; segment = segment->next) {
    if (segment->start > last_end + 1) {
      fprintf(fp, "%c %c %d- %d", flag ? ',' : ' ', s_flag ? 'S' : 'I', last_end, segment->start - 1);
      flag = 1;
    }
    last_end = segment->end + 1;
    s_flag = 0;
  }
  if (end > last_end) 
      fprintf(fp, "%c E %d- %d", flag ? ',' : ' ', last_end, end);
  else if (!flag) fprintf(fp, " None.");
}

int percent_contained(head, start, end) /* assumes start <= end */
     Segment *head;
     int start, end;
{
  Segment *segment;
  int start0, end0;
  int covered_bases;

  covered_bases = 0;
  for (segment = head; segment && segment->start <= end; segment = segment->next) {
    if (segment->end < start) continue;
    start0 = start < segment->start ? segment->start : start;
    end0 = end > segment->end ? segment->end : end;
    if (end0 > start0) covered_bases += end0 - start0 + 1;
  }
  return (int)((100.0 * covered_bases) / (end - start + 1));
}

/* following not possible with faster segment allocation above
free_segments(head)
     Segment *head;
{
  Segment *segment;

  for (segment = head; segment; segment = segment->next)
    free_segment(segment);
}

free_segment(segment)
  Segment *segment;
{
  our_free(segment);
}
*/


#define SEG_BLOCK_SIZE 5000

typedef struct seg_block {
  Segment *head_segment;
  struct seg_block *next;
} Seg_block;

/* head_block is current (most recently created) block; save_block (and all following blocks) should
   be saved */
static Seg_block *head_block, *save_block;
static Segment *head_segment, *save_head_segment;

static int n_segments, save_n_segments;

Segment *alloc_segment()
{
  Segment *segment;
  char *our_alloc();
  Seg_block *temp_block;

/*   return (Segment *)our_alloc(sizeof(Segment)); delete to restore original (faster) version */

  if (!n_segments) {
    head_segment = (Segment *)our_alloc(SEG_BLOCK_SIZE * sizeof(Segment));
    temp_block = (Seg_block *)our_alloc(sizeof(Seg_block));
    temp_block->next = head_block;
    temp_block->head_segment = head_segment;
    head_block = temp_block;
  }
  segment = head_segment + n_segments;
  n_segments = n_segments < SEG_BLOCK_SIZE - 1 ? n_segments + 1 : 0;
  return segment;
}

free_seg_blocks()
{
  Seg_block *next;

  for (; head_block != save_block; head_block = next) {
    our_free(head_block->head_segment);
    next = head_block->next;
    our_free(head_block);
  }
  n_segments = save_n_segments;
  head_segment = save_head_segment;  
}

mark_save_block()
{
  save_block = head_block;
  save_n_segments = n_segments;
  save_head_segment = head_segment;
} 

print_t_n_segments()
{
  int n_blocks;
  Seg_block *t_block;

  for (n_blocks = 0, t_block = head_block; t_block; t_block = t_block->next) n_blocks++;

  fprintf(stderr, "Total # segment blocks: %d, size: %.3f Mbytes\n",
	  n_blocks, n_blocks * SEG_BLOCK_SIZE *sizeof(Segment) / 1000000.);
}
 
