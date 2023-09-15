/*****************************************************************************
#   Copyright (C) 1995-2000 by Phil Green.                          
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

#define PAIR_BLOCK_SIZE 10000

typedef struct cand_pair_block {
  Cand_pair *head_pair;
  struct cand_pair_block *next;
} Cand_pair_block;

extern Parameters *parameters; /* from get_parameters() */

static Cand_pair_block *head_pair_block;
static int n_cand_pairs;

static int num_cand_pairs, num_dups, num_repeat_rejected_pairs;

cluster_pairs(entry1, entry2, start1, start2, reverse, rev_store_flag)
     int entry1, entry2, start1, start2;
     int reverse, rev_store_flag;
{
  int get_equiv_class();
  int entry, class, class1, class2, final_class;

  class1  = get_equiv_class(entry1);
  class2  = get_equiv_class(entry2);
  if (class1 == class2) return;

  for (class = class1, entry = entry1; 
       class != entry; 
       entry = class, class = get_equiv_class(entry));

  final_class = class;

  for (class = class2, entry = entry2; 
       class != entry; 
       entry = class, class = get_equiv_class(entry));

  if (class < final_class) final_class = class;

  for (class = class1, entry = entry1; 
       class != final_class; 
       entry = class, class = get_equiv_class(entry))
    set_equiv_class(entry, final_class);

  for (class = class2, entry = entry2; 
       class != final_class; 
       entry = class, class = get_equiv_class(entry))
    set_equiv_class(entry, final_class);

  set_equiv_class(entry, final_class);
}

make_new_cand_pair(entry1, entry2, start1, start2, reverse, rev_store_flag)
     int entry1, entry2, start1, start2;
     int reverse, rev_store_flag;
     /* offset currently not used -- but should be! */
{
  Cand_pair *pair, *parent;
  Cand_pair *get_cand_pairs();
  int off, off_min, off_max, temp;
  Segment *insert_segment();
  Seq_entry *get_seq_entry();
  Tag *tag;
  int i, j, t_start, t_end, store_entry, test_entry;
/*
  if (start1 < -1) { implies call is from make_full_pairs(); 
    off_min = 1;
    off_max = 0;
  }
*/

  off = start1 - start2;
  if (entry1 == entry2 && !reverse) {
    if (off < 0) off = -off;
    off_min = off <= parameters->bandwidth ? 1 : off - parameters->bandwidth;
  }
  else {
    off_min = off - parameters->bandwidth;
  }
  off_max = off + parameters->bandwidth;

  if (rev_store_flag) {
    for (parent = pair = get_cand_pairs(entry2); pair; ) {
      parent = pair;
      if (pair->entry1 == entry1 && pair->reverse == reverse) {
	/* readjust band size, if necessary, to encompass the diagonal containing
	   start1 and start2 */
	pair->band_segments = insert_segment(pair->band_segments, off_min, off_max);
	num_dups++;
	return;
      }
      pair = entry1 < pair->entry1 ? pair->left : pair->right;
    }
  }
  else {
    for (parent = pair = get_cand_pairs(entry1); pair; ) {
      parent = pair;
      if (pair->entry2 == entry2 && pair->reverse == reverse) {
	/* readjust band size, if necessary, to encompass the diagonal containing
	   start1 and start2 */
	pair->band_segments = insert_segment(pair->band_segments, off_min, off_max);
	num_dups++;
	return;
      }
      pair = entry2 < pair->entry2 ? pair->left : pair->right;
    }
  }
  num_cand_pairs++;
  
  append_cand_pair(entry1, entry2, reverse, off_min, off_max, parent, rev_store_flag);
}

append_cand_pair(entry1, entry2, reverse, left, right, parent, rev_store_flag)
     int entry1, entry2, reverse, left, right, rev_store_flag;
     Cand_pair *parent;
{
  char *our_alloc();
  Cand_pair *pair;
  static Cand_pair *head_pair;
  Cand_pair_block *temp_block;
  Segment *insert_segment();
  Seq_entry *seq_entry;
  Seq_entry *get_seq_entry();

  /* allocate pairs in blocks of PAIR_BLOCK_SIZE */
  if (!n_cand_pairs) {
    head_pair = (Cand_pair *)our_alloc(PAIR_BLOCK_SIZE * sizeof(Cand_pair));
    temp_block = (Cand_pair_block *)our_alloc(sizeof(Cand_pair_block));
    temp_block->head_pair = head_pair;
    temp_block->next = head_pair_block;
    head_pair_block = temp_block;
  }
  pair = head_pair + n_cand_pairs;
  n_cand_pairs = (n_cand_pairs + 1) % PAIR_BLOCK_SIZE;

  pair->entry1 = entry1;
  pair->entry2 = entry2;
  pair->reverse = reverse;
  pair->band_segments = insert_segment(0, left, right);
  pair->left = pair->right = 0;

  if (rev_store_flag) {
    if (!parent) {
      seq_entry = get_seq_entry(entry2);
      seq_entry->cand_pairs = pair;
    }
    else {
      if (entry1 < parent->entry1) parent->left = pair;
      else parent->right = pair;
    }
  }
  else {
    if (!parent) {
      seq_entry = get_seq_entry(entry1);
      seq_entry->cand_pairs = pair;
    }
    else {
      if (entry2 < parent->entry2) parent->left = pair;
      else parent->right = pair;
    }
  }
}

free_cand_pair_blocks()
{    
  Cand_pair_block *next;

  for (; head_pair_block; head_pair_block = next) {
    our_free(head_pair_block->head_pair);
    next = head_pair_block->next;
    our_free(head_pair_block);
  }
  n_cand_pairs = 0;
}

print_num_cand_pairs()
{
  fprintf(stderr, "\nnum_cand_pairs: %d, num add'l hits to pre-existing pair: %d\n", num_cand_pairs, num_dups);
}

