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

extern Parameters *parameters;
extern int num_pairs, t_num_entries;

#define MAX_CHIMERA_GAP 30  /* maximum gap between distinct confirmed segments for 
			       putative chimeras */
compare_pair_entries(pair_1, pair_2)
     Aligned_pair **pair_1, **pair_2;
{
  Aligned_pair *pair1, *pair2;
  int d;

  pair1 = *pair_1;
  pair2 = *pair_2;
  if (d = pair1->entry1 - pair2->entry1) return d;
  if (d = pair1->entry2 - pair2->entry2) return d;  
  if (d = is_reverse(pair1) - is_reverse(pair2)) return d;
  if ((d = !!pair1->score - !!pair2->score) || !pair1->score) return d;
  return (pair1->start1 - pair2->start1);
}

find_deletions()
{
  int i_ptr, entry1, entry2, prev_del, d_size, n_deletions;
  Aligned_pair *pair, *pair2, *pair3;
  Aligned_pair *get_aligned_pairs();
  Aligned_pair **sort_pairs();
  Aligned_pair **pair_pointers;
  Align_info *align_entry1;
  Align_info *get_align_entry();
  unsigned char *diff;
  char *get_id();
  int n_dels;
  int site1, t;

  notify("Finding deletions ... ");
  pair_pointers = sort_pairs(compare_pair_entries, 1);
  printf("\n\nProbable deletion reads (excluded from assembly):");
  n_deletions = 0;
  entry1 = entry2 = prev_del = -1;
  for (i_ptr = 0; i_ptr < 2 * num_pairs; i_ptr++) {
    pair = pair_pointers[i_ptr];
    if (!pair->score) continue;
    if (pair->entry1 == entry1 && pair->entry2 == entry2) { 
      pair2 = pair_pointers[i_ptr-1];
      if (is_reverse(pair) != is_reverse(pair2)) continue;
      if (entry1 == entry2) continue;
      d_size = pair->start2 - pair->start1 - (pair2->end2 - pair2->end1);
      if (d_size < 10 || d_size > 5000) continue; /* deletion smaller than 10 bases -- ignore */
      if (pair->start1 > pair2->end1 + 20) continue; /* segments are separated by more than
							20 bases; this is most likely
							due to a repeat of some sort, not a
							deletion */
      align_entry1 = get_align_entry(entry1);
/* use strong penalty alignments, rather than weak ones, to check for confirming
   read */
      for (pair3 = get_aligned_pairs(entry1); pair3; pair3 = pair3->next) {
	if (!pair3->score) continue;
	if (pair3->start1 < pair2->end1 - 5 && pair3->end1 > pair->start1 + 5) {
	  n_dels = 0;
	  site1 = pair3->start1 - 1;
	  for (diff = pair3->diffs; *(diff + 1); diff++) {
	    site1 += diff_gap1(*diff);
	    if (site1 > pair2->end1 - 5 
		&& site1 < pair->start1 + 5) {
	      if ((t = diff_type(*diff)) == 'I') n_dels++;
	      else if (t == 'D') n_dels--;
	    }
	  }
	  if (n_dels < 10) break;
	}
      }
      if (pair3) continue;
	    
/*
      for (segment = align_entry1->segments; segment; segment = segment->next) 
	if (segment->start < pair2->end1 - 5 && segment->end > pair->start1 + 5) 
	  break;
      if (segment) continue;
      if (align_entry1->segments->start < pair2->end1 - 5 && 
	  align_entry1->segments->end > pair->start1 + 5) continue;
*/
      /* there are spanning reads; should do a better test for this! */
      if (entry1 != prev_del) {
	printf("\n");
	prev_del = entry1;
	n_deletions++;
      }
      printf("\n%-12s  %2d   %3d-%3d  (%s %-12s  %3d-%3d)", 
	     get_id(entry1), d_size, 
	     pair2->end1 + 1, pair->start1 + 1, 
	     is_reverse(pair) ? "C" : " ", get_id(entry2), 
	     pair2->end2 + 1, pair->start2 + 1);
      set_deleted(align_entry1);
/* exists spanning read; note this may not be complete check! */
/* if (!align_entry1->segments || !align_entry1->segments->next) continue; 
   only consider candidate chimeras; note that chimera condition only uses pairs 
   with ->best == 1, so it is not as general as possible */
    }
    entry1 = pair->entry1;
    entry2 = pair->entry2;
  }
  if (n_deletions) printf("\n\n%d probable deletion reads.\n", n_deletions);
  else printf(" None.");
  our_free(pair_pointers);
  notify(" Done\n");
}

/*
reset_hq_mismatch()
{
  int entry1;
  Aligned_pair *pair;

  for (entry1 = 0; entry1 < t_num_entries; entry1++) 
    for (pair = align_array[entry1].first_pair; pair; pair = pair->next) 
      pair->hq_mismatch = pair->ss_mismatch = 0; 
}
*/

find_truncated_pairs()
{
  int entry1;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs(); 

  notify("Finding truncated pairs ... ");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (pair->entry1 >= pair->entry2) continue; 
      set_left_trunc_flag(pair, reject_left(pair));
      set_right_trunc_flag(pair, reject_right(pair));
      set_left_trunc_flag(pair->reversed_pair, reject_left(pair->reversed_pair));
      set_right_trunc_flag(pair->reversed_pair, reject_right(pair->reversed_pair));
/*      if (reject_left(pair) || reject_right(pair)) notify("*"); */
    }
  }
  notify(" Done\n");
}
/*
find_rejects()
{
  Align_info *align_entry, *align_entry2;
  int entry1;
  int i, j, length1, length2, n_hq_discrep, n_ss_discrep, hq_discrep, ss_discrep;
  int s1, s2;
  Aligned_pair *pair;
  unsigned char *diff;
  int q1, q2;
  int cum_score, max_score;
  int site1, site2, l_site1, l_site2;
  int max_sum_hq, max_sum_ss, sum_hq, sum_ss, d_flag;
  int LLR_hist[2][50];
  int LLR_score;
  int llr_score;

  notify("Finding rejects ... ");

  for (j = 0; j < 2; j++)
    for (i = 0; i < 50; i++) 
      LLR_hist[j][i] = 0;
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    align_entry = align_array + entry1;
    if (is_anomalous(align_entry)) continue;
    length1 = align_entry->db_entry->length - 1;
    for (pair = align_entry->first_pair; pair; pair = pair->next) {
      if (!is_best(pair)) continue; 
      if (pair->entry1 >= pair->entry2) continue; 
      if (is_reject_self(pair)) continue;
      set_left_trunc_flag(pair, reject_left(pair));
      set_right_trunc_flag(pair, reject_right(pair));
      set_left_trunc_flag(pair->reversed_pair, reject_left(pair->reversed_pair));
      set_right_trunc_flag(pair->reversed_pair, reject_right(pair->reversed_pair));
      align_entry2 = align_array + pair->entry2;
      if (is_anomalous(align_entry2)) continue;
      length2 = align_entry2->db_entry->length - 1;

      n_hq_discrep = n_ss_discrep = 0;

      site1 = pair->start1 - 1;
      site2 = pair->start2 - 1;
      
      for (diff = pair->diffs; *(diff + 1); diff++) {
	site1 += diff_gap1(*diff);
	site2 += diff_gap2(*diff);
	if (diff_type(*diff) == 'M') continue;
	serious_discrep(pair, site1, site2, &hq_discrep, &ss_discrep);
	n_hq_discrep += hq_discrep; 
	n_ss_discrep += ss_discrep;
      }

      d_flag = n_hq_discrep > parameters->max_discrep 
	|| n_ss_discrep > parameters->ss_max_discrep;
      set_reject_mismatch_flag(pair, d_flag);
      set_reject_mismatch_flag(pair->reversed_pair, d_flag);

 now count number of hq or ss mismatches (if any) that would have to exist past ends 
   of alignment to account for non-positive score there. Note assumes
   -2 is mismatch penalty. Does not take N's or possible indels into account (so is 
   overly conservative) 

      max_sum_hq = max_sum_ss = sum_hq = sum_ss = 0;
      for (site1 = pair->end1 + 1, site2 = pair->end2 + 1; 
	   site1 <= length1 && site2 <= length2 && (sum_hq > -20 || sum_ss > -20); 
	   site1++, site2++) {
	serious_discrep(pair, site1, site2, &hq_discrep, &ss_discrep);
	sum_hq += hq_discrep ? 1 : -2;
	sum_ss += ss_discrep ? 1 : -2;
	if (sum_hq > max_sum_hq) max_sum_hq = sum_hq;
	if (sum_ss > max_sum_ss) max_sum_ss = sum_ss;
      }
      if (max_sum_hq > 0) n_hq_discrep += (max_sum_hq + 2) / 3;
      if (max_sum_ss > 0) n_ss_discrep += (max_sum_ss + 2) / 3;

      max_sum_hq = max_sum_ss = sum_hq = sum_ss = 0;
      for (site1 = pair->start1 - 1, site2 = pair->start2 - 1; 
	   site1 >= 0 && site2 >= 0 && (sum_hq > -20 || sum_ss > -20); 
	   site1--, site2--) {
	serious_discrep(pair, site1, site2, &hq_discrep, &ss_discrep);
	sum_hq += hq_discrep ? 1 : -2;
	sum_ss += ss_discrep ? 1 : -2;
	if (sum_hq > max_sum_hq) max_sum_hq = sum_hq;
	if (sum_ss > max_sum_ss) max_sum_ss = sum_ss;
      }
      if (max_sum_hq > 0) n_hq_discrep += (max_sum_hq + 2) / 3;
      if (max_sum_ss > 0) n_ss_discrep += (max_sum_ss + 2) / 3;

      d_flag = n_hq_discrep > parameters->max_discrep || n_ss_discrep > parameters->ss_max_discrep;
      LLR_score = pair->LLR_score;
      llr_score = LLR_score / 10.0 + 10;
      if (llr_score < 0) llr_score = 0;
      if (llr_score > 20) llr_score = 20;
      LLR_hist[d_flag][llr_score] += 1;
      set_reject_total_flag(pair, d_flag);
      set_reject_total_flag(pair->reversed_pair, d_flag);

      pair->reversed_pair->hq_discrep = pair->hq_discrep = n_hq_discrep > 99 ? 99 : n_hq_discrep;
      pair->reversed_pair->ss_discrep = pair->ss_discrep = n_ss_discrep > 99 ? 99 : n_ss_discrep;

      l_site1 = site1 = pair->start1 - 1;
      l_site2 = site2 = pair->start2 - 1;
      
      for (diff = pair->diffs; !d_flag && *diff; diff++) {
	site1 += diff_gap1(*diff);
	site2 += diff_gap2(*diff);
	if (diff_type(*diff) == 'M') continue;
	for (s1 = l_site1 + 1, s2 = l_site2 + 1; 
	     s1 < site1 || s2 < site2; s1++, s2++) {
	  serious_discrep(pair, s1, s2, &hq_discrep, &ss_discrep);
	  n_hq_discrep += hq_discrep; 
	  n_ss_discrep += ss_discrep;
	  if (n_hq_discrep > 5 * (1 + parameters->max_discrep)
	      || n_ss_discrep > 5 * (1 + parameters->ss_max_discrep)) {
	    d_flag = 1;
	    break;
	  }
	}
	l_site1 = site1;
	l_site2 = site2;
      }
      set_rejectable_flag(pair, d_flag);
      set_rejectable_flag(pair->reversed_pair, d_flag);
    }
  }
  printf("\n\nLLR histogram (nonrejected, rejected pairs): ");
  for (i = 0; i < 50; i++)
    if (LLR_hist[0][i] || LLR_hist[1][i])
      printf("\n%5.1f  %4d    %4d", i - 10.0, LLR_hist[0][i], LLR_hist[1][i]);
  notify(" Done\n");
}
*/

/* find alignments that appear to be completely vector, and eliminate them 
*/

elim_complete_vector_matches(db) 
   Database *db;
{
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  int entry1, j;
  Aligned_pair *pair, *pair_prev;
  Aligned_pair *get_aligned_pairs();
  Seq_entry *get_seq_entry();
  unsigned char *get_seq();
  unsigned char *seq;
  int v_bound, end1, start1, end2, start2, length2;

  notify("Finding vector matches ... ");

  v_bound = parameters->vector_bound;
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry1 = get_align_entry(entry1);
    for (pair = get_aligned_pairs(entry1), pair_prev = 0; pair; pair = pair->next) {
      start1 = pair->start1;
      end1 = pair->end1;
      if (is_reverse(pair)) {
	length2 = get_seq_length(pair->entry2) - 1;
	start2 = length2 - pair->end2;
	end2 = length2 - pair->start2;
      }
      else {
	start2 = pair->start2;
	end2 = pair->end2;
      }
      if (!is_reverse(pair) && end1 < v_bound && end2 < v_bound
	  || is_reverse(pair) && (end1 < v_bound && start2 > v_bound + 100
				  || end2 < v_bound && start1 > v_bound + 100) ) {
/*
	       && (!align_entry->segments || align_entry->segments->end < 55 
	       || align_entry->segments->start > pair->end1
	       || !align_entry2->segments || align_entry2->segments->end < 55 
	       || align_entry2->segments->start > pair->end2 ? ' ' : '*')
*/
/*
	if (end1 < v_bound) {
	  if (end1 > align_entry1->last_vec) 
	    align_entry1->last_vec = end1;
	  if (start1 < align_entry1->first_vec) 
	    align_entry1->first_vec = start1;
	}
	if (end2 < v_bound) {
	  align_entry2 = get_align_entry(pair->entry2);
	  if (end2 > align_entry2->last_vec) 
	    align_entry2->last_vec = end2;
	  if (start2 < align_entry2->first_vec) 
	    align_entry2->first_vec = start2;
	}
*/
	if (pair_prev) pair_prev->next = pair->next;
	else get_seq_entry(entry1)->aligned_pairs = pair->next;
/* notify("decrement"); */
	num_pairs--;
/*
	set_reject_vector_flag(pair, 1);
	set_reject_vector_flag(pair->reversed_pair, 1);
*/
      }
      else pair_prev = pair;
    }
    if (parameters->screen) {
      seq = get_seq(entry1);
      for (j = 0; j <= align_entry1->last_vec; j++) seq[j] = 'X';
/* N.B. THIS ASSUMES SEQUENCE IS IN MEMORY -- IF NOT IT WILL FAIL !! */
    }
  }
/* look for other matches involving putative vector segment ? */
  notify(" Done\n");
}

find_more_vector(db)
     Database *db;
{
  Align_info *align_entry1;
  Align_info *get_align_entry();
  int entry1;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Segment *insert_segment();
  Segment *find_gap_segment_list(), *insert_segment();
  Segment *novec_segs, *diff_segs, *segment;
  int v_bound;
  
  notify("Finding more vector ... "); 

  v_bound = parameters->vector_bound;
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    mark_save_block();
    align_entry1 = get_align_entry(entry1);
    novec_segs = 0;
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (!same_subclone(pair->entry1, pair->entry2) &&
	  (is_reverse(pair) || pair->start2 >= v_bound)) 
	novec_segs = insert_segment(novec_segs, pair->start1, pair->end1);
    }
    diff_segs = find_gap_segment_list(novec_segs, align_entry1->segments);
    if (diff_segs && diff_segs->end < parameters->vector_bound) {
      if (align_entry1->first_vec > diff_segs->start)
	align_entry1->first_vec = diff_segs->start;
      for (segment = diff_segs; 
	   segment && segment->end < parameters->vector_bound; 
	   segment = segment->next)
	if (align_entry1->last_vec < segment->end)
	  align_entry1->last_vec = segment->end > parameters->vector_bound ? parameters->vector_bound : segment->end;
    }
    free_seg_blocks();
  }
/* look for other matches involving putative vector segment ? */
  notify(" Done\n");
}

/*  */
print_vector(db)
     Database *db;
{
  Align_info *align_entry1;
  Align_info *get_align_entry();
  int entry1, i, orig_i, vector_flag, start, n_X;
  char *get_id(), *get_orig_qual(); 
  unsigned char *get_seq();
  unsigned char *seq;
  char *o_qual;
  int length;

  printf("\n\nProbable unremoved sequencing vector (matches excluded from assembly, quality reduced to 0): ");
  vector_flag = 0;
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry1 = get_align_entry(entry1);
    seq = get_seq(entry1);
    if (align_entry1->last_vec > 0) {
      vector_flag = 1;
      printf("\n%-12s %2d-%2d   ", get_id(entry1), 
	     align_entry1->first_vec + 1, align_entry1->last_vec + 1);
      for (i = align_entry1->first_vec; i <= align_entry1->last_vec; i++) 
	printf("%c", seq[i]);
/*
      for (pair = align_entry1->first_pair; pair; pair = pair->next) 
	if (pair->start1 < align_entry1->last_vec - 5) 
	  printf("\n  Unremoved overlapping match: %d-%d  with  %c %s %d-%d",
		 1 + pair->start1, 1 + pair->end1, is_reverse(pair) ? 'C' : ' ',
		 align_array[pair->entry2].db_entry->id, 
		 1 + pair->start2, 1 + pair->end2);
*/
    }
/* also set orig_qual to 0 in initial part of read, before first_start */
    start = get_seq_length(entry1);
    length = start - 1;
    if (start > parameters->vector_bound) start = parameters->vector_bound;
    
    if (align_entry1->segments && align_entry1->segments->start < start)
      start = align_entry1->segments->start;

    if (start <= align_entry1->last_vec) start = align_entry1->last_vec + 1;
/*
    if (start < align_entry1->qual_start && align_entry1->qual_start < length) 
      start = align_entry1->qual_start;
*/
    if (align_entry1->last_vec > length) fatalError("length");
    o_qual = get_orig_qual(entry1);
    for (i = 0; i < start && o_qual[i] != 99; i++); 

    orig_i = i;
    for (n_X = 0; seq[i] && n_X < 20; i++) {
      if (seq[i] != 'X') n_X++;
      else n_X = 0;
    }
    if (seq[i]) i -= 20;
     
    align_entry1->last_vec = i - 1; 

/*
    if (i > orig_i) printf("\n%-12s  %d add'l X's", get_id(entry1), i - orig_i);   
*/
  }
  if (!vector_flag) printf(" None.");
}

/*
serious_discrep(pair, site1, site2, add_hq_discrep, add_ss_discrep)
     Aligned_pair *pair;
     int site1, site2;
     int *add_hq_discrep, *add_ss_discrep;
{
  int q1, q2;

  *add_hq_discrep = *add_ss_discrep = 0;
  if (is_reverse(pair)) {
    site2 = align_array[pair->entry2].db_entry->length - 1 - site2;
  }
  q1 = align_array[pair->entry1].adj_qual[site1];
  q2 = align_array[pair->entry2].adj_qual[site2];
  if (q1 >= parameters->qual_cut && q2 >= parameters->qual_cut) *add_hq_discrep = 1;
  if (is_reverse(pair) || different_chemistry(pair->entry1, pair->entry2)) { 
    if (q1 >= parameters->ss_qual_cut && q2 >= parameters->qual_cut
	|| q1 >= parameters->qual_cut && q2 >= parameters->ss_qual_cut)
      *add_ss_discrep = 1;
  }
  else 
  if (q1 >= parameters->ss_qual_cut && q2 >= parameters->ss_qual_cut)
    *add_ss_discrep = 1;
}
*/
/*
print_reject_summary()
{
  Align_info *align_entry, *align_entry2;
  Aligned_pair *pair;
  int entry1;
  int reject_hist[2][100], max_score[2][100];
  int ss_reject_hist[2][100], ss_max_score[2][100];
  int nonrejectable_hist[2][300];
  int i, j, cum0, cum1;
  char has_pair, has_non_rej, set_flag, has_left_right;

  for (j = 0; j < 2; j++)
    for (i = 0; i < 100; i++) 
      reject_hist[j][i] = max_score[j][i] = ss_reject_hist[j][i] = ss_max_score[j][i] = 0;
  for (j = 0; j < 2; j++)
    for (i = 0; i < 300; i++) nonrejectable_hist[j][i] = 0;

  set_flag = 0;
  printf("\n\nOther reads with all matches rejected: ");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    align_entry = align_array + entry1;
    if (is_anomalous(align_entry)) continue;
    has_pair = has_non_rej = 0;
    for (pair = align_entry->first_pair; pair; pair = pair->next) {
      if (!is_best(pair)) continue;
      if (is_reject_self(pair)) continue;
      has_pair = 1;
      if (!is_reject_any(pair)) has_non_rej = 1;
      if (pair->entry1 >= pair->entry2) continue; 
      if (is_anomalous(align_array + pair->entry2)) continue;
      if (!is_rejectable(pair))
	nonrejectable_hist[is_reverse(pair)][pair->score > 299 ? 299 : pair->score] += 1;
      reject_hist[is_reverse(pair)][pair->hq_discrep] += 1;
      if (pair->score > max_score[is_reverse(pair)][pair->hq_discrep])
	max_score[is_reverse(pair)][pair->hq_discrep] = pair->score;
      if (pair->hq_discrep <= parameters->max_discrep) {
	ss_reject_hist[is_reverse(pair)][pair->ss_discrep] += 1;
	if (pair->score > ss_max_score[is_reverse(pair)][pair->ss_discrep])
	  ss_max_score[is_reverse(pair)][pair->ss_discrep] = pair->score;
      }
    }
    if (has_pair && !has_non_rej) {
      set_flag = 1;
      printf("\n%s", align_entry->db_entry->id);
    }
  }
  if (!set_flag) printf("None.");

  printf("\n\nHQ mismatch histogram\nhq_mismatches, # fwd pairs (max_score), cum, # rev pairs (max_score), cum");
  cum0 = cum1 = 0;
  
  for (i = 99; i >= 0; i--) {
    if (i == parameters->max_discrep) printf("\n****************CUTOFF*******************");
    if (reject_hist[0][i] || reject_hist[1][i]) {
      cum0 += reject_hist[0][i];
      cum1 += reject_hist[1][i];
      printf("\n%2d  %4d (%3d)  %4d    %4d (%3d)  %4d", i, 
	     reject_hist[0][i], max_score[0][i], cum0, 
	     reject_hist[1][i], max_score[1][i], cum1);
    }
  }
  printf("\n\nSS mismatch histogram (non-hq-rejected pairs)\nss_mismatches, # fwd pairs (max_score), cum, # rev pairs (max_score), cum");
  cum0 = cum1 = 0;
  
  for (i = 99; i >= 0; i--) {
    if (i == parameters->ss_max_discrep) printf("\n****************CUTOFF*******************");
    if (ss_reject_hist[0][i] || ss_reject_hist[1][i]) {
      cum0 += ss_reject_hist[0][i];
      cum1 += ss_reject_hist[1][i];
      printf("\n%2d  %4d (%3d)  %4d    %4d (%3d)  %4d", i, 
	     ss_reject_hist[0][i], ss_max_score[0][i], cum0, 
	     ss_reject_hist[1][i], ss_max_score[1][i], cum1);
    }
  }

  printf("\n\nUnrejectable pair histogram\nscore # fwd pairs, cum, # rev pairs , cum");
  cum0 = cum1 = 0;
  
  for (i = 299; i >= 0; i--) {
    if (nonrejectable_hist[0][i] || nonrejectable_hist[1][i]) {
      cum0 += nonrejectable_hist[0][i];
      cum1 += nonrejectable_hist[1][i];
      printf("\n%2d  %4d   %4d     %4d   %4d", i, 
	     nonrejectable_hist[0][i], cum0, 
	     nonrejectable_hist[1][i], cum1);
    }
  }


}
*/
/* reject alignment due to incomplete extent to left */
reject_left(pair)
  Aligned_pair *pair;
{
  int length2;
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();

  align_entry = get_align_entry(pair->entry1);
  align_entry2 = get_align_entry(pair->entry2);

/* for reverse sense alignments, allow more limited extent on one of the two strands */
  if (is_reverse(pair)) {
    length2 = get_seq_length(pair->entry2) - 1;
    return 
      pair->start1 - align_entry->first_start > parameters->maxgap
	&& pair->start2 - (length2 - align_entry2->last_end) > parameters->maxgap
	  || pair->start1 - align_entry->first_start > parameters->maxgap
	    && pair->start2 - (length2 - align_entry2->last_end) > parameters->maxgap;
  }
  else
    return 
      pair->start1 - align_entry->first_start > parameters->maxgap
	&& pair->start2 - align_entry2->first_start > parameters->maxgap;
}

/* reject alignment due to incomplete extent to right */
reject_right(pair)
  Aligned_pair *pair;
{
  int length2;
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();

  align_entry = get_align_entry(pair->entry1);
  align_entry2 = get_align_entry(pair->entry2);

  if (is_reverse(pair)) {
    length2 = get_seq_length(pair->entry2) - 1;

    return align_entry->last_end - pair->end1 > parameters->maxgap
      && (length2 - align_entry2->first_start) - pair->end2 > parameters->maxgap
    || align_entry->last_end - pair->end1 > parameters->maxgap
      && (length2 - align_entry2->first_start) - pair->end2 > parameters->maxgap;
  }
  else
    return align_entry->last_end - pair->end1 > parameters->maxgap
      && align_entry2->last_end - pair->end2 > parameters->maxgap;
}

/*
debug_alignments()
{
  Aligned_pair *pair;
  Profile *make_profile_from_seq();
  Profile *q_profile;
  int temp, off1, off2, entry1;
  Align_info *align_entry1, *align_entry2;
  int i, first_start2, last_end2, start_gap, end_gap, length2;
  int orig_score;
  unsigned char *seq2;
  Db_entry *db_entry;

  printf("\n\n");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    align_entry1 = align_array + entry1;
    if (is_anomalous(align_entry1)) continue;
    for (pair = align_entry1->first_pair; pair; pair = pair->next) {
      if (pair->score < 250) continue;
**      if (pair->hq_mismatch != 32 || pair->ss_mismatch != 32) continue; **
**      if (pair->hq_mismatch & 31) > 0 || pair->ss_mismatch < 5) continue; **
      if (pair->entry2 <= pair->entry1) continue;
      
      align_entry2 = align_array + pair->entry2;
      if (is_anomalous(align_entry2)) continue;

      length2 = align_entry2->db_entry->length - 1;

      first_start2 = is_reverse(pair) ? length2 - align_entry2->last_end : align_entry2->first_start;
      last_end2 = is_reverse(pair) ? length2 - align_entry2->first_start : align_entry2->last_end;
      start_gap = pair->start1 - align_entry1->first_start;
      if (start_gap > pair->start2 - first_start2) start_gap = pair->start2 - first_start2;
      end_gap = align_entry1->last_end - pair->end1;
      if (end_gap > last_end2 - pair->end2) end_gap = last_end2 - pair->end2;

      printf("\n\n %d  hq: %d ss:%d    %d %d", pair->score, (int)pair->hq_discrep,  (int)pair->ss_discrep, start_gap, end_gap);
      set_score_mat(0, -2, -2, 0, -1);
      q_profile = make_profile_from_seq((Profile *)0, align_entry1->db_entry->seq, align_entry1->db_entry->length, 0);  
      printf("\n%s (Query)", align_entry1->db_entry->id);
      if (start_gap > 0) {
	printf("\nLeading: ");
	for (i = 0; i < pair->start1; i++) 
	  printf("%c", align_entry1->db_entry->seq[i]);
      }
      if (end_gap > 0) {
	printf("\nTrailing: ");
	for (i = pair->end1 + 1; i < align_entry1->db_entry->length; i++) 
	  printf("%c", align_entry1->db_entry->seq[i]);
      }
      seq2 = (align_entry2->db_entry + (is_reverse(pair) ? num_query_entries : 0))->seq;
 
      printf("\n%s %c", align_entry2->db_entry->id, is_reverse(pair) ? 'C' : ' ');
      if (start_gap > 0) {
	printf("\nLeading: ");
	for (i = 0; i < pair->start2; i++) printf("%c", seq2[i]);
      }
      if (end_gap > 0) {
	printf("\nTrailing: ");
	for (i = pair->end2 + 1; i < align_entry2->db_entry->length; i++) 
	  printf("%c", seq2[i]);
      }

      off1 = pair->start1 - pair->start2;
      off2 = pair->end1 - pair->end2;
      if (off1 > off2) {
	temp = off1;
	off1 = off2;
	off2 = temp;
      }
    db_entry = align_entry2->db_entry +
			  (is_reverse(pair) ? num_query_entries : 0);
      full_smith_waterman(q_profile, db_entry->seq, db_entry->length,
			  off1 - parameters->bandwidth, 
			  off2 + parameters->bandwidth, 0, 0, 0, 0, parameters->minscore, &orig_score);
      print_alignment();
    }
  }
}
*/
/* return 1 if pair is rejected on basis of mismatches (rather than incomplete alignment) */

/*
mismatch_reject(pair)
     Aligned_pair *pair;
{
  return (pair->hq_mismatch & 31) > parameters->max_discrep
    || (pair->ss_mismatch & 31) > parameters->ss_max_discrep;
}

any_reject(pair)
     Aligned_pair *pair;
{
  return pair->hq_mismatch > parameters->max_discrep
    || pair->ss_mismatch > parameters->ss_max_discrep;
}
*/

/* find exact duplicate reads, mark them, and remove pairs involving them. (no longer removes pairs -- because
   now using tree structure rather than linked list. So duplicate entries are marked only)
   Requires reads to be stored in memory */
int dup_flag;
find_duplicates(db)
     Database *db;
{
  Cand_pair *pair, *pair_prev, *cand_pairs;
  Cand_pair *get_cand_pairs();
  int entry1, entry2, elim_flag;
  unsigned char *seq, *seq2;
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  Seq_entry *get_seq_entry();
  unsigned char *get_seq(), *get_comp_seq();
  char *get_id();
  
  Segment *segment;

  dup_flag = 0;
  notify("Finding duplicates ... ");
  printf("\n\nExact duplicate reads: ");
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry1 = get_align_entry(entry1);
    if (!(cand_pairs = get_cand_pairs(entry1))) continue;
    visit_cand_pairs_dups(entry1, cand_pairs);
    /*
    for (pair = cand_pairs, pair_prev = 0; pair; pair = pair->next) { 
      entry2 = pair->entry2;
      align_entry2 = get_align_entry(entry2);
      elim_flag = 0;
      if (entry1 == entry2) goto elim;
      if (is_duplicate(align_entry2)) {
	elim_flag = 1;
	goto elim;
      }
      /* look for segment spanning 0 -- which must be present if there is to be a word match
	 on the main diagonal */
    /*
      for (segment = pair->band_segments; segment; segment = segment->next) 
	if (segment->start < 0 && segment->end > 0) break;

      if (!segment) goto elim;
      seq = get_seq(entry1);
      seq2 = pair->reverse ? get_comp_seq(entry2) : get_seq(entry2);
      for ( ; *seq && *seq == *seq2; seq++, seq2++);
      if (*seq == *seq2) {
	printf("\n%-15s   %c  %-15s (perfect)", 
	       get_id(entry1), pair->reverse ? 'C' : ' ', get_id(entry2));
	if (!parameters->retain_duplicates) {
	  elim_flag = 1;
	  set_duplicate(align_entry2);
	  get_seq_entry(entry2)->cand_pairs = 0;
	}

      }
    elim:
      if (elim_flag) {
	if (pair_prev) pair_prev->next = pair->next;
	else get_seq_entry(entry1)->cand_pairs = pair->next;
	dup_flag = 1;
      }
      else pair_prev = pair;
    }
    */
  }
  printf(dup_flag ? "\n\n2d read in each pair excluded from assembly." : " None.");
  /* now remove pairs from previous entries */
  /*
  if (dup_flag) {
    for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
      for (pair = get_cand_pairs(entry1), pair_prev = 0; 
	 pair; pair = pair->next) {
	if (is_duplicate(get_align_entry(pair->entry2))) {
	  if (pair_prev) pair_prev->next = pair->next;
	  else get_seq_entry(entry1)->cand_pairs = pair->next;
	}
	else pair_prev = pair;
      }
    }
  }
  */
  notify(" Done\n");
}

visit_cand_pairs_dups(entry1, pair)
     int entry1;
     Cand_pair *pair;
{
  int entry2, elim_flag;
  Align_info *align_entry2;
  Align_info *get_align_entry();
  unsigned char *seq, *seq2;
  Segment *segment;
  unsigned char *get_seq(), *get_comp_seq();
  char *get_id();
  Seq_entry *get_seq_entry();

  if (!pair) return;

  entry2 = pair->entry2;
  align_entry2 = get_align_entry(entry2);
  elim_flag = 0;
  if (entry1 == entry2) goto elim;
  if (is_duplicate(align_entry2)) {
    elim_flag = 1;
    goto elim;
  }
  /* look for segment spanning 0 -- which must be present if there is to be a word match
     on the main diagonal */
  for (segment = pair->band_segments; segment; segment = segment->next) 
    if (segment->start < 0 && segment->end > 0) break;

  if (!segment) goto elim;
  seq = get_seq(entry1);
  seq2 = pair->reverse ? get_comp_seq(entry2) : get_seq(entry2);
  for ( ; *seq && *seq == *seq2; seq++, seq2++);
  if (*seq == *seq2) {
    printf("\n%-15s   %c  %-15s (perfect)", 
	   get_id(entry1), pair->reverse ? 'C' : ' ', get_id(entry2));
    if (!parameters->retain_duplicates) {
      elim_flag = 1;
      set_duplicate(align_entry2);
      get_seq_entry(entry2)->cand_pairs = 0;
    }
  }
 elim:
  if (elim_flag) {
    dup_flag = 1;
  }
  visit_cand_pairs_dups(entry1, pair->left);
  visit_cand_pairs_dups(entry1, pair->right);
}

find_near_duplicates()
{
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int entry1, dup_flag;
  char *get_id();

  dup_flag = 0;
  notify("Finding near duplicates ... ");
  printf("\n\nNear duplicate reads: ");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (pair->entry2 <= entry1) continue;
      if (same_subclone(pair->entry1, pair->entry2) || is_reverse(pair)) continue; 
      if (abs(pair->start1 - pair->start2) < 10 
	  && pair->start1 < 50 && pair->end1 > get_seq_length(entry1) - 50
	  && pair->start2 < 50 && pair->end2 > get_seq_length(pair->entry2) - 50) {
	printf("\n%-15s   %c   %-15s (imperfect: %d-%d (%d)   %d-%d (%d) )", 
	       get_id(entry1), is_reverse(pair) ? 'C' : ' ', get_id(pair->entry2),
	       pair->start1, pair->end1, get_seq_length(entry1) - 1 - pair->end1, 
	       pair->start2, pair->end2, get_seq_length(pair->entry2) - 1 - pair->end2);
	dup_flag = 1;
	/*
	  set_duplicate(align_array + pair->entry2);
	  */
      }
    }
  }
  printf(dup_flag ? "" : " None.");
  notify(" Done\n");
}

find_blocked_reads(db)
     Database *db;
{
  int entry1, n_blocked_left, n_blocked_right, n_blocked_both, n_0_qual;
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  char *get_id();
  int left_flag, right_flag;
  int length2, qual_start1, qual_end1, temp, qual_start2, qual_end2;

  notify("Finding blocked reads ... ");
  printf("\n\nBlocked reads: ");
  n_blocked_left = n_blocked_right = n_blocked_both = n_0_qual = 0;
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    if (!get_aligned_pairs(entry1)) continue;
    align_entry1 = get_align_entry(entry1);
    left_flag = right_flag = 1;
    
    find_trimmed_quals(align_entry1, &qual_start1, &qual_end1);
    if (qual_start1 >= qual_end1) {
      set_lowqual(align_entry1);
      n_0_qual++; 
      for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
	set_reject_qual_flag(pair, 1);
	set_reject_qual_flag(pair->reversed_pair, 1);
      }
      continue;
    }

    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (is_reject_qual(pair)) continue;
      align_entry2 = get_align_entry(pair->entry2);
      find_trimmed_quals(align_entry2, &qual_start2, &qual_end2);
      if (qual_start2 >= qual_end2)  {
	set_reject_qual_flag(pair, 1);
	set_reject_qual_flag(pair->reversed_pair, 1);
	continue;
      }
      if (is_reverse(pair)) {
	length2 = get_seq_length(pair->entry2) - 1;
	temp = length2 - qual_start2;
	qual_start2 = length2 - qual_end2;
	qual_end2 = temp;
      }
      qual_start2 += pair->offset;
      qual_end2 += pair->offset;
      if (qual_start2 > qual_end1 - 10 || qual_end2 < qual_start1 + 10)  {
	set_reject_qual_flag(pair, 1);
	set_reject_qual_flag(pair->reversed_pair, 1);
	continue;
      }
      if (left_flag && pair->start1 < qual_start1 + 15 && qual_start2 < qual_start1) 
	left_flag = 0;
      if (right_flag && pair->end1 > qual_end1 - 15 && qual_end2 > qual_end1) 
	right_flag = 0;
    }
    
    if (left_flag || right_flag) {
      set_blocked(align_entry1);
      if (left_flag && !right_flag) n_blocked_left++;
      else if (!left_flag && right_flag) n_blocked_right++;
      else n_blocked_both++;
      printf("\n%s %d %d  %s %s", get_id(entry1), qual_start1, qual_end1, 
	     left_flag ? "left" : "", right_flag ? "right" : "");
    }
  }
  printf("\n\n%d blocked reads: %d left only, %d right only, %d both.",
	 n_blocked_left + n_blocked_right + n_blocked_both, 
	 n_blocked_left, n_blocked_right, n_blocked_both);
  printf("\n%d reads (not shown) lack a high-quality segment.", n_0_qual);
  notify(" Done\n");
}

find_chimeras(db)
     Database *db;
{
  Align_info *align_entry;
  Align_info *get_align_entry();
  Segment *segment;
  Aligned_pair *pair, *pairs;
  Aligned_pair *get_aligned_pairs();
  char left_partial_match, right_partial_match, x_flag;
  int n_chimeras, entry1, n_doubles;
  int max_seg_size, max_score, seg_start, seg_end, i_bit;
  int f_max_score, f_seg_start, f_seg_end, i;
  unsigned char *seq;
  unsigned char *get_seq();
  char *get_id();

  notify("Finding multi-segment reads ... ");
  printf("\n\nMulti-segment reads (initially rejected segments in parentheses) -- XXX means segments flank X'd region: ");
  n_chimeras = n_doubles = 0;

  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
/* require two non-overlapping confirmed segments */
    if (!align_entry->segments || !align_entry->segments->next) continue;
    seq = get_seq(entry1);
    pairs = get_aligned_pairs(entry1);
    x_flag = 0;
    for (segment = align_entry->segments; segment->next; segment = segment->next) 
      for (i = segment->end + 1; i < segment->next->start; i++)
	if (seq[i] == 'X') {
	  x_flag = 1; /* there are x's between two matching segments -- so reject pairs in least favorable
			 one */
	  break; 
	}
    n_doubles++;
    max_seg_size = 0;
    max_score = f_max_score = 0;
    seg_start = seg_end = f_seg_start = f_seg_end = -1;
    for (pair = pairs; pair; pair = pair->next) {
      if (pair->entry1 == pair->entry2 || is_reject_chimeric(pair) 
	  || is_reject_node(pair)) continue;
      if (same_subclone(pair->entry1, pair->entry2) 
/* FOLLOWING TAKEN OUT 7/20/98 */
	  /* && !different_chemistry(pair->entry1, pair->entry2) */) continue; 
      if (is_reverse(pair)) {
	if (pair->score > max_score) {
	  max_score = pair->score;
	  seg_start = pair->start1;
	  seg_end = pair->end1;
	}
      }
      else {
	if (pair->score > f_max_score) {
	  f_max_score = pair->score;
	  f_seg_start = pair->start1;
	  f_seg_end = pair->end1;
	}
      }
    }
    if (!max_score) {
      seg_start = f_seg_start;
      seg_end = f_seg_end;
    }
    for (segment = align_entry->segments; segment; segment = segment->next) 
      if (seg_start >= segment->start && seg_end <= segment->end) {
	seg_start = segment->start;
	seg_end = segment->end;
	break;
      }
    for (segment = align_entry->segments, i_bit = 1; segment; 
	 segment = segment->next, i_bit *= 2) 
      if (segment->start != seg_start) {
	align_entry->chimera_bits |= i_bit;
	for (pair = pairs; pair; pair = pair->next) 
	  if (pair->start1 >= segment->start && pair->end1 <= segment->end) {
	    set_reject_chimeric_flag(pair, 1);
	    set_reject_chimeric_flag(pair->reversed_pair, 1);
	    if (x_flag) {
	      set_reject_vector_flag(pair, 1);
	      set_reject_vector_flag(pair->reversed_pair, 1);
	    }
	  }
      }
    printf("\n%-12s  %s  ", get_id(entry1), x_flag ? "XXX" : "   ");
    print_segs(align_entry);
/*      if (segment->next->start > segment->end + MAX_CHIMERA_GAP) continue; */

      left_partial_match = right_partial_match = 1;
/*
      for (pair = pairs; pair; pair = pair->next) {
	if (pair->entry1 == pair->entry2 || is_reject_self(pair)) continue;
	if (same_subclone(pair->entry1, pair->entry2) 
	   && !different_chemistry(pair->entry1, pair->entry2)) continue; 

	if (!left_partial_match  && is_right_trunc(pair) 
	    && pair->start1 >= segment->start && pair->end1 <= segment->end) {
	      left_partial_match = 1;
	}

	if (!right_partial_match  && is_left_trunc(pair) 
	    && pair->start1 >= segment->next->start && pair->end1 <= segment->next->end) 
	  right_partial_match = 1;
      } 
*/
/*
      n_chimeras++;
      set_chimera(align_entry); 
      else n_doubles++;
*/
  }
/*
  if (n_chimeras) printf("\n\n%d probable chimeric reads.", n_chimeras);
  else printf("None.");
*/
  printf("\n\n%d reads with multiple segments.", n_doubles);
  notify(" Done\n");
}
 
print_segs(align_entry)
  Align_info *align_entry;
{
  Segment *segment;
  int i_bit;

  for (segment = align_entry->segments, i_bit = 1; segment; 
       segment = segment->next, i_bit *= 2) 
    printf(i_bit & align_entry->chimera_bits ? " (%d %d) " : " %d %d ", 
	   segment->start + 1, segment->end + 1);

}

find_self_matches()
{
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int entry1, i;
  char flag;
  char *get_id();
  unsigned char *get_seq();
  unsigned char *seq;

  notify("Finding self matches ... ");
  printf("\n\nInternal read matches (same orientation) : ");
  flag = 0;
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (!pair->score) continue;
      if (entry1 != pair->entry2) continue;
      set_reject_self_flag(pair, 1);
      if (is_reverse(pair)) continue;
      if (pair->start1 < pair->start2) {
	seq = get_seq(entry1);
	flag = 1;
	printf("\n%3d   %-12s ", pair->score, get_id(entry1));
	if (pair->start2 < pair->end1 + 5) {
	  printf(" tandem (");
	  if (pair->start2 - pair->start1 < 10) 
	    for (i = pair->start1; i < pair->start2; i++) 
	      printf("%c", seq[i]);
	  else printf("%d-mer", pair->start2 - pair->start1);
	  printf(")_%d    %3d-%3d ", 
		 (pair->end2 - pair->start1)/(pair->start2 - pair->start1), 
		 pair->start1, pair->end2);
	}
	else printf(" disjoint %d-mers  %d-%d / %d-%d ", 
		    pair->end1 - pair->start1 + 1,
		    pair->start1, pair->end1, pair->start2, pair->end2);
       }
    }
  }
  if (!flag) printf(" None.");
  notify(" Done\n");
} 


/* test chimeras for deletions; N.B. this will not find all deleted subclones!
   (many of them may not show up as chimeras) */
/* since doing full smith_waterman search, if either seq is over 5000 it
   is skipped */
old_find_deletions(db)
     Database *db;
{
  Align_info *align_entry;
  Align_info *get_align_entry();
  int entry1;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Profile *make_profile_from_seq();
  Profile *q_profile;
  int l_span, r_span, best_score, best_entry, n_tested;
  int score, start1, end1, start2, end2, mismatches, insertions, deletions;
  int orig_score;
  unsigned char *get_seq(), *get_comp_seq();
  char *get_id();
  unsigned char *seq2;

  set_score_mat(); /* need to revise call + test what is 
						    appropriate here;  */
  set_gap_penalties(2 * parameters->penalty, 0, 0, 0, 0);
  printf("\n\nDeletion reads:\n");
  n_tested = 0;
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    if (get_seq_length(entry1) > 5000) continue;
    if (!align_entry->segments || !align_entry->segments->next) continue;
    l_span = align_entry->segments->end - 20;
    r_span = align_entry->segments->next->start + 20;
    best_score = 0;
    best_entry = -1;
    q_profile = make_profile_from_seq((Profile *)0, get_seq(entry1), get_seq_length(entry1), 0); 
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
 /*     if (!is_best(pair)) continue;  check this */
      if (get_seq_length(pair->entry2) > 5000) continue;
     /* check whether entries have same name -- if so, don't use in test */
      if (same_subclone(pair->entry1, pair->entry2)) continue;
      seq2 = is_reverse(pair) ? get_comp_seq(pair->entry2) : get_seq(pair->entry2);
      score = full_smith_waterman(q_profile, seq2, get_seq_length(pair->entry2),
			 1, 0, /* full search -- so l_edge > r_edge */
			 0,0,0,0, parameters->minscore, &orig_score);
      get_stats(&start1, &end1, &start2, &end2,
		&mismatches, &insertions, &deletions);
      n_tested++;
      if (score >= parameters->minscore && start1 < l_span && end1 > r_span &&
	  score > best_score) {
	best_score = score;
	best_entry = pair->entry2;
      }
    }
    if (best_score) {
      printf("\n%-12s  (%d vs.  %s)", 
	     get_id(entry1), best_score, get_id(best_entry));
      set_deleted(align_entry);
    }
  }
  printf("\n\n%d alignments were tested", n_tested);
}

/* test sequence for signif. alignment with itself (to right of diagonal)
   note would not want to do this with long sequences */

test_self()
{
  Profile *make_profile_from_seq();
  Profile *q_profile;
  int entry1, score, length;
  int start1, end1, start2, end2;
  int mismatches, insertions, deletions;
  int full_smith_waterman();
  int orig_score;
  unsigned char *get_seq();
  char *get_id();
  unsigned char *seq;

  set_score_mat(); /* need to revise call */
  set_gap_penalties(parameters->penalty - 2, parameters->penalty - 1, parameters->penalty - 1, parameters->penalty - 1, 0);
/* N.B. NEED TO RESET PENALTIES FOLLOWING THIS!!! */
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    seq = get_seq(entry1);
    length = get_seq_length(entry1);
    q_profile = make_profile_from_seq((Profile *)0, seq, length, 0);
    score = full_smith_waterman(q_profile, seq, length,
		       1, length, 0,0,0,0, parameters->minscore, &orig_score);
    get_stats(&start1, &end1, &start2, &end2,
		       &mismatches, &insertions, &deletions);
    if (score >= parameters->minscore) {
      printf("\n%-12s %3d  %5d-%5d matches %5d-%5d",
	     get_id(entry1), score, start2, end2, start1, end1);
    }
    free_profile(q_profile);
  }
}

set_deleted(align_entry)
     Align_info *align_entry;
{
  align_entry->anomalies |= 1;
}

is_deleted(align_entry)
     Align_info *align_entry;
{
  return (int)(align_entry->anomalies & 1);
}

set_chimera(align_entry)
     Align_info *align_entry;
{
  align_entry->anomalies |= 2;
}

is_chimera(align_entry)
     Align_info *align_entry;
{
  return (int)(align_entry->anomalies & 2);
}

set_duplicate(align_entry)
     Align_info *align_entry;
{
  align_entry->anomalies |= 4;
}

is_duplicate(align_entry)
     Align_info *align_entry;
{
  return (int)(align_entry->anomalies & 4);
}

set_blocked(align_entry)
     Align_info *align_entry;
{
  align_entry->blocked |= 1;
}

is_blocked(align_entry)
     Align_info *align_entry;
{
  return (int)(align_entry->blocked & 1);
}

set_lowqual(align_entry)
     Align_info *align_entry;
{
  align_entry->blocked |= 2;
}

is_lowqual(align_entry)
     Align_info *align_entry;
{
  return (int)(align_entry->blocked & 2);
}



is_anomalous(align_entry)
     Align_info *align_entry;
{
  return (int)align_entry->anomalies;
}

