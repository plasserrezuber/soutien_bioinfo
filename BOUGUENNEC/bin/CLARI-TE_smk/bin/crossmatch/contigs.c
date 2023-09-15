/*****************************************************************************
#   Copyright (C) 1994-2000 by Phil Green.                          
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

#define LOCATION_FUDGE1 25 /* allows for uncertainty in exact position of reads in layout */
#define LOCATION_FUDGE2 25 /* larger -- for purposes of confirmation test */
 
extern Parameters *parameters;
extern int num_pairs, t_num_entries;
extern FILE *fp_log; /* log file */
extern Database *query_db;
Contig *contig_array;

static int num_contigs, num_singletons, num_duplicates; 

/* number of unordered aligned pairs
			 (= half the number of ordered pairs) */

/* contig initialization -- needs to be done after align_array has been created and
   pairs have been set */

alloc_contigs()
{
  char *our_alloc();

  contig_array = (Contig *)our_alloc(t_num_entries * sizeof(Contig));
}

init_contigs(strong_trim)
     int strong_trim;
{
  int i, j, j_start, hi_q;
  Align_info *align_entry;
  Align_info *get_align_entry();
  Contig *contig;
  char *get_adj_qual(), *get_id();
  char *adj_qual;
  unsigned char *seq;
  int qual_start, qual_end, length1;

  for (i = 0; i < t_num_entries; i++) {
    contig = contig_array + i;
    align_entry = get_align_entry(i);
    length1 = get_seq_length(i);
    adj_qual = get_adj_qual(i);

    if (strong_trim) {
      best_qual_seg(adj_qual, length1, &qual_start, &qual_end);
/* Following uses input trace qualities instead -- less useful when those are not
   available; but may have advantages in other situations (e.g. helps prevent
   incorporation of contaminants, whose adjusted quals may be significantly
   downgraded.
      find_trimmed_quals(align_entry, &qual_start, &qual_end);
*/
/* following may be too restrictive; unconfirmed reads will be trimmed too much */

/*
      if (qual_start < align_entry->first_start) 
	qual_start = align_entry->first_start;
*/

/* removed 980626
      if (!parameters->forcelevel) {
	if (qual_start > align_entry->rev_first_start) 
	  qual_start = align_entry->rev_first_start;
      }
*/
      contig->first_start = qual_start;

/* OLD VERSION:
      j_start = align_entry->last_vec + 1;
      if (j_start < align_entry->first_start) j_start = align_entry->first_start;
      for (j = j_start, hi_q = 0; j < align_entry->rev_first_start; j++) {
	hi_q += align_entry->adj_qual[j] >= 20;
	if (j - 10 >= j_start)
	  hi_q -= align_entry->adj_qual[j - 10] >= 20;
	if (hi_q >= 5) {
	  j -= 4;
	  break;
	}
      }
      contig->first_start = j; 
*/
      /* align_array[i].db_entry->length - 1;  */

/*
      if (qual_end > align_entry->last_end) qual_end = align_entry->last_end; 
*/

/* removed 980626
      if (!parameters->forcelevel) {
	if (qual_end < align_entry->rev_last_end) 
	  qual_end = align_entry->rev_last_end;
      }
*/
      contig->last_end = qual_end;

/* OLD VERSION
      j_start = align_entry->last_end;
      for (j = j_start, hi_q = 0; j > align_entry->rev_last_end; j--) {
 	hi_q += align_entry->adj_qual[j] >= 20;
	if (j + 10 <= j_start)
	  hi_q -= align_entry->adj_qual[j + 10] >= 20;
	if (hi_q >= 5) {
	  j -= 4;
	  break;
	}
      }
      contig->last_end = j; 
*/
    }
    else {
      contig->first_start = align_entry->first_start; 
      contig->last_end = align_entry->last_end; 
    }
    contig->tig_node = 0;
    contig->num_entries = 1;
    contig->num_matches = 0;
    contig->first = contig->last = align_entry;
    contig->index = 0;
    contig->n_pads = 0;
    contig->base_segment = 0;
    contig->comp_status = 0;
    contig->merge_reject = 0;
    contig->length = get_seq_length(align_entry->seq_entry);

    align_entry->contig = contig;
    align_entry->next = 0;
    align_entry->reverse = 0;
    align_entry->start = 0;
    align_entry->end = get_seq_length(align_entry->seq_entry) - 1;

    align_entry->score = align_entry->LLR_score = 0;
    align_entry->diffs = 0;
    align_entry->m_start = 0;
    align_entry->m_end = align_entry->m_length = align_entry->end;
    align_entry->bypassed = 0;
  }
}

compare_scores(pair1, pair2)
     Aligned_pair **pair1, **pair2;
{
  int d;

  if (d = (*pair2)->LLR_score - (*pair1)->LLR_score) return d;
  return (*pair2)->score - (*pair1)->score;
/*
  if (d = (*pair2)->score - (*pair1)->score) return d;
  return (*pair2)->score - (*pair1)->score;
*/
}

#define MAX_GAP_INDEX 1000
#define MAX_HIST_OFFSET 1000

Aligned_pair **LLR_pair_pointers;
int gap_hist1[MAX_GAP_INDEX], gap_hist2[MAX_GAP_INDEX];
int reject_hist1[MAX_GAP_INDEX], reject_hist2[MAX_GAP_INDEX];
int offset_hist[MAX_HIST_OFFSET];

old_merge_master()
{
  Aligned_pair *reject_pair;
  int i, i_ptr, pass;
  Aligned_pair *pair;
  Aligned_pair **sort_pairs();
  int LLR_join_cutoff, LLR_reject_cutoff, gap_cutoff;
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  Contig *contig1, *contig2;
  int offset, gap, discrep;
  int left_flag, right_flag, first_file_only, chim_flag;
  int w_offset;
  int entry1, entry2;

  notify("Merging ...");
  mark_save_block(); /* so segments can be freed later */
  alloc_contigs();
  init_contigs(1);
  LLR_pair_pointers = sort_pairs(compare_scores, 0);
  first_file_only = parameters->subject_files != 0;
  
  for (pass = 1; pass <= (first_file_only ? 1 : 4); pass++) {

    LLR_join_cutoff = (pass - 1) * -20;
    LLR_reject_cutoff = (pass - 1) * -20;
    gap_cutoff = pass * parameters->maxgap;

    for (i = 0; i < MAX_GAP_INDEX; i++) gap_hist1[i] = gap_hist2[i] = reject_hist1[i] = reject_hist2[i] = 0;
    for (i = 0; i < MAX_HIST_OFFSET; i++) offset_hist[i] = 0;
    fprintf(fp_log,"\n\nPass %d\n             Gap  Score    ", pass);

    for (i_ptr = 0; 
	 i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	 i_ptr++) {
   /* if (!is_best(pair)) continue;  only use original (weak) alignment for two entries */
      chim_flag = is_reject_chimeric(pair);
      if (is_used(pair) || is_repeat(pair) || pass < 4 && chim_flag) continue; 
/*      if (is_triple_reject(pair)) continue; */
      if (first_file_only) pair = pair->reversed_pair; 
      align_entry1 = get_align_entry(pair->entry1);
      align_entry2 = get_align_entry(pair->entry2);
      if (pass < 3 && (align_entry1->segments && align_entry1->segments->next 
		       || align_entry2->segments && align_entry2->segments->next))
	continue;
      if (first_file_only && (/* pair->entry2 >= num_query_entries
			      || */ align_entry2->contig->num_entries > 1)) continue;

      if (!first_file_only && pass < 4 && (align_entry1->first_start > align_entry1->last_end 
		       || align_entry2->first_start > align_entry2->last_end)) continue;
      /* keep out low-quality data during early passes -- since it tends to block ends of contigs */
      /* ignore other possible chimeras on first pass */
/*
      if (align_entry1->segments && align_entry1->segments->next
		    || align_entry2->segments && align_entry2->segments->next)
	continue; 
*/
      pair_merge(pair, chim_flag ? parameters->maxgap : gap_cutoff, 
		 LLR_reject_cutoff  - pair->LLR_score / 4 , LLR_join_cutoff, 1, 0, 0);
    }
    
    fprintf(fp_log, "\n\nN.B. Following not based on all pairs!!");
    fprintf(fp_log, "\n\nLowest   # merges  # failures");
    fprintf(fp_log, "\nLLR score");
    for (i = 0; i < MAX_GAP_INDEX; i++)
      if (reject_hist1[i] || reject_hist2[i])
	fprintf(fp_log, "\n%.1f    %4d     %4d", -i / 10.0, reject_hist1[i], reject_hist2[i]);
    fprintf(fp_log, "\n\nGap   # merges  # failures");
    fprintf(fp_log, "\n size");
    for (i = 0; i < MAX_GAP_INDEX; i++)
      if (gap_hist1[i] || gap_hist2[i])
	fprintf(fp_log, "\n%2d    %4d     %4d", i, gap_hist1[i], gap_hist2[i]);
    fprintf(fp_log, "\n\nOffset");
    for (i = 0; i < MAX_HIST_OFFSET; i++)
      if (offset_hist[i])
	fprintf(fp_log, "\n%3d    %4d ", i, offset_hist[i]);

    free_seg_blocks();
  } 

  if (parameters->revise_greedy && !first_file_only) {
    for (i_ptr = 0; i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > 10; 
	 i_ptr++) {
      if (is_used(pair)) continue;
      if (!is_best(pair)) continue;  
      if (is_reject_self(pair) || is_reject_node(pair) /*  || is_reject_vector(pair) */ ) continue;
      if (pair->entry1 == pair->entry2) continue;
      
      align_entry1 = get_align_entry(pair->entry1);
      align_entry2 = get_align_entry(pair->entry2);
      if (is_anomalous(align_entry1) || is_anomalous(align_entry2)) continue;
      
      contig1 = align_entry1->contig;
      contig2 = align_entry2->contig;

      if (contig1 == contig2) continue; /* NEED TO DEAL WITH THIS CASE 
					   (drop this condition eventually */
      
      if (is_reverse(pair) != (align_entry1->reverse != align_entry2->reverse)) {
	if (contig1 == contig2) continue; /* NEED TO DEAL WITH THIS CASE (inverted repeat) */
	complement(contig2); /* to ensure contigs are compatible */
      }
      
      offset = find_pair_offset(pair);
      
      test_merge(contig1, contig2, offset, gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff, 
		 &gap, &discrep, &left_flag, &right_flag, &w_offset, &reject_pair, &entry1, &entry2);
      
 /*     if (left_flag || right_flag) */
	split_contigs(contig1, contig2, offset, gap_cutoff, -60, 
		      10 /* LLR_join_cutoff */, &gap, &discrep, &left_flag, &right_flag);
    }
    /*    
      for (entry1 = 0; entry1 < t_num_entries; entry1++) {
      align_entry1 = align_array + entry1;
      if (unused_pos_match(align_entry1, 10)) 
      for (pair = align_entry1->first_pair; pair; pair = pair->next) {
      set_used_flag(pair, 0);
      set_repeat_flag(pair, 0);
      }
      else {
      for (pair = align_entry1->first_pair; pair; pair = pair->next) {
      set_repeat_flag(pair, 0);
      }
      }
      }
      for (entry1 = 0; entry1 < t_num_entries; entry1++) {
      align_entry1 = align_array + entry1;
      for (pair = align_entry1->first_pair; pair; pair = pair->next) 
      if (is_used(pair)) set_used_flag(pair->reversed_pair, 1);
      }
      */
    init_contigs(1);
    
    for (pass = 1; pass <= 2; pass++) {
      fprintf(fp_log,"\n\nReconstruction2 Pass %d\n             Gap  Score    ", pass);
      for (i_ptr = 0; i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_used(pair)) 
	  pair_merge(pair, gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff, 1, 0, 0);
      }
    }
    
    LLR_join_cutoff = 10; 
    LLR_reject_cutoff = -500;

    find_best_paths(gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff);
    
    init_contigs(1);
    
    LLR_join_cutoff = -50;
    LLR_reject_cutoff = -500; /* 750 */
    gap_cutoff = 100; /* 200 */

    for (pass = 1; pass <= 2; pass++) {
      fprintf(fp_log,"\n\nReconstruction3 Pass %d\n             Gap  Score    ", pass);
      for (i_ptr = 0; 
	   i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_used(pair)) pair_merge(pair, gap_cutoff,  LLR_reject_cutoff, LLR_join_cutoff, 1, 0, 0);
      }
    }
    for (i_ptr = 0; i_ptr < num_pairs; i_ptr++) {
      pair = LLR_pair_pointers[i_ptr];
      if (is_used(pair)) {
/*
	if (align_array[pair->entry1].contig != align_array[pair->entry2].contig) 
	  fprintf(stderr, "\nFAILED MERGE");
*/
	set_used_flag(pair, 0);
	set_used_flag(pair->reversed_pair, 0);
      }
      set_repeat_flag(pair, 0);
    }

    fprintf(fp_log,"\n\nReconstruction4    ");
    for (pass = 1; pass <= 2; pass++)
      for (i_ptr = 0; i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_used(pair) || is_repeat(pair) || pass == 1 && is_reject_chimeric(pair)) continue; 

	pair_merge(pair, gap_cutoff,  LLR_reject_cutoff, LLR_join_cutoff, 1, 0, 0);
      }
  }
  notify(" Done\n");
}
  
free_LLR_pair_pointers()
{
  our_free(LLR_pair_pointers);
}

static int num_repeats, num_used, num_split;
static int LLR_join_cutoff, LLR_reject_cutoff, gap_cutoff, reject_cutoff;

merge_master()
{
  Aligned_pair *reject_pair;
  int i, d, i_ptr, pass;
  Aligned_pair *pair;
  Aligned_pair **sort_pairs();
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  Contig *contig1, *contig2;
  int offset, gap, discrep;
  int left_flag, right_flag, first_file_only, weak_chim_flag, chim_flag;
  int w_offset;
  int entry1, entry2;

  notify("Merging ...");
  mark_save_block(); /* so segments can be freed later */
  alloc_contigs();
  init_contigs(1);
  LLR_pair_pointers = sort_pairs(compare_scores, 0);
  first_file_only = parameters->subject_files != 0;
  
  for (pass = 0; pass <= (first_file_only ? 1 : 5); pass++) {

    if (!pass && !parameters->preassemble) continue;
    if (pass == 2 /* && !parameters->revise_greedy */) continue;

    if (pass < 3) {
      LLR_join_cutoff = LLR_reject_cutoff = 0;
      gap_cutoff = parameters->maxgap;
    }
    else if (pass < 5) {
      LLR_join_cutoff = (pass - 1) * -20;
      LLR_reject_cutoff = (pass - 1) * -20; /* -150; */
      gap_cutoff = pass * parameters->maxgap;
    } 
    else {
      find_read_clusters();
      if (!parameters->forcelevel && !parameters->revise_greedy && !parameters->force_high) 
	break;
      /* leave LLR_join_cutoff the same; instead make 0? */
      LLR_reject_cutoff += ((FORCE_REJECT_SCORE + 1 - LLR_reject_cutoff) * parameters->forcelevel) / 10;
      gap_cutoff += ((1000 - gap_cutoff)  * parameters->forcelevel) / 10;
    }

/* force (essentially) all positive LLR scores in final pass -- probably need to
 be able to turn this off */

    for (i = 0; i < MAX_GAP_INDEX; i++) gap_hist1[i] = gap_hist2[i] = reject_hist1[i] = reject_hist2[i] = 0;
    for (i = 0; i < MAX_HIST_OFFSET; i++) offset_hist[i] = 0;
    fprintf(fp_log,"\n\nPass %d\n             Gap  Score    ", pass);

    for (i_ptr = 0; 
	 i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	 i_ptr++) {
   /* if (!is_best(pair)) continue;  only use original (weak) alignment for two entries */
      if (is_used(pair)) continue; 
      
/*      if (is_triple_reject(pair)) continue; */
      if (first_file_only) pair = pair->reversed_pair; 
      align_entry1 = get_align_entry(pair->entry1);
      align_entry2 = get_align_entry(pair->entry2);

      if (!pass && !same_group(pair->entry1, pair->entry2)) continue;
      if (pass < 3 && (is_blocked(align_entry1) || is_blocked(align_entry2))) continue;
/*
      if (pass < 4 && (is_lowqual(align_entry1) || is_lowqual(align_entry2))) continue;
*/
      if (pass < 4 && is_reject_chimeric(pair)) continue; 
      if (pass < 4 && is_reject_qual(pair)) continue;
      if (pass < 4 && is_repeat(pair) && 
	  align_entry1->contig == align_entry2->contig) 
	continue;
/* ignore other possible chimeras on first pass */

      if (pass < 3 && (align_entry1->segments && align_entry1->segments->next 
	  || align_entry2->segments && align_entry2->segments->next)) continue;

      if (first_file_only && (/* pair->entry2 >= num_query_entries
			      || */ align_entry2->contig->num_entries > 1)) 
	continue;

      if (pass < 4 && !first_file_only && (align_entry1->first_start > align_entry1->last_end 
		       || align_entry2->first_start > align_entry2->last_end)) 
	continue;
      /* keep out low-quality data during early passes -- since it tends to block ends of contigs */
      reject_cutoff = LLR_reject_cutoff  - pair->LLR_score / 4;
      if (reject_cutoff <= FORCE_REJECT_SCORE) reject_cutoff = FORCE_REJECT_SCORE + 1;

      chim_flag = align_entry1->segments && align_entry1->segments->next &&
	align_entry1->contig->num_entries > 1
	  || align_entry2->segments && align_entry2->segments->next &&
	align_entry2->contig->num_entries > 1;

      pair_merge(pair, chim_flag ? parameters->maxgap : gap_cutoff, 
		 reject_cutoff, LLR_join_cutoff, 1, 0, 
		 (pass >= 3) + 2 * (parameters->force_high && pass >= 5)); 

    }
    count_contigs(pass);

    if (pass != 2) {
      fprintf(fp_log, "\n\nN.B. Following not based on all pairs!!");
      fprintf(fp_log, "\n\nLowest   # merges  # failures");
      fprintf(fp_log, "\nLLR score");
      for (i = 0; i < MAX_GAP_INDEX; i++)
	if (reject_hist1[i] || reject_hist2[i])
	  fprintf(fp_log, "\n%.1f    %4d     %4d", 
		-i / 10.0, reject_hist1[i], reject_hist2[i]);
      fprintf(fp_log, "\n\nGap   # merges  # failures");
      fprintf(fp_log, "\n size");
      for (i = 0; i < MAX_GAP_INDEX; i++)
	if (gap_hist1[i] || gap_hist2[i])
	  fprintf(fp_log, "\n%2d    %4d     %4d", i, gap_hist1[i], gap_hist2[i]);
      fprintf(fp_log, "\n\nOffset");
      for (i = 0; i < MAX_HIST_OFFSET; i++)
	if (offset_hist[i])
	  fprintf(fp_log, "\n%3d    %4d ", i, offset_hist[i]);
    }
    else {
 /*     reassemble contigs */

      for (i_ptr = 0; 
	   i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	chim_flag = is_reject_chimeric(pair);
	if (is_repeat(pair)) num_repeats++;
	if (is_used(pair)) num_used++;

	if (is_used(pair) || is_repeat(pair) || chim_flag) continue; 
	align_entry1 = get_align_entry(pair->entry1);
	align_entry2 = get_align_entry(pair->entry2);
	if (align_entry1->segments && align_entry1->segments->next 
	    || align_entry2->segments && align_entry2->segments->next)
	  continue;
	if (align_entry1->first_start > align_entry1->last_end 
	    || align_entry2->first_start > align_entry2->last_end) 
	  continue;
/* now pair_merge tests */

	if (!is_best(pair)) continue;  
	if (is_reject_self(pair) || is_reject_node(pair))  continue;
	if (pair->entry1 == pair->entry2)  continue;
	if (is_anomalous(align_entry1) || is_anomalous(align_entry2)) continue;

	  /* what is left are unused, non-repeat-rejected pairs -- split
	     contigs at corresponding reads. Note that there will be redundant
	     calculations here! -- since used_pairs further down list are being reset to 0.*/
	num_split++;
	set_split_flag(pair, 1);
      }
      for (i_ptr = 0; 
	   i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_split(pair)) {
	  split_contigs(pair);
	  split_contigs(pair->reversed_pair);
	}
      }
      fprintf(fp_log, "\n\nNum repeat pairs: %d, num_used (orig): %d, num_split: %d", 
	      num_repeats, num_used, num_split);
      print_num_split();
      init_contigs(1);
      fprintf(fp_log,"\n\nReconstruction1 Pass\n             Gap  Score    ");
      for (i_ptr = 0; 
	   i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_used(pair)) {
	  reject_cutoff = LLR_reject_cutoff  - pair->LLR_score / 4;
	  if (reject_cutoff <= FORCE_REJECT_SCORE) reject_cutoff = FORCE_REJECT_SCORE + 1;
	  if (!pair_merge(pair, gap_cutoff, reject_cutoff, LLR_join_cutoff, 1, 0, 0)) {
	    notify("\nWARNING: Failed obligate remerge 1");
	    fprintf(fp_log,"\nWARNING: Failed obligate remerge 1");
	  }
	}
      }
      count_contigs(pass);
 
      find_best_paths(gap_cutoff, LLR_reject_cutoff, 10);
      
      init_contigs(1);
      fprintf(fp_log,"\n\nReconstruction2 Pass\n             Gap  Score    ");
      for (i_ptr = 0; 
	   i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_used(pair)) { /* missing brace pointed out by G. Montry 1/15/98
				may account for poor revise_greedy performance ? */
	  reject_cutoff = LLR_reject_cutoff  - pair->LLR_score / 4;
	  if (reject_cutoff < FORCE_REJECT_SCORE) reject_cutoff = FORCE_REJECT_SCORE + 1;
	  if (!pair_merge(pair, gap_cutoff, reject_cutoff, LLR_join_cutoff, 1, 1, 0)) {
	    notify("\nWARNING: Forced obligate remerge 2");
	    fprintf(fp_log,"\nWARNING: Forced obligate remerge 2");
	    set_used_flag(pair, 0);
	    set_used_flag(pair->reversed_pair, 0);
	  }
	}
      }
      count_contigs(pass);
    }
    free_seg_blocks();
  } 

/* following now obsolete */
  if (0 && parameters->revise_greedy && !first_file_only) {
    for (i_ptr = 0; i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > 10; 
	 i_ptr++) {
      if (is_used(pair)) continue;
      if (!is_best(pair)) continue;  
      if (is_reject_self(pair) || is_reject_node(pair)  || is_reject_vector(pair) ) continue;
      if (pair->entry1 == pair->entry2) continue;
      
      align_entry1 = get_align_entry(pair->entry1);
      align_entry2 = get_align_entry(pair->entry2);
      if (is_anomalous(align_entry1) || is_anomalous(align_entry2)) continue;
      
      contig1 = align_entry1->contig;
      contig2 = align_entry2->contig;

      if (contig1 == contig2) continue; /* NEED TO DEAL WITH THIS CASE 
					   (drop this condition eventually */
      
      if (is_reverse(pair) != (align_entry1->reverse != align_entry2->reverse)) {
	if (contig1 == contig2) continue; /* NEED TO DEAL WITH THIS CASE (inverted repeat) */
	complement(contig2); /* to ensure contigs are compatible */
      }
      
      offset = find_pair_offset(pair);
      
      test_merge(contig1, contig2, offset, gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff, 
		 &gap, &discrep, &left_flag, &right_flag, &w_offset, &reject_pair, &entry1, &entry2);
      
 /*     if (left_flag || right_flag) */
	split_contigs(contig1, contig2, offset, gap_cutoff, -60, 
		      10 /* LLR_join_cutoff */, &gap, &discrep, &left_flag, &right_flag);
    }
    /*    
      for (entry1 = 0; entry1 < t_num_entries; entry1++) {
      align_entry1 = align_array + entry1;
      if (unused_pos_match(align_entry1, 10)) 
      for (pair = align_entry1->first_pair; pair; pair = pair->next) {
      set_used_flag(pair, 0);
      set_repeat_flag(pair, 0);
      }
      else {
      for (pair = align_entry1->first_pair; pair; pair = pair->next) {
      set_repeat_flag(pair, 0);
      }
      }
      }
      for (entry1 = 0; entry1 < t_num_entries; entry1++) {
      align_entry1 = align_array + entry1;
      for (pair = align_entry1->first_pair; pair; pair = pair->next) 
      if (is_used(pair)) set_used_flag(pair->reversed_pair, 1);
      }
      */
    init_contigs(1);
    
    for (pass = 1; pass <= 2; pass++) {
      fprintf(fp_log,"\n\nReconstruction2 Pass %d\n             Gap  Score    ", pass);
      for (i_ptr = 0; i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_used(pair)) 
	  pair_merge(pair, gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff, 1, 0, 0);
      }
    }
    
    LLR_join_cutoff = 10; 
    LLR_reject_cutoff = -500;

    find_best_paths(gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff);
    
    init_contigs(1);
    
    LLR_join_cutoff = -50;
    LLR_reject_cutoff = -500; /* 750 */
    gap_cutoff = 100; /* 200 */

    for (pass = 1; pass <= 2; pass++) {
      fprintf(fp_log,"\n\nReconstruction3 Pass %d\n             Gap  Score    ", pass);
      for (i_ptr = 0; 
	   i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_used(pair)) pair_merge(pair, gap_cutoff,  LLR_reject_cutoff, LLR_join_cutoff, 1, 0, 0);
      }
    }
    for (i_ptr = 0; i_ptr < num_pairs; i_ptr++) {
      pair = LLR_pair_pointers[i_ptr];
      if (is_used(pair)) {
/*
	if (align_array[pair->entry1].contig != align_array[pair->entry2].contig) 
	  fprintf(stderr, "\nFAILED MERGE");
*/
	set_used_flag(pair, 0);
	set_used_flag(pair->reversed_pair, 0);
      }
      set_repeat_flag(pair, 0);
    }

    fprintf(fp_log,"\n\nReconstruction4    ");
    for (pass = 1; pass <= 2; pass++)
      for (i_ptr = 0; i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
	   i_ptr++) {
	if (is_used(pair) || is_repeat(pair) || pass == 1 && is_reject_chimeric(pair)) continue; 


	pair_merge(pair, gap_cutoff,  LLR_reject_cutoff, LLR_join_cutoff, 1, 0, 0);
      }
  }
  notify(" Done\n");
}

find_read_clusters()
{
  int i, j, k, m, length1, i_ptr;
  char *our_alloc();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Align_info *get_align_entry();
  Align_info *align_entry, *align_entry2;
  int *entry_classes, *entry_counts;
  int hist[10001];
  int pass;

  pass = 4;
  entry_classes = (int *)our_alloc(t_num_entries * sizeof(int));
  entry_counts = (int *)our_alloc(t_num_entries * sizeof(int));

  for (i = 0; i <= 10000; i++) hist[i] = 0;
  for (i = 0; i < t_num_entries; i++) {
    entry_classes[i] = i;
    entry_counts[i] = 0;
  }

  for (i = 0; i < t_num_entries; i++) {
    for (j = i; j != entry_classes[j]; j = entry_classes[j]);
    k = j;
    for (j = i; k != entry_classes[j];
	 m = j, j = entry_classes[j], entry_classes[m] = k);

    align_entry = get_align_entry(i);
    if (align_entry->contig->num_entries == 1) continue;
    length1 = get_seq_length(i);
    for (pair = get_aligned_pairs(i); pair; pair = pair->next) {
/* need better test here -- for chimera situation */
      if (!is_used(pair) && pair->LLR_score > 0 && !pair_merge_reject(pair)
	  && !is_reject_chimeric(pair)
	  && (align_entry2 = get_align_entry(pair->entry2))->contig->num_entries > 1) {
	if (align_entry2->contig != align_entry->contig
	    || is_reverse(pair) != (align_entry->reverse != align_entry2->reverse)  
	    || abs(find_pair_offset(pair)) > LOCATION_FUDGE1) {
	
	  break;
	}
      }
    }
    if (pair) continue;
    for (pair = get_aligned_pairs(i); pair; pair = pair->next) 
      if (is_used(pair) && pair->LLR_score > 0) {
	for (j = pair->entry2; k != entry_classes[j]; 
	     m = j, j = entry_classes[j], entry_classes[m] = k);
      }
  }

  for (i = 0; i < t_num_entries; i++) {
    for (j = i; j != entry_classes[j]; j = entry_classes[j]);
    entry_classes[i] = j;
    entry_counts[j] += 1;
  }    

  for (i = 0; i < t_num_entries; i++) {
    if (entry_counts[i]) 
      hist[entry_counts[i] < 10000 ? entry_counts[i] : 10000] += 1;
  }

  fprintf(stderr, "\n\nRead equivalence class histogram:");
  fprintf(fp_log, "\n\nRead equivalence class histogram:");

  for (i = 0; i <= 10000; i++)
    if (hist[i]) {
      fprintf(stderr, "\n%5d  %3d", i, hist[i]);
      fprintf(fp_log, "\n%5d  %3d", i, hist[i]);
    }
  fprintf(stderr, "\n");
  fprintf(fp_log, "\n");

  if (!parameters->revise_greedy && !parameters->shatter_greedy) 
    goto freeup;
/* Now -- make unused all pairs between different classes, and remerge;
   then do search of potential merges */

  for (i_ptr = 0; 
       i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
       i_ptr++) {
    if (is_used(pair) && 
	entry_classes[pair->entry1] != entry_classes[pair->entry2]
	&& (entry_counts[pair->entry1] != 1 || entry_counts[pair->entry2] != 1)) {
      set_used_flag(pair, 0);
      set_used_flag(pair->reversed_pair, 0);
    }
  }
  init_contigs(1);

  fprintf(fp_log,"\n\nReconstruction1 Pass\n             Gap  Score    ");
  for (i_ptr = 0; 
       i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
       i_ptr++) {
    if (is_used(pair)) {
      reject_cutoff = LLR_reject_cutoff  - pair->LLR_score / 4;
      if (reject_cutoff <= FORCE_REJECT_SCORE) reject_cutoff = FORCE_REJECT_SCORE + 1;
      if (!pair_merge(pair, gap_cutoff, reject_cutoff, LLR_join_cutoff, 1, 0, 1)) {
	notify("\nWARNING: Failed obligate remerge 1");
	fprintf(fp_log,"\nWARNING: Failed obligate remerge 1");
	set_used_flag(pair, 0);
	set_used_flag(pair->reversed_pair, 0);
      }
    }
  }
  count_contigs(pass);

  if (parameters->shatter_greedy) goto freeup;

/* remainder is executed only if revise_greedy option is selected */
  find_best_paths(gap_cutoff, LLR_reject_cutoff, 10);
      
  init_contigs(1);

  fprintf(fp_log,"\n\nReconstruction2 Pass\n             Gap  Score    ");
  for (i_ptr = 0; 
       i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
       i_ptr++) {
    if (is_used(pair)) {	
      reject_cutoff = LLR_reject_cutoff - pair->LLR_score / 4;
      if (reject_cutoff < FORCE_REJECT_SCORE) reject_cutoff = FORCE_REJECT_SCORE + 1;
      if (!pair_merge(pair, gap_cutoff, reject_cutoff, LLR_join_cutoff, 1, 0, 1)) {
	notify("\nWARNING: Failed obligate remerge 2");
	fprintf(fp_log,"\nWARNING: Failed obligate remerge 2");
	set_used_flag(pair, 0);
	set_used_flag(pair->reversed_pair, 0);
      }
    }
  }
  free_seg_blocks();
  count_contigs(pass);

 freeup:
  our_free(entry_classes);
  our_free(entry_counts);
}



merge_chimeras()
{
  int entry1;
  Align_info *align_entry;
  Align_info *get_align_entry();
  Aligned_pair *get_aligned_pairs();
  Aligned_pair *pair, *max_pair;
  char *o_qual1, *a_qual1;
  char *get_adj_qual(), *get_orig_qual();
  Segment *segments1, *segments2;
  int score, max_score;

  notify("Merging chimeras ...");
  fprintf(fp_log,"\nChimera merges: \n");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    align_entry = get_align_entry(entry1);
    if (align_entry->contig->num_entries > 1
	|| !align_entry->segments || !align_entry->segments->next) continue;
    o_qual1 = get_orig_qual(entry1);
    a_qual1 = get_adj_qual(entry1);
    segments1 = align_entry->segments;
/* recompute LLR scores, ignoring discrepancies outside alignment */
    max_score = 0;
    max_pair = 0;
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      segments2 = get_align_entry(pair->entry2)->segments;
      score = pair->LLR_score = pair->reversed_pair->LLR_score = 
	parameters->subject_files ? pair->score :
	  get_LLR(pair->diffs, o_qual1, get_orig_qual(pair->entry2), 
		  a_qual1, get_adj_qual(pair->entry2), 
		  get_seq_length(entry1), pair->start1, pair->end1, 
		  get_seq_length(pair->entry2), pair->start2, pair->end2, 
		  is_reverse(pair), segments1, segments2, 0, (FILE *)0, 0, 0, 0, 0, 1);
      if (score > max_score) {
	max_score = score;
	max_pair = pair;
      }
    }
    if (max_pair) {
      pair_merge(max_pair, 1000, -20, 0, 1, 0, 0);
      for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
	if (pair->LLR_score > 0 && pair != max_pair)
	  pair_merge(pair, 1000, -20, 0, 1, 0, 0);
      }
    }
  }
  notify(" Done\n");
}

merge_other_singletons()
{
  Align_info *get_align_entry();
  Aligned_pair *pair;
  int i_ptr;

  notify("Merging other singletons ...");
  fprintf(stderr,"\nLLR_join_cutoff: %d\n", LLR_join_cutoff);
  for (i_ptr = 0; 
       i_ptr < num_pairs && (pair = LLR_pair_pointers[i_ptr])->LLR_score > LLR_join_cutoff; 
       i_ptr++) {
    if ((get_align_entry(pair->entry1))->contig->num_entries > 1
      && (get_align_entry(pair->entry2))->contig->num_entries > 1) continue;
    reject_cutoff = LLR_reject_cutoff - pair->LLR_score / 4;
    if (reject_cutoff < FORCE_REJECT_SCORE) 
      reject_cutoff = FORCE_REJECT_SCORE + 1;
    pair_merge(pair, gap_cutoff, reject_cutoff, LLR_join_cutoff, 1, 0, 2);
  }
  notify(" Done\n");
}

pair_merge_reject(pair)
  Aligned_pair *pair;
{
  Align_info *get_align_entry();

  return !is_best(pair) || is_reject_self(pair) || is_reject_node(pair)
    || is_reject_vector(pair) || pair->entry1 == pair->entry2
      || is_anomalous(get_align_entry(pair->entry1)) 
      || is_anomalous(get_align_entry(pair->entry2));
}

int pair_merge(pair, gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff, print_flag, force_flag, allow_bad_entry)
  Aligned_pair *pair;
  int LLR_reject_cutoff, LLR_join_cutoff, gap_cutoff;
  int print_flag, force_flag, allow_bad_entry;
{
  int offset, gap, comp_flag;
  Contig *contig1, *contig2;
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  char *get_adj_qual(), *get_orig_qual();
  int score, success, discrep, t_print_flag, left_flag, right_flag, entries1, entries2;
  Merge_reject *find_merge_reject();
  Merge_reject *merge_reject;
  int reject_flag, reason;
  int w_offset;
  Aligned_pair *reject_pair;
  int r_entry1, r_entry2;
  int q1, q2, q3, q4, length2;
  int entry1, entry2, allow_flag;
  char *get_id();

  if (!parameters->bypasslevel) allow_bad_entry = 0; /* no bypasses allowed in this case */
  success = allow_flag = 0;
  if (pair_merge_reject(pair)) return 0;

  align_entry1 = get_align_entry(pair->entry1);
  align_entry2 = get_align_entry(pair->entry2);

  if (align_entry1->contig->num_entries < align_entry2->contig->num_entries) {
    pair = pair->reversed_pair;
    align_entry1 = get_align_entry(pair->entry1);
    align_entry2 = get_align_entry(pair->entry2);
  }
  /* ensure second contig has fewer entries, to reduce # tests */
  
  contig1 = align_entry1->contig;
  contig2 = align_entry2->contig;
  entries1 = contig1->num_entries;
  entries2 = contig2->num_entries;
  t_print_flag = 1;
    
  reject_flag = 0;

  if (contig1 == contig2) {
    if (is_reverse(pair) != (align_entry1->reverse != align_entry2->reverse)) {
/*
      if (print_flag) fprintf(fp_log,"\nIIR:            ");
*/
      offset = find_pair_offset(pair);
      /* ideally should also test pair->reversed_pair (here and in other
	 calls to test_read_merge) and set_repeat_flag should be invoked if
	 either meets condition (note that test_read_merge not symmetrical) */
      if (test_read_merge(contig1, align_entry2, offset, LLR_reject_cutoff, 1)) {
	set_repeat_flag(pair, 1);
	set_repeat_flag(pair->reversed_pair, 1);
      }
    }
    else {
      offset = find_pair_offset(pair);
      if (abs(offset) > LOCATION_FUDGE1) {
/*
	if (print_flag) fprintf(fp_log,"\nIDR:       %5d ", abs(offset)); 
*/
	offset = find_pair_offset(pair);
	if (test_read_merge(contig1, align_entry2, offset, LLR_reject_cutoff, 0)) {
	  set_repeat_flag(pair, 1);
	  set_repeat_flag(pair->reversed_pair, 1);
	}
      }
      /*
	else if (score = test_tandem(pair)) 
	fprintf(fp_log,"\nTAN:       %5d %3d >? ", abs(offset), score); 
	*/
      else {
	t_print_flag = 0; /* confirming matchups within contig
			   need not be printed */
	set_used_flag(pair, 1);
	set_used_flag(pair->reversed_pair, 1);
	contig1->num_matches += 1;
	success = 1;
      }
    }
  }
  else {
    /*	if (align_entry1->reverse) complement(contig1); ensure align_entry1 is
	in forward orientation */

    comp_flag = 0;
    if (is_reverse(pair) != (align_entry1->reverse != align_entry2->reverse)) { 
      complement(contig2); /* to ensure contigs are compatible */
      comp_flag = 1;
    }
    offset = find_pair_offset(pair);

    if (merge_reject = find_merge_reject(contig1, contig2, offset)) {
      reason = merge_reject->reject_reason;
      reject_flag = reason == 4
	|| merge_reject->lowest_LLR_score < LLR_reject_cutoff
	  || merge_reject->join_score <= LLR_join_cutoff 
	    && (merge_reject->gap > gap_cutoff || reason == 3); 
      if (!reject_flag) {
	delete_merge_reject(merge_reject->reverse);
	delete_merge_reject(merge_reject);
      }
/*
      else notify(".");
*/
    }
    w_offset = -1;
    if (!reject_flag && 
	(!(discrep = test_merge(contig1, contig2, offset, gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff, &gap, &score, &left_flag, &right_flag, &w_offset, &reject_pair, &entry1, &entry2)) 
	 || force_flag
	 || allow_bad_entry && (allow_flag = (entry1 > -2 || entry2 > -2) 
	    && (entries1 > 1 && entries2 > 1 || entries1 != entries2 && allow_bad_entry == 2)
	    && (score > FORCE_REJECT_SCORE || allow_bad_entry == 3)) && discrep == 2)) {
      /* what about similar upgrades to first_start, last_end when pair is within an already merged contig? */
      merge(contig1, contig2, offset);
      success = 1;
      gap_hist1[gap >= MAX_GAP_INDEX ? MAX_GAP_INDEX - 1 : gap] += 1;
      reject_hist1[-score >= MAX_GAP_INDEX ? MAX_GAP_INDEX - 1 : -score] += 1;
      set_used_flag(pair, 1);
      set_used_flag(pair->reversed_pair, 1);
    }
    
    else {
      if (!reject_flag) {
	append_merge_reject(contig1, contig2, offset, discrep,
			    pair->LLR_score, score, gap, pair, LLR_join_cutoff);
	
	reject_hist2[-score >= MAX_GAP_INDEX ? MAX_GAP_INDEX - 1 : -score] += 1;
	gap_hist2[gap >= MAX_GAP_INDEX ? MAX_GAP_INDEX - 1 : gap] += 1;
      }
      if (test_read_merge(contig1, align_entry2, offset, LLR_reject_cutoff, 0)) {
	set_repeat_flag(pair, 1);
	set_repeat_flag(pair->reversed_pair, 1);
      }
      if (comp_flag) complement(contig2);
    }
    if (!reject_flag && print_flag) {
      fprintf(fp_log,"\n%s%s(%d)  %3d  %5.1f", 
	      discrep && force_flag ? "FORCED " : "",
	      discrep && !(allow_bad_entry && allow_flag && discrep == 2) ? (discrep == 1 ? "Fail: " : "Rjct: ") : "Merge:", 
	      discrep, gap, score / 10.0);
    }
  } 
  if (!reject_flag && print_flag && t_print_flag && contig1 != contig2) {
    fprintf(fp_log," %.1f (%d,%d:%d) ", pair->LLR_score / 10.0, entries1, entries2, offset);
    print_pair(pair, fp_log, 1);
    if (allow_flag && discrep == 2 || (discrep == 2 || !discrep || discrep == 1 && gap < 500) && reject_pair && score < -30) {
      if (allow_flag) {
	fprintf(fp_log,"\n  Bypassed reads: contig1: %s; contig2: %s", 
		entry1 >= 0 ? get_id(entry1) : "none",
		entry2 >= 0 ? get_id(entry2) : "none");
	printf("\nBypassed reads: %s %s", 
		entry1 >= 0 ? get_id(entry1) : "",
		entry2 >= 0 ? get_id(entry2) : "");
	if (allow_bad_entry && entry1 >= 0 && discrep == 2) 
	  get_align_entry(entry1)->bypassed = 1;
	if (allow_bad_entry && entry2 >= 0 && discrep == 2) 
	  get_align_entry(entry2)->bypassed = 1;
      }
      fprintf(fp_log,"\n  Worst pair: offset %d ", w_offset);
      r_entry1 = reject_pair->entry1;
      r_entry2 = reject_pair->entry2;
      q1 = get_align_entry(r_entry1)->qual_start;
      q2 = get_align_entry(r_entry1)->qual_end;
      if (is_reverse(reject_pair)) {
	length2 = get_seq_length(r_entry2);
	q3 = length2 - 1 - get_align_entry(r_entry2)->qual_end;
	q4 = length2 - 1 - get_align_entry(r_entry2)->qual_start;
      }
      else {
	q3 = get_align_entry(r_entry2)->qual_start;
	q4 = get_align_entry(r_entry2)->qual_end;
      }
      get_LLR(reject_pair->diffs, get_orig_qual(r_entry1), get_orig_qual(r_entry2), 
	      get_adj_qual(r_entry1), get_adj_qual(r_entry2), get_seq_length(r_entry1),
	      reject_pair->start1, reject_pair->end1, get_seq_length(r_entry2), 
	      reject_pair->start2, reject_pair->end2, is_reverse(reject_pair), 
	      get_align_entry(r_entry1)->segments, get_align_entry(r_entry2)->segments, 
	      1, fp_log, q1, q2, q3, q4, 0);
      print_pair(reject_pair, fp_log, 1);
    }
  }
  return success;
}

delete_merge_reject(merge_reject)
     Merge_reject *merge_reject;
{
  if (merge_reject->next) merge_reject->next->prev = merge_reject->prev;
  if (merge_reject->prev) merge_reject->prev->next = merge_reject->next;
  else merge_reject->contig1->merge_reject = merge_reject->next;
  our_free(merge_reject);
}

Merge_reject *find_merge_reject(contig1, contig2, offset)
     Contig *contig1, *contig2;
     int offset;
{
  Merge_reject *merge_reject;

  for (merge_reject = contig1->merge_reject; merge_reject; 
       merge_reject = merge_reject->next) {
    if (contig2 == merge_reject->contig2 
	&& offset < merge_reject->offset + 5
	&& offset > merge_reject->offset - 5
	&& contig2->comp_status == merge_reject->comp_status2
	&& contig1->comp_status == merge_reject->comp_status1)
      return merge_reject;
  }
  return (Merge_reject *)0;
}

append_merge_reject(contig1, contig2, offset, reject_reason,
		    highest_LLR_score, lowest_LLR_score, gap, reject_pair, join_score)
     Contig *contig1, *contig2;
     int offset, reject_reason, highest_LLR_score, lowest_LLR_score, gap, join_score;
     Aligned_pair *reject_pair;
{
  char *our_alloc();
  Merge_reject *merge_reject1, *merge_reject2;

  merge_reject1 = (Merge_reject *)our_alloc(sizeof(Merge_reject));
  merge_reject2 = (Merge_reject *)our_alloc(sizeof(Merge_reject));
  merge_reject1->contig1 = merge_reject2->contig2 = contig1;
  merge_reject1->contig2 = merge_reject2->contig1 = contig2;
  merge_reject1->comp_status1 = merge_reject2->comp_status2 = contig1->comp_status;
  merge_reject1->comp_status2 = merge_reject2->comp_status1 = contig2->comp_status;
  merge_reject1->offset = merge_reject2->offset = offset;
  merge_reject1->reject_reason = merge_reject2->reject_reason = reject_reason;
  merge_reject1->highest_LLR_score = merge_reject2->highest_LLR_score = highest_LLR_score;
  merge_reject1->lowest_LLR_score = merge_reject2->lowest_LLR_score = lowest_LLR_score;
  merge_reject1->gap = merge_reject2->gap = gap;
  merge_reject1->reject_pair = merge_reject2->reject_pair = reject_pair;
  merge_reject1->join_score = merge_reject2->join_score = join_score;
 merge_reject1->next = contig1->merge_reject;
  merge_reject1->prev = 0;
  if (contig1->merge_reject) contig1->merge_reject->prev = merge_reject1;
  contig1->merge_reject = merge_reject1;
  merge_reject1->reverse = merge_reject2;

  merge_reject2->next = contig2->merge_reject;
  merge_reject2->prev = 0;
  if (contig2->merge_reject) contig2->merge_reject->prev = merge_reject2;
  contig2->merge_reject = merge_reject2;
  merge_reject2->reverse = merge_reject1;
}


/*
debug_contig(contig)
     Contig *contig;
{
  Align_info *align_entry;
  char *get_id();
  int flag1, flag2;

  flag1 = flag2 = 0;
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    if (!strcmp(align_entry->db_entry->id, "zb60g3.s1")
	|| !strcmp(align_entry->db_entry->id, "zb56c3.s1")
	|| !strcmp(align_entry->db_entry->id, "zc10c5.s1"))
      flag1 = 1;
    if (!strcmp(align_entry->db_entry->id, "zb42e10.s1")
	|| !strcmp(align_entry->db_entry->id, "zb44e8.s1")
	|| !strcmp(align_entry->db_entry->id, "zb45a7.s1")
	|| !strcmp(align_entry->db_entry->id, "zb52e8.f1c")
	|| !strcmp(align_entry->db_entry->id, "zb48f10.s1"))
      flag2 = 1;
  }
  if (flag1 && flag2) exit(1);
}
*/
/* find implied offset for two contigs implied by a pair relating them (assumes
   contigs are either equal, or have been properly oriented with respect 
   to each other */

int find_pair_offset(pair)
     Aligned_pair *pair;
{
  int offset1, offset2;
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();

  align_entry1 = get_align_entry(pair->entry1);
  align_entry2 = get_align_entry(pair->entry2);
  
/*
  offset = (pair_start1 + pair_end1) / 2.0 - (pair_start2 + pair_end2) / 2.0;
  return align_entry1->reverse ? align_entry1->end - align_entry2->end - offset
    : align_entry1->start - align_entry2->start + offset;
*/
  if (align_entry1->contig == align_entry2->contig 
      && is_reverse(pair) != (align_entry1->reverse != align_entry2->reverse)) {

    if (align_entry1->reverse) {
      offset1 = align_entry1->end - (-align_entry2->start) - pair->offset; /* pair->start1 + pair->start2; */
      offset2 = align_entry1->start + get_seq_length(pair->entry1)
	- (-align_entry2->end + get_seq_length(pair->entry2)) - pair->offset; /*- (pair->end1   pair->end2); */
    }
    else {
      offset1 = align_entry1->start - (-align_entry2->end) + pair->offset; /* pair->start1 - pair->start2; */
      offset2 = align_entry1->end - get_seq_length(pair->entry1)
	- (-align_entry2->start - get_seq_length(pair->entry2)) + pair->offset; /* pair->end1 - + pair->end2  ;  */
    }
  }
  else {
    if (align_entry1->reverse) {
      offset1 = align_entry1->end - align_entry2->end - pair->offset; /* pair->start1 + pair->start2; */
      offset2 = align_entry1->start + get_seq_length(pair->entry1)
	- (align_entry2->start + get_seq_length(pair->entry2)) - pair->offset; /*- (pair->end1   pair->end2); */
    }
    else {
      offset1 = align_entry1->start - align_entry2->start + pair->offset; /* pair->start1 - pair->start2; */
      offset2 = align_entry1->end - get_seq_length(pair->entry1)
	- (align_entry2->end - get_seq_length(pair->entry2)) + pair->offset; /* pair->end1 - + pair->end2  ;  */
    }
  }

  return abs(offset1) < abs(offset2) ? offset1 : offset2;
}

count_contigs(pass)
     int pass;
{
  int contig_hist[10001];
  Contig *contig;
  int i, n;
  Aligned_pair *get_aligned_pairs();
  Aligned_pair *pair;

  for (i = 0; i <= 10000; i++) contig_hist[i] = 0;

  for (i = 0; i < t_num_entries; i++) {
    contig = contig_array + i;
    if (n = contig->num_entries) {
      if (n == 1) { /* check no pairwise matches to other entries */
	for (pair = get_aligned_pairs(i); pair && i == pair->entry2; pair = pair->next);
	if (!pair) continue;
      }
      contig_hist[n > 10000 ? 10000 : n] += 1;
    }
  }

  fprintf(fp_log, "\n\nPass: %d\n#reads  #contigs (not counting singlets)", pass);
  fprintf(stderr, "\n\nPass: %d\n#reads  #contigs (not counting singlets)", pass);
  for (i = 0; i <= 10000; i++) 
    if (contig_hist[i]) {
      fprintf(stderr, "\n%5d  %5d", i, contig_hist[i]);
      fprintf(fp_log, "\n%5d  %5d", i, contig_hist[i]);
    }
}

/* check whether reads are from same subclone or not */

/* revise first_start, last_end of contig -- so that only determined from overlaps already
   incorporated into contig (in particular single unincorporated reads, e.g. undetected chimeras, 
   will basically have no restrictions -- is this a good idea??) */

reset_starts_ends()
{
  int i, start, end, start1, start2, end1, end2;
  Contig *contig;
  Align_info *align_entry;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();

  for (i = 0; i < t_num_entries; i++) {
    contig = contig_array + i;
    if (!contig->num_entries) continue;
    contig->last_end = 0;
    contig->first_start = get_seq_length(contig->first->seq_entry) - 1;
    for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
      for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
	if (!is_used(pair)) continue;
	start1 = pair->start1;
	end1 = pair->end1;
	start2 = pair->start2;
	end2 = pair->end2;
	start = align_entry->reverse ? align_entry->end - end1 
	  : align_entry->start + start1;
	end = align_entry->reverse ? align_entry->end - start1
	  : align_entry->start + end1;
	if (start < contig->first_start) contig->first_start = start;
	if (end > contig->last_end) contig->last_end = end;
      }
    }
  }
}

complement(contig)
     Contig *contig;
{
  Align_info *align_entry;
  int temp;

  temp = contig->first_start;
  contig->first_start = -contig->last_end;
  contig->last_end = -temp;
  contig->comp_status = !contig->comp_status;
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    temp = align_entry->start; 
    align_entry->start = -align_entry->end;
    align_entry->end = -temp;
    align_entry->reverse = !align_entry->reverse;
  }
}
     
merge(contig1, contig2, offset)
     Contig *contig1, *contig2;
     int offset; /* proposed offset of contig2 relative to contig1 */
{
  Align_info *align_entry;
  Merge_reject *merge_reject, *next;

  contig1->last->next = contig2->first;
  contig1->last = contig2->last;
  contig1->num_entries += contig2->num_entries;
  for (align_entry = contig2->first; align_entry; align_entry = align_entry->next) {
    align_entry->contig = contig1;
    align_entry->start += offset;
    align_entry->end += offset;
  }
  if (contig1->first_start > contig2->first_start + offset)
/*      && contig2->first_start < contig2->last_end) */
    contig1->first_start = contig2->first_start + offset;
  if (contig1->last_end < contig2->last_end + offset)
/*       && contig2->first_start < contig2->last_end)  */
    contig1->last_end = contig2->last_end + offset;
  
  contig2->first = contig2->last = 0;
  contig2->num_entries = 0;
  contig1->num_matches += 1 + contig2->num_matches;
  for (merge_reject = contig1->merge_reject; merge_reject; 
       merge_reject = next) {
    next = merge_reject->next;
    delete_merge_reject(merge_reject->reverse);
    delete_merge_reject(merge_reject);
  }
  for (merge_reject = contig2->merge_reject; merge_reject; 
       merge_reject = next) {
    next = merge_reject->next;
    delete_merge_reject(merge_reject->reverse);
    delete_merge_reject(merge_reject);
  }
}

/* for contigs in correct orientation, with specified offset, identifies clumps of pos. scoring matches
   in each contig and splits off those from rest of contig */

old_split_contigs(contig1, contig2, offset, gap_cutoff, reject_score, join_score, gap, min_score, left_flag, right_flag)
     Contig *contig1, *contig2;
     int offset; /* proposed offset of contig2 relative to contig1 */
     int gap_cutoff, reject_score, join_score; /* max gap size, LLR_scores for rejecting, accepting */
     int *gap, *min_score;
     int *left_flag, *right_flag; /* indicate that contig2 extends to left, to right */
{
  Align_info *align_entry, *align_entry2, *align_entry3;
  Align_info *get_align_entry();
  Aligned_pair *pair, *pair2, *pair3;
  Aligned_pair *get_aligned_pairs();

  for (align_entry = contig2->first; align_entry; align_entry = align_entry->next) 
    align_entry->temp = 0;

  for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) 
    align_entry->temp = 0;

  for (align_entry = contig2->first; align_entry; align_entry = align_entry->next) 
    for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
      if (!is_best(pair)) continue;  
      if (is_reject_self(pair)) continue;
      align_entry2 = get_align_entry(pair->entry2);
      if (align_entry2->contig != contig1
	  || is_reverse(pair) != (align_entry->reverse != align_entry2->reverse)) 
	continue;
      if (abs(offset + find_pair_offset(pair)) > LOCATION_FUDGE1) continue;
/*
      if (pair->LLR_score < reject_score || is_left_trunc(pair) || is_right_trunc(pair)) {
	align_entry->temp |= 2;
	align_entry2->temp |= 2;
      }
*/
      if (pair->LLR_score <= join_score) continue;
      for (pair2 = get_aligned_pairs(align_entry2->seq_entry); pair2; pair2 = pair2->next) {
	if (!is_used(pair2)) continue;
	for (pair3 = get_aligned_pairs(align_entry->seq_entry); pair3; pair3 = pair3->next) {
	  if (pair3->entry2 != pair2->entry2) continue;
	  if (pair3->LLR_score >= 0 /* reject_score */  && !is_left_trunc(pair3) && !is_right_trunc(pair3)) 
	    continue;
	  align_entry3 = get_align_entry(pair3->entry2);;
	  if (align_entry3->contig != contig1
	      || is_reverse(pair3) != (align_entry->reverse != align_entry3->reverse)) 
	    continue;
	  if (abs(offset + find_pair_offset(pair3)) > LOCATION_FUDGE1) continue;
	  set_used_flag(pair2, 0);
	  set_used_flag(pair2->reversed_pair, 0);
  
	}
      }
/*
      align_entry->temp |= 1;
      align_entry2->temp |= 1;
*/
    }

/*
  for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) { 
    if (align_entry->temp & 1) 
      for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) 
	if (is_used(pair) && (get_align_entry(pair->entry2)->temp & 2) || !get_align_entry(pair->entry2)->temp) {
	  set_used_flag(pair, 0);
	  set_used_flag(pair->reversed_pair, 0);
	}
  }
  for (align_entry = contig2->first; align_entry; align_entry = align_entry->next) { 
    if (align_entry->temp & 1) 
      for (pair = align_entry->first_pair; pair; pair = pair->next) 
	if (is_used(pair) && (align_array[pair->entry2].temp & 2 || !align_array[pair->entry2].temp)) {
	  set_used_flag(pair, 0);
	  set_used_flag(pair->reversed_pair, 0);
	}
  }
*/
}

/* N.B. This is very different from old version!!! */
static int num_split_secondary, num_split_tertiary;
static int num_split_primary[5];

split_contigs(pair)
     Aligned_pair *pair;
{
  int direction, direction2;
  Aligned_pair *pair1, *pair2;
  Aligned_pair *get_aligned_pairs();

  direction = get_direction(pair);
  num_split_primary[direction] += 1;

  for (pair1 = get_aligned_pairs(pair->entry1); pair1; 
       pair1 = pair1->next)
    if (is_used(pair1) && (direction & get_direction(pair1))) {
      set_used_flag(pair1, 0);
      set_used_flag(pair1->reversed_pair, 0);
      num_split_secondary++;
      direction2 = get_direction(pair1->reversed_pair);
      for (pair2 = get_aligned_pairs(pair1->entry2); pair2; 
	   pair2 = pair2->next) 
	if (is_used(pair2) && (direction2 & get_direction(pair2))) {
	  set_used_flag(pair2, 0);
	  set_used_flag(pair2->reversed_pair, 0);
	  num_split_tertiary++;
	}
    }
}

print_num_split()
{
  int i;

  fprintf(fp_log, "\n\nNo. of split pairs:\n");
  for (i = 0; i < 5; i++)
    if (num_split_primary[i])
      fprintf(fp_log, "\n%d  %d", i, num_split_primary[i]);
  fprintf(fp_log, "\n\nSecondary: %d, tertiary: %d\n", 
	  num_split_secondary, num_split_tertiary);
}

/* determines whether entry2 extends to left or right (or both, or neither)
   of entry1 */

int get_direction(pair)
     Aligned_pair *pair;
{
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  int start1, end1, start2, end2, length2;
  int d;

  align_entry1 = get_align_entry(pair->entry1);
  align_entry2 = get_align_entry(pair->entry2);

  start1 = align_entry1->first_start;
  end1 = align_entry1->last_end;

  if (is_reverse(pair)) {
    length2 = get_seq_length(pair->entry2) - 1;
    start2 = length2 - align_entry2->last_end;
    end2 = length2 - align_entry2->first_start;
  }
  else {
    start2 = align_entry2->first_start;
    end2 = align_entry2->last_end;
  }
  
  d = 0;
  if (start2 < start1 - pair->offset) d |= 1;
  if (end2 > end1 - pair->offset) d |= 2;

  return d;
}

/* following applies to ordered linked list of non-overlapping "segments" of sequence --
   i.e. segment->start < segment->end < segment->next->start
*/

/* in deciding to merge, always merge the shorter into the longer */

test_merge(contig1, contig2, offset, gap_cutoff, reject_score, join_score, gap, min_score, left_flag, right_flag, w_offset, reject_pair, add_entry1, add_entry2)
     Contig *contig1, *contig2;
     int offset; /* proposed offset of contig2 relative to contig1 */
     int gap_cutoff, reject_score, join_score; /* max gap size, LLR_scores for rejecting, accepting */
     int *gap, *min_score;
     int *left_flag, *right_flag; /* indicate that contig2 extends to left, to right */
     int *w_offset;
     Aligned_pair **reject_pair;
     int *add_entry1, *add_entry2;
{
  /* create segments */
  Segment *head1, *head2, *segment;
  Segment *insert_segment();
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int start, end, start1, end1, start2, end2;
  int pair_start1, pair_start2, pair_end1, pair_end2;
  int t_offset, pair_offset; 
  int temp, length2, gap1, gap2;
  char vector_reject;
  int entry1, entry2;
  char *get_id();

  head2 = head1 = 0;
  *min_score = 0;
  vector_reject = 0;
  *gap = 0;

  *reject_pair = 0;
  *add_entry1 = *add_entry2 = -1;
  for (align_entry = contig2->first; align_entry; align_entry = align_entry->next) {
    for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
/* N.B. This condition inappropriate for tandem repeats! Is it otherwise nec???
      if (!is_best(pair)) continue;  
*/

/*      if (is_triple_reject(pair)) continue; */
/*      if (pair->LLR_score < score) continue; */
/*      if (pair->score < score && 
	  (!is_reject_total(pair) || pair->score < parameters->confirm_score)) 
	continue; 2d condition -- to ensure no significant rejected alignments 
		    are implied 
*/
      align_entry2 = get_align_entry(pair->entry2);
      if (align_entry2->contig != contig1
	  || is_reverse(pair) != (align_entry->reverse != align_entry2->reverse)) 
	continue;
      pair_offset = -find_pair_offset(pair);
      t_offset = abs(offset - pair_offset); 

      if (t_offset > LOCATION_FUDGE2 && (pair->LLR_score != FORCE_REJECT_SCORE || t_offset > 50)) {
/*
	   if (pair->LLR_score > 10 && t_offset < 100) 
	  fprintf(stderr, "%d ", t_offset); 
	  fprintf(stderr, "\n%d %d  %d", offset, pair_offset, is_reverse(pair)); 
*/
	continue;
      }
      if (t_offset < 10 /* LOCATION_FUDGE2 */|| pair->LLR_score == FORCE_REJECT_SCORE) { /* MAKE MORE STRINGENT */
	if (pair->LLR_score < *min_score) {
	  *min_score = pair->LLR_score;
	  *reject_pair = pair->reversed_pair;
	}
	if (pair->LLR_score < reject_score) {
	  if (*add_entry1 == -1) *add_entry1 = pair->entry2;
	  else if (*add_entry1 != -2 && !same_subclone(*add_entry1, pair->entry2)) *add_entry1 = -2;
	  if (*add_entry2 == -1) *add_entry2 = pair->entry1;
	  else if (*add_entry2 != -2 && !same_subclone(*add_entry2, pair->entry1)) *add_entry2 = -2;
	}
/* 	if (*min_score < reject_score) return 2; N.B. THIS PREVENTS FULL CALCULATION !!! */
	*w_offset = t_offset;
      }
      if (is_reject_vector(pair)) vector_reject = 1;
      if (pair->LLR_score <= join_score) continue;
      if (pair->LLR_score > join_score) 
	offset_hist[t_offset >= MAX_HIST_OFFSET ? MAX_HIST_OFFSET - 1 : t_offset] += 1;
      pair_start1 = pair->start1;
      pair_start2 = pair->start2;
      pair_end1 = pair->end1;
      pair_end2 = pair->end2;
      
      if (is_reverse(pair)) {
	length2 = get_seq_length(pair->entry2) - 1;
	temp = pair_start2;
	pair_start2 = length2 - pair_end2;
	pair_end2 = length2 - temp;
      }
      head2 = align_entry->reverse ? insert_segment(head2, align_entry->end - pair_end1,
						    align_entry->end - pair_start1)
	: insert_segment(head2, align_entry->start + pair_start1, 
			 align_entry->start + pair_end1);
      
      head1 = align_entry2->reverse ? insert_segment(head1, align_entry2->end - pair_end2, 
						     align_entry2->end - pair_start2)
	: insert_segment(head1, align_entry2->start + pair_start2, 
			 align_entry2->start + pair_end2);
    }
  }
  /* define start and end here -- based on implied overlap */
  if (!head2 || !head1) {
    fprintf(stderr, "\nMISSING SEGMENT %d %d", offset, (int)(contig2->first->reverse));
    return 3;
  }

  start = contig2->first_start;
  if (start < contig1->first_start - offset) start = contig1->first_start - offset;
  end = contig2->last_end;
  if (end > contig1->last_end - offset) end = contig1->last_end - offset;
  gap2 = find_max_gap(head2, start, end);
/*  free_segments(head2);  */
  if (*min_score >= reject_score && gap2 > gap_cutoff) {
    fprintf(fp_log,"\n %d  %d  %d %d    %d   ", contig2->first_start, contig1->first_start - offset,
	    contig2->last_end, contig1->last_end - offset, gap2);
    print_segment_gaps(fp_log, head2, start, end);
  }

  start = contig1->first_start;
  if (start < contig2->first_start + offset) start = contig2->first_start + offset;
  end = contig1->last_end;
  if (end > contig2->last_end + offset) end = contig2->last_end + offset;
  gap1 = find_max_gap(head1, start, end);

  *gap = gap1 < gap2 ? gap1 : gap2;


   
  *left_flag = *right_flag = 0;

  start1 = head1->start;
  for (segment = head1; segment && segment->next; segment = segment->next);
  end1 = segment->end;

  start2 = head2->start;
  for (segment = head2; segment && segment->next; segment = segment->next);
  end2 = segment->end;

  if (start1 < contig1->first_start + gap_cutoff) /* match includes left end of contig */
    *left_flag = *gap <= gap_cutoff && end2 > contig2->last_end - gap_cutoff ? 2 : 1;

  if (end1 > contig1->last_end  - gap_cutoff) /* match includes right end of contig */
    *right_flag = *gap <= gap_cutoff && start2 < contig2->first_start + gap_cutoff ? 2 : 1;

/* TAKEN OUT 7/19/98
  if (vector_reject) return 4;
*/

  return *gap > gap_cutoff ? 1 : (*min_score < reject_score ? 2 : 0) ;

/*  free_segments(head1);  */

}

test_read_merge(contig1, align_entry, offset, reject_score, self_comp)
     Contig *contig1;
     Align_info *align_entry;
     int offset; /* proposed offset of contig2 relative to contig1 */
     int reject_score;
     int self_comp;
{
  Align_info *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int t_offset, pair_offset; 
  int n_pairs, min_score;
 
  n_pairs = 0;
  min_score = 100;
  for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
    align_entry2 = get_align_entry(pair->entry2);
    if (align_entry2->contig != contig1
	|| self_comp != (is_reverse(pair) != (align_entry->reverse != align_entry2->reverse))) 
      continue;
    pair_offset = self_comp ? find_pair_offset(pair) : -find_pair_offset(pair);
    t_offset = abs(offset - pair_offset); 
    if (t_offset > LOCATION_FUDGE2) continue;
    if (pair->LLR_score < min_score) {
      min_score = pair->LLR_score;
      if (min_score < reject_score) return 2;
    }
    n_pairs++;
  }
/*
  fprintf(stderr,"\n%d tested pairs, %d", n_pairs, min_score);
*/
  return 0;
}

compare_contigs(contig1, contig2)
     Contig **contig1, **contig2;
{
  int d;
/*
  if ((*contig1)->num_contigs != (*contig2)->num_contigs)
    return (*contig1)->num_contigs - (*contig2)->num_contigs;
  if ((*contig1)->t_num_entries != (*contig2)->t_num_entries)
   return (*contig1)->t_num_entries - (*contig2)->t_num_entries;
  if ((*contig1)->parent != (*contig2)->parent)
    return (*contig1)->parent - (*contig2)->parent;
*/
  if (d = is_duplicate((*contig2)->first) - is_duplicate((*contig1)->first))
    return d;
  /* put duplicates before other singlet contigs */

  if (d = (*contig1)->num_entries - (*contig2)->num_entries)
    return d; /* sort by no. of subclones within contig */
  if ((*contig1)->t_num_entries == 1 || (*contig2)->t_num_entries == 1)
    return (*contig1)->t_num_entries - (*contig2)->t_num_entries;
  /* put singletons before everything else */
  return (*contig2)->first->anomalies - (*contig1)->first->anomalies;
 /* sort anomalous contigs by type of anomaly */
}

compare_starts(entry1, entry2)
     Align_info **entry1, **entry2;
{
  return (*entry1)->start - (*entry2)->start;
}
  
static Contig **contig_ptrs;

clean_contigs()
{
  int i, j, base_offset, n_entries, index;
  Align_info *align_entry;
  Align_info *get_align_entry();
  char *our_alloc();
  Align_info **entry_ptrs;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Contig *contig1, *contig2;

  for (i = 0; i < t_num_entries; i++) {
    contig1 = contig_array + i;
    if (contig1->num_entries) {
      contig1->parent = contig1;
      contig1->num_contigs = contig1->t_num_entries = 0;
      /* ensure no more than half the reads are in reverse orientation 
	 -- for consistency */
      if (!parameters->subject_files) {
	j = 0;
	for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) 
	  if (align_entry->reverse) j++;
	if (j + j > contig1->num_entries) complement(contig1);
      }
    }
  }

/* can following be commented out */
  for (i = 0; i < t_num_entries; i++) {
    contig1 = contig_array + i;
    if (contig1->num_entries) {
      while (contig1->parent != contig1->parent->parent)
	contig1->parent = contig1->parent->parent;
      for (align_entry = contig1->first; align_entry;
	   align_entry = align_entry->next) {
	for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
	  if (is_reject_self(pair)) continue;
	  contig2 = get_align_entry(pair->entry2)->contig;
	  if (contig2 > contig1) { /* only merge in one direction */
	    while (contig2->parent != contig2->parent->parent)
	      contig2->parent = contig2->parent->parent;
	    if (contig2->parent != contig1->parent)
	      contig2->parent = contig2->parent->parent = contig1->parent;
	  }
	}
      }
    }
  }
  
  for (i = num_contigs = 0; i < t_num_entries; i++) {
    contig1 = contig_array + i;
    if (contig1->num_entries) {
      num_contigs++;
      while (contig1->parent != contig1->parent->parent)
	contig1->parent = contig1->parent->parent;
      contig1->parent->num_contigs += 1;
      contig1->parent->t_num_entries += contig1->num_entries;
    }
  }

  contig_ptrs = (Contig **)our_alloc(num_contigs * sizeof(Contig *));
  for (i = j = 0; i < t_num_entries; i++) {
    contig1 = contig_array + i;
    if (contig1->num_entries) {
      contig_ptrs[j++] = contig1;
      contig1->t_num_entries = contig1->parent->t_num_entries;
      contig1->num_contigs = contig1->parent->num_contigs;
    }
  }
  qsort(contig_ptrs, num_contigs, sizeof(Contig *), compare_contigs);
  
  index = 1;
  for (i = num_duplicates = 0; i < num_contigs && is_duplicate(contig_ptrs[i]->first); i++) 
    num_duplicates++;

  for (num_singletons = 0; i < num_contigs && contig_ptrs[i]->t_num_entries == 1; i++) 
    num_singletons++;

  for (; i < num_contigs; i++) {
    contig1 = contig_ptrs[i];
    if (n_entries = contig1->num_entries) {
      entry_ptrs = (Align_info **)our_alloc(n_entries * sizeof(Align_info *));
      for (j = 0, align_entry = contig1->first; align_entry; align_entry = align_entry->next, j++) {
	entry_ptrs[j] = align_entry;
      }
      qsort(entry_ptrs, n_entries, sizeof(Align_info *), compare_starts);
      base_offset = entry_ptrs[0]->start - 1;
      contig1->index = index++;
      contig1->first_start -= base_offset;
      contig1->last_end -= base_offset;
      contig1->first = entry_ptrs[0];
      contig1->last = entry_ptrs[n_entries - 1];
      contig1->length = 0;
      for (j = 0; j < n_entries; j++) {
	align_entry = entry_ptrs[j];
	align_entry->start -= base_offset;
	align_entry->end -= base_offset;
	if (contig1->length < align_entry->end) contig1->length = align_entry->end;
	align_entry->next = j < n_entries - 1 ? entry_ptrs[j+1] : 0;
      }
      our_free(entry_ptrs); 
    }
  }
}

analyze_slack(contig1, fp)
     Contig *contig1;
     FILE *fp;
{
  int i, j, offset1, offset2, score;
  int unused_hist1[100][2], used_hist1[100][2];
  int diff_offset_hist[100];
  char *get_id();
  Align_info *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int LLR_hist[2][50];
  int LLR_score;
  char *get_adj_qual(), *get_orig_qual();
  int r_entry1, r_entry2;
  int q1, q2, q3, q4, length2;

  for (j = 0; j < 2; j++)
    for (i = 0; i < 100; i++) 
      unused_hist1[i][j] = used_hist1[i][j] = 0;

  for (j = 0; j < 2; j++)
    for (i = 0; i < 50; i++) 
      LLR_hist[j][i] = 0;

  for (i = 0; i < 100; i++) diff_offset_hist[i] = 0;
  for (align_entry1 = contig1->first; align_entry1; align_entry1 = align_entry1->next) {
    for (pair = get_aligned_pairs(align_entry1->seq_entry); pair; pair = pair->next) {

      if (!is_best(pair)) continue;
      if (is_reject_self(pair)) continue;
      if (pair->entry1 >= pair->entry2) continue;  
      align_entry2 = get_align_entry(pair->entry2);
      if (is_anomalous(align_entry2)) continue;
      
      score = pair->score;
      LLR_score = pair->LLR_score / 10.0 + 10;
      if (LLR_score < 0) LLR_score = 0;
      if (LLR_score > 20) LLR_score = 20;
      LLR_hist[!is_used(pair)][LLR_score] += 1;
      if (contig1 == align_entry2->contig 
	  && is_reverse(pair) == (align_entry1->reverse != align_entry2->reverse)) {
	offset1 = abs(find_pair_offset(pair));
	offset2 = abs(find_pair_offset(pair->reversed_pair));
	offset2 = abs(offset1 - offset2);
	if (offset2 > 99) offset2 = 99;
	diff_offset_hist[offset2] += 1;

	if (offset1 > 99) offset1 = 99;
	if (is_used(pair)) {
	  if (!used_hist1[offset1][0] || pair->LLR_score > used_hist1[offset1][1]) 
	    used_hist1[offset1][1] = pair->LLR_score;
	  used_hist1[offset1][0] += 1;
	  /* rejected pairs may have been used if rejection is on basis
	     of alignment extent, rather than high-quality mismatches */
	}
	else {
	  if (!unused_hist1[offset1][0] || pair->LLR_score > unused_hist1[offset1][1]) 
	    unused_hist1[offset1][1] = pair->LLR_score;
	  unused_hist1[offset1][0] += 1;
	  if (offset1 < 20 && !is_reject_node(pair)) {
	    fprintf(fp, "\n%s Unused pair: %s %s  %.1f   %d",
		    is_reject_chimeric(pair) ? "Chimeric" : "",
		    get_id(pair->entry1), get_id(pair->entry2),
		    pair->LLR_score / 10.0, offset1);
	    r_entry1 = pair->entry1;
	    r_entry2 = pair->entry2;
	    q1 = get_align_entry(r_entry1)->qual_start;
	    q2 = get_align_entry(r_entry1)->qual_end;
	    if (is_reverse(pair)) {
	      length2 = get_seq_length(r_entry2);
	      q3 = length2 - 1 - get_align_entry(r_entry2)->qual_end;
	      q4 = length2 - 1 - get_align_entry(r_entry2)->qual_start;
	    }
	    else {
	      q3 = get_align_entry(r_entry2)->qual_start;
	      q4 = get_align_entry(r_entry2)->qual_end;
	    }
	    get_LLR(pair->diffs, get_orig_qual(r_entry1), get_orig_qual(r_entry2), 
		    get_adj_qual(r_entry1), get_adj_qual(r_entry2), get_seq_length(r_entry1),
		    pair->start1, pair->end1, get_seq_length(r_entry2), 
		    pair->start2, pair->end2, is_reverse(pair), 
		    get_align_entry(r_entry1)->segments, get_align_entry(r_entry2)->segments, 
		    1, fp, q1, q2, q3, q4, 0);
	    print_pair(pair, fp, 1);


	  }
	}
      }
    }
  }
/*  fprintf(fp, "\n\nContig %d", contig1->index); */
  fprintf(fp, "\n\nSlack, # used pairs (max_score), unused");
  for (i = 0; i < 100; i++) {
    if (used_hist1[i][0] || unused_hist1[i][0] 
	|| diff_offset_hist[i]) 
      fprintf(fp, "\n%2d  %4d  (%4.1f)  %4d (%4.1f)     %4d", 
	     i, used_hist1[i][0], used_hist1[i][1] / 10.0,
	      unused_hist1[i][0], unused_hist1[i][1] / 10.0,
	      diff_offset_hist[i]);
  }
  printf("\n\nLLR histograms (used, unused pairs): ");
  for (i = 0; i < 50; i++)
    if (LLR_hist[0][i] || LLR_hist[1][i])
      printf("\n%5.1f  %4d    %4d", i - 10.0, LLR_hist[0][i], LLR_hist[1][i]);
}

make_contig_db()
{
  int i;
  char *our_alloc();
  Contig *contig1;
  char *db_file;
  
  notify("Making contig sequences ...");
  db_file = parameters->query_files->name;
  set_contig_graph_weights();

/*  reset_extents(query_db);  Check whether this is appropriate! -- will distort output */
/* set other structure members here? set file_name that will be written to
   -- then put routine in db.c */
/* process in reverse order so that memory requests for largest contig will go in first
   to reduce risk of wasted allocation */
  for (i = num_contigs - 1; i >= num_singletons + num_duplicates; i--) {
    contig1 = contig_ptrs[i];
    contig1->id = (char *)our_alloc(strlen(db_file) + 15);
    sprintf(contig1->id, "Contig%d", contig1->index);
    contig1->descrip = 0;
    contig1->seq = 0;
    contig1->orig_qual = contig1->adj_qual = 0;
/*    find_depths(contig1); */
    node_master(contig1, query_db);
    continue;
 
/*
   align_entry = contig1->first;

    db_entry = align_entry->db_entry + (align_entry->reverse ? num_query_entries : 0);
    i_seq = i_seq1 = 0;
    do {
      max_final_end = length1 = db_entry->length - 1;
      max_end = align_entry->reverse ? 
	length1 - align_entry->first_start : align_entry->last_end;
      best_pair = best_pair2 = best_pair3 = 0;
      best_score = best_score2 = best_score3 = 0;
      for (pair = align_entry->first_pair; pair; pair = pair->next) {
	if (!is_used(pair)) continue; 
	if (pair->score <= best_score) continue;
	align_entry2 = align_array + pair->entry2; 
	if (align_entry2->used) continue;
	if (align_entry2->contig != align_entry->contig 
	    || is_reverse(pair) !=
	    (align_entry->reverse != align_entry2->reverse)) {
	  printf("\nERROR: used pair in wrong contig or orientation\n");
	  continue;
	}
	diff = align_entry->reverse ?
	  (align_entry->end - pair_end1) - (align_entry2->end - pair_end2)
	    : align_entry->start + pair_start1 -
	      (align_entry2->start + pair_start2);
	if (abs(diff) > LOCATION_FUDGE1) {
	  printf("\nERROR: used pair with wrong displacement\n");
	  continue;
	}
	length2 = align_entry2->db_entry->length - 1;
	diff1 = align_entry->reverse ?
	  length1 - pair_start1 - (length2 - pair_start2)
	    : pair_end1 - pair_end2;
	end2 = diff1 + (align_entry2->reverse ?
			length2 - align_entry2->first_start : align_entry2->last_end);
	final_end2 = diff1 + length2;
	
	start1 = align_entry->reverse ? length1 - pair_end1 : pair_start1;
	
	if (start1 <= i_seq1) {
	  start1 = align_entry->reverse ? length1 - pair_start1 : pair_end1;
	  pos_flag = 1;
	}
	else pos_flag = 0;
	
	if (start1 <= i_seq1) continue;
	
	if (end2 > max_end) {
	  if (pos_flag) {
	    if (best_score2 < pair->score) {
	      best_pair2 = pair;
	      best_score2 = pair->score;
	    }
	  }
	  else {
	    best_pair = pair;
	    best_score = pair->score;
	  }
	}
	else if (end2 == max_end && final_end2 > max_final_end && pair->score > best_score3) {
	  best_pair3 = pair;
	  best_score3 = pair->score;
	}
      }
      if (!best_pair) best_pair = best_pair2 ? best_pair2 : best_pair3;
      if (best_pair) { 
	align_entry2 = align_array + best_pair->entry2;
	db_entry2 = align_entry2->db_entry +
	  (align_entry2->reverse ? num_query_entries : 0);
	length2 = db_entry2->length - 1;
	start1 = align_entry->reverse ?
	  length1 - best_pair_end1 : best_pair_start1;
	start2 = align_entry->reverse ?
	  length2 - best_pair_end2 : best_pair_start2;
	if (start1 <= i_seq1) {
	  start1 = align_entry->reverse ?
	    length1 - best_pair_start1 : best_pair_end1;
	  start2 = align_entry->reverse ?
	    length2 - best_pair_start2 : best_pair_end2;
	}
	if (start1 <= i_seq1) {
	  printf("\n%s  %s  %d  %d  %d  %d\n  %d %d %d %d %d %d %d %d %d %d %d",
		 db_entry->id, db_entry2->id,
		 i_seq1, start1, length1, length2,
		 align_entry->reverse, align_entry2->reverse,
		 align_entry->first_start, align_entry->last_end,
		 align_entry2->first_start, align_entry2->last_end,
		 best_pair_start1, best_pair_start2,
		 best_pair_end1, best_pair_end2,
		 best_pair->reverse);
	}
      }
      else start1 = length1 + 1;
      while (i_seq1 < start1) 
	contig_entry->seq[i_seq++] = db_entry->seq[i_seq1++];
      if (i_seq >= contig1->length + 2000) {
	fprintf(stderr,"\n%s %d %d",
		db_entry->id, i_seq, contig1->length);
	fatalError("Contig sequence allocation exceeded");
      }
      align_entry->used = 1;
      i_seq1 = start2;
      align_entry = align_entry2;
      db_entry = db_entry2;
    } while (best_pair);
    contig_entry->seq[i_seq] = 0;
    contig_entry->length = contig1->length = i_seq;
*/
  }
  notify(" Done\n");
}

reverse_base_segments( contig1, num_base_segs_ptr ) 
  Contig *contig1;
  int *num_base_segs_ptr;
{
   Base_segment *previous_base_segment, *current_base_segment, *next_base_segment;

   *num_base_segs_ptr = 0;
   previous_base_segment = 0;
   current_base_segment = contig1->base_segment;
   while ( current_base_segment ) {
      ++(*num_base_segs_ptr);
      next_base_segment = current_base_segment->next;
      current_base_segment->next = previous_base_segment;
      previous_base_segment = current_base_segment;
      current_base_segment = next_base_segment;
   }

   /* now make the new head of the list what was the tail of the list */
   contig1->base_segment = previous_base_segment;
}

write_contigs()
{
  char *our_alloc();
  char *contig_file, *contig_name, *contig_qual_file, *ace_file;
  FILE *fp, *fr, *fa, *fv;
  FILE *fopenWrite();
  int i, j, i_size, i_seq, i_ins, contig_j, contig_i_seq;
  char read_flag, contig_flag;
  int *insert_size, *translate;
  Contig *contig1;
  unsigned char *diff;
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *get_aligned_pairs();
  Aligned_pair *pair;
  char *db_file,*view_file;
  char c;
  int l_trim, r_trim, l_pad_trim, r_pad_trim, length, length1, length2;
  int l_qual_trim, r_qual_trim, l_pad_qual_trim, r_pad_qual_trim;
  int l_align_trim, r_align_trim, l_pad_align_trim, r_pad_align_trim;
  Base_segment *base_segment;
  int c_start, c_end, r_start, r_end;
  int site1, site2, type, seq_num;
  char *get_id(), *get_descrip(), *get_adj_qual();
  unsigned char *get_seq(), *get_comp_seq();
  unsigned char *seq;
  char *adj_qual, *id, *descrip;
  Align_info *get_align_entry();
  int r_start1, r_end1, c_start1, c_end1, c_start2, c_end2, c_diff, orient, unalign, displace;
  int n_chimeras, num_reads;
  Segment *segment;
  char new_flag;
  int whole_read_infos, read_tags, num_base_segs;
  int padded_align_start, padded_align_end, padded_bases, n;
  int qual_start, qual_end, temp;
  Tag *tag;

  notify("Writing contigs ...");
  db_file = parameters->query_files->name;
  contig_file = (char *)our_alloc(strlen(db_file) + 10);
  contig_name = (char *)our_alloc(strlen(db_file) + 15);
  contig_qual_file = (char *)our_alloc(strlen(db_file) + 15);
  strcpy(contig_file, db_file);
  strcat(contig_file,".singlets");
  fp = fopenWrite(contig_file);
  for (i = num_singletons + num_duplicates, num_reads = 0; i < num_contigs; i++) 
    num_reads += contig_ptrs[i]->num_entries;
  if (parameters->view) {
    view_file = (char *)our_alloc(strlen(db_file) + 10);
    strcpy(view_file, db_file);
    strcat(view_file, ".view");
    fv = fopenWrite(view_file);
    fprintf(fv,"HEADER FORMAT1 %s %d %d %d %d %d\n", db_file, 
	    num_reads + num_singletons + num_duplicates,
	    num_duplicates, num_singletons, num_reads, 
	    num_contigs - num_singletons - num_duplicates);
  }
  
  for (i = num_duplicates; i < num_singletons + num_duplicates; i++) {
    seq_num = contig_ptrs[i]->first->seq_entry;
    id = get_id(seq_num);
    write_entry(fp, id, get_descrip(seq_num), get_seq(seq_num));
    if (parameters->view) fprintf(fv,"READ %s %d\n", id, get_seq_length(seq_num));
  }
  fclose(fp);

  n_chimeras = 0;
  strcpy(contig_file, db_file);
  strcat(contig_file,".contigs");
  fp = fopenWrite(contig_file);
  for (i = num_singletons + num_duplicates; i < num_contigs; i++) {
    contig1= contig_ptrs[i];
    sprintf(contig_name, "%s.%s", db_file, contig1->id);
    write_entry(fp, contig_name, contig1->descrip, contig1->seq);
    if (!parameters->view) continue;
    fprintf(fv,"CONTIG %s %d\nCONTIG_QUAL %s", contig1->id, contig1->length, contig1->id);
    for (j = 0; j < contig1->length; j++) fprintf(fv," %d", (int) contig1->adj_qual[j]);
    fprintf(fv,"\nDISCREP_QUAL %s", contig1->id);
    if (contig1->num_entries > 1) 
      for (j = 0; j < contig1->length; j++) 
	fprintf(fv," %d", (int) contig1->discrep_qual[j]);
    else 
      for (j = 0; j < contig1->length; j++) fprintf(fv," 0");
    fprintf(fv,"\n");
    for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) {
      length1 = get_seq_length(align_entry->seq_entry);
      id = get_id(align_entry->seq_entry);
      r_start1 = align_entry->reverse ? length1 - align_entry->m_end : align_entry->m_start - 1;
      r_end1 = align_entry->reverse ? length1 - align_entry->m_start : align_entry->m_end - 1;
      c_start1 = align_entry->start + align_entry->m_start - 2;
      c_end1 = align_entry->end - (length1 - align_entry->m_end) - 1;
      fprintf(fv,"READ %s %d %.1f %d %d %s %d %d", id, length1, align_entry->LLR_score / 10.0,
	      r_start1, r_end1, contig1->id, 
	      align_entry->reverse ? c_end1 : c_start1, 
	      align_entry->reverse ? c_start1 : c_end1);
      for (segment = align_entry->segments; segment; segment = segment->next) {
	if (segment->start > r_end1 - 10 || segment->end < r_start1 + 10) {
	  n_chimeras++;
	  fprintf(fv," Chimera%d %d %d", n_chimeras, segment->start, segment->end);
	}
      }
      fprintf(fv,"\n");
      for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
	if (pair->entry2 < pair->entry1) continue;
	length2 = get_seq_length(pair->entry2) - 1;
	align_entry2 =  get_align_entry(pair->entry2);
	c_start2 = align_entry2->start + align_entry2->m_start - 2;
	c_end2 = align_entry2->end - (length2 + 1 - align_entry2->m_end) - 1;
	c_diff = align_entry2->contig != contig1;
	orient =  is_reverse(pair) != (align_entry->reverse != align_entry2->reverse);
	unalign = 0;
	site1 = find_position(align_entry, pair->start1) - 1;
	if (unalign < c_start1 - site1) unalign = c_start1 - site1;
	if (unalign < site1 - c_end1) unalign = site1 - c_end1;
	site2 = find_position(align_entry2, 
			      is_reverse(pair) ? length2 - pair->start2 : pair->start2) - 1;
	if (unalign < c_start2 - site2) unalign = c_start2 - site2;
	if (unalign < site2 - c_end2) unalign = site2 - c_end2;
	displace = abs(site1 - site2);
	site1 = find_position(align_entry, pair->end1) - 1;
	if (unalign < c_start1 - site1) unalign = c_start1 - site1;
	if (unalign < site1 - c_end1) unalign = site1 - c_end1;
	
	site2 = find_position(align_entry2, 
			      is_reverse(pair) ? length2 - pair->end2 : pair->end2) - 1;
	if (unalign < c_start2 - site2) unalign = c_start2 - site2;
	if (unalign < site2 - c_end2) unalign = site2 - c_end2;
	if (displace > abs(site1 - site2)) displace = abs(site1 - site2);
	
	if (pair->LLR_score < 0 || c_diff || orient || displace > 30 || unalign > 25) 
	  fprintf(fv,"MATCH %.1f %s %d %d %s %d %d %d %d\n", pair->LLR_score / 10.0, 
		  id, pair->start1, pair->end1,
		  get_id(pair->entry2), 
		  is_reverse(pair) ? length2 - pair->start2 : pair->start2,
		  is_reverse(pair) ? length2 - pair->end2 : pair->end2,
		  c_diff || orient || displace > 30, unalign);
	
      }
    }
  } /* write contig length (contig1->length) as descrip? */
  fclose(fp);
  if (parameters->view) fclose(fv);

  strcpy(contig_qual_file, contig_file);
  strcat(contig_qual_file, ".qual");
  fr = fopenWrite(contig_qual_file);
  our_free(contig_qual_file);
  our_free(contig_file);
  for (i = num_singletons + num_duplicates; i < num_contigs; i++) {
    contig1 = contig_ptrs[i];
    fprintf(fr,">%s.%s ", db_file, contig1->id);
    for (j = 0; j < contig1->length; j++) {
      if (!(j % 50)) fprintf(fr,"\n");
      fprintf(fr,"%d ", (int)contig1->adj_qual[j]);
    }
    fprintf(fr,"\n");
  }

  fclose(fr);
  if (parameters->exp_output) {
    call_write_exp_files(num_singletons + num_duplicates, num_contigs, contig_ptrs);
    notify(" Done\n");
    return;
  }

  if (!parameters->old_ace && !parameters->new_ace) { 
    notify(" Done\n");
    return; /* the rest (which finds padded sequences) is only need for .ace file */
  }

  new_flag = parameters->new_ace;

  ace_file = (char *)our_alloc(strlen(db_file) + 10);
  strcpy(ace_file, db_file);
  strcat(ace_file, ".ace");
  fa = fopenWrite(ace_file);
  our_free(ace_file);

  if (new_flag) {
    fprintf(fa,"AS %d %d\n",
            num_contigs - (num_singletons + num_duplicates),
            num_reads);
  }
  for (i = num_singletons + num_duplicates; i < num_contigs; i++) {
    contig1 = contig_ptrs[i];
    if (contig1->num_entries > 1) 
      fprintf(fp_log, "\n\nContig %d unpadded => padded conversion:", contig1->index);
    insert_size = (int *)our_alloc((contig1->length + 1) * sizeof(int));
/* translate converts between sequence without pads and sequence with pads */
    translate = (int *)our_alloc((contig1->length + 1) * sizeof(int));
    for (j = 0; j <= contig1->length; j++) insert_size[j] = 0;
    for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) {
      site1 =  align_entry->start + align_entry->m_start - 2; 
      if (diff = align_entry->diffs)
	for (; *(diff + 1); diff++) {
	  site1 += diff_gap1(*diff);
	  if (diff_type(*diff) == 'I') {
	    i_size = 1;
	    while (diff_type(*(diff + 1)) == 'I' && !diff_gap1(*(diff + 1))) {
	      diff++;
	      i_size++;
	    }
	    if (i_size > insert_size[site1])
	      insert_size[site1] = i_size;
	  }
	}
    }

    if (new_flag) {
       /* figure out number of padded bases */
       padded_bases = contig1->length;
       for( n = 0; n< contig1->length; ++n )
         padded_bases += insert_size[n];

       /* put base segments in left to right order for printing out */
       reverse_base_segments( contig1, &num_base_segs );
       fprintf(fa,"\n\nCO %s %d %d %d U", contig1->id, padded_bases, contig1->num_entries, num_base_segs );
    }
    else fprintf(fa,"\n\nDNA %s", contig1->id);

    for (j = i_seq = i_ins = 0; j < contig1->length; i_seq++) {
      if (!(i_seq % 50)) 
	fprintf(fa,"\n");
      if (contig1->num_entries > 1 && j && !(j % 200)) 
	fprintf(fp_log, "\n%5d => %5d", j, i_seq);

      if (i_ins < insert_size[j]) {
	fprintf(fa,"*");
	i_ins++;
      }
      else {
	fprintf(fa, "%c", contig1->adj_qual[j] >= parameters->qual_show ? 
		contig1->seq[j] : tolower(contig1->seq[j]));
	j++;
	translate[j] = i_seq + 1;
	i_ins = 0;
      }
    }
    if (new_flag)
      fprintf(fa,"\n\nBQ");
    else
      fprintf(fa,"\n\n\nBaseQuality %s", contig1->id);
    for (j = 0; j < contig1->length; j++) {
      if (!(j % 50)) fprintf(fa,"\n");
      fprintf(fa, " %d", contig1->adj_qual[j]);
    }
    contig1->n_pads = translate[j] - j;

    if (!new_flag) fprintf(fa,"\n\n\nSequence %s", contig1->id);
    else fprintf(fa, "\n" );
    for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) {
      if (!new_flag) {
	fprintf(fa,"\nAssembled_from  %s%s  %d  %d", 
		get_id(align_entry->seq_entry), align_entry->reverse ? ".comp" : "", 
		align_entry->start, align_entry->end);
      fprintf(fa,"\nAssembled_from* %s%s  %d  %d",     
	      get_id(align_entry->seq_entry), align_entry->reverse ? ".comp" : "", 
	      align_entry->start > 0 && align_entry->start <= contig1->length ? 
	      translate[align_entry->start] : align_entry->start, 
	      align_entry->end <= contig1->length ? translate[align_entry->end]
	      : align_entry->end + contig1->n_pads);
      }
      else
         fprintf(fa, "\nAF %s %c %d", get_id(align_entry->seq_entry), 
          align_entry->reverse ? 'C' : 'U',       
	      align_entry->start > 0 && align_entry->start <= contig1->length ? 
	      translate[align_entry->start] : align_entry->start );
    }

    for (base_segment = contig1->base_segment; base_segment; base_segment = base_segment->next) {
      align_entry = get_align_entry(base_segment->entry);
      if (!new_flag) {
	fprintf(fa,"\nBase_segment  %d %d  %s%s  %d %d",
		base_segment->contig_start, base_segment->contig_end,
		get_id(base_segment->entry), align_entry->reverse ? ".comp" : "", 
		base_segment->read_start, base_segment->read_end);
      }
      c_start = translate[base_segment->contig_start];
      c_end = translate[base_segment->contig_end];
      r_start = c_start + 1 - (align_entry->start > 0 ? 
			       translate[align_entry->start] : align_entry->start);
      r_end = r_start + (c_end - c_start); 
      if (new_flag)
        fprintf(fa,"\nBS %d %d %s",
  	      c_start, c_end,
	      get_id(base_segment->entry) );
      else
        fprintf(fa,"\nBase_segment* %d %d  %s%s  %d %d",
  	      c_start, c_end,
	      get_id(base_segment->entry), align_entry->reverse ? ".comp" : "", 
	      r_start, r_end);
    }

    for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) {
      seq_num = align_entry->seq_entry;
      length1 = get_seq_length(seq_num);

      if (parameters->screen) {
	seq = get_seq(seq_num);
	for (j = 0; j <= align_entry->last_vec; j++) seq[j] = 'X';
	
/* N.B. THIS ASSUMES SEQUENCE IS IN MEMORY -- IF NOT IT WILL FAIL !! */
	r_start1 = align_entry->reverse ? length1 - align_entry->m_end : align_entry->m_start - 1;
	r_end1 = align_entry->reverse ? length1 - align_entry->m_start : align_entry->m_end - 1;
	for (segment = align_entry->segments; segment; segment = segment->next) {
	  if (segment->start > r_end1 - 10 || segment->end < r_start1 + 10) 
	    for (j = segment->start; j <= segment->end; j++) seq[j] = 'X';
	    
	}
      }
      seq = align_entry->reverse ? get_comp_seq(seq_num) : get_seq(seq_num);
      id = get_id(seq_num);
      descrip = get_descrip(seq_num);
      /* in following contig_j and contig_i_seq (positions in original and padded
	 contig sequences, resp.) have origin 1, while 
	 j and i_seq (positions in original and padded read) have origin 0 (Sorry!) */
      contig_j = align_entry->start;
      contig_i_seq = (contig_j > 0 && contig_j <= contig1->length) ?
	translate[contig_j] : contig_j;
      diff = align_entry->diffs;
      site2 = align_entry->m_start - 1 + diff_gap2(*diff);


      if (new_flag) {

	if (align_entry->m_end < align_entry->m_start + 10)
	   append_tag(align_entry, "unaligned", -1, -1);

         padded_align_start = 
	      align_entry->start > 0 && align_entry->start <= contig1->length ? 
	      translate[align_entry->start] 
            : align_entry->start;

         padded_align_end =
	      align_entry->end <= contig1->length ? 
            translate[align_entry->end]
              : align_entry->end + contig1->n_pads;


         padded_bases = padded_align_end - padded_align_start + 1;

         count_tags( align_entry, &whole_read_infos, &read_tags );

         fprintf(fa, "\n\nRD %s %d %d %d", 
	      id, padded_bases, whole_read_infos, read_tags ); 

      }
      else 
        fprintf(fa, "\n\n%s %s%s", "DNA",
	      id, align_entry->reverse ? ".comp" : ""); 

      length = length1 - 1;
      adj_qual = get_adj_qual(seq_num);
      if (align_entry->reverse) {
	for (l_trim = 0; l_trim <= length && l_trim < align_entry->m_start - 1
	     && adj_qual[length - l_trim] < parameters->qual_show; l_trim++);
	for (r_trim = length; r_trim > l_trim && r_trim > align_entry->m_end - 1
	     && adj_qual[length - r_trim] < parameters->qual_show; r_trim--);
      }
      else {
	for (l_trim = 0; l_trim <= length && l_trim < align_entry->m_start - 1
	     && adj_qual[l_trim] < parameters->qual_show; l_trim++);
	for (r_trim = length; r_trim > l_trim && r_trim > align_entry->m_end - 1
	     && adj_qual[r_trim] < parameters->qual_show; r_trim--);
      }

      find_trimmed_quals(align_entry, &qual_start, &qual_end);
      l_qual_trim = align_entry->reverse ? length - qual_end : qual_start;
      r_qual_trim = align_entry->reverse ? length - qual_start : qual_end;
      l_align_trim = align_entry->m_start - 1;
      r_align_trim = align_entry->m_end - 1;
      if (align_entry->reverse) {
	for (tag = align_entry->tags; tag; tag = tag->next) {
	  temp = tag->start;
	  tag->start = length - tag->end;
	  tag->end = length - temp;
	}
      }
      /* just to check that these values are set */

      l_pad_qual_trim = -666;
      r_pad_qual_trim = -666;


      for (j = i_seq = 0; j < length1; i_seq++, contig_i_seq++) {
	if (!(i_seq % 50)) {
	  fprintf(fa,"\n");
	}
/*
	read_flag = contig_flag = 1; 
	if (diff && diff->type != 'E' && diff->site2 <= j + 1) {
	  if (diff->type == 'D') {
	    if (diff->site2 == j) read_flag = 0;
	    else diff--; 
	  }
	  else if (diff->type == 'I') contig_flag = 0;
	  diff++;
	}
	else if (contig_j > 0 && contig_j <= contig1->length 
		 && contig_i_seq < translate[contig_j]) 
	  read_flag = contig_flag = 0;
*/
	if (contig_j > 0 && contig_j <= contig1->length 
	    && contig_i_seq < translate[contig_j]) 
	  read_flag = contig_flag = 0;
	else read_flag = contig_flag = 1; 

	if (*(diff + 1) && site2 <= j + 1) {
	  if ((type = diff_type(*diff)) == 'D') {
	    if (site2 <= j && contig_flag) read_flag = 0;
	    else {
	      site2 -= diff_gap2(*diff);
	      diff--; 
	    }
	  }
	  else if (type == 'I') {
	    if (!contig_flag) read_flag = 1;
	    else notify("\nWARNING: Alignment error I");
	  }
	  diff++;
	  site2 += diff_gap2(*diff);
	}

	if (!read_flag) c = '*';
	else {
	  if (adj_qual[align_entry->reverse ? 
				length - j : j] >= parameters->qual_show) 
	    c = seq[j];
	  else c = tolower(seq[j]);
	  if (j == l_trim) l_pad_trim = i_seq;
	  if (j == r_trim) r_pad_trim = i_seq;
	  if (j == l_qual_trim) l_pad_qual_trim = i_seq;
	  if (j == r_qual_trim) r_pad_qual_trim = i_seq;
	  if (j == l_align_trim) l_pad_align_trim = i_seq;
	  if (j == r_align_trim) r_pad_align_trim = i_seq;
	  for (tag = align_entry->tags; tag; tag = tag->next) {
	    if (j == tag->start) tag->start = -2 - i_seq;
	    if (j == tag->end) tag->end = -2 - i_seq;
	  }
	  j++;
	}
	fprintf(fa,"%c", c);
	if (contig_flag) contig_j++;
      }
      if (new_flag) {
         fprintf(fa, "\n" );
	 if (qual_start > qual_end ) {
	   /* signal consed that the whole read is low quality */
	   l_pad_qual_trim = -1;
	   r_pad_qual_trim = -1;
	 }
	 else {
	   /* convert to 1-based coordinates */

	   ++l_pad_qual_trim;
	   ++r_pad_qual_trim;
	 }
	 if (r_align_trim < l_align_trim + 10) {
	   /* unaligned read */
	   l_pad_align_trim = -2;
	   r_pad_align_trim = -2;
	 }
	 
	 fprintf(fa, "\nQA %d %d %d %d",
		 l_pad_qual_trim, r_pad_qual_trim,
		 l_pad_align_trim + 1, r_pad_align_trim + 1 );
	 fprintf(fa, "\nDS");
         if (descrip && descrip[0])
           fprintf(fa, " %s", descrip);
	 write_tags(fa, align_entry);
      }
      else {
         fprintf(fa, "\n\nSequence %s%s",
                 id, align_entry->reverse ? ".comp" : ""); 

         fprintf(fa, "\nClipping  %d %d", l_trim + 1, r_trim + 1);
         fprintf(fa, "\nClipping* %d %d", l_pad_trim + 1, r_pad_trim + 1); 
	 if (descrip && descrip[0]) 
	   fprintf(fa, "\nDescription  %s", descrip);
      }
    }
    our_free(insert_size);
    our_free(translate);
  }
  fprintf(fa, "\n\n");
  if (new_flag) {

    fprintf(fa,"\nWA{\nphrap_params phrap %s\n", parameters->date);
    for (i = 0; i < parameters->argc; i++) 
      fprintf(fa, "%s ", parameters->argv[i]);
    fprintf(fa,"\nphrap version %s\n}\n\n", parameters->version);

    fprintf(fa,"\nWA{\nsinglets phrap %s\n", parameters->date);
    fprintf(fa,"SN %d\n", num_singletons);
    for (i = num_duplicates; i < num_singletons + num_duplicates; i++) {
      seq_num = contig_ptrs[i]->first->seq_entry;
      fprintf(fa,"%s %s\n", get_id(seq_num), get_descrip(seq_num));
    }
    fprintf(fa,"}\n\n");
  }
  fclose(fa);
  notify(" Done\n");
}

revise_contigs()
{
  int i;
  Contig *contig1;

  notify("Revising contigs ...");
/* reverse order so biggest contigs done first */
  for (i = num_contigs - 1; i >= num_singletons + num_duplicates; i--) {
    contig1 = contig_ptrs[i];
    if (contig1->num_entries > 1) contig_revise(contig1);
  }
  notify(" Done\n");
  count_base_segments();
}


align_reads_to_contigs(pass)
     int pass;
{
  int i, j, score, slop, temp;
  int full_smith_waterman();
  Contig *contig1;
  Profile *make_profile_from_seq();
  Profile *q_profile;
  Align_info *align_entry;
  int start1, end1, start2, end2;
  int l_edge, r_edge, length1, seq_num;
  Segment *insert_segment();
  unsigned char *make_diffs();
  int mismatches, insertions, deletions;
  int orig_score;
  int n_successes, n_failures;
  char *get_id();
  unsigned char *get_seq(), *get_comp_seq();
  char *our_alloc();
  unsigned char *seq;
  char *qual, *qual_buffer;
  int qual_buffer_length;
  char *get_orig_qual();

  notify("Aligning reads to contigs ...");

  qual_buffer_length = 0;
  qual = qual_buffer = 0;

  set_score_mat();
  set_gap_penalties(parameters->gap_init, parameters->gap_ext, parameters->ins_gap_ext, parameters->del_gap_ext, 0);
/*
  if (pass == 1)
    set_score_mat(parameters->matrix, parameters->gap_init, parameters->gap_ext, 0, 
		parameters->penalty);
  else set_score_mat(parameters->matrix, -3, -2, 0, -1);
*/

  n_successes = n_failures = 0;
  for (i = num_singletons + num_duplicates; i < num_contigs; i++) {
    contig1 = contig_ptrs[i];
    q_profile = make_profile_from_seq((Profile *)0, contig1->seq, 0, 0); /* instead put length in */
    contig1->top_segments = contig1->bottom_segments = 0;
    contig1->first_start = contig1->length;
    contig1->last_end = 0;
    for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) {
      seq_num = align_entry->seq_entry;
      seq = align_entry->reverse ? get_comp_seq(seq_num) : get_seq(seq_num);

      length1 = get_seq_length(seq_num);
      goto skip_pairs;
/* SKIP FOLLOWING FOR NOW 
      l_edge = r_edge = align_entry->start;
      for (pair = align_entry->first_pair; pair; pair = pair->next) {
	if (!is_used(pair)) continue;
	align_entry2 = pair->entry2 + align_array;
	if (align_entry2->score < parameters->minscore  N.B. also covers case where
					       score = 0, i.e. entry not yet
					       processed in this contig 
	    || align_entry2->contig != align_entry->contig
	    || is_reverse(pair) != (align_entry->reverse != align_entry2->reverse))
	  continue;
	db_entry2 = align_entry2->db_entry;
	length2 = db_entry2->length;
	pair2 = pair->reversed_pair;
	start = align_entry2->reverse ?
	  (length2 - pair2->end1) - (length1 - pair2->end2)
	    : pair2->start1 - pair2->start2;
	start +=  align_entry2->start;  start now is predicted relative offset of
					  align_entry w.r.t contig 
	if (abs(start - align_entry->start) < 150) {
	  if (start < l_edge) l_edge = start;
	  if (start > r_edge) r_edge = start;
	}
	else printf("\nMISMATCH1: %s  %s  %d  %d  score: %d",
		    db_entry->id, db_entry2->id,
		    start, align_entry->start, pair->score);
	start = align_entry2->reverse ?
	  (length2 - pair2->start1) - (length1 - pair2->start2)
	    : pair2->end1 - pair2->end2;
	start +=  align_entry2->start; 
	if (abs(start - align_entry->start) < 150) {
	  if (start < l_edge) l_edge = start;
	  if (start > r_edge) r_edge = start;
	}
	else printf("\nMISMATCH2: %s  %s  %d  %d  score: %d",
		    db_entry->id, db_entry2->id,
		    start, align_entry->start, pair->score);
      }
*/
    skip_pairs:
      l_edge = align_entry->start;
      r_edge = align_entry->end - length1 - 1;
      if (r_edge < l_edge) {
	temp = r_edge;
	r_edge = l_edge;
	l_edge = temp;
      }
      if (align_entry->start > contig1->length || align_entry->end < 0) {
	fprintf(stderr, "\nALIGNMENT ERROR %s %d %d  %d", get_id(seq_num), align_entry->start, 
		align_entry->end, contig1->length);
	fprintf(stdout, "\nALIGNMENT ERROR %s %d %d  %d", get_id(seq_num), align_entry->start, 
		align_entry->end, contig1->length);
      }
      slop = 45;
/*
      score = quick_full_smith_waterman(q_profile, db_entry->seq, db_entry->length,
				  l_edge - slop, r_edge + slop, 
				  0, 0, 0, 0, 0, &orig_score, &success);
      success ? n_successes++ : n_failures++;
*/
      score = full_smith_waterman(q_profile, seq, length1, l_edge - slop, r_edge + slop, 
				  0, 0, 0, 0, 0, &orig_score);

      get_stats(&start1, &end1, &start2, &end2,
		&mismatches, &insertions, &deletions);
      if (!score) {
	start1 = end1 = start2 = end2 = 1;
      }
      align_entry->start = start1 - start2 + 1;
      align_entry->end = end1 + (length1 - end2);
      l_edge = align_entry->start + 
	(align_entry->reverse ? length1 - 1 - align_entry->last_end 
	 : align_entry->first_start);
      r_edge = align_entry->end - 
	(align_entry->reverse ? align_entry->first_start : length1 - 1 - align_entry->last_end );
      /* note above are still approximate! */
      if (score >= parameters->minscore && l_edge < contig1->first_start) 
	contig1->first_start = l_edge;
      if (score >= parameters->minscore && r_edge > contig1->last_end) 
	contig1->last_end = r_edge;
       
      if (score) {
	qual = get_orig_qual(align_entry->seq_entry);
	if (!align_entry->reverse) {
	  if (qual_buffer_length < length1) {
	    if (qual_buffer) our_free(qual_buffer);
	    qual_buffer = (char *)our_alloc(length1 * sizeof(char));
	    qual_buffer_length = length1;
	  }
	  for (j = 0; j < length1; j++)
	    qual_buffer[j] = qual[length1 - 1 - j];
	  qual = qual_buffer;
	}
	
	qual += length1 - end2 -1;
	if (pass) qual_slide_indels(qual, 1);

	align_entry->diffs = make_diffs(); 

	align_entry->LLR_score = 
	  get_LLR(align_entry->diffs, contig1->orig_qual, align_entry->orig_qual, 
		  contig1->adj_qual, align_entry->adj_qual, 
		  contig1->length, start1 - 1, end1 - 1, 
		  length1, start2 - 1, end2 - 1, align_entry->reverse, (Segment *)0, align_entry->segments, 0, (FILE *)0,
		  0, 0, 0, 0, 0);

      }
      else {
	align_entry->diffs = (unsigned char *)our_alloc(2 * sizeof(char));
	align_entry->diffs[0] = align_entry->diffs[1] = 0;
      }
      
      if (score >= parameters->minscore) {
	if (align_entry->reverse)
	  contig1->bottom_segments = 
	    insert_segment(contig1->bottom_segments, start1, end1);
	else contig1->top_segments = 
	  insert_segment(contig1->top_segments, start1, end1);
      }
      align_entry->m_start = start2;
      align_entry->m_end = end2;
      align_entry->m_length = end1 - start1 + 1;
      align_entry->score = score;
      /*
	if ((end1 = find_position(align_entry, align_entry->reverse ? 0 :db_entry->length - 1)) != 
	align_entry->end)
	fprintf(stderr, "\nMismatch: %d %d  %d", 
	end1, align_entry->end, align_entry->reverse); 
	*/
    }
    free_profile(q_profile);
  }
/*
  fprintf(stderr, "%d successes, %d failures", n_successes, n_failures);
*/
  if (qual_buffer) our_free(qual_buffer);
  notify(" Done\n");
}

/*
debug_print(i)
     int i;
{
  Diff *diff;
  Align_info *align_entry;
  Db_entry *db_entry;
  
  align_entry = align_array + 135;
  for (diff = align_entry->diffs; diff->type != 'E'; diff++) {
    if (diff->site2 > 600) 
      fprintf(stderr, "\n%d %c %d %d", i, diff->type, diff->site1, diff->site2);
  }
}
*/

/* print unattached singleton contigs */
print_singletons()
{
  int i, j, n, seq_num;
  Align_info *align_entry;
  char *get_id();
  unsigned char *get_seq();
  unsigned char *seq;

  printf("\n\n%d perfect duplicates", num_duplicates);
  i = 0;
  if (num_duplicates) {
    printf(": \n  Read         Length");
    for (; i < num_duplicates; i++) {
      align_entry = contig_ptrs[i]->first;
      printf("\n %-12s  %4d", get_id(align_entry->seq_entry), get_seq_length(align_entry->seq_entry));
    }
  }
  printf("\n\n%d isolated singlets (having no non-vector match to any other read)", 
	 num_singletons);
  if (num_singletons) {
    printf(": \n  Read         Length      (# trimmed non-X bases)");
    for (; i < num_singletons + num_duplicates; i++) {
      align_entry = contig_ptrs[i]->first;
      seq_num = align_entry->seq_entry;
      seq = get_seq(seq_num);
      j = align_entry->qual_start > align_entry->last_vec ? 
	align_entry->qual_start : align_entry->last_vec + 1;
      for (n = 0; j < align_entry->qual_end; j++) n += seq[j] != 'X';
      printf("\n %-12s  %4d   (%d)", 
	     get_id(align_entry->seq_entry), get_seq_length(align_entry->seq_entry), n);
    }
  }
}

compare_pair_offsets(pair_1, pair_2)
     Aligned_pair **pair_1, **pair_2;
{
  int i, offset_1, offset_2, start1, start2;
  int length1, length2, length3, length4;
  Aligned_pair *pair1, *pair2;
  char rel_orient1, rel_orient2;
  Align_info *align_entry1, *align_entry2, *align_entry3, *align_entry4;
  Align_info *get_align_entry();
  int unalign1, unalign2, unalign3, unalign4;
 
  pair1 = *pair_1;
  pair2 = *pair_2;
  align_entry1 = get_align_entry(pair1->entry1);
  align_entry2 = get_align_entry(pair1->entry2);
  align_entry3 = get_align_entry(pair2->entry1);
  align_entry4 = get_align_entry(pair2->entry2);
  
  if (i = align_entry1->contig->index - align_entry3->contig->index)
    return i;
  if (i = align_entry2->contig->index - align_entry4->contig->index)
    return i;
  rel_orient1 = (is_reverse(pair1) != (align_entry1->reverse != align_entry2->reverse));
  rel_orient2 = (is_reverse(pair2) != (align_entry3->reverse != align_entry4->reverse));
  if (i = rel_orient1 - rel_orient2) return i;

  unalign1 = is_unaligned(pair1);
  unalign2 = is_unaligned(pair2);
  unalign3 = is_unaligned(pair1->reversed_pair);
  unalign4 = is_unaligned(pair2->reversed_pair);
  if (i = unalign1 - unalign2) return i;
  if (unalign1 && (i = pair1->entry1 - pair2->entry1)) return i;
  if (i = unalign3 - unalign4) return i;
  if (unalign3 && (i = pair1->entry2 - pair2->entry2)) return i;

  length1 = get_seq_length(pair1->entry1);
  length2 = get_seq_length(pair1->entry2);
  length3 = get_seq_length(pair2->entry1);
  length4 = get_seq_length(pair2->entry2);
  
  if (align_entry1->reverse) {
    start1 = unalign1 ? pair1->end1 : find_position(align_entry1, pair1->end1);
    start2 = unalign3 ? (is_reverse(pair1) ? length2 - 1 - pair1->end2 : pair1->end2)
      : find_position(align_entry2, 
		      is_reverse(pair1) ? length2 - 1 - pair1->end2 : pair1->end2);
  }
  else {
    start1 = unalign1 ? pair1->start1 : find_position(align_entry1, pair1->start1);
    start2 = unalign3 ? (is_reverse(pair1) ? length2 - 1 - pair1->start2 : pair1->start2)
      : find_position(align_entry2, 
		      is_reverse(pair1) ? length2 - 1 - pair1->start2 : pair1->start2);
  } 
  if (rel_orient1) start2 = -start2;
  offset_1 = start1 - start2;
  
  if (align_entry3->reverse) {
    start1 = unalign2 ? pair2->end1 : find_position(align_entry3, pair2->end1);
    start2 = unalign4 ? (is_reverse(pair2) ? length4 - 1 - pair2->end2 : pair2->end2)
      : find_position(align_entry4, 
		      is_reverse(pair2) ? length4 - 1 - pair2->end2 : pair2->end2);
  }
  else {
    start1 = unalign2 ? pair2->start1 : find_position(align_entry3, pair2->start1);
    start2 = unalign4 ? (is_reverse(pair2) ? length4 - 1 - pair2->start2 : pair2->start2)
      : find_position(align_entry4, 
		      is_reverse(pair2) ? length4 - 1 - pair2->start2 : pair2->start2);
  }
  if (rel_orient2) start2 = -start2;
  offset_2 = start1 - start2;
  
  return offset_1 - offset_2;
}

print_contigs()
{
  int i, j, k, offset, n_entries, start1, start2, end1, end2;
  int i_ptr, i_ptr0, index1, index2, prev_index2, offset1, offset2, prev_offset;
  unsigned char *seq;
  Align_info *align_entry, *align_entry1, *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair **pair_pointers;
  Aligned_pair **sort_pairs();
  Segment *head1, *head2;
  Segment *insert_segment();
  char rel_orient, prev_rel_orient;
  int length1, length2, cross_score, n_lower, n_trimmed;
  Contig *contig1;
  int t_mismatches, t_insertions, t_deletions;
  double length, t_length;
  char *our_alloc();
  int *inserts[2], *deletes[2], *subs[2], *depth_trans[2];
  int strand, first_zer, zer_flag, n_ss, n_errs;
  unsigned char *diff;
  int temp;
  /* signed */ char *qual_flag;
  char *qual_type;
  /* signed */ char q2;
  int mismatches, insertions, deletions;
  int n_matches;
  int site1, site2, diff_t;
  int n_pos_LLR, max_pos_LLR, n_neg_LLR, min_neg_LLR;
  int seq_num;
  char *get_id();
  int unalign1, unalign2, prev_unalign1, prev_unalign2, 
        pair_entry1, pair_entry2, prev_pair_entry1, prev_pair_entry2;
  FILE *fp, *fq;
  FILE *fopenWrite();
  
  char *db_file, *prob_file;
  
  notify("Printing contigs ... ");

  db_file = parameters->query_files->name;
  prob_file = (char *)our_alloc(strlen(db_file) + 15);
  strcpy(prob_file, db_file);
  strcat(prob_file, ".problems");
  fp = fopenWrite(prob_file);
  strcat(prob_file, ".qual");
  fq = fopenWrite(prob_file);
  our_free(prob_file);

  if (parameters->print_extraneous_matches) {
    pair_pointers = sort_pairs(compare_pair_offsets, 1);
    i_ptr = 0;
    for (i_ptr = 0; i_ptr < 2 * num_pairs; i_ptr++) { /* jump through all pairs in which member is
						       in contig of index 0 (should be duplicates
						       only) */
      align_entry1 = get_align_entry(pair_pointers[i_ptr]->entry1);
      if ((int)align_entry1->contig->index > 0) break;
    }
  }
  
  for (i = num_singletons + num_duplicates; i < num_contigs; i++) {
    contig1 = contig_ptrs[i];
/*    transform_contig_quality(contig1); */
    /*
      if (i == num_singletons || contig1->parent != contig_ptrs[i-1]->parent){
      index_first = contig1->index;
      for (j = i + 1;
      j < num_contigs && contig_ptrs[j]->parent == contig1->parent;
      j++);
      index_last = contig_ptrs[j-1]->index;
      }
      */
    n_entries = contig1->num_entries;
    seq = contig1->seq;
    length1 = contig1->length;
     
    for (j = 0; seq[j] && contig1->adj_qual[j] == 0; j++);
    n_trimmed = j;
    for (j = length1 - 1; j > n_trimmed && contig1->adj_qual[j] == 0; j--);
    for (k = n_trimmed, n_lower = 0; k < j; k++) 
      if (contig1->adj_qual[k] == 0) n_lower++;
    n_trimmed += length1 - 1 - j;
    printf("\n\n%sContig %d.  %d read%s; %d bp (untrimmed), %d (trimmed).", 
	   parameters->tags ? "CONTIG    " : "",
	   contig1->index, n_entries, n_entries > 1 ? "s" : "", 
	   length1, length1 - n_trimmed);
    if (contig1->num_contigs <= 1) printf("  Isolated contig.");  
    /*
      length1 -= n_trimmed;
      printf("\n   Trimmed region is %.1f%% lower case; quality score: %d. ", 
      length1 ? (100.0 * n_lower) / (float)length1 : 100.0, 
      contig1->score);
      */
    
    /*   if (contig1->num_contigs > 1)   
	 printf("In linked group, contigs %d - %d.", index_first, index_last); 
	 */
    
    t_mismatches = t_insertions = t_deletions = t_length = 0;

    mark_save_block(); /* so can later free segments created in match_elsewhere; */

    for (align_entry = contig1->first; align_entry;
	 align_entry = align_entry->next) {
      seq_num = align_entry->seq_entry;
      length1 = get_seq_length(seq_num) - 1;

      cross_score = match_elsewhere(align_entry);
      
      get_discreps(align_entry->diffs, &mismatches, &insertions, &deletions);
      length =  align_entry->m_length;
      t_mismatches += mismatches;
      t_insertions += insertions;
      t_deletions += deletions;
      t_length += align_entry->m_length;
      /* align_entry->m_end - align_entry->m_start + 1; */
      
      if (is_deleted(align_entry)) 
	printf("\n ****  PROBABLE DELETION READ");
      
      else if (is_chimera(align_entry) /* align_entry->segments && align_entry->segments->next */) {
	/* potential chimeras (unsupported links) */
	printf("\n ****  PROBABLE CHIMERA  (segments ");
	if (align_entry->reverse)
	  printf("%d-%d, %d-%d) :",
		 1 + length1 - align_entry->segments->next->end,
		 1 + length1 - align_entry->segments->next->start,
		 1 + length1 - align_entry->segments->end,
		 1 + length1 - align_entry->segments->start);
	else printf("%d-%d, %d-%d) :",
		    1 + align_entry->segments->start,
		    1 + align_entry->segments->end,
		    1 + align_entry->segments->next->start,
		    1 + align_entry->segments->next->end);
      }
      
      printf("\n%s%c %5d %5d ", 
	     parameters->tags ? "READ      " : "", align_entry->reverse ? 'C' : ' ',
	     align_entry->start, align_entry->end);
      /* Following is old version -- i.e. positions with respect to padded sequence. No longer
	 can be used here because pad_translate and n_pads are set later -- in write_contigs 
	 
	 align_entry->start > 0 && align_entry->start <= contig1->length ? 
	 contig1->pad_translate[align_entry->start] : align_entry->start, 
	 align_entry->end <= contig1->length ? contig1->pad_translate[align_entry->end]
	 : align_entry->end + contig1->n_pads);
	 */ 
      printf("%-12s %4d (%3d)  %4.2f %4.2f %4.2f %4d (%3d) %4d (%3d) ",
	     get_id(seq_num), align_entry->score, cross_score,
	     length ? 100 * mismatches / length : 0.0,
	     length ? 100 * insertions / length : 0.0,
	     length ? 100 * deletions / length : 0.0,
	     align_entry->m_start - 1,
	     align_entry->reverse ?
	     length1 - align_entry->last_end : align_entry->first_start,
	     get_seq_length(seq_num) - align_entry->m_end,
	     align_entry->reverse ?
	     align_entry->first_start : length1 - align_entry->last_end);
      /*      
	if (align_entry->next) {
	for (pair = align_entry->first_pair; pair; pair = pair->next) {
	if (pair->entry2 == align_entry->next - align_array
	&& is_reverse(pair) ==
	(align_entry->reverse != align_entry->next->reverse))
	break;
	}
	printf(" %3d", pair ? pair->score : 0);
	}
	else printf("    ");    
	*/
      /* following 4 lines are special to artificially constructed reads */
      if (sscanf(get_id(seq_num),"read%d",&offset)) {
	offset1 = offset - (align_entry->reverse ? align_entry->end : align_entry->start);
	offset2 = offset + (align_entry->reverse ? align_entry->end : align_entry->start);
	printf("  %6d %6d", offset1, offset2);
      }
    }
    free_seg_blocks();

    if (n_entries > 1 && t_length)
      printf("\n\nOverall discrep rates (%%):             %.2f %.2f %.2f",
	     100 * t_mismatches / t_length,
	     100 * t_insertions / t_length,
	     100 * t_deletions / t_length);
    
    /*    check_qual(contig1, seq); */
    
    init_qual_hist();
    incr_qual_hist(contig1->adj_qual, contig1->length);
    print_qual_hist("Contig quality");
    print_qual_segments(contig1->adj_qual, contig1->length, parameters->qual_show);
    
    printf("\n\nFirst_start: %d, last_end: %d", contig1->first_start, 
	   contig1->last_end);
    if (contig1->num_entries > 1) {
      analyze_slack(contig1, stdout);
      if (contig1->first_start < 1)
	printf("\nWARNING: %d bases of confirmed read sequence omitted from left end of contig",
	       1 - contig1->first_start);
      if (contig1->last_end > contig1->length)
	printf("\nWARNING: %d bases of confirmed read sequence omitted from right end of contig",
	       contig1->last_end - contig1->length);
      print_gaps(contig1);
      
      read_contig_diff_hist(contig1, 0);
      read_contig_diff_hist(contig1, 1);
      
      qual_flag = (/* signed */ char *)our_alloc((contig1->length + 1) * sizeof(char));
      qual_type = (char *)our_alloc((contig1->length + 1) * sizeof(char));
      for (j = 0; j <= contig1->length; j++) {
	qual_flag[j] = 0;
	qual_type[j] = ' ';
      }
      
      for (strand = 0; strand < 2; strand++) {
	inserts[strand] = (int *)our_alloc((contig1->length + 1) * sizeof(int));
	deletes[strand] = (int *)our_alloc((contig1->length + 1) * sizeof(int));
	subs[strand] = (int *)our_alloc((contig1->length + 1) * sizeof(int));
	depth_trans[strand] = (int *)our_alloc((contig1->length + 2) * sizeof(int));
	
	for (j = 0; j <= contig1->length; j++) 
	  inserts[strand][j] = deletes[strand][j] = subs[strand][j] = depth_trans[strand][j] = 0;
	depth_trans[strand][j] = 0;
      }
      for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) {
	seq_num = align_entry->seq_entry;
	start1 = align_entry->start + align_entry->m_start - 1;
	strand = align_entry->reverse;
	depth_trans[strand][start1] += 1; 	
	depth_trans[strand][start1 + align_entry->m_length] -= 1; 
	if (strand) 
	  length1 = get_seq_length(seq_num);
	site2 = align_entry->m_start - 1;
	site1 = align_entry->start + align_entry->m_start - 2;
	for (diff = align_entry->diffs; *(diff + 1); diff++) {
	  site1 += diff_gap1(*diff);
	  site2 += diff_gap2(*diff);
	  if (diff_type(*diff) == 'M') continue;
	  q2 = align_entry->adj_qual[strand ? length1 - site2 : site2 - 1];
	  if (qual_flag[site1] < q2) {
	    qual_flag[site1] = q2;
	    qual_type[site1] = diff_type(*diff);
	  }
	}
	site2 = align_entry->m_start - 1;
	site1 = align_entry->start + align_entry->m_start - 2;
	for (diff = align_entry->diffs; *(diff + 1); diff++) {
	  site1 += diff_gap1(*diff);
	  site2 += diff_gap2(*diff);
	  if ((diff_t = diff_type(*diff)) == 'M') continue;
	  if (diff_t == 'I') {
	    inserts[strand][site1] += 1;
	    while (diff_type(*(diff + 1)) == 'I' && !diff_gap1(*(diff + 1))) {
	      diff++;
	      site1 += diff_gap1(*diff);
	      site2 += diff_gap2(*diff);
	    }
	  }
	  else if (diff_t == 'D') deletes[strand][site1] += 1;
	  else if (diff_t == 'S') subs[strand][site1] += 1;
	}
      }
      for (strand = 0; strand < 2; strand++) 
	for (j = 1; j <= contig1->length + 1; j++)
	  depth_trans[strand][j] += depth_trans[strand][j - 1];
      printf("\n\nDepth 0 regions:"); 
      first_zer = zer_flag = 0;
      for (j = 1; j <= contig1->length; j++) {
	if (!(depth_trans[0][j] || depth_trans[1][j])) {
	  if (!zer_flag) {
	    zer_flag = 1;
	    first_zer = j;
	  }
	}
	else if (zer_flag) {
	  zer_flag = 0;
	  printf("\n %d - %d", first_zer, j - 1);
	}
      }
      if (zer_flag) printf("\n %d - %d", first_zer, j - 1);
      
      /*	      
	for (start1 = 0; start1 <= contig1->length 
	&& !depth_trans[0][start1] && !depth_trans[1][start1]; start1++);
	if (start1 != 1) {
	printf("\n%d unaligned bases at contig start:\n ", start1 - 1);
	for (j = 0; j < start1 - 1; j++) printf("%c",seq[j]);
	}
	
	for (end1 = contig1->length + 1; end1 >= 1 && 
	!depth_trans[0][end1] && !depth_trans[1][end1]; end1--);
	if (end1 < contig1->length + 1) {
	printf("\n%d unaligned bases at contig end:\n ", contig1->length + 1 - end1);
	for (j = end1 - 1; j < contig1->length; j++) printf("%c",seq[j]);
	}
	*/
      /*
	printf("\n\nMinority calls\nSite  base(qual) #subs #ins #dels depth");
	*/
      count_blocks(contig1->adj_qual, contig1->length);
      
      /*      printf("\n\nLow quality suspects:"); */
      n_ss = n_errs = 0;
      for (j = 1; j <= contig1->length; j++) {
	if (!depth_trans[0][j] || !depth_trans[1][j]) {
	  n_ss++;
	  /*	  continue; */
	}
	/*
	  if (2 * (inserts[strand][j] + deletes[strand][j] + subs[strand][j]) >= depth[strand])
	  */
	/*
	  if (subs[0][j] + inserts[0][j] + deletes[0][j] >= .25 * depth_trans[0][j] 
	  || subs[1][j] + inserts[1][j] + deletes[1][j] >= .25 * depth_trans[1][j]) {
	  if (subs[0][j] + inserts[0][j] + deletes[0][j] > .5 * depth_trans[0][j] 
	  || subs[1][j] + inserts[1][j] + deletes[1][j] > .5 * depth_trans[1][j] || qual_flag[j] != ' ') {
	  n_errs++;
	  printf("\n%c%5d  %c(%d)  %3d  %3d  %3d  %3d     %3d  %3d  %3d  %3d ", 
	  qual_flag[j], j, 
	  seq[j - 1], contig1->adj_qual[j - 1],
	  subs[0][j], inserts[0][j], deletes[0][j], depth_trans[0][j], 
	  subs[1][j], inserts[1][j], deletes[1][j], depth_trans[1][j]);
	  }
	  */
	if (qual_flag[j] && qual_flag[j] < parameters->qual_show
	    && contig1->adj_qual[j - 1] <= qual_flag[j]) {
	  n_errs++;
	  /*
	    printf("\n%5d  %c  %c(%d)  %d", 
	    j, qual_type[j], seq[j - 1], contig1->adj_qual[j - 1], qual_flag[j]);
	    */
	}
      }
      /*      if (!n_errs) printf(" None."); */
      printf("\n\nSS region: %d (%.2f%%), flagged: %d (%.2f%%)",
	     n_ss, 100.0 * n_ss / (float)contig1->length,
	     n_errs, 100.0 * n_errs / (float)contig1->length);
      if (depth_trans[0][j] || depth_trans[1][j]) 
	printf("\n\nERROR in alignment positions:  remaining depth %d %d", 
	       depth_trans[0][j], depth_trans[1][j]);
      for (strand = 0; strand < 2; strand++) {
	our_free(inserts[strand]);
	our_free(deletes[strand]);
	our_free(subs[strand]);
	our_free(depth_trans[strand]);
      }
      our_free(qual_flag);
      our_free(qual_type);

      find_compressions(contig1);
      print_neg_LLR_sites(contig1);
      print_mismatch_reads(contig1, seq);
     }
    
    print_problem_reads(contig1, fp, fq);
    find_unique(contig1); /* may be using memory -- in segment list */
    if (!parameters->print_extraneous_matches) continue;
    printf("\n\nExtraneous matches: ");
    
    /* Simplify the following by making use of is_used(pair) */
    i_ptr0 = i_ptr;
    head1 = head2 = 0;
    
    index1 = contig1->index;
    n_matches = 0;
    n_pos_LLR = max_pos_LLR = n_neg_LLR = min_neg_LLR = 0;
    for (; i_ptr < 2 * num_pairs; i_ptr++) {
      pair = pair_pointers[i_ptr];
      /*       if (!is_best(pair)) continue;   CHECK THIS -- CURRENTLY NEC. BECAUSE
	       WEAK ALIGNMENTS NOT DONE IN GENERAL */
      if (is_reject_self(pair)) continue;
/*      if (is_reject_chimeric(pair)) continue; RE-EXAMINE THIS -- MAY MISS SOME MATCHES */
      align_entry1 = get_align_entry(pair->entry1);
      if ((int)align_entry1->contig->index > index1) break;
      else if ((int)align_entry1->contig->index < index1) fatalError("in extraneous matches");
      align_entry2 = get_align_entry(pair->entry2);
      if (is_anomalous(align_entry2) && align_entry1->contig->num_entries > 1) 
	continue; /* ignore matches to anomalous subclones -- unless 
		     both contigs are singlets */
      length1 = get_seq_length(pair->entry1);
      length2 = get_seq_length(pair->entry2);
      index2 = align_entry2->contig->index;
      rel_orient = is_reverse(pair) !=
	(align_entry1->reverse != align_entry2->reverse);

      unalign1 = is_unaligned(pair);
      unalign2 = is_unaligned(pair->reversed_pair);

      pair_entry1 = pair->entry1;
      pair_entry2 = pair->entry2;

      start1 = unalign1 ? pair->start1 : find_position(align_entry1, pair->start1);
      end1 = unalign1 ? pair->end1 : find_position(align_entry1, pair->end1);
      if (start1 > end1) {
	temp = start1;
	start1 = end1;
	end1 = temp;
      }
      start2 = unalign2 ? (is_reverse(pair) ? length2 - 1 - pair->end2 : pair->start2)
	: find_position(align_entry2, 
			is_reverse(pair) ? length2 - 1 - pair->end2 : pair->start2);
      end2 = unalign2 ? (is_reverse(pair) ? length2 - 1 - pair->start2 : pair->end2)
	: find_position(align_entry2, 
			is_reverse(pair) ? length2 - 1 - pair->start2 : pair->end2);
      if (rel_orient) {
	start2 = -start2;
	end2 = -end2;
      }
      if (start2 > end2) {
	temp = start2;
	start2 = end2;
	end2 = temp;
      }
      /*
	if (align_entry1->reverse) {
	start1 = align_entry1->start + (length1 - 1 - pair->end1);
	end1 = align_entry1->end - pair->start1;
	start2 = (rel_orient ? -align_entry2->end : align_entry2->start)
	+ (length2 - 1 - pair->end2);
	end2 = (rel_orient ? -align_entry2->start : align_entry2->end)
	- pair->start2;
	}
	else {
	start1 = align_entry1->start + pair->start1;
	end1 = align_entry1->end - (length1 - 1 - pair->end1);
	start2 = (rel_orient ? -align_entry2->end : align_entry2->start)
	+ pair->start2;
	end2 = (rel_orient ? -align_entry2->start : align_entry2->end)
	- (length2 - 1 - pair->end2);
	}
	*/
      offset = start1 - start2;

      if (i_ptr == i_ptr0 || index2 != prev_index2 || rel_orient != prev_rel_orient ||
	  unalign1 != prev_unalign1 || unalign1 && pair_entry1 != prev_pair_entry1 ||
	  unalign2 != prev_unalign2 || unalign2 && pair_entry2 != prev_pair_entry2 ||
	  abs(offset - prev_offset) > LOCATION_FUDGE1) {
	/*	  printf("\n%d %d %d",index2,(int)rel_orient,offset); */
	if (head1 || head2) {
	  print_match(&head1, &head2, index1, prev_index2, prev_rel_orient,
		      prev_unalign1, prev_unalign2, prev_pair_entry1, prev_pair_entry2,
		      contig1, &n_pos_LLR, &max_pos_LLR, &n_neg_LLR, &min_neg_LLR); 
	  n_matches++;
	}
	prev_index2 = index2;
	prev_rel_orient = rel_orient;
	prev_unalign1 = unalign1;
	prev_unalign2 = unalign2;
      }
      prev_pair_entry1 = pair_entry1;
      prev_pair_entry2 = pair_entry2;

      prev_offset = offset;
      if (index1 != index2 || !rel_orient && offset <= -LOCATION_FUDGE1 
	  || rel_orient && start1 < -start2 || unalign1 || unalign2) {
	if (outside_match(pair)) {
	  
	  printf("\n*%c  %s %d %d  %c  %d %d  | ", 
		 is_reject_chimeric(pair) ? 'c' : ' ',
		 get_id(align_entry1->seq_entry),
		 pair->start1 + 1, pair->end1 + 1, 
		 align_entry1->reverse ? 'C' : ' ',
		 align_entry1->reverse ? get_seq_length(align_entry1->seq_entry) - align_entry1->m_end + 1: 
		 align_entry1->m_start,
		 align_entry1->reverse ? get_seq_length(align_entry1->seq_entry) - align_entry1->m_start + 1: 
		 align_entry1->m_end);
	  print_segs(align_entry1);
	}
	else if (outside_match(pair->reversed_pair)) {
	  
	  printf("\n *%c %s %d %d  %c  %d %d  | ", 
		 is_reject_chimeric(pair) ? 'c' : ' ',
		 get_id(align_entry2->seq_entry),
		 pair->reversed_pair->start1, pair->reversed_pair->end1, 
		 align_entry2->reverse ? 'C' : ' ',
		 align_entry2->reverse ? get_seq_length(align_entry2->seq_entry) - align_entry2->m_end + 1: 
		 align_entry2->m_start,
		 align_entry2->reverse ? get_seq_length(align_entry2->seq_entry) - align_entry2->m_start + 1: 
		 align_entry2->m_end);
	  print_segs(align_entry2);
	}
	/* show one of two copies only; and not the match to itself */
	/*
	  printf(" %3d %3d %3d %.2f,", offset, pair->end1 - pair->start1 + 1, pair->score,
	  100.0 * (float)(pair->mismatches + 
	  pair->insertions + 
	  pair->deletions)/(float)(pair->end1 - pair->start1 + 1));
	  */
	head1 = insert_segment(head1, start1, end1);
	head2 = insert_segment(head2, start2, end2);
	if (pair->LLR_score >= 0) {
	  n_pos_LLR++;
	  if (pair->LLR_score > max_pos_LLR) max_pos_LLR = pair->LLR_score;
	}
	else {
	  n_neg_LLR++;
	  if (pair->LLR_score < min_neg_LLR) min_neg_LLR = pair->LLR_score;
	}
      }
    }
    if (head1 || head2) {
      print_match(&head1, &head2, index1, index2, rel_orient, 
		  unalign1, unalign2, pair_entry1, pair_entry2,
		  contig1, &n_pos_LLR, &max_pos_LLR, &n_neg_LLR, &min_neg_LLR); 
      n_matches++;
    }
    /* in case loop ends due to i_ptr = 2 * num_pairs */
    if (!n_matches) printf("  None.");
    printf("\n\n");
  }
  if (parameters->print_extraneous_matches) our_free(pair_pointers);
  notify(" Done\n");
}

read_contig_diff_hist(contig, type)
     Contig *contig;
     int type;
{
  Align_info *align_entry;
  int length1;
  unsigned char *seq;
  unsigned char *get_seq(), *get_comp_seq();

  init_diff_hist();
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    seq = align_entry->reverse ? get_comp_seq(align_entry->seq_entry) : get_seq(align_entry->seq_entry);
    length1 = get_seq_length(align_entry->seq_entry) - 1;
      
    set_int_qual(type ? align_entry->adj_qual : align_entry->orig_qual, 
		 length1 + 1, (int)align_entry->reverse);
    incr_diff_hist(seq,
		   align_entry->diffs, length1 + 1, contig->length, 
		   align_entry->m_start - 1, 
		   align_entry->start + align_entry->m_start - 2, 
		   align_entry->m_end - 1, 
		   align_entry->end + align_entry->m_end - length1 - 2, 1, 0);
    free_int_qual();
  }
      
  print_diff_hist(type ? "Read/contig alignment summary, by read base; adjusted qualities" :
"Read/contig alignment summary, by read base; trace qualities");
  free_diff_hist();
}


/* test whether matching part of first read extends significantly past part of
   read that aligns to contig */

print_match(add_head1, add_head2, index1, index2, rel_orient, 
	    unalign1, unalign2, pair_entry1, pair_entry2,
	    contig1, add_n_pos_LLR, add_max_pos_LLR, add_n_neg_LLR, add_min_neg_LLR)
     Segment **add_head1, **add_head2;
     int index1, index2;
     char rel_orient;
     int unalign1, unalign2, pair_entry1, pair_entry2;
     Contig *contig1;
     int *add_n_pos_LLR, *add_max_pos_LLR, *add_n_neg_LLR, *add_min_neg_LLR;
{ 
  printf("\n");
  print_segments(*add_head1, index1, 0); 
  printf(" | ");
  print_segments(*add_head2, index2, rel_orient); 
  if (*add_n_pos_LLR) {
    find_spans(*add_head1, contig1);
    printf("\n %d positive LLR pairs (max LLR score %.1f);  %d negative LLR pairs (min LLR score %.1f).", 
	   *add_n_pos_LLR, *add_max_pos_LLR / 10.0, *add_n_neg_LLR, *add_min_neg_LLR / 10.0);
  }
  printf("\n");
/*  free_segments(*add_head1); */
/*  free_segments(*add_head2); */
  *add_head1 = *add_head2 = 0;
  *add_n_pos_LLR = *add_max_pos_LLR = *add_n_neg_LLR = *add_min_neg_LLR = 0;
}

print_problem_reads(contig, fp, fq)
     Contig *contig;
     FILE *fp, *fq;
{
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Segment *segment, *local_segs, *local_aligned_segs, *local_or_aligned_segs, *aligned_segs;
  Segment *lu_segs, *da_segs, *du_segs;
  Segment *insert_segment(), *find_gap_segment_list();
  int length, start, end, length2, start2, end2;
  int local_flag, site1, site2, gap, start1, end1, slop, lu_flag, da_flag, du_flag, chimera_flag;
  int i, j, hi_q, hi_overlap, lo_overlap;
  int first_i, length1, seq_num, first_q, last_q;
  int n_X;
  char *get_id(), *get_descrip();
  unsigned char *get_seq();
  unsigned char *seq;
  char flag;
  int pre_unaligned_cf, pre_unaligned_hq, post_unaligned_cf, post_unaligned_hq;
  int q3, q4;

  flag = 0;
  slop = 20;

  mark_save_block();

  printf("\n\nReads with neg LLR score, or confirmed or high-qual unaligned seg > 20 bases, or other problem:");
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    seq_num = align_entry->seq_entry;
    length = get_seq_length(seq_num);
    length1 = length - 1;
    seq = get_seq(seq_num);

    first_i = align_entry->reverse ? length1 - align_entry->last_end 
      : align_entry->first_start;
    pre_unaligned_cf = align_entry->m_start - 1 - first_i;
    if (pre_unaligned_cf < 0) pre_unaligned_cf = 0;

    first_i = align_entry->reverse ? length1 - align_entry->first_start 
      : align_entry->last_end;
    post_unaligned_cf = first_i - align_entry->m_end + 1;
    if (post_unaligned_cf < 0) post_unaligned_cf = 0;

    find_trimmed_quals(align_entry, &first_q, &last_q);

    first_i = align_entry->reverse ? length1 - last_q : first_q;
    pre_unaligned_hq = first_q < last_q ? align_entry->m_start - 1 - first_i : 0;
    if (pre_unaligned_hq < 0) pre_unaligned_hq = 0;

    first_i = align_entry->reverse ? length1 - first_q : last_q;
    post_unaligned_hq = first_q < last_q ? first_i - align_entry->m_end + 1 : 0;
    if (post_unaligned_hq < 0) post_unaligned_hq = 0;

    if (align_entry->reverse) {
      start = length - align_entry->m_end;
      end = length - align_entry->m_start;
    }
    else {
      start = align_entry->m_start - 1;
      end = align_entry->m_end - 1;
    }
/*
    for (segment = align_entry->segments; segment; segment = segment->next) 
      if (segment->end > end + slop || segment->start < start - slop) break;
    if (!segment) continue;
*/
    local_segs = local_aligned_segs = aligned_segs = local_or_aligned_segs = 0;
 
/*
    aligned_segs = insert_segment(aligned_segs, start, end);
    local_or_aligned_segs = insert_segment(local_or_aligned_segs, start, end);
    local_segs = insert_segment(local_segs, start, end);
    local_aligned_segs = insert_segment(local_aligned_segs, start, end);
*/
   for (pair = get_aligned_pairs(seq_num); pair; pair = pair->next) {
      if (same_subclone(pair->entry1, pair->entry2)) continue; 
      start1 = pair->start1;
      end1 = pair->end1;
/*
      if (start1 >= start && end1 <= end) continue;
*/

      align_entry2 = get_align_entry(pair->entry2);
      length2 = get_seq_length(align_entry2->seq_entry);
      if (is_used(pair)) local_flag = 1;
      else {
	local_flag = 0;
	if (align_entry->contig->index == align_entry2->contig->index
	    && is_reverse(pair) == (align_entry->reverse != align_entry2->reverse)) {
	  
	  site1 = find_position(align_entry, pair->start1);
	  site2 = find_position(align_entry2, 
				is_reverse(pair) ? length2 - 1 - pair->start2 : pair->start2);
	  if (abs(site1 - site2) <= 30) 
	    local_flag = 1; /* match is local, i.e. compatible with assembly */
	}
      }
      if (local_flag) {
	local_segs = insert_segment(local_segs, start1, end1);
	local_or_aligned_segs = insert_segment(local_or_aligned_segs, start1, end1);
      }
      if (start1 < start) start1 = start;
      if (end1 > end) end1 = end;
      if (start1 < end1) {
	if (local_flag)
	  local_aligned_segs = insert_segment(local_aligned_segs, start1, end1);
	aligned_segs = insert_segment(aligned_segs, start1, end1);
	local_or_aligned_segs = insert_segment(local_or_aligned_segs, start1, end1);
      }
      if (align_entry2->reverse) {
	start2 = length2 - align_entry2->m_end;
	end2 = length2 - align_entry2->m_start;
      }
      else {
	start2 = align_entry2->m_start - 1;
	end2 = align_entry2->m_end - 1;
      }
      if (start2 < pair->reversed_pair->start1) start2 = pair->reversed_pair->start1;
      if (end2 > pair->reversed_pair->end1) end2 = pair->reversed_pair->end1;
      if (start2 < end2) { /* part of matching segment corresponds to contig */
	if (is_reverse(pair)) {
	  start1 = find_entry_position(pair, end2);
	  end1 = find_entry_position(pair, start2);
	}
	else {
	  start1 = find_entry_position(pair, start2);
	  end1 = find_entry_position(pair, end2);
	}
	if (local_flag)
	  local_aligned_segs = insert_segment(local_aligned_segs, start1, end1);
	aligned_segs = insert_segment(aligned_segs, start1, end1);
	local_or_aligned_segs = insert_segment(local_or_aligned_segs, start1, end1);
      }
    }
    gap = find_max_gap_list(local_aligned_segs, align_entry->segments);

    if (gap > slop) {
      
      chimera_flag = 0;
      if (align_entry->segments && align_entry->segments->next
	  && local_segs && !local_segs->next) {
	hi_overlap = lo_overlap = 0;
	for (segment = align_entry->segments; segment; segment = segment->next) {
	  if (segment->start >= local_segs->start && segment->end <= local_segs->end) 
	    hi_overlap++;
	  else if (segment->end < local_segs->start + slop || segment->start > local_segs->end - slop)
	    lo_overlap++;
	}
	if (hi_overlap == 1 && lo_overlap) {
	  chimera_flag = 1;
	  for (segment = align_entry->segments; segment; segment = segment->next) 
	    if (segment->end < local_segs->start + slop 
		|| segment->start > local_segs->end - slop) {
	      append_tag(align_entry, "chimera", segment->start, segment->end);
	    }
	}


      }
      
      lu_flag = 0;
      if (contig->num_entries > 1 && 
	  (!local_aligned_segs || find_max_gap_list(local_aligned_segs, local_segs) 
	  > local_aligned_segs->end - local_aligned_segs->start))
	lu_flag = 1;
      
      if (align_entry->last_vec > 0) { /* mask out vector from reports */
	aligned_segs = insert_segment(aligned_segs, 0, align_entry->last_vec);
	local_or_aligned_segs = insert_segment(local_or_aligned_segs, 0, align_entry->last_vec);
	local_segs = insert_segment(local_segs, 0, align_entry->last_vec);
	local_aligned_segs = insert_segment(local_aligned_segs, 0, align_entry->last_vec);
      }
      
      lu_segs = find_gap_segment_list(local_aligned_segs, local_segs);
      if (!lu_flag) {
	for (segment = lu_segs; segment; segment = segment->next) {
	  if (segment->end + 1 - segment->start > slop) {
	    for (i = segment->start, hi_q = 0; i <= segment->end; i++) 
	      if (align_entry->adj_qual[i] >= 20) hi_q++;
	    if (hi_q > .1 * (segment->end + 1 - segment->start)) {
	      lu_flag = 1;
	      break;
	    }
	  }
	}
      }
      da_flag = find_max_gap_list(local_aligned_segs, aligned_segs) > slop;    
      du_flag = find_max_gap_list(local_or_aligned_segs, align_entry->segments) > slop;    
      /*
	 if (!align_entry->segments->next && !da_flag && !du_flag 
	     && gap < .25 * (align_entry->segments->end - align_entry->segments->start)) 
	     continue;
	     */
    }    
    if ((gap <= slop || !da_flag && !du_flag && !lu_flag)
	&& align_entry->LLR_score >= 0 
	&& pre_unaligned_hq <= 20 && post_unaligned_hq  <= 20
/* 	&& pre_unaligned_cf <= 10 && post_unaligned_cf <= 10 RECONSIDER THIS!! -- BUT
 IT DOES ACCOUNT FOR MOST OF THE NOISE IN THE OUTPUT */
	&& !align_entry->bypassed) continue;

    if (!flag) {
      printf("\n\n(Lngths init HQ, CF unalgn segs); aligned contig pos [LLR score]; (terminal HQ, CF unalgn); read id; aligned read pos");
      printf("\n| confirmed read segments | unaligned read segments matching elsewhere");
      printf("\n (LU = local unaligned, DU = distant unaligned, DA = distant aligned, ** = overlaps trimmed region");
      printf("\n || Best local and distant LLR scores");
/*
      printf("\n ** = overlaps trimmed region; \n (read positions : contig positions : #hiqual bases / total) : ");
*/
      flag = 1;
    }

    printf("\n%s(%d, %d) %5d-%5d [%4.1f] (%d,%d)   %c %-18s %d-%d", 
	   align_entry->bypassed ? "Bypassed: " : "",
	   pre_unaligned_hq, pre_unaligned_cf,
	   align_entry->start + align_entry->m_start - 1, 
	   align_entry->end - (length - align_entry->m_end), 
	   align_entry->LLR_score / 10.0,
	   post_unaligned_hq, post_unaligned_cf,
	   align_entry->reverse ? 'C' : ' ', get_id(seq_num), 
	   align_entry->reverse ? length + 1 - align_entry->m_start 
	     : align_entry->m_start, 
	   align_entry->reverse ? length + 1 - align_entry->m_end 
	     : align_entry->m_end);

    write_entry(fp, get_id(seq_num), 0/*get_descrip(seq_num)*/, get_seq(seq_num));
    fprintf(fq, ">%s", get_id(seq_num));
    for (j = 0; j < length; j++) {
      if (!(j % 50)) fprintf(fq,"\n");
      fprintf(fq,"%d ", (int)align_entry->orig_qual[j]);
    }
    fprintf(fq,"\n");
    
    if (gap > slop && (da_flag || du_flag || lu_flag || chimera_flag)) {

      printf(" |");
      print_segs(align_entry); 
      /*
	 p_flag = m_flag = 0;
	 for (segment = align_entry->segments, i_bit = 1; segment; 
	      segment = segment->next, i_bit *= 2) {
	      printf(align_entry->chimera_bits & i_bit ? "  (%d %d) " : "  %d %d ", 
		     segment->start + 1, segment->end + 1);
		     if (find_max_gap(local_aligned_segs, segment->start, segment->end) <= slop) continue;
		     printf("[");
		     if ((gap = find_max_gap(local_segs, segment->start, segment->end)) > slop) {
		       printf("c");
		       if (gap < segment->end - segment->start + 1) printf(":%d ", gap);
		     }
		     if ((gap = find_max_gap(aligned_segs, segment->start, segment->end)) > slop)  {
		       printf("u");
		       if (gap < segment->end - segment->start + 1) printf(":%d ", gap);
		     }
		     printf("]");
		     }
		     */
      
      /* print chimeric segments that are partly included in the aligned read */
      /* non-chimeric segment that is partly missing */
      /*
	 if (align_entry->chimera_bits & i_bit) {
	   if (segment->end > start + 10 && segment->start < end - 10)
	     p_flag = 1;
	 }
	 else {
	   min_size =  end - start > segment->end - segment->start ?
	     segment->end - segment->start + 1 :  end - start + 1;
	   min_size *= .25;
	
	   min_size = 0; 
	   if (start - segment->start > min_size || segment->end - end  > min_size)
	     m_flag = 1;
	 }
	 */
      printf("|");
      if (lu_flag) {
	print_anomalous_segs("LU", align_entry, lu_segs);
      }
      if (da_flag) {
	da_segs = find_gap_segment_list(local_aligned_segs, aligned_segs);
	print_anomalous_segs("DA", align_entry, da_segs);
      }
      if (du_flag) {
	du_segs = find_gap_segment_list(local_or_aligned_segs, align_entry->segments);
	print_anomalous_segs("DU", align_entry, du_segs);
	for (segment = du_segs; segment; segment = segment->next) {
	  for (pair = get_aligned_pairs(seq_num); pair; pair = pair->next) {
	    if (same_subclone(pair->entry1, pair->entry2)) continue; 
	    if (pair->start1 > segment->end || pair->end1 < segment->start) continue;
	    align_entry2 = get_align_entry(pair->entry2);
	    printf(" [%d %d  with %c %s  %d %d", pair->start1 + 1, pair->end1 + 1,
		   is_reverse(pair) ? 'C' : ' ', get_id(pair->entry2),
		   pair->start2, pair->end2);
	    if (align_entry->contig->index == align_entry2->contig->index
		&& is_reverse(pair) == (align_entry->reverse != align_entry2->reverse)) {
	  
	      length2 = get_seq_length(align_entry2->seq_entry);
	      site1 = find_position(align_entry, pair->start1);
	      site2 = find_position(align_entry2, 
				    is_reverse(pair) ? length2 - 1 - pair->start2 : pair->start2);
	      printf("-- displ. %d", abs(site1 - site2));
	    }
	    printf("] ");
	  }
	}
      }
      if (chimera_flag) printf(" CHIMERIC");
    }
    print_local_distant_scores(align_entry);
    
    if (align_entry->LLR_score < 0) {
      if (align_entry->reverse) {
	q3 = length - 1 - align_entry->qual_end;
	q4 = length - 1 - align_entry->qual_start;
      }
      else {
	q3 = align_entry->qual_start;
	q4 = align_entry->qual_end;
      }
      get_LLR(align_entry->diffs, contig->orig_qual, align_entry->orig_qual, 
	      contig->adj_qual, align_entry->adj_qual, 
	      contig->length, align_entry->start + align_entry->m_start - 2, 
	      align_entry->end - length + align_entry->m_end - 1, 
	      length, align_entry->m_start - 1, align_entry->m_end - 1, 
	      align_entry->reverse, (Segment *)0, align_entry->segments, 
	      1, stdout,
	      0, contig->length - 1, q3, q4, 0);

    }
  }
  if (!flag) printf(" None.");
  free_seg_blocks();
}

print_anomalous_segs(label, align_entry, segments)
     char *label;
     Align_info *align_entry;
     Segment *segments;
{
  Segment *segment;
  int i, hi_q, hi_flag;

  printf(" %s:", label);
  hi_flag = 0;
  for (segment = segments; segment; segment = segment->next) {
    for (i = segment->start, hi_q = 0; i <= segment->end; i++) {
      if (align_entry->adj_qual[i] >= 20) hi_q++;
    }
    hi_flag = align_entry->qual_start < segment->end && align_entry->qual_end > segment->start;

    printf("(%s%d %d%s)", 
	   hi_flag ? "**" : "", segment->start + 1, segment->end + 1, hi_flag ? "**" : ""); 
/*
    printf(" (%s%d %d%s : %d %d : %d / %d)", 
	   hi_flag ? "**" : "", segment->start + 1, segment->end + 1, hi_flag ? "**" : "", 
	   find_position(align_entry, segment->start), find_position(align_entry, segment->end),
	   hi_q, segment->end + 1 - segment->start);
*/
  }
}

unused_pos_match(align_entry, min_LLR)
     Align_info *align_entry;
     int min_LLR;
{
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Align_info *get_align_entry();
 
  for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) 
    if (!is_used(pair) && pair->LLR_score > min_LLR /* && !is_reject_total(pair)  */ 
	&& !is_reject_self(pair) && !is_reject_chimeric(pair)
	&& !(is_anomalous(get_align_entry(pair->entry2)) && align_entry->contig->num_entries > 1)) return 1; 
  /* has matching alignment elsewhere -- but be careful about spuriously unused pairs */
  
  return 0;
}

/* find contig segments that are unique -- i.e. have no non-rejected matches to other reads
   in data set */
/* check that all used pairs are flagged; and make rejection condition more conservative? */
find_unique(contig)
     Contig *contig;
{
  Align_info *align_entry;
  Segment *head;
  Segment *insert_segment();
  int length1, e_start, e_end;
  
  mark_save_block();
  head = 0;
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    if (unused_pos_match(align_entry, 0)) continue;
    length1 = get_seq_length(align_entry->seq_entry) - 1;
    e_start = align_entry->reverse ? length1 - align_entry->last_end : align_entry->first_start;
    e_end = align_entry->reverse ? align_entry->first_start : length1 - align_entry->last_end;
    if (e_start < align_entry->m_start - 1) e_start = align_entry->m_start - 1;
    if (e_end < length1 - align_entry->m_end) e_end = length1 - align_entry->m_end;
    e_start += align_entry->start; 
    e_end = align_entry->end - e_end; 
    if (e_end > e_start)
      head = insert_segment(head, e_start, e_end);
  }
  check_segments(head);
  printf("\n\nGaps in unique-read coverage: ");
  print_segment_gaps(stdout, head, contig->first_start, contig->last_end);
  /*  free_segments(head); */
  free_seg_blocks();
}


/* check qualities in contig (forward/reverse confirmation), adjust contig qual values, and
   print out change.
   Note increases (i.e. bases now confirmed on both strands, which previously weren't
   generally reflect use of lower penalty for aligning reads against contigs; while
   decreases generally reflect fact that some reads do not assemble properly; this
   probably indicates a need to revise pairwise quality calcs. Some decreases may
   indicate problems with alignments (e.g. mononuc runs).
   */

/* N.B. THIS ROUTINE ASSUMES CURRENT DEFAULT QUALITY CONVENTIONS:
   REVERSE CONFIRMED <=> QUALITY >= 4. CHANGE THIS, USING CONTIG->ORIG_QUAL */
check_qual(contig, seq)
     Contig *contig;
     unsigned char *seq;
{
  Align_info *align_entry;
  Segment *insert_segment();
  Segment *for_segs, *rev_segs, *forward, *reverse;
  int last_disc, site1;
  unsigned char *diff;
  int i, j, q, q1, last_ch, first_ch;
  int d_site1, type;
  
  for_segs = rev_segs = 0;
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    last_disc = d_site1 = align_entry->start + align_entry->m_start - 2;
    for (diff = align_entry->diffs; *diff; diff++) {
      d_site1 += diff_gap1(*diff);
      if ((type = diff_type(*diff)) == 'M') continue;
      site1 = type == 'I' ? d_site1 + 1 : d_site1;
      if (site1 > last_disc + parameters->confirm_length) {
	if (align_entry->reverse) rev_segs = 
	  insert_segment(rev_segs, 
			 last_disc + parameters->confirm_trim + 1, 
			 site1 - parameters->confirm_trim - 1);
	else for_segs = 
	  insert_segment(for_segs, 
			 last_disc + parameters->confirm_trim + 1, 
			 site1 - parameters->confirm_trim - 1);
      }
      last_disc = d_site1;
    }
  }
  /*  printf("\n\nQuality increases:"); */
  forward = for_segs;
  reverse = rev_segs;
  last_ch = first_ch = -1;
  for (i = 1; i <= contig->length; i++) {
    while (forward && forward->end < i) forward = forward->next;
    while (reverse && reverse->end < i) reverse = reverse->next;
    q1 = (forward && forward->start <= i && reverse && reverse->start <= i);
    q = contig->adj_qual[i - 1];
    if (q1 /* && q < parameters->rev_match_credit */) {
      /*
	if (i > last_ch + 1) {
	if (last_ch != first_ch) printf("-%d, ", last_ch);
	else if (last_ch != -1) printf(", ");
	else printf("\n");
	printf("%d(%d)", i, q);
	first_ch = i;
	}
	last_ch = i;
	*/
      /* contig->adj_qual[i - 1] +=  parameters->rev_match_credit */;
								  }
  }
  /*
    if (last_ch != first_ch) printf("-%d", last_ch);
    else if (last_ch == -1) printf(" None.");
    */
  printf("\n\nQuality decreases:");
  forward = for_segs;
  reverse = rev_segs;
  last_ch = first_ch = -1;
  for (i = 1; i <= contig->length; i++) {
    while (forward && forward->end < i) forward = forward->next;
    while (reverse && reverse->end < i) reverse = reverse->next;
    q1 = (forward && forward->start <= i && reverse && reverse->start <= i);
    q = contig->adj_qual[i - 1];
    if (!q1 /* && q >= parameters->rev_match_credit */) {
      if (i > last_ch + 1) {
	if (last_ch != first_ch) {
	  printf("-%d", last_ch);
	}
	printf(" ");
	if (last_ch > -1 && last_ch < first_ch + 20)
	  for (j = first_ch > 5 ? first_ch - 5 : 0; 
	       j < last_ch + 4 && j < contig->length; j++) 
	    printf("%c", j < first_ch - 1 || j > last_ch - 1 ? tolower(seq[j]) : seq[j]);
	printf("\n");
	printf("%d(%d)", i, q);
	first_ch = i;
      }
      last_ch = i;
      /* contig->adj_qual[i - 1] -= parameters->rev_match_credit; */
      /* this is specific to current defaults */
    }
  }
  if (last_ch != first_ch) printf("-%d", last_ch);
  printf(" ");
  if (last_ch > -1 && last_ch < first_ch + 20)
    for (j = first_ch > 5 ? first_ch - 5 : 0; 
	 j < last_ch + 4 && j < contig->length; j++) 
      printf("%c", j < first_ch - 1 || j > last_ch - 1 ? tolower(seq[j]) : seq[j]);
  else if (last_ch == -1) printf(" None.");
  /*  free_segments(for_segs); */
  /*  free_segments(rev_segs); */
}  

print_mismatch_reads(contig, seq)
     Contig *contig;
     unsigned char *seq;
{
  Align_info *align_entry;
  char flag;
  int q, max_q, q_c,  site1, site2, n_disc, length1, q_flag, qual;
  unsigned char *diff, *diff2;
  int i, start, middle;
  int t_n_disc, t_n_reads, t_o_disc;
  int d_site1, d_site2;
  int d2_site1, d2_site2;
  int gap1, gap2, d_gap1, d_gap2, type;
  int n_q, n_q_c;
  char *our_alloc();
  int *lower_qual;
  char *discrep_qual;
  double err_q, err_q_c, rel_err;
  int min_seg, min_inset, seq_num;
  unsigned char *get_seq(), *get_comp_seq();
  char *get_id();
  unsigned char *seq1;
  int q_start, q_end;

  flag = 0;
  min_seg = 3; /* 3 minimum size of flanking perfect matching segments */
  min_inset = 1; /* 1 no. of bases to be used in the flanking segments */
  printf("\n\nRead/contig discrepancies (* = higher-quality):");
  t_n_disc = t_n_reads = t_o_disc = 0;
  lower_qual = (int *)our_alloc(contig->length * sizeof(int));
  if (parameters->view)
    discrep_qual = contig->discrep_qual = (char *)our_alloc(contig->length * sizeof(char));
  for (i = 0; i < contig->length; i++) {
    lower_qual[i] = contig->adj_qual[i];
    if (parameters->view) discrep_qual[i] = 0;
  }
  
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    seq_num = align_entry->seq_entry;
    seq1 = align_entry->reverse ? get_comp_seq(seq_num) : get_seq(seq_num);
    length1 = get_seq_length(seq_num) - 1;
    if (align_entry->reverse) {
      q_start = length1 - align_entry->qual_end;
      q_end = length1 - align_entry->qual_start;
    }
    else {
      q_start = align_entry->qual_start;
      q_end = align_entry->qual_end;
    }
    n_disc = 0;
    d_site2 = align_entry->m_start - 2;
    d_site1 = align_entry->start + align_entry->m_start - 3;
    diff = align_entry->diffs - 1;
    for (; *(diff + 1); ) {
      do {
	get_next_diff(&diff, &d_site1, &d_site2, &d_gap1, &d_gap2, &type);
      } while (*(diff + 1) && d_gap2 <= min_seg);
      if (!*(diff + 1)) break;
       
      d2_site2 = d_site2;
      d2_site1 = d_site1; 
      diff2 = diff;
      do {
	get_next_diff(&diff2, &d2_site1, &d2_site2, &gap1, &gap2, &type);
      } while (*(diff2 + 1) && gap2 <= min_seg);
      if (gap2 <= min_seg) break;
      get_prev_diff(&diff2, &d2_site1, &d2_site2, &type); 
      
      for (site2 = d_site2 - min_inset, q_flag = q = max_q = 0; site2 <= d2_site2 + min_inset; site2++) {
	qual = align_entry->adj_qual[align_entry->reverse ? length1 - site2 : site2]; 
	if (max_q < qual) max_q = qual;
	q += qual;
	if (qual >= parameters->qual_show) q_flag = 1;
      } 
      n_q = d2_site2 - d_site2 + 1 + 2 * min_inset;
      
      for (site1 = d_site1 + (d_gap2 - d_gap1) - min_inset, q_c = 0; site1 <= d2_site1 + min_inset; site1++)
	q_c += contig->adj_qual[site1];
      n_q_c = d2_site1 - (d_site1 + (d_gap2 - d_gap1)) + 1 + 2 * min_inset;
      /*       if (q < q_c) continue; favors inserted bases in one or other */
      /*
	
	for (site2 = diff2->site2, site1 = diff2->site1; 
	q == q_c && site2 < (diff2+1)->site2; site2++, site1++){
	q += align_entry->adj_qual[align_entry->reverse ? length1 - site2 : site2]; 
	q_c += contig->adj_qual[site1];
	}
	for (site2 = diff->site2 - 2, site1 = diff->site1 - 2; 
	q == q_c && site2 > (diff2 - 1)->site2; site2--, site1--){
	q += align_entry->adj_qual[align_entry->reverse ? length1 - site2 : site2]; 
	q_c += contig->adj_qual[site1];
	}
	*/      
      err_q = pow(10.0, q * -.1 / n_q);
      err_q_c = pow(10.0, q_c * -.1 / n_q_c);
      if (err_q > .95) err_q = .95;
      if (err_q_c > .95) err_q_c = .95;
      rel_err = err_q_c * (1 - err_q) / (err_q_c * (1 - err_q) + err_q * (1 - err_q_c));
      rel_err = 10.0 * -log10(rel_err);
      if (rel_err < 1) rel_err = 1;
      start = d_site1 + (d_gap2 - d_gap1) - 1;
      middle = (start + d2_site1 + 1) / 2;
      if (parameters->view)
	if (discrep_qual[middle] < max_q) 
	  discrep_qual[middle] = max_q;

      for (site1 = start; site1 <= d2_site1 + 1; site1++) {
	if (rel_err < lower_qual[site1] && d2_site2 >= q_start && d2_site2 <= q_end) lower_qual[site1] = rel_err; 
      }
      if (/* q && */ q * n_q_c >= q_c * n_q) { /* discrepant read has higher avg. quality */
	if (q_flag) n_disc++; /* no longer meaningful */
	else t_o_disc++;
	flag = 1;
	printf("\n%c%5d  %c   %c %-12s   ", q * n_q_c > q_c * n_q /* q_flag */ ? '*' : ' ',
	       d2_site1 + 1, type,	       
	       align_entry->reverse ? 'C' : ' ', get_id(seq_num));
	start = d_site1 + (d_gap2 - d_gap1) - 1;
	printf(" (%d)/(%d) ", q_c, q);
	printf(" %d ", start + 1);
	for (site1 = start; site1 <= d2_site1 + 1; site1++) 
	  printf("%c", seq[site1]);
	printf(" / ");
	for (site2 = d_site2 - 1; site2 <= d2_site2 + 1; site2++) 
	  printf("%c", seq1[site2]);
      }
    }
    /*
      for (diff = align_entry->diffs + 1; diff->type != 'E'; diff++) {
      site = diff->site2 - 1;
      q = align_entry->adj_qual[align_entry->reverse ? length1 - site : site]; 
      q_c = contig->adj_qual[diff->site1 - 1];
      if (q >= parameters->qual_show || q && q >= q_c) {
      if (q >= parameters->qual_show) n_disc++;
      else t_o_disc++;
      flag = 1;
      printf("\n%c%5d  %c  %c(%d)  %c(%d)  %c %-12s   %d  %c", 
      q >= parameters->qual_show ? '*' : ' ',
      diff->site1, diff->type,
      seq[diff->site1 - 1], q_c, db_entry->seq[diff->site2 - 1], q, 
      align_entry->reverse ? 'C' : ' ', db_entry->id, 
      align_entry->reverse ? length1 + 1 - diff->site2 : diff->site2, 
      n_disc > parameters->max_discrep ? '*' : ' ');
      start = diff->site1 - 1 - 1;
      if (start < 0) start = 0;
      end = diff->site1 + 1;
      if (end > contig->length) end = contig->length;
      for (i = start; i < end; i++) 
      contig->adj_qual[i] = 0;
      }
      }
      */
    
    
    if (n_disc) {
      t_n_reads++;
      t_n_disc += n_disc;
    }
  }
  if (!flag) printf(" None.");
  else printf("\n\n%d HQ discrepancies in %d reads.", t_n_disc, t_n_reads);
  printf("\n%d lower quality discrepant sites.", t_o_disc);
  /*
    for (cum_qual = i = 0; cum_qual < 7 && i < contig->length; i++) {
    cum_qual += (contig->adj_qual[i] > 0);
    if (i >= 10) 
    cum_qual -= (contig->adj_qual[i - 10] > 0);
    }
    i = i < 7 ? 0 : i - 7;
    
    for (cum_qual = 0, last_i = contig->length - 1; cum_qual < 7 && last_i > i; last_i--) {
    cum_qual += (contig->adj_qual[last_i] > 0);
    if (last_i < contig->length - 10) 
    cum_qual -= (contig->adj_qual[last_i + 10] > 0);
    }
    last_i = last_i >= contig->length - 7 ? contig->length - 1 : last_i + 7;
    
    for (; i <= last_i; i++) {
    if (lower_qual[i] || !contig->adj_qual[i])
    contig->adj_qual[i] = 3 - lower_qual[i];
    }
    */
/* FOLLOWING SHOULD BE UNNECESSARY -- GIVEN CHANGES TO FIND_EXTENTS 

  for (i = 0; i < contig->length; i++)
    contig->adj_qual[i] = lower_qual[i];
*/

  our_free(lower_qual);
}

/* histogram of no. of blocks of quality <= a given level */
count_blocks(qual, length)
     char *qual;
     int length;
{
  int n_blocks[256], n_bases[256];
  int i, j, cum;
  
  printf("\n\nBlock histogram:\nQual   bases    cum    blocks");
  for (i = 0; i < 256; i++) n_blocks[i] = n_bases[i] = 0;
  n_bases[qual[0]] = 1;
  for (i = 1; i < length; i++) {
    n_bases[qual[i]] += 1;
    for (j = qual[i - 1]; j < qual[i]; j++) n_blocks[j] += 1;
  }
  for (j = qual[i - 1]; j < 256; j++) n_blocks[j] += 1;
  
  for (i = cum = 0; i < 256; i++)
    if (n_bases[i]) {
      cum += n_bases[i]; 
      printf("\n%3d   %5d    %5d    %5d", i, n_bases[i], cum, n_blocks[i]);
    }
  
}


find_compressions(contig)
     Contig *contig;
{
  Align_info *align_entry;
  int **delete_score;
  int **cand_compression_motifs();
  int d_site1, d_site2, length1, type, d_gap1, d_gap2;
  unsigned char *diff;
  unsigned char *seq;
  unsigned char *get_seq(), *get_comp_seq();
  int chemistry, q_start, q_end;

/* SHOULD CHECK FOR DYE PRIMER CHEMISTRY */
  delete_score = cand_compression_motifs(contig->seq, contig->length);
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    chemistry = align_entry->chemistry;
    seq = align_entry->reverse ? 
      get_comp_seq(align_entry->seq_entry) : get_seq(align_entry->seq_entry);
    length1 = get_seq_length(align_entry->seq_entry) - 1;
    q_start = align_entry->reverse ? length1 - align_entry->qual_end : align_entry->qual_start;
    q_end = align_entry->reverse ? length1 - align_entry->qual_start : align_entry->qual_end;
    d_site2 = align_entry->m_start - 2;
    d_site1 = align_entry->start + align_entry->m_start - 3;
    diff = align_entry->diffs - 1;
    for (; *(diff + 1); ) {
      get_next_diff(&diff, &d_site1, &d_site2, &d_gap1, &d_gap2, &type);
      if (chemistry == 0 && d_site1 < contig->length 
	  && delete_score[align_entry->reverse][d_site1] >= 4 
	  && (diff_type(*diff) == 'D' 
	      || diff_type(*diff) == 'S' && contig->seq[d_site1] == 'G'
	         && contig->seq[d_site1 + 1] == 'C' 
	        && seq[d_site2] == 'C' && seq[d_site2 + 1] == 'G' )) {
	append_tag(align_entry, "compression", 
		   align_entry->reverse ? length1 - d_site2 - 1: d_site2,
		   align_entry->reverse ? length1 - d_site2 : d_site2 + 1);
      }
      if (chemistry == 1 && d_site1 < contig->length 
	  && diff_type(*diff) == 'S' && 
	  d_site2 > q_start && d_site2 < q_end &&
	  (seq[d_site2] == contig->seq[d_site1 - 1] || seq[d_site2] == contig->seq[d_site1 + 1]) &&
	      (contig->seq[d_site1] == 'G' && contig->seq[d_site1 - 1] == 'A' && !align_entry->reverse
	       || contig->seq[d_site1] == 'C' && contig->seq[d_site1 + 1] == 'T' && align_entry->reverse)) {
	append_tag(align_entry, "G_dropout", 
		   align_entry->reverse ? length1 - d_site2: d_site2,
		   align_entry->reverse ? length1 - d_site2 :d_site2);
      }
      
    }
  }  
  free_vectors(); /* added per Bonfield suggestions */
}

print_neg_LLR_sites(contig)
     Contig *contig;
{
  Align_info *align_entry;
  char flag;
  int i, d_site1, d_site2, d_loc, q, length1, type, d_gap1, d_gap2;
  int *sites, *max_pos, *max_neg, *top, *bottom;
  char *our_alloc();
  unsigned char *get_seq();
  unsigned char *diff;
  unsigned char *seq;
  
  sites = (int *)our_alloc(contig->length * sizeof(int));
  max_pos = (int *)our_alloc(contig->length * sizeof(int));
  max_neg = (int *)our_alloc(contig->length * sizeof(int));
  top = (int *)our_alloc(contig->length * sizeof(int));
  bottom = (int *)our_alloc(contig->length * sizeof(int));

  for (i = 0; i < contig->length; i++) 
    sites[i] = max_pos[i] = max_neg[i] = top[i] = bottom[i] = 0;
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {

    if (align_entry->LLR_score < 0) continue;

    length1 = get_seq_length(align_entry->seq_entry) - 1;
    seq = get_seq(align_entry->seq_entry);
    d_site2 = align_entry->m_start - 2;
    d_site1 = align_entry->start + align_entry->m_start - 3;
    diff = align_entry->diffs - 1;
    for (; *(diff + 1); ) {
      get_next_diff(&diff, &d_site1, &d_site2, &d_gap1, &d_gap2, &type);
      if (d_site2 <= length1) {
	q = align_entry->orig_qual[align_entry->reverse ? length1 - d_site2 : d_site2]; 
	if (q > 15 && 
	    'X' != seq[align_entry->reverse ? length1 - d_site2 : d_site2]) {
	  d_loc = d_site1 < contig->length ? d_site1 : contig->length - 1;
	  sites[d_loc] += q;
	  if (align_entry->LLR_score < 0) { 
	    if (max_neg[d_loc] < q) max_neg[d_loc] = q;
	  }
	  else { 
	    if (max_pos[d_loc] < q) max_pos[d_loc] = q;
	  }
	  if (align_entry->reverse) bottom[d_loc] += 1;
	  else top[d_loc] += 1;
	}
      }
    }
  }
  
  flag = 0;
  printf("\n\nSites with total LLR scores < -3.0  [max pos LLR read, max neg LLR read]  (#discrep top reads, #discrep bottom reads):");
  for (i = 0; i < contig->length; i++) {
    if (sites[i] >= 30 && (max_pos[i] >= 20 || max_neg[i] >= 20)) {
      flag = 1;
      printf("\n%5d   %6.1f  [%4.1f, %4.1f]  (%d, %d)", 
	     i + 1, -sites[i] / 10.0, -max_pos[i] / 10.0, -max_neg[i] / 10.0, top[i], bottom[i]);
    }
  }
  if (!flag) printf(" None.");
  our_free(sites);
  our_free(max_pos);
  our_free(max_neg);
  our_free(top);
  our_free(bottom);
}


/* return score of best matching read that isn't assembled here */

int match_elsewhere(align_entry)
     Align_info *align_entry;
{    
  Align_info *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int cross_score, first_start, last_end;
  Segment *HQ_head_segment, *LQ_head_segment, *segment, *segment2;
  Segment *weak_insert_segment();
  unsigned char *diff;
  int last_disc, last_disc2, diff_site, diff_site2, d, g1, g2, q_overlap;
  int q_start, q_end, flag, length2;

  cross_score = 0;

  HQ_head_segment = LQ_head_segment = 0;
  for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
    if (is_used(pair) || pair->LLR_score <= 0 || is_reject_self(pair)) continue;

    align_entry2 = get_align_entry(pair->entry2);
    if (is_anomalous(align_entry2)) continue;
    if (align_entry2->contig != align_entry->contig
	|| is_reverse(pair) != (align_entry->reverse != align_entry2->reverse)  
	|| abs(find_pair_offset(pair)) > LOCATION_FUDGE1) {
      length2 = get_seq_length(pair->entry2);
      q_start = align_entry2->m_start - 1;
      q_end = align_entry2->m_end - 1;
      if (align_entry2->reverse) {
	q_start = length2 - 1 - q_end;
	q_end = length2 - 1 - q_start;
      }
      if (q_start < align_entry2->qual_start) 
	q_start = align_entry2->qual_start;
      if (q_end > align_entry2->qual_end) 
	q_end = align_entry2->qual_end;
      if (is_reverse(pair)) {
	q_start = length2 - 1 - q_end;
	q_end = length2 - 1 - q_start;
      }
      q_start += 5; /* to avoid trivial overlaps */
      q_end -= 5;
      if (cross_score < pair->score) cross_score = pair->score;
/*
      q_overlap = is_qual_overlap(pair);
*/
      last_disc = diff_site = pair->start1 - 1;
      last_disc2 = diff_site2 = pair->start2 - 1;
      for (diff = pair->diffs; *diff; diff++) {
	g1 = diff_gap1(*diff);
	g2 = diff_gap2(*diff);
	d = diff_type(*diff);
	diff_site += g1;
	diff_site2 += g2;
	if (d == 'M') continue;
	if (diff_site > last_disc + 15) {
	  if (last_disc2 + 1 <= q_end && diff_site2 - 1 >= q_start)
	    HQ_head_segment = weak_insert_segment(HQ_head_segment, last_disc + 1, diff_site - 1);
	  else
	    LQ_head_segment = weak_insert_segment(LQ_head_segment, last_disc + 1, diff_site - 1);
	}
	last_disc = diff_site;
	last_disc2 = diff_site2;
      }      
/*
      head_segment = insert_segment(head_segment, pair->start1, pair->end1);
*/
    }
  }
  for (segment = HQ_head_segment; segment; segment = segment->next)
    append_tag(align_entry, "matchElsewhereHighQual", segment->start, segment->end);
  for (segment = LQ_head_segment; segment; segment = segment->next) {
    flag = 0;
    for (segment2 = HQ_head_segment; 
	 segment2 && segment2->start <= segment->end; 
	 segment2 = segment2->next) {
      if (segment->start <= segment2->end && segment->end >= segment2->start) {
	flag = 1;
	break;
      }
    }
    if (!flag)
      append_tag(align_entry, "matchElsewhereLowQual", segment->start, segment->end);
  }
  return cross_score;
} 

print_local_distant_scores(align_entry)
     Align_info *align_entry;
{    
  Align_info *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int local_score[2], distant_score[2];
  int i;
  
  for (i = 0; i < 2; i++)
    local_score[i] = distant_score[i] = 0;
  for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
    if (pair->LLR_score <= 0 || is_reject_self(pair)) continue;
    if (is_used(pair)) {
      i = get_align_entry(pair->entry2)->LLR_score > 0;
      if (local_score[i] < pair->LLR_score) 
	  local_score[i] = pair->LLR_score;
      continue;
    }

    align_entry2 = get_align_entry(pair->entry2);
    if (is_anomalous(align_entry2)) continue;
    i = get_align_entry(pair->entry2)->LLR_score > 0;
    if (align_entry2->contig != align_entry->contig
	|| is_reverse(pair) != (align_entry->reverse != align_entry2->reverse)  
	|| abs(find_pair_offset(pair)) > LOCATION_FUDGE1) {
      if (distant_score[i] < pair->LLR_score) 
	distant_score[i] = pair->LLR_score;
    }
    else if (local_score[i] < pair->LLR_score) 
      local_score[i] = pair->LLR_score;
  }
  printf(" || local(+/-) (%.1f,%.1f), distant (%.1f,%.1f)", 
	local_score[1] / 10.0, local_score[0] / 10.0, distant_score[1] / 10.0, distant_score[0] / 10.0);
}

/* need to check whether ->used pairs are still consistent, following final
   contig sequence generation --- */

print_gaps(contig)
     Contig *contig;
{
  Segment *segment;
  Align_info *align_entry, *best_entry;
  int best_start, prev_end;
  char *get_id();
  
  printf("\n\n   DS Gap         Size      Closest read (Start)   Covers   Read length required");
  printf("\n                                                    now?        to cover");
  printf("\nTop strand: \n left - ");
  
  best_entry = 0;
  prev_end = 0;
  
  for (segment = contig->top_segments; segment; segment = segment->next) {
    printf("%5d    %5d%c", 
	   segment->start - 1, segment->start - prev_end - 1, prev_end ? ' ' : '+');
    if (best_entry) {
      printf("      %-12s (%4d)    %s         %5d",
	     get_id(best_entry->seq_entry), best_entry->start, 
	     best_entry->end >= segment->start - 1 ? "Yes" : "No ",
	     segment->start - best_start);
    }
    printf("\n%5d - ", segment->end + 1);
    best_entry = 0;
    for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
      if (align_entry->reverse) continue;
      if (align_entry->start > segment->end) continue; 
      if (!best_entry || align_entry->start > best_start) {
	best_entry = align_entry;
	best_start = align_entry->start;
      }
    }
    prev_end = segment->end;
  }
  printf("right    %5d+", contig->length - prev_end);
  if (best_entry) 
    printf("      %-12s (%4d)    %s         %5d+",
	   get_id(best_entry->seq_entry),  best_entry->start, 
	   /* best_entry->end >= contig->length - 1 ? "Yes" : */ "No ",
	   contig->length - best_start);
  
  printf("\n\nBottom strand: \n left - ");
  prev_end = 0;
  for (segment = contig->bottom_segments; segment; segment = segment->next) {
    printf("%5d    %5d%c", 
	   segment->start - 1, segment->start - prev_end - 1, prev_end ? ' ' : '+');
    best_entry = 0;
    for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
      if (align_entry->reverse && align_entry->end > segment->start && 
	  (!best_entry || align_entry->end < best_start)) {
	best_entry = align_entry;
	best_start = align_entry->end;
      }
    }
    if (best_entry) 
      printf("      %-12s (%4d)    %s         %5d%c",
	     get_id(best_entry->seq_entry),  best_entry->end, 
	     best_entry->start < prev_end ? "Yes" : "No ",
	     best_start - prev_end, prev_end ? ' ' : '+');
    printf("\n%5d - ", segment->end + 1);
    prev_end = segment->end;
  }
  printf("right    %5d+", contig->length - prev_end);
}

/* find contig position (in origin 1) to which site (origin 0 in align_entry)
   corresponds) */

int find_position(align_entry, site)
     Align_info *align_entry;
     int site;
{
  unsigned char *diff;
  int position, site2;
  
  if (align_entry->reverse) 
    site = get_seq_length(align_entry->seq_entry) - 1 - site;
  position = site + align_entry->start;
  site++;
  if (diff = align_entry->diffs) {
    site2 = align_entry->m_start - 1 + diff_gap2(*diff);
    for (; *(diff + 1) && site2 < site; diff++, site2 += diff_gap2(*diff)) {
      if (diff_type(*diff) == 'I') position--;
      else if (diff_type(*diff) == 'D') position++;
    }
  }
  return position;
}   

/* find position (in pair->entry1) to which site (in pair->entry2)
   corresponds); all origin 0 */

int find_entry_position(pair, site)
     Aligned_pair *pair;
     int site;
{
  unsigned char *diff;
  int position, site2;
  
  if (is_reverse(pair)) 
    site = get_seq_length(pair->entry2) - 1 - site;
  position = site + pair->start1 - pair->start2;
  if (diff = pair->diffs) {
    site2 = pair->start2 - 1 + diff_gap2(*diff);
    for (; *(diff + 1) && site2 < site; diff++, site2 += diff_gap2(*diff)) {
      if (diff_type(*diff) == 'I') position--;
      else if (diff_type(*diff) == 'D') position++;
    }
  }
  return position;
}   

find_spans(head, contig)
     Segment *head;
     Contig *contig;
{
  Segment *segment;
  int start, end, e_start, e_end, length1, e_diff, e_max, s_min;
  Align_info *align_entry, *best_entry;
  int best_start, best_end;
  char *get_id();
  
  for (segment = head; segment; segment = segment->next) {
    start = segment->start;
    end = segment->end;
    e_max = -contig->length;
    s_min = contig->length;
    for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
      if (align_entry->start > end) break; /* assumes contig entries are sorted */
      if (align_entry->end < start) continue;
      length1 = get_seq_length(align_entry->seq_entry) - 1;
      e_start = align_entry->reverse ? length1 - align_entry->last_end : align_entry->first_start;
      e_end = align_entry->reverse ? align_entry->first_start : length1 - align_entry->last_end;
      if (e_start < align_entry->m_start - 1) 
	e_start = align_entry->m_start - 1;
      if (e_end < length1 - align_entry->m_end) 
	e_end = length1 - align_entry->m_end;
      e_start += align_entry->start;
      e_end = align_entry->end - e_end;
      e_diff = start - e_start > e_end - end ? e_end - end : start - e_start;
      e_diff -= parameters->maxgap;
      
      if (e_diff > e_max && 
	  (align_entry->reverse ? e_end >= end + parameters->maxgap
	   : e_start <= start - parameters->maxgap)) {
	e_max = e_diff;
	best_entry = align_entry;
	best_start = e_start;
	best_end = e_end;
      }
    }
    printf("\n %d %d %s", start, end, 
	   e_max >= 0 ? " spanned by:" : "  ***UNSPANNED***   Best candidate:");
    if (e_max > -contig->length) 
      printf(" %c %s  %d  %d   (%d)", 
	     best_entry->reverse ? 'C' : ' ',
	     get_id(best_entry->seq_entry), best_start, best_end, e_max);
    else printf(" None.");
  }
}
