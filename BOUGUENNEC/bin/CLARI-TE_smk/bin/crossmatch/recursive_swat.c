/*****************************************************************************
#   Copyright (C) 1994-1998, 2006 by Phil Green.                          
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

/* MAKE TREATMENT OF SHORT READS MORE EFFICIENT -- AVOID RECURSIVE CALLS IF REMAINDER IS TOO SHORT */

extern Parameters *parameters;
extern int num_pairs;

static int num_swat, num_good_scores;
static int n_successes, n_failures;

reset_n_successes()
{
  n_successes = n_failures = 0;
}

reset_swat_statics()
{
  num_swat = num_good_scores = 0;
}

get_num_swat()
{
  return num_swat;
}

get_num_good_scores()
{
  return num_good_scores;
}

print_n_swats()
{
  fprintf(stderr, "Quickalign: %d successes, %d failures\n", n_successes, n_failures);
  fprintf(stderr, "%d SWAT alignments performed. %d pairs have score >= %d\n",
	  num_swat, num_pairs, parameters->minscore);
}

Aligned_pair *sig_pair, *match_pair;

reset_match_pairs()
{
  sig_pair = match_pair = 0;
}

recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
	       q_left, q_right, s_left, s_right)
     int minscore, nextscore;
     Cand_pair *cand_pair;
     Profile *q_profile;
     unsigned char *seq;
     int seq_length;
     int l_edge, r_edge, q_left, q_right, s_left, s_right;
{
  int quick_full_smith_waterman();
  Aligned_pair *append_pair();
  int start1, end1, start2, end2, score;
  int mismatches, deletions, insertions;
  unsigned char *make_diffs();
  int orig_score;
  int success, stat_flag;
  Seq_entry *get_seq_entry(), *seq_entry2;
  Score_hist *append_score_hist();
  Query_domain *query_domain, *update_query_data();

  score = quick_full_smith_waterman(q_profile, seq, seq_length, l_edge, r_edge, 
				    q_left, q_right, s_left, s_right, minscore, &orig_score, &success);

/* NEED TO UPDATE QUERY_DOMAIN->BEST_PAIR BELOW IF APPROPRIATE */
  stat_flag = 0;
  if (parameters->score_flag && (score >= minscore || score >= parameters->min_record_score)) {
    get_stats(&start1, &end1, &start2, &end2, &mismatches, &insertions, &deletions); 
    stat_flag = 1;
    query_domain = update_query_data(cand_pair->entry2, score, cand_pair->reverse, start2 - 1, end2 - 1);

    if (parameters->score_hist)
      append_score_hist(query_domain, score, 0, 1);

    /*
    seq_entry2 = get_seq_entry(cand_pair->entry2);
    if (seq_entry2->score < score) {
      seq_entry2->score = score;
    }
    else if (parameters->masklevel <= 0 || (parameters->minmargin >= 0 && orig_score < seq_entry2->score - parameters->minmargin)) {
      keep_flag = 0;
    }
    if (parameters->score_hist && score >= parameters->min_record_score) {
      append_score_hist(query_domain, score, 0, 1);
    }
    */
  }

  num_swat++;
  if (!(num_swat % 1000)) notify(".");
  if (orig_score >= minscore) {
    if (!stat_flag)
      get_stats(&start1, &end1, &start2, &end2, &mismatches, &insertions, &deletions);

    success ? n_successes++ : n_failures++;

    if (start1 < q_left || end1 > q_right || start1 >= end1 ||
	start2 < s_left || end2 > s_right || start2 >= end2)
      fprintf(stderr, "\nWARNING: %d %d %d start/end out of bounds %d %d %d %d / %d %d %d %d",
	      orig_score, score, minscore, start1, end1, q_left, q_right, start2, end2, s_left, s_right);
    if (score >= minscore &&
      (!parameters->score_flag || parameters->masklevel >= 101 || 
       score >= query_domain->best_score + (parameters->minmargin >= 0 ? 0.0 : parameters->minmargin) 
       && (parameters->minmargin < 1 || query_domain->n_best <= 1)  /* note: incomplete filter in this case -- but need filtering at printout anyway */
       && (parameters->minmargin != 0.5 || !(random() % query_domain->n_best)) )) {


      /* reuse previous pair, if it is isolated */

      if (!match_pair 
	  || sig_pair && match_pair->start1 - match_pair->start2 < sig_pair->start1 - sig_pair->start2 + parameters->max_intron_length
	  || start1 - start2 < match_pair->start1 - match_pair->start2 + parameters->max_intron_length) {

	num_good_scores++;
	num_pairs++; 
	match_pair = append_pair(cand_pair->entry1, cand_pair->entry2, 0);
	set_reverse_flag(match_pair, cand_pair->reverse);
      }
      match_pair->score = score;
      match_pair->diffs = make_diffs();
      match_pair->start1 = start1 - 1;
      match_pair->start2 = start2 - 1;
      match_pair->end1 = end1 - 1;
      match_pair->end2 = end2 - 1;
      if (parameters->score_flag) {
	match_pair->query_data->query_domain = query_domain;
	query_domain->recent_pair = match_pair;
	if (score == query_domain->best_score)
	  query_domain->best_pair = match_pair;
      }
      if (parameters->splice_edge_length) {
	get_edges(match_pair);
      }
      if (score >= parameters->minscore) {
	sig_pair = match_pair;
      }
    } 
    minscore = nextscore; /* next pass done at (potentially) more sensitive score */
/* N.B. FOLLOWING ASSUMES REQUIRE EACH SEQUENCE NEEDS TO BE AT LEAST MINSCORE
   IN LENGTH TO HAVE POSSIBLE SCORE OF MINSCORE -- NOT TRUE WITH MATRIX */

/*    fprintf(stderr, "*%d %d %d %d\n", start1, end1, start2, end2); */
#if defined(OLDVERSION)
    if (l_edge <= r_edge) {
      if (q_left <= start1 - minscore) {
	if (s_left <= start2 - minscore)
	  recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
			 q_left, start1 - 1, s_left, start2 - 1);
	if (s_right >= end2 + minscore)
	  recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
			 q_left, start1 - 1, end2 + 1, s_right);
      }
      if (q_right >= end1 + minscore) {
	if (s_left <= start2 - minscore)
	  recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
			 end1 + 1, q_right, s_left, start2 - 1);
	if (s_right >= end2 + minscore)
	  recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
			 end1 + 1, q_right, end2 + 1, s_right);
      }
    }
    else { /* perhaps this should be the generic case */
      if (s_left < start2)
	recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
		       q_left, q_right, s_left, start2 - 1);
      if (s_right > end2)
	recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
		       q_left, q_right, end2 + 1, s_right);
      if (!parameters->subject_files) {
/* ignore overlapping matches in first input file --- */
	if (q_left < start1) 
	  recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
			 q_left, start1 - 1, s_left, s_right);
	if (q_right > end1) 
	  recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
			 end1 + 1, q_right, s_left, s_right);
      }
    }
#else
    if (s_left < start2)
      recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
		     q_left, q_right, s_left, start2 - 1);
    if (s_right > end2)
      recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
		     q_left, q_right, end2 + 1, s_right);

/* NOTE NON-SYMMETRIC CONDITIONS -- BUT WANT TO AVOID REDUNDANCY */

/*
    if (q_left < start1) 
      recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
		     q_left, start1 - 1, s_left, s_right);
    if (q_right > end1) 
      recursive_swat(minscore, nextscore, cand_pair, q_profile, seq, seq_length, l_edge, r_edge, 
		     end1 + 1, q_right, s_left, s_right);
*/
#endif
  }
}
