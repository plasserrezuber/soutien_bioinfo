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

#define MAX_QUAL 99

#define EXTENT_FUDGE 7 /* adjustment to extent segments to allow for spurious extensions
			  to extents -- also in pairs.c */

#define MIN_QUAL_SEG_SIZE 30 /* minimum size of high-quality region -- if less than this, then
				is reset to size 0 */

extern Parameters *parameters;
extern int t_num_entries;

static double **LLR_discrep[2], **LLR_match[2];
static int **disc_qual_matrix[2], **match_qual_matrix[2]; /* used for computing log_likelihoods of
						     discrepancies, cond'l on qualities */
static int *qual_LLR;
static double *err_rate, *qual_freq;
static double ref_conserve, match_credit;
static int mismatch_LLR[256];
double *scaled_err_probs; /* used for finding high-quality segments */

init_qual_arrays()
{
  int i;
  double x, log_ref, log95;
  char *our_alloc();

  ref_conserve = parameters->repeat_stringency;
  log_ref = log10(1. / ref_conserve);
  log95 = log10(1. / .95); /* this is scaling factor -- so that negative cutoffs retain
			      similar meaning */
  match_credit = 10 * log95; /* log_ref; */
  x = -log(10.0) / 10.0;
/*
  fprintf(stderr, "\nDiscrepancy penalties, by min quality -- repeat_stringency = %.3f: ",
	  ref_conserve);
*/
  for (i = 0; i < 99; i++) {
    mismatch_LLR[i] = (log95 / log_ref) * (-10 * log10(1.0 - ref_conserve + ref_conserve * exp(i * x)) - i);
/*
    fprintf(stderr, "\n%2d  %d ", i, mismatch_LLR[i]);
*/
  }
  for (; i < 256; i++) 
    mismatch_LLR[i] = FORCE_REJECT_SCORE; /* to force marked discrepant reads (having qualities of 99) 
					     not to overlap */

/* following is temporarily inactivated 
  qual_LLR  = (int *)our_alloc((MAX_QUAL + 1) * sizeof(int));
  err_rate  = (double *)our_alloc((MAX_QUAL + 1) * sizeof(double));
  qual_freq  = (double *)our_alloc((MAX_QUAL + 1) * sizeof(double));
  for (k = 0; k < 2; k++) {
    disc_qual_matrix[k] = (int **)our_alloc((MAX_QUAL + 1) * sizeof(int *));
    match_qual_matrix[k] = (int **)our_alloc((MAX_QUAL + 1) * sizeof(int *));
    LLR_discrep[k] = (double **)our_alloc((MAX_QUAL + 1) * sizeof(double *));
    LLR_match[k] = (double **)our_alloc((MAX_QUAL + 1) * sizeof(double *));
    for (i = 0; i <= MAX_QUAL; i++) {
      disc_qual_matrix[k][i] = (int *)our_alloc((MAX_QUAL + 1) * sizeof(int));
      match_qual_matrix[k][i] = (int *)our_alloc((MAX_QUAL + 1) * sizeof(int));
      LLR_discrep[k][i] = (double *)our_alloc((MAX_QUAL + 1) * sizeof(double));
      LLR_match[k][i] = (double *)our_alloc((MAX_QUAL + 1) * sizeof(double));
    }
  }    
*/
} 

set_qual_arrays(db)
     Database *db;
{
  int i, j, k, entry1, reverse, t_n_disc, t_n_match, n_disc, n_match;
  int site1, site2, q1, q2, last_disc1, last_disc2, length2;
  unsigned char *diff;
  char *qual1, *qual2;
  double quot, eij;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int t_pairs[2];
  char *get_adj_qual();

/* check that qualities are in bounds -- this should be a separate routine */
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    qual1 = get_adj_qual(entry1);
    for (i = 0; i < get_seq_length(entry1); i++)
      if (qual1[i] > MAX_QUAL || qual1[i] < 0) {
	fprintf(stderr, "\nQuality out of bounds: %d, max_qual: %d ", 
		qual1[i], MAX_QUAL);
	exit(1);
      }
  }

 /* reference repeat */
  for (i = 0; i <= MAX_QUAL; i++) 
    qual_LLR[i] = qual_freq[i] = err_rate[i] = 0;

  for (k = 0; k < 2; k++) 
    for (i = 0; i <= MAX_QUAL; i++) 
      for (j = 0; j <= MAX_QUAL; j++)
	 disc_qual_matrix[k][i][j] = match_qual_matrix[k][i][j] = 
	   LLR_discrep[k][i][j] = LLR_match[k][i][j] = 0;

  for (entry1 = db->first_entry; entry1 <=  db->last_entry; entry1++) {
/*    if (is_anomalous(align_entry1)) continue; */
    qual1 = get_adj_qual(entry1);
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (pair->entry1 >= pair->entry2) continue;
      if (pair->LLR_score < 0) continue;
/*
      if (pair->score < parameters->confirm_score) continue;
      if (!is_best(pair)) continue;
      if (is_reject_any(pair)) continue;
*/
/*      if (is_anomalous(align_entry2)) continue; */
      qual2 = get_adj_qual(pair->entry2);
      reverse = is_reverse(pair);
      length2 = get_seq_length(pair->entry2) - 1;
      last_disc1 = site1 = pair->start1 - 1;
      last_disc2 = site2 = pair->start2 - 1;
      for (diff = pair->diffs; *diff; diff++) {
	site1 += diff_gap1(*diff);
	site2 += diff_gap2(*diff);
	if (diff_type(*diff) == 'M') continue;
	for (i = last_disc1 + 1, j = last_disc2 + 1; i < site1 || j < site2; i++, j++) {
	  q1 = qual1[i];
	  q2 = qual2[reverse ? length2 - j : j];
	  match_qual_matrix[reverse][q1][q2] += 1;
	  match_qual_matrix[reverse][q2][q1] += 1;
	}
	if (!*(diff + 1)) break;
	last_disc1 = site1;
	last_disc2 = site2;
	q1 = qual1[site1];
	q2 = qual2[reverse ? length2 - site2 : site2];
	disc_qual_matrix[reverse][q1][q2] += 1;
	disc_qual_matrix[reverse][q2][q1] += 1;
      }
    }
  }

  t_n_disc = t_n_match = 0;
  t_pairs[0] = t_pairs[1] = 0;

  for (i = 1; i <= MAX_QUAL; i++) {
    n_disc = n_match = 0;
    for (j = 1; j <= MAX_QUAL; j++) {
      n_disc += disc_qual_matrix[0][i][j] + disc_qual_matrix[1][i][j];
      n_match += match_qual_matrix[0][i][j] + match_qual_matrix[1][i][j];
      t_pairs[0] += disc_qual_matrix[0][i][j] + match_qual_matrix[0][i][j];
      t_pairs[1] += disc_qual_matrix[1][i][j] + match_qual_matrix[1][i][j];
    }
    t_n_disc += n_disc;
    t_n_match += n_match;
    qual_freq[i] = n_disc + n_match;
    err_rate[i] = (n_disc + 1.0) / (n_disc + n_match + 1.0);
  } 

  quot = .5 * (t_n_disc + 1.0) / (t_n_disc + t_n_match + 1.0);
  printf("\n\nQual freq (%%), err_rate (%%), LLR:");
  for (i = 1; i <= MAX_QUAL; i++) {
    if (!qual_freq[i]) continue;
    qual_freq[i] /= t_n_disc + t_n_match;
    err_rate[i] -= quot;
    if (err_rate[i] < .0001) err_rate[i] = .0001;
    qual_LLR[i] = 10 * log10(err_rate[i]);
    printf("\n%2d   %6.2f  %6.2f  %6.1f", i, 100.0 * qual_freq[i], 
	   100.0 * err_rate[i], qual_LLR[i] / 10.0);
  }

  for (k = 0; k < 2; k++)  
    for (i = 1; i <= MAX_QUAL; i++) 
      for (j = 1; j <= MAX_QUAL; j++) {
/* version 950810
	quot = (1.0 + disc_qual_matrix[k][i][j]) / 
	  (1.0 + disc_qual_matrix[k][i][j] + match_qual_matrix[k][i][j]);
	LLR_discrep[k][i][j] = 10 * log10(quot / (ref_conserve * quot + 1 - ref_conserve));
*/
	eij = ref_conserve * (1.0 - err_rate[i]) * (1.0 - err_rate[j]);
	if (disc_qual_matrix[k][i][j] + match_qual_matrix[k][i][j]) {
	  LLR_discrep[k][i][j] = 
	    10 * log10(((1.0 + disc_qual_matrix[k][i][j]) / (1.0 + t_pairs[k]))
		       / (qual_freq[i] * qual_freq[j] * (1 - eij)));
	  if (LLR_discrep[k][i][j] > 0) LLR_discrep[k][i][j] = 0; 
/*
	  LLR_match[k][i][j] = 
	    10 * log10(((1.0 + match_qual_matrix[k][i][j]) / (1.0 + t_pairs[k]))
		       / (qual_freq[i] * qual_freq[j] * eij));
	  if (LLR_match[k][i][j] < 0) LLR_match[k][i][j] = 0; 
*/
	  LLR_match[k][i][j] = match_credit;
	}
      }
 
  printf("\n\nRead/read discrepancies, by quality:");
  printf("\n  ");
  for (i = 0; i <= MAX_QUAL; i++) 
    if (qual_freq[i]) printf("  %5d", i);
  
  for (k = 0; k < 2; k++) {
    printf("\n          %s strand", k ? "Opposite" : "Same");
    for (i = 0; i <= MAX_QUAL; i++) {
      if (qual_freq[i]) {
	printf("\n%2d", i);
	for (j = 0; j <= i; j++) 
	  if (qual_freq[j]) printf("  %5d", disc_qual_matrix[k][i][j]);
      }
    }
  }

  printf("\n\nPer cent:");
  for (k = 0; k < 2; k++) {
    printf("\n          %s strand", k ? "Opposite" : "Same");
    for (i = 0; i <= MAX_QUAL; i++) {
      if (!qual_freq[i]) continue;
      printf("\n%2d", i);
      for (j = 0; j <= i; j++) {
	if (!qual_freq[j]) continue;
	quot = disc_qual_matrix[k][i][j] + match_qual_matrix[k][i][j];
	if (!quot) quot = 1.0;
	printf(" %6.3f", (100.0 * disc_qual_matrix[k][i][j]) / quot);
      }
    }
  }

  printf("\n\nLLR_discrep matrix (for %.2f repeat):", ref_conserve);
  for (k = 0; k < 2; k++) {
    printf("\n          %s strand", k ? "Opposite" : "Same");
    for (i = 0; i <= MAX_QUAL; i++) {
      if (!qual_freq[i]) continue;
      printf("\n%2d", i);
      for (j = 0; j <= i; j++) 
	if (qual_freq[j]) printf(" %6.3f", LLR_discrep[k][i][j] / 10.0);
    }
  }

  printf("\n\nLLR_match matrix (for %.2f repeat):", ref_conserve);
  for (k = 0; k < 2; k++) {
    printf("\n          %s strand", k ? "Opposite" : "Same");
    for (i = 0; i <= MAX_QUAL; i++) {
      if (!qual_freq[i])  continue;
      printf("\n%2d", i);
      for (j = 0; j <= i; j++) 
	if (qual_freq[j]) printf(" %6.3f", LLR_match[k][i][j] / 10.0);
    }
  }
}

/* compute LLR scores; BUT if first_file_only, just make them equal to raw swat score
   (this should be changed, at some point) */

set_LLR_scores()
{
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int entry1, i, cum, score;
  int LLR_hist[400];
  char *o_qual1, *a_qual1;
  char *get_adj_qual(), *get_orig_qual();
  Align_info *get_align_entry();
  Segment *segments1, *segments2;
  char *get_id();

  notify("Computing LLR scores ...");
  for (i = 0; i < 400; i++) LLR_hist[i] = 0;

  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    o_qual1 = get_orig_qual(entry1);
    a_qual1 = get_adj_qual(entry1);
    segments1 = get_align_entry(entry1)->segments;
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (entry1 >= pair->entry2) continue;
      segments2 = get_align_entry(pair->entry2)->segments;
      score = pair->LLR_score = pair->reversed_pair->LLR_score = 
	parameters->subject_files ? pair->score :
	  get_LLR(pair->diffs, o_qual1, get_orig_qual(pair->entry2), 
		  a_qual1, get_adj_qual(pair->entry2), 
		  get_seq_length(entry1), pair->start1, pair->end1, 
		  get_seq_length(pair->entry2), pair->start2, pair->end2, 
		  is_reverse(pair), segments1, segments2, 0, (FILE *)0, 0, 0, 0, 0, 0);
      score = (score / 50) + 300;
      score = score < 0 ? 0 : (score > 399 ? 399 : score);
      LLR_hist[score] += 1;
    }
  }
  printf("\n\nLLR score histogram:\nScore    #   cum # ");
  for (i = cum = 0; i < 400; i++) {
    if (LLR_hist[i]) {
      cum += LLR_hist[i];
      printf("\n%5.1f  %4d  %4d", (i - 300) * 5.0, LLR_hist[i], cum);
    }
  }
  notify(" Done\n");
} 

get_site_LLR(q_val)
     int q_val;
{
  return qual_LLR[q_val];
}

#define MAX_UNALIGNED 49
#define CELL_SIZE 3

get_LLR(diffs, orig_qual1, orig_qual2, adj_qual1, adj_qual2, length1, start1, end1, length2, start2, end2, reverse, segments1, segments2, print_flag, fp, q_start1, q_end1, q_start2, q_end2, ignore_ends)
     unsigned char *diffs;
     char *orig_qual1, *orig_qual2, *adj_qual1, *adj_qual2;
     int length1, start1, end1, length2, start2, end2, reverse;
     Segment *segments1, *segments2;
     int print_flag;
     FILE *fp;
     int q_start1, q_end1, q_start2, q_end2;
     int ignore_ends;
{
  int site1, site2, q1, q2, last_disc1, last_disc2;
  int i, j, unaligned;
  int n_matches;
  double score, score_lt_20, score_gt_20, score_in_HQ, score_out_HQ;
  int score_vector[MAX_UNALIGNED];
  int get_unaligned_LLR();
  unsigned char *diff;
  int d_score;
  Segment *segment;
  int num_lt_20, num_gt_20;

/* following are necessary, since init_qual_arrays no longer used */

  n_matches = score = score_lt_20 = score_gt_20 = score_in_HQ = score_out_HQ = 0;
  num_lt_20 = num_gt_20 = 0;
  length1--;
  length2--;
  last_disc1 = site1 = start1 - 1;
  last_disc2 = site2 = start2 - 1;
  for (diff = diffs; *diff; diff++) {
    site1 += diff_gap1(*diff);
    site2 += diff_gap2(*diff);
    if (diff_type(*diff) == 'M') continue;
    for (i = last_disc1 + 1, j = last_disc2 + 1; i < site1 || j < site2; i++, j++) {
/*
      if ((q1 = adj_qual1[i]) && (q2 = adj_qual2[reverse ? length2 - j : j])) 
	score += LLR_match[reverse][q1][q2]; 
      if ((q1 = orig_qual1[i]) && (q2 = orig_qual2[reverse ? length2 - j : j])) 
	score += match_credit; 
*/
      if (adj_qual1[i] && adj_qual2[reverse ? length2 - j : j]) 
	n_matches++;

    }
						
/* 950810	n_matches++; */
    if (!*(diff + 1)) break;
    last_disc1 = site1;
    last_disc2 = site2;
    q1 = adj_qual1[site1];
    q2 = adj_qual2[reverse ? length2 - site2 : site2];
    if ((d_score = discrep_score(q1, q2)) == FORCE_REJECT_SCORE)
      return FORCE_REJECT_SCORE;

    score += d_score; 
    if (print_flag) {
      if (d_score > -20) {
	score_lt_20 += d_score;
	num_lt_20++;
      }
      else {
	score_gt_20 += d_score;
	num_gt_20++;
      }
      if (site1 >= q_start1 && site1 <= q_end1 && 
	  site2 >= q_start2 && site2 <= q_end2) 
	score_in_HQ += d_score;
      else score_out_HQ += d_score;
    }
/* this is approximation  */
/*
    score += LLR_discrep[reverse][q1][q2];
    q1 = orig_qual1[site1];
    q2 = orig_qual2[reverse ? length2 - site2 : site2];
    score -= q1 < q2 ? q1 : q2; 
*/

  }

  if (print_flag) {
    fprintf(fp,"\nLLR breakdown: discreps: %.1f (<20 part: %.1f (#=%d), >20:%.1f (#=%d); in HQ: %.1f, out HQ %.1f), match: %.1f ", 
	    score/10.0, score_lt_20/10.0, num_lt_20, score_gt_20/10.0, num_gt_20, score_in_HQ/10.0, score_out_HQ/10.0, (double)n_matches * match_credit/10.0);
  }
  score += n_matches * match_credit;

/* 950810	  score += n_matches * 10 * log_ref; */
/* find minimum possible penalty for unaligned parts of reads */
  unaligned = MAX_UNALIGNED; /* look at maximum of 33 bases at each end (8 "cells") */
  if (unaligned > length1 - end1) unaligned = length1 - end1;
/* adjust for chimeric reads */
  if (segments1 && segments1->next) {
    for (segment = segments1->next; segment; segment = segment->next) {
      if (segment->start + EXTENT_FUDGE > end1) {
	if (unaligned > segment->start - end1) 
	  unaligned = segment->start - end1;
      }
    }
  }
  if (unaligned > length2 - end2) unaligned = length2 - end2;
  if (segments2 && segments2->next) {
    if (reverse) {
      for (segment = segments2; segment->next; segment = segment->next) {
	if ((length2 - segment->end) + EXTENT_FUDGE > end2) {
	  if (unaligned > (length2 - segment->end) - end2) 
	    unaligned = (length2 - segment->end) - end2;
	}
      }
    }
    else {
      for (segment = segments2->next; segment; segment = segment->next) {
	if (segment->start + EXTENT_FUDGE > end2) {
	  if (unaligned > segment->start - end2) 
	    unaligned = segment->start - end2;
	}
      }
    }
  }

  if (unaligned > 0) {
    for (i = 0; i < unaligned; i++) {
      q1 = adj_qual1[i + end1 + 1];
      site2 = i + end2 + 1;
      q2 = adj_qual2[reverse ? length2 - site2 : site2];
      score_vector[i] = discrep_score(q1, q2);
/*
      score_vector[i] = LLR_discrep[reverse][q1][q2];
      q1 = orig_qual1[i + end1 + 1];
      q2 = orig_qual2[reverse ? length2 - site2 : site2];
      score_vector[i] -= q1 < q2 ? q1 : q2; 
*/
    }
    if ((d_score = get_unaligned_LLR(score_vector, unaligned)) == FORCE_REJECT_SCORE)
      return FORCE_REJECT_SCORE;
    if (print_flag) 
      fprintf(fp," trail: %.1f ", d_score/10.0);
    if (!ignore_ends) score += d_score;
  }

  unaligned = MAX_UNALIGNED;

  if (unaligned > start1) unaligned = start1;
  if (segments1 && segments1->next) {
    for (segment = segments1; segment->next; segment = segment->next) {
      if (segment->end - EXTENT_FUDGE < start1) {
	if (unaligned > start1 - segment->end) 
	  unaligned = start1 - segment->end;
      }
    }
  }
  if (unaligned > start2) unaligned = start2;
  if (segments2 && segments2->next) {
    if (reverse) {
      for (segment = segments2->next; segment; segment = segment->next) {
	if ((length2 - segment->start) - EXTENT_FUDGE < start2) {
	  if (unaligned > start2 - (length2 - segment->start)) 
	    unaligned = start2 - (length2 - segment->start);
	}
      }
    }
    else {
      for (segment = segments2; segment->next; segment = segment->next) {
	if (segment->end - EXTENT_FUDGE < start2) {
	  if (unaligned > start2 - segment->end) 
	    unaligned = start2 - segment->end;
	}
      }
    }
  }



  if (unaligned > 0) {
    for (i = 0; i < unaligned; i++) {
      q1 = adj_qual1[start1 - 1 - i];
      site2 = start2 - 1 - i;
      q2 = adj_qual2[reverse ? length2 - site2 : site2];
      score_vector[i] = discrep_score(q1, q2);
/*
      score_vector[i] = LLR_discrep[reverse][q1][q2];
      q1 = orig_qual1[start1 - 1 - i];
      q2 = orig_qual2[reverse ? length2 - site2 : site2];
      score_vector[i] -= q1 < q2 ? q1 : q2; 
*/
    }
    if ((d_score = get_unaligned_LLR(score_vector, unaligned)) == FORCE_REJECT_SCORE)
      return FORCE_REJECT_SCORE;
    if (!ignore_ends) score += d_score;
    if (print_flag) 
      fprintf(fp," lead: %.1f ", d_score/10.0);
  }
  if (print_flag) 
    fprintf(fp," total: %.1f ", score/10.0);
  return (int)score < FORCE_REJECT_SCORE + 1 ? FORCE_REJECT_SCORE + 1 : (int)score;
}

int discrep_score(q1, q2)
     int q1, q2;
{
  q1 = q1 < q2 ? q1 : q2; /* this is approximate! Should actually be slightly smaller */

  return mismatch_LLR[q1];
}


/*  get any negative score in unaligned region (keeping in mind that there must
    be enough mismatches there to make score non-positive; and assuming those
    mismatches are in most favorable positions with respect to LLR score */
get_unaligned_LLR(score_vector, unaligned)
     int *score_vector;
     int unaligned;
{
  int n_cells, i_cell, i;
  int score, temp;

  score = score_vector[0];
  n_cells = (unaligned - 1) / CELL_SIZE; /* here cell size = mismatch penalty + 1 (now avg. of indel and 
				    mismatch penalty is no. of positions
				    required to include one mismatch + compensating matches, 
				    assuming a match reward of 1. */
  for (i_cell = 1; i_cell <= n_cells; i_cell++) {
    /* since maximum score is 0, if attain it don't need to go any further */
/* do not allow marked regions to be displaced */
/* following removed 8/5/98 -- causes incorrect rejection of some merges.
   should be unnec with new consed.
    if (score_vector[i_cell] == FORCE_REJECT_SCORE) 
      return FORCE_REJECT_SCORE;
*/
    for (i = i_cell + 1; score_vector[i_cell] && i <= i_cell * CELL_SIZE; i++) 
      if (score_vector[i_cell] < (temp = score_vector[i])) { 
	score_vector[i] = score_vector[i_cell];
	score_vector[i_cell] = temp;
      }

    if (score_vector[i_cell] == FORCE_REJECT_SCORE) 
      return FORCE_REJECT_SCORE;
    score += score_vector[i_cell];
  }
  return score < FORCE_REJECT_SCORE + 1 ? FORCE_REJECT_SCORE + 1 : score;
}

set_scaled_err_probs()
{
  double trim_err, x;
  int i;
  char *our_alloc();

  trim_err = pow(10.0, -parameters->trim_qual / 10.0);
  scaled_err_probs = (double *)our_alloc(256 * sizeof(double));
  for (i = 0; i < 256; i++) {
    x = i < 100 ? pow(10.0, -i / 10.0) : 0.0;
    scaled_err_probs[i] = trim_err - x;
  }
  
}

/* read in quality data */
read_qual(db)
     Database *db;
{
  FILE *fp;
  FILE *fopen();
  File *file;
  int i, j, k, length, entry1, c;
  char *our_alloc();
  char *qf_name, *orig_qual, *id, *idq;
  unsigned char *seq;
  unsigned char *get_seq();
  char *get_orig_qual();
  char *get_id();
  Align_info *align_entry;
  Align_info *get_align_entry();
  long int n_high;
/*
  char dig_buffer[1000];
*/
  int i_string;

  file = db->file; 
  parameters->qual_flag = 1;

  init_qual_hist();

  if (file->type) { /* CALF file -- qualities already known */
    for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++)
      incr_qual_hist(get_orig_qual(entry1), get_seq_length(entry1));

    print_qual_hist("Input quality");
    return; 	      
  }
  notify("Reading quality files ...");

  entry1 = db->first_entry;
  qf_name = (char *)our_alloc((strlen(file->name) + 6) * sizeof(char));
  strcpy(qf_name, file->name);
  strcat(qf_name, ".qual");
  fp = fopen(qf_name, "r");
  if (fp) {
    n_high = 0;
    printf("\n\nQuality file: %s", qf_name);
    fgetc(fp);
    do {
      id = get_id(entry1);
      for (c = fgetc(fp), idq = id; c == *idq; 
	   c = fgetc(fp), idq++);
      if (!isspace(c) && c != '.') { /* 2d condition is to allow for trimming of .seq */
	fprintf(stderr,"\nFile %s, entry %s", qf_name, id);
	fatalError("Inconsistency between quality file and sequence file");
      }
      for (; c != '\n' && c != EOF; c = fgetc(fp));
      
      seq = get_seq(entry1);
      length = get_seq_length(entry1); 
      
      orig_qual = get_orig_qual(entry1);
      align_entry = get_align_entry(entry1);

      for (j = 0; j < length; j++) {
	for (; isspace(c); c = fgetc(fp));
	for (i_string = k = 0; isdigit(c); c = fgetc(fp), i_string++) 
	  k = 10 * k + (c - '0'); /* assumes ASCII convention */
	
	if (!i_string || k > MAX_QUAL) {
	  if (!i_string && !isspace(c))
	    fprintf(stderr, 
		    "Error: ASCII character %d found where digit expected in %s, entry %s: position %d out of %d", 
		    c, qf_name, get_id(entry1), j + 1, length);
	  else 
	    fprintf(stderr, 
		  "Error or quality value out of bounds (< 0 or > %d) in file %s, entry %s", 
		  MAX_QUAL, qf_name, get_id(entry1));
	  fatalError("");
	}
	if (parameters->DNA_flag == 3) { /* quality is stored in same bytes as seq */
	  if (seq[j] & 63) { /* non ACGT chars will get qual 0 */
	    if (k > 60) {
	      n_high++;
	      k = 60;
	    }
	    seq[j] = (seq[j] & 192) + k + 1;
	  }
	}
	else {
	  orig_qual[j] = k != 98 ? k : 0; /* change to 0 for 'N' ?? */
	}
      }
      
      while (isspace(c = fgetc(fp)));
      if (c != '>' && c != EOF) {
	fprintf(stderr,"\nFile %s, entry %s\n%c\n", qf_name, id, c);
	fatalError("Error or file inconsistency: quality and sequence file");
      }
      incr_qual_hist(orig_qual, length);

/* note if qualities are rescaled, this may not be what is desired */
/* Moved to set_qual_segs, to deal with case where no input qualities are provided 
     best_qual_seg(orig_qual, length, &(align_entry->qual_start), 
		    &(align_entry->qual_end));
*/
      entry1++;
    } while (c != EOF); 
    fclose(fp);
    print_qual_hist("Input quality");
    if (parameters->DNA_flag == 3) {
      printf("\n(%ld quals > 60 were truncated to 60)", n_high);
    }
  }
  else { 
    parameters->qual_flag = 0;
    printf("\n\nNO QUALITY FILE %s WAS FOUND. REMAINING INPUT QUALITIES SET TO %d.", 
	   qf_name, parameters->default_qual);
    fprintf(stderr, "\n\nNO QUALITY FILE %s WAS FOUND. REMAINING INPUT QUALITIES SET TO %d.", 
	    qf_name, parameters->default_qual);
  }
/* rescale quality if necessary (done only if 0 < max_qual <= 20; rescaled values have
 maximum of 30 */
/*
  rescale_qual(db);
*/
  our_free(qf_name);
  notify(" Done\n");    
  if (parameters->DNA_flag == 3) { 
    /* assuming that complements are stored!! */
    store_complements(db); /* should really be done outside */ 
  }
}

set_qual_segs(db)
     Database *db;
{
  int entry1;
  char *get_orig_qual();
  Align_info *get_align_entry();
  Align_info *align_entry;

  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    best_qual_seg(get_orig_qual(entry1), get_seq_length(entry1), &(align_entry->qual_start), 
		    &(align_entry->qual_end));
  }
}

/* NOT CORRECT WITH COMPRESSED QUAL */
trim_qual(db, window, min_pos)
     Database *db;
  int window, min_pos;
{
  int entry1, length, cum_qual, j, k, j_min;
  char *qual;
  char *get_orig_qual();

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage1");

  init_qual_hist();
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    qual = get_orig_qual(entry1);
    length = get_seq_length(entry1);
    cum_qual = 0;
    for (j = 0; j < window; j++) cum_qual += (qual[j] > 0);
    for (; cum_qual < min_pos && j < length; j++) 
      cum_qual += ((qual[j] > 0) - (qual[j - window] > 0));
    for (k = 0; k < j - window; k++) qual[k] = 0; /* NOT CORRECT WITH COMPRESSED QUAL */

    j_min = j;
    cum_qual = 0;
    for (j = length - 1; j >= length - window; j--) cum_qual += (qual[j] > 0);
    for (; cum_qual < min_pos && j > j_min; j--)
      cum_qual += (qual[j] > 0) - (qual[j + window] > 0);
    for (k = length - 1; k > j + window; k--) qual[k] = 0;  /* NOT CORRECT WITH COMPRESSED QUAL */
    incr_qual_hist(qual, length);
  }
  print_qual_hist("Trimmed quality");
}

/* NOT CORRECT WITH COMPRESSED QUAL */
trim_quals(db, qual_type)
     Database *db;
     int qual_type;
{
  Align_info *get_align_entry();
  char *get_adj_qual(), *get_orig_qual(), *get_id();
  Align_info *align_entry;
  char *qual;
  int i, length1, start, entry1;
  int cum, max_cum, best_start, best_end;

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage2");
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    length1 = get_seq_length(entry1);
    qual = qual_type ? get_adj_qual(entry1) : get_orig_qual(entry1);

    if (!qual_type) {
      best_start = align_entry->qual_start;
      best_end = align_entry->qual_end;
    }
    else
      best_qual_seg(qual, length1, &best_start, &best_end);

    if (best_start > align_entry->first_start) best_start = align_entry->first_start;
    if (best_end < align_entry->last_end) best_end = align_entry->last_end;

    for (i = 0; i < best_start; i++) qual[i] = 0;
    for (i = best_end + 1; i < length1; i++) qual[i] = 0;

  }
}
 
/* NOT CORRECT WITH COMPRESSED QUAL */
trim_X_quals(db)
     Database *db;
{
  char *get_adj_qual(), *get_orig_qual(), *get_id();
  unsigned char *get_seq();
  char *qual;
  unsigned char *seq;
  int i, length1, start, entry1;
  int start1, end1, start2, end2;

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage3");
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    qual = get_orig_qual(entry1);
    seq = get_seq(entry1);
    for (i = 0; seq[i] && seq[i] != 'X'; i++);
    start1 = i;
    while (seq[i] == 'X') i++;
    end1 = i;
    for (; seq[i] && seq[i] != 'X'; i++);
    start2 = i;
    while (seq[i] == 'X') i++;
    end2 = i;
    if (start1 > parameters->vector_bound) start1 += 15;
    for (i = start1; i < end1 - 15; i++) qual[i] = 0;
    for (i = start2 + 15; i < end2; i++) qual[i] = 0;
  }
}

/* NOT CORRECT WITH COMPRESSED QUAL */
best_qual_seg(qual_array, length, add_start, add_end)
     char *qual_array;
     int length;
     int *add_start, *add_end;
{
  int i, start, q_start, q_end, q, qual_mask;
  double max_cum, cum;

  qual_mask = parameters->DNA_flag >= 3 ? 63 : 127;
  q_start = length;
  q_end = -1;
  max_cum = cum = 0;
  start = 0; 
  for (i = 0; i < length; i++) {
    cum += scaled_err_probs[qual_mask & qual_array[i]];      /* qual_array[i] - trim_level; */
    if (cum > max_cum) {
      max_cum = cum;
      q_end = i;
      q_start = start;
    }
    else if (cum < 0) {
      cum = 0; 
      start = i + 1;
    }
  }
  *add_start = q_end > q_start + MIN_QUAL_SEG_SIZE ? q_start : length;
  *add_end = q_end > q_start + MIN_QUAL_SEG_SIZE ? q_end : -1;
}

/* NOT CORRECT WITH COMPRESSED QUAL */
find_trimmed_quals(align_entry, add_qual_start, add_qual_end)
     Align_info *align_entry;
     int *add_qual_start, *add_qual_end;
{
  int n_x, qual_start, qual_end;
  unsigned char *seq;
  unsigned char *get_seq();

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage5");
  qual_start = align_entry->qual_start;
  qual_end = align_entry->qual_end;
  seq = get_seq(align_entry->seq_entry);

  for (n_x = 0; qual_end >= 0 && n_x < 10; qual_end--) 
    n_x += seq[qual_end] != 'X';
  if (qual_end >= 0) qual_end += 10;
  
  if (qual_start < align_entry->last_vec + 1) 
    qual_start = align_entry->last_vec + 1;

  *add_qual_end = qual_end;
  *add_qual_start = qual_start;
}

/* NOT CORRECT WITH COMPRESSED QUAL */
rescale_qual(db)
     Database *db;
{
  SEQ_AREA i;
  int max_qual;
  double scale_param;

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage6");
  max_qual = get_max_qual();
  if (!max_qual || max_qual > 20) return;
  scale_param = 30.0 / max_qual;

  for (i = 0; i < db->t_length; i++) 
    db->orig_qual_area[i] *= scale_param;

}

/* sets adj_qual equal to orig_qual */
init_adj_qual(db)
     Database *db;
{
  SEQ_AREA i;

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage7");
  for (i = 0; i < db->t_length; i++) 
    db->adj_qual_area[i] = db->orig_qual_area[i];
}


/* convert (mostly) mononucleotide, quality < qual_cut, ends of sequence to a single char 
(usually 'N'); also converts trim_start letters at beginning to 'N' */

convert_ends(db) 
     Database *db;
{   
  int min_score, penalty, qual_cut;
  int trim_start, j, max, cum, max_pos;
  char c, end, end_char;
  char test_chars[30];
  char *chars;
  unsigned char *seq;
  int entry1, length, qual_start, qual_end;
  unsigned char *get_seq();
  char *get_id();
  Align_info *align_entry;
  Align_info *get_align_entry();
  
  notify("Testing ends ...");
  min_score = parameters->trim_score;
  penalty = parameters->trim_penalty;
  trim_start = parameters->trim_start;
  end_char = 'N';
  qual_cut = 20;
  strcpy(test_chars,"ACGT");

  printf("\n\nFollowing regions converted to %c's\n", end_char);

  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    qual_start = align_entry->qual_start;
    qual_end = align_entry->qual_end;
/*
    qual = get_orig_qual(entry1);
*/
    seq = get_seq(entry1);
    length = get_seq_length(entry1);
    for (end = 0; end <= 1; end++) {
      for (chars = test_chars; *chars; chars++) {
	cum = max = 0;
	for (j = end ? 0 : length - 1; 
/*
	     (end ? j < length : j >= 0) && qual[j] <= qual_cut; 
*/
	     end ? j < qual_start : j > qual_end;
	     j += (end ? 1 : -1)) {
	  c = seq[j];
	  cum += (c == *chars) ? 1 : (c == 'N' || c == 'X' ? 0 : penalty);
	  if (cum >= max) {
	    max = cum;
	    max_pos = j;
	  }
	}
	if (max >= min_score) {
	  printf("\n%s %d-%d ", get_id(entry1), end ? 1 : max_pos + 1, 
		 end ? max_pos + 1 : length);
	  for (j = end ? 0 : max_pos; j <= (end ? max_pos : length - 1); j++) {
	    printf("%c", seq[j]);
	    seq[j] = end_char;
	  }
	}      
      } 
    }
    for (j = 0; j < trim_start && seq[j]; j++) seq[j] = end_char;
  }     
  notify(" Done\n");
}

transform_contig_quality(contig1)
  Contig *contig1;
{
  int j;

  for (j = 0; j < contig1->length; j++)
    contig1->adj_qual[j] = -qual_LLR[contig1->adj_qual[j]];
}

print_quality(label, db)
     char *label;
     Database *db;
{
  int entry1;
  char *get_adj_qual();

  notify("Printing qualities ...");
  init_qual_hist();
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    incr_qual_hist(get_adj_qual(entry1), get_seq_length(entry1));
  }
  print_qual_hist(label);
  notify(" Done\n");
}

print_qual_segments(qual, length, cut)
     /* signed */ char *qual;
     int length, cut;
{
  int i, j, n_segs, t_size;

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage8");
  for (i = 0; i < length && !qual[i]; i++);
  for (j = length - 1; j > i && !qual[j]; j--);
  printf("\n\nInitial, terminal qual 0 segments: ");
  if (i) printf(" 1-%d,", i);
  else printf(" (None),");
  if (j < length - 1) printf(" %d-%d", j + 2, length);
  else printf(" (None)");

  printf("\n\nRegions of LLR- adjusted quality < %.1f:\n", cut / 10.0);
  for (i = n_segs = t_size = 0; i < length; i = j) {
    for (; i < length && (int)qual[i] >= cut; i++);
    for (j = i; j < length && (int)qual[j] < cut; j++);
    if (j == i + 1) printf("%d, ", j);
    else if (j > i) printf("%d-%d, ", i + 1, j);
    t_size += j - i;
    n_segs++;
    if (!(n_segs % 8)) printf("\n");
  } 
  if (!n_segs) printf(" (None).");
  else printf("\n\n%d regions, avg size %.1f, avg spacing %.1f",
	      n_segs, t_size / (float)n_segs, length / (float)n_segs);
}

static int qual_hist_c[257];
static int *qual_hist;
static int num_entries;

init_qual_hist()
{
  int i;

  qual_hist = qual_hist_c + 1;
  num_entries = 0;
  for (i = -1; i < 256; i++) qual_hist[i] = 0;
}

incr_qual_hist(qual, length)
     /* signed */ char *qual;
     int length;
{
  int i, q;

  num_entries++;
  if (parameters->compact_qual) {
    for (i = 0; i < length && (qual[i] & 63) <= 1; i++) qual_hist[-1] += 1;
    for (length--; (qual[length] & 63) <= 1 && length > i; length--) qual_hist[-1] += 1;
    for (; i <= length; i++) {
      q = qual[i] & 63;
      qual_hist[q ? q - 1 : 0] += 1;
    }
  }
  else {
    for (i = 0; i < length && !qual[i]; i++) qual_hist[-1] += 1;
    for (length--; !qual[length] && length > i; length--) qual_hist[-1] += 1;
    for (; i <= length; i++) qual_hist[qual[i]] += 1;
  }
}

get_max_qual()
{
  int i;

  for (i = 255; i >= -1; i--) if (qual_hist[i]) return i;
  return 0;
}

print_qual_hist(label)
     char *label;
{
  int i, cum;
  double mean, total, div, cum_errs;

  printf("\n\n%s (quality, n_residues, %%,  cum, cum %%,  cum expected errs):", label);
  total = qual_hist[-1] + qual_hist[0];
  mean = 0;
  for (i = 1; i < 256; i++) {
    if (qual_hist[i]) {
      total += qual_hist[i];
      mean += (double)i * qual_hist[i];
    }
  }
  div = total ? total : 1;
  cum_errs = 0;
  for (i = 255, cum = 0; i >= -1; i--) 
    if (qual_hist[i]) {
      cum += qual_hist[i];
      cum_errs += qual_hist[i] * (i > 0 ? pow(10.0, i * -.1) : 1);
      printf("\n%3d  %6d %5.1f  %6d %5.1f  %6.2f", 
	     i, qual_hist[i], 100.0 * qual_hist[i] / div, cum, 100.0 * cum / div, cum_errs);
    }
  printf("   (quality -1 = terminal quality 0)");
  printf("\n\nAvg. full length: %.1f, trimmed (qual > -1): %.1f\nAvg. quality: %.1f per base", 
	 total / num_entries, (total - qual_hist[-1]) / num_entries, mean / div); 
}

static int *qual;

set_int_qual(c_qual, length, reverse)
     /* signed */ char *c_qual;
     int length, reverse;
{
  int j, q;
  char *our_alloc();

  qual = (int *)our_alloc(length * sizeof(int));

  if (parameters->compact_qual) {
    for (j = 0; j < length; j++) {
      q = c_qual[reverse ? length - 1 - j : j] & 63;
      qual[j] = q ? q - 1 : 0;
      if (qual[j] > MAX_QUAL) fatalError("Quality out of bounds");
    }
  }
  else {
    for (j = 0; j < length; j++) {
      qual[j] = c_qual[reverse ? length - 1 - j : j];
      if (qual[j] > MAX_QUAL) fatalError("Quality out of bounds");
    }
  }
  for (j = 0; j < length && !qual[j]; j++) qual[j] = -1; 
  for (j = length - 1; !qual[j]; j--) qual[j] = -1;
}

free_int_qual()
{
  if (qual) our_free(qual);
}
 

/* structures for recording difference/confirmation info for subject sequences, 
   in resequencing project
*/

typedef struct diffdata {
  struct diffdata *next;
  unsigned char *seq; /* not needed for query del */
  int length, count[2], max_score, max_margin;
  char type, max_qual; /* type: 1: del, 2: sub, 3: ins, 4: end (to left), 5: end (to right)
			note that for types 3,5 ins is to right of current base */
} Diffdata;

typedef struct diffsegnode {
  Diffdata *diffdata;
  struct diffsegnode *child[2]; /* 0 = left, 1 = right */
  int entry, seg_start, seg_end; /* subject entry, and segment ends */
  int conf_count[2], conf_max_score, conf_max_margin; /* # confirming reads */
} Diffsegnode;

int data_type, data_length, data_qual, data_reverse, data_score, data_margin;
unsigned char *data_seq;

Diffsegnode *head_diffsegnode, data_diffsegnode;

/* Note that upper/lowercase distinction used in last printf --
   may lead to confusion if sequence itself is in lowercase
   -- but currently this applies only to query sequence so should be OK
*/
extern int qseq_trans[256];
extern double *match_score_adjust, *mismatch_score_adjust; /* score adjustments for each quality value, to
								    correct rounding simplifications */

print_diffs(pair, length2)
     Aligned_pair *pair;
     int length2;
{
  unsigned char *diff, *diff1;
  int i, j, length, print_flag, entry2, s2, last_s2, output, qual_threshold;
  unsigned char *seq;
  int site1, site2, type, strand;
  int d_site1, d_site2, large_flag;
  unsigned char *get_seq();
  Query_domain *query_domain;

  seq = get_seq(pair->entry1);
  length = get_seq_length(pair->entry1);
  output = parameters->output_bcdsites;
  qual_threshold = parameters->bcdsites_qual_threshold;
  entry2 = pair->entry2;
  data_reverse = is_reverse(pair) ? 1 : 0;
  site1 = pair->start1 - 1;
  site2 = pair->start2 - 1;
  strand = pair->spl3 & 64;
  if (output) {
    data_diffsegnode.conf_count[data_reverse] = 1;
    data_diffsegnode.conf_count[!data_reverse] = 0;
    data_diffsegnode.conf_max_score = data_score = pair->score;
    query_domain = pair->reversed_pair->query_data->query_domain;
    while (query_domain != query_domain->parent) query_domain = query_domain->parent; /* may be unnec */
    data_diffsegnode.conf_max_margin = data_margin = data_score - 
      (data_score < query_domain->best_score || query_domain->n_best > 1 ?  query_domain->best_score :
      query_domain->next_best ? query_domain->next_best : parameters->min_record_score - 1);
    s2 = data_reverse ? length2 - site2 : site2 + 1;
    if (site1 >= 0 && site2 >= 0) { /* discrepant site!! */
      data_qual = qual[site1];
      printf("\n%s E   %5d  %c(%d)  %5d  ", 
	     parameters->tags ? "DISCREPANCY  " : "",
	     site1 + 1, qseq_trans[seq[site1]], data_qual, s2);
      for (i = site1 > 6 ? site1 - 6 : 0; i < length && i < site1 + 7; i++) 
	printf("%c", i != site1 ? tolower(qseq_trans[seq[i]]) : qseq_trans[seq[i]]);
      if (data_qual >= qual_threshold) {
	data_type = data_reverse? 5 : 4;
	data_length = 1 + (site1 < site2 ? site1 : site2);
	if (data_length > 20) data_length = 20;
	data_seq = seq + site1 - (data_length - 1);
	i = data_reverse ? s2 - 1 : s2 + 1;
	append_diffsegnode(head_diffsegnode, entry2, i, i, 0, &data_diffsegnode);
      }
    }
    last_s2 = s2;
  }
  for (diff = pair->diffs; *(diff + 1); diff++) {
    d_site1 = site1 += diff_gap1(*diff);
    d_site2 = site2 += diff_gap2(*diff);
    if ((type = diff_type(*diff)) == 'M') continue;
    /* note reversal of type : so that refers to change in 1st sequence with respect to the 2d */

    print_flag = 0;
    data_qual = qual[site1];
    diff1 = diff;
    large_flag = 0;
    if (type == 'D') {
      for (diff1 = diff + 1; diff_type(*diff1) == 'D' && !diff_gap2(*diff1); diff1++);
      i = diff1 - diff; 
      if (i > 1) {
	if (diff_type(*diff1) == 'm') {
	  large_flag = 1;
	  for (diff1++; diff_type(*diff1) != 'm'; diff1++) i += diff_gap1(*diff1);
	  for (diff1++; diff_type(*diff1) == 'D' && !diff_gap2(*diff1); diff1++) i++;
	}

	print_flag = 1;
	printf("\n%s %c-%d %5d  ", 
	       parameters->tags ? "DISCREPANCY  " : "", large_flag && strand ? 'i' : 'I',
	       i, site1 + 1);
	data_length = i;
	if (!large_flag) {
	  for (j = data_qual = 0; j < i; j++) {
	    printf("%c", qseq_trans[seq[site1 + j]]);
	    if (qual[site1 + j] > data_qual) data_qual = qual[site1 + j];
	  }
	}
      }
      d_site1 += i - 1;
      diff1--;
    }  
    else if (type == 'I') {
      for (diff1 = diff + 1; diff_type(*diff1) == 'I' && !diff_gap1(*diff1); diff1++);
      i = diff1 - diff; 
      if (i > 1) {
	print_flag = 1;
	if (diff_type(*diff1) == 'm') {
	  large_flag = 1;
	  for (diff1++; diff_type(*diff1) != 'm'; diff1++) i += diff_gap1(*diff1);
	  for (diff1++; diff_type(*diff1) == 'I' && !diff_gap1(*diff1); diff1++) i++;
	}

	printf("\n%s %c-%d %5d  %c", 
	  parameters->tags ? "DISCREPANCY  " : "", large_flag && strand ? 'd' : 'D',
	       i, site1 + 1, qseq_trans[seq[site1]]);
	data_length = i;
      }
      d_site2 += i - 1;
      diff1--;
    }	  
    if (!print_flag) {
      printf("\n%s %c   %5d  %c", 
	  parameters->tags ? "DISCREPANCY  " : "",
	     type == 'D' ? 'I' : (type == 'I' ? 'D' : type),
	     site1 + 1, qseq_trans[seq[site1]]);
      data_length = 1;
    }
    s2 = data_reverse ? length2 - site2 : site2 + 1;
    printf("(%d)  %5d  ", data_qual, s2);
    if (output) {

      if (data_reverse) {
	append_diffsegnode(head_diffsegnode, entry2, type == 'D' ? s2 : s2 + 1, last_s2 - 1, 1, &data_diffsegnode); 
	last_s2 = type == 'I' ? s2 - data_length + 1: s2;
      }
      else {
	append_diffsegnode(head_diffsegnode, entry2, last_s2 + 1, type == 'D' ? s2 : s2 - 1, 1, &data_diffsegnode); 
	last_s2 = type == 'I' ? s2 + data_length - 1: s2;
      }

      if (large_flag || data_qual >= qual_threshold) {
	data_type = type == 'I' ? (large_flag && strand ? 6 : 1) : type == 'D' ? (large_flag && strand ? 7 : 3) : 2;
	data_seq = seq + site1;

	if (data_reverse)
	  i = type == 'I' ? s2 - data_length + 1 : type == 'D' ? s2 - 1 : s2;
	else i = s2;

	append_diffsegnode(head_diffsegnode, entry2, i, i, 0, &data_diffsegnode);
      }
    }
    /*
    if (diff < diff1)
      for (diff++; diff <= diff1; diff++) {
	d_site1 += diff_gap1(*diff);
	d_site2 += diff_gap2(*diff);
      }
    */
    for (i = site1 > 6 ? site1 - 6 : 0; i < length && i < d_site1 + 7; i++) 
      printf("%c", i < site1 || i > d_site1 ? tolower(qseq_trans[seq[i]]) : qseq_trans[seq[i]]);

    diff = diff1;
    site1 = d_site1;
    site2 = d_site2;
  }
  
  if (output) {
    site1 = pair->end1 + 1;
    site2 = pair->end2 + 1;
    s2 = data_reverse ? length2 - site2 : site2 + 1;
    if (data_reverse) {
      append_diffsegnode(head_diffsegnode, entry2, s2 + 1, last_s2 - 1, 1, &data_diffsegnode); 
    }
    else {
      append_diffsegnode(head_diffsegnode, entry2, last_s2 + 1, s2 - 1, 1, &data_diffsegnode); 
    }
    if (site1 < length && site2 < length2) { /* this is a discrepant site!! */
      data_qual = qual[site1];
      printf("\n%s E   %5d  %c(%d)  %5d  ", 
	     parameters->tags ? "DISCREPANCY  " : "",
	     site1 + 1, qseq_trans[seq[site1]], data_qual, s2);
      for (i = site1 > 6 ? site1 - 6 : 0; i < length && i < site1 + 7; i++) 
	printf("%c", i != site1 ? tolower(qseq_trans[seq[i]]) : qseq_trans[seq[i]]);
      if (data_qual >= qual_threshold) { /* extend right */
	data_type = data_reverse? 4 : 5;
	data_length = length - site1 < length2 - site2 ? length - site1 : length2 - site2;
	if (data_length > 20) data_length = 20;
	data_seq = seq + site1;
	i = data_reverse ? s2 + 1 : s2 - 1;
	append_diffsegnode(head_diffsegnode, entry2, i, i, 0, &data_diffsegnode);
      }
    }
  }
}

/* following is slightly wrong -- because omits adjustment for 1 matching base when there is deletion in query */
revise_score(pair)
     Aligned_pair *pair;
{
  double match_tot, mismatch_tot;
  unsigned char *diff;
  int site1, d1, i;
  unsigned char *get_seq(), *seq;

  seq = get_seq(pair->entry1);

  match_tot = mismatch_tot = 0.0;

  site1 = pair->start1 - 1;
  for (diff = pair->diffs; *diff; diff++) {
    d1 = diff_gap1(*diff);
    if (d1 > 1) {
      /* printf("\nd1: %d ", d1 - 1); */
      for (i = site1 + 1; i < site1 + d1; i++) {
	/* printf("%d:%d ", i, qual[i]); */
	match_tot += match_score_adjust[qual[i]];
      }
    }
    site1 += d1;
    if (*(diff + 1) && diff_type(*diff) == 'S' && qseq_trans[seq[site1]] != 'N') {
      /* printf("\nM%d:%d ", site1, qual[site1]); */
      mismatch_tot += mismatch_score_adjust[qual[site1]];
    }
  }
  /* printf("\n%d + match_tot: %.2f + mismatch_tot: %.2f", pair->score, match_tot, mismatch_tot); */
  pair->score += match_tot + mismatch_tot;
  if (pair->score < 0) pair->score = 0;
}

/* N.B. DIFF_HIST SHOULD REFLECT DISCREPANCIES AT BEGINNINGS & ENDS -- BUT DOESN'T!! */

static int *n_align, *n_unalign, *n_p_unalign, *n_subs, *n_dels, *n_ins, *n_Ns, *n_Xs; 

init_diff_hist()
{
  int i;
  char *our_alloc();

  n_Ns = 1 + (int *)our_alloc((MAX_QUAL + 2) * sizeof(int));
  n_Xs = 1 + (int *)our_alloc((MAX_QUAL + 2) * sizeof(int));
  n_subs = 1 + (int *)our_alloc((MAX_QUAL + 2) * sizeof(int));
  n_dels = 1 + (int *)our_alloc((MAX_QUAL + 2) * sizeof(int));
  n_ins = 1 + (int *)our_alloc((MAX_QUAL + 2) * sizeof(int));
  n_align = 1 + (int *)our_alloc((MAX_QUAL + 2) * sizeof(int));
  n_unalign = 1 + (int *)our_alloc((MAX_QUAL + 2) * sizeof(int));
  n_p_unalign = 1 + (int *)our_alloc((MAX_QUAL + 2) * sizeof(int));
  for (i = -1; i <= MAX_QUAL; i++)
    n_subs[i] = n_Ns[i] = n_Xs[i] = n_dels[i] = n_ins[i] 
      = n_align[i] = n_unalign[i] = n_p_unalign[i] = 0;
}
 
incr_diff_hist(seq, diff, length1, length2, start1, start2, end1, end2, reverse, ignore_ends)
     unsigned char *seq;
     unsigned char *diff;
     int length1, length2, start1, start2, end1, end2, reverse, ignore_ends;
{
  int j, q, q1, offset;
  unsigned char *diff1, *orig_diff;
  char c;
  int site1, site2, type, largegap_flag;

  if (reverse) {
    site2 = start1 - 1;
    site1 = start2 - 1;
  }
  else {
    site1 = start1 - 1;
    site2 = start2 - 1;
  }
  largegap_flag = 0;
  orig_diff = diff;
  for (; *(diff + 1); diff++) {
    type = diff_type(*diff);
    if (type == 'm') {
      largegap_flag = largegap_flag ? 0 : diff_type(*(diff - 1));
      continue;
    }
    if (largegap_flag != 'I')
      site1 += diff_gap1(*diff);
    if (largegap_flag != 'D')
      site2 += diff_gap2(*diff);
    if (type == 'M') continue;
    q = qual[reverse ? site2 : site1];
/*
    if (reverse && q >= parameters->qual_show) 
      printf("\n%d: quality %d", site1, q);
*/
    if (type == 'S') {
      c = toupper(seq[reverse ? site2 : site1]);
      if (c == 'N') n_Ns[q] += 1;
      else if (c == 'X') n_Xs[q] += 1;
      else n_subs[q] += 1; 
    }
      
    else if (type == 'D') {
      for (diff1 = diff + 1; diff_type(*diff1) == 'D' && !diff_gap2(*diff1); diff1++)
	if ((q1 = qual[reverse ? site2 : site1 + (diff1 - diff)]) > q) q = q1;
      if ((diff == orig_diff || diff_type(*(diff - 1)) != 'm') && diff_type(*diff1) != 'm')
	n_ins[q] += 1;
      site1 += (diff1 - diff) - 1;
      diff = diff1 - 1;
    } 
    /* note reversal of type : so that refers to change
       if 1st sequence with respect to the 2d */
      
    else if (type == 'I') {
      for (diff1 = diff + 1; diff_type(*diff1) == 'I' && !diff_gap1(*diff1); diff1++)
	if ((q1 = qual[reverse ? site2 + (diff1 - diff) : site1]) > q) q = q1;
      if ((diff == orig_diff || diff_type(*(diff - 1)) != 'm') && diff_type(*diff1) != 'm')
	n_dels[q] += 1;
      site2 += (diff1 - diff) - 1;
      diff = diff1 - 1;
    }
  }

  if (ignore_ends) return;

  offset = start2 - start1;
  for (j = 0; j < start1; j++) {
    if (qual[j] > MAX_QUAL || qual[j] < -1) {
      fprintf(stderr, "\n%d %d %d", j, qual[j], MAX_QUAL);
      fatalError("start1");
    }
    n_unalign[qual[j]] += 1;
    if (j + offset >= 0) n_p_unalign[qual[j]] += 1;
  }
  
  offset = end2 - end1;
  for (j = end1 + 1; j < length1; j++) {
    n_unalign[qual[j]] += 1;
    if (j + offset < length2) n_p_unalign[qual[j]] += 1;
  }
  
  for (j = start1; j <= end1; j++) n_align[qual[j]] += 1;
}

print_diff_hist(label)
     char *label;
{
  int i, t_errs, cum_align, rcum_align, rcum_errs, cum_errs, t_n_align, t_n_errs;

  printf("\n\n%s\nQual algn  cum    rcum    (%%)    unalgn X    N  sub del ins  total (%%)   cum  rcum (%%)", label);
  for (t_n_errs = t_n_align = i = 0; i <= MAX_QUAL; i++) {
    t_n_align += n_align[i];
    t_n_errs += n_subs[i] + n_dels[i] + n_ins[i] + n_Ns[i];
  }

  cum_align = cum_errs = 0;
  for (i = MAX_QUAL; i >= -1; i--) {
    if (!n_align[i] && !n_unalign[i]) continue;
    rcum_align = t_n_align - cum_align;
    cum_align += n_align[i];
    rcum_errs = t_n_errs - cum_errs;
    t_errs = n_subs[i] + n_dels[i] + n_ins[i] + n_Ns[i];
    cum_errs += t_errs;
    printf("\n%2d %6d %6d %6d (%6.2f)  %4d %2d   %2d %3d %3d %3d   %3d (%.2f)  %3d  %3d (%.2f)", 
	   i, /* n_align[i] + n_unalign[i], */
/*	   length ? (100.0 * (n_align[i] + n_unalign[i])) / length : 0.0, */
	   n_align[i], cum_align, 
/*	   t_n_align && i >= 0 ? 100.0 * cum_align / t_n_align : 0.0, */
	   rcum_align,
	   t_n_align && i >= 0 ? 100.0 * rcum_align / t_n_align : 0.0,
	   n_p_unalign[i], 
	   n_Xs[i], n_Ns[i], n_subs[i], n_dels[i], n_ins[i],
	   t_errs, 
	   n_align[i] ? (100.0 * t_errs) / n_align[i] : 0.0, 
	   cum_errs, 
/*	   cum_align && i >= 0 ? (100.0 * cum_errs) / cum_align : 0.0, */
	   rcum_errs, 
	   rcum_align && i >= 0 ? (100.0 * rcum_errs) / rcum_align : 0.0);
  }
  printf("\n");
}

free_diff_hist()
{
  our_free(n_Ns - 1);
  our_free(n_Xs - 1);
  our_free(n_subs - 1);
  our_free(n_dels - 1);
  our_free(n_ins - 1);
  our_free(n_align - 1);
  our_free(n_unalign - 1);
  our_free(n_p_unalign - 1);
}

/* routines for recording difference/confirmation info for subject sequences, 
   in resequencing project
*/

Diffsegnode *node_block;
int node_index;
extern unsigned char area_comp_mat[];

append_diffdata(diffsegnode)
     Diffsegnode *diffsegnode;
{
  char *our_alloc();
  static Diffdata *head_diffdata;
  Diffdata *diffdata;
  static int diffdata_index;
  static unsigned char *head_chardata;
  static int chardata_index;
  int i;

  for (diffdata = diffsegnode->diffdata; diffdata; diffdata = diffdata->next) {
    if (diffdata->type != data_type || diffdata->length != data_length) continue; 
    if (data_type != 1 && data_type != 6) {
      if (data_reverse) 
	for (i = 0; i < data_length && diffdata->seq[i] == area_comp_mat[data_seq[data_length - 1 - i]]; i++);
      else
	for (i = 0; i < data_length && diffdata->seq[i] == data_seq[i]; i++);
    }
    if (data_type == 1 || data_type == 6 || i == data_length) {
      if (diffdata->max_qual < data_qual) diffdata->max_qual = data_qual;
      if (diffdata->max_score < data_score) diffdata->max_score = data_score;
      if (diffdata->max_margin < data_margin) diffdata->max_margin = data_margin;
      diffdata->count[data_reverse] += 1;
      return;
    }
  }

  if (!diffdata_index) {
    head_diffdata = (Diffdata *)our_alloc(1000 * sizeof(Diffdata));
  }
  diffdata = head_diffdata + diffdata_index;
  diffdata_index = (diffdata_index + 1) % 1000;
  diffdata->next = diffsegnode->diffdata;
  diffsegnode->diffdata = diffdata;

  /* allocate sequence & copy */
  
  if (data_type != 1 && data_type != 6) {
    if (!chardata_index || chardata_index + data_length > 998) {
      head_chardata = (unsigned char *)our_alloc(1000 * sizeof(unsigned char));
      head_chardata[0] = 0;
      chardata_index = 1;
    }
    diffdata->seq = head_chardata + chardata_index;
    if (data_reverse) {
      for (i = 0; i < data_length; i++) diffdata->seq[i] = area_comp_mat[data_seq[data_length - 1 - i]];
    }
    else
      for (i = 0; i < data_length; i++) diffdata->seq[i] = data_seq[i];
    diffdata->seq[i] = 0;
    chardata_index += data_length + 1;
  }
  else diffdata->seq = 0;

  diffdata->type = data_type;
  diffdata->length = data_length;
  diffdata->max_qual = data_qual;
  diffdata->max_score = data_score;
  diffdata->max_margin = data_margin;
  diffdata->count[data_reverse] = 1;
  diffdata->count[!data_reverse] = 0;
}

append_diffsegnode(node, entry, seg_start, seg_end, conf_flag, data_node)
     Diffsegnode *node, *data_node; /* latter is 0 if reading from global vars */
     int entry, seg_start, seg_end, conf_flag;
{
  char *our_alloc();
  int i, s, direc, min_start, min_end, max_start, max_end;
  Diffsegnode *new_node;

  if (seg_start > seg_end) return;

  /* printf("\n%d-%d(%d)", seg_start, seg_end, conf_flag); */
  
  if (!head_diffsegnode || entry < node->entry || entry == node->entry && seg_end < node->seg_start) 
    direc = 0;
  else if (entry > node->entry || entry == node->entry && seg_start > node->seg_end)
    direc = 1;
  else { /* new segment overlaps current node; so take intersection */
    max_start = seg_start > node->seg_start ? seg_start : node->seg_start;
    min_end = seg_end < node->seg_end ? seg_end : node->seg_end;
    min_start = seg_start < node->seg_start ? seg_start : node->seg_start;
    max_end = seg_end > node->seg_end ? seg_end : node->seg_end;
    /* 1st segment: from min_start to max_start - 1; 2d from max_start to min_end; last: from min_end + 1 to max_end */

/* cases: 1) adding diff site (always size 1)
               i) coincides with existing site
               ii) new site, not contained
               iii) lies within confirmed segment
          2) adding confirmed segment (more than one below can apply
                i) overlaps confirmed segment
                ii) overlaps site

*/ 
    
    node->seg_start = max_start;
    node->seg_end = min_end;

    /* intersection: must be single site, if either new one is, or old one was; this
       is only node to receive diffdata
    */

    if (min_start < max_start) { 
      append_diffsegnode(node, entry, min_start, max_start - 1, 1,
			 min_start == seg_start ? data_node : node);
    }
    if (min_end < max_end) {
      append_diffsegnode(node, entry, min_end + 1, max_end, 1,
			 max_end == seg_end ? data_node : node);
    }

    if (conf_flag) {
      for (i = 0; i < 2; i++)
	node->conf_count[i] += data_node->conf_count[i];
      if (node->conf_max_score < data_node->conf_max_score) 
	node->conf_max_score = data_node->conf_max_score;
      if (node->conf_max_margin < data_node->conf_max_margin
	  || !(node->conf_count[0] + node->conf_count[1])) 
	node->conf_max_margin = data_node->conf_max_margin;
    }
    else 
      append_diffdata(node);
    return;
  }
  if (!head_diffsegnode || !node->child[direc]) {
    if (!(node_index % 100)) {
      node_block = (Diffsegnode *)our_alloc(100 * sizeof(Diffsegnode));
      node_index = 0;
    }
    new_node = node_block + node_index++;
    if (!head_diffsegnode) {
      head_diffsegnode = new_node;
    }
    else node->child[direc] = new_node;
    new_node->entry = entry;
    new_node->seg_start = seg_start;
    new_node->seg_end = seg_end;
    new_node->conf_count[0] = new_node->conf_count[1] = 0;
    new_node->conf_max_score = new_node->conf_max_margin = 0;
    new_node->diffdata = 0;
    if (conf_flag) {
      for (i = 0; i < 2; i++)
	new_node->conf_count[i] = data_node->conf_count[i];
      new_node->conf_max_score = data_node->conf_max_score;
      new_node->conf_max_margin = data_node->conf_max_margin;
    }
    else
      append_diffdata(new_node);
    new_node->child[0] = new_node->child[1] = 0;
  }
  else {
    append_diffsegnode(node->child[direc], entry, seg_start, seg_end, conf_flag, data_node);
  }
  return;
}

int last_entry, last_seg_end;
FILE *fp_diffseg;
char cdb[10];
int diff_seg_hist[101];

open_diffseg_file()
{
  FILE *fopenWrite();
  char filename[200];
  int i;

  sprintf(filename, "%s.bcdsites", parameters->query_files->name);
  fp_diffseg = fopenWrite(filename);

  fprintf(fp_diffseg, "#");
  for (i = 0; i < parameters->argc; i++) 
    fprintf(fp_diffseg, " %s", parameters->argv[i]);
  fprintf(fp_diffseg, "\n# %s version %s", parameters->calling_program, parameters->version);
  fprintf(fp_diffseg, "\n# Run date:time  %s\n", parameters->date);

  strcpy(cdb, "CDSILRdi");
}

write_diffsegnode(node)
     Diffsegnode *node;
{
  char *get_id();
  int i, j, n, merge_flag;
  Diffdata *diffdata, *diffdata2;

  if (!node) return;
  write_diffsegnode(node->child[0]);
  /*
  if (last_entry != node->entry) {
    last_entry = node->entry;
    fprintf(fp_diffseg, "\n\n%s:", get_id(last_entry));
    last_seg_end = -1;
  }
  */
  /*
  if (node->seg_start > last_seg_end + 10)
    fprintf(fp_diffseg, " (gap %d)\n", node->seg_start - (last_seg_end + 1));
  */
  fprintf(fp_diffseg, "%s %d", get_id(node->entry), node->seg_start);
  if (node->conf_count[0] + node->conf_count[1]) {
    fprintf(fp_diffseg, " C");
    if (node->seg_start != node->seg_end)
      fprintf(fp_diffseg, "-%d", node->seg_end - node->seg_start + 1);
    fprintf(fp_diffseg, ",%d,%d,%d,%d,0", 
	    node->conf_count[0], node->conf_count[1],
	    node->conf_max_score, node->conf_max_margin);
  }
  for (diffdata = node->diffdata; diffdata; diffdata = diffdata->next) {
    if (diffdata->type > 3 && diffdata->type < 6) { /* must consider potential merges */
      merge_flag = 0;
      for (diffdata2 = diffdata->next; diffdata2; diffdata2 = diffdata2->next) {
	if (diffdata2->type != diffdata->type || diffdata2->length == diffdata->length) continue; /* cannot be perfectly overlapping */
	if (diffdata->type == 5)
	  for (i = j = 0; diffdata->seq[i] == diffdata2->seq[j]; i++, j++);
	else 
	  for (i = diffdata->length - 1, j = diffdata2->length - 1; diffdata->seq[i] == diffdata2->seq[j]; i--, j--);
	if (!diffdata->seq[i] || !diffdata2->seq[j]) { /* one is contained in other -- so merge into second one, & skip current one */
	  merge_flag = 1;
	  if (diffdata2->max_qual < diffdata->max_qual) diffdata2->max_qual = diffdata->max_qual;
	  if (diffdata2->max_score < diffdata->max_score) diffdata2->max_score = diffdata->max_score;
	  if (diffdata2->max_margin < diffdata->max_margin) diffdata2->max_margin = diffdata->max_margin;
	  for (i = 0; i < 2; i++)
	    diffdata2->count[i] += diffdata->count[i];
	  if (diffdata2->length < diffdata->length) {
	    diffdata2->length = diffdata->length;
	    diffdata2->seq = diffdata->seq;
	  }
	  break;
	}
      }
      if (merge_flag) continue;
    }
    fprintf(fp_diffseg, " %c", cdb[diffdata->type]);
    if (diffdata->type == 1 || diffdata->type == 6) {
      if (diffdata->length != 1) {
	fprintf(fp_diffseg, "-%d", diffdata->length);
	if (diffdata->length >= parameters->min_intron_length) {
	  n = diffdata->length / 100;
	  diff_seg_hist[n > 100 ? 100 : n] += 1;
	}
      }
    }
    else
      fprintf(fp_diffseg, "-%s", diffdata->seq); 
    fprintf(fp_diffseg, ",%d,%d,%d,%d,%d", 
	    diffdata->count[0], diffdata->count[1],
	    diffdata->max_score, diffdata->max_margin, diffdata->max_qual);
  }
  fprintf(fp_diffseg, "\n");

  last_seg_end = node->seg_end;
  write_diffsegnode(node->child[1]);
}

write_diffsegnodes()
{
  int i;
  double cum, tot;

  for (i = 0; i < 101; i++) diff_seg_hist[i] = 0;
  write_diffsegnode(head_diffsegnode);

  printf("\n\nIntron size histogram");
  for (i = tot = 0; i < 101; i++) tot += diff_seg_hist[i];
  for (i = 100, cum = 0; i >= 0; i--)
    if (diff_seg_hist[i]) {
      cum += diff_seg_hist[i];
      printf("\n%5d %5d %5.0f %.3f", i * 100, diff_seg_hist[i], cum, cum / tot);
    }
}

