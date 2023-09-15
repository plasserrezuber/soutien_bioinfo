/*****************************************************************************
#   Copyright (C) 1994-1999, 2006-2009 by Phil Green.                          
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

/* diff representation of alignments: single byte with high order 2 bits indicating
   type ('S', 'I', 'D', or 'M'), low order bits indicating spacing since previous diff.
   Special case: diff == 63 (assigned type 'm') (which does not occur with actual 'M's) indicates
   boundaries of 'large' indel (e.g. intron). This diff must occur after 1 or more diffs
   with type 'I' or 'D' at beginning at large gap and again preceding diffs of same
   type ('I' or 'D') at end of indel; between the diffs of type 'm' there occur diffs
   of type 'M' whose collective low order bits indicate the (remaining) size of the large gap.
*/

extern Parameters *parameters;
extern int n_c_rows, row_convert[];
extern double lambda, *log_freqs[2];

/* the following retrieve additional information -- after the initial call */
/*BUG 0 SCORE ALIGNMENTS DO NOT HAVE INFORMATION APPROPRIATELY SET */
static unsigned char *string1, *string2, *string3;
static int string_end;
static int query_start, query_end, subject_start, subject_end;
static int mismatches, insertions, deletions;
static char warnings, complexity_adjust;
static double complexity_factor,  complexity_lambda, *complexity_log_freqs;
static int complexity_n_freq_indices, complexity_counts[256], complexity_freq_type, *complexity_convert;
static int score;
static Profile *q_profile;
static unsigned char *db_entry_seq;
static int db_entry_length;
static int q_left, s_left; /* boundaries of the query and 
			      subject sequences to be used -- origin 0 */
static int band_width, r_edge_orig, l_edge_orig;
static int gap_init, ins_gap_ext, del_gap_ext;
static int string_len2;
static double nlogn[101], nlogf_n[101], log_factor;

/* should recheck complexity adjustment procedure: it is probably not optimal! */
/* N.B. must be called AFTER set_score_mat, not before!! */ 
set_complexity_adjust(factor, type) /* type = 0 for row freqs, 1 for col freqs */
     double factor;
     int type;
{
  int n;

  complexity_adjust = 1;
  complexity_factor = factor;
  complexity_freq_type = type;
  nlogn[0] = nlogf_n[0] = 0;
  for (n = 1; n < 101; n++) {
    nlogn[n] = n * log((double)n);
    nlogf_n[n] = n * log(factor / n);
  }
  log_factor = factor;

  /* defaults -- if type == 1,  reset with q_profile specific data */
  complexity_n_freq_indices = n_c_rows;
  complexity_lambda = lambda;
  complexity_log_freqs = log_freqs[0];
  complexity_convert = row_convert;
}

unset_complexity_adjust()
{
  complexity_adjust = 0;
}

warnings_on()
{
  warnings = 1;
}

warnings_off()
{
  warnings = 0;
}

init_alignment_globals(q_p, db_e_seq, db_e_len, b_w, r_e, l_e, q_l, s_l)
     Profile *q_p;
     unsigned char *db_e_seq;
     int db_e_len;
     int q_l, s_l; 
     int b_w, r_e, l_e;
{
  subject_end = query_end = subject_start = query_start = 0;
  mismatches = insertions = deletions = 0;
  
  q_profile = q_p;
  db_entry_seq = db_e_seq;
  db_entry_length = db_e_len;
  band_width = b_w;
  r_edge_orig = r_e;
  l_edge_orig = l_e;
  q_left = q_l;
  s_left = s_l;
  gap_init = q_profile->gap_init;
  ins_gap_ext = q_profile->ins_gap_ext;
  del_gap_ext = q_profile->del_gap_ext;

  if (complexity_freq_type) {
    complexity_n_freq_indices = q_profile->n_c_cols;
    complexity_lambda = q_profile->lambda;
    complexity_log_freqs = q_profile->log_col_freqs;
    complexity_convert = q_profile->convert;
  }
}

alloc_strings(restart, t_len)
     int restart, t_len;
{
  char *our_alloc();
  int i, k, i1, new_i1;
  unsigned char *new_string1;
  static int post_first, string_len, align_extend, align_extend2;

  if (post_first) {
    align_extend = parameters->splice_edge_length;
    align_extend2 = 2 * align_extend;
    post_first = 1;
  }

  t_len += align_extend2;
  if (t_len >= string_len) {
    fprintf(stderr, " Allocating alignment string: %d %d ", t_len, string_len);
    t_len = 2 * t_len + 1000;
    new_string1 = (unsigned char *)our_alloc(3 * t_len * sizeof(unsigned char));
    if (string_len) {
     /* only free output strings when need a new one -- so
	 are available for use between initial comparison and
	 subsequent one */

      string1 -= align_extend;
      if (!restart) { /* copy */
	for (k = 0; k < 3; k++) {
	  /* for (i = 0; i < 3 * string_len; i++) */
	  new_i1 = k * t_len;
	  i1 = k * string_len;
	  for (i = 0; i < string_len; i++)
	    new_string1[i + new_i1] = string1[i + i1];
	}
      }
      our_free(string1);
    }
    string1 = new_string1 + align_extend;
    string2 = string1 + t_len;
    string3 = string2 + t_len;
    string_len = t_len;
    string_len2 = string_len - align_extend2 - 1;
  }

  if (restart) {
    string_end = -1;

    if (complexity_adjust) {
      for (i = 0; i < complexity_n_freq_indices; i++) complexity_counts[i] = 0;
    }
  }
}

alignment_from_direc(mdir, mdiri, direc, add_max_score)
     char **direc, **mdir;
     char *mdiri;
     int *add_max_score;
{
  char edge_flag, d;
  int i, j, k;
  int s, score_check, cum_neg, min_cum_neg;
  char u_ext, t_ext; /* flags for whether current u and t are
			gap extending */
  int pos_score;
  
  edge_flag = 0;
  
  alloc_strings(1, 10); /* db_entry_length + q_profile->length); /* 10);*/
  
  u_ext = t_ext = 0;
  score_check = pos_score = cum_neg = min_cum_neg = 0;
  subject_end = (mdir - direc) + 1;
  query_end = (mdiri - *mdir) + 1;
  
  for (i = mdir - direc, j =  mdiri - *mdir, k = 0;
       i >= s_left && j >= q_left && (d = direc[i][j]); k++) {
    
    if (band_width && (j - i < l_edge_orig + 2 || j - i > r_edge_orig - 2))
      edge_flag = 1;
    
    string2[k] = string1[k] = '-';
    
    if (t_ext) {
      d = 1 + (d & 12);
    }
    if (u_ext) {
      d = 2 + (d & 12);
    }
    if (!(d & 2)) t_ext = d & 4;
    if (!(d & 1)) u_ext = d & 8;
    
    if (d & 1) { 
      string2[k] = db_entry_seq[i];
      i--;
    }
    if (d & 2) {
      string1[k] = q_profile->seq[j];
      j--;
    }
    string3[k] = ' ';
    if (string1[k] == '-') {
      s = string1[k - 1] == '-' ? ins_gap_ext : gap_init;
      score_check += s;
      string3[k] = '-';
      insertions++;
    }
    else if (string2[k] == '-') {
      s = string2[k - 1] == '-' ? del_gap_ext : gap_init;
      score_check += s;
      string3[k] = '-';
      deletions++;
    }
    else /* if (string2[k] != '-' && string1[k] != '-') */ {
      s = q_profile->scores[q_profile->convert[string2[k]]][j + 1];
      score_check += s;
      if (complexity_adjust) {
	pos_score += s;
/*
	complexity_counts[q_profile->convert[string1[k]]] += 1;
*/
	complexity_counts[complexity_convert[complexity_freq_type ? string2[k] : string1[k]]] += 1;
      }
      if (s > 0) {
	string3[k] = string2[k] == string1[k] ? string1[k] : '+';
/*
	if (complexity_adjust) {
	  pos_score += s;
	  complexity_counts[complexity_convert[string2[k]]] += 1;
	}
*/
      }
      else {
	mismatches++; /* previously only was set if s < 0 */
	/*	fprintf(stderr, " 1:%c;%d;%d ", string2[k], s, mismatches); */
	/*
	if (string3[k] != ' ') fprintf(stderr, "%d-%c ", string3[k], string3[k]);
	fprintf(stderr, "(%d)", mismatches);
	*/
	/* fprintf(stderr, " %d:%d ", mismatches, k); */
      }
    }
    cum_neg += s;
    if (cum_neg > 0) cum_neg = 0;
    else if (min_cum_neg > cum_neg) min_cum_neg = cum_neg;
    /* notify("1"); /* */
    if (k >= string_len2)
      alloc_strings(0, k + 1); 
  }
  if (score_check != *add_max_score) 
    fprintf(stderr, "\nSCORE DISCREPANCY: initial %d, from alignment %d  %d %d",
	    *add_max_score, score_check, q_left, s_left);
  
  i += 2;
  j += 2;
  subject_start = i;
  query_start = j;
  string_end = k - 1;
  
  if (complexity_adjust) {
    *add_max_score = complexity_correct(*add_max_score, complexity_counts);
/*
    fprintf(stderr, "\nmin_cum_neg = %d", min_cum_neg);
*/
  }
  if ((band_width || q_left) && edge_flag && warnings)
    fprintf(stderr, "\nWARNING: Possible band edge overflow: score %d",
	    *add_max_score);
  score = *add_max_score;
  /*
  for (i = j = 0; i <= string_end; i++) 
    if (string3[i] == ' ') fprintf(stderr, " %d:%d ", ++j, i);
  */
}

int alignment_from_max_list(n_maxes, max_list, add_score, l_edge_orig, r_edge_orig)
     int n_maxes;
     Max_list *max_list;
     int *add_score; 
     int l_edge_orig, r_edge_orig;
{
  int i, j, k, t_score, target_score, s, pos_score;
  int s1, s2, s3, j1, k1;
  char indel_flag, active_flag1, active_flag2, active_flag3, block_succeed;

  if (!n_maxes) return 0;
  i = n_maxes - 1;
  j = max_list[i].q_pos; 
  k = max_list[i].s_pos;
  query_end = j + 1;
  subject_end = k + 1;
  block_succeed = indel_flag = 0;

  pos_score = 0;
  target_score = max_list[i].score;
  score = 0; 
  alloc_strings(1, 10); /* db_entry_length + q_profile->length); /* 10); */

  /*
  if (k == 50) {
    for (k1 = 0; k1 < 30; k1++) {
      printf("\nDEBUG: %d %d %d %c (%d) %c %d %d %d %d %d %d %d %d", 
	     q_profile->scores[q_profile->convert[db_entry_seq[k - k1]]][j - k1], j + 1 - k1, k + 1- k1,
	     db_entry_seq[k - k1], q_profile->convert[db_entry_seq[k - k1]], q_profile->seq[j - k1],
	     q_profile->scores[q_profile->convert['A']][j - k1],
	     q_profile->scores[q_profile->convert['C']][j - k1],
	     q_profile->scores[q_profile->convert['G']][j - k1],
	     q_profile->scores[q_profile->convert['T']][j - k1],
	     q_profile->scores[q_profile->convert['k']][j - k1],
	     q_profile->scores[q_profile->convert['m']][j - k1],
	     q_profile->scores[q_profile->convert['r']][j - k1],
	     q_profile->scores[q_profile->convert['s']][j - k1]
	     );
    }
  }
  */

  for (; ;) {
    s = q_profile->scores[q_profile->convert[db_entry_seq[k]]][j];
    if (s <= 0 && !block_succeed && j && k 
	&& (q_profile->scores[q_profile->convert[db_entry_seq[k - 1]]][j] > s
	    || q_profile->scores[q_profile->convert[db_entry_seq[k]]][j - 1] > s) 
	&& score + gap_init > 0) {
      /* try looking at entire block first */
      if (i) {
	for (j1 = j - 1, k1 = k - 1, t_score = score + s; t_score > 0 &&
	     j1 > max_list[i - 1].q_pos && k1 > max_list[i - 1].s_pos; j1--, k1--) 
	  t_score += q_profile->scores[q_profile->convert[db_entry_seq[k1]]][j1];
	if (t_score > 0) {
	  if (j1 > max_list[i - 1].q_pos) {
	    t_score += gap_init + del_gap_ext * (j1 - max_list[i - 1].q_pos - 1);
	  }
	  else if (k1 > max_list[i - 1].s_pos) {
	    t_score += gap_init + ins_gap_ext * (k1 - max_list[i - 1].s_pos - 1);
	  }
	  if (t_score + max_list[i - 1].score == target_score) {
	    block_succeed = 1;
	  }
	}
	if (!block_succeed && q_profile->scores[q_profile->convert[db_entry_seq[k]]][j - 1] > s) {
	  for (j1 = j - 1, k1 = k, t_score = score + gap_init; t_score > 0 &&
	       j1 > max_list[i - 1].q_pos && k1 > max_list[i - 1].s_pos; j1--, k1--) 
	    t_score += q_profile->scores[q_profile->convert[db_entry_seq[k1]]][j1];
	  if (t_score > 0) {
	    if (j1 > max_list[i - 1].q_pos) {
	      t_score += gap_init + del_gap_ext * (j1 - max_list[i - 1].q_pos - 1);
	    }
	    else if (k1 > max_list[i - 1].s_pos) {
	      t_score += gap_init + ins_gap_ext * (k1 - max_list[i - 1].s_pos - 1);
	    }
	    if (t_score + max_list[i - 1].score == target_score) {
	      block_succeed = 1;
	      k++;
	      indel_flag = 3;
	    }
	  }
	}
	if (!block_succeed && q_profile->scores[q_profile->convert[db_entry_seq[k - 1]]][j] > s) {
	  for (j1 = j, k1 = k - 1, t_score = score + gap_init; t_score > 0 &&
	       j1 > max_list[i - 1].q_pos && k1 > max_list[i - 1].s_pos; j1--, k1--) 
	    t_score += q_profile->scores[q_profile->convert[db_entry_seq[k1]]][j1];
	  if (t_score > 0) {
	    if (j1 > max_list[i - 1].q_pos) {
	      t_score += gap_init + del_gap_ext * (j1 - max_list[i - 1].q_pos - 1);
	    }
	    else if (k1 > max_list[i - 1].s_pos) {
	      t_score += gap_init + ins_gap_ext * (k1 - max_list[i - 1].s_pos - 1);
	    }
	    if (t_score + max_list[i - 1].score == target_score) {
	      block_succeed = 1;
	      j++;
	      indel_flag = 2;
	    }
	  }
	}

/*	  notify("success"); */
      }

      if (!block_succeed) {
	active_flag1 = active_flag2 = active_flag3 = 1;
	
	for (j1 = j - 1, k1 = k - 1, s1 = s, s2 = s3 = gap_init; 
	     (active_flag2 || active_flag3) && j1 >= 0 && k1 >= 0; j1--, k1--) {
	  if (active_flag1) s1 += q_profile->scores[q_profile->convert[db_entry_seq[k1]]][j1];
	  if (active_flag2) {
	    s2 += q_profile->scores[q_profile->convert[db_entry_seq[k1 + 1]]][j1];
	  }
	  if (active_flag3) {
	    s3 += q_profile->scores[q_profile->convert[db_entry_seq[k1]]][j1 + 1];
	  }
	  if (active_flag2 && 
	      (active_flag1 && s1 + gap_init > s2 || active_flag3 && s3 + gap_init > s2)) 
	    active_flag2 = 0;
	  if (active_flag3 && 
	      (active_flag1 && s1 + gap_init > s3 || active_flag2 && s2 + gap_init > s3)) 
	    active_flag3 = 0;
	  if (active_flag1 && 
	      (active_flag2 && s2 + gap_init > s1 || active_flag3 && s3 + gap_init > s1)) 
	    active_flag1 = 0;
	  
	  if (active_flag2 && !active_flag1 && s2 > s3) {
	    k++;
	    indel_flag = 3;
	    break;
	  }
	  if (active_flag3 && !active_flag1 && s3 > s2) {
	    j++;
	    indel_flag = 2;
	    break;
	  }
	}
      }
    }
    if (indel_flag) {
      string_extend((int)indel_flag, q_profile->seq[j], db_entry_seq[k]);
      score += gap_init;
      indel_flag = 0;
    }
    else {
      score += s;
      string_extend(s <= 0, q_profile->seq[j], db_entry_seq[k]);
      if (complexity_adjust) {
	pos_score += s; 
/*
	complexity_counts[complexity_convert[q_profile->seq[j]]] += 1;
*/
	complexity_counts[complexity_convert[complexity_freq_type ? db_entry_seq[k] : q_profile->seq[j]]] += 1;

      }
/*
      if (s > 0 && complexity_adjust) {
	pos_score += s;
	complexity_counts[complexity_convert[q_profile->seq[j]]] += 1;
      }
*/
    }
    if (score >= target_score) { /* success */
      if (score > target_score) notify("\nScore exceeded");
      subject_start = k + 1;
      query_start = j + 1;
      score = *add_score = max_list[n_maxes - 1].score;
      if (complexity_adjust) {
	*add_score = complexity_correct(score, complexity_counts);
      }
      return 1;
    }
    else if (score <= 0 || j <= 0 || k <= 0 || j - k < l_edge_orig || j - k > r_edge_orig) { /* failure */
/*
      subject_start = s_left + 1;
      query_start = q_left + 1;
*/
      subject_start = k + 1;
      query_start = j + 1;

      score = max_list[n_maxes - 1].score;
      return 0;
    }
    j--;
    k--;
    if (i && (j == max_list[i - 1].q_pos || k == max_list[i - 1].s_pos)) {
      if (j > max_list[i - 1].q_pos) {
	s = gap_init + del_gap_ext * (j - max_list[i - 1].q_pos - 1);
      }
      else if (k > max_list[i - 1].s_pos) {
	s = gap_init + ins_gap_ext * (k - max_list[i - 1].s_pos - 1);
      }
      else s = 0;
      if (score + s + max_list[i - 1].score == target_score) {
	for (; j > max_list[i - 1].q_pos; j--) string_extend(3, q_profile->seq[j], '-');
	for (; k > max_list[i - 1].s_pos; k--) string_extend(2, '-', db_entry_seq[k]);
	target_score = max_list[i - 1].score;
	block_succeed = 0;
	score = 0;
      }
      i--;
    }
  }
}

string_extend(type, q_res, s_res)
     int type; /* 0 -- matching pair of residues; 1 -- mismatch;
		  2 -- insertion; 3 -- deletion */
     char q_res, s_res;
{
  string_end++;
  /*  notify("2"); /* */
  if (string_end >= string_len2)
    alloc_strings(0, string_end + 1);
  if (type < 2) {
    string2[string_end] = s_res;
    string1[string_end] = q_res;
    if (type) {
      mismatches++;
      /*      fprintf(stderr, " 2:%d ", mismatches); */
      string3[string_end] = ' ';
      /* fprintf(stderr, "(%d)* ", mismatches); */
    }
    else {
      string3[string_end] = (s_res == q_res ? s_res : '+');
    }
  }
  else if (type == 2) {
    string3[string_end] = string1[string_end] = '-';
    string2[string_end] = s_res;
    insertions++;
  }
  else {
    string3[string_end] = string2[string_end] = '-';
    string1[string_end] = q_res;
    deletions++;
  }
}

/* must set complexity_lambda, *complexity_log_freqs, complexity_n_freq_indices before calling */

int complexity_correct(orig_score, counts) /* previously passed pos_score also */
     int orig_score, *counts;
{
  int i, adj_score, old_adj_score, t_counts, n_letters;
  double t_factor, t_sum;
  
  for (i = t_counts = t_factor = t_sum = /* n_letters = */ 0; i < complexity_n_freq_indices; i++) {

    if (counts[i]) {
      if (complexity_log_freqs[i]) {
	t_factor += counts[i] < 101 ? nlogn[counts[i]] : counts[i] * log((double)counts[i]);
	t_sum += counts[i] * complexity_log_freqs[i];
	t_counts += counts[i];
	/*	n_letters++; */
      }
    }
  }
 
  if (t_counts) 
    t_factor -= t_counts < 101 ? nlogn[t_counts] : t_counts * log((double)t_counts);
  t_sum -= t_factor; /* could instead compute using nlogf_n -- but need to verify it is correct */
  /*  printf(" %.2f %.3f\n", t_sum / c_lambda, c_lambda); */

  /* out 4/25/08
  t_factor /=  t_counts * log(n_letters > 1 / complexity_factor ? 1. / n_letters : complexity_factor);
  old_adj_score = orig_score + t_factor * pos_score - pos_score + .5;
  */
  adj_score = orig_score + t_sum / complexity_lambda + .999;

  if (adj_score < 0) adj_score = 0;
/*
  fprintf(stderr, "\nORIG, old/new ADJUSTED SCORES: %d, %d / %d  t_factor %f  pos_score %d  n_letters %d", 
	  orig_score, old_adj_score, adj_score, t_factor, pos_score, n_letters);  
  fprintf(stderr, "\nc_lambda %f, t_sum %f", c_lambda, t_sum);
*/
/*
    fprintf(stderr, "\nScore increase: %d %d", orig_score, adj_score);
    for (i = 0; i < 256; i++) 
      if (counts[i])
	fprintf(stderr, "\n %d %d", i, counts[i]);
*/
  if (adj_score > orig_score) {
    fprintf(stderr, "\nScore increase: %d %d", orig_score, adj_score);
    for (i = 0; i < complexity_n_freq_indices; i++) 
      if (counts[i])
	fprintf(stderr, "\n %c %d", (char)i, counts[i]);
  }
  return adj_score;
}

print_alignment(profile)
     Profile *profile;
{
  int i, j, k, k0, alen, hist;

  if (!score) return; /* variable values may be undefined in this case */
  i = subject_start;
  j = query_start;
  hist = parameters->query_histograms;
  for (alen = string_end; alen >= 0; alen -= 50) {
    printf("\n\nSubject %5d ",i);
    for (k = alen; k > alen - 50 && k >= 0; k--) {
      printf("%c",string2[k]);
      if (string2[k] != '-') i++;
    }
    printf(" %d\n              ", i - 1);
    for (k = alen; k > alen - 50 && k >= 0; k--) printf("%c",string3[k]);
    printf("\nQuery   %5d ", j);
    for (k = alen; k > alen - 50 && k >= 0; k--) {
      printf("%c",string1[k]);
      if (string1[k] != '-') {
	if (hist) {
	  profile->hist[j - 1][(string2[k] == '-') + profile->convert[string2[k]]] += 1;
	}
	j++;
      }
      else {
	if (hist && string1[k + 1] != '-') {
	  for (k0 = k - 1; string1[k0] == '-'; k0--);
	  k0 = k - k0;
	  if (k0 > 10) k0 = 10;
	  profile->hist[j - 1][1 + profile->convert['-'] + k0] += 1;
	}
      }
    }
    printf(" %d", j - 1);
  }
}

get_stats(add_query_start, add_query_end, 
	  add_subject_start, add_subject_end, 
	  add_mismatches, add_insertions, add_deletions)
  int *add_query_start, *add_query_end, *add_subject_start, *add_subject_end,
      *add_mismatches, *add_insertions, *add_deletions;
{
  *add_query_start = query_start;
  *add_query_end = query_end;
  *add_subject_start = subject_start;
  *add_subject_end = subject_end;
  *add_mismatches = mismatches;
  *add_insertions = insertions;
  *add_deletions = deletions;
}

/*
Diff *make_diffs()
{
  char *our_alloc();
  int k, n_diffs;
  int n_i, n_d, n_s;

  n_diffs = mismatches + insertions + deletions;
  diffs = (Diff *)our_alloc((n_diffs + 2) * sizeof(Diff));

  i_diffs = n_diffs + 1;
  set_diff(query_end + 1, subject_end + 1, 'E');

  i = subject_end;
  j = query_end;
  n_i = n_d = n_s = 0;

  for (k = 0; i_diffs > 0; k++) {
    if (string1[k] == '-') {
      set_diff(j, i, 'I'); 
      i--;
      n_i++;
    }
    else if (string2[k] == '-') {
      set_diff(j, i, 'D'); 
      j--;
      n_d++;
    }
    else {
      if (string3[k] == ' ') {
	set_diff(j, i, 'S'); 
	n_s++;
      }
      i--;
      j--;
    }
  }
  set_diff(query_start - 1, subject_start - 1, 'B');

  if (n_i != insertions || n_d != deletions || n_s != mismatches)
    fatalError("Mutation count inconsistency");
  return diffs;
}

set_diff(site1, site2, type)
     int site1, site2;
     char type;
{
  Diff *diff;

  diff = diffs + i_diffs--;
  diff->site1 = site1;
  diff->site2 = site2;
  diff->type = type;
  if (diff->site1 != site1 || diff->site2 != site2)
    fatalError("Diff structure needs redefinition -- ints required");
}

*/

static unsigned int t_n_diffs, t_n_diff_lists;
static unsigned int diff_block_size;
static unsigned char *diff_block;
#define DIFF_BLOCK_SIZE 100000

/* NEED TO CHECK AND COMPLETE THE FOLLOWING BEFORE IMPLEMENTING --
   PURPOSE IS TO PLACE THE LOCATION OF THE DELETION OPTIMALLY.
*/

qual_slide_indels(qual, qual_index)
     int qual_index; /* qual_index = 0 for string1, 1 for string2 */
     char *qual; /* must point to last value preceding alignment */
{
  x_slide_indels(string1, string2, qual, qual_index);
  x_slide_indels(string2, string1, qual, !qual_index);
}

x_slide_indels(string_a, string_b, qual, qual_index)
     unsigned char *string_a, *string_b;
     int qual_index; /* qual_index = 0 for string_a, 1 for string_b */
     char *qual; /* must point to last value preceding alignment */
{
  int j, k, j1;
  char res;

  for (k = 0; k <= string_end; k++) {
    if (!qual_index && string_a[k] != '-') qual++;
    if (qual_index && string_b[k] != '-') qual++;

    if (string_b[k] == '-') {
      res = string_a[k];
      for (j = k - 1; j >= 0 && string_b[j] == res && string_a[j] == res; j--);
      if (j < k - 1 && j >= 0) {
	if (qual[1 + j - k] < qual[0]) {
	  j++;
	  for (j1 = k; j1 > j; j1--) {
	    string_b[j1] = string_b[j1 - 1];
	    string3[j1] = string3[j1 - 1];
	  }
	  string_b[j] = string3[j] = '-';
/*
	  for (j++; j <= k; j++) {
	    string_b[j] = string3[j] = res;
	  }
*/
	}
      }
      else {
	for (j = k + 1; 
	     j <= string_end && string_b[j] == res && string_a[j] == res; 
	     j++);
	if (j > k + 1 && j <= string_end) {
	  if (qual[j - k - 1] < qual[0]) {
	    j--;
	    if (qual_index) qual++;
	    for (j1 = k; j1 < j; j1++) {
	      string_b[j1] = string_b[j1 + 1];
	      string3[j1] = string3[j1 + 1];
	    }
	    string_b[j] = string3[j] = '-';
/*
	    for (j--; j >= k; j--) string_b[j] = string3[j] = res;
*/
	  }
	}
      }
    }
  }
}

unsigned char *make_diffs()
{
  int k, n_diffs;
  int n_i, n_d, n_s, n_m;
  int i_diffs, gap;
  unsigned char *diffs;
  unsigned char *set_diff_block();

  n_i = n_d = n_s = n_m = 0;

  for (k = string_end, gap = 0; k >= 0; k--) {
    if (string1[k] == '-') {
      gap = 0;
      n_i++;
    }
    else {
      gap++;
      if (string2[k] == '-') {
	gap = 0;
	n_d++;
      }
      else {
	if (string3[k] == ' ') {
	  gap = 0;
	  n_s++;
	}
	else if (gap > 60) {
	  gap = 0;
	  n_m++;
	}
      }
    }
  }
  if (n_i != insertions || n_d != deletions || n_s != mismatches) {
    /*
    for (k = string_end; k >= 0; k--) {
      if (string1[k] != string2[k] || string1[k] != string3[k])
	fprintf(stderr, "\n %c %c %c", string1[k], string2[k], string3[k]); 

    }
    */
    fprintf(stderr, "\n\nn_i: %d, n_d: %d, n_s: %d; insertions: %d, deletions: %d, mismatches: %d; string_end: %d\n",
	    n_i, n_d, n_s, insertions, deletions, mismatches, string_end);
    fatalError("Mutation count inconsistency"); 
  }

  n_diffs = n_i + n_d + n_s + n_m;

  diffs = set_diff_block(n_diffs);

/*  diffs = (unsigned char *)our_alloc((n_diffs + 2) * sizeof(char)); */
  t_n_diffs += n_diffs + 2;
  t_n_diff_lists += 1;
  i_diffs = 0;

  for (k = string_end, gap = 0; k >= 0; k--) {
    if (string1[k] == '-') {
      diffs[i_diffs++] = set_diff(gap, 'I');
      gap = 0;
    }
    else {
      gap++;
      if (string2[k] == '-') {
	diffs[i_diffs++] = set_diff(gap, 'D');
	gap = 0;
      }
      else {
	if (string3[k] == ' ') {
	  diffs[i_diffs++] = set_diff(gap, 'S');
	  gap = 0;
	}
	else if (gap > 60) {
	  diffs[i_diffs++] = set_diff(gap, 'M');
	  gap = 0;
	}
      }
    }
  }
  gap++;
  diffs[i_diffs++] = set_diff(gap, 'S'); /* final gap (to end of alignment) */
  /* diffs[i_diffs++] = set_diff(0, 'M'); /* string end */
  return diffs;
}

unsigned char *set_diff_block(n_diffs)
     int n_diffs;
{
  char *our_alloc();
  unsigned char *diffs;

  if (n_diffs + 2 > diff_block_size) { /* new block needed */
    diff_block_size = n_diffs + 2 > DIFF_BLOCK_SIZE ? n_diffs + 2 : DIFF_BLOCK_SIZE;
    diff_block = (unsigned char *)our_alloc(diff_block_size * sizeof(char));
  }
  diffs = diff_block;
  diff_block += n_diffs + 2;
  diff_block_size -= n_diffs + 2;
  diffs[n_diffs + 1] = 0; /* string end */
  return diffs;
}

set_diff(gap, type)
     int gap;
     char type;
{

  if (type == 'm') return 63; /* multiplicity flag -- since gap size 63 never assigned */
  if (type == 'I') return gap + 128;
  if (type == 'D') return gap + 64;
  if (type == 'S') return gap + 192;
  if (type == 'M') return gap;
}

set_rev_diff(diff)
     unsigned char diff;
{
  int t, g;

  t = diff & 192;
  if (!t || t == 192) return (int)diff; /* 'S' and 'M' are unchanged */
  g = diff & 63;
  g = t == 128 ? g + 1 : g - 1;
  t = 192 - t; /* switch 'D' and 'I' */
  return t + g;
}

diff_type(diff)
     unsigned char diff;
{
  int d;

  d = diff;
  if (d == 63) return (int)'m';
  d &= 192;
  if (d == 0) return (int)'M';
  if (d == 192) return (int)'S';
  if (d == 128) return (int)'I';
  if (d == 64) return (int)'D';
  fprintf(stderr, "\nUndefined diff character: %d", d);
  exit(1);
}

rev_diff_type(diff)
     unsigned char diff;
{
  int d;

  d = diff;
  if (d == 63) return (int)'m';
  d &= 192;
  if (d == 0) return (int)'M';
  if (d == 192) return (int)'S';
  if (d == 64) return (int)'I';
  if (d == 128) return (int)'D';
  fprintf(stderr, "\nUndefined diff character: %d", d);
  exit(1);
}

diff_gap1(diff)
     unsigned char diff;
{
  return diff == 63 ? 0 : (int)(diff & 63);
}

diff_gap2(diff)
     unsigned char diff;
{
  int g, t;

  if (diff == 63) return 0;

  g = diff & 63;
  t = diff_type(diff);
  if (t == 'I') g++;
  else if (t == 'D') g--;
  return g;
}
 
get_next_diff(diff, site1, site2, gap1, gap2, type)
     unsigned char **diff;
     int *site1, *site2, *gap1, *gap2, *type;
{
  int g1, g2; 

  *gap1 = *gap2 = 0;
  if (*(*diff + 1)) {
    do {
      (*diff) += 1;
      *site1 += (g1 = diff_gap1(**diff));
      *site2 += (g2 = diff_gap2(**diff));
      *gap1 += g1;
      *gap2 += g2;
      *type = diff_type(**diff);
    } while (*type == 'M');
    return 1;
  }
  else return 0;
}
 
get_prev_diff(diff, site1, site2, type)
     unsigned char **diff;
     int *site1, *site2, *type;
{
  do {
    *site1 -= diff_gap1(**diff);
    *site2 -= diff_gap2(**diff);
    (*diff) -= 1;
    *type = diff_type(**diff);
  } while (*type == 'M');
}
 
/*
get_starts_ends(diffs, n_diffs, pair_start1, pair_start2, pair_end1, pair_end2)
     Diff *diffs;
     int n_diffs;
     int *pair_start1, *pair_start2, *pair_end1, *pair_end2;
{
  *pair_start1 = diffs->site1 + 1;
  *pair_start2 = diffs->site2 + 1;
  *pair_end1 = diffs[n_diffs - 1].site1 - 1;
  *pair_end2 = diffs[n_diffs - 1].site2 - 1;
}
*/
/* return differences to origin 0 
adjust_diffs(diff)
     Diff *diff;
{
    for (; ; diff++) {
      diff->site1 -= 1; 
      diff->site2 -= 1;
      if (diff->type == 'E') break;
    }
}

*/

print_t_n_diffs()
{
  fprintf(stderr, "Total # diffs: %d, in %d lists, size: %.3f Mbytes\n", 
	  t_n_diffs, t_n_diff_lists, t_n_diffs / 1000000.);
}

get_edges(match_pair)
   Aligned_pair *match_pair;
{
  int s_e_8, s_e_4, s_e_2, s_e_3, s_e, k, i_q, i_s;
  char *edge;

  s_e = parameters->splice_edge_length;
  s_e_2 = 2 * s_e;
  s_e_3 = 3 * s_e;
  s_e_4 = 4 * s_e;
  s_e_8 = 8 * s_e;
  edge = match_pair->query_data->edges12; /* = (char *)our_alloc((s_e_8 + 1) * sizeof(char)); */

  for (k = 0; k < s_e; k++) {
    i_q = query_start - 1 - s_e + k;
    i_s = subject_start - 1 - s_e + k;
    edge[k] = i_q >= 0 ? q_profile->seq[i_q] : '-';
    edge[k + s_e_4] = i_s >= 0 ? db_entry_seq[i_s] : '-';
  }

  /* strings are in reverse order!! */
  for (; k < s_e_2; k++) {
    edge[k] = string1[string_end - k + s_e];
    edge[k + s_e_4] = string2[string_end - k + s_e];
  }

  for (; k < s_e_3; k++) {
    edge[k] = string1[s_e_3 - 1 - k];
    edge[k + s_e_4] = string2[s_e_3 - 1 - k];
  }

  for (; k < s_e_4; k++) {
    i_q =  query_end + k - s_e_3;
    i_s = subject_end + k - s_e_3;
    edge[k] = i_q < q_profile->length ? q_profile->seq[i_q] : '-';
    edge[k + s_e_4] = i_s < db_entry_length ? db_entry_seq[i_s] : '-';
  }

  edge[s_e_8] = 0;

}    
