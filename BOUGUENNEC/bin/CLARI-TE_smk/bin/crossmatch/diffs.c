/*****************************************************************************
#   Copyright (C) 1994-1996 by Phil Green.                          
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
/*
extern int **c_score_mat;
extern int row_convert[], col_convert[];
*/
extern int merge_residues[256];


unsigned char *make_reversed_diffs(diffs, rev)
     unsigned char *diffs;
     int rev;
{
  char *our_alloc();     
  unsigned char *diff, *diff2, *diffs2;
  int n_diffs;
  int d, d_next, g, last_g, largegap_flag, next_diff, last_diff;
  unsigned char *set_diff_block();

  if (!diffs) return (unsigned char *)0;
  n_diffs = strlen((char *)diffs) + 1; 
  diffs2 = set_diff_block(n_diffs - 2);
  /*
  diffs2 = (unsigned char *)our_alloc(n_diffs * sizeof(char));
  */
  if (rev) {
    /* diffs2[n_diffs - 1] = 0; */
    d = 'S';
    for (diff = diffs, diff2 = diffs2 + n_diffs - 2, largegap_flag = 0; *diff; diff2--, diff++) {
      next_diff = *diff;
      d_next = rev_diff_type(next_diff); /* d and d_next are the types of discrepancies
				    which flank the interval */
      if (d == 'm') {
	largegap_flag = !largegap_flag;
	*diff2 = last_diff;
      }
      else if (d == 'M' && largegap_flag) {
	*diff2 = last_diff;
      }
      else {
	if (d_next == 'm') {
	  g = d == 'D'; /* d == 'I' only other relevant case */
	}
	else {
	  g = diff_gap2(next_diff);
	  if (d_next == 'I') g++;
	  if (d == 'I') g--;
	}
	*diff2 = set_diff(g, d);
      }
      d = d_next;
      last_diff = next_diff;
    }
  }
  else {
    for (diff = diffs, diff2 = diffs2; *diff; diff2++, diff++) {
      *diff2 = set_rev_diff(*diff);
      /*
      d = diff_type(*diff);
      g = diff_gap2(*diff);
      if (d == 'I') *diff2 = set_diff(g, 'D');
      else if (d == 'D') *diff2 = set_diff(g, 'I');
      else *diff2 = *diff;
      */
    }
    /* *diff2 = 0; */
  }
  return diffs2;
}
/*
unsigned char *make_reversed_diffs(pair, diffs)
     Aligned_pair *pair;
     unsigned char *diffs;
{
  char *our_alloc();     
  Diff *diff, *diff2, *pair2_diffs;
  int i, length1, length2;
  int n_diffs;

  n_diffs = strlen(diffs) + 1;
  diff2 = pair2_diffs = (Diff *)our_alloc(n_diffs * sizeof(Diff));
  if (is_reverse(pair)) {
    length1 = get_seq_length(pair->entry1) - 1;
    length2 = get_seq_length(pair->entry2) - 1;
    diff = diffs + n_diffs - 1;
    for (i = 0; i < n_diffs; i++, diff2++, diff--) {
      diff2->site1 = length2 - diff->site2;
      diff2->site2 = length1 - diff->site1;
      if (diff->type == 'I') {
	diff2->type = 'D';
	diff2->site2 -= 1;
      }
      else if (diff->type == 'D') {
	diff2->type = 'I';
	diff2->site1 -= 1;
      }
      else diff2->type = diff->type;
    }
    if ((diff2 - 1)->type != 'B') fatalError("Diff count");
    (diff2 - 1)->type = 'E';
    if (pair2_diffs->type != 'E') fatalError("Diff count");
    pair2_diffs->type = 'B';
  }
  else {
    diff = diffs;
    for (i = 0; i < n_diffs; i++, diff2++, diff++) {
      diff2->site1 = diff->site2;
      diff2->site2 = diff->site1;
      if (diff->type == 'I') diff2->type = 'D';
      else if (diff->type == 'D') diff2->type = 'I';
      else diff2->type = diff->type;
    }
    if ((diff2 - 1)->type != 'E') fatalError("Diff count");
    if (pair2_diffs->type != 'B') fatalError("Diff count");
  }
  return pair2_diffs;
}
*/
extern int qseq_trans[256];
print_alignment_from_diffs(pair)
     Aligned_pair *pair;
{
  int i, j, n, s, cc, start_i, start_j;
  int g_init, g_tot, last_i_init, last_j_init;
  double transitions, transversions, length;
  unsigned char *diff_i, *diff_j, *start_diff_i, *start_diff_j;
  char print_char;
  int site_j, site_i, start_site_j, start_site_i, c_1, c_2;
  char *get_id();
  unsigned char *get_seq(), *get_comp_seq();
  char *id1, *id2;
  unsigned char *seq1, *seq2;
  int length1, length2, mark_mode, type_j, type_i, largegap_flag_j, largegap_flag_i, largegap_flag_j_start, largegap_flag_i_start;

  int left_extend, right1, right2, right_extend, start1, start2, end1, end2;

  /* N.B. This could be a problem for cross_match -- may need to construct complement! */

  mark_mode = parameters->DNA_flag; /* mode = 0: indicating matching residues; 1 : mark transitions, transversions
				       2: merged base reads -- indicate agreement/disagreement according to score */

  if (mark_mode >= 3) mark_mode = 1;

  id1 = get_id(pair->entry1);
  id2 = get_id(pair->entry2);
  seq1 = get_seq(pair->entry1);
  seq2 = is_reverse(pair) ? get_comp_seq(pair->entry2) : get_seq(pair->entry2);
  length1 = get_seq_length(pair->entry1);
  length2 = get_seq_length(pair->entry2);

  g_init = g_tot = 0;
  last_i_init = last_j_init = -1;
  transitions = transversions = 0;


  /* fprintf(stderr, " %s %s", id1, id2); */

  /*
  printf("\n");
  for (diff_i = pair->diffs; *diff_i; diff_i++)
    printf(" %c-%d", diff_type(*diff_i), diff_gap1(*diff_i));
  printf("\n");
  */

  if (parameters->align_extend) {
    left_extend = pair->start1 < pair->start2 ? pair->start1 : pair->start2;
    left_extend = parameters->align_extend < left_extend ? parameters->align_extend : left_extend;
    right2 = length2 - 1 - pair->end2;
    right1 = length1 - 1 - pair->end1;
    right_extend = right1 < right2 ? right1 : right2;
    right_extend = parameters->align_extend < right_extend ? parameters->align_extend : right_extend;
  }
  else left_extend = right_extend = 0;

  start1 = pair->start1 - left_extend;
  start2 = pair->start2 - left_extend;

  end1 = pair->end1 + right_extend;
  end2 = pair->end2 + right_extend;

  i = start2;
  j = start1; 
  diff_i = diff_j = pair->diffs;
  site_j = pair->start1 - 1 + diff_gap1(*diff_j);
  site_i = pair->start2 - 1 + diff_gap2(*diff_i);
  largegap_flag_j = largegap_flag_i = 0;  /* 0: not in large gap; 'D' in large D gap, 'I' in large I gap */

  for ( ; i <= end2 && j <= end1; ) {
    start_i = i;
    start_j = j;
    start_diff_i = diff_i;
    start_diff_j = diff_j;
    start_site_j = site_j;
    start_site_i = site_i;
    largegap_flag_j_start = largegap_flag_j;
    largegap_flag_i_start = largegap_flag_i;
    type_j = diff_type(*diff_j);
    type_i = diff_type(*diff_i);

    printf("\n\n%c %-15.15s %9d ", is_reverse(pair) ? 'C' : ' ', id2, 
	   is_reverse(pair) ? length2 - i : i + 1);
    for (cc = 0; cc < 50 && i <= end2; cc++) {
      while (*(diff_i + 1) && (site_i < i - 1 || type_i != 'D' && !largegap_flag_i && type_i != 'm')) {
	diff_i++;
	type_i = diff_type(*diff_i);
	if (type_i == 'm') 
	  largegap_flag_i = diff_type(*(diff_i - 1)); /* must be entering gap */
	else 
	  site_i += diff_gap2(*diff_i);
      }
      if (i <= site_i || i > pair->end2) printf("%c", i < pair->start2 || i > pair->end2 ? tolower(qseq_trans[seq2[i++]]) : qseq_trans[seq2[i++]]);
      else {
	if (largegap_flag_i || type_i == 'm') {
	  printf(".");
	}
	else {
	  g_tot++;
	  if (i != last_i_init) {
	    g_init++;
	    last_i_init = i;
	  }
	  printf("-");
	}
	diff_i++;
	type_i = diff_type(*diff_i);
	if (type_i == 'm') {
	  largegap_flag_i = largegap_flag_i ? 0 : diff_type(*(diff_i - 1));
	}
	else {
	  if (!largegap_flag_i)
	    site_i += diff_gap2(*diff_i);
	  else if (largegap_flag_i == 'I') {
	    site_i += diff_gap2(*diff_i);
	    i = site_i + 1;
	  }
	}
      }
    }
    printf(" %d\n%15.15s             ", is_reverse(pair) ? length2 - i + 1 : i, "");

    i = start_i;
    diff_i = start_diff_i;
    largegap_flag_i = largegap_flag_i_start;
    type_i = diff_type(*diff_i);
    site_i = start_site_i;

    for (cc = 0; cc < 50 && i <= end2; cc++) {
      while (*(diff_i + 1) && (site_i < i - 1 || type_i != 'D' && !largegap_flag_i && type_i != 'm'))  {
	diff_i++;
	type_i = diff_type(*diff_i);
	if (type_i == 'm') 
	  largegap_flag_i = diff_type(*(diff_i - 1)); /* must be entering gap */
	else 
	  site_i += diff_gap2(*diff_i);
      }
      while (*(diff_j + 1) && (site_j < j - 1 || type_j != 'I' && !largegap_flag_j && type_j != 'm'))  {
	diff_j++;
	type_j = diff_type(*diff_j);
	if (type_j == 'm') 
	  largegap_flag_j = diff_type(*(diff_j - 1));
	else 
	  site_j += diff_gap1(*diff_j);
      }

      if (i <= site_i && j <= site_j || i > pair->end2) {
	c_2 = toupper(qseq_trans[seq2[i]]);
	c_1 = toupper(seq1[j]);
	if (mark_mode == 0) 
	  print_char = c_2 == c_1 ? qseq_trans[seq2[i]] : ' '; 
	else if (mark_mode == 1) {
	  print_char = ' ';
	  
	  if (c_2 != c_1) {
	    if (!strchr("ACGTYR", c_1) || !strchr("ACGTYR", c_2))
	      print_char = '?';
	    else if (strchr("CTY", c_1) && strchr("CTY", c_2)
		     || strchr("AGR", c_1) && strchr("AGR", c_2)) {
	      if (i >= pair->start2 && i <= pair->end2)
		transitions++;
	      print_char = 'i';
	    }
	    else {
	      if (i >= pair->start2 && i <= pair->end2)
		transversions++;
	      print_char = 'v';
	    }
	  }
	}
	else {
	  print_char = ' ';
	  if (c_2 != c_1) {
	    /*
	    s = c_score_mat[row_convert[seq1[j]]][col_convert[seq2[i]]];
	    print_char = s > 0 ? '+' : s == 0 ? '0' : 'x';
	    */
	    s = merge_residues[seq1[j]] & merge_residues[qseq_trans[seq2[i]]];
	    print_char = !s ? 'X' : s & 15 ? '2' : '1';

	  }
	}
	i++;
	j++;
      }
      else {
	print_char = largegap_flag_i || type_i == 'm' ? '.' : '-';
	if (j > site_j) {
	  diff_j++;
	  type_j = diff_type(*diff_j);
	  if (type_j == 'm') {
	    largegap_flag_j = largegap_flag_j ? 0 : diff_type(*(diff_j - 1));
	  }
	  else {
	    if (!largegap_flag_j)
	      site_j += diff_gap1(*diff_j);
	    else if (largegap_flag_j == 'D') {
	      site_j += diff_gap1(*diff_j);
	      j = site_j + 1;
	    }
	  }
	}
	else j++;
	if (i > site_i) {
	  diff_i++;
	  type_i = diff_type(*diff_i);
	  if (type_i == 'm') {
	    largegap_flag_i = largegap_flag_i ? 0 : diff_type(*(diff_i - 1));
	  }
	  else {
	    if (!largegap_flag_i)
	      site_i += diff_gap2(*diff_i);
	    else if (largegap_flag_i == 'I') {
	      site_i += diff_gap2(*diff_i);
	      i = site_i + 1;
	    }
	  }
	}
	else i++;
      }
      printf("%c", print_char);
    }

    j = start_j;
    diff_j = start_diff_j;
    largegap_flag_j = largegap_flag_j_start;
    type_j = diff_type(*diff_j);
    site_j = start_site_j;

    /*    printf(" %d %d %c", j, site_j, type_j); */

    printf("\n  %-15.15s %9d ", id1, j + 1);
    for (cc = 0; cc < 50 && j <= end1; cc++) {
      while (*(diff_j + 1) && (site_j < j - 1 || type_j != 'I' && !largegap_flag_j && type_j != 'm'))  {
	diff_j++;
	type_j = diff_type(*diff_j);
	if (type_j == 'm') 
	  largegap_flag_j = diff_type(*(diff_j - 1));
	else 
	  site_j += diff_gap1(*diff_j);
      }
      if (j <= site_j || j > pair->end1) printf("%c", j < pair->start1 || j > pair->end1 ? tolower(seq1[j++]) : seq1[j++]); 
      else {
	if (largegap_flag_j || type_j == 'm') {
	  printf(".");
	}
	else {
	  g_tot++;
	  if (j != last_j_init) {
	    g_init++;
	    last_j_init = j;
	  }
	  printf("-");
	}
	diff_j++;
	type_j = diff_type(*diff_j);
	if (type_j == 'm') {
	  largegap_flag_j = largegap_flag_j ? 0 : diff_type(*(diff_j - 1));
	}
	else {
	  if (!largegap_flag_j)
	    site_j += diff_gap1(*diff_j);
	  else if (largegap_flag_j == 'D') {
	    site_j += diff_gap1(*diff_j);
	    j = site_j + 1;
	  }
	}
      }
    }
    printf(" %d", j);
  }
  if (mark_mode == 1) {
    length = pair->end2 - pair->start2 + 1;
    printf("\n\nTransitions / transversions = %.2f (%.0f / %.0f); ", transversions ? transitions / transversions : 1.0,
	   transitions, transversions);
    printf("Gap_init rate = %.2f (%d / %.0f), avg. gap size = %.2f (%d / %d)", 
	   g_init / length, g_init, length, g_init ? g_tot / (double)g_init : 0.0, g_tot, g_init);
  }
}

