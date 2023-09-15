/*****************************************************************************
#   Copyright (C) 1993-1996 by Phil Green.                          
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

/* Needleman-Wunsch version of smith_wat.c.  */

/*
   If function compiled with the flag -DFINDALIGN, alignments will be printed out;
   otherwise, only the score is computed and returned (latter version is
   appropriate for database searches, as it is faster).
*/

#include "swat.h"
#define BIGNEGATIVE -32000

int *maxstu_vec;
int maxstu_alloc_size;

#ifdef FINDALIGN

int full_nw(q_profile, db_entry_seq, db_entry_length, print_alignment, query_start, query_end,
			subject_start, subject_end)
     int print_alignment;
     int *query_start, *query_end, *subject_start, *subject_end;
#else
int nw(q_profile, db_entry_seq, db_entry_length)
#endif
     Profile *q_profile;
     char *db_entry_seq;
     int db_entry_length;
{
  int gap_init;
  char *seq;
  register int maxstu_above_el, u, t, maxstu, gap_ext,
     temp, max_score;
  register SCORE *score_vec;
  register int *maxstu_above, *t_vec, *t_above;
  int q_length, s_length;
  int qbegin_gap_init,qbegin_gap_ext,qend_gap_init,qend_gap_ext,
      sbegin_gap_init,sbegin_gap_ext,send_gap_init,send_gap_ext;
  int lead_penalty, save_t, min_entry;

#ifdef FINDALIGN
  char *string1, *string2, *string3;
  char **direc, **dir, **mdir;
  int i, j, k;
  int s, score_check;
  register char *diri,*mdiri;
  register char d;
  char u_ext, t_ext; /* flags for whether current u and t are
			     gap extending */
  int alen; /* alignment length */
  char *our_alloc();
  int last_q, last_s, g_ext, g_init;
  int first_s, first_q, first_alen;
#endif

  seq = db_entry_seq;
  q_length = q_profile->length;

  s_length = db_entry_length;

#ifdef FINDALIGN
  dir = direc = (char **)our_alloc(s_length * sizeof(char *));
  string1 = (char *)our_alloc((q_length + s_length) * sizeof(char));
  string2 = (char *)our_alloc((q_length + s_length) * sizeof(char));
  string3 = (char *)our_alloc((q_length + s_length) * sizeof(char));
#endif
/* Assume score matrix has been adjusted by subtracting 2 * gap_ext from each
   entry; in following algorithm then can replace gap_init by gap_init - gap_ext,
   and gap_ext by 0, and readjust total score at end */

  gap_ext = q_profile->gap_ext;
  gap_init = q_profile->gap_init - gap_ext;

  qbegin_gap_init = q_profile->qbegin_gap_init - gap_ext;
  qbegin_gap_ext = q_profile->qbegin_gap_ext - gap_ext;
  qend_gap_init = q_profile->qend_gap_init - gap_ext;
  qend_gap_ext = q_profile->qend_gap_ext - gap_ext;

  sbegin_gap_init = q_profile->sbegin_gap_init - gap_ext;
  sbegin_gap_ext = q_profile->sbegin_gap_ext - gap_ext;
  send_gap_init = q_profile->send_gap_init - gap_ext;
  send_gap_ext = q_profile->send_gap_ext - gap_ext;

  /*  maxstu_vec = q_profile->maxstu_vec; NEED TO CHECK ALLOCATION ETC. -- LIKE IN smith_wat.c */
  t_vec = q_profile->t_vec;

  min_entry = q_profile->min_entry - 1; /* maximum negative value in score matrix */
  lead_penalty = sbegin_gap_init;

  /* FOLLOWING PROB. NOW INCORRECT -- SINCE maxstu_vec NO LONGER INITIALIZED IN profile.c */
  for (maxstu_above = maxstu_vec, t_above = t_vec;
       *maxstu_above >= *t_above; 
      maxstu_above++, t_above++) {
    *maxstu_above = lead_penalty;
    *t_above = BIGNEGATIVE; /* to ensure is always too small to
						affect first row calculations */
    lead_penalty += sbegin_gap_ext;
  }
  save_t = BIGNEGATIVE;
  maxstu_above_el = 0;
  u = BIGNEGATIVE;
  lead_penalty = qbegin_gap_init;
  for (; *seq; seq++) {
#ifdef FINDALIGN
    diri = *dir = (char *)our_alloc(q_length * sizeof(char));
    u_ext = 0;
#endif
    score_vec = q_profile->scores[q_profile->convert[*seq]];
    maxstu_above = maxstu_vec;
    t_above = t_vec;
    if (*(seq + 1)) {
      for(;;) {
	maxstu = maxstu_above_el + *score_vec++;
	if ((t = *t_above) > (maxstu_above_el = *maxstu_above)) {
     /* end of row has been reached; recalculate t values for last column entry */
	  temp = *(--maxstu_above) + qend_gap_init;
	  t = save_t + qend_gap_ext;
#ifdef FINDALIGN
	  d = *(--diri);
#endif
	  if (t <= temp) {
	    t = temp;
#ifdef FINDALIGN
	    if (d & 4) *diri = d -= 4;
#endif
	  }
#ifdef FINDALIGN
	  else {
	    if (!(d & 4)) *diri = d += 4;
	  }
#endif
	  *(--t_above) = save_t = t;
	  *maxstu_above = t; /* to retain sentinel */
	  
	  break;
	}
#ifdef FINDALIGN
	d = 3; 
#endif
	if (maxstu < t) { 
	  maxstu = t; 
#ifdef FINDALIGN
	  d = 1;
#endif
	} 
	if (maxstu < u) {
	  maxstu = u; 
#ifdef FINDALIGN
	  d = 2; 
#endif
	}
#ifdef FINDALIGN
	if (maxstu_above_el + gap_init != t) 
	  d += 4; /* set 4's bit -- to indicate that extending gap */
	if (u_ext) d += 8; /* set 8's bit -- to indicate that extending gap */
	u_ext = 1;
#endif
	temp = maxstu + gap_init; 
	if (u < temp) {
	  u = temp;
#ifdef FINDALIGN
	  u_ext = 0; 
#endif
	}
	if (t < temp) *t_above = temp; 
	
#ifdef FINDALIGN
	*diri++ = d; 
#endif
	*maxstu_above++ = maxstu;
	t_above++;
      }
      maxstu_above_el = lead_penalty;
      u = BIGNEGATIVE;
      lead_penalty += qbegin_gap_ext;
#ifdef FINDALIGN
      dir++;
#endif
    }
    else {
      for(;;) {
	maxstu = maxstu_above_el + *score_vec++;
	if ((t = *t_above) > (maxstu_above_el = *maxstu_above)) {
     /* end of row has been reached; recalculate t values for last column entry */
	  maxstu_above--;
#ifdef FINDALIGN
	  temp = *maxstu_above + qend_gap_init;
	  t = save_t + qend_gap_ext;
	  d = *(--diri);
	  if (t > temp) {
	    if (!(d & 4)) *diri = d += 4;
	  }
	  else {
	    if (d & 4) *diri = d -= 4;
	  }
#endif
	  break;
	}
#ifdef FINDALIGN
	d = 3; 
#endif
	if (maxstu < t) { 
	  maxstu = t; 
#ifdef FINDALIGN
	  d = 1;
#endif
	} 
	if (maxstu < u) {
	  maxstu = u; 
#ifdef FINDALIGN
	  d = 2; 
#endif
	}
#ifdef FINDALIGN
	if (maxstu_above_el + gap_init != t) 
	  d += 4; /* set 4's bit -- to indicate that extending gap */
	if (u_ext) d += 8; /* set 8's bit -- to indicate that extending gap */
	u_ext = 1;
#endif
	temp = maxstu + send_gap_init; 
	u += send_gap_ext;

	if (u < temp) {
	  u = temp;
#ifdef FINDALIGN
	  u_ext = 0; 
#endif
	}
#ifdef FINDALIGN
	*diri++ = d; 
#endif
	*maxstu_above++ = maxstu;
	t_above++;
      }
    }
  }
  /* Readjust score to compensate for prev use of gap_ext = 0. */
  max_score = *maxstu_above + (s_length + q_length) * gap_ext;
/*  printf("\n MAXSCORE %d %d %d %d %d %d\n",max_score, *maxstu_above,
	 *(maxstu_above - 1), s_length, q_length, gap_ext);
*/						  
#ifdef FINDALIGN
  u_ext = t_ext = 0;
  score_check = 0;

  last_q = last_s = first_q = first_s = 0;
  *subject_end = s_length;
  *query_end = q_length;
  for (i = s_length - 1, j =  q_length - 1, k = 0;
       i >= 0 || j >= 0; k++) {
    
    string2[k] = string1[k] = '-';
    if (i < 0) d = 2;
    else if (j < 0) d = 1;
    else {
      d = direc[i][j]; 
      
      if (t_ext) {
	d = 1 + (d & 28);
      }
      if (u_ext) {
	d = 2 + (d & 28);
      }
      if (!(d & 2)) t_ext = d & 4;
      if (!(d & 1)) u_ext = d & 8;
    }

    if (d & 1) { 
      string2[k] = db_entry_seq[i];
      last_s = k;
      if (!first_s) first_s = k + 1;
      i--;
    }
    if (d & 2) {
      string1[k] = q_profile->seq[j];
      last_q = k;
      if (!first_q) first_q = k + 1;
      j--;
    }
    s = 0;

    if (string1[k] == '-') {
      score_check += (!k || string1[k - 1] != '-') ? gap_init : 0;
    }
    else if (string2[k] == '-') {
      score_check += (!k || string2[k - 1] != '-') ? gap_init : 0;
    }
    else /* (string2[k] != '-' && string1[k] != '-') */ {
      s = q_profile->scores[q_profile->convert[string2[k]]][j + 1];
      score_check += s;
      s += 2 * gap_ext;
    }
    string3[k] = string2[k] == string1[k] ? string1[k] : (s > 0 ? '+' : ' ');
  }

  score_check += (q_length + s_length) * gap_ext;
  k--;
  first_q--;
  first_s--;

  if (last_q < last_s) {
    alen = last_q;
    i += k - alen;
    score_check +=
      qbegin_gap_init - gap_init + (k - last_q - 1) * qbegin_gap_ext;
  }
  else {
    alen = last_s;
    j += k - alen;
    if (last_s < last_q) score_check +=
      sbegin_gap_init - gap_init + (k - last_s - 1) * sbegin_gap_ext;
  }
  if (first_q < first_s) {
    first_alen = first_s;
    score_check +=
      send_gap_init - gap_init + (first_s - 1) * send_gap_ext;
  }
  else {
    first_alen = first_q;
    if (first_q) score_check +=
      qend_gap_init - gap_init + (first_q - 1) * qend_gap_ext;
  }

  if (score_check != max_score)
    printf("\nSCORE DISCREPANCY: initial %d, from alignment %d\n",
	   max_score, score_check);
  i += 2;
  j += 2;
  *subject_start = i;
  *query_start = j;

  if (print_alignment) {
    for (; alen >= first_alen; alen -= 50) {
      printf("\n\n%5d ",i);
      for (k = alen; k > alen - 50 && k >= first_alen; k--) {
	printf("%c",string2[k]);
	if (string2[k] != '-') i++;
      }
      printf(" %d\n      ", i - 1);
      for (k = alen; k > alen - 50 && k >= first_alen; k--) printf("%c",string3[k]);
      printf("\n%5d ", j);
      for (k = alen; k > alen - 50 && k >= first_alen; k--) {
	printf("%c",string1[k]);
	if (string1[k] != '-') j++;
      }
      printf(" %d", j - 1);
    }
  }
  for (i = 0; i < s_length; i++) our_free(direc[i]);
  our_free(direc);
  our_free(string1);
  our_free(string2);
  our_free(string3);
#endif

/* now subtract approximate end gap penalty due to discrepancy in
   sequence lengths */

  if (q_length > s_length) {
    max_score -= (q_length - s_length) *
      (sbegin_gap_ext < send_gap_ext ? q_profile->sbegin_gap_ext : q_profile->send_gap_ext);
  }
  else
    max_score -= (s_length - q_length) *
      (qbegin_gap_ext < qend_gap_ext ? q_profile->qbegin_gap_ext : q_profile->qend_gap_ext);
    
  return max_score; 
  printf("\n MAXSCORE %d\n",max_score);
}
