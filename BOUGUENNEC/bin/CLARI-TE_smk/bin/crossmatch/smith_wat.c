/*****************************************************************************
#   Copyright (C) 1993-1998, 2006 by Phil Green.                          
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

/* modified 060725  -- to remove inefficiencies in banded case 
   these change origins of query and subject */


/* Efficient implementation of Smith-Waterman algorithm */

/*
   If function compiled with the flag -DFINDALIGN, alignments will be printed out;
   otherwise, only the score is computed and returned (latter version is
   appropriate for database searches, as it is faster).
   The defs of the key variables (using matrix notation, altho the matrix
   itself is not actually stored) are as follows. We imagine the columns
   of the matrix (indexed by j) as corresponding to the elements of the
   query sequence (or profile), while the rows correspond to the elements
   of the subject (database) sequence.
   s_ij is the score of the best alignment that ends at positions i (in subject)
   and j (in the query), AND that aligns i with j:  y_j
                                                    x_i
   t_ij is score of best alignment ending at i and j that has i at the
     end of the alignment, opposite a gap in the query.    y_j   -
                                                                x_i
   u_ij is score of best alignment ending at i and j that has j at the end of the
     alignment, opposite a gap in the subject.       y_j
                                                 x_i  -

   Then s_ij = max( max (s_{i-1,j-1},t_{i-1,j-1},u_{i-1,j-1}) + score(x_i,y_j), 0)
     t_ij = max(s_{i-1,j} + gap_init, t_{i-1,j} + gap_ext, u_{i-1,j} + gap_init)
     u_ij = max(s_{i,j-1} + gap_init, t_{i,j-1} + gap_init, u_{i,j-1} + gap_ext) 
   (Note: don't need to max t and u with 0 (?) )
   The following efficiency simplifications are made:
     (i) 
         store max (s,t,u) + gap_init, rather than s or max (s,t,u).
         Formula for t_ij is then max(maxstu_{i-1,j}, t_{i-1,j} + gap_ext);
	 similarly for u_ij. This saves an addition or two. -- NB: in current version
	 this no longer valid. Instead only initiate new t and u when max + gap_init > 0.
     (ii) only store a single row at a time of maxstu; maxstu_ij in the
         rightmost 16 bits of each entry in this vector. t_i+1 j is stored
         in the leftmost 16 bits and only extracted when necessary (since t is
	 usually 0, this is inexpensive).
         Only need a single value of u at a time, so don't need to
	 store it in a vector. 
     (iii) entries in score(x_i, y_j) are accessed as vectors (in j) which
         have all scores for a fixed residue (the current x_i) against residues
	 in j. This reduces matrix access overhead, and has the side benefit that
	 one can use exactly the same function when the query is a profile
	 instead of a sequence (altho in that case one might want to
	 force a global alignment rather than a local one); the variable names
	 reflects this. 
      (iv) the inner loop is coded quite efficiently, so that only 8
         instructions are carried out in most iterations (at least when gap_init is
	 reasonably large, relative to score matrix entries). These are:
	   2 memory retrieves  (maxstu and score).
	   2 pointer increments (""       "").
	   1 addition: maxstu + score.
	   3 tests: (next stored) maxstu < 0, and (current) maxstu > 0.
             1 additional test is carried out when maxstu > 0.
	     the first test subsumes the tests for the end of the current row and
	     for t being positive.
	   1 memory store (maxstu).
*/
/* 5/11/94: added changes to allow banded smith-waterman search.
   band_width: the width of the band. Generally an odd no. (so that band
   is symmetric around a diagonal), or 0 (original non-banded search).
   offset defines the central diagonal; it is
   (index of position in query) - (index of corresponding position in
   subject).
   5/17/94: allow for scores exceeding 2^15 - 1, by checking for max-score
     exceeding 32000, and (in quick version) returning immediately, or (in full
     version) reducing all scores in row, and keeping track of amount of
     reduction, and adding back later.
   5/23/94: in full version, return # mismatches, insertions, deletions (relative
     to query -- i.e. an insertion is an insertion in the subject, or a deletion
     in the query); "mismatch" is a negative score.
   6/1/94: Changed treatment of banded searches: now pass l_edge and r_edge
     instead of band_width and offset. For non-banded search, put r_edge
     < l_edge.
 */
#include "swat.h"

int *maxstu_vec;
int maxstu_alloc_size;

#ifdef QUICKALIGN
static int n_maxes, max_list_len;
static Max_list *max_list;
#endif

#if defined(FINDALIGN)
int full_smith_waterman(q_profile, db_entry_seq, db_entry_length, l_edge, r_edge,
			q_left, q_right, s_left, s_right, minscore, add_orig_score)
     int *add_orig_score;
     int minscore;
#elif defined(QUICKALIGN)
int quick_full_smith_waterman(q_profile, db_entry_seq, db_entry_length, l_edge, r_edge, 
		   q_left, q_right, s_left, s_right, minscore, add_orig_score, add_success)
     int *add_orig_score, *add_success;
     int minscore;
#elif defined(COUNTS)

double n1, n2, n3, n4, n5, n6, n7, n8, n9, 
   n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, 
   n20, n21, n22, n23, n24, n25, n26, n27, n28, n29;

int smith_waterman_counts(q_profile, db_entry_seq, db_entry_length, l_edge, r_edge, 
		   q_left, q_right, s_left, s_right)
#else 
int smith_waterman(q_profile, db_entry_seq, db_entry_length, l_edge, r_edge, 
		   q_left, q_right, s_left, s_right)
#endif



     Profile *q_profile;
     unsigned char *db_entry_seq;
     int db_entry_length;
     int l_edge, r_edge; /* left and right edges of band (in first row) --
			    origin 0 */
     int q_left, q_right, s_left, s_right; /* boundaries of the query and 
					      subject sequences to be used 
					      in this search (origin 1) */
{
  int band_width, r_edge_orig, l_edge_orig, saved_pos, base_l_edge;
  int i, gap_init, q_len, s_len;
  char s_replace;
  unsigned char *seq;
  register int maxstu_above_el, u, t, maxstu, t_gap_ext, u_gap_ext,
     temp, max_score, n_gap_init;
  register SCORE *score_vec;
  register int *maxstu_above;
  char *our_alloc();

#ifdef QUICKALIGN
  int q_pos, s_pos;
  int q_left_orig, s_left_orig, start1, start2, end1, end2, temp2;
#endif

#ifdef FINDALIGN
  char **direc, **dir, **mdir;
  register char *diri,*mdiri;
  register char d;
  char u_ext; /* flags for whether current u and t are
			     gap extending */
  int alloc_size;
#endif

#if defined(FINDALIGN) || defined(QUICKALIGN)
  char excess_flag;
  int excess, score;
  char *our_alloc();

  *add_orig_score = 0;
#endif

#ifdef QUICKALIGN
  q_left_orig = q_left;
  s_left_orig = s_left;

  *add_success = 1;
#endif

  l_edge_orig = l_edge;
  r_edge_orig = r_edge;
  s_len = (s_right <= 0 || s_right > db_entry_length) ? 
    db_entry_length : s_right; 

  q_len = (q_right <= 0 || q_right > q_profile->length) ?
    q_profile->length : q_right; 

  /* return s_left, q_left to origin 0 */

  s_left = (s_left > 1 && s_left <= s_len) ? s_left - 1 : 0;
  q_left = (q_left > 1 && q_left <= q_len) ? q_left - 1 : 0;
    
  band_width = r_edge - l_edge + 1;
  if (band_width <= 0) {
    band_width = 0;
    l_edge = q_left;
    r_edge = q_len - 1;
  }
  else {
    if (l_edge + s_len > q_len) s_len = q_len - l_edge;
    if (r_edge + s_left < q_left) s_left = q_left - r_edge;

    l_edge += s_left;
    if (q_left < l_edge) q_left = l_edge;
    if (q_len > r_edge + s_len) q_len = r_edge + s_len;

    r_edge += s_left;

    if (r_edge > q_len - 1) r_edge = q_len - 1;
  }

  if (s_left >= s_len || q_left >= q_len)
      return 0; /* band does not intersect rectangle */

  seq = db_entry_seq;
  saved_pos = s_len;
  s_replace = seq[saved_pos];
  seq[saved_pos] = 0; /* create artificial end for seq; replace later */
  seq += s_left;

  max_score = 0; 

  gap_init = q_profile->gap_init;
  t_gap_ext = q_profile->ins_gap_ext;
  u_gap_ext = q_profile->del_gap_ext;
  n_gap_init = -gap_init;
  /* maxstu_vec = q_profile->maxstu_vec; */

  init_alignment_globals(q_profile, db_entry_seq, db_entry_length, band_width, r_edge_orig, l_edge_orig, 
			 q_left, s_left);

  /* reset so q_left, s_left are at 0 */

  l_edge -= q_left;
  r_edge -= q_left;
  q_len -= q_left;
  s_len -= s_left;

  base_l_edge = l_edge;

#ifdef QUICKALIGN
  if (s_len > max_list_len) {
    if (max_list_len) our_free(max_list);
    max_list = (Max_list *)our_alloc(s_len * sizeof(Max_list));
    max_list_len = s_len;
    /* fprintf(stderr, " s_len alloc: %d", s_len); */
  }
  n_maxes = 0;
#endif

#ifdef FINDALIGN
  dir = direc = (char **)our_alloc(s_len * sizeof(char *));
  /* fprintf(stderr, " s_len alloc: %d", s_len); */
  direc -= s_left; 
  alloc_size = band_width ? band_width + 2 : q_len;
#endif

#if defined(FINDALIGN) || defined(QUICKALIGN)
  excess = excess_flag = 0;
#endif

  if (maxstu_alloc_size < q_len + 2) {
    /*    fprintf(stderr," maxstu_vec realloc %d => %d (%d,%d,%d)", maxstu_alloc_size, i + 2, q_len, s_len, r_edge); */
    if (maxstu_vec) our_free(maxstu_vec);
    maxstu_vec = (int *)our_alloc((q_len + 10) * sizeof(int));
    /* fprintf(stderr, " q_len alloc: %d", q_len); */
    maxstu_alloc_size = q_len + 10;
  }

  for (i = 0; i <= r_edge; i++) maxstu_vec[i] = 0;
  maxstu_vec[i] = -1;

  /*  fprintf(stderr, " stu: %d", r_edge + 1); */

/* Orig. version -- prior to band_width changes 
    for (maxstu_above = maxstu_vec; *maxstu_above != -1; maxstu_above++)
      *maxstu_above = 0;
*/

  for (; *seq; seq++) {

#if defined(FINDALIGN) || defined(QUICKALIGN)
    if (excess_flag) { /* reduce scores to prevent overflow -- should not
			  interfere with alignment information */
      max_score -= 16000;
      excess += 16000; 
      excess_flag = 0;
      for (maxstu_above = maxstu_vec; *maxstu_above != -1; maxstu_above++) {
	maxstu_above_el = *maxstu_above;
	t = 0;
	if (maxstu_above_el < 0) {
	  maxstu_above_el = -maxstu_above_el;
	  t = maxstu_above_el >> 16;
	  maxstu_above_el &= 65535;
	}
	t = t < 16000 ? 0 : t - 16000;
	maxstu_above_el = maxstu_above_el < 16000 ? 0 : maxstu_above_el - 16000;
	*maxstu_above = t > 0 ? -((t << 16) + maxstu_above_el) : maxstu_above_el;
      }
    }
#endif

    score_vec = q_profile->scores[q_profile->convert[*seq]] + q_left; 
    maxstu_above = maxstu_vec;
    maxstu_above_el = 0; 

#ifdef FINDALIGN
    diri = (char *)our_alloc(alloc_size * sizeof(char));
    /* fprintf(stderr, " alloc: %d", alloc_size); */
    *dir = diri - l_edge - q_left; 
    if (l_edge < 0) diri = *dir + q_left; 
/* N.B. *dir always points to 0 position; diri to l_edge or q_left, whichever
   is greater */
#endif

    if (band_width)  {
      if (l_edge > 0) {
	score_vec += l_edge;
	maxstu_above += l_edge;

	maxstu_above_el = *(maxstu_above - 1);
	if (maxstu_above_el < 0) {
	  maxstu_above_el = -maxstu_above_el;
	  maxstu_above_el &= 65535;
	}
      }
    
      maxstu_vec[r_edge + 1] = -1;
      /*      fprintf(stderr, "%d,", r_edge + 1);  */
      if (maxstu_vec[r_edge] == -1) maxstu_vec[r_edge] = 0;
      if (r_edge < q_len - 1) r_edge++;
      l_edge++;
    }

  zero_u:
    for(;;) {
      maxstu = maxstu_above_el;
      maxstu_above_el = *maxstu_above;

#ifdef COUNTS
      n1++;
#endif

      if (maxstu_above_el >= 0) { /* maxstu_above_el >= 0  -- so t is 0 */
#ifdef COUNTS
	n6++;
#endif
	maxstu += *score_vec++;
	if (maxstu <= 0) {
#ifdef FINDALIGN
	  d = 0; 
#endif
#ifdef COUNTS
	    n10++;
#endif
	  *maxstu_above++ = 0;
	}
	else { 
#ifdef FINDALIGN
	  d = 3; 
#endif
#ifdef COUNTS
	    n7++;
#endif
	  if (maxstu > n_gap_init) {
#ifdef COUNTS
	    n8++;
#endif
	    t = u = maxstu + gap_init;
	    goto transition0;
	  }
	  else {

#ifdef COUNTS
	    n9++;
#endif
	    *maxstu_above++ = maxstu;
	  }
	}
      }
      else {

#ifdef COUNTS
	n2++;
#endif
	
	maxstu_above_el = -maxstu_above_el;
	if (maxstu_above_el == 1)
	  goto next_row; /* -1 is sentinel for end of row */

/* decompose maxstu_above_el into t and maxstu parts */
	maxstu += *score_vec++;
	t = maxstu_above_el >> 16;
	maxstu_above_el &= 65535;

#ifdef FINDALIGN
	d = maxstu_above_el + gap_init != t ? 7 : 3;
	  /* set 4's bit -- to indicate that extending gap */
#endif
	if (maxstu > t) {
	  if (maxstu > n_gap_init) {

#ifdef COUNTS
	    n3++;
#endif

	    t += t_gap_ext;
	    u = maxstu + gap_init;
	    if (t < u) t = u;
	    goto transition0;
	  }
#ifdef COUNTS
	    n4++;
#endif
	}
	else { 
#ifdef FINDALIGN
	  d -= 2; /* remove two's bit */
#endif

#ifdef COUNTS
	    n5++;
#endif
	  maxstu = t; 
/* note no new maximum in score is possible, and no transition to positive
   u is necessary (because would imply adjacent gaps on opposite strands);
   so test for maxstu > n_gap_init is unnecessary */
	} 
	t += t_gap_ext;
	*maxstu_above++ = t > 0 ? -((t << 16) + maxstu) : maxstu;
      }
#ifdef FINDALIGN
      *diri++ = d; 
#endif
    }
  transition0:
#ifdef FINDALIGN
    u_ext = 0; 
#endif

#ifdef COUNTS
	    n11++;
#endif

  transition: /* transition (to case u positive) has occurred */
#ifdef COUNTS
	    n12++;
#endif
    if (max_score < maxstu) { /* N.B. This placement only allows detection of
				 scores exceeding -gap_init */
#ifdef COUNTS
	    n13++;
#endif
      max_score = maxstu;  
#ifdef QUICKALIGN
      s_pos = seq - db_entry_seq;
      q_pos = maxstu_above - maxstu_vec + q_left; 

      if (n_maxes && max_list[n_maxes - 1].s_pos == s_pos) 
	n_maxes--;
      max_list[n_maxes].q_pos = q_pos;
      max_list[n_maxes].s_pos = s_pos;
      max_list[n_maxes].score = max_score + excess;
      n_maxes++;
#endif
#ifdef FINDALIGN
      mdir = dir; 
      mdiri = diri; 
#endif
      if (max_score > 32000) {
#if defined(FINDALIGN) || defined(QUICKALIGN)
	excess_flag = 1;
#else
	goto return_code; /* this sufficient for rapid database searches */
#endif
/*	printf("\nWARNING: Score exceeds 32000; alignment truncated\n"); 
#ifdef FINDALIGN
	*diri++ = d; 
#endif
	goto finish;
*/
      }
    }

    *maxstu_above++ = -((t << 16) + maxstu);

#ifdef FINDALIGN
    *diri++ = d; 
#endif

/* loop for case u positive */
    do {
#ifdef COUNTS
	    n14++;
#endif
      maxstu = maxstu_above_el + *score_vec++;
      maxstu_above_el = *maxstu_above;
      if (maxstu_above_el < 0) {
#ifdef COUNTS
	    n15++;
#endif
	maxstu_above_el = -maxstu_above_el;
	if (maxstu_above_el == 1)
	  goto next_row; /* -1 is sentinel for end of row */
	
	/* decompose maxstu_above_el into t and maxstu parts */
	t = maxstu_above_el >> 16;
	maxstu_above_el &= 65535;
	
#ifdef FINDALIGN
	d = u_ext ? 11 : 3; /* set 8's bit -- to indicate that extending u gap */
	if (maxstu_above_el + gap_init != t) d += 4;
	  /* set 4's bit -- to indicate that extending t gap */
#endif
	if (maxstu > t) {
#ifdef COUNTS
	  n16++;
#endif
	  if (maxstu > u) {
#ifdef COUNTS
	    n17++;
#endif
	    if (maxstu > n_gap_init) {
#ifdef COUNTS
	      n18++;
#endif
	      temp = maxstu + gap_init; 
	      t += t_gap_ext;
	      if (t < temp) t = temp; 
	      u += u_gap_ext;
	      if (u < temp) {
#ifdef COUNTS
		n19++;
#endif
		u = temp;
		goto transition0;
	      }
#ifdef FINDALIGN
	      u_ext = 1;
#endif
#ifdef COUNTS
	      n20++;
#endif
	      goto transition;
	    }
	  }
	  else {
#ifdef COUNTS
	    n21++;
#endif
	    maxstu = u; 
#ifdef FINDALIGN
	    d -= 1; 
#endif
	  }
	}
	else { 
#ifdef COUNTS
	    n22++;
#endif
	  if (u > t) {
#ifdef COUNTS
	    n23++;
#endif
	    maxstu = u; 
#ifdef FINDALIGN
	    d -= 1; 
#endif
	  }
	  else {
#ifdef COUNTS
	    n24++;
#endif
	    maxstu = t; 
#ifdef FINDALIGN
	    d -= 2; 
#endif
	  }
	} 
	t += t_gap_ext;
	*maxstu_above++ = t > 0 ? -((t << 16) + maxstu) : maxstu;
      }
      else { /* maxstu_above_el >= 0  -- so t is 0; and maxstu >= u, which
	      is positive */
#ifdef FINDALIGN
	d = u_ext ? 11 : 3;
#endif
#ifdef COUNTS
	    n25++;
#endif
	if (maxstu > u) {
#ifdef COUNTS
	  n26++;
#endif
	  if (maxstu > n_gap_init) {
#ifdef COUNTS
	    n27++;
#endif
	    t = maxstu + gap_init; 
	    u += u_gap_ext;
	    if (u < t) {
#ifdef COUNTS
	      n28++;
#endif
	      u = t;
	      goto transition0;
	    }
#ifdef FINDALIGN
	    u_ext = 1;
#endif
	    goto transition;
	  }
	  *maxstu_above++ = maxstu;
	}
	else {
#ifdef COUNTS
	    n29++;
#endif
	  *maxstu_above++ = u;
#ifdef FINDALIGN
	  d -= 1; 
#endif
	}
      }
#ifdef FINDALIGN
      u_ext = 1;
      *diri++ = d; 
#endif
      u += u_gap_ext;
    } while (u > 0);
    goto zero_u;

  next_row:
#ifdef FINDALIGN
    dir++;
#endif
    ; 
  }

#if defined(FINDALIGN) || defined(QUICKALIGN)
  max_score += excess;
  *add_orig_score = score = max_score;
#endif

#ifdef QUICKALIGN
  db_entry_seq[saved_pos] = s_replace; 
  if (max_score && max_score >= minscore) {
    *add_success = alignment_from_max_list(n_maxes, max_list, &score, l_edge_orig, r_edge_orig);
    if (!*add_success) {
/*
      printf("\n\nQuick:");
      print_alignment();
*/
      get_stats(&start1, &end1, &start2, &end2, &temp2, &temp2, &temp2);
      max_score = full_smith_waterman(q_profile, db_entry_seq, db_entry_length, l_edge_orig, r_edge_orig,
			  q_left_orig, end1, s_left_orig, end2, 
			  minscore, &temp2);
/*
      printf("\n\nFull:");
      print_alignment();
*/
      if (temp2 != *add_orig_score) fatalError("Score discrepancy");
    }
    else max_score = score;

/*
    if (success) {
      get_stats(&start1, &end1, &start2, &end2,
		&mismatches, &insertions, &deletions);
      if (score1 != score || orig_score1 != orig_score || start1a != start1 || end1a != end1
	  || start2a != start2 || end2a != end2 
	  || mismatchesa != mismatches || insertionsa != insertions || deletionsa != deletions) 
      fprintf(stderr,"\nDiscrepancy %d-%d : %d-%d : %d-%d : %d-%d : %d-%d : %d-%d : %d-%d : %d-%d : %d-%d", 
	      score1, score, orig_score1, orig_score, start1a, start1, end1a, end1, 
	      start2a, start2, end2a, end2 , mismatchesa, mismatches, 
	      insertionsa, insertions, deletionsa, deletions);
    }
*/
  }
  return max_score;
#endif

#ifdef FINDALIGN
  if (max_score && max_score >= minscore) {
    alignment_from_direc(mdir, mdiri, direc, &score);
    max_score = score;
  }

  direc += s_left; 
  l_edge = base_l_edge; 
  for (i = 0; i < s_len; i++) {
    our_free(direc[i] + l_edge + q_left); 
    if (band_width) l_edge++;
  }

  our_free(direc);
#endif

 return_code:
  db_entry_seq[saved_pos] = s_replace;
  return max_score;
}

#if defined(COUNTS)

print_counts()
{
  fprintf(stderr, "\nn1: %.0lf, n2: %.0lf, n3: %.0lf, n4: %.0lf, n5: %.0lf\nn6: %.0lf, n7: %.0lf, n8: %.0lf, n9: %.0lf\nn10: %.0lf, n11: %.0lf, n12: %.0lf, n13: %.0lf, n14: %.0lf, n15: %.0lf\nn16: %.0lf, n17: %.0lf, n18: %.0lf, n19: %.0lf\nn20: %.0lf, n21: %.0lf, n22: %.0lf, n23: %.0lf, n24: %.0lf, n25: %.0lf\nn26: %.0lf, n27: %.0lf, n28: %.0lf, n29: %.0lf\n", 
	  n1, n2, n3, n4, n5, n6, n7, n8, n9, 
	  n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, 
	  n20, n21, n22, n23, n24, n25, n26, n27, n28, n29);
}
#endif
