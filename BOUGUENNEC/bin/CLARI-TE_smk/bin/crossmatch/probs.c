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


/* Probabilites for Smith-Waterman algorithm */

#include "swat.h"

double aafq[20]
	/* Amino acid residue frequencies used by S. Altschul */
	= {.081, .057, .045, .054, .015, .039, .061, .068, .022, .057,
 		  .093, .056, .025, .040, .049, .068, .058, .013, .032, .067 }
	;

main(argc, argv)
     int argc;
     char *argv[];
{
  char rand_query[301];
  double uniform_probs[20];
  char aa_list[26], matrix[50];
  int i, j, gap_init, gap_ext, length, target_score, end_gap;
  Profile *make_profile_from_seq();
  Profile *q_profile;
  double prob, fudge;
  double smith_waterman_prob();
  FILE *fp;
  FILE *fopenWrite();
  int q_len, s_len;

  fudge = .30; /* fudge factor -- should eliminate when understand 
		  the dependencies better */
  srandom(1971);
  s_len = 100;
  q_len = 300;
  strcpy(aa_list, "ARNDCQEGHILKMFPSTWYV"); /* "ACDEFGHIKLMNPQRSTVWY"); */
  strcpy(matrix, "BLOSUM62");
  gap_init = -12;
  gap_ext = -2;
  end_gap = -1;
  printf("Matrix: %s, gap_init: %d, gap_ext: %d\nq_len: %d, s_len: %d, fudge: %.2f\n", 
	 matrix, gap_init, gap_ext, q_len, s_len, fudge);
  fp = fopenWrite("temp.query");
  fprintf(fp, ">Rand.query\n");
  for (i = 0; i < q_len; i++) {
    rand_query[i] = aa_list[random() % 20];
    printf("%c",rand_query[i]);
    fprintf(fp, "%c",rand_query[i]);
  }
  rand_query[i] = 0;
  fprintf(fp,"\n");
  fclose(fp);
  fp = fopenWrite("temp.random");
  for (j = 0; j < 10000; j++) {
    fprintf(fp, ">Rand%d\n", j);
    for (i = 0; i < s_len; i++) fprintf(fp, "%c", aa_list[random() % 20]);
    fprintf(fp,"\n");
  }
  fclose(fp);
  for (i = 0; i < 20; i++) uniform_probs[i] = .05;
/* N.B. WILL NEED CODE TO CORRECTLY MATCH UP PROBS WITH LETTERS */

  set_score_mat(matrix, gap_init, gap_ext, end_gap, 0);
  q_profile = make_profile_from_seq((Profile *)0, rand_query, 0, 0); /* instead put in length */
  printf("\n");
  for (i = 35; i <= 100; i += 5) {
    prob = smith_waterman_prob(q_profile, aa_list, uniform_probs, i, s_len, fudge);
    printf("%d %g\n", i, prob);
  }
}

  
double smith_waterman_prob(q_profile, letters, lett_probs, target_score, length, fudge)
     Profile *q_profile;
     int target_score, length;
     double *lett_probs;
     char *letters;
     double fudge;

{
  int i, j, k, targ_plus_gap, targ_plus_score, q_len;
  char *seq;
  int score, gap_init, gap_ext;
  SCORE *score_vec;
  char *our_alloc();
  double **s_prob, **t_prob, **maxstu_prob, **new_maxstu_prob, **temp_vec;
  double *maxstu_probi, *t_probi, *s_probi, *new_maxstu_probi;
  double *u_prob, *new_u_prob, *u_ptr, *new_u_ptr, *temp_ptr;
  double aa_prob, target_prob, row_target_prob, aa_target_prob, cell_target_prob;
  double x, y, t, u, s;

  q_len = q_profile->length;
/*  for (i = 0; i < q_len; i++)    printf("%c",q_profile->seq[i]); */
  gap_init = q_profile->gap_init;
  gap_ext = q_profile->gap_ext;
  targ_plus_gap = target_score + gap_init;
/*  printf("%d %d %d %d %d\n", q_len, gap_init, gap_ext, length, target_score); */

  u_prob = (double *)our_alloc((target_score - gap_init) * sizeof(double));
  new_u_prob = (double *)our_alloc((target_score - gap_init) * sizeof(double));
  u_prob -= gap_init;
  new_u_prob -= gap_init;

  /* assumes gap_init is negative */
  for (j = gap_init; j <= 0; j++) new_u_prob[j] = u_prob[j] = 1;
  for (j = 1; j < target_score; j++) new_u_prob[j] = u_prob[j] = 0;

  maxstu_prob = (double **)our_alloc(q_len * sizeof(double *));
  t_prob = (double **)our_alloc(q_len * sizeof(double *));
  s_prob = (double **)our_alloc(q_len * sizeof(double *));
  new_maxstu_prob = (double **)our_alloc(q_len * sizeof(double *));

  for (i = 0; i < q_len; i++) {
    maxstu_prob[i] = (double *)our_alloc(target_score * sizeof(double));
    new_maxstu_prob[i] = (double *)our_alloc(target_score * sizeof(double));
    t_prob[i] = (double *)our_alloc(target_score * sizeof(double));
    s_prob[i] = (double *)our_alloc(target_score * sizeof(double));
    new_maxstu_prob[i][0] = maxstu_prob[i][0] = t_prob[i][0] = s_prob[i][0] = 1.0;
    for (j = 1; j < target_score; j++) 
      new_maxstu_prob[i][j] = maxstu_prob[i][j] = t_prob[i][j] = s_prob[i][j] = 0;
  }
/* note: following assumes that gap_init is larger (in abs value) than largest score 
 --- which need not be correct! */  
  for (seq = letters; *seq; seq++) { /* now seq plays the role of list of possible aas */
    aa_prob = lett_probs[q_profile->convert[*seq]];
    score_vec = q_profile->scores[q_profile->convert[*seq]];
    score = score_vec[0];
    for (j = 1; j <= score; j++) {
      new_maxstu_prob[0][j] += aa_prob;
    }
  }
  for (j = 1; j < target_score; j++) {
    maxstu_prob[0][j] = s_prob[0][j] = new_maxstu_prob[0][j];
  }

  target_prob = 0;
  
  for (k = 0; k < length; k++) {
/*    printf("%d\n",k); */
    row_target_prob = 0;
    for (i = 1; i < q_len; i++) 
      for (j = 1; j < target_score; j++) 
	new_maxstu_prob[i][j] = s_prob[i][j] = 0;   /* first column (j = 0) never changes */
    
    for (seq = letters; *seq; seq++) { /* now seq plays the role of list of possible aas */
      aa_target_prob = 0;
      aa_prob = lett_probs[q_profile->convert[*seq]];
      score_vec = q_profile->scores[q_profile->convert[*seq]];
      
/*      for (j = gap_init; j <= 0; j++) u_prob[j] = 1;   gap_init + score */
      for (j = 1; j < targ_plus_gap; j++) u_prob[j] = 0;
/* target probs not set -- assuming target_score > largest pos score in matrix */

      for (i = 1; i < q_len; i++) {
	u_ptr = u_prob + gap_init - gap_ext;
	new_u_ptr = new_u_prob + gap_init;
	score = score_vec[i];
	maxstu_probi = maxstu_prob[i - 1] - score;
	t_probi = t_prob[i];
	s_probi = s_prob[i];
	new_maxstu_probi = new_maxstu_prob[i];
	if (score > 0) {
	  cell_target_prob = maxstu_probi[target_score]; /* require target_score 
							    >= largest score */
	  for (j = 1; j <= score; j++) {
	    new_maxstu_probi[j] += aa_prob;
	    s_probi[j] += aa_prob;
	    new_u_ptr[j] = 1;
	  }
	  if (cell_target_prob) {
	    aa_target_prob += (1 - aa_target_prob) * cell_target_prob;
	    if (cell_target_prob >= 1.0) 
	      printf("\n %.6g cell_target_prob ERROR\n", cell_target_prob);
	    y = 1. / (1 - cell_target_prob);
	    for (; j < targ_plus_gap; j++) {
	      u = fudge * u_prob[j];
	      t = fudge * t_probi[j];
	      s = maxstu_probi[j];
	      s = (s - cell_target_prob) * y; 
	      x = u + (1 - u) * (t + (1 - t) * s);
/*
	      if (s && (x-s)/s > .1 && k > 10) 
		printf("%c %d %d %d %g %g %g %g\n", *seq, k, i, j, (x-s)/s, u, t, s);
*/
	      new_maxstu_probi[j] += aa_prob * x; /* */
	      s_probi[j] += aa_prob * s;
	      new_u_ptr[j] = s + (1 - s) * u_ptr[j];
/*
	      if (new_u_ptr[j] == .05 && j > -gap_init ) 
		printf(" %d %d %g %g\n",i,j,s,u_ptr[j]);
*/
	    }
	    for (; j < target_score; j++) {
	      s = maxstu_probi[j];
	      s = (s - cell_target_prob) * y; 
	      new_maxstu_probi[j] += aa_prob * s; 
	      s_probi[j] += aa_prob * s;
	      new_u_ptr[j] = s + (1 - s) * u_ptr[j];
/*
	      if (new_u_ptr[j] == .05 && j > -gap_init ) 
		printf(" %d %d %g %g\n",i,j,s,u_ptr[j]);
*/

	    }
	  }
	  else {
	    for (; j < targ_plus_gap; j++) {
	      u = fudge * u_prob[j];
	      t = fudge * t_probi[j];
	      s = maxstu_probi[j];
	      x = u + (1 - u) * (t + (1 - t) * s);
/*
	      if (s && (x-s)/s > .1 && k > 10) 
		printf("%c %d %d %d %g %g %g %g\n", *seq, k, i, j, (x-s)/s, u, t, s);
*/
	      new_maxstu_probi[j] += aa_prob * x; /* */
	      s_probi[j] += aa_prob * s;
	      new_u_ptr[j] = s + (1 - s) * u_ptr[j];
/*
	      if (new_u_ptr[j] == .05 && j > -gap_init ) 
		printf(" %d %d %g %g\n",i,j,s,u_ptr[j]);
*/

	    }
	    for (; j < target_score; j++) {
	      s = maxstu_probi[j];
	      new_maxstu_probi[j] += aa_prob * s; 
	      s_probi[j] += aa_prob * s;
	      new_u_ptr[j] = s + (1 - s) * u_ptr[j];
/*
	      if (new_u_ptr[j] == .05 && j > -gap_init ) 
		printf(" %d %d %g %g\n",i,j,s,u_ptr[j]);
*/

	    }
	  }
	}
	else {
	  targ_plus_score = target_score + score;
	  for (j = 1; j < targ_plus_gap; j++) { 
                /* require gap_init <= smallest neg substitution score */
	    u = fudge * u_prob[j];
	    t = fudge * t_probi[j];
	    s = maxstu_probi[j];
	    x = u + (1 - u) * (t + (1 - t) * s);
/*
	      if (s && (x-s)/s > .1 && k > 10) 
		printf("%c %d %d %d %g %g %g %g\n", *seq, k, i, j, (x-s)/s, u, t, s);
*/
	    new_maxstu_probi[j] += aa_prob * x; /* */
	    s_probi[j] += aa_prob * s;
	    new_u_ptr[j] = s + (1 - s) * u_ptr[j];
/*
	      if (new_u_ptr[j] == .05 && j > -gap_init ) 
		printf(" %d %d %g %g\n",i,j,s,u_ptr[j]);
*/

	  }
	  for (; j < targ_plus_score; j++) {
	    s = maxstu_probi[j];
	    new_maxstu_probi[j] += aa_prob * s; 
	    s_probi[j] += aa_prob * s;
	    new_u_ptr[j] = s + (1 - s) * u_ptr[j];
/*
	      if (new_u_ptr[j] == .05 && j > -gap_init ) 
		printf(" %d %d %g %g\n",i,j,s,u_ptr[j]);
*/

	  }
	  for (; j < target_score; j++) {
	    new_u_ptr[j] = u_ptr[j];
/*
	      if (new_u_ptr[j] == .05 && j > -gap_init ) 
		printf(" %d %d %g %g\n",i,j,s,u_ptr[j]);
*/

	  }
	}

/* need to define zero scores for maxstu_prob[-1] or [i][j - score] when j - score > target_score, i.e. gap_init smaller than score */


/* j = 0 unnec. -- since will always be 1 */
	
	/* assumes gap_init >= gap_ext */
	temp_ptr = new_u_prob;
	new_u_prob = u_prob;
	u_prob = temp_ptr;
	/* remaining u_probs are 0 -- have already been set */
      }
      row_target_prob += aa_prob * aa_target_prob; 
    }
/*    printf("%g\n",row_target_prob); */
    temp_vec = new_maxstu_prob;
    new_maxstu_prob = maxstu_prob;
    maxstu_prob = temp_vec;

    for (i = 1; i < q_len; i++) { 
      t_probi = t_prob[i];
      s_probi = s_prob[i] - gap_init;
      for (j = 1; j < targ_plus_gap; j++) {
	x = t_probi[j - gap_ext];
        t_probi[j] = x + (1 - x) * s_probi[j];
      }
    }
    target_prob += (1 - target_prob) * row_target_prob;
/*  printf("%d  %.5g  %.5g\n", k, row_target_prob, target_prob); */
  }
/*
  printf("%10.5g\n", row_target_prob);
  for (j = 0; j < target_score; j++)
    printf("%2d %10.5g %10.5g %10.5g %10.5g\n", 
	   j, u_prob[j], t_prob[q_len - 1][j], s_prob[q_len - 1][j], 
	   maxstu_prob[q_len - 1][j]);
*/
/*  printf("%.5g  ", row_target_prob); */
  our_free(u_prob + gap_init);
  our_free(new_u_prob + gap_init);

  for (i = 0; i < q_len; i++) {
    our_free(maxstu_prob[i]);
    our_free(new_maxstu_prob[i]);
    our_free(t_prob[i]);
    our_free(s_prob[i]);
  }

  our_free(maxstu_prob);
  our_free(new_maxstu_prob);
  our_free(t_prob);
  our_free(s_prob);
  return target_prob;
}
