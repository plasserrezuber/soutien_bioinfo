/*****************************************************************************
#   Copyright (C) 1993-1998 by Phil Green.                          
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

#define MAX_HISTSCORE 200 /* maximum score allowed in score histogram */
#define MAX_HISTLENGTH 20000 /* maximum length in length histogram */

#define LAMBDA_TOL .0001  /* controls accuracy of lambda estimation algorithm */
#define H_TOL .001  /* controls accuracy of h estimation algorithm */
#define CENSOR .01  /* defines amount (at top) of score distribution to be
		       censored for purposes of estimating lambda and K */
#define LENGTH_CUTOFF 30 /* minimum database sequence length allowed, for fitting */
#define Z_CUTOFF 7.0

extern Parameters *parameters;

static int *score_hist, *length_hist;
static double *score_sums, *score2_sums, *length_sums;
static double logK, lambda, rho, mu, rho2, mu2, h;
static int max_histscore, max_histlength,
  max_score, min_score, max_length, min_length, zero_scores, z_rejects;
static int num_db_entries, total_db_length, n_pruned;
static float avg_entry_length, var_cutoff, z_cutoff;

reject_entry(score_entry)
     Score_entry *score_entry;
{
  return (!parameters->nw_flag && !score_entry->score || score_entry->length < LENGTH_CUTOFF);
}

alloc_hist()
{
  char *our_alloc();

  max_histlength = MAX_HISTLENGTH;
  max_histscore = MAX_HISTSCORE;

  score_hist = max_histscore + (int *)our_alloc((2 * max_histscore + 1) * sizeof(int));
  length_sums = max_histscore + (double *)our_alloc((2 * max_histscore + 1) * sizeof(double));
  length_hist = (int *)our_alloc((max_histlength + 1) * sizeof(int));
  score_sums = (double *)our_alloc((max_histlength + 1) * sizeof(double));
  score2_sums = (double *)our_alloc((max_histlength + 1) * sizeof(double));
}
  
initialize_hist()
{
  int i;

  for (i = -max_histscore; i <= max_histscore; i++)
    score_hist[i] = length_sums[i] = 0;
  for (i = 0; i <= max_histlength; i++)
    length_hist[i] = score_sums[i] = score2_sums[i] = 0;
  zero_scores = 0;
  z_rejects = 0;
  z_cutoff = Z_CUTOFF;
}

update_hist(score_entry, z_flag)
     Score_entry *score_entry;
     int z_flag;
{
  int score, length;

  if (reject_entry(score_entry)) {
    zero_scores++;
    return;
  }
  if (z_flag && score_entry->z > z_cutoff) {
    z_rejects++;
    return;
  }
  score = score_entry->score;
  length = score_entry->length;
  if (score > max_histscore) score = max_histscore;
  if (score < -max_histscore) score = -max_histscore;
  length_sums[score] += length;
  score_hist[score] += 1;
  if (length > max_histlength) length = max_histlength;
  length_hist[length] += 1;
  score_sums[length] += score;
  score2_sums[length] += score * score;
}  

prune_hist(score_entry)
     Score_entry *score_entry;
{
  int score, length;

  if (reject_entry(score_entry)) return;

  score = score_entry->score;
  length = score_entry->length;

  if (score > max_histscore) score = max_histscore;
  if (score < -max_histscore) score = -max_histscore;
  length_sums[score] -= length;
  score_hist[score] -= 1;
  if (length > max_histlength) length = max_histlength;
  length_hist[length] -= 1;
  score_sums[length] -= score;
  score2_sums[length] -= score * score;
  n_pruned++;
}  
 
process_hist()
{
  int j, cum;

  for (j = -max_histscore; j <= max_histscore && !score_hist[j]; j++);
  min_score = j; 

  for (j = max_histscore; j >= -max_histscore && !score_hist[j]; j--);
  max_score = j; 

  for (j = 0; j <= max_histlength && !length_hist[j]; j++);
  min_length = j; 
  
  for (j = max_histlength; j >= 0 && !length_hist[j]; j--);
  max_length = j;

  num_db_entries =  total_db_length = cum = 0;
  for (j = max_score; j >= min_score; j--) {
    if (score_hist[j]) {
      num_db_entries += score_hist[j];
      total_db_length += length_sums[j];
      cum += score_hist[j];
      printf("\n%3d   %3d   %3d", j, score_hist[j], cum);
    }
  }


  if (!num_db_entries)
    fatalError("No entries have scores exceeding -gap_init");
/*   printf("\n%d %d %d \n",min_score, max_score, num_db_entries); */
  avg_entry_length = total_db_length / (float)num_db_entries;
  printf("\n\n%d entries with score 0 or length < %d are excluded from the statistical analyses,",
	 zero_scores, (int)LENGTH_CUTOFF);
  printf("\n  leaving %d entries.", num_db_entries);
  printf("\nAverage length: %.1f, (truncated) range %d - %d.",
       avg_entry_length, min_length, max_length);
  printf("\n(Truncated) score range: %d - %d. ", min_score, max_score);
  if (max_score == max_histscore)
    printf("(%d scores are %d or greater)", score_hist[max_score], max_score);
}

fit_log_n(q_length)
     int q_length;
{
  int j;
  int n, last_j, cell_size;
  double mean, x, w, y, y2, t, t2, u, u2, v, p, z, s, s2;
  double mean_x, mean_y, var_x, var_y, covar_xy;
  double mean_y2, covar_xy2;
  double mean_w, var_w, covar_wy, covar_wy2, covar_xw;
  
  /* now fit scores to best linear function of log(n), using ordinary regression */
  
  /*  for (p = -.55; p <= .55; p += .1) {  
      printf("\nPower: %.2f",p);
      */
  mean_x = mean_w = mean_y = mean_y2 = 0.0;
  var_x = var_y = var_w = 0.0;
  covar_xw = covar_wy = covar_xy = covar_xy2 = covar_wy2 = 0.0;
  n = 0;
  for (j = min_length; j <= max_length; j++) 
    if (length_hist[j]) {
      n += length_hist[j];
      x = log((float)j);
      /* w = log((float)(j < q_length ? j : q_length));  sqrt((float)j + p); */
      w = j; /* exp(p * x); j; (j < q_length ? (q_length - j) : 0); */
      mean_x += length_hist[j] * x;
      mean_w += length_hist[j] * w;
      mean_y += score_sums[j];
      var_x += length_hist[j] * x * x;
      var_w += length_hist[j] * w * w;
      var_y += score2_sums[j];
      covar_xy += x * score_sums[j];
      covar_wy += w * score_sums[j];
      covar_xw += length_hist[j] * w * x;
    }
  mean_x /= n;
  mean_w /= n;
  mean_y /= n;
  var_x = var_x / n - mean_x * mean_x;
  var_y = var_y / n - mean_y * mean_y;
  var_w = var_w / n - mean_w * mean_w;
  if (var_x < .0001) var_x = 0.0;
  if (var_y < .0001) var_y = 0.0;
  if (var_w < .0001) var_w = 0.0;
  
  covar_xy = covar_xy / n - mean_x * mean_y;
  covar_xw = covar_xw / n - mean_x * mean_w;
  covar_wy = covar_wy / n - mean_w * mean_y;
  rho = var_x ? covar_xy / var_x : 0.0;
  mu = mean_y - rho * mean_x;
  z = covar_wy - rho * covar_xw;
  printf("\n\nZ-statistic analysis, excluding %d entries with z-scores above %.1f:",
	 z_rejects, z_cutoff);
  printf("\nProportion of score variance explained by log(n): %.2f, by n: %.2f (adds %.2f)",
	 var_x && var_y ? covar_xy * covar_xy / (var_x * var_y) : 0.0,
	 var_w && var_y ? covar_wy * covar_wy / (var_w * var_y) : 0.0,
	 var_x && var_w && var_y ? z * z / ((var_w - covar_xw * covar_xw / var_x) * var_y) : 0.0);
  
  for (j = min_length; j <= max_length; j++) 
    if (length_hist[j]) {
      x = log((float)j);
      /* w = log((float)(j < q_length ? j : q_length)); sqrt((float)j + p); */
      w = j;   /* exp(p * x); (j < q_length ? (q_length - j) : 0); */
      u = rho * x + mu;
      if (u < 0) u = 0; /* Assumes Smith-Waterman */
      y2 = score2_sums[j] - 2 * score_sums[j] * u + length_hist[j] * u * u;
      mean_y2 += y2;
      covar_xy2 += x * y2;
      covar_wy2 += w * y2;
    }
  
  mean_y2 /= n;
  covar_xy2 = covar_xy2 / n - mean_x * mean_y2;
  covar_wy2 = covar_wy2 / n - mean_w * mean_y2;
  rho2 = var_x ? covar_xy2 / var_x : 0.0;
  if (rho2 < 0) rho2 = 0;
  mu2 = mean_y2 - rho2 * mean_x;
  
  if (rho2 < 0) 
    z = (rho2 * log(5000.0) + mu2 > 0) ? 5000 : exp(-1.0 - mu2 / rho2);
  else z =  rho2 ? exp(1.0 - mu2 / rho2) : LENGTH_CUTOFF;
  if (z < LENGTH_CUTOFF) z = LENGTH_CUTOFF;
  printf("\nMinimum allowed predicted variance at n = %.0f",z);
  var_cutoff = rho2 * log(z) + mu2;
  z = covar_wy2 - rho2 * covar_xw;
  printf("\nVariance of residuals explained by log(n): %.2f, by n: %.2f (adds %.2f)",
	 var_x ? covar_xy2 * covar_xy2 / var_x : 0.0, var_w ? covar_wy2 * covar_wy2 / var_w : 0.0,
	 var_x && var_w ? z * z / (var_w - covar_xw * covar_xw / var_x) : 0.0);
  /* } */
  printf("\n\nBest fits to log(n):");
  /*
    printf("\n  mean = %.2f * log(n) + %.2f,   var = %.2f * log(n) + %.2f\n",
    rho, mu, rho2, mu2);
    */
  printf("\n  mean = %.2f * log(n) + %.2f,   var = %.2f * log(n) + %.2f\n",
	 rho, mu, rho2, mu2);
  
  cell_size = n / 10;
  printf("\n\n                          S             S ln(300)/ln(n)         z");
  printf("\n\nLengths (n)   No.     Mean    Std        Mean    Std        Mean    Std\n");
  for (j = last_j = min_length; j <= max_length; ) {
    for (s = s2 = t = t2 = u = u2 = n = 0; n <= cell_size && j <= max_length; j++) 
      if (length_hist[j]) {
	n += length_hist[j];
	s += score_sums[j];
	s2 += score2_sums[j];
	
	x = log(300.0) / log((float)j);
	t += score_sums[j] * x;
	t2 += score2_sums[j] * x * x;
	
	y = log((float)j);
	w = y * y;
	x = rho * y + mu;
	if (x < 0) x = 0; /* Assumes Smith-Waterman */
	v = rho2 * y + mu2; /* rho2 * y + mu2; */
	if (v < var_cutoff) v = var_cutoff;
	u += (score_sums[j] - length_hist[j] * x) / sqrt(v);
	u2 += (score2_sums[j] - 2 * x * score_sums[j] + length_hist[j] * x * x) /
	  v;
      }
    mean = s / (float)n;
    printf("\n%3d - %4d:   %4d  %6.1f  %6.2f", last_j, j - 1, n, mean,
	   sqrt(s2 / (float)n - mean * mean));
    mean = t / n;
    printf("     %6.1f  %6.2f", mean, sqrt(t2 / n - mean * mean));
    mean = u / n;
    printf("     %6.1f  %6.2f", mean, sqrt(u2 / n - mean * mean));
    last_j = j;
  }  
}  

est_lambda_K(q_length) 
     int q_length;
{
  int n_scores, i, j;
  double av_score, av_score2, av_length, t_length, new_lambda, sum1, sum2,
         cum, prev_cum, expected, logL, cutoff;
  double high_lambda, low_lambda, d_lambda, increment, oldLogL, oldoldLogL;
  double temp1, temp2, temp;
  int cum_hist;
  double find_score_Evalue();
  float censor; /* fraction of scores at top of dist'n to be censored (not
			used in estimating parameters) */


  if (!min_score) min_score = 1; /* for consistency */
  censor = CENSOR;
  i = 0;
  if (!censor) {
    cutoff = max_score + 100;
    j = max_score + 1;
  }
  else 
    for (j = max_score; j >= min_score; j--)
      if (score_hist[j]) {
	i += score_hist[j];
	if (i >= censor * num_db_entries) {
	  cutoff = j - .5;
	  break;
	}
      }

  n_scores = av_score = av_score2 = t_length = 0;
  for (j--; j >= min_score; j--)
    if (score_hist[j]) {
      t_length += length_sums[j];
      n_scores += score_hist[j];
      av_score += j * score_hist[j];
      av_score2 += j * j * score_hist[j];

    }
  av_length = t_length / n_scores;
  av_score /= n_scores;
  av_score2 = sqrt(av_score2 / n_scores - av_score * av_score);

  printf("\n\n%d scores below %.1f (average: %.1f, std: %.1f) were used in estimating \ndistribution parameters. \nTotal length of entries used: %.0f, avg. length: %.1f\n",
	 n_scores, cutoff, av_score, av_score2, t_length, av_length);

  high_lambda = 2.0;
  low_lambda = 0.0;
  new_lambda = 1.0;
  do {
    lambda = new_lambda;
    sum1 = sum2 = 0;
    for (i = min_score; i <= cutoff; i++) {
      if (score_hist[i]) {
	temp1 = length_sums[i] * exp(-lambda * (i - min_score));
	sum1 += temp1;
	sum2 += i * temp1;
      }
    }
    temp2 = censor ? t_length * exp(-lambda * (cutoff - min_score)) : 0.0;
    sum1 -= temp2;
    sum2 -= cutoff * temp2;
    logK = log((double)n_scores / q_length) - log(sum1) + lambda * min_score;
    logL = logK +log(lambda) - lambda * av_score; 
    d_lambda = 1.0/lambda - av_score + sum2 / sum1;
    if (d_lambda > 0) low_lambda = lambda;
    else high_lambda = lambda;
    new_lambda = (low_lambda + high_lambda) * .5;
/*    new_lambda = 1.0 / (av_score - sum2 / sum1);  
    if (new_lambda < 0.0)
      fatalError("lambda estimation failed");
*/
  } while (fabs(new_lambda - lambda) > LAMBDA_TOL);
  printf("\nLambda: %.4f,  logK: %.4f,  logL: %.4f, d_lambda: %10.4f",
	 lambda, logK, logL, d_lambda);
  printf("\n\n");
  prev_cum = num_db_entries - find_score_Evalue(max_score + 1, q_length);
  cum_hist = 0;
  printf("\nObserved and expected score distribution:\n");
  printf("\n                        Cumulative");
  printf("\nScore  Obs.   Exp.     Obs.      Exp.");
  for (i = max_score; i >= min_score; i--) {
    expected = find_score_Evalue(i, q_length);
    cum = num_db_entries - expected;
    cum_hist += score_hist[i];
    if (score_hist[i]) {
      printf(expected < 1 ? "\n%3d   %4d %6.1f    %5d %10.3g" : "\n%3d   %4d %6.1f    %5d %10.1f",
	   i, score_hist[i], (prev_cum - cum),
			cum_hist, expected);
    }
    prev_cum = cum;
  }
  printf("\n");
  printf("\n\nPredicted fits to seq. length (n) :\n  mean = %.2f * log(n) + %.2f, var = %.2f\n",
	 1.0 / lambda, (logK + log((float)q_length) + .58) / lambda,
	 1.28 * 1.28 / (lambda * lambda));
}

new_est_lambda_K(q_length, score_entries, last_score_entry) /* using non-extreme-value dist'n */
     int q_length;
     Score_entry *score_entries, *last_score_entry;
{
  Score_entry *score_entry;
  int n_scores, i, j;
  double av_score, av_score2, av_length, t_length, new_lambda, sum1, sum2, sum3,
         cum, prev_cum, expected, logL, cutoff;
  double hx, x, n, m, explx, h;
  double high_lambda, low_lambda, d_lambda, d_h, increment, oldLogL, oldoldLogL;
  double high_h, low_h, new_h;
  double temp1, temp2, temp3, temp;
  int cum_hist;
  double new_find_score_Evalue();
  float censor; /* fraction of scores at top of dist'n to be censored (not
			used in estimating parameters) */


  if (!min_score) min_score = 1; /* for consistency */
  censor = CENSOR;
  i = 0;
  if (!censor) {
    cutoff = max_score + 100;
    j = max_score + 1;
  }
  else 
    for (j = max_score; j >= min_score; j--)
      if (score_hist[j]) {
	i += score_hist[j];
	if (i >= censor * num_db_entries) {
	  cutoff = j - .5;
	  break;
	}
      }

  n_scores = av_score = av_score2 = t_length = 0;
  for (j--; j >= min_score; j--)
    if (score_hist[j]) {
      t_length += length_sums[j];
      n_scores += score_hist[j];
      av_score += j * score_hist[j];
      av_score2 += j * j * score_hist[j];

    }
  av_length = t_length / n_scores;
  av_score /= n_scores;
  av_score2 = sqrt(av_score2 / n_scores - av_score * av_score);

  printf("\n\n%d scores below %.1f (average: %.1f, std: %.1f) were used in estimating \ndistribution parameters. \nTotal length of entries used: %.0f, avg. length: %.1f\n",
	 n_scores, cutoff, av_score, av_score2, t_length, av_length);

  m = q_length;
  high_h = 2.9;
  low_h = 0.0;
  new_h = 1.0;

  do {
    h = new_h;
    high_lambda = 2.0;
    low_lambda = 0.0;
    new_lambda = 1.0;
    do {
      lambda = new_lambda;
      sum1 = sum2 = sum3 = 0;
      for (i = min_score; i <= cutoff; i++) {
	if (score_hist[i]) {
	  explx = exp(-lambda * (i - min_score));
	  temp1 = (m - h * i) * (length_sums[i] - h * score_hist[i] * i) * explx;
	  sum1 += temp1;
	  sum2 += i * temp1;
	}
      }
      if (censor) {
	explx = exp(-lambda * (cutoff - min_score));
	temp1 = (m - h * cutoff) * (t_length - h * n_scores * cutoff) * explx;
	sum1 -= temp1;
	sum2 -= cutoff * temp1;
      }
      logK = log((double)n_scores) - log(sum1) + lambda * min_score;

      d_lambda =  sum2 / sum1 - av_score;

      sum2 = 0.0;
      for (score_entry = score_entries; score_entry < last_score_entry; score_entry++) {
	if (reject_entry(score_entry)) continue;
	x = score_entry->score;
	if (censor && x > cutoff) continue;
	n = score_entry->length;
	if (h * x >= m || h * x >= n) continue; /* NEED TO CORRECT TERMS ABOVE ALSO */
	hx = h * x;
	temp1 = (m - hx) * (n - hx);
	temp3 = m + n - (hx + hx);
	temp2 = h * temp3 + lambda * temp1;
	sum2 += temp1 / temp2;
      }
      d_lambda += sum2 / n_scores;

      /* now sum over all entries -- but multiply by 1 / n_scores */

      if (d_lambda > 0) low_lambda = lambda;
      else high_lambda = lambda;
      new_lambda = (low_lambda + high_lambda) * .5;

      printf("\nLambda: %.4f, K: %.4f, h: %.2f, d_lambda: %.4f",
	     lambda, exp(logK), h, d_lambda);
      /*    new_lambda = 1.0 / (av_score - sum2 / sum1);  
	    if (new_lambda < 0.0)
	    fatalError("lambda estimation failed");
	    */
    } while (fabs(new_lambda - lambda) > LAMBDA_TOL);

    sum1 = sum3 = 0;
    for (i = min_score; i <= cutoff; i++) {
      if (score_hist[i]) {
	explx = exp(-lambda * (i - min_score));
	temp1 = (m - h * i) * (length_sums[i] - h * score_hist[i] * i) * explx;
	sum1 += temp1;
	sum3 += i * (length_sums[i] + (m - 2 * h * i)  * score_hist[i]) * explx; 
      }
    }
    if (censor) {
      explx = exp(-lambda * (cutoff - min_score));
      temp1 = (m - h * cutoff) * (t_length - h * n_scores * cutoff) * explx;
      sum1 -= temp1;
      temp1 = cutoff * (t_length + (m - 2 * h * cutoff) * n_scores) * explx;
      sum3 -= temp1;
    }
    logL = logK - 1.0 - lambda * av_score; 
    d_h = sum3 / sum1;

    sum1 = sum3 = 0.0;
    for (score_entry = score_entries; score_entry < last_score_entry; score_entry++) {
      if (reject_entry(score_entry)) continue;
      x = score_entry->score;
      if (censor && x > cutoff) continue;
      n = score_entry->length;
      if (h * x >= m || h * x >= n) continue; /* NEED TO CORRECT TERMS ABOVE ALSO */
      hx = h * x;
      temp1 = (m - hx) * (n - hx);
      temp3 = m + n - (hx + hx);
      temp2 = h * temp3 + lambda * temp1;
      sum3 += ((1.0 - lambda * x) * temp3 - hx - hx) / temp2;
      sum1 += log(temp2);
    }
    d_h += sum3 / n_scores;
    logL += sum1 / n_scores;

    printf("\n\nLambda: %.4f, h: %.2f, logL: %.4f, d_h: %.4f",
	   lambda, h, logL, d_h);
    printf("\n\n");

    prev_cum = num_db_entries - new_find_score_Evalue(max_score + 1, q_length, h);
    cum_hist = 0;
    printf("\nObserved and expected score distribution:\n");
    printf("\n                        Cumulative");
    printf("\nScore  Obs.   Exp.     Obs.      Exp.");
    for (i = max_score; i >= min_score; i--) {
      expected = new_find_score_Evalue(i, q_length, h); /* CHANGE h!!! */
      cum = num_db_entries - expected;
      cum_hist += score_hist[i];
      if (score_hist[i]) {
	printf(expected < 1 ? "\n%3d   %4d %6.1f    %5d %10.3g" : "\n%3d   %4d %6.1f    %5d %10.1f",
	       i, score_hist[i], (prev_cum - cum),
	       cum_hist, expected);
      }
      prev_cum = cum;
    }
    if (d_h > 0) low_h = h;
    else high_h = h;
    new_h = (low_h + high_h) * .5;
  } while (fabs(new_h - h) > H_TOL);
  printf("\n");
  printf("\n\nPredicted fits to seq. length (n) :\n  mean = %.2f * log(n) + %.2f, var = %.2f\n",
	 1.0 / lambda, (logK + log((float)q_length) + .577216) / lambda,
	 1.28255 * 1.28255 / (lambda * lambda));
}

/* assumed global variables: K, lambda, num_db_entries, total_db_length,
   min_length, max_length, length_hist, h */
double new_find_Evalue(score, q_length, s_length)
     int score, q_length, s_length;
{
  double E, temp, cum, s;
  int j;

  s = score + 1.0;

  if (q_length <= h * s || s_length <= h * s) return 0.0;

  temp = (q_length - h * s) * (s_length - h * s) * exp(-lambda * s + logK);
  E = num_db_entries * (temp < .01 ? temp : 1.0 - exp(-temp));

  return E;
}

/* assumed global variables: K, lambda, num_db_entries, total_db_length,
   min_length, max_length, length_hist */
/* finds expected no. of entries with a given score or higher */
double new_find_score_Evalue(score, q_length, h)
     int score, q_length;
     double h;
{
  double E, temp, cum, s;
  int j;

  s = score - .5;

  if (q_length <= h * s || total_db_length < num_db_entries * h * s) return 0.0;

  temp = (q_length - h * s) * exp(-lambda * s + logK);
  if (temp * max_length < .05)  E = (total_db_length - num_db_entries * h * s) * temp;
  else {
    cum = 0.0;
    for (j = min_length; j <= max_length; j++)
      if (length_hist[j] && (h * s < j)) 
	cum += length_hist[j] * exp(-temp * (j - h * s));
    E = num_db_entries - cum; 
  }

  return E;
}

/* assumed global variables: K, lambda, num_db_entries, total_db_length,
   min_length, max_length, length_hist */
double find_score_Evalue(score, q_length)
     int score, q_length;
{
  double E, temp, cum;
  int j;

  temp = q_length * exp(-lambda * (score - .5) + logK);
  if (temp * max_length < .05)  E = total_db_length * temp;
  else {
    cum = 0.0;
    for (j = min_length; j <= max_length; j++)
      if (length_hist[j]) 
	cum += length_hist[j] * exp(-temp * j);
    E = num_db_entries - cum; 
  }

  return E;
}

static double *log_table;

find_z(entry)
     Score_entry *entry;
{
  float log_len, var, mean;
  int length;

  length = entry->length > LENGTH_CUTOFF ? entry->length : LENGTH_CUTOFF;
  if (length <= MAX_HISTLENGTH) {
    if (!log_table) make_log_table();
    
    if (log_table[length] < 0) log_table[length] = log((float)length);
    log_len = log_table[length];
  }
  else log_len = log((float)length);
  var = rho2 * log_len + mu2;
  if (var < var_cutoff) var = var_cutoff;
  mean = rho * log_len + mu;
  if (mean < 0) mean = 0; /* Assumes Smith-Waterman */
  entry->z = (entry->score - mean) / sqrt(var);
}

make_log_table()
{
  char *our_alloc();
  int i;

  log_table = (double *)our_alloc((MAX_HISTLENGTH + 1) * sizeof(double));
  for (i = 0; i <= MAX_HISTLENGTH; i++) log_table[i] = -1;
}


double z_to_E(z) /* computes E value for a given z value, assuming (Gumbel) extreme
		   value distribution with mean 0 and variance 1 -- should be valid
		   when behavior is local */
     float z;
{
  double e;

  e = exp(-1.28255 * z - .577216); /* pi / sqrt(6); Euler-Mascheroni */
  return num_db_entries * (e > .01 ? 1.0 - exp(-e) : e);
}


