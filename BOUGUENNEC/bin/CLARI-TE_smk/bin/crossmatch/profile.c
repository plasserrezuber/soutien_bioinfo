/*****************************************************************************
#   Copyright (C) 1993-2009 by Phil Green.                          
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

#define BIGNEG -1000000

extern Parameters *parameters;

static char alphabet[256];
static int alphsize, uc_alphsize;
static char *row_labels;
static int n_rows, n_cols, n_places, min_entry, max_entry;
static int **score_mat, **comp_c_score_mat;
static int *score_area, *c_score_area, *comp_c_score_area, *full_score_area;
static int gap_init, gap_ext, ins_gap_ext, del_gap_ext, end_gap;
static int profile_flag;
static double *freqs[2];
static int freq_flag;
static int merged_flag;

double *log_freqs[2], lambda;
int **c_score_mat, **full_score_mat, n_c_rows, n_c_cols;
int row_convert[256], col_convert[256];
double *match_score_adjust, *mismatch_score_adjust; /* score adjustments for each quality value, to
								    correct rounding simplifications */
/* Changes made 9/9/06 to compress sizes of profiles, by eliminating duplicate scoring rows and columns from score matrix;
   also allowing lower case. Need to check whether this affects lambda, log_freqs, query_histogram
*/

/* Following routine reads in score matrix, and sets gap penalties; if file_name is 0 or the empty
 string, it is assumed that the matrix is for nucleotides (upper or lower
 case) and mismatch penalties are set to mismatch, except for N and X (which
 have penalties equal to 0 and -1 respectively). mismatch is ignored when
 file_name is not empty */

int get_DNA_flag()
{
  /*  n < 10 ? 1 : (n < 20 ? 2 : 0); /* NOT A GOOD CRITERION !!! */
  return merged_flag ? 2 : uc_alphsize < 20 ? 1 : 0;
}

int get_alphsize()
{
  return alphsize;
}

int get_uc_alphsize()
{
  return uc_alphsize;
}

int get_gap_init()
{
  return gap_init;
}

int get_gap_ext()
{
  return gap_ext;
}

int get_ins_gap_ext()
{
  return ins_gap_ext;
}

int get_del_gap_ext()
{
  return del_gap_ext;
}

set_gap_penalties(gap_init0, gap_ext0, ins_gap_ext0, del_gap_ext0, end_gap0)
     int gap_init0, gap_ext0, ins_gap_ext0, del_gap_ext0, end_gap0;
{
  gap_init = gap_init0;
  gap_ext = gap_ext0;
  ins_gap_ext = ins_gap_ext0;
  del_gap_ext = del_gap_ext0;
  end_gap = end_gap0;
}

int set_score_mat()
{
  FILE *fp;
  FILE *fopenRead();
  char buffer[501];
  char *our_alloc();
  int i, j, k, c, d, n_lines, alph_flag, file_flag, gap_flag, i_row, i_col, i_c_row, i_c_col, s, q, orig_qual;
  int full_n_cols;
  char letter;
  float frequency;
  double total, x;
  double pre_freqs[256];
  int uc_residues[256];
  int row_to_c_row[256], c_row_to_row[256], col_to_c_col[256], c_col_to_col[256];
  int *convert;

  file_flag = parameters->matrix && parameters->matrix[0];
  if (score_mat) {
    our_free(score_area);
    our_free(score_mat);
    our_free(c_score_area);
    our_free(c_score_mat);
    our_free(comp_c_score_area);
    our_free(comp_c_score_mat);
  }
  gap_flag = 0;
  if (row_labels) our_free(row_labels);
  if (file_flag) {
    freqs[0] = freqs[1] = 0;
    freq_flag = 0;
    merged_flag = 0;
    for (i = 0; i < 256; i++) pre_freqs[i] = 0;
    fp = fopenRead(parameters->matrix);
    n_rows = alph_flag = n_lines = 0;
    while (fgets(buffer, 500, fp)) {
      if (!alph_flag) n_lines++; /* count lines prior to matrix itself */
      for (i = 0; buffer[i] && isspace(buffer[i]); i++);
      if (!buffer[i] || buffer[i] == '#') continue; /* ignore empty or commented lines */
      if (!strncmp(buffer + i, "GAP", 3)) {
	gap_flag = 1;
	if (2 != sscanf(buffer + i + 3, "%d %d", &gap_init, &gap_ext))
	  fatalError("Score matrix format: GAP line incorrect");
      }
      else if (!strncmp(buffer + i, "MERGED BASE DNA", 15)) {
	merged_flag = 1;
      }
      else if (!strncmp(buffer + i, "FREQS", 5)) {
	freq_flag = 1;
	k = i + 5; 
	do {
	  while (buffer[k] && isspace(buffer[k])) k++;
	  if (!buffer[k]) break;
	  if (1 != sscanf(buffer + k, "%c ", &letter))
	    fatalError("Incorrect score matrix format -- FREQS line");
	  for (k++; buffer[k] && isspace(buffer[k]); k++);
	  if (1 != sscanf(buffer + k, "%f", &frequency))
	    fatalError("Incorrect score matrix format -- FREQS line");
	  while (buffer[k] && !isspace(buffer[k])) k++;
	  pre_freqs[letter] = frequency; /* store backgd frequencies assoc to given letters */
	} while (buffer[k]);
      }
      else if (!alph_flag) {
	for (alphsize = 0; buffer[i]; i++) 
	  if (!isspace(buffer[i])) {
	    if ((isalpha(buffer[i]) || buffer[i] == '*') && alphsize < 256)
	      alphabet[alphsize++] = buffer[i]; /* toupper(buffer[i]); */
	    else fatalError("Incorrect score matrix format-- reading column labels");
	  }
	alph_flag = 1;
      }
      else n_rows++;
    }
    rewind(fp);
    for (i = 0; i < n_lines; i++) fgets(buffer, 500, fp);
  }
  else {
    alphsize = n_rows = 6;
    strcpy(alphabet,"ACGTNX");
    gap_flag = 1;
    if (parameters->qual_scores) {
      /* use quality-based scoring: i is a quality scored nuc (upper 2 bits: nuc, lower 6 bits: quality)
	 s(i,j) = 6 if i,j bases match, = 6 - (q + 5) if they mismatch and q is base quality of i
            here 6 = -10log10(1/4), 5 = -10log10(1/3)
         effectively s is to give scores ~ -10log10(x/1+x) where x/(1+x) is posterior prob of contam
           (contam has prior of .5, given location of read has prior of 1/2G, probs of read are 
           (1/4)^L for contam, and prods of (1-e_i) and (e_i/3) where e_i are base error prob.
           Ignore 1+x since usually is ~1.
	 reserved chars: those with i & 63 == 0: 0 for N, large neg
	 gap_init ?? should not make 2 indels cheaper than subst
	 gap_ext : 
       
	 N.B. Can also have compact_qual without using quality-based scoring; in this case need appropriate translation
         matrix
      */
      match_score_adjust = 1 + (double *)our_alloc(62 * sizeof(double));
      mismatch_score_adjust = 1 + (double *)our_alloc(62 * sizeof(double));
      gap_init = -30; /* with score_offset of 60, used -120 and -60 -- why? */
      gap_ext = -20;
    }
    else {
      gap_init = parameters->penalty - 2;
      gap_ext = parameters->penalty - 1;
    }
  }

  for (i = 0; i < 256; i++) uc_residues[i] = 0;
  for (i = 0; i < alphsize; i++) uc_residues[toupper(alphabet[i])] += 1;
  for (i = uc_alphsize = 0; i < 256; i++) if (uc_residues[i]) uc_alphsize++;

  for (i = 1; i < alphsize; i++)
    for (j = 0; j < i; j++)
      if (alphabet[i] == alphabet[j])
	fatalError("Incorrect score matrix format -- duplicate column labels");


/* N.B. NEED CODE FOR PROTEINS; AND FOR PROFILES !!!  */

  n_cols = parameters->compact_qual ? 256 : alphsize;

  score_mat = (int **)our_alloc(n_rows * sizeof(int *));
  score_area = (int *)our_alloc(n_rows * n_cols * sizeof(int));
  row_labels = (char *)our_alloc(n_rows * sizeof(char));
  for (i = 0; i < n_rows; i++) {
    row_labels[i] = alphabet[i]; /* default */
    score_mat[i] = score_area + i * n_cols; /* (int *)our_alloc(n_cols * sizeof(int)); */
    if (file_flag) {
      do {
	fgets(buffer, 500, fp);
	for (k = 0; buffer[k] && isspace(buffer[k]); k++);
      } while (!buffer[k] || buffer[k] == '#');
      if (isalpha(buffer[k]) || buffer[k] == '*') {
	row_labels[i] = buffer[k];
	k++;
      }
      else if (i >= n_cols) 
	fatalError("Incorrect score matrix format -- too many rows");
      for (j = 0; j < n_cols; j++) {
	while (buffer[k] && isspace(buffer[k])) k++;
	if (!buffer[k] || 1 != sscanf(buffer + k, "%d", &score_mat[i][j]))
	  fatalError("Incorrect score matrix format -- while reading rows");
	while (buffer[k] && !isspace(buffer[k])) k++;
      }
    }
    else {
/*      printf("\n"); */
      c = alphabet[i];
      for (j = 0; j < n_cols; j++) {
	if (parameters->compact_qual) { /* holds true if parameters->qual_scores is set */
	  q = j & 63;
	  if (q) {
	    d = "ACGT"[j >> 6];
	    if (parameters->qual_scores) { /* otherwise q is irrelevant */
	      q--; /* converts back to orig qual scale */
	      if (q < 2) q = 2; /* qualities <= 1 are not sensible (implied error rate > 3/4) */
	      orig_qual = q;
	      q = q / parameters->qual_resolution + .5;  
	      q *= parameters->qual_resolution; /* now rounded to nearest mult of qual_resolution */
	      if (!q) q = parameters->qual_resolution; /* round up lowest qualities -- so there is some penalty */
	      mismatch_score_adjust[orig_qual] = q - orig_qual;
	      q += 5; /* allows for factor of 1/3 in prob of a particular erroneous base occurring */
	    }
	  }
	  else d = j == 64 ? 'N' : 'X';
	}
	else {
	  d = alphabet[j];
	}
	score_mat[i][j] = toupper(d) == 'N' || toupper(c) == 'N' ? 0
	  : d == 'X' || c == 'X' ? parameters->X_penalty 
	  : toupper(c) == toupper(d) ?  parameters->match_reward 
	  : parameters->qual_scores ? parameters->match_reward - q : parameters->penalty;

/*	else if (islower(c) || islower(d)) score_mat[i][j] = gap_ext / 2; */
/*	if (score_mat[i][j] >= 60) fprintf(stderr, "%d:%d:%d ",i, j, score_mat[i][j]);  */
      }
    }
  }
  if (parameters->qual_scores) { /* otherwise q is irrelevant */
    for (q = 2; q < 61; q++) { /* qualities of 1 or less have error probs > .75 -- i.e. worse than random guessing! */
      x = pow(10.0, -q / 10.0);
      match_score_adjust[q] = 10.0 * log10(1 - x);
    }
    for (q = -1; q < 2; q++)
      match_score_adjust[q] = match_score_adjust[2];
    mismatch_score_adjust[-1] = mismatch_score_adjust[0];
    /* 
    for (q = -1; q < 61; q++) {
      fprintf(stderr, "\n%d %.2f %.2f", q, match_score_adjust[q], mismatch_score_adjust[q]); 
    }
    /* */
  }
  
  if (file_flag) fclose(fp); 
  profile_flag = 0;
  for (i = 1; i < n_rows && !profile_flag; i++)
    for (j = 0; j < i && !profile_flag; j++)
      if (row_labels[i] == row_labels[j]) profile_flag = 1;

  if (!profile_flag) {

    row_to_c_row[0] = 0; /* for each i_row <_row, gives compressed index i_c_row */
    c_row_to_row[0] = 0; /* for each i_c_row up to n_c_row, gives 1st corresponding i_row in score matrix */ 
    n_c_rows = 1; /* # of distinct rows */

    for (i_row = 1; i_row < n_rows; i_row++) {
      for (i_c_row = 0; i_c_row < n_c_rows; i_c_row++) {
	for (j = 0; j < n_cols && score_mat[i_row][j] == score_mat[c_row_to_row[i_c_row]][j]; j++);
	if (j == n_cols) break;
      }
      row_to_c_row[i_row] = i_c_row;
      if (i_c_row == n_c_rows) {
	c_row_to_row[n_c_rows++] = i_row;
      }
    }
    if (n_rows > n_c_rows) printf("\n%d distinct rows out of %d in score matrix", n_c_rows, n_rows);

    for (i = 0; i < 256; i++) row_convert[i] = -1; 
    for (i_row = 0; i_row < n_rows; i_row++) row_convert[row_labels[i_row]] = row_to_c_row[i_row];
    for (i = 0; i < 256; i++) {
      if (row_convert[i] == -1) {
	if (row_convert[toupper(i)] != -1) row_convert[i] = row_convert[toupper(i)];
	else if (row_convert[tolower(i)] != -1) row_convert[i] = row_convert[tolower(i)];
	else row_convert[i] = n_c_rows - 1; 
      }
    }
    for (i = 0; i < 256; i++) 
      if (row_convert[i] < 0 || row_convert[i] >= n_c_rows)
	fatalError("row conversion");

  }   
  else {
    n_c_rows = n_rows;
    for (i_row = 0; i_row < n_rows; i_row++) {
      row_to_c_row[i_row] = c_row_to_row[i_row] = i_row;
    }
  }

  n_places = n_cols + 1 + 10;

  col_to_c_col[0] = 0; /* for each i_col <_col, gives compressed index i_c_col */
  c_col_to_col[0] = 0; /* for each i_c_col up to n_c_col, gives i_col in score matrix */ 
  n_c_cols = 1; /* # of distinct cols */

  for (i_col = 1; i_col < n_cols; i_col++) {
    for (i_c_col = 0; i_c_col < n_c_cols; i_c_col++) {
      for (j = 0; j < n_rows && score_mat[j][i_col] == score_mat[j][c_col_to_col[i_c_col]]; j++);
      if (j == n_rows) break;
    }
    col_to_c_col[i_col] = i_c_col;
    if (i_c_col == n_c_cols) {
      c_col_to_col[n_c_cols++] = i_col;
    }
  }
  if (n_cols > n_c_cols) printf("\n%d distinct cols out of %d in score matrix", n_c_cols, n_cols);

  for (i = 0; i < 256; i++) col_convert[i] = -1; /* n_c_cols - 1; */
  for (i_col = 0; i_col < n_cols; i_col++) col_convert[parameters->compact_qual ? i_col : alphabet[i_col]] = col_to_c_col[i_col];
  for (i = 0; i < 256; i++) {
    if (col_convert[i] == -1) {
      if (col_convert[toupper(i)] != -1) col_convert[i] = col_convert[toupper(i)];
      else if (col_convert[tolower(i)] != -1) col_convert[i] = col_convert[tolower(i)];
      else col_convert[i] = n_c_cols - 1; 
    }
  }

  for (i = 0; i < 256; i++) 
    if (col_convert[i] < 0 || col_convert[i] >= n_c_cols)
      fatalError("column conversion");

  c_score_mat = (int **)our_alloc(n_c_rows * sizeof(int *));
  c_score_area = (int *)our_alloc(n_c_rows * n_c_cols * sizeof(int));
  comp_c_score_mat = (int **)our_alloc(n_c_cols * sizeof(int *));
  comp_c_score_area = (int *)our_alloc(n_c_rows * n_c_cols * sizeof(int));
  max_entry = min_entry = 0;
  for (i_c_row = 0; i_c_row < n_c_rows; i_c_row++) 
    c_score_mat[i_c_row] = c_score_area + i_c_row * n_c_cols; 

  for (i_c_col = 0; i_c_col < n_c_cols; i_c_col++) 
    comp_c_score_mat[i_c_col] = comp_c_score_area + i_c_col * n_c_rows;

  for (i_c_row = 0; i_c_row < n_c_rows; i_c_row++) {
    for (i_c_col = 0; i_c_col < n_c_cols; i_c_col++) {
      c_score_mat[i_c_row][i_c_col] = s = score_mat[c_row_to_row[i_c_row]][c_col_to_col[i_c_col]];
      comp_c_score_mat[i_c_col][i_c_row] = s;
      if (max_entry < s) max_entry = s;
      if (min_entry > s) min_entry = s;
    }
  } 

  full_n_cols = parameters->compact_qual ? 256 : 128;
  full_score_mat = (int **)our_alloc(128 * sizeof(int *));
  full_score_area = (int *)our_alloc(55 * full_n_cols * sizeof(int));
  for (i_c_row = i = 0; i_c_row < 128; i_c_row++) {
    if (i_c_row == 0 || isalpha(i_c_row)) {
      full_score_mat[i_c_row] = full_score_area + i;
      for (i_c_col = 0; i_c_col < full_n_cols; i_c_col++)
	full_score_mat[i_c_row][i_c_col] = !i_c_row || !i_c_col ? BIGNEG : c_score_mat[row_convert[i_c_row]][col_convert[i_c_col]];
      i += full_n_cols;
    }
    else full_score_mat[i_c_row] = 0;
  }

  /*
  for (i_c_row = 0; i_c_row < 128; i_c_row++) {
    if (!full_score_mat[i_c_row]) continue;
    fprintf(stderr, "\n%c ", (char)i_c_row);
    for (i_c_col = 0; i_c_col < full_n_cols; i_c_col++) {
      fprintf(stderr, "%d ", full_score_mat[i_c_row][i_c_col]);
    }
  }
  exit(1);
  */
  /*
    for (i_c_col = 0; i_c_col < n_c_cols; i_c_col++) {
      printf("\nDEBUG: row %d ", i_c_col);
      for (i_c_row = 0; i_c_row < n_c_rows; i_c_row++) {
	printf(" %d", comp_c_score_mat[i_c_col][i_c_row]);
      }
    }
  */
  /*
    for (i = 0; i < 256; i++) col_convert[i] = n_cols - 1;
    for (i = 0; i < n_cols; i++) col_convert[alphabet[i]] = i;
  */

  if (!freqs[0]) {
    for (j = 0; j < 2; j++) {
      freqs[j] = (double *)our_alloc(256 * sizeof(double));
      log_freqs[j] = (double *)our_alloc(256 * sizeof(double));
      for (i = 0; i < 256; i++) freqs[j][i] = 0;
    }
    if (freq_flag) {
      for (j = 0; j < 2; j++) {
	convert = j ? col_convert : row_convert;
	for (i = 0; i < 256; i++) freqs[j][convert[i]] += pre_freqs[i];   /* IS SUMMATION CORRECT?? */
      }
    }
    else if (uc_alphsize < 20) {
      freqs[0][row_convert['A']] += .25;
      freqs[0][row_convert['C']] += .25;
      freqs[0][row_convert['G']] += .25;
      freqs[0][row_convert['T']] += .25;
      if (parameters->compact_qual) {
	for (i = 0; i < 256; i++) {
	  if (i & 63)
	    freqs[1][col_convert[i]] += 1.0;
	}
      }
      else {
	freqs[1][col_convert['A']] += .25;
	freqs[1][col_convert['C']] += .25;
	freqs[1][col_convert['G']] += .25;
	freqs[1][col_convert['T']] += .25;
      }
    }
    else {
      for (j = 0; j < 2; j++) {
	convert = j ? col_convert : row_convert;
	freqs[j][convert['A']] += .0756;
	freqs[j][convert['C']] += .0171;
	freqs[j][convert['D']] += .0531;
	freqs[j][convert['E']] += .0633;
	freqs[j][convert['F']] += .0408;
	freqs[j][convert['G']] += .0688;
	freqs[j][convert['H']] += .0224;
	freqs[j][convert['I']] += .0574;
	freqs[j][convert['K']] += .0596;
	freqs[j][convert['L']] += .0934;
	freqs[j][convert['M']] += .0214;
	freqs[j][convert['N']] += .0456;
	freqs[j][convert['P']] += .0493;
	freqs[j][convert['Q']] += .0404;
	freqs[j][convert['R']] += .0517;
	freqs[j][convert['S']] += .0721;
	freqs[j][convert['T']] += .0578;
	freqs[j][convert['V']] += .0654;
	freqs[j][convert['W']] += .0127;
	freqs[j][convert['Y']] += .0322;
      }
    }
    for (j = 0; j < 2; j++) {
      for (i = total = 0; i < 256; i++) total += freqs[j][i];
      if (total) for (i = 0; i < 256; i++) freqs[j][i] /= total;
      for (i = 0; i < 256; i++) log_freqs[j][i] = freqs[j][i] ? log(freqs[j][i]) : 0.0;
    }
  }

  get_lambda();

/*  printf("\n%s\n",alphabet); */
  return gap_flag;
}

print_score_mat()
{
  int i, j;

  printf("\n ");
  for (i = 0; i < alphsize; i++) printf("   %c",  alphabet[i]);
  for (i = 0; i < n_rows; i++) {
    printf("\n%c",  row_labels[i]);
    for (j = 0; j < n_cols; j++) {
      printf(" %3d", score_mat[i][j]);
    }
  } 
}

/* do row frequencies instead? following is from when there was only one freq set */
print_background_freqs(complexity_freq_type)
     int complexity_freq_type;
{
  int i, j, n_indices, *convert;

  n_indices = complexity_freq_type ? n_c_cols : n_c_rows;
  convert = complexity_freq_type ? col_convert : row_convert;
  printf(" Assumed background frequencies:\n ");
  for (i = 0; i < n_indices; i++) {
    if (i < n_indices - 1)
      for (j = 0; j < 256; j++)
	if (convert[j] == i) printf("%c", j);
    printf(": %.3f  ", freqs[complexity_freq_type][i]);
  }
}


/* N.B. NO LONGER COPY seq -- SO REQUIRE IT STAYS CONSTANT DURING USAGE OF profile!!! */

Profile *make_profile_from_seq(current_profile, seq, length, nw_flag)
     Profile *current_profile;
     unsigned char *seq;
     int length, nw_flag;
{
  char *our_alloc();
  Profile *profile;
  int i, j, k, alloc_flag;
  int offset, max_score_cutoff; /* , packed_length; */

  if (!length)
    for (length = 0; seq[length]; length++); 

  if (profile_flag) {
    if (length != n_rows)
	fatalError("Query sequence does not match profile matrix row labels");
    for (i = 0; i < n_rows; i++)
      if (seq[i] != row_labels[i])
	fatalError("Query sequence does not match profile matrix row labels");
  }

  alloc_flag = 1;
  if (current_profile) {
    if (current_profile->alloc_length >= length) {
      profile = current_profile;
      alloc_flag = 0;
    }
    else {
      free_profile(current_profile);
    }
  }
  if (alloc_flag) {
    profile = (Profile *)our_alloc(sizeof(Profile));
    profile->alloc_length = length;
    /*
    profile->seq = (char *)our_alloc((length + 1) * sizeof(char));
    */
    /*   profile->maxstu_vec = (int *)our_alloc((length + 1) * sizeof(int)); */
/*    profile->t_vec = (int *)our_alloc((length + 1) * sizeof(int)); CURRENTLY NOT NEEDED 
*/
    profile->convert = (int *)our_alloc(256 * sizeof(int));
    profile->n_c_cols = n_c_cols;
    profile->scores = (SCORE **)our_alloc(n_c_cols * sizeof(SCORE *));
    profile->score_area = (SCORE *)our_alloc(n_c_cols * (length + 1) * sizeof(SCORE));
    for (i = 0; i < n_c_cols; i++) {
      profile->scores[i] = profile->score_area + i * (length + 1); /* (int *)our_alloc((1 + length) * sizeof(int)); */
    }

    /* fprintf(stderr, "\nProfile allocation: %d (cols) x %d (int size) x %d (len) bytes", n_c_cols, sizeof(SCORE), length); */
    if (parameters->query_histograms) {
      profile->hist = (int **)our_alloc(length * sizeof(int *));
      for (i = 0; i < length; i++) {
	profile->hist[i] = (int *)our_alloc(n_places * sizeof(int));
      }
    }
  }
  profile->length = length;
  profile->seq = seq;

 /* strcpy(profile->seq, seq); */
  /*  for (j = 0; j < length; j++) profile->maxstu_vec[j] = 0; PRECEDING MAY STILL BE NEC FOR NW.C /* profile->t_vec[j] = 0; */
  /*  profile->maxstu_vec[length] = -1; /* sentinel for end of vector in smith-waterman
				  algorithm -- ; can never occur
				  with actual scores */
  /* profile->t_vec[length] = 0; /* sentinel (in comb w/ maxstu_vec[j] in n-w algorthm,
			    since never should exceed maxstu_vec with actual
			    scores */
  for (i = 0; i < 256; i++) profile->convert[i] = col_convert[i];

/*
  for (i = 1, j = 0; i < 256; i++) N.B. strchr(alphabet,0) is non-null !
    profile->convert[i] = strchr(alphabet, i) ?
      strchr(alphabet, i) - alphabet : alphsize - 1;
*/
/* previously used toupper(i) instead of i 
*/
/*  profile->score_pos = (int **)our_alloc(alphsize * sizeof(int *)); 
*/
  profile->max_entry = max_entry;
  profile->min_entry = min_entry; /* ACTUALLY THIS COULD BE SLIGHTLY DIFFERENT FROM MAX & MIN OVER PROFILE ITSELF!! */
  for (i = 0; i < n_c_cols; i++) {
/*    profile->score_pos[i] = (int *)our_alloc((1 + length) * sizeof(int)); 
*/
    if (!profile_flag && !nw_flag)
      for (j = 0; j < length; j++) 
	profile->scores[i][j] = comp_c_score_mat[i][row_convert[seq[j]]];
    else {
      for (j = k = 0; j < length; j++) {
	profile->scores[i][j] = 
	  profile_flag ? comp_c_score_mat[i][j] : comp_c_score_mat[i][row_convert[seq[j]]];

	if (nw_flag) profile->scores[i][j] -= gap_ext + gap_ext;

	/* last term is adjustment for needleman-wunsch algorithm, to allow assumption that gap_ext = 0 
	 */
	/* N.B. ADJUSTMENTS NEEDED FOR INS_GAP_EXT, DEL_GAP_EXT 
	 */
	/*      if (profile->scores[i][j] > 0) profile->score_pos[i][k++] = j; 
	 */
      }
    }
    profile->scores[i][length] = 0;
/*    profile->score_pos[i][k] = length; 
*/
  }
  if (parameters->query_histograms) {
    for (i = 0; i < length; i++) 
      for (j = 0; j < n_places; j++) profile->hist[i][j] = 0;
  }

  profile->lambda = lambda;
  profile->log_col_freqs = log_freqs[1];
  profile->gap_init = gap_init;
  profile->gap_ext = gap_ext;
  profile->ins_gap_ext = ins_gap_ext;
  profile->del_gap_ext = del_gap_ext;
  /* temporary settings */
  profile->qbegin_gap_init = end_gap;
  profile->qbegin_gap_ext = end_gap;
  profile->qend_gap_init = end_gap;
  profile->qend_gap_ext = end_gap;

  profile->sbegin_gap_init = end_gap;
  profile->sbegin_gap_ext = end_gap;
  profile->send_gap_init = end_gap;
  profile->send_gap_ext = end_gap;
/*  profile->score_mat = score_mat; 
*/
 
/*  profile->packed_length = packed_length = (length + 8) / 8; /* +8, rather than +7, to allow
								for shift of last byte to next
								iteration (in fast_smith_wat) 
*/
  offset = 1;
  if (profile->min_entry < 0 && offset < -profile->min_entry) offset = -profile->min_entry;
  if (profile->gap_init < 0 && offset < -profile->gap_init) offset = -profile->gap_init;
  if (profile->gap_ext < 0 && offset < -8 * profile->gap_ext) offset = -8 * profile->gap_ext;
  profile->poly_offset = profile->poly_one = 0;
  profile->max_score_cutoff = 0;
  profile->poly_gap_init = 0;
  profile->poly_gap_ext = 0;
  max_score_cutoff = 127 - profile->max_entry;
/* N.B. COULD BE 255 - profile->max_entry -- IF WEREN'T USING 1 BIT FOR A POSITIVE T
   INDICATOR 
*/
  for (i = 0; i < 8; i++) {
    profile->poly_offset <<= 8;
    profile->poly_offset += offset;
    profile->poly_one <<= 8;
    profile->poly_one += 1;
    profile->max_score_cutoff <<= 8;
    profile->max_score_cutoff += max_score_cutoff;
    profile->poly_gap_init <<= 8;
    profile->poly_gap_init += -gap_init;
    profile->poly_gap_ext <<= 8;
    profile->poly_gap_ext += -gap_ext;
  }
/* N.B. ASSUMING gap_init, gap_ext are NON-POSITIVE 
*/
/*  profile->packed_maxstu_vec = (unsigned long int *)our_alloc((2 * packed_length + 2) * sizeof(unsigned long int));
  profile->packed_scores = (unsigned long int **)our_alloc(n_cols * sizeof(unsigned long int *));
  for (i = 0; i < n_cols; i++) {
    profile->packed_scores[i] = (unsigned long int *)our_alloc((1 + packed_length) * sizeof(unsigned long int));

    for (j = 0; j < packed_length * 8; j++) {
      profile->packed_scores[i][j / 8] <<= 8;
      profile->packed_scores[i][j / 8] += offset + (j < length ? profile->scores[i][j] : 0);
    }
    profile->packed_scores[i][packed_length] = 0;
  }
*/
/*
  fprintf(stderr, "\n%d %d %d %d %d\n", gap_init, gap_ext, offset, packed_length, max_score_cutoff);
  fprintf(stderr, "\n%ld %ld %ld %d %ld\n", profile->poly_gap_init & 255, profile->poly_gap_ext & 255, profile->poly_offset & 255, profile->packed_length & 255, profile->max_score_cutoff & 255);
*/
  return profile;
}

Profile *maximal_profile_segments(profile, seq, length, minscore, id, screen)
     Profile *profile;
     unsigned char *seq;
     char *id;
     int minscore, screen, length;
{
  int i;
  Profile *make_profile_from_seq();

  profile = make_profile_from_seq(profile, seq, length, 0); 

  for (i = 0; i < n_c_cols; i++)
    find_segs(profile, i, 0, profile->length - 1, minscore, seq, id, screen);

  return profile;
  /* free_profile(profile); */
}

find_segs(profile, i_alph, j_start, j_end, minscore, seq, id, screen)
     Profile *profile;
     int i_alph, j_start, j_end, minscore, screen;
     unsigned char *seq;
     char *id;
{
  int j;
  int max, cum, min, j_min, j_max, j_min_best;

  max = min = cum = 0;
  j_min = j_start - 1;
  for (j = j_start; j <= j_end; j++) {
    cum += profile->scores[i_alph][j];
    if (cum < min) {
      min = cum;
      j_min = j;
    }
    else if (cum - min > max) {
      max = cum - min;
      j_min_best = j_min + 1;
      j_max = j;
    }
  }
  if (max >= minscore) {
    printf("\n%s   ", id);
    for (j = 0; j < 256; j++) 
      if (col_convert[j] == i_alph) printf("%c", j);
    printf("   Score: %4d   Residues: %5d - %5d", max, j_min_best + 1, j_max + 1);
    if (screen) {
      for (j = j_min_best; j <= j_max; j++) seq[j] = 'X';
    }
    find_segs(profile, i_alph, j_start, j_min_best, minscore, seq, id, screen);
    find_segs(profile, i_alph, j_max + 1, j_end, minscore, seq, id, screen);
  }
}

/* THIS MAY BE WRONG */
show_query_hist(profile)
  Profile *profile;
{
  int i, j, r_total, n_match, match_pos;

  printf("\n\nQuery discrepancy histogram (MAY BE WRONG due to changes):\n              ");
  for (i = 0; i < alphsize; i++) printf(" %c ", alphabet[i]);
  printf("    -  I1 I2 ...");
  for (i = 0; i < profile->length; i++) {
    for (j = r_total = 0; j < n_places; j++) 
      r_total += profile->hist[i][j];
    match_pos = profile->convert[profile->seq[i]]; /* CONVERSION SHOULD NOT APPLY TO QUERY SEQ!! */
    n_match = profile->hist[i][match_pos];
    profile->hist[i][match_pos] = 0;
    printf("\n%c %4d %4d  ", profile->seq[i], r_total, n_match);
    for (j = 0; j < n_cols; j++) printf(" %2d", profile->hist[i][j]);
    printf("   ");
    for (; j < n_places; j++) printf(" %2d", profile->hist[i][j]);
  }
}  


free_profile(profile)
  Profile *profile;
{
  int i;

  if (parameters->query_histograms) {
    for (i = 0; i < profile->length; i++) our_free(profile->hist[i]);
    our_free(profile->hist);
  }
  /*
  for (i = 0; i < n_c_cols; i++) {
    our_free(profile->scores[i]);
    /*    our_free(profile->packed_scores[i]); 
  }
*/
  our_free(profile->score_area); 
  our_free(profile->scores); 
  our_free(profile->convert); 
  /*  our_free(profile->t_vec); 
  our_free(profile->maxstu_vec); 
  */
  /* our_free(profile->seq);  */
  /*
  our_free(profile->packed_maxstu_vec);
  our_free(profile->packed_scores);
  */
  our_free(profile); 
}

/* generate random sequence -- not needed at present 
unsigned char *get_random_seq(length)
     int length;
{
  char *our_alloc();
  unsigned char *seq;
  int i;

  seq = (unsigned char *)our_alloc((length + 1) * sizeof(unsigned char));
  for (i = 0; i < length; i++)
    seq[i] = alphabet[random() % alphsize];
  seq[length] = 0;
  return seq;
}
*/
      
get_lambda()
{
  double get_sum();
  double sum, lambda_lower, lambda_upper;

  lambda_lower = 0.0;
  lambda = 0.5;
  for (;;) {
    sum = get_sum();
    if (sum >= 1.0) break;
    lambda_lower = lambda;
    lambda *= 2.0; 
  } 
  lambda_upper = lambda;
  while (lambda_upper - lambda_lower > .000001) {
    lambda = (lambda_lower + lambda_upper) / 2.0;
    sum = get_sum();
    if (sum >= 1.0) lambda_upper = lambda;
    else lambda_lower = lambda;
  } 
  /* 
  fprintf(stderr,"\n lambda: %f  1/lambda: %f  sum: %f", lambda, 1/ lambda, sum);
  /* */
}

/* MAY BE WRONG -- SINCE ASSUMES SAME FREQS APPLY TO BOTH ROWS AND COLUMNS */
double get_sum()
{
  int i, j;
  double total;
  double check;

  check = 0;
  for (i = total = 0; i < n_c_rows; i++)
    for (j = 0; j < n_c_cols; j++)
      if (freqs[0][i] && freqs[1][j]) {
	total += freqs[0][i] * freqs[1][j] * exp(lambda * c_score_mat[i][j]);
	check += freqs[0][i] * freqs[1][j];
      }
  if (check > 1.001 || check < .999) fatalError("frequency sums"); /* RESTORE THIS */
  return total;
}

/* penalties for spliced alignments: idea is LLR type adjustment. rel likelihood of partic splicing model =
   P(intron start at x) P(intron length L) P(seq | splice model) / P(seq | bkgd model)
   P(intron start at x) = 1 / avg exon size
   P(intron length L) = (1 / u) exp(-(L - min) / u) where u = mean length, min = min size [exponential model]
                      = (1 / (L - min + .5)) (1 / sqrt(2 pi var)) exp(-(log(L - min + .5) - mu)^2 / 2 var) [lognormal model with
		                                                                    mean mu and variance var in logspace]
    note that lognormal fits better than single gamma dist'n, but not as well as 'spliced' gamma (distinct gammas going
                       in opp directions from a given cutpoint). A lognormal mixture might be best. Param choices below
                       are from C. elegans intron analysis (see ubc_solexa directory on solduc, align.perl runs)
   P(seq) ratio includes aligned bases, but also factor of 256/4 = 64 (since 2 poss dinuc pair combins x 2 strands)
   LLR = LLR(seq excluding ss) - log(ex_size * u) + log(64) - (L-min)/u; [exponential]
   score = LLR/lambda = score(seq excluding ss) - (1/lambda)log(ex_size * u / 64) - (1/lambda)((L-min)/u);

   in lognormal case: 
        score = score(seq excluding ss) - (1/lambda)log(ex_size / 64) - (1/lambda)(log(L - min + .5) + .5 * log(2 pi var) + (log(L - min + .5) - mu)^2 / 2 var);
*/

double ipenalty, imultiplier, imultiplier2, lognormal_mean;
int min_intron_size;

set_intron_penalties()
{
  double avg_exon_size, mean_intron_size, ss_seq, mult_fac, lognormal_var, fac1, fac2;
  int i, ipen, last_ipen;

  avg_exon_size = 150;
  ss_seq = 64;
  min_intron_size = parameters->min_intron_length;
  mean_intron_size = 1000; /* arbitrary -- for exponential use only */

  lognormal_mean = 4.562; /* from C. elegans -- but perhaps OK for other organisms if min size is adjusted */
  lognormal_var = 2.821;

  /* for exponential
  fac1 = fac2 = mean_intron_size;
  /* */
  /* */
  fac1 = sqrt(lognormal_var * 2 * 3.14159265);
  fac2 = 1.0;
  /*  */

  ipenalty = log(avg_exon_size * fac1 / ss_seq) / lambda + .5; /*.5 is for rounding */
  imultiplier = 1.0 / (lambda * fac2);
  imultiplier2 = 1.0 / (2 * lognormal_var);
  /* 
  for (i = 30; i < 20000; i++) {
    ipen = get_intron_penalty(i);
    if (ipen != last_ipen) {
      printf("\n%d %d", i, ipen);
      last_ipen = ipen;
    }
  }
  exit(1);
  /* */
}

get_intron_penalty(intron_size)
     int intron_size;
{
  double x;

  if (intron_size < min_intron_size) fatalError("intron size");

  x = log(intron_size - min_intron_size + .5);
  /* 
  return (int)(ipenalty + imultiplier * (intron_size - min_intron_size + 1));
  /* */
  /* */
  return (int)(ipenalty + imultiplier * (x + (x - lognormal_mean) * (x - lognormal_mean) * imultiplier2));
  /* */
}
