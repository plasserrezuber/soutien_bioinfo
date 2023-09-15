/*****************************************************************************
#   Copyright (C) 1994-2003 by Phil Green.                          
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

#define HIST_WORD_SIZE 25

extern FILE *fp_log; /* log file */
extern Parameters *parameters; 

static Database *db;
static char *seq_area; /* array containing sequence data for all entries (see db.c) */
static char *seq_area_plus; /* seq_area + index_word_size */
static int minmatch, maxmatch, max_group_size; /* minimum length of matching words */
static int *s_array; /* array of positions of word starts */
static int s_length; /* length of s_array */
static char d_char1, d_char2; /* degenerate characters -- words containing
				 these are ignored */
static int num_entries, start_entry;
static int num_uncomp_entries;
static int total_length;
static int vector_matches;
static int residues[256], m_residues[256];
static int alphabet_size, num_index_words, index_word_size, mod_value;
static int *index_words;
static int *position_table;
static char *word_block;
static unsigned int mask;
static int word_int2;
static int (*pair_routine)();
static int complexity_flag;
static int num_pass_words, num_le_fail_words, num_repeat_fail_words;
static int repeat_screen;

/*
static int alloc_size;
*/

set_word_db(db0, complements)
     Database *db0;
     int complements;
{
  char *our_alloc();
  int i;
  int make_new_cand_pair(), cluster_pairs();

  db = db0;
  seq_area = db->seq_area;
  num_uncomp_entries = db->num_entries;
  num_entries = complements ? 2 * num_uncomp_entries: num_uncomp_entries;
  total_length = complements ? 2 * db->t_length : db->t_length;
  position_table = (int *)our_alloc((num_entries + 1) * sizeof(int));
  start_entry = db->first_entry;
  position_table[0] = 1;
  for (i = 0; i < num_uncomp_entries; i++) 
    position_table[i + 1] = position_table[i] + get_seq_length(start_entry + i) + 1;
  for (; i < num_entries; i++) 
    position_table[i + 1] = position_table[i] + 
      get_seq_length(start_entry + num_entries - 1 - i) + 1;
  /* check */
  for (i = 0; i < num_entries; i++) 
    if (seq_area[position_table[i] - 1]) fatalError("position_table");

  index_word_size = parameters->indexwordsize;
  minmatch = parameters->minmatch;
  maxmatch = parameters->maxmatch;
  max_group_size = parameters->max_group_size;

  set_residue_arrays();

  if (alphabet_size == 4) {
    mask = 0;
    for (i = 0; i < index_word_size; i++)
      mask = (mask << 2) + 3;
  }

  seq_area_plus = seq_area + index_word_size;
  for (i = 0, mod_value = 1; i < index_word_size - 1; i++) 
    mod_value *= alphabet_size;

  num_index_words = mod_value * alphabet_size;

  for (i = 0; i < 256; i++) 
    m_residues[i] = residues[i] > -1 ? residues[i] * mod_value : -1;

  index_words = (int *)our_alloc((num_index_words + 1) * sizeof(int));
  for (i = 0; i <= num_index_words; i++) index_words[i] = 0;

  if (!strcmp(parameters->calling_program, "cluster"))
    pair_routine = cluster_pairs;
  else pair_routine = make_new_cand_pair;

  repeat_screen = parameters->repeat_screen;

  make_nlogn();
  complexity_flag = !parameters->word_raw;
}

set_residue_arrays()
{
  int i;

  for (i = 0; i < 256; i++) residues[i] = -1;

  if (parameters->DNA_flag) {
    alphabet_size = 4;
    residues['A'] = residues['a'] = 0;
    residues['C'] = residues['c'] = 1;
    residues['G'] = residues['g'] = 2;
    residues['T'] = residues['t'] = 3;
  }
  else {
    alphabet_size = 20;
    residues['A'] = residues['a'] = 0;
    residues['C'] = residues['c'] = 1;
    residues['D'] = residues['d'] = 2;
    residues['E'] = residues['e'] = 3;
    residues['F'] = residues['f'] = 4;
    residues['G'] = residues['g'] = 5;
    residues['H'] = residues['h'] = 6;
    residues['I'] = residues['i'] = 7;
    residues['K'] = residues['k'] = 8;
    residues['L'] = residues['l'] = 9;
    residues['M'] = residues['m'] = 10;
    residues['N'] = residues['n'] = 11;
    residues['P'] = residues['p'] = 12;
    residues['Q'] = residues['q'] = 13;
    residues['R'] = residues['r'] = 14;
    residues['S'] = residues['s'] = 15;
    residues['T'] = residues['t'] = 16;
    residues['V'] = residues['v'] = 17;
    residues['W'] = residues['W'] = 18;
    residues['Y'] = residues['y'] = 19;
  }
}

/* find and sort all words of length minmatch, not containing N, X, 
   or a lower case letter */

new_sort_words()
{
  char *our_alloc();
  int i, j, last_j, w, i_lower, i_upper, word_int;
  char *word_end, *seqj, *seqm;
  int t_length;
  int d, min, max, i_min, i_max, gsize;
  int histogram[201], hist20[201];
  int group_histogram[3][201];
  int *next_lower_length, *next_lower;
  int next_ll[128]; 
  int max_index_block, k, last_k, left_bdry_length, w_length;
  char *seq1, *seq2;

  notify("Finding words ...");

  s_length = t_length = total_length;
  last_j = t_length + num_entries + 1 - minmatch;

/* count no. of words -- to allow minimum allocation */

  seqm = seq_area;
 reinit1:
  seqj = seqm + 1;
  if (seqj >= seq_area + last_j) goto out1; /* always reach this point -- when c is final 0 */
  word_int = 0;
  word_end = seqj + index_word_size;

  for (seqm = seqj; seqm < seqj + minmatch; seqm++) {
    w = residues[*seqm];
    if (w < 0) goto reinit1;
    if (seqm < word_end) word_int = word_int * alphabet_size + w;
  }

  index_words[word_int] += 1;

  if (alphabet_size == 4) {
    for (; ; seqm++, word_end++) {
      if (residues[*seqm] < 0) goto reinit1;
      word_int = ((word_int << 2) & mask) + residues[*word_end];
/*
      word_int = (word_int % mod_value) * alphabet_size + residues[*word_end];
*/
      index_words[word_int] += 1;
    }
  }
  else {
    for (; ; seqm++, word_end++) {
      if (residues[*seqm] < 0) goto reinit1;
      word_int = (word_int - m_residues[*(word_end - index_word_size)]) * alphabet_size + residues[*word_end];
/*
    word_int = (word_int % mod_value) * alphabet_size + residues[*word_end];
*/
      index_words[word_int] += 1;
    }
  }

 out1:

/* temp */
  max = 0;
  min = 1000000;
  max_index_block = 0;
  for (i = 0; i <= 200; i++) histogram[i] = 0;
  for (i = 0; i < num_index_words; i++) {
    d = index_words[i];
    if (max_index_block < d) max_index_block = d;

    if (min > d) {min = d; i_min = i; }
    if (max < d) {max = d; i_max = i; }
    d = (d + 9) / 10;
    histogram[d > 100 ? 100 : d] += 1;
  }
  for (i = 0; i <= 200; i++) {
    if (histogram[i]) fprintf(stderr, "\n %d  %d", i * 10, histogram[i]);
  }
  fprintf(stderr, "\nmin no. words: %d at %d", min, i_min);
  fprintf(stderr, "\nmax no. words: %d at %d", max, i_max);
  next_lower_length = (int *)our_alloc(max_index_block * sizeof(int));

  for (i = 1; i < num_index_words; i++)
    index_words[i] += index_words[i - 1];

  s_length = index_words[num_index_words] = index_words[num_index_words - 1];
  s_array = (int *)our_alloc((s_length + 1) * sizeof(int));

  fprintf(stderr, "\navg no. words: %.1f\n", index_words[num_index_words] / (float)num_index_words);

/*
  alloc_size = s_length + 1;
*/
  s_array[0] = 0; /* create sentinel -- word starting with 0, so less than
		     every other word */
  s_array++;
  word_block = (char *)our_alloc(s_length * sizeof(char));
  
  seqm = seq_area;
 reinit2:
  seqj = seqm + 1;
  if (seqj >= seq_area + last_j) goto out2; /* always reach this point -- when c is final 0 */
  word_int = 0;
  word_end = seqj + index_word_size;
  for (seqm = seqj; seqm < seqj + minmatch; seqm++) {
    w = residues[*seqm];
    if (w < 0) goto reinit2;
    if (seqm < word_end) word_int = word_int * alphabet_size + w;
  }
  index_words[word_int] -= 1;
  s_array[index_words[word_int]] = j = seqj - seq_area;

  if (alphabet_size == 4) {
    for (; ; seqm++, word_end++) {
      if (residues[*seqm] < 0) goto reinit2;
      word_int = ((word_int << 2) & mask) + residues[*word_end];
/*
    word_int = (word_int % mod_value) * alphabet_size + residues[*word_end];
*/
      index_words[word_int] -= 1;
      s_array[index_words[word_int]] = ++j;
    }
  }
  else {
    for (; ; seqm++, word_end++) {
      if (residues[*seqm] < 0) goto reinit2;
      word_int = (word_int - m_residues[*(word_end - index_word_size)]) * alphabet_size + residues[*word_end];
/*
    word_int = (word_int % mod_value) * alphabet_size + residues[*word_end];
*/
      index_words[word_int] -= 1;
      s_array[index_words[word_int]] = ++j;
    }
  }

 out2:

  notify(" Done\n");

  fprintf(fp_log, "No. words: %d; after pruning: %d\n", 
	  t_length, s_length);
  fprintf(stderr, "%d words; after pruning: %d\n", t_length, s_length);

  notify("Sorting words ...");

  for (j = 0; j < 3; j++)
    for (i = 0; i <= 200; i++) 
      group_histogram[j][i] = 0;

  for (i = 0; i <= 200; i++) 
      histogram[i] = hist20[i] = 0;

  for (i = 0; i < num_index_words; i++) {
    i_lower = index_words[i];
    i_upper = index_words[i + 1] - 1;
/* correction 1/14/99: the following statement (i.e. use of an offset pointer) apparently causes
   problems, e.g. segmentation faults on HP machines, apparently because for
   large values of i_lower, next_lower may point to something outside the dataspace;
   this should not be a problem because it is only used with a positive displacement
   that results in pointing to a valid element, and it is not a problem on most machines.
   anyway have redefined next_lower. New version slightly less computationally efficient 
    next_lower = next_lower_length - i_lower;
*/
    next_lower = next_lower_length;

    if (i_upper < i_lower) continue;
    if (i_upper > i_lower)
      quicksort(i_lower, i_upper, index_word_size);
    /* Find length of match between each word and its following one */

    gsize = 1;
    for (j = i_lower; j < i_upper; j++) {
      seq1 = seq_area + s_array[j];
      seq2 = seq_area + s_array[j + 1];
      for (k = index_word_size; seq1[k] == seq2[k] && seq1[k] && k < maxmatch; k++);
      if (k >= HIST_WORD_SIZE) gsize++;
      else {
	hist20[gsize > 200 ? 200 : gsize] += 1;
	gsize = 1;
      }
      word_block[j] = k;
    }
    hist20[gsize > 200 ? 200 : gsize] += 1;
    word_block[i_upper] = minmatch - 1;
    /* Now work back, finding distance to next occurrence of a smaller matching word size */
    for (k = 0; k <= maxmatch; k++)
      next_ll[k] = i_upper;
    next_lower[i_upper - i_lower] = i_upper;
    for (j = i_upper - 1; j >= i_lower; j--) {
      next_lower[j - i_lower] = next_ll[word_block[j] - 1];
      /* Could probably improve efficiency of following loop -- which is potential
	 time hog */
      for (k = word_block[j]; k <= maxmatch; k++)
	next_ll[k] = j;
    }
    left_bdry_length = minmatch - 1; /* index_word_size - 1; */
    for (j = i_lower; j <= i_upper; j++) {
      
      for (last_k = j, k = next_lower[j - i_lower]; 
	   word_block[last_k] > left_bdry_length
	   && (k - j < max_group_size || word_block[last_k] >= maxmatch); 
	   last_k = k, k = next_lower[k - i_lower]);
      
      w_length = word_block[last_k] > left_bdry_length ? word_block[last_k] : left_bdry_length;
      w_length++;
      group_histogram[last_k > j + 1 ? 2 : (last_k > j ? 1 : 0)][w_length] += 1;

      histogram[last_k >= j + 200 ? 200 : last_k - j + 1] += 1;
      left_bdry_length = minmatch - 1 > word_block[last_k] ? 
	minmatch - 1 : word_block[last_k] ;
      
      for (k = j; k <= last_k; k++) {
	word_block[k] = w_length;
      }
      word_block[j] |= 128; /* left bdry flag */
      j = last_k;
    }
  }
  notify(" Done\n\n");

  fprintf(stderr, "\nHistogram of no. of %d-mers with given no. of occurrences", 
	  HIST_WORD_SIZE);
  for (i = 0; i <= 200; i++) 
    if (hist20[i]) 
      fprintf(stderr, "\n%3d  %6d", i, hist20[i]);

  fprintf(stderr, "\n\nHistogram of no. of word groups of given length (group sizes 1, 2, and > 2)\n");
  for (i = 0; i <= 200; i++) {
    if (group_histogram[0][i] || group_histogram[1][i] || group_histogram[2][i]) 
      fprintf(stderr, "%3d  %6d  %6d  %6d\n", 
	      i, group_histogram[0][i], group_histogram[1][i], group_histogram[2][i]);
  }

  fprintf(stderr, "\n\nHistogram of no. of word groups of given size\n");
  for (i = 0; i <= 200; i++) {
    if (histogram[i]) 
      fprintf(stderr, "%3d  %6d\n", i, histogram[i]);
  }



  our_free(next_lower_length);
/* Check whether sort is correct
  for (i = 1; i < s_length; i++)
    if (compare_words(seq_area + s_array[i], seq_area + s_array[i - 1], seq_area + s_array[i] + maxmatch, &last_j) < 0) {
      fprintf(stderr, "\n%d\n", i);
      fprintf(stderr,"%s",seq_area + s_array[i - 1]);
      fprintf(stderr, "\n");
      fprintf(stderr,"%s",seq_area + s_array[i]);
      fatalError("Sort");
    }
*/
/*
  print_n_compares();
  exit(1);
*/
}     

compare_positions(position1, position2)
     int *position1, *position2;
{
  char *seq1, *seq2;
  int k;

  seq1 = seq_area + *position1;
  seq2 = seq_area + *position2;
  for (k = 0; k < minmatch; k++, seq1++, seq2++)
    if (*seq1 != *seq2) return *seq1 - *seq2;
  return 0;
}


/* find and sort all words of length minmatch, not containing N, X, 
   or a lower case letter */

sort_words(t_length, minmatch0, degen_chars)
     int t_length;
     int minmatch0;
     char *degen_chars; /* degenerate characters -- normally "NX" for nucleotides
			   and "X" for proteins */
{
  char *our_alloc();
  int j, k, m, last_j;
  char c;

  notify("Finding words ...");

  if (strlen(degen_chars) > 2) 
    fatalError("Only two degenerate characters permitted in call to sort_words");
  d_char1 = degen_chars[0];
  d_char2 = degen_chars[1]; /* may be 0 */

  minmatch = minmatch0;

  s_length = t_length;
  last_j = t_length + num_entries + 1 - minmatch;
/*  fprintf(stderr, "\n%s\n", seq_area + last_j); */
/* count no. of words -- to allow minimum allocation */
/* N.B. FOLLOWING DOES NOT EXCLUDE WORDS IN TRIMMED PART OF SEQUENCE! -- CHANGE */
  for (m = j = k = 0; ; j++) {
    for (; m < j + minmatch; m++) {
      c = seq_area[m];
      if (c == d_char1 || c == d_char2 || !c) {
	j = m + 1; /*  || islower(c) */
	if (j >= last_j) goto out1; /* always reach this point -- when c is final 0 */
      }
    }
    k++;
  }

 out1:
  s_length = k;
  s_array = (int *)our_alloc((k + 1) * sizeof(int));
  s_array[0] = 0; /* create sentinel -- word starting with 0, so less than
		     every other word */
  s_array++;

  for (m = j = k = 0; ; j++) {
    for (; m < j + minmatch; m++) {
      c = seq_area[m];
      if (c == d_char1 || c == d_char2 || !c) {
	j = m + 1; /*  || islower(c) */
      }
    }
    s_array[k++] = j;
    if (k >= s_length) break;
  }

  notify(" Done\n");

  fprintf(fp_log, "\n\nNo. words: %d; after pruning: %d\n", 
	  t_length, s_length);

  notify("Sorting words ...");

/*  qsort(s_array, s_length, sizeof(int), compare_positions);
*/
  quicksort(0, s_length - 1, 0); 
/*
  for (i = 1; i < s_length; i++)
    if (compare_words(seq_area + s_array[i], seq_area + s_array[i - 1], seq_area + s_array[i] + minmatch, &last_j) < 0) {
      printf("\nError: %d", i);
      fatalError("Sort");
    }
*/
  notify(" Done\n");
/*
  print_n_compares();
  exit(1);
*/
}     

/*
static int n_compares, t_length;

print_n_compares()
{
  printf("\n%d comparisons, t_length %d, %.1f\n", 
	 n_compares, t_length, t_length/(float)n_compares);
}
*/

compare_words(seq1, seq2, w_stop, k_end)
     char *seq1, *seq2;
     char *w_stop;
     int *k_end;
{
  char *w_start;

/*
  n_compares++;
*/
  w_start = seq1;
  for (; seq1 < w_stop; seq1++, seq2++)
    if (*seq1 != *seq2) {
      *k_end = seq1 - w_start;
/*
      t_length += *k_end + 1;
*/
      return *seq1 - *seq2;
    }

  *k_end = seq1 - w_start;
/*
  t_length += *k_end + 1;
*/
  return 0;
}

/* quicksort (based in part on generic recursive version with insertion sort
   refinement in Sedgwick, but using bookkeeping of matching leading parts 
   of words in a block to substantially improve
   efficiency) */

quicksort(left, right, k_start) 
     int left, right, k_start;
{
  char *w_right, *w_stop, *s_plus_k;
  int i, j, temp;
  int k_end, k_stop, k_left_start, k_right_start;

/* use insertion sort for small pieces */
  k_stop = maxmatch - k_start;
  s_plus_k = seq_area + k_start;
  if (right < left + 12) {
    for (j = left + 1; j <= right; j++) {
      temp = s_array[j];
      w_right = s_plus_k + temp;
      w_stop = w_right + k_stop;
      for (i = j; i > left &&
	   compare_words(w_right, s_plus_k + s_array[i - 1], w_stop, &k_end) < 0; 
	   i--) 
	s_array[i] = s_array[i - 1];
      s_array[i] = temp;
    }
    return;
  }
      
  w_right = s_plus_k + s_array[right];
  w_stop = w_right + k_stop;
  i = left - 1;
  j = right;
  k_left_start = k_right_start = k_stop;
  for (;;) {
    for (i++; ; i++) {
      if (compare_words(w_right, s_plus_k + s_array[i], w_stop, &k_end) <= 0) {
	if (k_end < k_right_start) k_right_start = k_end;
	break;
      }
      else if (k_end < k_left_start) k_left_start = k_end;
    }
    for (j--; j >= left; j--) {
      if (compare_words(w_right, s_plus_k + s_array[j], w_stop, &k_end) >= 0) {
	if (k_end < k_left_start) k_left_start = k_end;
	break;
      }
      else if (k_end < k_right_start) k_right_start = k_end;
    }
    if (j < i) break;
    temp = s_array[i];
    s_array[i] = s_array[j];
    s_array[j] = temp;
  } 
  temp = s_array[i];
  s_array[i] = s_array[right];
  s_array[right] = temp;
  /* leading k_left_start characters match for all words in left part of this block,
     similarly for k_right_start */
  k_left_start += k_start;
  if (k_left_start < maxmatch && left < j) 
    quicksort(left, j, k_left_start);
  k_right_start += k_start;
  if (k_right_start < maxmatch && i + 1 < right) 
    quicksort(i + 1, right, k_right_start);
}

static int n_cross_matches, n_same_matches;

find_subject_matches(files)
     File *files;
{
  File *file;
  int entry1;
  Aligned_pair *get_aligned_pairs();
  char *seq;
  Database *sdb;
  Segment *seg_list;
  Segment *repeat_parse_descrip();

  notify("Finding subject matches ...");


  num_pass_words = num_le_fail_words = num_repeat_fail_words = 0;
  for (file = files; file; file = file->next) {
    sdb = file->db;
/* Note -- this assumes file is not disturbed!! */
    while (get_next_file_entry(sdb)) {
      entry1 = append_seq_entry(sdb);
      seg_list = repeat_screen ? repeat_parse_descrip(sdb->descrip_buffer) : (Segment *)0;
      find_external_word_matches(entry1, seg_list, sdb->seq_buffer);
/* in this case, will not be doing comp or internal matches; need to generate list
   of cand. pairs, and do swat alignments immediately */
      find_scores(entry1, sdb->seq_buffer);
      if (!get_aligned_pairs(entry1)) remove_seq_entry(sdb);
      free_cand_pair_blocks();
      free_seg_blocks();
    }
  }
  notify(" Done\n");
  fprintf(stderr, "# pass words: %d, #low-entropy-fail words: %d, #repeat_fail_words: %d\n", 
	  num_pass_words, num_le_fail_words, num_repeat_fail_words);

}

find_comp_matches()
{
  int entry1;
  char *seq;
  char *get_comp_seq();

  notify("Finding complement word matches ...");

  num_pass_words = num_le_fail_words = num_repeat_fail_words = 0;
  for (entry1 = start_entry; entry1 < start_entry + num_uncomp_entries; entry1++) {
    seq = get_comp_seq(entry1);
    find_external_word_matches(entry1, (Segment *)0, seq);

/* in this case, will not be doing subject matches; create candidate pairs to process
   in the usual fashion */

  }
  notify(" Done\n");
  fprintf(stderr, "\n %d opp_sense matches:  %d pass words,  %d low-entropy-fail words,  %d repeat-fail words\n", 
	  n_cross_matches, num_pass_words, num_le_fail_words, num_repeat_fail_words);
  print_num_cand_pairs();
}

find_external_word_matches(entry1, seg_list, seq)
     int entry1;
     Segment *seg_list;
     char *seq;
{
  char *seqm, *seqj, *word_end;
  int length, last_j, j, word_int, w;

  length = get_seq_length(entry1);
  last_j = length - minmatch;
  seqm = seq - 1;

 reinit2:
  seqj = seqm + 1;
  if (seqj > seq + last_j) {
    return; /* always reach this point -- when c is final 0 */
  }
  word_int = 0;
  word_end = seqj + index_word_size;
  for (seqm = seqj; seqm < seqj + minmatch; seqm++) {
    w = residues[*seqm];
    if (w < 0) goto reinit2;
    if (seqm < word_end) word_int = word_int * alphabet_size + w;
  }
  j = seqj - seq;
  new_lookup_words(entry1, seg_list, j, word_int, word_end, seqm);
  
  if (alphabet_size == 4) {
    for (; ; ) {
      if (residues[*seqm] < 0) goto reinit2;
      word_int = ((word_int << 2) & mask) + residues[*word_end];
/*
    word_int = (word_int % mod_value) * alphabet_size + residues[*word_end];
*/
      j++;
      seqm++;
      word_end++;
      new_lookup_words(entry1, seg_list, j, word_int, word_end, seqm);
    }
  }
  else {
    for (; ; ) {
      if (residues[*seqm] < 0) goto reinit2;
      word_int = (word_int - m_residues[*(word_end - index_word_size)]) * alphabet_size + residues[*word_end];
/*
    word_int = (word_int % mod_value) * alphabet_size + residues[*word_end];
*/
      j++;
      seqm++;
      word_end++;
      new_lookup_words(entry1, seg_list, j, word_int, word_end, seqm);
    }
  }
}

new_lookup_words(entry1, seg_list, offset1, word_int, seqw, seqm)
     int entry1, offset1, word_int;
     Segment *seg_list;
     char *seqw, *seqm;
{
  int i, k, k_end, l_lower, u_upper, mid, c, found;
  int entry2, offset2, complement, msize;
  char *seq1, *seq2;

  l_lower = index_words[word_int];
  u_upper = index_words[word_int + 1] - 1;
  if (u_upper < l_lower) return;

  seq1 = seqw - index_word_size;
  found = 0;
  do {
    mid = (u_upper + l_lower + 1) / 2;
	    seq2 = seq_area + s_array[mid];
    for (k = index_word_size; 
	 seq1[k] == seq2[k] && seq1[k] && k < maxmatch; k++);
    if (k >= (word_block[mid] & 127)) {
      found = 1;
      break;
    }
    if (seq1[k] < seq2[k]) u_upper = mid - 1;
    else l_lower = mid + 1;
  } while (u_upper >= l_lower);
  if (!found) return;

/* following is inefficient; instead should mimic procedure in old
  lookup_words */

  for (l_lower = mid; 
       !(word_block[l_lower] & 128); 
       l_lower--);
  for (u_upper = mid + 1; 
       u_upper < index_words[word_int + 1]  && 
       !(word_block[u_upper] & 128); 
       u_upper++);
  u_upper--;

  c = residues[seqm[-minmatch - 1]];
/*
  if (compare_words(seqm - minmatch, seq_area + s_array[l_lower], seqm, &k_end)
      || l_lower >= 1 && !compare_words(seqm - minmatch, seq_area + s_array[l_lower - 1], seqm, &k_end)
      || compare_words(seqm - minmatch, seq_area + s_array[u_upper], seqm, &k_end)
      || !compare_words(seqm - minmatch, seq_area + s_array[u_upper + 1], seqm, &k_end))
    fatalError("wordmatches");
*/
  for (i = l_lower; i <= u_upper; i++) {

    if (c < 0 || c != residues[seq_area[s_array[i] - 1]]) {

      find_entry_and_offset(s_array[i], &entry2, &offset2, &complement);
      if (entry2 <= entry1 && !parameters->subject_files) continue; 
/* use only one of two equivalent (complementary) word matches, and
				     don't allow matches to reverse complement
				     of self for now; fix this later (requires
				     attention to num_pairs) */

      n_cross_matches++;
      msize = get_matchsize(seqm - minmatch, seq_area + s_array[i]);
      if (complexity_flag && get_corrected_matchsize() < minmatch) {
	num_le_fail_words++;
	continue;
      }

      if (parameters->subject_files) {
	if (repeat_screen &&
	    (contained_in_tag(entry2, "repeat", offset2 + 1, offset2 + msize, complement)
	     || seg_list && find_max_gap(seg_list, offset1 + 1, offset1 + msize) < 15)
	    ) {
	  num_repeat_fail_words++;
	  continue;
	}
	pair_routine(entry1, entry2, offset1, offset2, complement);
	num_pass_words++;
      }
      else { /* no subject file, so must be complement of query file sequence */
	if (repeat_screen && 
	    (contained_in_tag(entry1, "repeat", offset1 + 1, offset1 + msize, 1)
	    || contained_in_tag(entry2, "repeat", offset2 + 1, offset2 + msize, 0))) {
	  num_repeat_fail_words++;
	  continue;
	}
	pair_routine(entry2, entry1, offset2, offset1, 1);
	num_pass_words++;
      }
    }
  }
}

old_lookup_words(entry1, offset1, word_int, seqw, seqm)
     int entry1, offset1, word_int;
     char *seqw, *seqm;
{
  int i, k_end, l_lower, l_upper, u_lower, u_upper, mid, c, found;
  int entry2, offset2, complement;

  l_lower = u_lower = index_words[word_int];
  l_upper = u_upper = index_words[word_int + 1] - 1;
  if (l_upper < l_lower) return;
  found = 0;
  c = compare_words(seqw, seq_area_plus + s_array[l_lower], seqm, &k_end); 
  if (c < 0) return;
  if (c == 0) {
    l_upper = l_lower;
    found = 1;
  }
  else l_lower++;

  if (u_lower < u_upper) {
    c = compare_words(seqw, seq_area_plus + s_array[u_upper], seqm, &k_end); 
    if (c > 0) return;
    if (c == 0) {
      u_lower = u_upper;
      found = 1;
    }
    else {
      u_upper--;
      if (l_upper > u_upper) l_upper = u_upper;
    }
  }


  while (l_upper > l_lower || l_upper == l_lower && !found) {
    mid = (l_upper + l_lower) / 2;
    c = compare_words(seqw, seq_area_plus + s_array[mid], seqm, &k_end); 
    if (c > 0) l_lower = mid + 1;
    else if (c < 0) l_upper = u_upper = mid - 1;
    else {
      l_upper = mid;
      if (!found) {
	found = 1;
	u_lower = mid;
      }
    }
  }
  if (!found) return;
  if (u_lower < l_lower) u_lower = l_lower;

  while (u_upper > u_lower) {
    mid = (u_upper + u_lower + 1) / 2;
    c = compare_words(seqw, seq_area_plus + s_array[mid], seqm, &k_end); 
    if (c < 0) u_upper = mid - 1;
    else u_lower = mid;
  } 

  c = residues[seqm[-minmatch - 1]];
/*
  if (compare_words(seqm - minmatch, seq_area + s_array[l_lower], seqm, &k_end)
      || l_lower >= 1 && !compare_words(seqm - minmatch, seq_area + s_array[l_lower - 1], seqm, &k_end)
      || compare_words(seqm - minmatch, seq_area + s_array[u_upper], seqm, &k_end)
      || !compare_words(seqm - minmatch, seq_area + s_array[u_upper + 1], seqm, &k_end))
    fatalError("wordmatches");
*/
  for (i = l_lower; i <= u_upper; i++) {
    if (c < 0 || c != residues[seq_area[s_array[i] - 1]]) {
      find_entry_and_offset(s_array[i], &entry2, &offset2, &complement);
      get_matchsize(seqm - minmatch, seq_area + s_array[i]);
   if (!complexity_flag || get_corrected_matchsize() >= minmatch) {
	num_pass_words++;
	if (parameters->subject_files) {
	  pair_routine(entry1, entry2, offset1, offset2, complement);
	  n_cross_matches++;
	}
	else if (entry2 > entry1) { /* use only one of two equivalent (complementary) word matches, and
				       don't allow matches to reverse complement
				       of self for now; fix this later (requires
				       attention to num_pairs) */
	  pair_routine(entry2, entry1, offset2, offset1, 1);
	  n_cross_matches++;
	}
      } 
      else num_le_fail_words++;
    }
  }
}

old_find_internal_word_matches()
{
  char *start_gp, *start2;
  int offset1, offset2, entry1, entry2, msize;
  int i, j, k, m, c, i_upper, group_start, k_begin, k_end;
  int phrap_flag, complement;

  notify("Finding internal word matches ...");

  k_begin = minmatch - 1;
  k_end = index_word_size; 

  phrap_flag =  !strcmp(parameters->calling_program, "phrap") || !strcmp(parameters->calling_program, "gcphrap") ;
  for (m = 0; m < num_index_words; m++) { 
    i_upper = index_words[m + 1];
    group_start = index_words[m];
    start_gp = seq_area + s_array[group_start];
    for (i = group_start + 1; i < i_upper; i++) {
      start2 = seq_area + s_array[i];
      for (k = k_begin; start_gp[k] == start2[k] && k >= k_end; k--);
/* start at end, where mismatches are most likely; and only check some of the
   bases */
      if (k < k_end) {
	entry2 = -1;
	c = residues[start2[-1]];
	for (j = group_start; j < i; j++) {
	  if (c < 0 || residues[seq_area[s_array[j] - 1]] != c) {
	    /* only consider full-length matches */
	    if (entry2 < 0) find_entry_and_offset(s_array[i], &entry2, &offset2, &complement);
	    find_entry_and_offset(s_array[j], &entry1, &offset1, &complement);
/* ignore nucleotide word matches entirely contained within first 
   parameters->vector_bound residues of both reads  */
	    msize = get_matchsize(seq_area + s_array[i], seq_area + s_array[j]);
	    if (phrap_flag && offset1 + minmatch  
		< parameters->vector_bound && offset2 + minmatch < parameters->vector_bound) {
	      for (k = minmatch; start2[k] == start_gp[k]; k++);
	      if (k + offset1 < parameters->vector_bound) {
		vector_matches++;
		continue; 
	      }
	    }
	    if (!complexity_flag 
		|| get_corrected_matchsize() >= minmatch) {
	      if (entry1 > entry2) pair_routine(entry1, entry2, offset1, offset2, 0);
	      else pair_routine(entry2, entry1, offset2, offset1, 0);
	      n_same_matches++;
	    }
	  }
	} 
      }
      else {
	group_start = i;
	start_gp = start2;
      }
    }
  }
  fprintf(stderr,"  n_same_matches: %d\n", n_same_matches);
  print_num_cand_pairs();
  notify(" Done\n");
}

new_find_internal_word_matches()
{
  char *start_gp, *start2, *start3;
  int offset1, offset2, entry1, entry2, msize;
  int i, j, k, m, c, i_upper, group_start;
  int vector_bound, complement;


  notify("Finding internal word matches ...");

  num_pass_words = num_le_fail_words = num_repeat_fail_words = 0;

  vector_bound = parameters->vector_bound;
  for (m = 0; m < num_index_words; m++) { 

    i_upper = index_words[m + 1];
    group_start = index_words[m];
/*
    if (group_start > alloc_size - 1)
      fatalError("illegal s_array read");
*/
    start_gp = seq_area + s_array[group_start];
    for (i = group_start + 1; i < i_upper; i++) {
      start2 = seq_area + s_array[i];
      if (!(word_block[i] & 128)) {
	entry2 = -1;

	c = residues[start2[-1]];
	for (j = group_start; j < i; j++) {

	  start3 = seq_area + s_array[j];
	  if (c < 0 || residues[start3[-1]] != c) {
	    /* only consider full-length matches */
	    if (entry2 < 0) 
	      find_entry_and_offset(s_array[i], &entry2, &offset2, &complement);
	    find_entry_and_offset(s_array[j], &entry1, &offset1, &complement);
/* ignore nucleotide word matches entirely contained within first 
   parameters->vector_bound residues of both reads  */
	    n_same_matches++;
	    msize = get_matchsize(start2, start3);
	    /* note assymetry in condition -- should this be changed? */
	    if (offset1 + msize < vector_bound && offset2 + minmatch < vector_bound) {
	      vector_matches++;
	      continue; 
/* old version -- gives slightly different results because "reads thru"
		 X's, N's etc.
	      for (k = minmatch; start2[k] == start3[k]; k++);
	      if (k + offset1 < vector_bound) {
		vector_matches++;
		continue; 
	      }
*/
	    }

	    if (complexity_flag && get_corrected_matchsize() < minmatch) {
	      num_le_fail_words++;
	      continue;
	    }
	    if (repeat_screen && 
		(contained_in_tag(entry1, "repeat", offset1 + 1, offset1 + msize, 0)
		 || contained_in_tag(entry2, "repeat", offset2 + 1, offset2 + msize, 0))) {
	      num_repeat_fail_words++;
	      continue;
	    }

	    num_pass_words++;
	    if (entry1 > entry2) 
	      pair_routine(entry1, entry2, offset1, offset2, 0);
	    else 
	      pair_routine(entry2, entry1, offset2, offset1, 0);
	  }
	} 
      }
      else {
	group_start = i;
	start_gp = start2;
      }
    }
  }
  fprintf(stderr,"\n  %d same-sense word matches:  %d vector_matches,  %d pass words, %d low-entropy-fail words, %d repeat-fail words\n", 
	  n_same_matches, vector_matches, num_pass_words, num_le_fail_words, num_repeat_fail_words);
  print_num_cand_pairs();
  notify(" Done\n");
}

find_entry_and_offset(position, add_entry, add_offset, add_complement)
     int position;
     int *add_entry, *add_offset, *add_complement;
{
  unsigned int upper, lower, mid;

  upper = num_entries;
  lower = 0;
  do {
    mid = (upper + lower) >> 1;
    if (position_table[mid] > position)
      upper = mid;
    else lower = mid;
  } while (upper > lower + 1);
  if (lower >= num_uncomp_entries) {
    *add_entry = num_entries - 1 - lower;
    *add_complement = 1;
  }
  else {
    *add_entry = lower;
    *add_complement = 0;
  }
  *add_offset = position - position_table[lower];
  *add_entry += start_entry;
}

free_word_arrays()
{

  our_free(s_array - 1); /* N.B. -- s_array no longer available! */
  our_free(position_table);
/* temp */
  our_free(index_words);
  our_free(word_block);
}

static double nlogn[500];
static int res_counts[256];
static double p;

make_nlogn()
{
  int n;

  p = -log((double)alphabet_size);
  nlogn[0] = 0;
  for (n = 1; n < 500; n++)
    nlogn[n] = n * log((double)n) / p;
}

int get_matchsize(seq1, seq2)
     char *seq1, *seq2;
{
  int i, r1, r2, r;

  for (i = 0; i < alphabet_size; i++) res_counts[i] = 0;

  for (r = 0;;) {
    r1 = residues[*seq1];
    r2 = residues[*seq2];
    if (r1 == r2 && r1 >= 0) {
      res_counts[r1] += 1;
      seq1++;
      seq2++;
      r++;
    }
    else break;
  }
  if (r < minmatch) fatalError("matchsize < minmatch");
  return r;
}

int get_corrected_matchsize()
{
  double size, x;
  int i, c, n;

  size = n = 0;
  for (i = 0; i < alphabet_size; i++) {
    c = res_counts[i];
    n += c;
    x = c < 500 ? nlogn[c] : c * log((double)c) / p; 
    size += x;
  }
  size -= n < 500 ? nlogn[n] : n * log((double)n) / p; 
  return (int)(size + .999); /* round up to next highest integer */
}
