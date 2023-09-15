/*****************************************************************************
#   Copyright (C) 1994-2009 by Phil Green.                          
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

int merge_residues[256]; /* used in print_alignments (for parameters->DNA_flag == 2) -- nowhere else (outside this file) at present */

static FILE *fp_debug; /* normally stderr */
static unsigned char *seq_area, *comp_area; /* array containing sequence data for all entries (see db.c) */
static SEQ_AREA seq_area_size;
static int minmatch, maxmatch, word_offset; /* from parameters */
static double x_minmatch;
static int num_entries, start_entry, num_uncomp_entries, merged_flag;
static int max_length;
static SEQ_AREA total_length;
static SEQ_AREA *position_table;
static int residues[256], *trans_index_residues[2], *mod_val_trans_index_residues[2], *trans_sort_residues[2], *trans_index_res, *mod_val_trans_index_res, *trans_sort_res, *mod_val_residues[2], *mod_val_res; /* residues gives numerical conversion of genomic, trans_residues of query */
static int alphabet_size, num_index_words[2], index_word_size[2], index_shift[2], mod_value[2];
static SEQ_AREA s_array_length[2];
static SEQ_AREA *index_words[2];
static SEQ_AREA *s_array[2]; /* arrays of index words, positions of word starts */
static SEQ_AREA *index_lower, *index_upper;
static unsigned char *word_block[2];
static int (*pair_routine)();
static int complexity_flag, q_repeat_screen, s_repeat_screen;
static double num_pass_words, num_le_fail_words, num_le_test_words, num_s_repeat_fail_words, num_q_repeat_fail_words, vector_matches, n_cross_matches, n_same_matches, num_merge_words[5];

static double nlogn[500], nlogn_diff[500];
static int res_counts[256], base_res_counts[256];
static double clump_hist[2][101], match_size_hist[101];
static int *trans_sort_vec[2], *trans_sort_v;
static int gap_init, dropoff, long_word_shift;
static unsigned char *long_word_array;
static int n_skips[3];


typedef struct sw_edge {
  int curr0, curr1, curr2;
  int max_score;
} Sw_edge;


/*
  long_word_array -- allows quick lookup whether word of length minmatch is present in query 
    (saves time when processing genome -- since looking up query sequences is expensive)

  rel probs of matches in merged base data. True base corresponding to stronger signal should match stronger peak.
   True base corresponding to weaker signal should match
   (i) single called base if there is only one. 
   (ii) 2d base -- if both are 'real'     1 - .25 z
   (iii) 1st base -- if true signals are same, 2d peak is noise  .25 z (z = noise prob)
  so (cond'l on observed no. of peaks) likelihood ratio rel to random base = 
    (i) 1/.25 = 4
    (ii) (1 - .25 z) / .25 = 4 - z
    (iii) .25 z / .25 = z

  have option for only taking strongest signal, in merged read data??  also, what about merged reads
   with single input file (i.e. no subject files?)

  repeat markup: might be preferable to do only for near-exact repeats (not distinguishable).

 check usage of complexity_flag -- where is it only appropriate for primary peaks?

 for merged base reads, want to find word matches as follows:
   ordinary matches to top peaks
   matches that start with several bottom peaks, then extend by matching either top or bottom
   representation: in general a single subject word will match multiple ambiguity codes. 
   Order codes as follows (for this purpose, not distinguishing case, i.e. which peak is larger)
     A  A
     B  AC
     C  AG
     D  AT
     E  T
     F  CT
     G  GT    
     H  G
     I  CG
     J  C 

     Then A entries are 0-3 (so one entry point, at 0)
          T are 3-6         (1 entry)
          C are 1, 5, 8-9   (3 entries)
          G are 2, 6-8      (2 entries)
      so 7 different entry points
    so fewer entries for A+T rich words, since these are most common in most genomes

    sort pointers in two ways: 
           1. only considering n top peaks
           2. bottom peaks for first k positions, then above order on ambig codes for m subsequent positions (k + m = n).

    for a given word, need list of corresponding merged read words to start search at (then end when have mismatch at some position)
     so have list of 4^n word indices; each having multiple entry points. total no. of entry points is 4^k 7^m. e.g. for k = m = 5
     this is ~1000 50^2 7 ~ 17.5 million. Not correct!!
     
     simpler: first occurrence of index word, and last occurrence of index word, allowing some intervening non-occurrences
       Do sequence lookups to find 1st, last occurrences (starting with preword of length k)
       A >= A, <= D.
       C >= B, <= J
       G >= C, <= I
       T >= D, <= G

      compare word at each position, 
*/

typedef struct ssite { /* potential splice site locations */
  int pos, keep_length;
  unsigned int word_int, long_word_int;
} Ssite;

static int n_ssites, first_ssite[2], last_ssite[2], max_intron, max_intron2, min_intron, indwd2, indwd22, left_ss_wdsize, left_ss_wdshift, left_ss_wdmask, right_ss_wdsize, right_ss_wdshift, right_ss_wdmask, right_ss_long_wdshift, right_ss_long_wdmask;
static double n_splice_calls;
static Ssite *ssites[2];
static int *off1s, i_off1, n_off1, spliced_match_left, spliced_match_right;

set_word_db(db, complements)
     Database *db;
     int complements;
{
  char *our_alloc();
  int i, j, n_long_words, long_word_size, right_ss_long_wdsize;
  int make_new_cand_pair(), cluster_pairs();
  char alphabet[30];

  fp_debug = parameters->print_word_data ? stdout : 0; /* stderr; */

  for (i = 0; i < 3; i++)
    n_skips[i] = 0;
  for (i = 0; i < 101; i++) clump_hist[0][i] = clump_hist[1][i] = match_size_hist[i] = 0;

  seq_area = db->seq_area;
  comp_area = db->comp_area;
  seq_area_size = db->seq_area_size;
  num_uncomp_entries = db->num_entries;
  num_entries = complements ? 2 * num_uncomp_entries: num_uncomp_entries;
  total_length = complements ? 2 * db->t_length : db->t_length;

  max_length = 1 << 30;

  position_table = (SEQ_AREA *)our_alloc((num_entries + 1) * sizeof(SEQ_AREA));
  start_entry = db->first_entry;
  position_table[0] = 1;
  for (i = 0; i < num_uncomp_entries; i++) 
    position_table[i + 1] = position_table[i] + get_seq_length(start_entry + i) + 1;

  for (; i < num_entries; i++) 
    position_table[i + 1] = position_table[i] + 
      get_seq_length(start_entry + num_entries - 1 - i) + 1;
  /* check */
  for (i = 0; i < num_entries; i++) 
    if (seq_area[position_table[i] - 1]) {
      fprintf(stderr, "\n%d %d", i, num_uncomp_entries);
      fatalError("position_table");
    }

  /* following not correct -- for phrap assemblies
  if (seq_area_size != position_table[num_entries]) {
    fprintf(stderr, "\n%lu %lu %d %d", (unsigned long int)seq_area_size, (unsigned long int)position_table[num_entries], num_entries, complements);
    fatalError("position_table");
  }
  */

  index_word_size[0] = parameters->indexwordsize;
  index_word_size[1] = parameters->indexwordsize2; 

  gap_init = parameters->gap_init;
  dropoff = parameters->gap1_dropoff;
  minmatch = parameters->minmatch;
  x_minmatch = minmatch - .999; /* in complexity adjustment -- to round up */
  maxmatch = parameters->maxmatch;
  word_offset = parameters->word_offset;

  for (i = 0; i < 256; i++) residues[i] = -2; /* -2 is minmatch stopper; -1 is indexword stopper */

  strcpy(alphabet, parameters->DNA_flag ? "ACGT" : "ACDEFGHIKLMNPQRSTVWY");
  alphabet_size = strlen(alphabet);

  for (i = 0; i < alphabet_size; i++)
    residues[alphabet[i]] = residues[tolower(alphabet[i])] = i;

  residues[' '] = -1;

  long_word_size = alphabet_size == 4 && minmatch > index_word_size[0] ? minmatch : 0;
  if (long_word_size > 14) long_word_size = 14;
  if (long_word_size) {
    long_word_shift = 2 * (long_word_size - index_word_size[0]); 
    if (long_word_shift < 0) long_word_size = long_word_shift = 0;
    else {
      for (n_long_words = 1, i = 0; i < long_word_size; n_long_words *= 4, i++);
      n_long_words /= 8;
      long_word_array = (unsigned char *)our_alloc(n_long_words * sizeof(unsigned char));
      for (i = 0; i < n_long_words; i++) long_word_array[i] = 0;
    }
  }
  else long_word_shift = 0;

  merged_flag = parameters->DNA_flag == 2; 

  if (merged_flag) {
    for (i = 0; i < 256; i++) merge_residues[i] = 0;

    set_merge_residues('A', 1, 1);
    set_merge_residues('C', 2, 2);
    set_merge_residues('G', 4, 4);
    set_merge_residues('T', 8, 8);
    set_merge_residues('M', 1, 2);
    set_merge_residues('R', 1, 4);
    set_merge_residues('W', 1, 8);
    set_merge_residues('S', 2, 4);
    set_merge_residues('Y', 2, 8);
    set_merge_residues('K', 4, 8);
    set_merge_residues('N', 15, 15);
  }
  max_intron = parameters->spliced_word_gapsize;
  max_intron2 = parameters->spliced_word_gapsize2;

  if (max_intron || max_intron2) {
    n_ssites = .5 * (max_intron2 > max_intron ? max_intron2 : max_intron) + 1;
    ssites[0] = (Ssite *)our_alloc(n_ssites * sizeof(Ssite));
    ssites[1] = (Ssite *)our_alloc(n_ssites * sizeof(Ssite));

    spliced_match_left = parameters->spliced_match_left;/* 0; */
    spliced_match_right = parameters->spliced_match_right;/* 16; */

    if (n_off1 = .5 * max_intron2)
      off1s = (int *)our_alloc(n_off1 * sizeof(int));

    min_intron = parameters->min_intron_length;

    indwd2 = index_word_size[0] - 2;
    indwd22 = 2 * indwd2;
    left_ss_wdsize = minmatch / 2; /* round up? */
    right_ss_wdsize = index_word_size[0] - left_ss_wdsize;
    if (right_ss_wdsize <= 0) fatalError("incompatible splice site word sizes"); /* instead put this in parameters.c */

    left_ss_wdshift = 2 * right_ss_wdsize;
    left_ss_wdmask = pow(2.0, 2.0 * left_ss_wdsize) - 1;

    right_ss_wdshift = 2 * (index_word_size[0] - (right_ss_wdsize + 2));
    right_ss_wdmask = pow(2.0, 2.0 * right_ss_wdsize) - 1;

    if (long_word_size) {
      right_ss_long_wdsize = long_word_size - left_ss_wdsize;
      right_ss_long_wdshift = 2 * (index_word_size[0] - (right_ss_long_wdsize + 2));
      if (right_ss_long_wdshift < 0) fatalError("negative shift");
      right_ss_long_wdmask = pow(2.0, 2.0 * right_ss_long_wdsize) - 1;
    }
  }

  for (i = 0; i < 2; i++) {
    trans_sort_residues[i] = (int *)our_alloc(256 * sizeof(int));
    trans_index_residues[i] = (int *)our_alloc(256 * sizeof(int));
    mod_val_trans_index_residues[i] = (int *)our_alloc(256 * sizeof(int));
    trans_sort_vec[i] = (int *)our_alloc(256 * sizeof(int));
  }

  for (j = 0; j < 2; j++) {
    index_shift[j] = 2 * (index_word_size[j] - 1);
    for (i = 0, mod_value[j] = 1; i < index_word_size[j] - 1; i++, mod_value[j] *= alphabet_size); 
    num_index_words[j] = mod_value[j] * alphabet_size;
    mod_val_residues[j] = (int *)our_alloc(256 * sizeof(int));
    for (i = 0; i < 256; i++) mod_val_residues[j][i] = residues[i] < 0 ? residues[i] : mod_value[j] * residues[i];
  }

  if (!strcmp(parameters->calling_program, "cluster"))
    pair_routine = cluster_pairs;
  else pair_routine = make_new_cand_pair;  

  q_repeat_screen = parameters->repeat_screen & 1;
  s_repeat_screen = parameters->repeat_screen & 2;

  complexity_flag = !parameters->word_raw;

  make_nlogn();

  if (merged_flag) { 
    new_sort_words(0, 1, 1, "AaCcGgTtMmRrWwYyKkSs", "AACCGGTTACAGATCTGTCG", "AACCGGTTACAGATCTGTCG");
    new_sort_words(1, 0, 0, "AaCcGgTtMmRrWwYyKkSs", "        CAGATATCTGGC", "AAJJHHEEBBCCDDFFGGII"); /* spaces mean index_word is prevented */
    /* simple sort, using qsort but for comparison translating each suffix sequence in two pieces -- using transalph2 and transalph3 */ 
    set_lower_upper();
  }
  else { /* this case also includes DNA_flag >= 3 -- in which case alphabet is ignored */
    new_sort_words(0, 1, 1, alphabet, alphabet, alphabet);
  }
}

/* find and sort all words of length minmatch or greater -- allowing for maxmatch grouping and complexity 

   ISSUE OF MAXMATCH, MAXGROUPSIZE AND COMPLEXITY IS NONTRIVIAL!! 
*/

set_merge_residues(code, num1, num2)
     int code, num1, num2;
{
  merge_residues[code] = (num1 << 4) + num2;
  merge_residues[tolower(code)] = (num2 << 4) + num1;
}

/* only used when DNA_flag == 2 */

set_lower_upper()
{
  unsigned int word_int, word_reduced, last_word_reduced;
  int i, j, k, c, ind_wd_size, ind_wd_size0, w, dup_flag;
  SEQ_AREA *seq_offsets;
  SEQ_AREA *ind_wds;
  long int l_upper, l_lower, u_upper, u_lower, mid, orig;
  unsigned char *seq2, min[256], max[256];
  int shift, n_change;
  double total_word_matches, restricted_word_matches;
  char *our_alloc(); 

  seq_offsets = s_array[1];
  ind_wds = index_words[1];
  ind_wd_size = index_word_size[1];
  ind_wd_size0 = index_word_size[0];
  trans_sort_res = trans_sort_residues[1];
  trans_sort_v = trans_sort_vec[1];
  index_lower = (SEQ_AREA *)our_alloc(num_index_words[0] * sizeof(SEQ_AREA)); /* now index with original word size, but point to
									     lower, upper bounds */
  index_upper = (SEQ_AREA *)our_alloc(num_index_words[0] * sizeof(SEQ_AREA));

  for (i = 0; i < 256; i++) min[i] = max[i] = 0;

  min['A'] = 'A';
  min['C'] = 'B';
  min['G'] = 'C';
  min['T'] = 'D';
  max['A'] = 'D';
  max['C'] = 'J';
  max['G'] = 'I';
  max['T'] = 'G';
  /* find lower bound -- 1st word in sorted seq area that can match word_int, i.e. anything lower cannot match;
   implies no letter is strictly less than min of corresponding letter in word_int ? */
  /* should do some sort of binary search instead!!! */

  if (ind_wd_size == ind_wd_size0) {
    for (word_int = 0; word_int < num_index_words[0]; word_int++) {
      index_lower[word_int] = ind_wds[word_int];
      index_upper[word_int] = ind_wds[word_int + 1]; /*  - 1 -- but could cause neg result if unsigned */
    }
    return;
  }

  shift = index_shift[0] - index_shift[1] - 2; 

  total_word_matches = restricted_word_matches = 0;
  for (word_int = 0; word_int < num_index_words[0]; word_int++) {
    word_reduced = word_int >> shift;
    dup_flag = 0;
    if (word_int && word_reduced == last_word_reduced) {
      index_lower[word_int] = index_lower[word_int - 1];
      index_upper[word_int] = index_upper[word_int - 1];
      dup_flag = 1;
    }

    last_word_reduced = word_reduced;

    word_reduced >>= 2;

    l_lower = ind_wds[word_reduced];
    u_upper = ind_wds[word_reduced + 1] - 1;
    total_word_matches += u_upper - l_lower + 1;
    if (dup_flag) continue;

    if (ind_wd_size < ind_wd_size0) {

      w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - ind_wd_size)) & 3];
      seq2 = seq_area + ind_wd_size;
      l_upper = u_upper;
      do {
	mid = (l_lower + l_upper) >> 1;
	c = trans_sort_v[seq2[seq_offsets[mid]]];
	if (c >= min[w]) l_upper = mid;
	else l_lower = mid + 1;
      } while (l_lower < l_upper);

      u_lower = l_lower;
      do {
	mid = (u_lower + u_upper + 1) >> 1;
	c = trans_sort_v[seq2[seq_offsets[mid]]];
	if (c <= max[w]) u_lower = mid;
	else u_upper = mid - 1;
      } while (u_lower < u_upper);
      
    }

    index_lower[word_int] = l_lower;
    index_upper[word_int] = u_upper + 1; /* was + 0 */
  }


  if (ind_wd_size >= ind_wd_size0 - 1) return;

  for (word_int = 0; word_int < num_index_words[0]; word_int++) {
    l_lower = index_lower[word_int];
    u_upper = index_upper[word_int] - 1; /* was + 0 */
    for (k = ind_wd_size; k < ind_wd_size0; ) {
      w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - k)) & 3];
      w = merge_residues[w];
      seq2 = seq_area + k;
      for (orig = l_lower; l_lower <= u_upper && !(w & merge_residues[seq2[seq_offsets[l_lower]]]); l_lower++);
      /*     fprintf(stdout, "\n%d %c%c %4d", k, w, c, l_lower - ind_wds[word_reduced]); */
      if (orig < l_lower) k = ind_wd_size;
      else k++;
    }

    for (k = ind_wd_size; k < ind_wd_size0; ) {
      w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - k)) & 3];
      w = merge_residues[w];
      seq2 = seq_area + k;
      for (orig = u_upper; u_upper >= l_lower && !(w & merge_residues[seq2[seq_offsets[u_upper]]]); u_upper--);
      if (orig > u_upper) k = ind_wd_size; /* restart */
      else k++;
    }
    restricted_word_matches += u_upper - l_lower + 1;
    index_lower[word_int] = l_lower;
    index_upper[word_int] = u_upper + 1; /* was + 0 */

    /* NEEDS TYPE CHANGES TO REVIVE
    for (i = 0; i < ind_wds[num_index_words[1]]; i++) {
      seq2 = seq_area + seq_offsets[i];
      for (k = 0; k < ind_wd_size; k++) {
	w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - k)) & 3];
	w = merge_residues[w];
	if (!(w & merge_residues[seq2[k]] & 15)) break;
      }
      if (k < ind_wd_size) continue;
      for (k = ind_wd_size; k < ind_wd_size0; k++) {
	w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - k)) & 3];
	w = merge_residues[w];
	if (!(w & merge_residues[seq2[k]])) break;
      }
      if (k == ind_wd_size0) break;
    }
    if (i < ind_wds[num_index_words[1]] && i != l_lower || i == ind_wds[num_index_words[1]] && l_lower <= u_upper) {
      fprintf(stderr, "\n%d %d %d %d %d", i, l_lower, u_upper, k, word_int);
      for (k = 0; k < ind_wd_size0; k++) {
	w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - k)) & 3];
	fprintf(stderr, " %c:%c:%c", w, seq2[k], seq_area[seq_offsets[l_lower] + k]);
      }
      fatalError("lword");
    }
    for (i = ind_wds[num_index_words[1]] - 1; i >= 0; i--) {
      seq2 = seq_area + seq_offsets[i];
      for (k = 0; k < ind_wd_size; k++) {
	w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - k)) & 3];
	w = merge_residues[w];
	if (!(w & merge_residues[seq2[k]] & 15)) break;
      }
      if (k < ind_wd_size) continue;
      for (k = ind_wd_size; k < ind_wd_size0; k++) {
	w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - k)) & 3];
	w = merge_residues[w];
	if (!(w & merge_residues[seq2[k]])) break;
      }
      if (k == ind_wd_size0) break;
    }
    if (i >= 0 && i != u_upper) {
      fprintf(stderr, "\n%d %d %d %d %d", i, l_lower, u_upper, k, word_int);
      for (k = 0; k < ind_wd_size0; k++) {
	w = "ACGT"[(word_int >> 2 * (ind_wd_size0 - 1 - k)) & 3];
	fprintf(stderr, " %c:%c:%c", w, seq2[k], seq_area[seq_offsets[l_lower] + k]);
      }
      fatalError("uword");
    }
    if (!(word_int % 1024)) notify(".");
    */
  }
  
  fprintf(stderr, "\n\nTotal, restricted word matches: %.0f  %.0f  %.3f\n", 
	  total_word_matches, restricted_word_matches, restricted_word_matches / total_word_matches);
}


translate_seq(alphabet, transalph, trans_type, sort_flag)
     char *alphabet, *transalph;
     int trans_type, sort_flag;
{
  int i;
  int trans_index_v[256];

  if (sort_flag) {
    trans_sort_v = trans_sort_vec[trans_type];
    if (parameters->DNA_flag >= 3) {
      trans_sort_v[0] = 0;
      for (i = 1; i < 256; i++) trans_sort_v[i] = i & 63 ? "ACGT"[i >> 6] : " NX*"[i >> 6]; 
    }
    else {
      for (i = 0; i < 256; i++) trans_sort_v[i] = i;
      for (i = 0; i < strlen(alphabet); i++) trans_sort_v[alphabet[i]] = transalph[i];
    }
    trans_sort_res = trans_sort_residues[trans_type];
    for (i = 0; i < 256; i++) trans_sort_res[i] = residues[trans_sort_v[i]];
  }
  else {
    /* trans_index_v = trans_index_vec[trans_type]; */
    if (parameters->DNA_flag >= 3) {
      for (i = 0; i < 256; i++) trans_index_v[i] = i & 63 ? "ACGT"[i >> 6] : i == 64 ? 'N' : 0;
    }
    else {
      for (i = 0; i < 256; i++) trans_index_v[i] = i;
      for (i = 0; i < strlen(alphabet); i++) trans_index_v[alphabet[i]] = transalph[i];
    }
    trans_index_res = trans_index_residues[trans_type];
    mod_val_trans_index_res = mod_val_trans_index_residues[trans_type];

    for (i = 0; i < 256; i++) {
      trans_index_res[i] = residues[trans_index_v[i]];
      mod_val_trans_index_res[i] = trans_index_res[i] < 0 ? trans_index_res[i] : mod_value[trans_type] * trans_index_res[i];
    }
  }
}

/*
 other input: total_length, num_entries, complexity_flag, alphabet_size, residues, minmatch, nlogn_diff, mod_value, index_shift */

set_index_words(pass, trans_type, index_type)
     int pass, trans_type, index_type;
{
  unsigned char *seqm;
  unsigned int word_int, long_word_int;
  size_t alloc_size;
  int c, w, i, k, incr, n_i_w, lw_flag, phase;
  SEQ_AREA j, n;
  double size_ind, size_min;
  char *our_alloc();
  int res_counts[256];
  static int ind_ws;
  static SEQ_AREA *ind_wds;
  static SEQ_AREA *seq_offsets; 

/* in pass 0, count no. of words; -- to allow minimum allocation;  in pass 1, populate seq_offsets */

  if (!pass) {
    ind_ws = index_word_size[trans_type];
    n_i_w = num_index_words[trans_type];
    index_words[trans_type] = ind_wds = (SEQ_AREA *)our_alloc((n_i_w + 1) * sizeof(SEQ_AREA));
    for (i = 0; i <= n_i_w; i++) ind_wds[i] = 0;
  }
  else {
    n_i_w = num_index_words[trans_type];
    for (i = 1; i < n_i_w; i++) ind_wds[i] += ind_wds[i - 1];
    ind_wds[n_i_w] = ind_wds[n_i_w - 1];

    alloc_size = ind_wds[n_i_w] + 2; /* corrected from 1 per James Bonfield's comment */
    alloc_size *= sizeof(SEQ_AREA);
    seq_offsets = (SEQ_AREA *)our_alloc(alloc_size);
    seq_offsets[0] = 0; /* sentinel -- word starting with 0, so less than every other word */
    seq_offsets++;
    s_array[trans_type] = seq_offsets;
    s_array_length[trans_type] = ind_wds[n_i_w];
  }
  /*
  fprintf(stderr, "\n sizeof(size_t), sizeof(int), n_i_w, ind_wds[n_i_w], alloc_size, total_length, num_entries: %d %d %d %d %lu %ld %d\n", 
	   (int)sizeof(size_t), (int)sizeof(int), n_i_w, ind_wds[n_i_w], alloc_size, (long int)total_length, num_entries);
  */
  /*  work backwards! */

  lw_flag = long_word_shift && !pass && !trans_type;

  seqm = seq_area + total_length + num_entries;

  size_min = 0;

  incr = pass ? -1 : 1;

 reinit:

  while (seqm >= seq_area && trans_index_res[*seqm] < 0) {
    if (trans_index_res[*seqm] < -1) size_min = 0; /* stop char for minmatch */
    else size_min++;
    seqm--;
  }
  if (seqm <= seq_area) {
    return; /* always reach this point */
  }

  long_word_int = word_int = 0;
  if (alphabet_size == 4) {
    if (complexity_flag && index_type) { /*  keep track of complexity, and only exit loop when exceed cutoff */
      for (c = 0; c < alphabet_size; c++) res_counts[c] = 0; 
      for (size_min = k = 0; size_min < x_minmatch && k < maxmatch; seqm--, k++) { /* maxmatch must be low enough that nlogn defined */
	if ((w = mod_val_trans_index_res[*seqm]) < 0) goto reinit;
	word_int = (word_int >> 2) + w; 
	if (lw_flag) 
	  long_word_int = (long_word_int >> 2) + (w << long_word_shift); 
	i = res_counts[trans_index_res[*seqm]] += 1;
	size_min += nlogn_diff[i] - nlogn_diff[k + 1]; 
      }
    }
    else {
      for (size_ind = 0; size_ind < ind_ws || size_min < minmatch /* take out preceding when need all indexed words */; seqm--, size_ind++, size_min++) {
	if ((w = mod_val_trans_index_res[*seqm]) < 0) goto reinit;
	word_int = (word_int >> 2) + w; 
	if (lw_flag) 
	  long_word_int = (long_word_int >> 2) + (w << long_word_shift); 
      }
    }

    n = ind_wds[word_int] += incr;
    if (pass) {
      seq_offsets[n] = j = (SEQ_AREA)(seqm + 1 - seq_area);
    }
    else if (lw_flag) {
      long_word_array[long_word_int >> 3] |= 1 << (long_word_int & 7);
    }

    for (phase = 0; ; seqm--) {
      if ((w = mod_val_trans_index_res[*seqm]) < 0) goto reinit;
      word_int = (word_int >> 2) + w; 
      phase++;
      if (lw_flag) {
	  long_word_int = (long_word_int >> 2) + (w << long_word_shift); 
      }
      if (phase == word_offset) {
	n = ind_wds[word_int] += incr;
	if (pass)	{
	  seq_offsets[n] = j -= phase;
	}
	else if (lw_flag) {
	  long_word_array[long_word_int >> 3] |= 1 << (long_word_int & 7);
	}
	phase = 0;
      }
    }
  }
  else {
    if (complexity_flag && index_type) { 
      for (c = 0; c < alphabet_size; c++) res_counts[c] = 0; 
      for (size_min = k = 0; size_min < x_minmatch && k < maxmatch; seqm--, k++) { 
	if ((w = mod_val_trans_index_res[*seqm]) < 0) goto reinit;
	word_int = word_int / alphabet_size + w; 
	i = res_counts[trans_index_res[*seqm]] += 1;
	size_min += nlogn_diff[i] - nlogn_diff[k + 1]; 
      }
    }
    else {
      for (size_ind = 0; size_ind < ind_ws || size_min < minmatch; seqm--, size_ind++, size_min++) {
	if ((w = mod_val_trans_index_res[*seqm]) < 0) goto reinit;
	word_int = word_int / alphabet_size + w; 
      }
    }

    n = ind_wds[word_int] += incr;
    if (pass) seq_offsets[n] = j = (SEQ_AREA)(seqm + 1 - seq_area);

    for (phase = 0; ; seqm--) {
      if ((w = mod_val_trans_index_res[*seqm]) < 0) goto reinit;
      word_int = word_int / alphabet_size + w; 
      phase++;
      if (phase == word_offset) {
	n = ind_wds[word_int] += incr;
	if (pass) seq_offsets[n] = j -= phase;
	phase = 0;
      }
    }
  }
}

/*
typedef struct comp_results {
  SEQ_AREA seq_offset;
  int d, k_end;
} Comp_results;

Comp_results *comp_results_base, *comp_results_vec;
SEQ_AREA n_comp_results;
*/

new_sort_words(trans_type, block_flag, index_type, alphabet, index_transalph, sort_transalph)
     int trans_type, block_flag, index_type;
     char *alphabet, *index_transalph, *sort_transalph;
{
  long int max_index_block;
  long int i;
  int j, n, c;
  long int print_indexword_histogram();
  SEQ_AREA *seq_offsets;
  char *our_alloc();

  notify("Indexing words ...");

  translate_seq(alphabet, index_transalph, trans_type, 0);

  set_index_words(0, trans_type, index_type);

  max_index_block = print_indexword_histogram(trans_type);

  set_index_words(1, trans_type, index_type);

  notify(" Done\n");

  notify("Sorting words ...");

  translate_seq(alphabet, sort_transalph, trans_type, 1);

  /* return; /* use if sorting/blocking not needed */
  /*
  if (n_comp_results < max_index_block) {
    if (comp_results_base) our_free(comp_results_base); 
    n_comp_results = max_index_block;
    comp_results_base = (Comp_results *)our_alloc(n_comp_results * sizeof(Comp_results));
  }
  */
  for (i = 0; i < num_index_words[trans_type]; i++) {
    if (index_words[trans_type][i] + 1 < index_words[trans_type][i + 1]) {
      /* 
     comp_results_vec = comp_results_base - index_words[trans_type][i];
      */
      quicksort(index_words[trans_type][i], index_words[trans_type][i + 1] - 1, index_word_size[trans_type], max_length /* to ensure complete sorting; was maxmatch */, s_array[trans_type]);
    }
  }
  notify(" Done\n");

  fprintf(fp_log, " Total no. words: %lu; after pruning: %lu; avg no. per indexword: %.1f\n", 
	  (unsigned long int)total_length, (unsigned long int)s_array_length[trans_type], s_array_length[trans_type] / (double)num_index_words[trans_type]);
  if (fp_debug)
    fprintf(fp_debug, " Total no. words: %lu; after pruning: %lu; avg no. per indexword: %.1f\n", 
	  (unsigned long int)total_length, (unsigned long int)s_array_length[trans_type], s_array_length[trans_type] / (double)num_index_words[trans_type]);

  /* Check whether sort is correct 
  seq_offsets = s_array[trans_type];
  notify("testing sort...");
  for (i = 1; i < s_array_length[trans_type]; i++)
    if ((n = compare_words(seq_area + seq_offsets[i], seq_area + seq_offsets[i - 1], (unsigned char *)0, &j)) < 0) {
      notify(" sort failure "); 
      fp_debug = stderr;
      fprintf(fp_debug, "\n%ld %d\n", i, n);
      fprintf(fp_debug,"%s",seq_area + seq_offsets[i - 1]);
      fprintf(fp_debug, "\n");
      fprintf(fp_debug,"%s",seq_area + seq_offsets[i]);
      fprintf(fp_debug, "\n");
      for (j = 0; ; j++) {
	c = *(seq_area + seq_offsets[i] + j);
	fprintf(fp_debug, "%d ", trans_sort_v[c]);
	if (!c) break;
      }
      fatalError("Sort");
    }
  notify(" Done\n");
  exit(1);
    /* */
  if (block_flag) 
    get_word_block(max_index_block, trans_type);

}     

/* Mark regions, and min length, of matching words, to meet maxgroupsize and complexity requirements.
   other input: num_index_words[trans_type], complexity_flag, alphabet_size, word_block, ind_wds */

get_word_block(max_index_block, trans_type)
     long int max_index_block;
     int trans_type;
{
  char *our_alloc();
  int gsize;
  int histogram[201], hist20[201];
  int group_histogram[3][201];
  SEQ_AREA *next_lower_length, *next_lower;
  SEQ_AREA next_ll[128]; 
  int left_bdry_length, w_length, k1;
  unsigned char *seq1, *seq2;
  double base_size, size;
  int i, c, k;
  long int i_lower, i_upper, j, last_jk, jk, max_group_size;
  int res_counts[256], base_res_counts[256];
  unsigned char *word_bl;
  int ind_wd_size;
  SEQ_AREA *seq_offsets;

  notify("Finding word boundaries ...");

  for (j = 0; j < 3; j++)
    for (i = 0; i <= 200; i++) 
      group_histogram[j][i] = 0;

  for (i = 0; i <= 200; i++) 
      histogram[i] = hist20[i] = 0;

  ind_wd_size = index_word_size[trans_type];
  seq_offsets = s_array[trans_type];
  word_block[trans_type] = word_bl = (unsigned char *)our_alloc(s_array_length[trans_type] * sizeof(unsigned char));
  next_lower_length = (SEQ_AREA *)our_alloc(max_index_block * sizeof(SEQ_AREA));
  max_group_size = parameters->max_group_size;

  if (max_group_size <= 0) max_group_size = total_length; /* index_words[trans_type][num_index_words[trans_type]] - 1; largest value of i_upper */

  for (i = 0; i < num_index_words[trans_type]; i++) {
    i_lower = index_words[trans_type][i];
    i_upper = index_words[trans_type][i + 1] - 1;
    if (i_upper < i_lower) continue;

/* correction 1/14/99: the following statement (i.e. use of an offset pointer) apparently causes
   problems, e.g. segmentation faults on HP machines, apparently because for
   large values of i_lower, next_lower may point to something outside the dataspace;
   this should not be a problem because it is only used with a positive displacement
   that results in pointing to a valid element, and it is not a problem on most machines.
   anyway have redefined next_lower. New version slightly less computationally efficient 
    next_lower = next_lower_length - i_lower;
*/
    next_lower = next_lower_length;

    /* Find length of match between each word and its following one */

    if (complexity_flag) {
      for (c = 0; c < alphabet_size; c++) base_res_counts[c] = 0;
      seq1 = seq_area + seq_offsets[i_lower];
      for (k = base_size = 0; k < ind_wd_size; k++) {
	c = trans_sort_res[seq1[k]];
	base_res_counts[c] += 1;
	base_size += nlogn_diff[base_res_counts[c]] - nlogn_diff[k + 1];
      }
      /*      fprintf(fp_debug, " %.1f", base_size); */
    }
    gsize = 1;
    for (j = i_lower; j < i_upper; j++) {
      seq1 = seq_area + seq_offsets[j];
      seq2 = seq_area + seq_offsets[j + 1];
      for (k = ind_wd_size; trans_sort_res[seq1[k]] == trans_sort_res[seq2[k]] && seq1[k] && k < maxmatch; k++);
      if (k >= HIST_WORD_SIZE) gsize++;
      else {
	hist20[gsize > 200 ? 200 : gsize] += 1;
	gsize = 1;
      }
      word_bl[j] = k; /* no. matching chars, or index of 1st mismatch */
    }
    hist20[gsize > 200 ? 200 : gsize] += 1;
    word_bl[i_upper] = minmatch - 1;
    /* Now work back, finding distance to next occurrence of a smaller matching word size */
    for (k = 0; k <= maxmatch; k++)
      next_ll[k] = i_upper;
    next_lower[i_upper - i_lower] = i_upper;
    for (j = i_upper - 1; j >= i_lower; j--) {
      next_lower[j - i_lower] = next_ll[word_bl[j] - 1];
      /* Could probably improve efficiency of following loop -- which is potential time hog */
      for (k = word_bl[j]; k <= maxmatch; k++)
	next_ll[k] = j;
    }
    left_bdry_length = minmatch - 1; /* ind_wd_size - 1; */
    for (j = i_lower; j <= i_upper; j++) {
      if (complexity_flag) {
	seq1 = seq_area + seq_offsets[j];
	for (c = 0; c < alphabet_size; c++) 
	  res_counts[c] = base_res_counts[c];
	for (size = base_size, k = ind_wd_size; size < x_minmatch && k < maxmatch; k++) {
	  c = trans_sort_res[seq1[k]];
	  if (c < 0) break;
	  res_counts[c] += 1;
	  size += nlogn_diff[res_counts[c]] - nlogn_diff[k + 1];
 	}
	if (left_bdry_length < k - 1) left_bdry_length = k - 1;
	/*	fprintf(fp_debug, "\n%.1f %d %d", size, k, c); */
      }
      for (last_jk = j, jk = next_lower[j - i_lower]; 
	   word_bl[last_jk] > left_bdry_length
	   && (jk - j < max_group_size || word_bl[last_jk] >= maxmatch); 
	   last_jk = jk, jk = next_lower[jk - i_lower]);
      
      w_length = word_bl[last_jk] > left_bdry_length ? word_bl[last_jk] : left_bdry_length;
      w_length++;
      group_histogram[last_jk > j + 1 ? 2 : (last_jk > j ? 1 : 0)][w_length] += 1;

      histogram[last_jk >= j + 200 ? 200 : last_jk - j + 1] += 1;
      left_bdry_length = minmatch - 1 > word_bl[last_jk] ? 
	minmatch - 1 : word_bl[last_jk] ;
      
      for (jk = j; jk <= last_jk; jk++) {
	word_bl[jk] = w_length;
	/*
	for (k1 = 0; k1 < w_length; k1++)
	  if (!seq_area[seq_offsets[k] + k1]) fatalError("word in 0");
	*/
      }
      word_bl[j] |= 128; /* left bdry flag */
      j = last_jk;
    }
  }
  notify(" Done\n");

  if (fp_debug) {
    fprintf(fp_debug, "\nHistogram of no. of %d-mers with given no. of occurrences", HIST_WORD_SIZE);
    for (i = 0; i <= 200; i++) 
      if (hist20[i]) 
	fprintf(fp_debug, "\n%3d  %6d", i, hist20[i]);

    fprintf(fp_debug, "\n\nHistogram of no. of word groups of given length (group sizes 1, 2, and > 2)\n");
    for (i = 0; i <= 200; i++) {
      if (group_histogram[0][i] || group_histogram[1][i] || group_histogram[2][i]) 
	fprintf(fp_debug, "%3d  %6d  %6d  %6d\n", 
		i, group_histogram[0][i], group_histogram[1][i], group_histogram[2][i]);
    }

    fprintf(fp_debug, "\n\nHistogram of no. of word groups of given size\n");
    for (i = 0; i <= 200; i++) {
      if (histogram[i]) 
	fprintf(fp_debug, "%3d  %6d\n", i, histogram[i]);
    }
  }
  our_free(next_lower_length);
}

long int print_indexword_histogram(trans_type)
     int trans_type;
{
  long int histogram[201];
  int i, i_min, i_max;
  long int min, max, max_index_block, d;

  max = 0;
  min = 100000000;
  max_index_block = 0;
  for (i = 0; i <= 200; i++) histogram[i] = 0;
  for (i = 0; i < num_index_words[trans_type]; i++) {
    d = index_words[trans_type][i];
    if (max_index_block < d) max_index_block = d;

    if (min > d) {min = d; i_min = i; }
    if (max < d) {max = d; i_max = i; }
    d = (d + 9) / 10;
    histogram[d > 100 ? 100 : d] += 1;
  }
  if (fp_debug) {
    for (i = 0; i <= 200; i++) {
      if (histogram[i]) fprintf(fp_debug, "\n %d  %ld", i * 10, histogram[i]);
    }
    fprintf(fp_debug, "\nmin no. words: %ld at %d", min, i_min);
    fprintf(fp_debug, "\nmax no. words: %ld at %d", max, i_max);
  }
  return max_index_block;
}

/* following no longer used
compare_positions(position1, position2)
     SEQ_AREA *position1, *position2;
{
  unsigned char *seq1, *seq2;
  int k;

  seq1 = seq_area + *position1;
  seq2 = seq_area + *position2;
  for (k = 0; k < minmatch; k++, seq1++, seq2++)
    if (*seq1 != *seq2) return *seq1 - *seq2;
  return 0;
}
*/

compare_words(seq1, seq2, w_stop, k_end)
     unsigned char *seq1, *seq2, *w_stop;
     int *k_end;
{
  unsigned char *w_start;
  int c, d;
/*
  n_compares++;
*/
  w_start = seq1;
  for (; /* seq1 < w_stop */; seq1++, seq2++) { /* test inactivated -- to ensure full sort & deal with addressing error 5/5/08 */
    c = trans_sort_v[*seq1]; /* was trans_sort_res */
    d = trans_sort_v[*seq2];
    if (c != d || !c) {
      *k_end = seq1 - w_start;
/*
      t_length += *k_end + 1;
*/
      return c - d;
    }
  }
}

/*
compare_comp_results(comp_ptr1, comp_ptr2)
     Comp_results *comp_ptr1, *comp_ptr2;
{
  if (comp_ptr1->d >= 0 && comp_ptr2->d <= 0
      || comp_ptr1->d <= 0 && comp_ptr2->d >= 0 
      || comp_ptr1->k_end == comp_ptr2->k_end) /* at this point know both d's positive, or both negative */ /*
    return comp_ptr1->d - comp_ptr2->d;

  return comp_ptr1->d > 0 ? comp_ptr2->k_end - comp_ptr1->k_end : comp_ptr1->k_end - comp_ptr2->k_end;
}
*/

/* new sorting strategy that minimizes sequence accesses -- by taking advantage of fact that sequence comparison yields
useful information that was being largely ignored
-- however it is slower!!! */

/*
newsort(left, right, k_start, maxcheck, seq_offsets)
     SEQ_AREA left, right;
     int k_start, maxcheck;
     SEQ_AREA *seq_offsets;
{
  unsigned char *w_right, *s_plus_k;
  SEQ_AREA i, j, temp;
  int k_end, start_i, start_d, start_k_end;
  Comp_results *comp_results, *base_comp_res;

  s_plus_k = seq_area + k_start;
  if (right < left + 4) {
    for (j = left + 1; j <= right; j++) {
      temp = seq_offsets[j];
      w_right = s_plus_k + temp;
      for (i = j; i > left &&
	   compare_words(w_right, s_plus_k + seq_offsets[i - 1], (unsigned char *)0, &k_end) < 0; i--) 
	seq_offsets[i] = seq_offsets[i - 1];
      seq_offsets[i] = temp;
    }
    return;
  }
      
  w_right = s_plus_k + seq_offsets[right];

  base_comp_res = comp_results_vec + left;
  for (i = left, comp_results = base_comp_res; i < right; i++, comp_results++) { 
    comp_results->seq_offset = seq_offsets[i];
    comp_results->d = compare_words(s_plus_k + seq_offsets[i], w_right, (unsigned char *)0, &k_end);
    comp_results->k_end = comp_results->d ? k_end : 0;
  }
  comp_results->d = comp_results->k_end = 0; /* w_right to itself; don't need k_end */ /*
  comp_results->seq_offset = seq_offsets[i];

  qsort(base_comp_res, (size_t)(right - left + 1), sizeof(Comp_results), compare_comp_results);

  for (i = left, comp_results = base_comp_res; i <= right; i++, comp_results++) { 
    seq_offsets[i] = comp_results->seq_offset;
  }
  comp_results = base_comp_res;
  start_i = left;
  start_d = comp_results->d;
  start_k_end = comp_results->k_end;
  for (i = left + 1, comp_results++; i <= right + 1; i++, comp_results++) { 
    if (i == right + 1 || comp_results->d != start_d || comp_results->k_end != start_k_end) {
      if (start_d && i - 1 > start_i)
	newsort(start_i, i - 1, k_start + start_k_end, maxcheck, seq_offsets);
      start_i = i;
      start_d = comp_results->d;
      start_k_end = comp_results->k_end;
    }
  }
}
*/

/* quicksort (based in part on generic recursive version with insertion sort
   refinement in Sedgwick, but using bookkeeping of matching leading parts 
   of words in a block to substantially improve
   efficiency) */

quicksort(left, right, k_start, maxcheck, seq_offsets)
     SEQ_AREA left, right;
     int k_start, maxcheck;
     SEQ_AREA *seq_offsets;
{
  unsigned char *w_right, *w_stop, *s_plus_k;
  SEQ_AREA i, j, temp;
  int k_end, k_stop, k_left_start, k_right_start;
  int i_lower, i_upper, c;

/* use insertion sort for small pieces */
  k_stop = maxcheck - k_start;
  s_plus_k = seq_area + k_start;
  if (right < left + 4) {
    for (j = left + 1; j <= right; j++) {
      temp = seq_offsets[j];
      w_right = s_plus_k + temp;
      w_stop = 0; /* w_right + k_stop; test inactivated -- to ensure full sort & deal with addressing error 5/5/08 */
      /* following involves fewer comparisons than standard insertion sort -- but less efficient !! 
      i_lower = left;
      i_upper = j;
      do {
	i = (i_lower + i_upper) >> 1;
	c = compare_words(w_right, s_plus_k + seq_offsets[i], w_stop, &k_end);
	if (c < 0) i_upper = i;
	else i_lower = i + 1;
      } while (i_lower < i_upper);
      for (i = j; i > i_upper; i--)
	seq_offsets[i] = seq_offsets[i - 1];
      /* */
      /* standard insertion sort */
      for (i = j; i > left &&
	   compare_words(w_right, s_plus_k + seq_offsets[i - 1], w_stop, &k_end) < 0; i--) 
	seq_offsets[i] = seq_offsets[i - 1];
	/* */
      seq_offsets[i] = temp;
    }
    return;
  }
      
  w_right = s_plus_k + seq_offsets[right];
  w_stop = 0; /* w_right + k_stop; test inactivated -- to ensure full sort & deal with addressing error 5/5/08 */ 
  i = left; /* left - 1; */
  j = right;
  k_left_start = k_right_start = k_stop;
  for (;; i++) {
    for (; i < j; i++) { 
      if (compare_words(w_right, s_plus_k + seq_offsets[i], w_stop, &k_end) <= 0) {
	if (k_end < k_right_start) k_right_start = k_end;
	break;
      }
      else if (k_end < k_left_start) k_left_start = k_end;
    }
    for (j--; j > i; j--) {
      if (compare_words(w_right, s_plus_k + seq_offsets[j], w_stop, &k_end) >= 0) {
	if (k_end < k_left_start) k_left_start = k_end;
	break;
      }
      else if (k_end < k_right_start) k_right_start = k_end;
    }
    if (j <= i) break;
    temp = seq_offsets[i];
    seq_offsets[i] = seq_offsets[j];
    seq_offsets[j] = temp;
  } 
  if (j == i && j > left) j--;
  temp = seq_offsets[i];
  seq_offsets[i] = seq_offsets[right];
  seq_offsets[right] = temp;
  /* leading k_left_start characters match for all words in left part of this block,
     similarly for k_right_start */
  k_left_start += k_start; /* += k_start; */
  if (k_left_start < maxcheck && left < j) 
    quicksort(left, j, k_left_start, maxcheck, seq_offsets);
  k_right_start += k_start; /* += k_start; */
  if (k_right_start < maxcheck && i + 1 < right) 
    quicksort(i + 1, right, k_right_start, maxcheck, seq_offsets);
}

find_subject_matches(files)
     File *files;
{
  File *file;
  int entry1, i;
  Aligned_pair *get_aligned_pairs();
  unsigned char *seq;
  Database *sdb;
  Segment *seg_list;
  Segment *repeat_parse_descrip();

  notify("Finding subject matches ...");

  init_word_counts();

  for (file = files; file; file = file->next) {
    sdb = file->db;
/* Note -- this assumes file is not disturbed!! */
    while (get_next_file_entry(sdb)) {
      entry1 = append_seq_entry(sdb);
      add_subject_length(get_seq_length(entry1));

      /*
      seg_list = s_repeat_screen ? repeat_parse_descrip(sdb->descrip_buffer) : (Segment *)0;
      */
      /* notify("Finding external word matches..."); */
      find_external_word_matches(entry1, /* seg_list, */ sdb->seq_buffer);
      /* notify(" Done\n"); */
/* in this case, will not be doing comp or internal matches; need to generate list
   of cand. pairs, and do swat alignments immediately */

      find_scores(entry1, sdb->seq_buffer, 1); /* have stored cand_pairs under entry2 rather than entry1 */

      if (!get_aligned_pairs(entry1)) remove_seq_entry(sdb);
      free_cand_pair_blocks();
      free_seg_blocks();
    }
    free_seq_buffers(sdb);
  }
  notify(" Done\n");
  print_word_counts();
}

find_comp_matches()
{
  int entry1;
  unsigned char *seq;
  unsigned char *get_comp_seq();

  notify("Finding complement word matches ...");

  init_word_counts();

  for (entry1 = start_entry; entry1 < start_entry + num_uncomp_entries; entry1++) {
    seq = get_comp_seq(entry1);
    find_external_word_matches(entry1, /* (Segment *)0, */ seq);

/* in this case, will not be doing subject matches; create candidate pairs to process
   in the usual fashion */

  }
  notify(" Done\n");
  print_word_counts();
}

new_find_internal_word_matches()
{
  unsigned char *start_gp, *start2, *start3;
  int offset1, offset2, entry1, entry2, msize;
  int k, m, c;
  int vector_bound, complement, vec_on;
  SEQ_AREA *seq_offsets;
  SEQ_AREA *ind_wds;
  SEQ_AREA i, j, i_upper, group_start;
  unsigned char *word_bl;

  notify("Finding internal word matches ...");

  init_word_counts();
  seq_offsets = s_array[0];
  ind_wds = index_words[0];
  word_bl = word_block[0];

  vector_bound = parameters->vector_bound;
  for (m = 0; m < num_index_words[0]; m++) { 

    i_upper = ind_wds[m + 1];
    group_start = ind_wds[m];

    start_gp = seq_area + seq_offsets[group_start];
    for (i = group_start + 1; i < i_upper; i++) {
      start2 = seq_area + seq_offsets[i];
      if (!(word_bl[i] & 128)) {
	entry2 = -1;
	vec_on = 0;
	c = trans_sort_res[start2[-1]]; /* PROBABLY SHOULD BE TRANS_INDEX_RES */
	for (j = group_start; j < i; j++) {

	  start3 = seq_area + seq_offsets[j];
	  if (c < 0 || trans_sort_res[start3[-1]] != c) {
	    /* only consider full-length matches */
	    if (entry2 < 0) {
	      find_entry_and_offset(seq_offsets[i], &entry2, &offset2, &complement);
	      if (offset2 + minmatch < vector_bound)
		vec_on = 1;
	    }
	    find_entry_and_offset(seq_offsets[j], &entry1, &offset1, &complement);
/* ignore nucleotide word matches entirely contained within first 
   parameters->vector_bound residues of both reads  */
	    n_same_matches++;
	    msize = 0;
	    if (vec_on) {
	      if (!msize) msize = get_matchsize_nocomplexity(start2, start3, 0); 
	      /* note asymmetry in condition -- should this be changed? */
	      if (offset1 + msize < vector_bound) {
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
	    }
	    /* FOLLOWING SHOULD ALREADY BE DEALT WITH BY WORD_BLOCK 
	    if (complexity_flag && get_corrected_matchsize() < minmatch) {
	      num_le_fail_words++;
	      continue;
	    }
	    */
	    if (q_repeat_screen) {
	      if (!msize) msize = get_matchsize_nocomplexity(start2, start3, 0); 
	      if (contained_in_tag(entry1, "repeat", offset1 + 1, offset1 + msize, 0, 0)
		 || contained_in_tag(entry2, "repeat", offset2 + 1, offset2 + msize, 0, 0)) {
		num_q_repeat_fail_words++;
		continue;
	      }
	    }

	    num_pass_words++;
	    if (entry1 > entry2) 
	      pair_routine(entry1, entry2, offset1, offset2, 0, 0);
	    else 
	      pair_routine(entry2, entry1, offset2, offset1, 0, 0);
	  }
	} 
      }
      else {
	group_start = i;
	start_gp = start2;
      }
    }
  }
  print_word_counts();
  notify(" Done\n");
}

find_entry_and_offset(position, add_entry, add_offset, add_complement)
     SEQ_AREA position;
     int *add_entry, *add_offset, *add_complement;
{
  unsigned int upper, lower, mid;
  double n_entries;
  SEQ_AREA n, upper_pos, lower_pos;

  upper = num_entries;
  lower = 0;
  upper_pos = seq_area_size;
  lower_pos = 1;

  do {
    n_entries = upper - lower;
    mid = lower + (n_entries * (position - lower_pos)) / (upper_pos - lower_pos);
    if (position_table[mid] <= position) {
      if ((n = position_table[mid + 1]) > position) {
	lower = mid;
	break;
      }
      else {
	lower = mid + 1;
	lower_pos = n;
      }
    }
    else {
      if ((n = position_table[mid - 1]) <= position) {
	lower = mid - 1;
	break;
      }
      else {
	upper = mid - 1;
	upper_pos = n;
      }
    }
    /* following by itself is old version */
    mid = (upper + lower) >> 1;
    if ((n = position_table[mid]) > position) {
      upper = mid;
      upper_pos = n;
    }
    else {
      lower = mid;
      lower_pos = n;
    }
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
  notify("Freeing word arrays ...");
  our_free(s_array[0] - 1); /* N.B. -- s_array no longer available! */
  if (s_array[1]) our_free(s_array[1] - 1); /* N.B. -- s_array no longer available! */
  our_free(position_table);
  our_free(index_words[0]);
  if (index_words[1]) our_free(index_words[1]);
  our_free(word_block[0]);
  if (word_block[1])our_free(word_block[1]);
  if (long_word_shift) our_free(long_word_array);
  notify(" Done\n");
}

static double p;

make_nlogn()
{
  int n;

  p = -log((double)alphabet_size);
  nlogn[0] = nlogn_diff[0] = 0;
  for (n = 1; n < 500; n++) {
    nlogn[n] = n * log((double)n) / p;
    nlogn_diff[n] = nlogn[n] - nlogn[n - 1];
  }
}

int get_matchsize_complexity(seq1, seq2, minsize)
     unsigned char *seq1, *seq2;
     int minsize;
{
  int i, r1, r2, r;
  unsigned char *init_seq1, *init_seq2;

  seq1 += minsize;
  seq2 += minsize;

  if (minsize) 
    for (i = 0; i < alphabet_size; i++) res_counts[i] = base_res_counts[i];
  else
    for (i = 0; i < alphabet_size; i++) res_counts[i] = 0;

  r = minsize;

  if (merged_flag) {
    for (;;) {
      r1 = merge_residues[*seq1]; /* residues[*seq1]; */  /* NOT CORRECT FOR PROTEIN */
      r2 = merge_residues[*seq2]; /* trans_sort_res[*seq2]; */
      if (r1 & r2 && r1 != 255) {
	res_counts[residues[*seq1]] += 1; 
	seq1++;
	seq2++;
	r++;
      }
      else break;
    }
  }
  else {
    for (;;) {
      r1 = residues[*seq1]; 
      r2 = trans_index_res[*seq2]; 
      if (r1 == r2 && r1 >= 0) {
	res_counts[r1] += 1; 
	seq1++;
	seq2++;
	r++;
      }
      else break;
    }
  }
  /*  if (r < minmatch) fatalError("matchsize < minmatch"); */
  return r;
}

int get_matchsize_nocomplexity(seq1, seq2, minsize)
     unsigned char *seq1, *seq2;
     int minsize;
{
  int i, r1, r2, r;
  unsigned char *init_seq1, *init_seq2;

  seq1 += minsize;
  seq2 += minsize;

  r = minsize;

  if (merged_flag) {
    for (;;) {
      r1 = merge_residues[*seq1]; /* residues[*seq1]; */  /* NOT CORRECT FOR PROTEIN */
      r2 = merge_residues[*seq2]; /* trans_sort_res[*seq2]; */
      if ((r1 & r2) && r1 != 255) { /* DON'T REWARD 'N' */
	seq1++;
	seq2++;
	r++; 
      }
      else break;
    }
  }
  else {
    for (;;) {
      r1 = residues[*seq1]; 
      r2 = residues[*seq2]; 
      if (r1 == r2 && r1 >= 0) { /* DON'T REWARD 'N' */
	seq1++;
	seq2++;
	r++;
      }
      else break;
    }
  }
  /*  if (r < minmatch) fatalError("matchsize < minmatch"); */
  return r;
}

int get_corrected_matchsize()
{
  double size, x;
  int i, c, n;


  /* fprintf(stderr, "\nLow entropy word %.8f:", p); */
  for (i = size = n = 0; i < alphabet_size; i++) {
    c = res_counts[i];
    n += c;
    x = c < 500 ? nlogn[c] : c * log((double)c) / p; 
    size += x;
    /* fprintf(stderr, " %d", c); */
  }
  size -= n < 500 ? nlogn[n] : n * log((double)n) / p; 
  

  return (int)(size + .999); /* round up to next highest integer */
}

double return_corrected_matchsize(counts)
     int counts[];
{
  double size, x;
  int i, c, n;

  size = n = 0;
  for (i = 0; i < alphabet_size; i++) {
    c = counts[i];
    n += c;
    x = c < 500 ? nlogn[c] : c * log((double)c) / p; 
    size += x;
  }
  size -= n < 500 ? nlogn[n] : n * log((double)n) / p; 
  return size;
}

static int entry1, strand; /* no. consec lower case letters from current point; for subject
			    files can only be positive if s_repeat_screen is on
			 */
			    
/* static Segment *seg_list; */
static unsigned char *base_seq, *word_bl;
/* static int word_stack[3]; 
static int last_ss, last_off1;

*/ 
 
find_external_word_matches(e1, seq) 
     int e1;
     unsigned char *seq;
{
  unsigned char *seqm;
  int length, i, k, w, c, mod_val, l_lower, u_upper, pos, n_lc, ind_wd_size2, found; 
  Ssite *ssite, *ssite2;
  unsigned int word_int, long_word_int, w2, w2_long, w5, w5_long;
  double size;
  SEQ_AREA *ind_wds;


  for (strand = 0; strand < 2; strand++)
    first_ssite[strand] = last_ssite[strand] = 0;

  i_off1 = 0;

  word_bl = word_block[0];
  trans_sort_v = trans_sort_vec[0];
  ind_wds = index_words[0];
  mod_val = mod_value[0];
  mod_val_res = mod_val_residues[0];
  ind_wd_size2 = index_word_size[0] - 2;

  entry1 = e1;
  base_seq = seq;

  length = get_seq_length(entry1);
  seqm = seq + length - 1;

  n_lc = 0;

 reinit:

  /* using sentinel -- first char in buffer is 0, so residues < 0 */

  while (seqm >= seq && residues[*seqm] < 0) seqm--;
  if (seqm < seq) return; /* sentinel */

  long_word_int = word_int = 0;
  /* word_stack[0] = word_stack[1] = word_stack[2] = 0; */

  /* FOLLOWING MAY BE INEFFICIENT -- IF LOTS OF SHORT SUBJECT SEQUENCES (E.G. READ COMPLEMENTS) */
  if (complexity_flag) { 
    for (c = 0; c < alphabet_size; c++) res_counts[c] = 0;
    for (size = k = 0; size < x_minmatch && k < maxmatch; seqm--, k++) { /* maxmatch must be low enough that nlogn defined */
      if ((w = residues[*seqm]) < 0) goto reinit;
      word_int = word_int / alphabet_size + w * mod_val; /* could use precomputed vals as below, but then w not appropr for res_counts? */
      if (long_word_shift)
	long_word_int = (long_word_int >> 2) + ((w * mod_val) << long_word_shift); 
      /*
      word_stack[2] = word_stack[1];
      word_stack[1] = word_stack[0];
      word_stack[0] = word_int;
      */
      n_lc = islower(*seqm) ? n_lc + 1 : 0;
      res_counts[w] += 1; 
      size += nlogn_diff[res_counts[w]] - nlogn_diff[k + 1];
    }
  }
  else {
    for (i = 0; i < minmatch; seqm--, i++) {
      if ((w = mod_val_res[*seqm]/*residues[*seqm]*/) < 0) goto reinit;
      word_int = word_int / alphabet_size + w; /* w * mod_val; */
      if (long_word_shift)
	long_word_int = (long_word_int >> 2) + (w << long_word_shift); 
      /*
      word_stack[2] = word_stack[1];
      word_stack[1] = word_stack[0];
      word_stack[0] = word_int;
      */
      n_lc = islower(*seqm) ? n_lc + 1 : 0;
    }
  }
  
  if (!long_word_shift || long_word_array[long_word_int >> 3] & (1 << (long_word_int & 7)))
    new_lookup_words(seqm + 1, 0, ind_wds[word_int], ind_wds[word_int + 1], 0, n_lc);

  if (merged_flag) {
    new_lookup_words(seqm + 1, 1, index_lower[word_int], index_upper[word_int], 0, n_lc);
  }
    
  if (alphabet_size == 4) {
    for (; ; seqm--) {
      if ((w = mod_val_res[*seqm]/*residues[*seqm]*/) < 0) goto reinit;
      /* if (!((seqm - seq) % 10000)) fprintf(stderr, " %d", seqm - seq); */
      n_lc = islower(*seqm) ? n_lc + 1 : 0;

      word_int = (word_int >> 2) + w; 
      if (long_word_shift) {
	long_word_int = (long_word_int >> 2) + (w << long_word_shift); 
	if (long_word_array[long_word_int >> 3] & (1 << (long_word_int & 7)))
	  new_lookup_words(seqm, 0, ind_wds[word_int], ind_wds[word_int + 1], 0, n_lc); 
      }
      /*
      word_stack[2] = word_stack[1];
      word_stack[1] = word_stack[0];
      word_stack[0] = word_int;
      */
      /*      if (word_int > num_index_words[0] - 1) fatalError("word out of range"); */

      else
	new_lookup_words(seqm, 0, ind_wds[word_int], ind_wds[word_int + 1], 0, n_lc); 

      if (merged_flag) {
	new_lookup_words(seqm, 1, index_lower[word_int], index_upper[word_int], 0, n_lc);
      }

      if (max_intron || max_intron2) {  
	w2 = word_int & 15;
	strand = w2 == 11 || w2 == 9 ? 0 : w2 == 7 ? 1 : 2;
	if (strand < 2 && n_lc < ind_wd_size2) { /* not quite right -- will allow some sites at intron edge */
	  pos = (seqm - seq) + indwd2; /* 1st base of intron */
	  w5 = ((word_int >> 4) & left_ss_wdmask) << left_ss_wdshift;
	  w5_long = w5 << long_word_shift;
	  if (strand >= 2 || strand < 0 || first_ssite[strand] < 0 || last_ssite[strand] > n_ssites || last_ssite[strand] < first_ssite[strand])
	    fatalError("ssite range");
	  for (i = first_ssite[strand], found = 0; i < last_ssite[strand]; i++) {
	    ssite = ssites[strand] + i;
	    if (ssite->pos - pos > ssite->keep_length) {
	      /*
	      if (i != first_ssite[strand]) {
		fprintf(stderr, "\n%d %d", i, first_ssite[strand]);
		fatalError("register shift");
	      }
	      first_ssite[strand] += 1; 
	      */
	      continue;
	    }
	    if (ssite->pos < pos + min_intron) break;
	    if (!found) {
	      found = 1;
	      first_ssite[strand] = i;
	    }
	    w2 = w5 + ssite->word_int;

	    if (long_word_shift) {
	      w2_long = w5_long + ssite->long_word_int;
	      if (long_word_array[w2_long >> 3] & (1 << (w2_long & 7)))
		new_lookup_words(seq + ssite->pos - left_ss_wdsize, 0, ind_wds[w2], ind_wds[w2 + 1], pos - ssite->pos, 0);
	    }
	    else 
	      new_lookup_words(seq + ssite->pos - left_ss_wdsize, 0, ind_wds[w2], ind_wds[w2 + 1], pos - ssite->pos, 0);

	    /*
	    printf("\n%d %d ", w2, strand);
	    for (k = 0; k < 12; k++) {
	      printf("%c", "acgt"[w2 % 4]);
	      w2 /= 4;
	    }
	    printf(" ");
	    for (k = pos - 7; k < ssite->pos + 5; k++)
	      printf("%c", seq[k]);
	    */

	    n_splice_calls++;
	  }
	}
	w2 = word_int >> indwd22;
	strand = w2 == 2 ? 0 : w2 == 1 || w2 == 9 ? 1 : 2;
	if (strand < 2 && n_lc < 2) {
	  ssite = ssites[strand] + last_ssite[strand];
	  pos = ssite->pos = (seqm - seq) + 2; /* 1st base of exon */
	  ssite->word_int = (word_int >> right_ss_wdshift) & right_ss_wdmask;
	  ssite->long_word_int = (word_int >> right_ss_long_wdshift) & right_ss_long_wdmask;
	  ssite->keep_length = max_intron;
	  last_ssite[strand] += 1;

	  if (max_intron2 > max_intron) {
	    for (i = i_off1 - 1; i >= 0 && pos >= off1s[i] - spliced_match_left; i--)
	      if (pos <= off1s[i] + spliced_match_right) {
		ssite->keep_length = max_intron2;
		break;
	      }
	  }
	  /*	  
	  last_ss = pos;
	  if (last_off1 < last_ss && last_off1)
	    fprintf(stderr, " %d", last_ss - last_off1);
	  */

	  for (i = first_ssite[strand]; i < last_ssite[strand]; i++) {
	    ssite = ssites[strand] + i;
	    if (ssite->pos - pos > ssite->keep_length) {
	      first_ssite[strand] += 1;
	    }
	    else break;
	  }

	  if (last_ssite[strand] >= n_ssites) {
	    for (i = first_ssite[strand], ssite2 = ssites[strand]; i < last_ssite[strand]; i++) {
	      ssite = ssites[strand] + i;
	      if (ssite->pos - pos > ssite->keep_length) continue;
	      ssite2->pos = ssite->pos;
	      ssite2->word_int = ssite->word_int;
	      ssite2->long_word_int = ssite->long_word_int;
	      ssite2->keep_length = ssite->keep_length;
	      ssite2++;
	    }
	    first_ssite[strand] = 0;
	    last_ssite[strand] = ssite2 - ssites[strand];
	    if (last_ssite[strand] >= n_ssites || last_ssite[strand] < 0) {
	      fprintf(stderr, "\n%d %d %d", strand, last_ssite[strand], n_ssites);
	      for (i = first_ssite[strand]; i < last_ssite[strand]; i++) {
		ssite = ssites[strand] + i;
		fprintf(stderr, "[%d %d]", ssite->pos, ssite->word_int);
	      }
		
	      fatalError("max no. ssites exceeded");
	    }
	  }
	}
	/*
	if (seqm - seq < length - 100) {
	  printf("\n%s", seqm);
	  for (strand = 0; strand < 2; strand++) {
	    printf("\n%d %d %d ", strand, first_ssite[strand], last_ssite[strand]);
	    for (i = 0; i < last_ssite[strand]; i++) {
	      ssite = ssites[strand] + i;
	      printf(" %d %d", ssite->pos - length, ssite->word_int);
	    }
	  }
	  exit(1);
	}
	*/
      }
    }
  }
  else {
    for (; ; seqm--) {
      if ((w = mod_val_res[*seqm]/*residues[*seqm]*/) < 0) goto reinit;
      word_int = word_int / alphabet_size + w; /* w * mod_val; */
      n_lc = islower(*seqm) ? n_lc + 1 : 0;
      new_lookup_words(seqm, 0, ind_wds[word_int], ind_wds[word_int + 1], 0, n_lc); 
      if (merged_flag) {
	new_lookup_words(seqm, 1, index_lower[word_int], index_upper[word_int], 0, n_lc);
      }
    }
  }
}

int known_match_size, target_match_size, offset1, sub_check, reuse;
extern int num_pairs;

new_lookup_words(seq1, trans_type, sa_l_lower, sa_u_upper, ss_displace, n_lc)
     int trans_type, ss_displace, n_lc;
     SEQ_AREA sa_l_lower, sa_u_upper;
     unsigned char *seq1;
{
  int k, k_end, b, c, d, min_k, u_k, l_k, r, s, u, v, k_target;
  int incr;
  int ind_wd_size, min_target;
  int last_n;
  unsigned char *seq2;
  double base_size, size;
  double return_corrected_matchsize();
  SEQ_AREA *seq_offsets;
  /* SEQ_AREA stack_l_lower, stack_u_upper; */
  int extend_match();
  long int l_lower, u_upper;
  SEQ_AREA mid, i;
  int c_test, test_left;

  if (sa_u_upper <= sa_l_lower) return;
  l_lower = sa_l_lower;
  u_upper = sa_u_upper - 1;

  seq_offsets = s_array[trans_type];
  ind_wd_size = index_word_size[trans_type];
  trans_index_res = trans_index_residues[trans_type];

  /* return; /* TESTING */

  if (!trans_type) {

    /* approach useful when indexing preserves order (& saves need to sort -- but not much time
       NEED TYPE CHANGES TO REVIVE THIS -- ALSO SET TEST_LEFT PARAMETER
    stack_l_lower = index_words[0][word_stack[2]];
    stack_u_upper = index_words[0][word_stack[2] + 1] - 1;
    last_n = 0;
    offset1 = seq1 - base_seq;
    c = residues[seq1[-1]];
    target_match_size = known_match_size = 14; 
    sub_check = known_match_size < target_match_size;
    for (i = l_lower, j = stack_l_lower; i <= u_upper && j <= stack_u_upper; i++) {
      n = seq_offsets[i] + 2;
      if (n <= last_n) fatalError("word order");
      last_n = n;
      for (; j <= stack_u_upper && seq_offsets[j] < n; j++);
      if (j <= stack_u_upper && seq_offsets[j] == n) {
	j++;
	seq2 = seq_area + seq_offsets[i];
	if (!ss_displace && c == trans_index_res[seq2[-1]] && c >= 0) continue; 
	if (4 == extend_match(seq1, seq2, ss_displace, test_left)) break; 
      }
    }
    return;
    */

    min_target = minmatch >= n_lc + 1 ? minmatch : n_lc + 1;

    incr = toupper(seq1[ind_wd_size]) == 'T'; /* */
    /* u_k = l_k = min_k = ind_wd_size; */
    for (;;) {
      mid = (u_upper + l_lower + incr ) >> 1; /* / 2 */
      /* following three lines && test vs k_target are to deal with long perfect match cases */
      k_target = word_bl[mid] & 127;
      if (k_target < min_target) k_target = min_target; /* is this necessary?? */

      /*
	if (mid > s_array_length[0]) {
	fprintf(stderr, "\n%ld %ld %ld\n", (long int)u_upper, (long int)l_lower, (long int)s_array_length[0]);
	fatalError("word range");
	}
      */
      seq2 = seq_area + seq_offsets[mid];

      /*
	for (k = 0; k <= target_match_size; k++)
	if (seq2[k] < 0) fatalError("mismatch");
      */
      for (k = ind_wd_size;
	   (u = toupper(seq1[k])) == (v = trans_sort_v[seq2[k]]) && u && k < k_target; k++); /* HAVE BETTER CHECK FOR SEQ END -- different chars */

      if (k < k_target /* k < min_target || k < (word_bl[mid] & 127) */) { /* 1st cond'n is to avoid word_bl access where unnecessary */

	if (u < v) { /* need to make sure this is how it's sorted!!! */ /* WHAT IF u = v = 0?? */
	  /* following saves some time, but relatively little -- & could be inefficient in some cases /* */
	  for (u_upper = mid; !(word_bl[u_upper] & 128); u_upper--);
	  u_upper--;
	  /* */
	  /* u_upper = mid - 1; /* so could go negative!! -- hence can't be unsigned */
	  /* u_k  = k; */
	  incr = 0;
	}
	else {
	  /* */
	     for (l_lower = mid + 1; l_lower <= u_upper && !(word_bl[l_lower] & 128); l_lower++); 
	     /* */
	     /* l_lower = mid + 1; /* */
	  /* l_k = k; */
	     incr = 1;
	}
	if (u_upper < l_lower) return; /* not found */
      }
      else break; /* found */
      /* min_k = u_k < l_k ? u_k : l_k; */
    }

    /* return; /* TESTING */
    target_match_size = k_target;
    /*
    target_match_size = word_bl[mid] & 127;
    if (target_match_size < min_target) target_match_size = min_target;
    */
    /*
      if (base_size < x_minmatch) {
      fprintf(fp_debug, "\nLow complexity: %.1f  %d  %.1f", base_size, target_match_size, return_corrected_matchsize(base_res_counts));
      }
    */
    /* following is inefficient; instead should mimic procedure in old
       lookup_words */

    offset1 = seq1 - base_seq;
    c = residues[seq1[-1]];

    c_test = word_offset == 1 && !ss_displace && c >= 0; 

    test_left = word_offset > 1 && !ss_displace && c >= 0;

    known_match_size = word_bl[mid] & 127;

    sub_check = known_match_size < target_match_size;
    /* */
    for (i = mid + 1; 
	 i <= u_upper && !(word_bl[i] & 128); 
	 i++);
    u_upper = i - 1;
    /* */
    reuse = 100;

    for (i-- /* i = mid */; ; i--) {
      seq2 = seq_area + seq_offsets[i];
      /* if (reuse > full_word_bl[i]) reuse = full_word_bl[i]; */
      if (c_test && c == trans_index_res[seq2[-1]]) {
	if (word_bl[i] & 128) break;
	continue; 
      }
      if (4 == extend_match(seq1, seq2, ss_displace, test_left) || (word_bl[i] & 128)) break; /* entire word is within subject lc region -- so remaining ones also */
    }
    l_lower = i;
    /* 
    for (i = mid + 1; 
	 i <= u_upper && !(word_bl[i] & 128); 
	 i++) {
      seq2 = seq_area + seq_offsets[i];
      if (c_test && c == trans_index_res[seq2[-1]]) continue; 
      if (4 == extend_match(seq1, seq2, ss_displace, test_left)) break; 
    }
    u_upper = i - 1;
    /* */
  }
  else { 
    /* N.B. IN THIS CASE MATCHES CAN EXTEND TO LEFT AS WELL AS RIGHT -- SO COMPLEXITY AND REPEAT FILTERING
       ARE TOO HARSH!!
    */
    known_match_size = ind_wd_size;
    sub_check = 0; /*  in this case target_match_size already met */
    
    if (complexity_flag) {
      for (c = 0; c < alphabet_size; c++) base_res_counts[c] = 0;
      for (k = base_size = 0; base_size < x_minmatch && k < maxmatch; k++) { 
	c = residues[seq1[k]];
	if (c < 0) return; 
	base_res_counts[c] += 1;
	base_size += nlogn_diff[base_res_counts[c]] - nlogn_diff[k + 1];  
      }
      target_match_size = k; 
    }
    else target_match_size = minmatch;
    if (target_match_size < n_lc + 1) target_match_size = n_lc + 1;

/*
  if (compare_words(seqm - minmatch, seq_area + seq_offsets[l_lower], seqm, &k_end)
      || l_lower >= 1 && !compare_words(seqm - minmatch, seq_area + seq_offsets[l_lower - 1], seqm, &k_end)
      || compare_words(seqm - minmatch, seq_area + seq_offsets[u_upper], seqm, &k_end)
      || !compare_words(seqm - minmatch, seq_area + seq_offsets[u_upper + 1], seqm, &k_end))
    fatalError("wordmatches");
*/

    offset1 = seq1 - base_seq;
    c = residues[seq1[-1]];

    b = merge_residues[seq1[target_match_size - 1]];
    for (i = l_lower; i <= u_upper; i++) {
      seq2 = seq_area + seq_offsets[i];
      if (!(b & merge_residues[seq2[target_match_size - 1]])) continue; /* last position is most likely to be different */
      for (r = known_match_size; (d = merge_residues[seq1[r]]) & merge_residues[seq2[r]]; r++); /* merge_residues[seq1[r]] */
      if (r < target_match_size) { 
	/*
	if (i == l_lower && r < index_word_size[0]) {
	  fprintf(stderr, "\n%d %d %d %d", r, index_word_size[0], l_lower, u_upper);
	}
	*/
	seq2 = seq_area + r; /* last position to fail is most likely to fail next time */
	for (i++; i <= u_upper && !(d & merge_residues[seq2[seq_offsets[i]]]); i++);
	i--;
	continue; 
      }
      /*
	for (r = 0; r < known_match_size; r++)
	if (!(merge_residues[seq1[r]] & merge_residues[seq2[r]] & 15)) fatalError("sorting");
      */
      num_merge_words[1] += 1;

      /* SHOULD BE BETTER SCORING SCHEME HERE */

      /*
      for (r = s = known_match_size; d = merge_residues[seq1[r]] & merge_residues[seq2[r]]; r++)
	if (d & 15) s++;
      
      for (r = -1; d = merge_residues[seq1[r]] & merge_residues[seq2[r]]; r--)
	if (d & 15) s++;

      if (s < minmatch2) continue;
      */

      /*
      fprintf(stderr, "\n%d ", s);
      for (k = 0; k < r; k++)
	fprintf(stderr, "%c%c ", seq1[k], seq2[k]);
      */

      seq2 = seq_area + seq_offsets[i];
      num_merge_words[2] += 1;
      if (c == trans_index_res[seq2[-1]] && c >= 0) continue; /* does extension always imply index found?? */
      num_merge_words[3] += 1;

      extend_match(seq1, seq2, ss_displace, 0); /* CURRENTLY NOT WORKING IF SS_DISPLACE != 0 */
    }
  }

  clump_hist[trans_type][u_upper - l_lower + 1 > 100 ? 100 : u_upper - l_lower + 1] += 1;
  if (trans_type) num_merge_words[0] += u_upper - l_lower + 1;
}

extern int **full_score_mat, n_c_rows, row_convert[];
unsigned char *last_seq2, *last_seq1;
int last_j, last_target_score, last_end_j, last_end_loc, last_rgap;
/* Sw_edge sw_edges[100]; */
#define BIGNEG -1000000

int t_score, base_score, n_segs, n_diffs;
int start_loc, end_loc, lgap, base_start_j, base_end_j, end1;
int end_types[5], end_jumps[5], ends[5];

int extend_match(seq1, seq2, ss_displace, test_left)
     unsigned char *seq1, *seq2;
     int ss_displace, test_left;
{
  int entry2, offset2, complement, msize, left_score, right_score, i_penalty;
  int j, k, m, score, max_drop, max_score, c, d, e, target_score, entry_flag, keep_flag, sh_flag;
  int prev0, prev1, prev2, curr0, curr1, curr2, prev_e, max, pair_type, start1, start2, end2, loc, spl5, spl3, lsize, gap;
  int pos, i, i_diff, n, gap1_only, gap_1, gap_2, rgap, span_new, span_old;
  int *c_score_vec, *prev_c_score_vec, match_size, end_j, start_j, match_flag, base_target, globality;
  int i_seg, s, t, pin_left, pin_right, orig_score, seg_flag, off1, reuse_flag;
  int counts[256];
  static int base_counts[256], base_end;
  static unsigned char *base_seq1, *base_seq3;
  unsigned char *seq3;
  Aligned_pair *match_pair, *append_pair();  
  unsigned char *set_diff_block(), *diffs;
  unsigned char *get_id();
  Score_hist *append_score_hist(), *s_h;
  Query_domain *update_query_data(), *query_domain;
  Ssite *ssite;

  if (test_left) {
    j = -word_offset;
    for (i = -1; i >= j && (c = residues[seq1[i]]) >= 0 && c == trans_index_res[seq2[i]]; i--);
    if (i < j) return 1; /* use a different code?? */
   }

  if (parameters->subject_files) entry_flag = 0;
  else {
    find_entry_and_offset((SEQ_AREA)(seq2 - seq_area), &entry2, &offset2, &complement);

    if (entry2 <= entry1) return 1; /* if this is only info needed, can do more quickly!! */
    /* use only one of two equivalent (complementary) word matches, and
       don't allow matches to reverse complement of self for now; 
       fix this later (requires attention to num_pairs) */
    entry_flag = 1;
  }

  /*
  print_flag = 0;
  if (offset1 > 8197000 && offset1 < 8198000) {
    print_flag = 1;
    fprintf(stderr, "\n%d ", offset1);
  }
  */

  gap1_only = ss_displace || parameters->gap1_only; 

  if (gap1_only && (globality = parameters->globality)) {
    if (seq2 > comp_area) { /* query is complemented */
      pin_left = globality & 2;
      pin_right = globality & 1;
    }
    else {
      pin_left = globality & 1;
      pin_right = globality & 2;
    }
  }
  else pin_left = pin_right = 0;

  match_size = known_match_size;

  base_target = gap1_only ? parameters->min_record_score : parameters->gap1_minscore;

  if (!parameters->gap1_minscore && !gap1_only) goto meets_gap1_filter; /* NOT CORRECT IF PIN_LEFT OR PIN_RIGHT IS SET!! */

  if (1 /* last_seq1 != seq1 || merged_flag */) { 
    match_flag = 0;
    last_j = last_end_j = 1000;
    last_seq1 = seq1;
    n_skips[0] += 1;
  }
  else { /* NOT CORRECT IN MIXED CASE!! ALSO -- NEED TO HANDLE */ /* NOT CORRECT IF PIN_LEFT OR PIN_RIGHT IS SET!! */
    for (j = match_size; seq2[j] == last_seq2[j] && seq2[j]; j++); 
    /* if (j != reuse) fprintf(stderr, "(%d,%d)", j, reuse); */
    /* j = reuse; */
    if (seq2[j] != last_seq2[j]) j--;
    match_flag = !seq2[j + 1] && j >= last_end_j || j >= last_j; 
    if (match_flag) {
      n_skips[j >= last_j ? 2 : 1] += 1;
    }
    else match_size_hist[j > 100 ? 100 : j] += 1;
    /*
    if (j >= last_end_j) {
      if (j >= last_j) match_flag = 1;
      else {
	for (max = max_vec[j], j++; seq2[j] && max <= 0; j++, max++);
	if (max <= 0) match_flag = 1;
      }
    }
    */
    /* if (!last_seq2[j + 1]) fprintf(stderr, " %d", j); */

    /* fprintf(stderr, "[%d**%d,%d]", j, last_j, last_end_j); /* */
    /* fprintf(stderr, "(%d%s,%d%s,%d)", j, seq2[j + 1] ? "" : "*", last_j, last_seq2[last_j] ? "" : "*", last_end_j);  /* */

  }

  /*
  for (j = match_size; ; j++) {
    if ((score = full_score_mat[seq1[j]][seq2[j]]) > 0)
      base_score += score;
    else break;
  }
  match_size = j;
  */

  /* 
     if (base_score != match_size) {
     fprintf(stderr, "\n%d %d\n", max_score, match_size);
     fatalError("scoring");
     }
  */
  /* printf("\n%d %d ", offset1, base_score); */

  /* following is both slower & less thorough!!
    end_j = match_size - 1;
    for (j = match_size, score = max_score = 0; score >= dropoff && (e = seq2[j]) && (d = seq1[j]); j++) {
    score += c_score_mat[row_convert[d]][col_convert[e]];
    if (max_score < score) {
    max_score = score;
    end_j = j;
    }
    }
    base_score += max_score;

    end_j++;
    max_score = 0;
    if (seq2[end_j])
    for (j = end_j, score = gap_init; score >= dropoff && (e = seq2[j + 1]) && (d = seq1[j]); j++) {
    score += c_score_mat[row_convert[d]][col_convert[e]];
    if (max_score < score) max_score = score;
    }
    if (seq1[end_j])
    for (j = end_j, score = gap_init; score >= dropoff && (e = seq2[j]) && (d = seq1[j + 1]); j++) {
    score += c_score_mat[row_convert[d]][col_convert[e]];
    if (max_score < score) max_score = score;
    }
    base_score += max_score;

    end_j = 0;
    for (j = -1, score = max_score = 0; score >= dropoff && (e = seq2[j]) && (d = seq1[j]); j--) {
    score += c_score_mat[row_convert[d]][col_convert[e]];
    if (max_score < score) {
    max_score = score;
    end_j = j;
    }
    }
    base_score += max_score;

    end_j--;
    max_score = 0;
    if (seq2[end_j])
    for (j = end_j, score = gap_init; score >= dropoff && (e = seq2[j - 1]) && (d = seq1[j]); j--) {
    score += c_score_mat[row_convert[d]][col_convert[e]];
    if (max_score < score) max_score = score;
    }
    if (seq1[end_j])
    for (j = end_j, score = gap_init; score >= dropoff && (e = seq2[j]) && (d = seq1[j - 1]); j--) {
    score += c_score_mat[row_convert[d]][col_convert[e]];
    if (max_score < score) max_score = score;
    }

    base_score += max_score;
  */

  /* search band around diagonal defined by word, but assuming at most one gap char allowed in alignment */

  seq3 = seq1 + ss_displace;

  if (match_flag) {
    target_score = last_target_score;
    if (target_score <= 0 && !gap1_only) goto meets_gap1_filter; 
  }
  else {
    target_score = base_target;

    reuse = 100;
    last_seq2 = seq2;

    if (parameters->matrix) { 
      for (j = left_score = 0; j < left_ss_wdsize; j++) {
	left_score += full_score_mat[seq3[j]][seq2[j]]; 
	/*	fprintf(stderr, "%d:%d:%d ", seq3[j], seq2[j], full_score_mat[seq3[j]][seq2[j]]);  */
      }
      for (right_score = 0; j < match_size; j++) {
	right_score += full_score_mat[seq1[j]][seq2[j]]; 
	/* fprintf(stderr, "%d:%d:%d ", seq1[j], seq2[j], full_score_mat[seq1[j]][seq2[j]]); */
      }
    }
    else {
      left_score = left_ss_wdsize * parameters->match_reward;
      right_score = (match_size - left_ss_wdsize) * parameters->match_reward;
    }

    base_score = left_score + right_score;

    /* fprintf(stderr, "Left:%d,right:%d,Base:%d ", left_ss_wdsize, match_size, base_score); */
    target_score -= base_score;	 
    /*
    if (print_flag) {
      for (j = 0 ; j < 4; j++) fprintf(stderr, "%c", seq1[j]);
      fprintf(stderr, " %d %d [%d: %d] ", match_size, left_ss_wdsize, t_score, target_score);  
      for (j = 0 ; j < 4; j++) fprintf(stderr, "%c", seq2[j]);
    }
    */
    if (target_score <= 0 && !gap1_only) goto meets_gap1_filter; 
    base_end_j = end_j = match_size - 1;
    end_loc = 0;

    for (j = match_size; (e = seq2[j]) && (d = seq1[j]); j++) { 
      prev_c_score_vec = full_score_mat[d]; /* c_score_mat[row_convert[d]]; */
      prev0 = prev_c_score_vec[e];
      /* max_vec[j] = 0; */
      if (0 < prev0) { /* NOT CORRECT WITH MATRIX!! */
	base_score += prev0; 
	right_score += prev0; /* avoid inefficiencies -- with more bookkeeping!! */
	target_score -= prev0;
	base_end_j = end_j = j;
	if (target_score <= 0 && !gap1_only) {
	  last_end_j = last_j = j;
	  last_target_score = target_score;
	  goto meets_gap1_filter; 
	}
      }
      else {
	break;
      }
    }

    if (e && d) {
      gap_1 = gap_2 = j;
      prev_e = e;
      if (pin_right) {
	max_score = seq2[j + 1] ? BIGNEG : prev0;
	max_drop = BIGNEG;
	end_j = j;
      }
      else {
	max_score = 0;
	max_drop = dropoff;
      }
      /* sw_edges[j].curr1 = sw_edges[j].curr2 = */ prev1 = prev2 = gap_init;
      /* sw_edges[j].curr0 = prev0; */ 

      for (j++; /* (e = seq2[j]) && (d = seq1[j]) */; j++) { 
	e = seq2[j];
	d = seq1[j];
	c_score_vec = full_score_mat[d]; /* c_score_mat[row_convert[d]]; */
	curr1 = curr2 = prev0 + gap_init;
	max = prev0 += c_score_vec[e];
	loc = 0;
	if (curr1 < (c = prev1 + c_score_vec[prev_e])) {
	  curr1 = c;
	  if (max < curr1 && (!pin_right || !e)) { /* unnecessary work here in pin_right case */
	    max = curr1; /* placement here assumes gap_init <= smallest mismatch penalty */
	    loc = 1;
	  }
	}
	else gap_1 = j;
	if (curr2 < (c = prev2 + prev_c_score_vec[e])) {
	  curr2 = c;
	  if (max < curr2) { /* unnecessary work here in pin_right case */
	    max = curr2;
	    loc = 2;
	  }
	}
	else gap_2 = j;
	/* max_vec[j] = max; */
	if (!pin_right || !e || !seq2[j + 1]) {
	  if (max < max_drop) break; /* check this in pin_right case!!! */
	  if (max_score <= max) {
	    max_score = max;
	    max_drop = max_score + dropoff;
	    end_j = j;
	    if (gap1_only) {
	      end_loc = loc;
	      rgap = loc == 1 ? gap_1 : gap_2;
	    }
	    else if (max_score >= target_score && !gap1_only) {
	      last_target_score = target_score -= max_score;
	      right_score += max_score;
	      last_end_j = last_j = j;
	      goto meets_gap1_filter; 
	    }
	  }
	}
	if (!e || !d) break;
	prev_c_score_vec = c_score_vec;
	prev_e = e;
	/* sw_edges[j].curr1 = */ prev1 = curr1;
	/* sw_edges[j].curr2 = */ prev2 = curr2;
	/* sw_edges[j].max_score = max_score; */ 
	/* sw_edges[j].curr0 = prev0; */ 
      }
      target_score -= max_score; 
      right_score += max_score;
    }
    last_end_j = end_j;
    last_end_loc = end_loc;
    last_rgap = rgap;
    last_j = j;
    last_target_score = target_score;
  }
  
  /* should restore original version */
  base_start_j = start_j = 0;
  start_loc = 0;

  for (j = -1; (e = seq2[j]) && (d = seq3[j]); j--) { 
    prev_c_score_vec = full_score_mat[d]; /* c_score_mat[row_convert[d]]; */
    prev0 = prev_c_score_vec[e];
    if (0 < prev0) { /* this unnecessary in many cases */ /* NOT CORRECT WITH MATRIX!! */
      base_score += prev0;
      left_score += prev0;
      base_start_j = start_j = j;
      target_score -= prev0;
      if (target_score <= 0 && !gap1_only) goto meets_gap1_filter;  
    }
    else {
      break;
    }
  }

  if (e && d) {
    gap_1 = gap_2 = j;
    prev_e = e;
    if (pin_left) {
      max_score = seq2[j - 1] ? BIGNEG : prev0;
      max_drop = BIGNEG;
      start_j = j;
    }
    else {
      max_score = 0;
      max_drop = dropoff;
    }

    prev1 = prev2 = gap_init;

    for (j--; /* (e = seq2[j]) && (d = seq3[j]) */; j--) { 
      e = seq2[j];
      d = seq3[j];
      c_score_vec = full_score_mat[d]; 
      curr1 = curr2 = prev0 + gap_init;
      max = prev0 += c_score_vec[e];
      loc = 0;
      if (curr1 < (c = prev1 + c_score_vec[prev_e])) {
	curr1 = c;
	if (max < curr1 && (!pin_left || !e)) {
	  max = curr1; 
	  loc = 1;
	}
      }
      else gap_1 = j;
      if (curr2 < (c = prev2 + prev_c_score_vec[e])) {
	curr2 = c;
	if (max < curr2) {
	  max = curr2;
	  loc = 2;
	}
      }
      else gap_2 = j;
      if (!pin_left || !e || !seq2[j - 1]) {
	if (max < max_drop) break;
	if (max_score <= max) {
	  max_score = max;
	  max_drop = max_score + dropoff;
	  if (gap1_only) {
	    start_j = j;
	    start_loc = loc;
	    lgap = loc == 1 ? gap_1 : gap_2;
	  }
	  else if (max_score >= target_score && !gap1_only) goto meets_gap1_filter; 
	}

      }
      if (!e || !d) break;
      prev_c_score_vec = c_score_vec;
      prev_e = e;
      prev1 = curr1;
      prev2 = curr2;
    }
    target_score -= max_score;
    left_score += max_score;
  }

  if (target_score > 0) return 5;

 meets_gap1_filter:


  if (gap1_only) {

    start1 = start2 = start_j;
    if (start_loc == 1) start2 += 1;
    else if (start_loc == 2) start1 += 1;
    start1 += ss_displace;
    i_penalty = 0;

    if (ss_displace) {
      i_penalty = get_intron_penalty(-ss_displace);
      if (!parameters->globality) {
	if (left_score <= i_penalty || right_score <= i_penalty) return 6; 
	/* FOLLOWING COULD BE IMPROVED -- CHECK LEFT SIDE AS WELL */
	for (i = left_ss_wdsize - 1, j = 0; (n = full_score_mat[seq1[i]][seq2[i]]) >= 0; i--) 
	  j += n;
	if (left_score <= i_penalty + j) return 6;
      }
      /* check for polyT -- THIS NEEDS IMPROVEMENT --CHECK LENGTH NOT SCORE!! */
      c = seq2[start2];
      for (i = start2 + 1; seq2[i] == c; i++);
      if (i >= start2 + left_score - 1) {
	/*	  notify("polyT reject"); */
	return 6;
      }
    }

    if (max_intron2) { /* allowing both spliced & unspliced; spliced */
      /*      fprintf(stderr, " *%d*", last_ss - (offset1 + start1)); */
      off1 = offset1 + start1;
      for (i = i_off1 - 1; i >= 0 && off1s[i] < off1; i--);
      if (i < 0 || off1s[i] > off1) { /* if ==, ignore */
	for (j = i_off1 - 1; j > i; j--) off1s[j + 1] = off1s[j];
	off1s[i + 1] = off1;
	i_off1++;
	if (i_off1 >= n_off1) {
	  for (j = i_off1 - 1; j >= 0 && off1s[j] <= offset1 + ss_displace; j--);
	  for (i = 0, j++; j < i_off1; i++, j++)
	    off1s[i] = off1s[j];
	  i_off1 = i;
	  if (i_off1 >= n_off1) fatalError("resizing failed");
	}

	for (j = 0; j < 2; j++)
	  for (i = last_ssite[j] - 1; i >= first_ssite[j]; i--) {
	    ssite = ssites[j] + i;
	    if (ssite->pos > off1 + spliced_match_right) break;
	    if (ssite->pos >= off1 - spliced_match_left) 
	      ssite->keep_length = max_intron2;
	  }

      }
    }
	  
    end1 = end2 = last_end_j;
    end_loc = last_end_loc;
    if (end_loc == 1) end2 -= 1;
    else if (end_loc == 2) end1 -= 1;

    s_h = 0;
    query_domain = 0;
    reuse_flag = 0;
    orig_score = base_target - target_score - i_penalty; 
    seg_flag = 0;

    if (!parameters->raw) { /*  COMPLEXITY ADJUST */

      if (seq1 != base_seq1 || seq3 != base_seq3 || base_end > end1) {
	for (i = 0; i < n_c_rows; i++) base_counts[i] = 0;
	for (j = 0; j < left_ss_wdsize; j++) 
	  base_counts[row_convert[seq3[j]]] += 1;
	base_end = end1 > 20 ? 20 : end1; /* may not be optimal */
	for (; j <= base_end; j++) 
	  base_counts[row_convert[seq1[j]]] += 1;
	base_seq1 = seq1;
	base_seq3 = seq3;

	if (start1 - ss_displace > 0) fatalError("start id");

      }

      for (i = 0; i < n_c_rows; i++) counts[i] = base_counts[i];

      for (j = start1 - ss_displace; j < 0; j++) 
	counts[row_convert[seq3[j]]] += 1;
      for (j = base_end + 1; j <= end1; j++) 
	counts[row_convert[seq1[j]]] += 1;

      /* should subtract counts opposite gap in query */

      if (start_loc == 1)
	counts[row_convert[seq3[lgap]]] -= 1;

      if (end_loc == 1)
	counts[row_convert[seq1[last_rgap]]] -= 1;
      /*
      printf("\n%d %d %c %c\n", start_loc, end_loc , seq3[lgap], seq1[last_rgap]);
      */
      score = complexity_correct(orig_score, counts);
    }
    else {
      score = orig_score;
    }
    keep_flag = score >= parameters->minscore;
    if (parameters->score_flag && score >= parameters->min_record_score) { /* latter cond'n redundant? */
      if (!entry_flag) {
	find_entry_and_offset((SEQ_AREA)(seq2 - seq_area), &entry2, &offset2, &complement);
	entry_flag = 1;
      }
      
      sh_flag = parameters->score_hist;
      query_domain = update_query_data(entry2, score, complement, offset2 + start2, offset2 + end2);
      match_pair = query_domain->recent_pair;   
      if (match_pair && score == match_pair->score && offset1 + start1 == match_pair->start1 && offset1 + end1 == match_pair->end1 && entry1 == match_pair->entry1
	  && offset2 + start2 == match_pair->start2 && offset2 + end2 == match_pair->end2 && complement == is_reverse(match_pair)) {
	keep_flag = sh_flag = 0;
	if (query_domain->n_best > 1 && score == query_domain->best_score)
	  query_domain->n_best -= 1;
      }

      else if (query_domain->n_best > 1 && score == query_domain->best_score) { /*  need to check if likely a duplicate match */
	match_pair = query_domain->best_pair;  
	/* following may not be necessary -- with recent pair test */
	if (match_pair && offset1 + start1 == match_pair->start1 && offset1 + end1 == match_pair->end1 && entry1 == match_pair->entry1
	    && offset2 + start2 == match_pair->start2 && offset2 + end2 == match_pair->end2 && complement == is_reverse(match_pair)) {
	  query_domain->n_best -= 1;
	  keep_flag = sh_flag = 0;
	}
	else if (match_pair && ss_displace) {
	  span_new = end1 - start1 + 1;
	  span_old = match_pair->end1 - match_pair->start1 + 1;
	  if (span_new < span_old) {
	    query_domain->n_best = 1; /* this will be new best pair */
	    query_domain->next_best = query_domain->best_score; 
	  }
	  else if (span_new > span_old) {
	    query_domain->n_best -= 1;
	    query_domain->next_best = query_domain->best_score; 
	    keep_flag = 0;
	  }
	}
      }
      if (keep_flag && parameters->masklevel < 101 && 
	  (
	   score < query_domain->best_score + (parameters->minmargin >= 0 ? 0.0 : parameters->minmargin)
	   || parameters->minmargin >= 1 && query_domain->n_best > 1  /* note: incomplete filter in this case -- but need filtering at printout anyway */
	   || parameters->minmargin == 0.5 && random() % query_domain->n_best)
	  )
	keep_flag = 0;

      if (sh_flag) {
	s_h = append_score_hist(query_domain, score, ss_displace ? 1 : 0, 1);
	if (s_h->count > 1) { /* need to check if was likely a duplicate match */
	  match_pair = s_h->best_pair;
	  if (match_pair && offset1 + start1 == match_pair->start1 && offset1 + end1 == match_pair->end1 && entry1 == match_pair->entry1
	      && offset2 + start2 == match_pair->start2 && offset2 + end2 == match_pair->end2 && complement == is_reverse(match_pair)) {
	    s_h->count -= 1;
	    keep_flag = 0;
	  }
	}
      }
      /* displace unspliced pair with higher-scoring one having same end */
      /*
      match_pair = query_domain->recent_pair;  
      if (parameters->fuse_gap1 && keep_flag && match_pair && 
	  (offset1 + start1 == match_pair->start1 && offset2 + start2 == match_pair->start2 
	   || offset1 + end1 == match_pair->end1 && offset2 + end2 == match_pair->end2)
	  && entry1 == match_pair->entry1 && complement == is_reverse(match_pair)) {
	if (!ss_displace && score < match_pair->score - 5) { 
	  keep_flag = 0;
	  notify(".");
	}
	else if (is_unaligned(match_pair) && core > match_pair->score + 5) {
	  reuse_flag = 1; 
	  notify("-");
	}
      }
      */
/* should revise bookkeeping here!! i.e. discard data for previous one */
    }
    if (keep_flag) {

      if (!entry_flag) {
	find_entry_and_offset((SEQ_AREA)(seq2 - seq_area), &entry2, &offset2, &complement);
	entry_flag = 1;
      }
      if (!seg_flag) {
	get_segs(ss_displace);
	seg_flag = 1;
      }
      for (i_seg = gap = 0, j = start1, k = start2; ; i_seg++) {
	for (; j < ends[i_seg]; j++, k++) {
	  t_score += s = full_score_mat[seq1[j]][seq2[k]];
	  gap++;
	  if (s <= 0 || gap > 60) {
	    n_diffs++;
	    gap = 0;
	  }
	}
	if (i_seg == n_segs - 1) break;
	m = end_types[i_seg];
	if (m == 1) {
	  j++;
	  gap = 0;
	}
	else if (m == 2) {
	  k++;
	  gap = 0;
	}
	else if (m) {
	  j += end_jumps[i_seg];
	  if (m == 4) {
	    k += gap = ss_displace + end_jumps[i_seg]; /* no intronic jump in query */
	  }
	  else {
	    k += end_jumps[i_seg]; /* no intronic jump in query */
	    gap += end_jumps[i_seg];
	  }
	  while (gap > 60) {
	    n_diffs++;
	    gap -= 60;
	  }
	}
      }

      t_score -= i_penalty;
      if (parameters->raw && t_score != score) {
	fprintf(stderr,"\n%d %d  %d %d %d %d %d %d %d %d\n", 
		t_score, score, start1, start2, end1, end2, j, k, lgap, rgap);
	for (j = 0; j < n_segs; j++) fprintf(stderr, "%d:%d:%d ", ends[j], end_types[j], end_jumps[j]);
	for (j = start1; j <= end1; j++) fprintf(stderr, "%c", seq1[j]);
	notify("\n");
	for (j = start2; j <= end2; j++) fprintf(stderr, "%c", seq2[j]);
	fatalError("score discrepancy");
      }
      
      if (reuse_flag) {
	match_pair->flags = match_pair->reject_flags = 0;
	match_pair->spl5 = match_pair->spl3 = 0;
      }
      else {
	num_pairs++; 
	match_pair = append_pair(entry1, entry2, 1);
      }
      set_reverse_flag(match_pair, complement);
      match_pair->score = score;
      match_pair->start1 = offset1 + start1; /* origin 0 */
      match_pair->start2 = offset2 + start2;
      match_pair->end1 = offset1 + end1;
      match_pair->end2 = offset2 + end2;
      if (s_h) s_h->best_pair = match_pair;
      if (query_domain) {
	match_pair->query_data->query_domain = query_domain;
	query_domain->recent_pair = match_pair;
	if (score == query_domain->best_score)
	  query_domain->best_pair = match_pair;
      }
      if (ss_displace) {
	set_unaligned_flag(match_pair, 1); /* flag as spliced */
	lsize = left_ss_wdsize;
	spl5 = strand ? 63 - (16 * residues[seq1[lsize - 1]] + 4 * residues[seq1[lsize - 2]] + residues[seq1[lsize - 3]])
	  : 16 * residues[seq3[lsize]] + 4 * residues[seq3[lsize + 1]] + residues[seq3[lsize + 2]];

	spl3 = strand ? 63 - (16 * residues[seq3[lsize + 2]] + 4 * residues[seq3[lsize + 1]] + residues[seq3[lsize]]) 
	  : 16 * residues[seq1[lsize - 3]] + 4 * residues[seq1[lsize - 2]] + residues[seq1[lsize - 1]];

	if (spl5 < 64 && spl3 < 64) {
	  match_pair->spl5 = spl5 + 128;
	  match_pair->spl3 = spl3 + 128;
	}
	else match_pair->spl5 = match_pair->spl3 = 128;
	if (strand) {
	  match_pair->spl5 += 64;
	  match_pair->spl3 += 64;
	}
      }
      /* if (n_diffs < 0) fatalError("n_diffs"); */
      match_pair->diffs = diffs = set_diff_block(n_diffs);

      for (i_seg = i_diff = gap = 0, j = start1, k = start2; ; i_seg++) {
	for (; j < ends[i_seg]; j++, k++) {
	  s = full_score_mat[seq1[j]][seq2[k]];
	  /*	  if (!strcmp(get_id(entry2), "8-1-180-879_1")) fprintf(stderr, "%s:%d ", get_id(entry2), s); */
	  gap++;
	  if (s <= 0) {
	    diffs[i_diff++] = set_diff(gap, 'S');
	    gap = 0; /* need 'M' case also */
	  }	
	  else if (gap > 60) {
	    diffs[i_diff++] = set_diff(gap, 'M');
	    gap = 0;
	  }
	}
	if (i_seg == n_segs - 1) break;
	m = end_types[i_seg];
	if (m == 1) {
	  j++;
	  gap++;
	  diffs[i_diff++] = set_diff(gap, 'D');
	  gap = 0;
	}
	else if (m == 2) {
	  k++;
	  diffs[i_diff++] = set_diff(gap, 'I');
	  gap = 0;
	}
	else if (m) {
	  j += end_jumps[i_seg];
	  if (m == 4) {
	    /*	  gap = -start1 + left_ss_wdsize + 1; /* NOT NEC CORRECT -- DEPENDS ON gaps, prev discreps */
	    diffs[i_diff++] = set_diff(gap + 1, 'D');   
	    n = parameters->word_intron_margin;
	    for (i = 0; i < n - 1; i++)
	      diffs[i_diff++] = set_diff(1, 'D');
	    diffs[i_diff++] = set_diff(0, 'm');
	    gap = 20;
	    for (pos = -ss_displace - 2 * n; pos > 0; pos -= gap) {
	      gap = pos < 20 ? pos : 20;
	      diffs[i_diff++] = set_diff(gap, 'M');
	    }
	    diffs[i_diff++] = set_diff(0, 'm');
	  
	    for (i = 0; i < n; i++)
	      diffs[i_diff++] = set_diff(1, 'D');
	    k += gap = end_jumps[i_seg] + ss_displace;;
	  }
	  else {
	    k += end_jumps[i_seg];
	    gap += end_jumps[i_seg];
	  }
	  while (gap > 60) {
	    diffs[i_diff++] = set_diff(60, 'M');
	    gap -= 60;
	  }
	}
      }

      diffs[i_diff++] = set_diff(gap + 1, 'S');
      if (i_diff != n_diffs + 1) {
	fprintf(stderr, "\n%d %d\n", i_diff, n_diffs);
	for (j = 0; j < n_segs; j++) fprintf(stderr, "(%d:%d:%d) ", ends[j], end_types[j], end_jumps[j]);
	fatalError("diff allocation");
      }
     /* 
      printf("\n%d %d %d  %d..%d %d..%d  %s\n", offset1, parameters->gap1_minscore - target_score, -ss_displace, 
	     match_pair->start1, match_pair->end1, match_pair->start2, match_pair->end2, get_id(entry2));
      for (j = 4; j < 35; j++) printf("%c", seq1[j]);

      printf("\n%s", seq2 + 4);
      printf("\n");
      for (j = lsize + 2; j >= -20; j--) printf("%c", seq3[j]);
      printf("\n");
      for (j = lsize + 2; j >= -20 && seq2[j]; j--) printf("%c", seq2[j]);

      /* */

      /* need to count these cases*/

      return 6;
    }
    else return 7;
  }

  if (!entry_flag) {
    find_entry_and_offset((SEQ_AREA)(seq2 - seq_area), &entry2, &offset2, &complement);
    entry_flag = 1;
  }

  if (!parameters->subject_files) { /* no subject file, so must be complement of query file sequence */
    n_cross_matches++;
    if (q_repeat_screen) {
      msize = get_matchsize_nocomplexity(seq1, seq2, match_size);
      match_size_hist[msize > 100 ? 100 : msize] += 1;

      if (msize < target_match_size) return 2; /* SHOULD NOT HAPPEN?? */
      if (contained_in_tag(entry1, "repeat", offset1 + 1, offset1 + msize, 1, 0)
         || contained_in_tag(entry2, "repeat", offset2 + 1, offset2 + msize, 0, 0)) {
	num_q_repeat_fail_words++;

	return 3;
      }
    }
    pair_routine(entry2, entry1, offset2, offset1, 1, 0); /* not reversing storage */
  }
  else {
    n_cross_matches++;
    msize = 0;

    if (sub_check) {
      msize = get_matchsize_nocomplexity(seq1, seq2, match_size);
      match_size_hist[msize > 100 ? 100 : msize] += 1;
      if (msize < target_match_size) {
	num_s_repeat_fail_words++;

	return 4; 
      }
    }

    if (q_repeat_screen) {
      if (!msize) {
	msize = get_matchsize_nocomplexity(seq1, seq2, match_size);
	match_size_hist[msize > 100 ? 100 : msize] += 1;
      }
      /*
	if (!complexity_pass) {
	num_le_test_words++;
	if (get_corrected_matchsize() < minmatch) {
	num_le_fail_words++;
	return;
	}
	}
      */
      if (contained_in_tag(entry2, "repeat", offset2 + 1, offset2 + msize, complement, 0)) {
	num_q_repeat_fail_words++;
	return 3;
      }
    }
    pair_routine(entry1, entry2, offset1, offset2, complement, 1); /* REVERSING STORAGE */
  }
  num_pass_words++;

  return 0;
}

init_word_counts()
{
  int i;

  num_pass_words = num_le_fail_words = num_le_test_words = num_q_repeat_fail_words = num_s_repeat_fail_words = 00;
  for (i = 0; i < 5; i++) num_merge_words[i] = 0;
}

print_word_counts()
{
  int i;
  double cum;

  if (!fp_debug) return;

  fprintf(fp_debug, "\n%.0f pass words, %.0f low-entropy-tested words, %.0f low-entropy-fail words, %.0f query repeat_fail_words, %.0f subject repeat_fail_words, %.0f vector matches\n", 
	  num_pass_words, num_le_test_words, num_le_fail_words, num_q_repeat_fail_words, num_s_repeat_fail_words, vector_matches);
  fprintf(fp_debug, "\n %.0f opp_sense (external) word matches: ", n_cross_matches);
  fprintf(fp_debug,"\n %.0f same-sense (internal) word matches: ", n_same_matches); 
  fprintf(fp_debug,"\n # merged word iterations (in millions) of different levels:");
  for (i = 0; i < 4; i++) fprintf(fp_debug," %.6f", num_merge_words[i] / 1000000.0);
  fprintf(fp_debug,"\n %.0f spliced word calls", n_splice_calls);
  fprintf(stderr,"\n %.0f spliced word calls", n_splice_calls);
  print_num_cand_pairs();
  fprintf(fp_debug, "\n\nClump size histogram:");
  for (i = cum = 0; i < 101; i++)
    if (clump_hist[0][i] + clump_hist[1][i]) {
      cum += clump_hist[0][i] + clump_hist[1][i];
      fprintf(fp_debug, "\n%3d %10.0f %10.0f  %10.0f", i, clump_hist[0][i], clump_hist[1][i], cum);
    }
  fprintf(fp_debug, "\n\nWord match size histogram:");
  if (1 /* complexity_flag || q_repeat_screen || s_repeat_screen */) {
    for (i = cum = 0; i < 101; i++)
      if (match_size_hist[i]) {
	cum += match_size_hist[i];
	fprintf(fp_debug, "\n%3d %10.0f   %10.0f", i, match_size_hist[i], cum);
      }
    fprintf(fp_debug, "\nNo. of skips: %d %d %d", n_skips[0], n_skips[1], n_skips[2]);
  }
}

get_segs(ss_displace)
     int ss_displace;
{
  int n;

  t_score = base_score;
  n_segs = 0;
  n_diffs = 0;

  if (start_loc) {
    end_types[n_segs] = start_loc;
    end_jumps[n_segs] = 0;
    ends[n_segs++] = ss_displace + (start_loc == 1 ? lgap : lgap + 1);
    t_score += gap_init;
    n_diffs++;
  }

  end_types[n_segs] = 3;
  end_jumps[n_segs] = (ss_displace ? left_ss_wdsize : base_end_j + 1) - base_start_j;
  ends[n_segs++] = ss_displace + base_start_j;
      
  if (ss_displace) {
    end_types[n_segs] = 4;
    end_jumps[n_segs] = -ss_displace + base_end_j + 1 - left_ss_wdsize;
    ends[n_segs++] = left_ss_wdsize + ss_displace;
    n = parameters->word_intron_margin;
    n_diffs += 2 * n + 2 + (-ss_displace - 2 * n + (20 - 1)) / 20;
  }
      
  if (end_loc) {
    end_types[n_segs] = end_loc;
    ends[n_segs++] = last_rgap;
    t_score += gap_init;
    n_diffs++;
  }
  end_types[n_segs] = 0;
  ends[n_segs++] = end1 + 1;
}

