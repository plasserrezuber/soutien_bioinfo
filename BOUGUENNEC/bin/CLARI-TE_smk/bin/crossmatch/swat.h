/*****************************************************************************
#   Copyright (C) 1993-2008 by Phil Green.                          
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
  
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int SEQ_AREA; /* options: int, unsigned int, long int; long int is slower & takes more memory */

typedef signed char SCORE; /* int may take less time than signed char (but more space -- for long seqs) */
/* NO LONGER HAVE SEPARATE MANYREADS VERSION
#if defined(MANYREADS)
*/
typedef unsigned int Entry_index; /* indices of reads -- int required if > 64,000 reads  */
/*
#else
typedef short unsigned int Entry_index; 
#endif
*/

/* NO LONGER HAVE SEPARATE LONGREADS VERSION
#if defined(LONGREADS)
*/
typedef unsigned int Read_position; /* read position (relevant to phrap 
					     node structure only) -- int required if any read has > 64,000 bases */
/*
#else
typedef short unsigned int Read_position; 
#endif
*/

typedef struct profile {
  double lambda; /* calculated lambda for score matrix */
  double *log_col_freqs; /* log residue freqs for score matrix */
  int gap_init, gap_ext, ins_gap_ext, del_gap_ext; /* gap initiating and gap extension penalties;
			    assumed constant across profile. Should be
			    negative numbers. gap_init is the penalty for
			    the first residue in a gap, gap_ext is the
			    penalty for each succeeding residue; ins_gap_ext
			    is for insertion (in subject rel to query),
			    del_gap_ext is for deletion (in subject rel to query)*/
  int qbegin_gap_init,qbegin_gap_ext,qend_gap_init,qend_gap_ext,
      sbegin_gap_init,sbegin_gap_ext,send_gap_init,send_gap_ext;
/* penalties at beginning, end of query and subject, for Needleman-Wunsch.
   e.g. qbegin penalties are for gaps in query preceding first residue in query */

  int alloc_length, length; /* length of profile; & allocated length (may be larger) */
  int n_c_cols; /* # of indices produced by convert (in score matrix) */
  int *convert; /* array to convert letters (e.g. in target aa seq being
		   compared to the profile) to indices for score matrix */
  SCORE **scores; /* array (of length 20, or size of alphabet)
		   of pointers to vectors of scores (corresponding to
		   rows in profile; each has length profile->length + 1) */
  SCORE *score_area;
  int **hist;
  int max_entry, min_entry; /* maximum positive entry, and minimum negative entries
			       (or 0, if none) in scores */
			       
/*  int **score_pos;  array of pointers to vectors of indices of positive
		      score positions in scores -- was used only with SPARSE */
  unsigned char *seq; /* characters in sequence which corresponds to the profile (if any) */
  /*  int *maxstu_vec; NOT NEEDED /* array of length length + 1, used to hold row of scores,
		      in smith-waterman alg; last entry must be -1 */
  int *t_vec; /* array of length length + 1, used to hold row of t scores,
		 in smith-waterman alg; last entry must be > -1 */
/*  int **work3, **work4;  working space, must be of length profile->length */
/*  int **score_mat;  score matrix */

/*  int packed_length; */
  unsigned long int poly_offset, poly_one, max_score_cutoff, poly_gap_init, poly_gap_ext;
/* N.B. NEED POLY_INS_GAP_EXT ETC. */
/*
  unsigned long int *packed_maxstu_vec;
  unsigned long int **packed_scores;
*/
} Profile;

typedef struct file {
  char *name;
  FILE *fp;
  struct file *next;
  struct database *db;
  int type; /* 0: fasta; 1: calf */
} File;

typedef struct database {
  char *id_buffer, *descrip_buffer;
  unsigned char *seq_buffer, *comp_seq_buffer; /* buffers containing current entry (when file is being read) */
  int seq_length, comp_seq_length; /* length of current entry */
  int seq_buffer_entry, comp_seq_buffer_entry; /* current contents of buffer and complement */
  int id_buffer_size, descrip_buffer_size, seq_buffer_size, comp_seq_buffer_size;
  int first_entry, last_entry; /* index (in list of all entries) of first and last entries
				  in this database */
  int num_entries; /* no. of entries in database 
			(before complementing) */
  SEQ_AREA t_length; /* total length of entries (before complementing) */
  int in_memory; /* flag indicating whether in memory */
  File *file; /* FASTA file used */
  char *id_area, *descrip_area;
  unsigned char *seq_area, *comp_area;
/* arrays containing ids, descriptions all residue sequences; first char is 0, 
 which is followed by each sequence in turn, each being terminated by 0; comp_area points to start of comp seqs in seq_area (or to end, if not stored) */
  SEQ_AREA id_area_size, descrip_area_size, seq_area_size; /* size (in bytes) of different areas */
  char *orig_qual_area, *adj_qual_area;
  struct align_info *align_array;
  int complements; /* flag indicating whether includes complements */
  struct database *next;
} Database;
 
/*
typedef struct db_entry {
  char *id;
  char *descrip;
  char *seq;
  long int fpos; ** offset in file stream of this entry **
  int score, length;
  float z; ** z-score **
} Db_entry;
*/

/* diff structure -- for recording positions at which two aligned sequences differ */

typedef int Seq_position; /* structure for indicating position within sequence */

/*
typedef struct diff {
  Seq_position site1, site2;
  char type;  = 'S', 'D', 'I' for Substitution, Deletion, Insertion in 2d
		sequence relative to first; 'B','E' for Beginning, End of alignment (
		(positions are one base before or after) 
} Diff; 
*/

typedef struct max_list {
  int q_pos, s_pos, score;
} Max_list;

/* segment structure -- used for keeping track of appropriate bands for SWAT, and of contig segments
   with similarity to other regions */

typedef struct segment {
  int start, end;
  struct segment *next;
} Segment;

typedef struct node {
  int score, best, link;
  Entry_index entry1; 
  Read_position site; /* N.B. THIS NOT ENOUGH IF READS ARE EXTREMELY LONG (> 64K) !!! */
/*  int position; position in contig sequence */
/*  int swat_score; */
/*  char used;  for detection of infinite loops -- altho there should be
		a better way! */
} Node; 

typedef struct splice_site {
  struct splice_site *next;
  /* char seq[10]; */
  float llr;
  int offset, intron_size;
  char strand, side, type; /* type: 5'ss: GT is 0, GC is 1 (penalty 2), ATATCC is 2;  for 3' ss: HAG is 3, GAG is 4 (penalty 2), AC is 5; 
				 allowed pairings: 5 only with 2, any other 5' with 3' OK */
  char penalty, cdna_overhang_adjust, genome_overhang_adjust, overhang_check[3], m_len, loc, j_offset; /* loc, j_offset must be signed -- others need not be, and some can be bits */
} Splice_site;

/* structure for linked list of aligned pairs */

typedef struct aligned_pair {
  struct aligned_pair *next, *reversed_pair; /* next is the next in the list assoc.
						to entry1; reversed_pair is the
						pair with entry1 and entry2 reversed */
  unsigned char *diffs; /* differences in initial  alignment */
  struct query_data *query_data;
  int start1, end1, start2, end2; 
/*   starting, ending positions of alignments in the two entries, initially,
     start1 and start2 are starting positions of alignments in the two
     sequences; if second sequence is complemented (i.e. pair->reverse = 1)
     then count from beginning of complement rather than original sequence.
     First sequence is always assumed to be uncomplemented (i.e. original
      database sequence).
     Always have start < end */
  int score; /* SWAT scores */
  int offset; /* average offset (entry1 position - entry2 position) for the alignment */

  Entry_index entry1, entry2; /* indices (in the database, db) 
					of the two entries */
  short int LLR_score; /* loglikelihood ratio score of alignment -- based on quality matrix */
  char flags; /* bit flags, including following:
             bit #         flag
              0       reverse    0 if both entries are in original orientation
		                  1 if 2d entry is complemented  
              1       used    1 if used in (or consistent with) contig merge; 
              2       best    1 if is highest scoring alignment (among all which overlap it) 
                                     for these particular entries
	      3       triple_reject    1 if there is another read with positive LLR_score against
	                       one read in pair and negative LLR_score against other read in pair
              4       repeat  1 if reads are from different repeats in same contig
	      5       chimeric 1 if match involves "non-retained" segment in an
	              apparently chimeric read.
              6       split   1 for pairs involved in splitting assembly in revise_greedy_algorithm
              7       unaligned  1 if substantial part of the matching region in entry1 is
                        not aligned to the contig. OR (cross_match only) is spliced match
*/
  char reject_flags; /* bit flags related to rejection, including following:
             bit #         flag
              0       rejectable    0 if pair does not have sufficient high_quality
	                              bases to allow rejection
	      1       self       1 if match is to self, or to perfect duplicate
	      2       vector     1 if match involves part of read on other side of vector (X'd out) stretch
              3       mismatch   1 if alignment is rejected for mismatches alone
	      4       total      1 if alignment is rejected for total discrepancies
              5       left_trunc   1 if alignment is prematurely truncated to left
              6       right_trunc  1 if alignment is prematurely truncated to right
	      7       node        1 if alignment is "node-rejected", i.e. there are
	                  no legal nodes (no solid regions of alignment)
*/
  unsigned char spl5, spl3; /* 128 bit indicates if spliced; 64 bit indicates strand; other bits indicate 3-base sequence at intron bdry */
} Aligned_pair;

/* candidate pair (having possible word match) */
typedef struct cand_pair {
  struct cand_pair *left, *right; /* left, right in tree next is the next in the list assoc. to entry1 */
  Segment *band_segments; /* linked list of segments (bands around offset) to
		       be searched by Smith-Waterman. */
  Entry_index entry1, entry2; /* indices (in the database, db) 
					of the two entries */
  char reverse;
} Cand_pair;

typedef struct tag {
  struct tag *next;
  char type[30]; /* type of tag */
  int start, end; /* starting, ending bases of tag */
} Tag;

typedef struct align_info {
  char *template_start, *template_end; /* beginning, end of template name */
  char *orig_qual, *adj_qual; /* original & adjusted quality measures: 
                                             (in range 0 to 100), for each base */  
  Segment *segments; /* (merged) segments that overlap other reads */
  struct contig *contig; /* contig containing this entry */
  struct align_info *next; /* next entry in the same contig */  
  Tag *tags; /* linked list of tags associated to this read */
  unsigned char *diffs; /* differences between contig sequence (site1)  
		  and current read (site2); origin 1 */
  int seq_entry; /* entry with read info (for the forward read only!) */
  int first_start, last_end; /* first, last aligned bases, w.r.t. all pairwise 
				alignments involving this entry */
  int rev_first_start, rev_last_end; /* first, last aligned bases, w.r.t. all reverse sense pairwise
				alignments involving this entry */
  int qual_start, qual_end; /* beginning, end of high quality segment */
  int first_vec, last_vec; /* apparent first, last vector bases (based on anomalous
			      matches) */
  int start, end; /* starting and ending positions relative to contig 0 base (which
		     is NOT necessarily the first base of the contig;
		   start < end regardless of sequence orientation */

/* following pertain to alignments against contig */
  int score, LLR_score; /* SWAT and LLR scores against contig */
  int m_start, m_end; /* positions of first, last matching bases relative to
			 contig */
  int m_length; /* length of matching segment (in consensus) -- is this nec? */
  int equiv_class; /* equivalence class -- for clustering purposes */
  char anomalies; /* +1 for deletion; +2 for chimera; +4 for duplicate read;
		     ; +16 for singlet
		   (no non-vector match to anything except duplicates) */
/*  int *depth;  (cumulative) depth of confirmed bases at each point */
  char chimera_bits; /* indicates chimeric segments */
  char temp; /* working character */

/* following pertain to placement of read in contig */
  char reverse; /* reverse = 0 if entry is in original orientation in contig,
		           = 1 if entry is complemented */
  char blocked; /* blocked (no extension past one end) */
  char bypassed; /* indicates read bypassed during assembly */
  char chemistry; /* 0 = dye primer; 1 = old dye terminator; 2 = bigDye term */
  char direction; /* 0 = forward read, 1 = reverse read */
} Align_info;

/* base_segment structure - to keep track of segments of read sequences from which 
   contig sequence is derived */

typedef struct base_segment {
  int entry; /* entry from which sequence is derived */
  int read_start, read_end, contig_start, contig_end; /* base positions (origin 1) of
							 segments */
  struct base_segment *next;
} Base_segment;

typedef struct contig {
  struct merge_reject *merge_reject;
  Segment *top_segments, *bottom_segments; /* segments covered by matching parts of reads */
/*   int *pad_translate; translate positions in unpadded seq to those in padded seq */
  /* signed */ char *orig_qual, *adj_qual, *discrep_qual; /* original and adjusted qualities -- inherited
     from appropriate reads (adj_qual may be modified); discrep_qual is highest adj_qual
     of a discrepant base in a read */
  Base_segment *base_segment;
  struct tig_node *tig_node;
  unsigned char *seq;
  char *id, *descrip;
  Align_info *first, *last; /* first, last members in linked list of entries for this
			       contig */
  struct contig *parent; /* parent contig class for linked group */
  int first_start, last_end; /* first, last aligned bases w.r.t. contig 0 */
  int num_entries; /* number of entries in contig */
  int num_matches; /* no. of matching pairs in contig */
  int length; /* best guess for length */
  int t_num_entries; /* total number of entries (reads)
				     in contig group (only defined for parent) */
  int score; /* Quality score of sequence */
  int n_pads; /* no. of pads inserted into sequence */
  Entry_index index; /* numerical index of final contigs */
  Entry_index num_contigs; /* no. of contigs in the class
				     (only defined for parent)*/
  char comp_status; /* complementation status */
} Contig;

typedef struct merge_reject {
  struct merge_reject *reverse, *next, *prev; /* next, prev refer to doubly linked list 
						 position for contig1;
						 reverse refers to contig2 list */
  Contig *contig1, *contig2;
  Aligned_pair *reject_pair; /* pair that caused rejection */
  int comp_status1, comp_status2;
  int offset;
  int reject_reason;
  int highest_LLR_score, lowest_LLR_score;  /* n_pairs?? */
  int gap; /* also want location; equiv want aligned part 
	      or biggest piece of it */
  int join_score;
} Merge_reject;

typedef struct seq_class {
  struct seq_class *parent;
  int offset, class_size;
} Seq_class;

typedef struct parameters {
  /* following are user-modifiable */
  int alignments; 
  int align_extend;
  int bandwidth;
  int bcdsites_qual_threshold;
  int bypasslevel; 
  int compact_qual;
  int confirm_length;
  int confirm_penalty;
  int confirm_score; 
  int confirm_trim;
  int contig_graph_weights; 
  int create_test_file;
  int default_qual; 
  int del_gap_ext;
  int discrep_lists;
  int discrep_tables; 
  float e_cutoff;
  int end_gap; 
  int exp_input;
  int exp_output;
  int file_flag;
  int force_high; 
  int forcelevel; 
  int fuse_gap1;
  int gap_ext;
  int gap_init; 
  int gap1_dropoff;
  int gap1_minscore;
  int gap1_only;
  int globality; 
  int group_delim; 
  int indexwordsize;
  int indexwordsize2;
  int ins_gap_ext; 
  char *logfile; /* name of .log file */
  int masklevel;
  char *matrix; /* score matrix for swat & cross_match alignments */
  int max_group_size;
  int max_intron_length;
  int max_num_alignments;
  int max_overlap; 
  int max_subclone_size; 
  int maxgap;
  int maxmatch;
  int min_exon_length; 
  int min_intron_length; 
  float minmargin; 
  int minmatch;
  int minscore; 
  int n_delim; 
  int near_minscore;
  int new_ace;
  int node_seg;
  int node_space;
  int nw_flag;
  int old_ace;
  int output_bcdsites; 
  int output_nonmatching_queries; 
  int penalty; 
  int preassemble; 
  int print_extraneous_matches; 
  int print_word_data;
  int qual_scores;
  int qual_show;
  int raw; 
  int repeat_screen; 
  float repeat_stringency;
  int retain_duplicates; 
  int revise_greedy; 
  int score_hist;
  int screen;
  int shatter_greedy; 
  int splice_edge_length; 
  int spliced_match_left;
  int spliced_match_right;
  int spliced_word_gapsize;
  int spliced_word_gapsize2;
  int subclone_delim;
  int tags; 
  int trim_penalty;
  int trim_qual;
  int trim_score;
  int trim_start;
  int truncatedb;
  int vector_bound;
  int view; 
  int word_intron_margin; 
  int word_offset;
  int word_raw; 
  float z_cutoff;

  /* following are non(directly)-user-modifiable) */
  int query_histograms;
  int (*align)(),(*full_align)();
  int splice_edge_alloc; /* length of edges of alignment boundaries to scan for splices (at each end, scan twice this length -- internal + external */
  int score_flag; /* depth of score monitoring to do */
  int min_record_score;
  int keep_query_data;
  int find_z_flag;
  int gap1_flag;
  int DNA_flag; /* 1 if ord DNA, 2 if merged DNA, 0 if protein */
  int use_e, use_z, use_n;
  int qual_flag; /* 1 if .qual file was provided, 0 otherwise */
  int X_penalty, match_reward; /* used in defining score matrix values */
  float qual_resolution; /* resolution in quality values */
  char calling_program[20], version[20], date[20];
  int argc;
  char **argv;
  File *query_files;
  File *subject_files;
/*   char **file_names; names of input sequence files 
  int num_files;
*/ 
/* following currently are for gcphrap only, & allow experiment file input/output */
  char *exp_dir;
} Parameters; 

typedef struct score_hist {
  int score, count;
  Aligned_pair *best_pair; /* best pair with current score -- currently only needed 
			      for dealing with gap1_only duplicate hits */
  struct score_hist *next;
} Score_hist;

typedef struct query_data {
  char *edges12; /* edges of alignments -- a single string with query, subject sequences at both ends of alignment, used in checking for splicing; */  
  Splice_site *splice_sites; 
  struct query_domain *query_domain; 
} Query_data;

typedef struct query_domain {
  Aligned_pair *best_pair; /* best (or currently selected) pair with current score */
  Aligned_pair *recent_pair; /* most recent pair in this domain */
  Score_hist *score_hist[2]; /* histogram of scores: 0 = unspliced, 1 = spliced */  
  struct query_domain *child[2], *parent; /* parent is NOT the tree parent -- instead, the equiv class parent */
  int start, end, trim_start, trim_end, best_score, next_best, n_best;
} Query_domain;

typedef struct seq_entry {
  Aligned_pair *aligned_pairs; /* (unordered) linked list of significant matches involving  
				this entry (all of them have have pair->entry1 =
				index of this entry) */
  Cand_pair *cand_pairs; /* (unordered) linked list of candidate pairs involving
				this entry (all of them have have pair->entry1 =
				index of this entry) */
  Query_domain *query_domains;
  long int id_pos, descrip_pos, seq_pos; /* long int, because must hold fseek position */
  int seq_length, score; /* Smith-Waterman score, for swat use; best_score for cross_match */
} Seq_entry;


typedef struct score_entry {
  int seq_entry, score, length;
  float z, E;
} Score_entry;

#define FORCE_REJECT_SCORE -1000  /* NEEDS TO BE CO-ORDINATED WITH MERGING CRITERIA; MUST BE STOREABLE
				   IN SHORT INT */
 
