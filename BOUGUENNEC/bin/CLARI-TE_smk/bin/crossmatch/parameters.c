/*****************************************************************************
#   Copyright (C) 1994-2009 by Phil Green.                          
#   All rights reserved.                           
#                                                                           
#   This software is part of a beta-test version of the swat/cross_match/phrap 
#   package.  It should not be redistributed, or
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

/* NEED TO CHECK THAT DEFAULTS AGREE WITH WHAT PHRAP.DOC SAYS */

  /* TAKE OUT FOLLOWING LINE IF time.h COMPILATION ERROR */
#include <time.h>

#include "swat.h"

#define VERSION "1.090518"

Parameters *parameters; /* exported */ 

FILE *fp_log; /* log file; exported */

static int i_arg, argc, prog_type, i_time;
static char **argv;
static time_t start_time;

get_parameters(t_argc, t_argv, calling_program)
     int t_argc;
     char **t_argv;
     char *calling_program;
{
  int c, gap_flag, n;
  int num_files, swat_flag, phrap_flag, gcphrap_flag, cross_flag, gccross_flag, cluster_flag;
  char *our_alloc();
  FILE *fp;
  FILE *fopenRead();
  File *file;
  int nw(), full_nw(), smith_waterman(), quick_full_smith_waterman(), 
      full_smith_waterman();
  SEQ_AREA sa_test;
  /* remove following if compile error */

  argc = t_argc;
  argv = t_argv;

  for (i_arg = 0; i_arg < argc; i_arg++) printf("%s ", argv[i_arg]);
  fprintf(stderr, "\n");
  for (i_arg = 0; i_arg < argc; i_arg++) fprintf(stderr, "%s ", argv[i_arg]);
  printf("\n%s version %s\n", calling_program, VERSION);
  fprintf(stderr, "\n%s version %s\n", calling_program, VERSION);
  notify("Reading parameters ... ");

  if (sizeof(int) < 4)
    fatalError("Ints are less than 4 bytes; use a different C compiler.");

  parameters = (Parameters *)our_alloc(sizeof(Parameters));
 
  strcpy(parameters->calling_program, calling_program);
  strcpy(parameters->version, VERSION);

  strcpy(parameters->date, "999999:999999");

  /* TAKE OUT FOLLOWING 4 LINES IF time.h COMPILATION ERROR */

  init_times();
  time(&start_time);
  strftime(parameters->date, 20, "%y%m%d:%H%M%S", localtime(&start_time));
  printf("\nRun date:time  %s", parameters->date);
  fprintf(stderr, "Run date:time  %s\n", parameters->date);

  parameters->argc = argc;
  parameters->argv = argv;

  swat_flag = !strcmp(calling_program, "swat");
  phrap_flag = !strcmp(calling_program, "phrap");
  gcphrap_flag = !strcmp(calling_program, "gcphrap");
  cross_flag = !strcmp(calling_program, "cross_match");
  gccross_flag = !strcmp(calling_program, "gccross_match");
  cluster_flag = !strcmp(calling_program, "cluster");

  prog_type = swat_flag ? 4 : cross_flag || gccross_flag ? 2 : phrap_flag || gcphrap_flag ? 1 : 0;

  parameters->query_histograms = swat_flag;
  parameters->query_files = 0;
  parameters->subject_files = 0;

  set_option((float *)0, &parameters->align_extend, "-align_extend", 0.0, 0, 0, 2);  /* # of additional residues to show at ends of displayed alignments */
  set_option((float *)0, &parameters->alignments, "-alignments", 0.0, 0, 1, 2);  /* display alignments */
  set_option((float *)0, &parameters->bandwidth, "-bandwidth", 0.0, 14, 0, 3); /* 1/2 band size for banded swat searches */
  set_option((float *)0, &parameters->bcdsites_qual_threshold, "-bcdsites_qual_threshold", 0.0, 25, 0, 2); /* quality threshold for discrepancy flagging */
  set_option((float *)0, &parameters->bypasslevel, "-bypasslevel", 0.0, 1, 0, 1); /* stringency for bypassing reads */
  set_option((float *)0, &parameters->compact_qual, "-compact_qual", 0.0, 0, 1, 2);  /* store quality values in same bytes as nucs */
  set_option((float *)0, &parameters->confirm_length, "-confirm_length", 0.0, 8, 0, 1);  /* minimum size of confirming segment (starts at first occurrence
									    of 3d distinct nuc following discrepancy) */
  set_option((float *)0, &parameters->confirm_penalty, "-confirm_penalty", 0.0, -5, 0, 1); /* penalty for confirming matches */
  set_option((float *)0, &parameters->confirm_score, "-confirm_score", 0.0, 30, 0, 1); /* minimum SWAT score (using confirm_penalty) for confirming matches */
  set_option((float *)0, &parameters->confirm_trim, "-confirm_trim", 0.0, 1, 0, 1); /* amount by which confirming segments are "trimmed" on
								       either side */
  set_option((float *)0, &parameters->contig_graph_weights, "-contig_graph_weights", 0.0, 0, 0, 1); /* controls weights in weighted directed graph used to find
										       contig sequence -- 0 = quality values, 1 = scaled error probs */
  set_option((float *)0, &parameters->create_test_file, "-create_test_file", 0.0, 0, 0, 3); /* create file for testing validity of qualities */ /* UNDOCUMENTED -- WHICH PROGRAMS?? */
  set_option((float *)0, &parameters->default_qual, "-default_qual", 0.0, 15, 0, 3); /* default quality value to be used when no .qual file is provided */
  set_option((float *)0, &parameters->del_gap_ext, "-del_gap_ext", 0.0, 0, 0, 7); /* default set later */
  set_option((float *)0, &parameters->discrep_lists, "-discrep_lists", 0.0, 0, 1, 2); /* list of discrepancies for each match */
  set_option((float *)0, &parameters->discrep_tables, "-discrep_tables", 0.0, 0, 1, 2); /* table of discrepancies for each match */
  set_option(&parameters->e_cutoff, (int *)0, "-E", 1.0, 0, 0, 4); /* default E value (upper) cutoff */
  set_option((float *)0, &parameters->end_gap, "-end_gap", 0.0, -1, 0, 4); /* gap penalty for gaps at ends of sequences (for Needleman-Wunsch algorithm) */
  set_option((float *)0, &parameters->exp_input, "-exp_input", 0.0, 0, 0, 3); /* flag indicating whether to take experiment files as input */
  set_option((float *)0, &parameters->exp_output, "-exp_output", 0.0, 0, 0, 3); /* flag indicating whether to produce experiment files as output */
  set_option((float *)0, &parameters->file_flag, "-file", 0.0, 0, 1, 4); /* swat only?? */
  set_option((float *)0, &parameters->force_high, "-force_high", 0.0, 0, 1, 1); /*  force isolated merges involving high-quality discrepancies: Allows ignoring high-quality discrepancies during final "forcing" merge pass */
  set_option((float *)0, &parameters->forcelevel, "-forcelevel", 0.0, 0, 0, 1); /* stringency for a final "forcing" merge pass */
  set_option((float *)0, &parameters->fuse_gap1, "-fuse_gap1", 0.0, 0, 1, 3); /* if set, fuse gap1 alignments */
  set_option((float *)0, &parameters->gap_ext, "-gap_ext", 0.0, 0, 0, 7); /* penalty for each succeeding residue in gap */ /* default set later */
  set_option((float *)0, &parameters->gap_init, "-gap_init", 0.0, 0, 0, 7); /* penalty for first residue in gap */ /* default set later */
  set_option((float *)0, &parameters->gap1_dropoff, "-gap1_dropoff", 0.0, -12, 0, 3);
  set_option((float *)0, &parameters->gap1_minscore, "-gap1_minscore", 0.0, 17, 0, 3); /* minscore for extension of word match */
  set_option((float *)0, &parameters->gap1_only, "-gap1_only", 0.0, 0, 1, 2); /* if set, only look for gap1 alignments */
  set_option((float *)0, &parameters->globality, "-globality", 0.0, 0, 0, 3); /* forces semiglobal or global (wrt query) alignments */ /* SHOULD THIS BE CROSS_MATCH ONLY?? */
  set_option((float *)0, &parameters->group_delim, "-group_delim", 0.0, '_', 2, 1); /* delimiter for part of name indicating group of reads to be preassembled */
  set_option((float *)0, &parameters->indexwordsize, "-indexwordsize", 0.0, 0, 0, 3); /* word size for indexing (partial hashing) seqs, in words.c */ /* default set later */
  set_option((float *)0, &parameters->indexwordsize2, "-indexwordsize2", 0.0, 0, 0, 3); /* word size for indexing 2d peaks in merged reads  */ /* default set later */
  set_option((float *)0, &parameters->ins_gap_ext, "-ins_gap_ext", 0.0, 0, 0, 7); /* default set later */
  parameters->logfile = 0;
  set_option((float *)0, &parameters->masklevel, "-masklevel", 0.0, 80, 0, 2); /* maximum percentage of aligned bases in query allowed to
								  overlap a higher scoring match */
  parameters->matrix = 0;
  set_option((float *)0, &parameters->max_group_size, "-max_group_size", 0.0, cross_flag || gccross_flag ? 0 : 10, 0, 3);
  set_option((float *)0, &parameters->max_intron_length, "-max_intron_length", 0.0, 10000, 0, 2);
  set_option((float *)0, &parameters->max_num_alignments, "-max_num_alignments", 0.0, 20, 0, 4); /* max no. displayed alignments */ 
  set_option((float *)0, &parameters->max_num_alignments, "-N", 0.0, 20, 0, 4); /* max no. displayed alignments */ 
  set_option((float *)0, &parameters->max_overlap, "-max_overlap", 0.0, 20, 0, 2); /* for cDNA to genomic alignments */
  set_option((float *)0, &parameters->max_subclone_size, "-max_subclone_size", 0.0, 5000, 0, 1); /* maximum possible size of a subclone - for consistency checks */
  set_option((float *)0, &parameters->maxgap, "-maxgap", 0.0, 30, 0, 1); /* maximum gap size */
  set_option((float *)0, &parameters->maxmatch, "-maxmatch", 0.0, 0, 0, 3); /* maximum length perfect matching word (with no ambiguity chars) */ /* default set later */
  set_option((float *)0, &parameters->min_exon_length, "-min_exon_length", 0.0, 6, 0, 2); /* for cDNA to genomic alignments */
  set_option((float *)0, &parameters->min_intron_length, "-min_intron_length", 0.0, 30, 0, 2); /* for cDNA to genome alignments */
  set_option(&parameters->minmargin, (int *)0, "-minmargin", 0.5, 0, 0, 2); /* only keep matches whose score is at least this amount above every other score */
  set_option((float *)0, &parameters->minmatch, "-minmatch", 0.0, 0, 0, 3); /* minimum length perfect matching word (with no ambiguity chars) */ /* default set later */
  set_option((float *)0, &parameters->minscore, "-minscore", 0.0, 30, 0, 3); /* minimum SWAT score (using PENALTY) for defining read matches */
  set_option((float *)0, &parameters->n_delim, "-n_delim", 0.0, 1, 0, 1); /* for alternative naming conventions in which the subclone name ends at the n-th occurrence of the delimiter */
  set_option((float *)0, &parameters->near_minscore, "-near_minscore", 0.0, 0, 0, 2);
  set_option((float *)0, &parameters->new_ace, "-new_ace", 0.0, 0, 1, 1); /* create new-format .ace file */
  set_option((float *)0, &parameters->new_ace, "-ace", 0.0, 0, 1, 1); 
  set_option((float *)0, &parameters->node_seg, "-node_seg", 0.0, 8, 0, 1);
  set_option((float *)0, &parameters->node_space, "-node_space", 0.0, 4, 0, 1);
  set_option((float *)0, &parameters->nw_flag, "-nw", 0.0, 0, 1, 4); 
  set_option((float *)0, &parameters->old_ace, "-old_ace", 0.0, 0, 1, 1); /* create old-format .ace file */
  set_option((float *)0, &parameters->output_bcdsites, "-output_bcdsites", 0.0, 0, 1, 2); 
  set_option((float *)0, &parameters->output_nonmatching_queries, "-output_nonmatching_queries", 0.0, 0, 1, 3); 
  set_option((float *)0, &parameters->penalty, "-penalty", 0.0, -2, 0, 3); /* mismatch penalty (match is +1, gap is penalty - 2) for primary matching */
  set_option((float *)0, &parameters->preassemble, "-preassemble", 0.0, 0, 1, 1); /* preassemble reads within groups, prior to merging with other groups */
  set_option((float *)0, &parameters->print_extraneous_matches, "-print_extraneous_matches", 0.0, 0, 1, 1); /* in phrap, print non-local matches between contig regions */ 
  set_option((float *)0, &parameters->print_word_data, "-print_word_data", 0.0, 0, 1, 3);
  set_option((float *)0, &parameters->qual_scores, "-qual_scores", 0.0, 0, 1, 2);  /* use quality-based alignment scores */
  set_option((float *)0, &parameters->qual_show, "-qual_show", 0.0, 20, 0, 1); /* LLR cutoff for displaying discrepancies, low-quality regions */
  set_option((float *)0, &parameters->raw, "-raw", 0.0, 0, 1, 7); /* use raw (rather than complexity-adjusted) Smith-Waterman scores */
  set_option((float *)0, &parameters->repeat_screen, "-repeat_screen", 0.0, 0, 0, 3); /* indicates which entries (query vs subject) should be screened for
			annotated repeats, in deciding which word matches to ignore:
			0 neither
			1 first only
			2 second only
			3 both
		     */
  set_option(&parameters->repeat_stringency, (int *)0, "-repeat_stringency", .95, 0, 0, 1); /* assumed frac'l identity of repeats in computing LLR scores;
										     must be < 1.0 -- the closer to 1, the higher the stringency */
  set_option((float *)0, &parameters->retain_duplicates, "-retain_duplicates", 0.0, 0, 1, 1); /*  retain exact duplicate reads, rather than eliminating them  */
  set_option((float *)0, &parameters->revise_greedy, "-revise_greedy", 0.0, 0, 1, 1); /* revise assembly following completion of greedy algorithm */
  set_option((float *)0, &parameters->score_hist, "-score_hist", 0.0, 0, 1, 2);
  set_option((float *)0, &parameters->screen, "-screen", 0.0, 0, 1, 3); /* create copy of first file with matching regions X'd out */
  set_option((float *)0, &parameters->shatter_greedy, "-shatter_greedy", 0.0, 0, 1, 1); /* shatter assembly into likely single-copy pieces */
  set_option((float *)0, &parameters->splice_edge_length, "-splice_edge_length", 0.0, 0, 0, 2); /* length of edge to scan for potential splice signals */
  set_option((float *)0, &parameters->spliced_match_left, "-spliced_match_left", 0.0, 0, 0, 2);
  set_option((float *)0, &parameters->spliced_match_right, "-spliced_match_right", 0.0, 20, 0, 2);
  set_option((float *)0, &parameters->spliced_word_gapsize, "-spliced_word_gapsize", 0.0, 0, 0, 2);
  set_option((float *)0, &parameters->spliced_word_gapsize2, "-spliced_word_gapsize2", 0.0, 0, 0, 2);
  set_option((float *)0, &parameters->subclone_delim, "-subclone_delim", 0.0, '.', 2, 1); /* delimiter for subclone name */
  set_option((float *)0, &parameters->tags, "-tags", 0.0, 0, 1, 3); /* tag selected output lines -- to facilitate output parsing */
  set_option((float *)0, &parameters->trim_penalty, "-trim_penalty", 0.0, -2, 0, 1); /* mismatch penalty for scoring degenerate end sequence */
  set_option((float *)0, &parameters->trim_qual, "-trim_qual", 0.0, 13, 0, 1); /* quality used in trimming reads (modified R. Mott method) */
  set_option((float *)0, &parameters->trim_score, "-trim_score", 0.0, 20, 0, 1); /* minimum score for converting degenerate end sequence to N */
  set_option((float *)0, &parameters->trim_start, "-trim_start", 0.0, 0, 0, 1); /* number of bases to be converted to 'N' at the start of each
								   read */
  set_option((float *)0, &parameters->truncatedb, "-truncatedb", 0.0, 0, 0, 4); 
  set_option((float *)0, &parameters->vector_bound, "-vector_bound", 0.0, phrap_flag || gcphrap_flag || cluster_flag ? 80 : 0, 0, 3); /* no. of bases at beginning of read, matches lying entirely within which should be considered vector */
  set_option((float *)0, &parameters->view, "-view", 0.0, 0, 1, 1);  /* create .view file */
  set_option((float *)0, &parameters->view, "-viewfile", 0.0, 0, 1, 1); 
  set_option((float *)0, &parameters->word_intron_margin, "-word_intron_margin", 0.0, 2, 0, 2); /* # bases to display at each end of spliced-word introns -- for diff list (& alignment display) purposes */
  set_option((float *)0, &parameters->word_offset, "-word_offset", 0.0, 1, 0, 2); /* # bases to offset successive query words that are indexed */

  set_option((float *)0, &parameters->word_raw, "-word_raw", 0.0, 0, 1, 3); /* use raw (rather than complexity-adjusted) word length, for
		    matching words (N.B. maxmatch and indexwordsize always
		    refer to raw word size) */
  set_option(&parameters->z_cutoff, (int *)0, "-z", 6.0, 0, 0, 4); /* default z-value (lower) cutoff */

  parameters->score_flag = parameters->keep_query_data = 0; 
  parameters->align = smith_waterman;
  parameters->full_align = swat_flag ? full_smith_waterman : quick_full_smith_waterman;

  num_files = 0;

  /* N.B. anything not a parameter is assumed to be a filename */
  for (i_arg = 1; i_arg < argc; i_arg++) {
    if (reset_options()) ;
    else if (!strcmp(argv[i_arg], "-exp") || !strcmp(argv[i_arg], "-expfile")) {
      if (!gcphrap_flag) 
	fatalError("option -exp is only available with gcphrap, not phrap -- see documentation");
      if (i_arg >= argc - 1) fatalError("Missing value for -expfile on command line"); 
      parameters->exp_dir = argv[++i_arg];
      parameters->exp_output = 1;
    }
    else if (!strcmp(argv[i_arg], "-M") || !strcmp(argv[i_arg], "-matrix")) { 
      if (i_arg >= argc - 1) fatalError("Missing value for -matrix on command line"); 
      parameters->matrix = argv[++i_arg];
    }
    else if (!strcmp(argv[i_arg], "-logfile")) { 
      if (i_arg >= argc - 1) fatalError("Missing value for -logfile on command line"); 
      parameters->logfile = argv[++i_arg];
    }
    else if (argv[i_arg][0] == '-') {
      fprintf(stderr, "\nFATAL ERROR: Command line option %s not recognized\n", argv[i_arg]);
      exit(1);
    }
    else {
      file = (File *)our_alloc(sizeof(File));
      file->name = argv[i_arg];
      file->fp = fopenRead(file->name); 
      n = strlen(file->name);
      file->type = n > 4 && !strcmp(file->name + n - 5, ".calf") ? 1 : 0;
      if (!file->type) {
	if (fgetc(file->fp) != '>') {
	  if (gcphrap_flag || gccross_flag) 
	    parameters->exp_input = 1;
	  else 
	    reset_file(file);
	}
      }
      else { 
	reset_file(file);
      }
      if (!parameters->query_files) {
	file->next = 0;
	parameters->query_files = file;
	if (file->type) parameters->compact_qual = 1;
      }
      else { /* note order is reversed */
	if (file->type) fatalError("subject files cannot be .calf files");
	file->next = parameters->subject_files;
	parameters->subject_files = file;
      }
    }
  }

  if (!parameters->query_files) 
    fatalError("Sequence files must be specified on command line. See documentation.");

  if ((phrap_flag || gcphrap_flag) && parameters->subject_files)
    fatalError("With current phrap version only one input sequence file may be specified");

  printf("\nQuery file(s): ");

  for (file = parameters->query_files; file; file = file->next)
    printf(" %s", file->name);

  if (parameters->logfile)
    init_log_file(parameters->logfile, 0);
  else 
    init_log_file(parameters->query_files->name, 1);

  if (parameters->subject_files) {
    printf("\nSubject file(s): ");
    for (file = parameters->subject_files; file; file = file->next)
      printf("  %s", file->name);
  }
  else {
    parameters->repeat_screen &= 1; /* won't scan subject files */
    if (is_modified("-minmargin") || parameters->score_hist || parameters->splice_edge_length
	|| parameters->output_bcdsites
	|| parameters->spliced_word_gapsize || parameters->spliced_word_gapsize2 || parameters->qual_scores || parameters->compact_qual)
      fatalError("Minmargin, score_hist, output_bcdsites, spliced_word_gapsize, spliced_word_gapsize2, splice_edge_length, map_qual, and compact_qual can currently only be used when there are 1 or more subject files.");
    if (parameters->masklevel != 101 && (cross_flag || gccross_flag)) {
      parameters->masklevel = 101;
      printf("\nMasklevel 101 required when there is no subject file: resetting");
    }
  }

  /* PRINT OUT FOLLOWING VALS; ALSO INDICATE HOW MINSCORE IS SET IN QUAL_SCORES CASE */
  if (parameters->qual_scores) { 
    printf("\nQuality-based scoring");
    if (parameters->matrix) fatalError("score matrix currently cannot be provided with quality scores");
    parameters->score_hist = 1;
    parameters->compact_qual = 1;
    /* parameters->raw = 1; */
    parameters->gap1_minscore = 100;
    parameters->near_minscore = 100;
    parameters->gap1_dropoff = -60;
    parameters->X_penalty = -6000;
    parameters->match_reward = 6; /* = -10 log_10(1/4) -- for random seq; should be adjusted
for biassed comp, and for low quality bases  */

    parameters->minscore = 0; /* will increase as subject sequences are read */

    parameters->qual_resolution = 5; /* resolution of quality values, for rounding purposes -- needed in qual_scores case to keep size of matrices and
			    (especially!) profile smaller by reducing no. of indep scores
			    e.g. if vals 1-40 are used, then 40 / 5 = 8 indep vals, x 4 nucs = 32

			    Note that should really solve profile issue generically!
			 */ 

    /* N.B. minscore will be reset as length of subject sequence increases */

    /* also need to define appropriate values for: (documentation
     should point out margin issue), minscore, gap_init, gap_ext,
     gap1_minscore ...  anything that is affected by changing the
     matrix */
    parameters->masklevel = 0;
  }
  else {
    parameters->X_penalty = parameters->penalty - 1; /* latter case was 1 -- what is is best? More generally,
							     what about scoring cases not in alphabet */
    parameters->match_reward = 1;
  }
  if (parameters->confirm_penalty > parameters->penalty)
    parameters->confirm_penalty = parameters->penalty;

  if (!parameters->matrix && !is_modified("-penalty")) {
    if (swat_flag) { 
      parameters->matrix = (char *)our_alloc(strlen("BLOSUM50") + 1);
      strcpy(parameters->matrix, "BLOSUM50"); /* default protein score matrix */
    }
  }

  gap_flag = set_score_mat();
  /* FOLLOWING MAY BE INCORRECT FOR QUAL_SCORES CASE */
  if (!is_modified("-gap_init")) parameters->gap_init = gap_flag ? get_gap_init() : -12;
  if (!is_modified("-gap_ext")) parameters->gap_ext = gap_flag ? get_gap_ext() : -2;
  if (!is_modified("-del_gap_ext")) parameters->del_gap_ext = parameters->gap_ext;
  if (!is_modified("-ins_gap_ext")) parameters->ins_gap_ext = parameters->gap_ext;

  set_gap_penalties(parameters->gap_init, parameters->gap_ext, 
		    parameters->ins_gap_ext, parameters->del_gap_ext, parameters->end_gap);


  parameters->DNA_flag = get_DNA_flag(); 
  printf("\nPresumed sequence type (from score matrix): %s", parameters->DNA_flag == 2 ? "merged trace DNA" : parameters->DNA_flag ? "DNA" : "protein");

  if (parameters->compact_qual) { 
    if (parameters->DNA_flag != 1) fatalError("-compact_qual can only be used with ordinary DNA queries");
    parameters->DNA_flag = parameters->query_files->type ? 4 : 3;
    printf("\nCompact quality storage used");
  }

  if (parameters->DNA_flag == 2) {
    printf(" ... forcing -masklevel 101");
    parameters->masklevel = 101;
    if (!(cross_flag || gccross_flag)) fatalError("Merged base DNA only usable with cross_match");
  }

  set_mats((cross_flag || gccross_flag) && (parameters->repeat_screen & 2), parameters->DNA_flag == 2); /* criteria for allowing lower case preservation in
													   subject & query seqs */

  if (parameters->forcelevel < 0) parameters->forcelevel = 0;
  if (parameters->forcelevel > 10) parameters->forcelevel = 10;
  if (parameters->bypasslevel < 0) parameters->bypasslevel = 0;
  if (parameters->bypasslevel > 10) parameters->bypasslevel = 10;

  printf("\n\nPairwise comparison algorithm: %s %s", swat_flag ? "" : parameters->gap1_only ? "gap1" : "banded",
	   parameters->nw_flag ? "Needleman-Wunsch" : "Smith-Waterman");
  /*
  if ((parameters->gap1_only || parameters->spliced_word_gapsize) && (!parameters->DNA_flag || parameters->matrix)) {
      fatalError("\ngap1_only and spliced_word_gapsize currently only usable with DNA sequences, & default penalties");
  }
  */
  parameters->gap1_flag = parameters->gap1_only || parameters->spliced_word_gapsize || parameters->spliced_word_gapsize2;
  if (parameters->globality && !parameters->gap1_flag)
    fatalError("globality can only be set with gap1_only or spliced_word alignments");
  if (parameters->fuse_gap1 && !parameters->gap1_flag)
    fatalError("fuse_gap1 can only be set with gap1_only or spliced_word alignments");
  if (parameters->gap1_only) {
    printf("\nOnly %sgap1 alignments will be found: settings of gap_ext, bandwidth, gap1_minscore, and near_minscore now irrelevant", parameters->fuse_gap1 ? "fused " : ""); 
    printf("\n alignments will be %s (with respect to query)",
	   !parameters->globality ? "local" : parameters->globality == 3 ? "global" : parameters->globality == 1 ? "left-global" : parameters->globality == 2 ? "right-global" : "***");
  }
  if (parameters->spliced_word_gapsize || parameters->spliced_word_gapsize2) {
    if (parameters->spliced_word_gapsize2 && parameters->spliced_word_gapsize2 <= parameters->spliced_word_gapsize)
      parameters->spliced_word_gapsize2 = 0; /* don't need -- since subsumed */
    printf("\nScanning for spliced words flanking gaps >= %d and <= %d (or <= %d near downstream match boundary) , with GT..AG or GC..AG bdries, using gap1 Smith-Waterman.",
	   parameters->min_intron_length, parameters->spliced_word_gapsize, parameters->spliced_word_gapsize2); 
    if (parameters->fuse_gap1) printf(" Alignments will be fused.");
    if (parameters->spliced_word_gapsize2) {
      printf(" Margins around downstream match left boundary: %d to left and %d to right", 
	     parameters->spliced_match_left, parameters->spliced_match_right);
      if (parameters->spliced_match_left < 0 || parameters->spliced_match_right < 0)
	fatalError("-spliced_match_left and -spliced_match_right must be non-negative");
    }
    printf("    word_intron_margin: %d", parameters->word_intron_margin);
    if (parameters->word_intron_margin < 2) 
      fatalError("parameter -word_intron_margin must be at least 2"); /* required for print_diffs */
    printf("\n spliced alignments will be %s (with respect to query)",
	   !parameters->globality ? "local" : parameters->globality == 3 ? "global" : parameters->globality == 1 ? "left-global" : parameters->globality == 2 ? "right-global" : "***");

    set_intron_penalties();
  }

  printf("\n\nScore matrix ");
  if (parameters->matrix) printf("%s", parameters->matrix);
  else if (parameters->qual_scores) printf("set for quality-based scoring");
  else {
    printf("(set by value of penalty: %d)", parameters->penalty);
    if (parameters->penalty >= 0)
      fatalError("Penalty must be negative");
  }
  if (parameters->DNA_flag && parameters->DNA_flag < 3) 
    print_score_mat();

  printf("\n\nGap penalties: gap_init: %d, gap_ext: %d, ins_gap_ext: %d, del_gap_ext: %d, ", 
	 parameters->gap_init, parameters->gap_ext, 
	 parameters->ins_gap_ext, parameters->del_gap_ext); 
  if (parameters->gap_init >= 0)
      fatalError("gap_init must be negative");

  if (swat_flag) {
    parameters->find_z_flag = parameters->raw ? 0 : 1;

    parameters->use_e = is_modified("-E");
    parameters->use_n = is_modified("-max_num_alignments") || is_modified("-N");
    parameters->use_z = is_modified("-z");

    if (parameters->nw_flag) {
      printf(", terminal: %d", parameters->end_gap);
      parameters->find_z_flag = 0; /* no valid E or z-values in this case */
      parameters->align = nw;
      parameters->full_align = full_nw;
    }
    if (!parameters->find_z_flag) {
      printf("\nz-scores and E-values will not be computed.");
      parameters->use_n = 1;
      parameters->use_z = parameters->use_e = 0;
    }
    else if (!(parameters->use_e || parameters->use_z || parameters->use_n)) 
      parameters->use_e = parameters->use_z = parameters->use_n = 1;
    printf("\n");

    if (parameters->use_z) printf("z cutoff: %.4g   ",  parameters->z_cutoff);
    if (parameters->use_e) printf("E cutoff: %.4g   ",  parameters->e_cutoff);
    if (parameters->use_n) printf("Max. no. displayed alignments: %d", parameters->max_num_alignments);  
    printf("\n");
  }
  if (phrap_flag || gcphrap_flag || cross_flag || gccross_flag || cluster_flag) {

    if (!is_modified("-indexwordsize")) 
      parameters->indexwordsize = parameters->DNA_flag ? 12 : 4;
    if (!is_modified("-indexwordsize2")) 
      parameters->indexwordsize2 = 4;
  
    if (!is_modified("-minmatch")) 
      parameters->minmatch = parameters->DNA_flag ? 14 : 4;
    
    if (!is_modified("-maxmatch")) {
      /*
	parameters->maxmatch = (cross_flag || gccross_flag) ? parameters->minmatch :
	(parameters->DNA_flag ? DNA_MAXMATCH : PROT_MAXMATCH);
      */
      parameters->maxmatch = parameters->DNA_flag ? 20 : 4;
    }

    if (parameters->maxmatch > 127) parameters->maxmatch = 127;
    if (parameters->minmatch > 127) parameters->minmatch = 127;

    /* necessary because of use of single byte (less 1 bit) to represent
       match length, in words.c */

    if (!parameters->near_minscore) 
      parameters->near_minscore = parameters->minscore;

    if (parameters->maxmatch < parameters->minmatch) 
      parameters->maxmatch = parameters->minmatch;
    if (parameters->indexwordsize > parameters->minmatch) 
      parameters->indexwordsize = parameters->minmatch;

    if (parameters->DNA_flag && parameters->indexwordsize > 4 * sizeof(int))
      parameters->indexwordsize = sizeof(int) * 4;
    /* necessary because of method used to compute indices rapidly in words.c
       (using the variable mask) */

    if (parameters->indexwordsize2 > parameters->indexwordsize) 
      parameters->indexwordsize2 = parameters->indexwordsize;

    if (!parameters->raw && parameters->DNA_flag) {
      printf("\nUsing complexity-adjusted scores.");
      print_background_freqs(parameters->compact_qual ? 0 : 1);
      set_complexity_adjust(.25, parameters->compact_qual ? 0 : 1); /* perhaps should always be 0? */
    }
    else printf("\nUsing raw scores.");
    if (parameters->maxmatch == parameters->minmatch) {
      printf("\nSince minmatch = maxmatch, setting max_group_size to 0 and turning off complexity-adjustment of word lengths");
      parameters->max_group_size = 0;
      parameters->word_raw = 1;
    }

    if (parameters->repeat_screen) {
      printf("\n\nWill ignore word matches lying within %s%s", 
	     parameters->repeat_screen & 1 ? "query-file-tagged repeats or " : "",
	     parameters->repeat_screen & 2 ? "subject file lower-case regions" : "");
    }
    parameters->min_record_score = parameters->minscore;

    if (!is_modified("-word_offset")) 
      parameters->word_offset = (parameters->DNA_flag == 1 || parameters->DNA_flag >= 3) && parameters->subject_files && !parameters->spliced_word_gapsize && !parameters->spliced_word_gapsize2 ? 2 : 1;

    if (parameters->word_offset < 1) fatalError("-word_offset must be at least 1");

    if (!parameters->subject_files && parameters->word_offset != 1) fatalError("with no subject files -word_offset must = 1");

    printf("\n\nminmatch: %d, maxmatch: %d, max_group_size: %d%s, minscore: %d, near_minscore: %d, bandwidth: %d, indexwordsize: %d, indexwordsize2: %d, word_offset: %d", 
	   parameters->minmatch, parameters->maxmatch, 
	   parameters->max_group_size, parameters->max_group_size > 0 ? "" : " (turned off)",
	   parameters->minscore, parameters->near_minscore, parameters->bandwidth,
	   parameters->indexwordsize, parameters->indexwordsize2, parameters->word_offset);
    printf("\nword_raw: %d", parameters->word_raw);
    printf("\nvector_bound: %d", parameters->vector_bound); 
    if (!parameters->qual_scores && parameters->gap1_minscore > parameters->minscore) {
      printf("\nGap1_minscore is greater than minscore: resetting");
      parameters->gap1_minscore = parameters->minscore;
    }
    printf("\ngap1_minscore: %d, gap1_dropoff: %d", parameters->gap1_minscore, parameters->gap1_dropoff);
    if (parameters->gap1_dropoff >= 0) fatalError("gap1_dropoff must be negative");
  }    
  if (phrap_flag || gcphrap_flag) {
    if (parameters->trim_start) 
      printf("\n%d characters trimmed from start of each read", parameters->trim_start);
    printf("\ntrim_penalty: %d, trim_score: %d, trim_qual: %d, maxgap: %d", 
	   parameters->trim_penalty, parameters->trim_score, parameters->trim_qual, parameters->maxgap);
    printf("\nrepeat_stringency: %f", parameters->repeat_stringency);
    if (parameters->repeat_stringency <= 0 || parameters->repeat_stringency >= 1)
      fatalError("repeat_stringency must be > 0 and < 1");
    printf("\nqual_show: %d", parameters->qual_show);
    
    printf("\nconfirm_length: %d, confirm_trim: %d, confirm_penalty: %d, confirm_score: %d",
	   parameters->confirm_length, parameters->confirm_trim, parameters->confirm_penalty, 
	   parameters->confirm_score);
    printf("\nnode_seg: %d, node_space: %d", 
	   parameters->node_seg, parameters->node_space); 
    printf("\nforcelevel: %d, bypasslevel: %d", parameters->forcelevel, parameters->bypasslevel);
    printf("\ncontig_graph_weights: %d", parameters->contig_graph_weights);
    printf("\nmax_subclone_size: %d", parameters->max_subclone_size);
  }
    
  if (cross_flag || gccross_flag) {
    /* cases:
        masklevel 101:  -minmargin irrelevant; score_hist needed only if set, in which case single query_domain covering entire length
                                   
                  < 101: -minmargin relevant  
                  0:     single query domain covering entire length
        minmargin <= 0: only need highest score for each qd
                  = 0.5: need highest score and # matches at that score (each qd)
                  = 1:          "                  "
                  >= 2: need partial score hist: scores within (minmargin - 1) of best; may as well keep all

	score_flag settings:
                  0: no need to check, or keep query domain (only applies to masklevel 101 with no score_hist)
                  1: update query domain (max score, n_max)
                  >1: update score_hist
    */


    

    printf("\nmasklevel: %d", parameters->masklevel);

    if (parameters->masklevel == 101) {
      printf(" (minmargin irrelevant)");
    }
    else if (parameters->masklevel < 0 || parameters->masklevel > 101)
      fatalError("Masklevel must be >= 0 and <= 101");
    else { /* FOLLOWING INCORRECT IN QUAL_SCORES CASE */
      printf("\nminmargin: %.1f", parameters->minmargin); 
      if (parameters->minmargin > 1) { 
	parameters->min_record_score = parameters->minscore - (int)parameters->minmargin + 1;
	if (parameters->min_record_score < 1) parameters->min_record_score = 1; 
      } 
    }

    if (parameters->gap1_flag || parameters->masklevel != 101) {
      parameters->score_flag |= 1; /* must turn on scorekeeping in gap1 case to avoid dup hits */
    }

    if (parameters->score_hist) {
      parameters->score_flag |= 2;
      printf("\nScore histograms printed: score threshold %d", parameters->min_record_score);
    }

    parameters->splice_edge_alloc = 8 * parameters->splice_edge_length + 1;
    printf("\nsplice_edge_length: %d, allocation: %d bytes", parameters->splice_edge_length, parameters->splice_edge_alloc);
    printf("\nmin_intron_length: %d, max_intron_length: %d, max_overlap: %d, min_exon_length: %d ", 
	   parameters->min_intron_length, parameters->max_intron_length, parameters->max_overlap, parameters->min_exon_length);
    if (parameters->output_bcdsites) {
      if (!parameters->discrep_lists) {
	printf("\n\nWarning: setting -discrep_lists (required for confirmed & discrepant subject segments)");
	parameters->discrep_lists = 1;
      }
      parameters->score_flag |= 1; /* need query domain accounting */
      printf("\n\nConfirmed & discrepant subject sites are output in files %s.bcdsites. Bcdsite quality threshold for discrepancies: %d", 
	     parameters->query_files->name, parameters->bcdsites_qual_threshold);
    }
    if (parameters->splice_edge_length || parameters->score_flag) {
      parameters->keep_query_data = 1;
    }
  }

  if (parameters->output_nonmatching_queries)
    printf("\n\nNonmatching queries & qualities are output in files %s.nonmatching, %s.nonmatching.qual", 
	   parameters->query_files->name, parameters->query_files->name);

  sa_test = -1;
  printf("\nVariable sizes: int: %d, long int: %d, SEQ_AREA: %d, signed: %d\n", (int)sizeof(int), (int)sizeof(long int), (int)sizeof(SEQ_AREA), sa_test < 0);

  notify(" Done\n");
}

typedef struct option {
  int *int_add;
  float *float_add;
  char name[30];
  int flag, reset, prog_scope;
} Option;

Option options[100];

int n_options;

/* flag: 0 if ordinary (numerical) param, 1 if flag, 2 if character-valued param; 
   prog_scope: 1 phrap
               2 cross_match
	       4 swat
       or sum of relevant values
*/

set_option(float_add, int_add, name, float_value, int_value, flag, prog_scope)
     int *int_add, int_value;
     float *float_add, float_value;
     char *name;
     int flag, prog_scope;
{
  Option *option;

  option = options + n_options++;
  if (n_options >= 100) fatalError("too many options");

  option->int_add = int_add;
  option->float_add = float_add;
  if (int_add) *int_add = int_value;
  else *float_add = float_value;
  if (strlen(name) >= 30) fatalError("option name too long");
  strcpy(option->name, name);
  option->flag = flag;
  option->reset = 0;
  option->prog_scope = prog_scope;
}

reset_options()
{
  Option *option;
  int i_option;
  double x;
  char c;

  for (i_option = 0; i_option < n_options; i_option++) {
    option = options + i_option;
    if (!strcmp(option->name, argv[i_arg])) {
      if (!(prog_type & option->prog_scope)) {
	fprintf(stderr, "\nFATAL ERROR: Option %s cannot be used with %s\n", 
		argv[i_arg], parameters->calling_program);
	exit(1);
      }
      option->reset = 1;
      if (option->flag == 1) *option->int_add = 1;
      else {
	if (i_arg >= argc - 1 
	    || !option->flag && !(sscanf(argv[i_arg + 1], "%lf", &x))
	    || option->flag == 2 && !(sscanf(argv[i_arg + 1], "%c", &c))

	    ) {
	  fprintf(stderr, "\nFATAL ERROR: Missing or misspecified value for %s on command line\n", 
		  argv[i_arg]);
	  exit(1);
	}
	if (option->int_add) *option->int_add = option->flag == 2 ? c : x;
	else *option->float_add = x;
	i_arg++;
      }
      return 1;
    }
  }
  return 0;
}

is_modified(option_name)
     char *option_name;
{
  Option *option;
  int i_option;
  char *name;

  for (i_option = 0; i_option < n_options; i_option++) {
    option = options + i_option;
    if (!strcmp(option->name, option_name)) {
      return option->reset;
    }
  }
  fatalError("option not found");
}

static time_t times[10];
static char **time_labels;

init_times()
{
  char *our_alloc();
  int i;

  time_labels = (char **)our_alloc(10 * sizeof(char *));

  for (i = 0; i < 10; i++) 
    time_labels[i] = (char *)our_alloc(20 * sizeof(char));

  set_time("start");
}

set_time(time_label)
     char *time_label;
{
  strcpy(time_labels[i_time], time_label);
  time(&times[i_time]);
  
  i_time++;
  if (i_time >= 10) fatalError("timing slots exceeded");
}

print_times()
{
  int i;
  double x, cum, difftime();

  printf("\n\nTimes in secs (cum)\n");
  for (i = 1, cum = 0; i < i_time; i++) {
    x = difftime(times[i], times[i - 1]);
    cum += x;
    printf("%12s %3.0f (%3.0f)\n", time_labels[i], x, cum);
  }
}

init_log_file(file, append_flag)
     char *file;
     int append_flag;
{
  char *our_alloc();
  char *log_file;
  FILE *fopenWrite();
  
  if (append_flag) {
    log_file = (char *)our_alloc(strlen(file) + 10);
    strcpy(log_file, file);
    strcat(log_file,".log");
  }
  else log_file = file;

  fp_log = fopenWrite(log_file);
  printf("   Log file: %s", log_file);
  if (append_flag) our_free(log_file);
}

static double total_subject_length;

add_subject_length(length)
     int length;
{
  total_subject_length += length;
  if (parameters->qual_scores) {
    parameters->min_record_score = parameters->minscore = 10.0 * log10(2 * total_subject_length); /* factor of 2 assumes DNA */
  }
}

double get_subject_length()
{
  return total_subject_length;
}
