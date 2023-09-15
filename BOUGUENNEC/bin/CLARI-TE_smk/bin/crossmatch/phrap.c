/*****************************************************************************
#   Copyright (C) 1994-1999 by Phil Green.                          
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

extern Parameters *parameters; /* from get_parameters() */
extern Database *query_db;

main(argc,argv)
     int argc;
     char *argv[];
{
/* get parameter values -- from defaults / command line */

#if defined(GCPHRAP)
  get_parameters(argc, argv, "gcphrap");
#else
  get_parameters(argc, argv, "phrap");
#endif

/*
  init_qual_arrays();
*/
/* read in sequence files, initialize array of aligned read info, and log file of ancillary information,
 and find pairwise matches */
  readin_and_match();

/* find confirmed starts and ends (using all matching reads that are not both same subclone and same chemistry */
  find_starts_ends(query_db, 0);

  trim_quals(query_db, 0);
  trim_X_quals(query_db);

  print_t_n_pairs();

  if (!parameters->subject_files) {

/* find parts of reads that match other reads (from different subclone or different chemistry) */
    find_segments(query_db);

/* Find additional likely vector parts of reads: occurring within 1st
   parameters->vector_bound bases, matching only similarly located parts
   of same-orientation reads. Set first_vec, last_vec.  Unlike
   elim_complete_vector_matches() now allow partial matches.  Done here
   (after make_reversed_pairs) in order to get more complete list of
   overlapping unremoved pairs. */

    find_more_vector(query_db); 

    /* print probable vector involved in matches. Also set last_vec equal to 
       min(first non-vector match, parameters->vector_bound, first marked base,
       extended by any X's) */

    print_vector(query_db); 
    
  }

/*  slide_indels(); */

/* find average offset of one read with respect to another, as implied
by a pairwise match */

  find_mean_offsets();

  if (!parameters->subject_files) {

/* print near duplicate matching reads that are not perfect
duplicates. No special processing is given to these */

    find_near_duplicates();

/* print matches of reads to themselves (i.e. short tandem
duplications); mark corresponding pair as "reject_self" */

    find_self_matches();

/* find and mark pairs for which there are no acceptable cross-linked
nodes (i.e. strongly aligned regions). These are not allowed in
joining contigs, because may result in inability to resolve layout */

    find_node_rejects(query_db); 

/* mark pairs for which alignment does not extend as far to left or to
right as it apparently should (this no longer used? In any case, needs
to be corrected, for vector) */

    find_truncated_pairs(); 

/* find multi-segment reads (which include potential chimeras), choose
favored segment */

    find_chimeras(query_db);

/* find confirmed starts and ends (of non-chimeric pairs) */
    find_starts_ends(query_db, 1);

/* find deletion reads */
    find_deletions();

/* find confirmed segments; compute adjusted quality. Affects:
adj_qual, first_start, last_end, rev_first_start, rev_last_end */

    find_extents(query_db, 0);

/*
    trim_quals(query_db, 1);
*/    
    print_quality("Revised quality", query_db);

/* initialize variables and arrays used in computing LLR scores */
    init_qual_arrays();

/*
  set_qual_arrays(query_db);
  set_LLR_scores();
  set_qual_arrays(query_db);
*/

/* find LLR score for each pair */
    set_LLR_scores();

/* 2d pass -- now using positive LLR-scoring pairs 
*/
    find_extents(query_db, 1);
/*
    trim_quals(query_db, 1);
*/
/*
  set_qual_arrays(query_db);
*/
  }
  else {
    init_adj_qual(query_db); /* CHECK THIS */
  }

/* compute LLR scores for each pair */
  set_LLR_scores();

/* in cases where there are several overlapping alignments involving
the same two reads, mark the one with the highest LLR score */

  find_best();

/*  find_rejects();  */

/*  find_triple_rejects(); */
/*  print_reject_summary(); perhaps restore this later */

/*  debug_alignments(); */

  if (!parameters->subject_files) {
/* histogram of adjusted quality */

    print_quality("2d revised quality", query_db);

    print_coverage(query_db);
    analyze_discreps(query_db);
  }

/* find blocked reads (to avoid in first merge passes) */
  find_blocked_reads(query_db);

/* make read layout, using greedy algorithm */
  
  merge_master();
  merge_chimeras();
  merge_other_singletons();
  free_LLR_pair_pointers();
/* sort contigs by no. reads, sort reads within contigs, etc. */
  clean_contigs();

 /* 3d pass -- now using pairs that are either used in merges, AND are
positive LLR-scoring */

  find_extents(query_db, 2); 

  trim_quals(query_db, 1);

/* make contig sequences */
  make_contig_db();

/* find best pairwise alignment of each read to contig */
  align_reads_to_contigs(0);

/* find local improvements to sequence -- short higher scoring regions
from particular reads, better resolved compressions */

  revise_contigs();

/* find best pairwise alignment of each read to contig */

  align_reads_to_contigs(1);

/* flag pairs which involve part of read not aligned to contig */
  find_unaligned_pairs();

  print_singletons();
  print_t_n_segments();
  print_contigs();
  print_t_n_segments();

  write_contigs(); /* done here rather than earlier, because quality values have 
		      been changed in print_contigs */

  map_reads_by_name(query_db);
  check_pairs();
  check_alloc_list();
  print_t_n_pairs();
  print_t_n_segments();
  print_t_n_diffs();
  print_t_n_tags();
  notify("\n");
  return 0;
}

/* temporary code to create assembly database, with 1% deletion, 1% insertion, and
   1% substitution errors 
  fq = fopenWrite("zk637.assembly.err.2");
  db_entry = read_db;
  get_next_entry(db_entry, fp);
  for (i = 0; i < 600; i++) {
    j = random() % db_entry->length;
    fprintf(fq,">read%d\n", j);
    if (random() % 2) {
      for (k = j; k < j + 350 && k < db_entry->length; ) {
	if (random() % 100) 
	  fprintf(fq,"%c", random() % 100 ? db_entry->seq[k] : "ACGT"[random() % 4]);
	if (random() % 100) k++; 
      }
    }
    else {
      for (k = j; k > j - 350 && k >= 0; ) {
	if (random() % 100) 
	  fprintf(fq,"%c", random() % 100 ? c_mat[db_entry->seq[k]]
	                                 : "ACGT"[random() % 4]);
	if (random() % 100) k--;
      }
    }
      
    fprintf(fq,"\n");
  }
  exit();
*/


