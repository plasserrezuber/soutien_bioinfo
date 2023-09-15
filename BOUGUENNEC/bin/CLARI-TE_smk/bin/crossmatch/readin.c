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

extern Parameters *parameters; /* from get_parameters() */

Database *query_db; 

/* read in sequence files, initialize array of aligned read info, and log file of ancillary information,
 and find pairwise matches */

readin_and_match()
{
  File *file;
  Database *append_db(), *call_append_db();
  int phrap_flag, gcphrap_flag, subject_flag, cluster_flag;
  /*
     fprintf(stderr, "\nStructure sizes (bytes): Segment %d, Db_entry %d, Align_info %d, Contig %d, Aligned_pair %d, Cand_pair %d, Entry_index %d\n", 
	     sizeof(Segment), sizeof(Db_entry), sizeof(Align_info), sizeof(Contig),
	     sizeof(Aligned_pair), sizeof(Cand_pair), sizeof(Entry_index));
	     */
  phrap_flag = !strcmp(parameters->calling_program, "phrap");
  gcphrap_flag = !strcmp(parameters->calling_program, "gcphrap");
  subject_flag = parameters->subject_files ? 1 : 0;
  cluster_flag = !strcmp(parameters->calling_program, "cluster");
  
  query_db = call_append_db(parameters->query_files, 1, parameters->DNA_flag > 0, 1); /* stores complements -- for speed;
										may want to use 1, 0, 1 for large
										assemblies */

  if (phrap_flag || gcphrap_flag) {
    name_check(query_db);
    set_templates(query_db); /* if experiment files are used, templates,
				  chemistries, directions should be set there */
    set_chemistries(query_db); 
    set_directions(query_db); 
  }
  if (parameters->repeat_screen & 1) set_repeat_tags(query_db); 
  /* may not want to do the following for cross_match */
  /* read quality files */
  read_qual(query_db);
  set_scaled_err_probs();
  set_qual_segs(query_db);
  set_time("read queries");
  for (file = parameters->subject_files; file; file = file->next) {
    if (phrap_flag) {
      append_db(file, 1, parameters->DNA_flag, 1);
      read_qual(file->db);
      set_qual_segs(query_db); /* WHY DUPLICATED?? */
      name_check(file->db);
      set_templates(file->db);
      set_chemistries(file->db);
      set_directions(file->db);
      set_repeat_tags(query_db); /* WHY DUPLICATED?? */
    }
    else append_db(file, 0, 0, 0);
  }
  
  /* reset quality at ends of reads -- require 7 pos qualities out of 10. This should be
     done in phred, instead */ 
  /*  trim_qual(query_db, 10, 7); No longer done */
  
  /* convert to N's the first trim_start bases at beginning of read, or
     regions at either end of read having quality <= 20 and a score of
     trim_score or greater when aligned against a mononucleotide sequence
     (using trim_penalty). Main effect is to reduce number of word matches
     due to mononucleotide runs involving spurious sequence, since use of
     complexity-adjusted Smith-Waterman scores tends to eliminate spurious
     swat matches */
  
  if (phrap_flag || cluster_flag) convert_ends(query_db);
  
  if (parameters->minmatch) {	/* also put this in phrap ? */

    /* find pairs of reads having matching words (not containing N or X)
       of length minmatch or greater.  Each such word match defines a band in
       the Smith-waterman matrix, centered on the diagonal containing the
       word match and of width 2 * bandwidth + 1. Overlapping bands (due to
       multiple word matches) are merged. */


    set_word_db(query_db, subject_flag && parameters->DNA_flag);
    set_time("sort queries");
    /* following for testing only 
    print_times();
    exit(1);
    /* */
    /*
       sort_words(query_db->t_length + query_db->lengths[0], parameters->minmatch, "NX");
       */
    if (subject_flag) {
      set_splice_params();
      find_subject_matches(parameters->subject_files);
      free_word_arrays();
    }
    else {
      if (parameters->DNA_flag /* && !cluster_flag */) find_comp_matches();
      new_find_internal_word_matches();
      /* find reads that are exact duplicates of other reads, mark them, and remove 
	 their cand_pairs; they are excluded from all subsequent analyses (in phrap).  */
      free_word_arrays();
      find_duplicates(query_db);
      if (cluster_flag) return;
      find_all_scores();
      free_cand_pair_blocks();
      free_seg_blocks();	/* assumes no segments except those in cand pairs have been allocated */
    }
  }
  else {
    fatalError("Currently require positive word size");
    make_full_pairs(query_db);
    find_all_scores();
    free_cand_pair_blocks();
    free_seg_blocks();		/* assumes no segments except those in cand pairs have been allocated */

  }

/* revise scores */
  if (parameters->qual_scores) {
    revise_scores(query_db);     
  }

  set_time("find matches");
  print_n_swats();
  
  /* reactivate later?
     test_self();
     */
  
  /* Do recursive banded swat search of bands in swat matrix that contain matching words */
  
  /* find matches that occur entirely within first parameters->vector_bound 
     bases of each read (with both reads in same orientation), and eliminate.  
     Also done in words.c */
  
  if (phrap_flag && !subject_flag) 
    elim_complete_vector_matches(query_db); 

  if (parameters->score_flag) 
    merge_qds(query_db, 0);

  if (parameters->fuse_gap1) {
    fuse_all_pairs();
    if (parameters->score_flag) 
      merge_qds(query_db, 1);
  }

  /* make pair information for reversed pairs (i.e. with order of entries reversed) */



  make_reversed_pairs();

  if (parameters->output_nonmatching_queries)
    output_nonmatching_queries();
/*  slide_indels();  may not really be necessary */
}

