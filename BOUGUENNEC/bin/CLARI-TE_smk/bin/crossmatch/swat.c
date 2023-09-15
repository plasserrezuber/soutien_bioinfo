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

extern Parameters *parameters; 
        
compare_scores(entry1, entry2)
     Score_entry *entry1, *entry2;
{
  return entry2->score - entry1->score;
}

compare_z_scores(entry1, entry2)
     Score_entry *entry1, *entry2;
{
  if (entry2->z == entry1->z) return 0;
  else return entry2->z > entry1->z ? 1 : -1;
}

compare_E_values(entry1, entry2)
     Score_entry *entry1, *entry2;
{
  if (entry2->E == entry1->E) return 0;
  else return entry2->E < entry1->E ? 1 : -1;
}

main(argc,argv)
     int argc;
     char *argv[];
{
  Profile *make_profile_from_seq();
  Profile *q_profile;
  FILE *fr;
  FILE *fopenWrite();
  double find_score_Evalue(), new_find_Evalue();
  double z_to_E();
  double e_value, zcut, e, prev_zcut, prev_e;
  char nameBuffer[100];
  Database *qdb, *sdb; 
  Database *append_db(); 
  Score_entry *score_entries, *score_entry, *last_score_entry;
  char *our_alloc();
  char *get_id(), *get_descrip();
  unsigned char *get_seq();
  Seq_entry *get_seq_entry();
  int print_flag, sdb_pop_flag;
  int i, s_entry;
  int nz, n_trimmed;
  int q_length, score, length, seq_num, orig_score, n_pos_scores, success;
  double t_cells;
  int max_score_cutoff;

  get_parameters(argc, argv, "swat"); 

  qdb = append_db(parameters->query_files, 0, 0, 0);

/*  (instead loop over subject_files) */
  sdb = append_db(parameters->subject_files, 0, 0, 0); 

  t_cells = 0;

  sdb_pop_flag = 0;
  alloc_hist();
  while (get_next_file_entry(qdb)) {
    printf("\n\n*********************************************************************************\n");
    q_length = qdb->seq_length;
    printf("\nQuery: %s  %s   Length: %d residues\n",
	   qdb->id_buffer, qdb->descrip_buffer, q_length);
    q_profile = make_profile_from_seq((Profile *)0, qdb->seq_buffer, q_length, parameters->nw_flag);

    max_score_cutoff = (q_profile->max_score_cutoff - q_profile->poly_offset) & 255;

    notify("Searching: ");
    if (!sdb_pop_flag) {
      while (get_next_file_entry(sdb)) {
	s_entry = append_seq_entry(sdb);
#if defined(ALPHA)
	score = fast_smith_waterman(q_profile, sdb->seq_buffer, sdb->seq_length, 1, 0, 0, 0, 0, 0);
	if (score >= max_score_cutoff)
	  score = smith_waterman(q_profile, sdb->seq_buffer, sdb->seq_length, 1, 0, 0, 0, 0, 0);

#elif defined(COUNTS)
	score = smith_waterman_counts(q_profile, sdb->seq_buffer, sdb->seq_length, 1, 0, 0, 0, 0, 0);
#else 
	score = parameters->align(q_profile, sdb->seq_buffer, sdb->seq_length, 1, 0, 0, 0, 0, 0);
#endif
	t_cells += q_length * sdb->seq_length;
	set_score(s_entry, score);
	if (!(sdb->num_entries % 1000)) notify(".");
	if (sdb->num_entries == parameters->truncatedb) break;  
      }
      score_entries = (Score_entry *)our_alloc(sdb->num_entries * sizeof(Score_entry));
      last_score_entry = score_entries + sdb->num_entries;
      sdb_pop_flag = 1;
    }
    else {
      for (s_entry = sdb->first_entry; s_entry <= sdb->last_entry; s_entry++) {
#if defined(ALPHA)
	score = fast_smith_waterman(q_profile, get_seq(s_entry), get_seq_length(s_entry), 1, 0, 0, 0, 0, 0);
	if (score >= max_score_cutoff)
	  score = smith_waterman(q_profile, get_seq(s_entry), get_seq_length(s_entry), 1, 0, 0, 0, 0, 0);
#elif defined(COUNTS)
	score = smith_waterman_counts(q_profile, get_seq(s_entry), get_seq_length(s_entry), 1, 0, 0, 0, 0, 0);
#else 
	score = parameters->align(q_profile, get_seq(s_entry), get_seq_length(s_entry), 1, 0, 0, 0, 0, 0);
#endif
	t_cells += q_length * get_seq_length(s_entry);
	set_score(s_entry, score);
	if (!((s_entry - sdb->first_entry) % 1000)) notify(".");
      }
    }

    fprintf(stderr, "\nDone. %d entries processed. Total no. cells: %.0f\n", 
	    sdb->num_entries, t_cells);
#if defined(ALPHA)
    print_cell_counts();
#endif
#if defined(COUNTS)
    print_counts();
#endif
    fprintf(stdout, "\n%d entries processed. Total no. cells: %.0f\n", 
	    sdb->num_entries, t_cells);
    fflush(stderr);
    initialize_hist(); 
    n_pos_scores = 0;
    for (s_entry = sdb->first_entry, score_entry = score_entries; 
	 s_entry <= sdb->last_entry; s_entry++, score_entry++) {
      score = get_seq_entry(s_entry)->score;
      length = get_seq_length(s_entry);
      score_entry->seq_entry = s_entry;
      score_entry->score = score;
      score_entry->length = length;
      update_hist(score_entry, 0);
      if (score) n_pos_scores++;
    }
    process_hist();

    /* sort database entries by (decreasing) score */
/*    if (sdb->num_entries > 10) est_lambda_K(q_length);
      No longer use lambda and K 
*/
    if (n_pos_scores < 100 && parameters->find_z_flag) {  /* no valid E or z-values in this case */
      printf("\nToo few positive scores to compute valid z-scores and E-values");
      parameters->find_z_flag = 0;
      parameters->use_n = 1;
      parameters->use_z = parameters->use_e = 0;
    }
    if (parameters->find_z_flag) {  
      for (i = 1; i <= 2; i++) {
	fit_log_n(q_length);
	initialize_hist();
	for (score_entry = score_entries; score_entry < last_score_entry; score_entry++) { 
	  find_z(score_entry);
	  update_hist(score_entry, 1);
	}
      }

      fit_log_n(q_length);
      for (score_entry = score_entries; score_entry < last_score_entry; score_entry++) {
	find_z(score_entry);
      }

      qsort((char *)score_entries, sdb->num_entries, sizeof(Score_entry), compare_z_scores);

      printf("\n\nObserved and expected z distribution (non-rejected entries):");
      printf("\n z range        Obs.     Exp.    Cum Exp.");
      prev_zcut = 99.0;
      prev_e = 0.0;
      score_entry = score_entries;
      for (zcut = 8.0; score_entry < last_score_entry; zcut -= 1.0) {
	for (nz = 0; score_entry < last_score_entry && score_entry->z >= zcut;
	     score_entry++)
	  if (!reject_entry(score_entry)) nz++; /* assumes Smith-Waterman */
	e = z_to_E(zcut);
	printf("\n%4.1f - %4.1f   %5d  %7.1f  %7.1f", zcut, prev_zcut, nz, e - prev_e, e);
	prev_e = e;
	prev_zcut = zcut;
      }

      new_est_lambda_K(q_length, score_entries, last_score_entry);
      for (score_entry = score_entries; score_entry < last_score_entry; score_entry++) {
	score_entry->E = new_find_Evalue(score_entry->score, q_profile->length, score_entry->length);
      }
      qsort((char *)score_entries, sdb->num_entries, sizeof(Score_entry), compare_E_values);

      printf("\n\nObserved and expected E distribution (non-rejected entries):");
      printf("\n     E range            Obs.     Exp.    Cum Exp.");
      prev_zcut = 0.0;
      prev_e = 0.0;
      score_entry = score_entries;
      for (zcut = 1.0; score_entry < last_score_entry; zcut *= 10.0) {
	for (nz = 0; score_entry < last_score_entry && score_entry->E <= zcut;
	     score_entry++)
	  if (!reject_entry(score_entry)) nz++; /* assumes Smith-Waterman */
	if (score_entry == last_score_entry) zcut = (score_entry - 1)->E;
	e = zcut;
	printf("\n%8.1f -%8.1f   %5d  %8.1f  %8.1f", prev_zcut, zcut, nz, e - prev_e, e);
	prev_e = e;
	prev_zcut = zcut;
      }
    }
    else {
      qsort((char *)score_entries, sdb->num_entries, sizeof(Score_entry), compare_scores);
    }

/* find alignments for all scores above cutoff */
    if (parameters->file_flag) {
      nameBuffer[0] = 0;
      strcpy(nameBuffer, qdb->id_buffer);
      strcat(nameBuffer,".allscores");
      fr = fopenWrite(nameBuffer);
    }
    printf("\n\n\nMatches ranked by decreasing %s:", 
	   parameters->find_z_flag ? "z-scores" : "raw alignment scores");

    print_flag = 1;
    for (score_entry = score_entries; score_entry < last_score_entry; score_entry++) {
      seq_num = score_entry->seq_entry;
      if (print_flag) {
	if (parameters->use_n && score_entry >= score_entries + parameters->max_num_alignments
	 || parameters->use_z && score_entry->z < parameters->z_cutoff || !score_entry->score) 
	  print_flag = 0;
	else if (parameters->find_z_flag) {
	  e_value = z_to_E(score_entry->z); 
/* find_score_Evalue(score_entry->score, q_length); */
	  if (parameters->use_e && e_value > parameters->e_cutoff) 
	    print_flag = 0;
	}
      }
      if (print_flag) {
	printf(parameters->find_z_flag ? "\n\n\n%s %s  Length: %d\n  Score: %d  z: %.2f  E: %.3g  new_E: %.3g"
	       : "\n\n%s %s  Length: %d\n  Score: %d",
	       get_id(seq_num), get_descrip(seq_num), get_seq_length(seq_num),
	       score_entry->score, score_entry->z, e_value, score_entry->E);
      
	/* find and print out alignment */
	parameters->full_align(q_profile, get_seq(seq_num), get_seq_length(seq_num), 1, 0, 0, 0, 0, 0, 0, &orig_score, &success);
	/* fprintf(stderr, "%d", success); */
	print_alignment(q_profile);
      }
      else if (!parameters->file_flag) break;
      if (parameters->file_flag)
	fprintf(fr, parameters->find_z_flag ? "%-25s %5d %5d  %6.2f\n" : "%-25s %5d %5d \n",
		get_id(seq_num), get_seq_length(seq_num), score_entry->score, score_entry->z);
    }
    if (parameters->file_flag) fclose(fr);
/* TEMPORARILY INACTIVATED
    show_query_hist(q_profile);
*/
    free_profile(q_profile);
  }
  printf("\n");
  return 0;
}

      
