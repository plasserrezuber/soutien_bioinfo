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

extern Parameters *parameters;
char subclone_delim, group_delim;
int n_delim, delim_count;

#define MAX_SUBCLONE_HIST 101

set_delims()
{
  subclone_delim = parameters->subclone_delim;
  group_delim = parameters->group_delim;
  n_delim = parameters->n_delim;
}

#ifdef SANGER

/* Suggestion from Richard Mott:
Another possibility would
be to modify the fasta-files so that the chemistry is stated
explicitly as the next word after the sequence name, eg:

>Xah90g12.s1t   dye_terminator
ACGGCCTT...

>Xah90g12.s1   dye_primer


One could also state if it is a forward or reverse strand  in a similar way
*/

/* hacked by rmott to cope with Sanger-style 
terminators */

/* Sanger rule for a terminator reaction is a 't' somewhere in the
extension, so this function returns true if exactly one of the entries
has a 't' */
 

set_chemistries(db)
     Database *db;
{
  int entry1;
  char *get_id();
  char c1;
  char *id1;
  
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    if (chem_parse_descrip(entry1)) continue; /* information found in descrip
					    field -- so don't parse name */
    for (id1 = get_id(entry1), delim_count = 0; 
	 *id1 && (delim_count += *id1 == subclone_delim) < n_delim; id1++);
    while (*id1 && *id1 != 't') id1++;
    if (*id1) set_chemistry(entry1, 1);
  }
}

#else
/* St Louis style 
Date: Mon, 1 Jun 1998 17:42:58 -0500
From: LaDeana Hillier <lhillier@alu.wustl.edu>
To: phg@u.washington.edu
Subject: phrap


Hi Phil -
        
        We are adding some new naming conventions and 
wondering about your willingness to incorporate this into phrap?
The ones below with asterisks are the new ones:


                  "s" for single strand read
                  "f" for forward read on double strand
                  "r" for reverse read on double strand
                  "t" for T7  (cDNAs)
                  "p" for SP6 (cDNAs)
                  "e" for T3  (cDNAs)
                  "d" for special
                  "x" forward direction SS standard dye terminators
                * "z" forward direction DS standard dye terminators
                  "y" reverse direction DS dye terminators
                * "i" forward direction SS big dye terminators
                * "b" forward direction DS big dye terminators
                * "g" reverse direction DS big dye terminators
                  "c" consensus pieces
                  "a" assembly pieces



Thanks,
LaDeana
*/
set_chemistries(db)
     Database *db;
{
  int i, entry1;
  char *get_id();
  char c1;
  char *id1;
  int descrip_chemistry[20], name_chemistry[20];
  int chem;

  for (i = 0; i < 20; i++) descrip_chemistry[i] = name_chemistry[i] = 0;

  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    if (chem_parse_descrip(entry1)) {
      descrip_chemistry[get_chemistry(entry1)] += 1;
      continue;  /*  information found in descrip
					    field -- so don't parse name */
    }
    for (id1 = get_id(entry1), delim_count = 0; 
	 *id1 && (delim_count += *id1 == subclone_delim) < n_delim; id1++);
    if (!*id1) {
      name_chemistry[0] += 1;
      continue;  /* no extension found -- leave chemistry as 0 */
    }
    c1 = *(id1 + 1);
    chem = (c1 == 'x' || c1 == 'y' || c1 == 'z') ? 1 : 
	(c1 == 'i' || c1 == 'b' || c1 == 'g' ? 2 : 
	 (c1 == 'c' || c1 == 'a' ? 3 : 0));
    set_chemistry(entry1, chem);
    name_chemistry[chem] += 1;
  }
  printf("\n\nChemistries inferred from description field:");
  printf("\n%5d  dye-primer", descrip_chemistry[0]);
  printf("\n%5d  old-dye-terminator", descrip_chemistry[1]);
  printf("\n%5d  big-dye-terminator", descrip_chemistry[2]);
  printf("\n%5d  other", descrip_chemistry[3]);
  printf("\n\nChemistries inferred from name:");
  printf("\n%5d  dye-primer", name_chemistry[0]);
  printf("\n%5d  old-dye-terminator", name_chemistry[1]);
  printf("\n%5d  big-dye-terminator", name_chemistry[2]);
  printf("\n%5d  other", name_chemistry[3]);
}

#endif

int opposite_direction(entry1, entry2)
     int entry1, entry2;
{
  return get_read_direction(entry1) != get_read_direction(entry2);
}

int different_chemistry(entry1, entry2)
     int entry1, entry2;
{
  char *id1, *id2;
  char *get_id();
  char c1, c2;
 
  return (get_chemistry(entry1) > 0) != (get_chemistry(entry2) > 0);

/*
  for (id1 = get_id(entry1); *id1 && *id1 != subclone_delim; id1++);
  if (!*id1) return 0;   
  for (id2 = get_id(entry2); *id2 && *id2 != subclone_delim; id2++);
  if(!*id2) return 0; 
  c1 = *(id1 + 1);
  c2 = *(id2 + 1);
  return ((c1 == 'x' || c1 == 'y' || c1 == 'z' || c1 == 'i' || c1 == 'b' || c1 == 'g') != 
	  (c2 == 'x' || c2 == 'y' || c2 == 'z' || c2 == 'i' || c2 == 'b' || c2 == 'g'));
*/
}

compare_templates(entry1, entry2)
     int *entry1, *entry2;
{
  Align_info *get_align_entry();
  Align_info *align_entry1, *align_entry2;
  int d;
  char *end1, *end2, *id1, *id2;

  align_entry1 = get_align_entry(*entry1);
  align_entry2 = get_align_entry(*entry2);

  end1 = align_entry1->template_end; 
  end2 = align_entry2->template_end; 

  for (id1 = align_entry1->template_start, 
       id2 = align_entry2->template_start; 
       *id1 == *id2 && id1 < end1 && id2 < end2; 
       id1++, id2++);

  if (d = (id2 >= end2) - (id1 >= end1)) return d;
  return id1 >= end1 ? 0 : *id1 - *id2;
}


set_templates(db)
     Database *db;
{
  int i, entry1;
  Align_info *get_align_entry();
  int descrip_template, name_template;
  char *id1;
  int subclone_hist[MAX_SUBCLONE_HIST];
  int *entry_ptrs;
  char *our_alloc();
  int read_count, new_flag;
  char *get_id();

  descrip_template = name_template = 0;

  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    if (template_parse_descrip(entry1)) {
      descrip_template += 1;
      continue;  /*  information found in descrip
		     field -- so don't parse name */
    }
    get_align_entry(entry1)->template_start = id1 = get_id(entry1);
    for (delim_count = 0; 
	 *id1 && (delim_count += *id1 == subclone_delim) < n_delim; id1++);
    get_align_entry(entry1)->template_end = id1;
    name_template += 1;
  }
  printf("\n\n\nTemplates inferred from description field: %5d", descrip_template );
  printf("\nTemplates inferred from name field:        %5d", name_template);

  for (i = 0; i < MAX_SUBCLONE_HIST; i++) subclone_hist[i] = 0;
  entry_ptrs = (int *)our_alloc((db->num_entries + 1) * sizeof(int));
  for (i = 0, entry1 = db->first_entry; entry1 <= db->last_entry; i++, entry1++) {
    entry_ptrs[i] = entry1;
  }
  if (i != db->num_entries) fatalError("database size");

  qsort(entry_ptrs, db->num_entries, sizeof(int), compare_templates);
  entry_ptrs[db->num_entries] = -1; /* dummy last id*/

  read_count = 1;
  for (i = 1; i <= db->num_entries; i++) {
    new_flag = 1;
    if (i < db->num_entries) {
      if (same_subclone(entry_ptrs[i], entry_ptrs[i - 1])) {
	read_count++;
	new_flag = 0;
      }
    }
      
    if (new_flag) {
      if (read_count >= MAX_SUBCLONE_HIST) {
	fprintf(stderr, "\nWARNING -- PROBABLE NOMENCLATURE VIOLATION; > %d READS PER SUBCLONE\n",
	       MAX_SUBCLONE_HIST - 1);
	printf("\nWARNING -- PROBABLE NOMENCLATURE VIOLATION; > %d READS PER SUBCLONE\n",
	       MAX_SUBCLONE_HIST - 1);
	subclone_hist[MAX_SUBCLONE_HIST - 1] += 1;
      }
      else subclone_hist[read_count] += 1;
      read_count = 1;
    }
  }
  printf("\n\nRead-template multiplicity analysis:\n # Reads      # templates");
  for (i = 0; i < MAX_SUBCLONE_HIST; i++)
    if (subclone_hist[i]) printf("\n %3d        %6d", i, subclone_hist[i]);
  our_free(entry_ptrs);
  
}

set_directions(db)
     Database *db;
{
  int i, entry1, direction;
  int descrip_direction[3], name_direction[3];
  int sense[256];
  char c1;
  char *id1;
  char *get_id();

  for (i = 0; i < 256; i++) sense[i] = 2;
  sense['f'] = sense['s'] = sense['x'] = sense['p'] = sense['z'] = sense['i'] = sense['b'] = sense['F'] = sense['H'] = 0;
  sense['r'] = sense['y'] = sense['q'] = sense['g'] = sense['R'] = 1;

  for (i = 0; i < 3; i++) descrip_direction[i] = name_direction[i] = 0;

  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    if (direction_parse_descrip(entry1)) {
      direction = get_read_direction(entry1);
      descrip_direction[direction] += 1;
      if (direction == 2) set_direction(entry1, 0);
      continue;  /*  information found in descrip
		     field -- so don't parse name */
    }
    for (id1 = get_id(entry1), delim_count = 0; 
	 *id1 && (delim_count += *id1 == subclone_delim) < n_delim; id1++);
    if (!*id1) {
      name_direction[2] += 1;
      continue;  /* no extension found -- leave direction as 0 */
    }
    c1 = *(id1 + 1);
    name_direction[sense[c1]] += 1;
    set_direction(entry1, sense[c1] == 2 ? 0 : sense[c1]);
  }
  printf("\n\nDirections inferred from description field:");
  printf("\n%5d  fwd", descrip_direction[0]);
  printf("\n%5d  rev", descrip_direction[1]);
  printf("\n%5d  unknown (set to fwd)", descrip_direction[2]);
  printf("\n\nDirections inferred from name:");
  printf("\n%5d  fwd", name_direction[0]);
  printf("\n%5d  rev", name_direction[1]);
  printf("\n%5d  unknown (set to fwd)", name_direction[2]);
}

/* N.B. FOLLOWING ASSUMES DESCRIPTION FIELD STORED IN MEMORY */
template_parse_descrip(entry1)
     int entry1;
{
  char *get_descrip();
  char *ptr;
  char *strstr();
  Align_info *get_align_entry();
  Align_info *align_entry;

  if (ptr = strstr(get_descrip(entry1), "TEMPLATE:")) {
    align_entry = get_align_entry(entry1);
    ptr += 9;
    for (; isspace(*ptr); ptr++);
    align_entry->template_start = ptr;
    for (; !isspace(*ptr); ptr++);
    align_entry->template_end = ptr;
    return 1;
  }
  return 0;
}


  
direction_parse_descrip(entry1)
     int entry1;
{
  char *get_descrip();
  char *ptr;
  char *strstr();
  char direction[30];

  if (ptr = strstr(get_descrip(entry1), "DIRECTION:")) {
    sscanf(ptr + 10, "%s", direction);
    set_direction(entry1, 
		  !strcmp(direction, "fwd") ? 0 : (!strcmp(direction, "rev") ? 1 : 2));
    return 1;
  }
  return 0;
}

/* descrip = get_descrip(entry1) */
  
Segment *repeat_parse_descrip(descrip)
     char *descrip;
{
  char *ptr;
  char *strstr();
  Segment *seg_list;
  Segment *insert_segment();
  int start, end;

  seg_list = 0;
  for (ptr = descrip; ptr = strstr(ptr, "REPEAT:"); ) {
    ptr += 7;
    sscanf(ptr, "%d %d", &start, &end);
    seg_list = insert_segment(seg_list, start, end); 
  }
  return seg_list;
}
  
chem_parse_descrip(entry1)
     int entry1;
{
  char *get_descrip();
  char *descrip, *ptr;
  char chem_type[30], dye_type[30];
  char *strstr();
  int chem, dye;
  
  descrip = get_descrip(entry1);
  chem = dye = 0;
  if (ptr = strstr(descrip, "CHEM:")) {
    sscanf(ptr + 5, "%s", chem_type);
    if (!strcmp(chem_type, "prim")) chem = 1;
    else if (!strcmp(chem_type, "term")) chem = 2;
  }
  if (ptr = strstr(descrip, "DYE:")) {
    sscanf(ptr + 4, "%s", dye_type);
    if (!strcmp(dye_type, "rhod")) dye = 1;
    else if (!strcmp(dye_type, "big")) dye = 2;
    else if (!strcmp(dye_type, "ET")) dye = 3;
    else if (!strcmp(dye_type, "d-rhod")) dye = 4;
  }
  
  if (chem == 1) {
    set_chemistry(entry1, 0);
    return 1;
  }
  else if (chem == 2) {
    set_chemistry(entry1, dye == 2 ? 2 : 1);
    return 1;
  }
  return 0;
}

int same_subclone(entry1, entry2)
     int entry1, entry2;
{
  char *id1, *id2, *end1, *end2;
  char *get_id();
  Align_info *get_align_entry();

/*  
  for (id1 = get_id(entry1), id2 = get_id(entry2), delim_count = 0; 
       *id1 == *id2 && (delim_count += *id1 == subclone_delim) < n_delim && *id1; 
       id1++, id2++);
  return *id1 == *id2;
*/

  end1 = get_align_entry(entry1)->template_end; 
  end2 = get_align_entry(entry2)->template_end; 

  for (id1 = get_align_entry(entry1)->template_start, 
       id2 = get_align_entry(entry2)->template_start; 
       *id1 == *id2 && id1 < end1; 
       id1++, id2++);
  
  return id1 >= end1 && id2 >= end2;
}

int same_group(entry1, entry2)
     int entry1, entry2;
{
  char *id1, *id2;
  char *get_id();
  
  for (id1 = get_id(entry1), id2 = get_id(entry2), delim_count = 0; 
       *id1 == *id2 && *id1 != group_delim && (delim_count += *id1 == subclone_delim) < n_delim && *id1; 
       id1++, id2++);
  
  return *id1 == group_delim && *id2 == group_delim;
}

compare_names(i,j)
     int *i, *j;
{
  char *get_id();

  return strcmp(get_id(*i), get_id(*j));
}

old_map_reads_by_name(db)
     Database *db;
{
  char *our_alloc();
  int *entry_ptrs;
  int i, j, k, i_orig, index1, index2, length, num_entries;
  char *id, *id_prev, *id0;
  char *get_id();
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();

  num_entries = db->num_entries;
  printf("\n\nSubclone/read consistency checks (* = inconsistency, L = contig link)");
  entry_ptrs = (int *)our_alloc(num_entries * sizeof(int));
  for (i = 0; i < num_entries; i++) entry_ptrs[i] = i + db->first_entry;
  qsort(entry_ptrs, num_entries, sizeof(int), compare_templates);
/*
  id_prev = get_id(entry_ptrs[0]);
*/
  i_orig = 0;

  for (i = 1; i < num_entries; i++) {
    id = get_id(entry_ptrs[i]);
/*
    for (j = 0; id[j] && id[j] != subclone_delim && id[j] == id_prev[j]; j++);
    if (id[j] && id[j] != subclone_delim) i_orig = i;
*/
    if (!same_subclone(entry_ptrs[i], entry_ptrs[i_orig])) i_orig = i;
    else for (k = i - 1; k >= i_orig; k--) {
      id0 = get_id(entry_ptrs[k]);
      if (1 /* id[j + 1] != id0[j + 1] */) {
	align_entry = get_align_entry(entry_ptrs[i]);
	align_entry2 = get_align_entry(entry_ptrs[k]);
	index1 = align_entry->contig->index;
	index2 = align_entry2->contig->index;
	printf("\n%-15s  %-15s ", id, id0);
	printf(" %5d  %5d    %5d  %5d",
	       index1, align_entry->reverse,
	       index2, align_entry2->reverse);
	
	if (!index1 && !index2) continue;
	if (index1 != index2) printf(" L");
	else if (id[j + 1] == 's' && id0[j + 1] == 'r'
#ifdef SANGER
               || id[j + 1] == 'q' && id0[j + 1] == 'p'
#endif /*SANGER*/
		 || id[j + 1] == 'r' && id0[j + 2] == 's') { /* 2d condition
								allows for .ls
								reads */
	  if (align_entry->reverse == align_entry2->reverse)
	    printf(" *");
	  else {
	    length = align_entry->reverse ?
	      align_entry->end - align_entry2->start :
		align_entry2->end - align_entry->start;
	    printf(" %c %6d", length < 0 ? '*' : ' ', length);
	  }
	}
	else {
	  if (align_entry->reverse != align_entry2->reverse)
	    printf(" *");
	  length = align_entry->reverse ?
	      align_entry2->end - align_entry->end :
		align_entry->start - align_entry2->start;
	    printf("          %6d", length);
	}
      }
    }
    id_prev = id;  
  }
  our_free(entry_ptrs);
}

typedef struct read_pair {
  int contig1, contig2, read1, read2, displace;
  int start, end; /* only for opp_sense reads, within a contig */
  char direction, complement, conflict, opp_sense;
} Read_pair;

compare_read_pairs(pair1, pair2)
     Read_pair *pair1, *pair2;
{
  int d;
  if (d = pair1->contig1 - pair2->contig1)
    return d;
  if (d = (pair1->contig1 == pair1->contig2) - (pair2->contig1 == pair2->contig2))
    return d;
  if (d = pair1->opp_sense - pair2->opp_sense)
    return d;
  if (d = pair1->direction - pair2->direction)
    return d;
  if (d = pair1->contig2 - pair2->contig2)
    return d;
  if (d = pair1->complement - pair2->complement)
    return d;
  if (d = pair1->conflict - pair2->conflict)
    return d;
  return pair1->start - pair2->start;
}

map_reads_by_name(db)
     Database *db;
{
  char *our_alloc();
  int *entry_ptrs;
  int i, j, k, i_orig, contig1, contig2, length, num_entries;
  char *id, *id_prev, *id0;
  char *get_id();
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();
  int num_read_pairs;
  int prev_contig1, prev_contig2, prev_direction, prev_complement, prev_opp_sense;
  int direction, complement, opp_sense;
  Read_pair *read_pairs, *read_pair;
  int sense[256];
  int start, start2, end, end2, length1, length2, contig_tail, contig_tail2, contig_lead, contig_lead2;
  int *size_hist;
  int hist_granularity, n_hist_entries, hist_entry;
  int link_flag;
  int *depth_shifts;
  int c_length, cum;

  notify("Mapping reads by name ...");
  /*
     for (i = 0; i < 256; i++) sense[i] = 0;
     sense['f'] = sense['s'] = sense['x'] = sense['p'] = sense['z'] = sense['i'] = sense['b'] = sense['F'] = sense['H'] = 1;
     sense['r'] = sense['y'] = sense['q'] = sense['g'] = sense['R'] = 2;
     */
  hist_granularity = 200;
  n_hist_entries = 1 + parameters->max_subclone_size / hist_granularity;
  size_hist = (int *)our_alloc(n_hist_entries  * sizeof(int));
  for (i = 0; i < n_hist_entries; i++) size_hist[i] = 0;
  num_entries = db->num_entries;
  printf("\n\nSubclone/read contig links and consistency checks (* = inconsistency; Contig 0 = singlets)");
  printf("\nMax subclone size: %d", parameters->max_subclone_size);
  entry_ptrs = (int *)our_alloc(num_entries * sizeof(int));
  for (i = 0; i < num_entries; i++) entry_ptrs[i] = i + db->first_entry;
  qsort(entry_ptrs, num_entries, sizeof(int), compare_templates);
  /*
     id_prev = get_id(entry_ptrs[0]);
     */
  i_orig = 0;
  num_read_pairs = 0;
  for (i = 1; i < num_entries; i++) {
    /*
       id = get_id(entry_ptrs[i]);
       for (j = delim_count = 0; id[j] && (delim_count += id[j] == subclone_delim) < n_delim && id[j] == id_prev[j]; j++);
       if (id[j] && delim_count < n_delim) i_orig = i;
       */
    if (!same_subclone(entry_ptrs[i], entry_ptrs[i_orig])) i_orig = i;
    else for (k = i - 1; k >= i_orig; k--) {
      /*
	 id0 = get_id(entry_ptrs[k]);
	 */
      align_entry = get_align_entry(entry_ptrs[i]);
      align_entry2 = get_align_entry(entry_ptrs[k]);
      if (align_entry->contig->index != align_entry2->contig->index)
	num_read_pairs += 2;
      else num_read_pairs += 1;
    }
    /*
       id_prev = id;  
       */
  }
  
  read_pairs = (Read_pair *)our_alloc(num_read_pairs * sizeof(Read_pair));
  /*
     id_prev = get_id(entry_ptrs[0]);
     */
  i_orig = 0;
  read_pair = read_pairs;
  for (i = 1; i < num_entries; i++) {
    /*
       id = get_id(entry_ptrs[i]);
       for (j = delim_count = 0; id[j] && (delim_count += id[j] == subclone_delim) < n_delim && id[j] == id_prev[j]; j++);
       if (id[j] && delim_count < n_delim) i_orig = i;
       */
    if (!same_subclone(entry_ptrs[i], entry_ptrs[i_orig])) i_orig = i;
    else for (k = i - 1; k >= i_orig; k--) {
      /*
	 id0 = get_id(entry_ptrs[k]);
	 */
      align_entry = get_align_entry(entry_ptrs[i]);
      align_entry2 = get_align_entry(entry_ptrs[k]);

      length1 = get_seq_length(align_entry->seq_entry);
      length2 = get_seq_length(align_entry2->seq_entry);
      start = align_entry->start + align_entry->m_start - 2;
      start2 =  align_entry2->start + align_entry2->m_start - 2;
      end = align_entry->end - length1 + align_entry->m_end - 1;
      end2 = align_entry2->end - length2 + align_entry2->m_end - 1;

      contig_tail = align_entry->reverse ? end : align_entry->contig->length - start;
      contig_tail2 = align_entry2->reverse ? end2 : align_entry2->contig->length - start2;
      contig_lead = align_entry->contig->length - contig_tail;
      contig_lead2 = align_entry2->contig->length - contig_tail2;
      /*      
	 if (!strcmp(get_id(entry_ptrs[i]), "XWW600030266.R1") 
	 || !strcmp(get_id(entry_ptrs[i]), "XWW600030266.F1")
	 || !strcmp(get_id(entry_ptrs[k]), "XWW600030266.R1") 
	 || !strcmp(get_id(entry_ptrs[k]), "XWW600030266.F1")) {
	 fprintf(stderr, "\n\n%s %d %d  %d %d %d %d ", 
	 get_id(entry_ptrs[i]), contig_tail, contig_lead, start, end, 
	 align_entry->contig->length, align_entry->reverse);
	 fprintf(stderr, "\n\n%s %d %d  %d %d %d %d ", 
	 get_id(entry_ptrs[k]), contig_tail2, contig_lead2, start2, end2, 
	 align_entry2->contig->length, align_entry2->reverse);
	 }
	 */

      read_pair->read1 = entry_ptrs[i];
      read_pair->read2 = entry_ptrs[k];
      read_pair->contig1 = contig1 = align_entry->contig->index;
      read_pair->contig2 = contig2 = align_entry2->contig->index;
      read_pair->direction = align_entry->reverse;
      read_pair->conflict = 0;
      read_pair->opp_sense = opp_sense = opposite_direction(entry_ptrs[i], entry_ptrs[k]);
      /*
	 (sense[id[j + 1]] && sense[id0[j + 1]] && 
	 sense[id[j + 1]] != sense[id0[j + 1]]);
	 */
      /*
	 (id[j + 1] == 's' && id0[j + 1] == 'r'
	 || id[j + 1] == 'q' && id0[j + 1] == 'p'
	 || id[j + 1] == 'y' && id0[j + 1] == 'x'
	 || id[j + 1] == 'y' && id0[j + 1] == 's'
	 || id[j + 1] == 'x' && id0[j + 1] == 'r'
	 || id[j + 1] == 'r' && id0[j + 2] == 's'));   allows for .ls reads 
	 */

      read_pair->complement = complement =
	opp_sense != (align_entry->reverse != align_entry2->reverse);

      if (contig1 == contig2) {
	if (align_entry->reverse) {
	  read_pair->start = start2;
	  read_pair->end = end;
	}
	else {
	  read_pair->start = start;
	  read_pair->end = end2;
	}
	if (!contig1) {
	  read_pair->displace = 0;
	  read_pair++;
	  continue;
	}
	read_pair->conflict = complement;

	if (opp_sense) {
	  read_pair->displace = read_pair->end - read_pair->start;

	  if (read_pair->displace < 0 || read_pair->displace > parameters->max_subclone_size || read_pair->conflict) {
	    read_pair->conflict = 1;
	    size_hist[n_hist_entries - 1] += 1;
	  }
	  else {
	    hist_entry = read_pair->displace / hist_granularity;
	    size_hist[hist_entry < n_hist_entries ? hist_entry : n_hist_entries - 1] += 1;
	  }
	}
	else {
	  read_pair->displace = align_entry->reverse ?
	    end2 - end : start - start2;
	}
	read_pair++;
      }
      else {
	read_pair->start = align_entry->reverse ? end : start;
	read_pair->end = align_entry2->reverse ? end2 : start2;
	
	read_pair->displace = opp_sense ? contig_tail + contig_tail2 : 
	(contig_tail + contig_lead2 < contig_tail2 + contig_lead ?
	 contig_tail + contig_lead2 : contig_tail2 + contig_lead);

	read_pair->conflict = read_pair->displace > parameters->max_subclone_size;
	read_pair++;
	read_pair->read1 = entry_ptrs[k];
	read_pair->read2 = entry_ptrs[i];
	read_pair->contig1 = contig2;
	read_pair->contig2 = contig1;
	read_pair->direction = align_entry2->reverse;
	read_pair->complement = complement;
	read_pair->displace = (read_pair - 1)->displace;
	read_pair->conflict = (read_pair - 1)->conflict;
	read_pair->start = (read_pair - 1)->end;
	read_pair->end = (read_pair - 1)->start;
	read_pair->opp_sense = opp_sense;
	read_pair++;
      }
    }
    id_prev = id;  
  }
  
  if (read_pairs + num_read_pairs != read_pair) 
    notify("\nWARNING -- READ PAIR DISCREPANCY\n");
  qsort(read_pairs, num_read_pairs, sizeof(Read_pair), compare_read_pairs);

  c_length = 0;
  depth_shifts = 0;
  for (i = 0; i < num_read_pairs; i++) {
    read_pair = read_pairs + i;
    contig1 = read_pair->contig1;
    direction = read_pair->direction;
    contig2 = read_pair->contig2;
    complement = read_pair->complement;
    opp_sense = read_pair->opp_sense;
    if (i == 0 || contig1 != prev_contig1 
	|| (contig1 != contig2 && direction != prev_direction) 
	|| contig2 != prev_contig2 || opp_sense != prev_opp_sense 
        || (contig1 != contig2 && complement != prev_complement)) {
      if (c_length) {
	print_covered_regions(depth_shifts, c_length, 1);
	c_length = 0;
	our_free(depth_shifts);
      }
      
      if (contig1) {
	link_flag = 0;
	if (contig1 == contig2) {
    
	  printf("\n\nINTERNAL Contig %d  %s  ", contig1,
		 opp_sense ? "opp sense" : "same sense");
	  if (opp_sense) {
	    c_length = get_align_entry(read_pair->read1)->contig->length;
	    depth_shifts = (int *)our_alloc((c_length + 1) * sizeof(int));
	    for (j = 0; j <= c_length; j++) depth_shifts[j] = 0;
	  }
	}
	else
	  printf("\n\nContig %d %s  %s    %s Contig %d", contig1,
		 opp_sense ? "opp sense" : "same sense",
		 direction ? "LEFT LINK:" : "RIGHT LINK:", 
		 complement ? "complement" : "", contig2);
      }
      prev_contig1 = contig1;
      prev_direction = direction;
      prev_contig2 = contig2;
      prev_complement = complement;
      prev_opp_sense = opp_sense;
    }
    if (contig1) {
      if (c_length && !read_pair->conflict) {
	depth_shifts[read_pair->start < 0 ? 0 :
		     (read_pair->start >= c_length ? c_length - 1 : read_pair->start)] += 1;
	depth_shifts[read_pair->end < 0 ? 0 :
		     (read_pair->end >= c_length ? c_length - 1 : read_pair->end)] -= 1;
      }
      if (opp_sense && contig1 != contig2 && !read_pair->conflict) {
	if (!link_flag) {
	  alloc_link(contig1, contig2, direction, complement);
	  link_flag = 1;
	}
	increment_n_pairs();
      }
      align_entry = get_align_entry(read_pair->read1);
      align_entry2 = get_align_entry(read_pair->read2);
      printf("\n %c  %c %-15s  %c %-15s", read_pair->conflict ? '*' : ' ',
	     align_entry->reverse ? 'C' : ' ',
	     get_id(read_pair->read1), 
	     align_entry2->reverse ? 'C' : ' ',
	     get_id(read_pair->read2));
      printf("  %6d  %6d  %6d", read_pair->displace, read_pair->start, read_pair->end);
    }
  }
  if (c_length) {
    print_covered_regions(depth_shifts, c_length, 1);
    c_length = 0;
    our_free(depth_shifts);
  }
  
  printf("\n\nSize histogram for consistent forward-reverse pairs (*** = inconsistent pairs)");
  for (i = 0; i < n_hist_entries; i++) {
    if (size_hist[i] || hist_granularity * i >= parameters->max_subclone_size) {
      if (hist_granularity * i >= parameters->max_subclone_size)
	printf("\n  ***  %4d", size_hist[i]);
      else printf("\n%5d  %4d", hist_granularity * i, size_hist[i]);
    }
  }

  print_contig_chains();
  print_links();

  our_free(read_pairs);
  our_free(entry_ptrs);
  our_free(size_hist);
  notify(" Done\n");
}

/* assumptions: depths has transition marks (converted to coverage by function);
   has length equal to length + 1 */
print_covered_regions(depths, length, depth_cutoff)
     int *depths;
     int length, depth_cutoff;
{
  int j, cum;

  for (j = cum = 0; j <= length; j++) cum = depths[j] += cum;
  printf("\nCovered regions: ");
  depths[0] = depths[length] = 0;
  for (j = 1; j <= length; j++) {
    if (depths[j] >= depth_cutoff && depths[j - 1] < depth_cutoff) printf("%d..", j);
    if (depths[j] < depth_cutoff && depths[j - 1] >= depth_cutoff) printf("%d ", j - 1);
  }
}


compare_ids(id1, id2)
     char **id1, **id2;
{
  return strcmp(*id1, *id2);
}
 
name_check(db)
     Database *db;
{
  char *our_alloc();
  char **name_ptrs;
  int entry1, i, j, c1, c2, read_count;
  char *get_id();
  int suffix[256];
  int subclone_hist[MAX_SUBCLONE_HIST];
  int new_flag;

  printf("\n\nRead name analysis:");
  for (i = 0; i < MAX_SUBCLONE_HIST; i++) subclone_hist[i] = 0;
  for (i = 0; i < 256; i++) suffix[i] = 0;

  set_delims();

  name_ptrs = (char **)our_alloc((db->num_entries + 1) * sizeof(char *));
  for (i = 0, entry1 = db->first_entry; entry1 <= db->last_entry; i++, entry1++) {
    name_ptrs[i] = get_id(entry1);
  }
  if (i != db->num_entries) fatalError("database size");

  qsort(name_ptrs, db->num_entries, sizeof(char *), compare_ids);
  name_ptrs[db->num_entries] = 0; /* dummy last id of 0 */

  read_count = 1;
  for (i = 1; i <= db->num_entries; i++) {
    new_flag = 1;
    delim_count = j = 0;
    c1 = name_ptrs[i - 1][j];
    if (i < db->num_entries) {
      for ( ; ; j++) {
	c1 = name_ptrs[i - 1][j];
	c2 = name_ptrs[i][j];
	delim_count += c1 == subclone_delim;
	if (c1 != c2 || delim_count >= n_delim || !c1) break;
      }
      if (c1 == c2) {
	read_count++;
	new_flag = 0;
      }
    }
      
    if (new_flag) {
      if (read_count >= MAX_SUBCLONE_HIST) {
	subclone_hist[MAX_SUBCLONE_HIST - 1] += 1;
      }
      else subclone_hist[read_count] += 1;
      read_count = 1;
    }
    if (delim_count < n_delim && c1) /* problem if i == db->num_entries */
      for (j++; 
	   (c1 = name_ptrs[i - 1][j]) 
	   && (delim_count += c1 == subclone_delim) < n_delim; 
	   j++);
    if (c1) suffix[name_ptrs[i - 1][j + 1]] += 1;
    else {
      suffix[0] += 1;
    }
  }

  printf("\n # Reads      # templates");
  for (i = 0; i < MAX_SUBCLONE_HIST; i++)
    if (subclone_hist[i]) printf("\n %3d        %6d", i, subclone_hist[i]);
  printf("\n\n Suffix counts:");
  if (suffix[0]) printf("\n(no suffix) %3d", suffix[0]);
  for (i = 1; i < 256; i++)
    if (suffix[i]) printf("\n  %c         %5d", (char)i, suffix[i]);

  our_free(name_ptrs);
}

typedef struct link {
  int contig1, contig2, n_pairs;
  char direction, complement, used_in_chain;
  struct link *next;
} Link;

Link *head_link;
int n_links;

alloc_link(contig1, contig2, direction, complement)
     int contig1, contig2;
     int direction, complement;
{
  Link *link;
  char *our_alloc();

  link = (Link *)our_alloc(sizeof(Link));
  link->contig1 = contig1;
  link->contig2 = contig2;
  link->direction = direction;
  link->complement = complement;
  link->n_pairs = 0;
  link->used_in_chain = 0;
  link->next = head_link;
  head_link = link;
  n_links++;
}

increment_n_pairs()
{
  head_link->n_pairs += 1;
}

compare_link_ptrs(link1, link2)
     Link **link1, **link2;
{
  int d;

  if (d = (*link2)->n_pairs - (*link1)->n_pairs)
    return d;
  return (*link2)->contig1 + (*link2)->contig2 - ((*link1)->contig1 + (*link1)->contig2);
/* sort in decreasing order first by number of links, then by sum of contig indices
   (which is surrogate for combined sizes) */
}

print_links()
{
  Link *link;

  printf("\n\n Consistent opp sense links (* = not used in chain, ** = multiple non-zero):");
  for (link = head_link; link; link = link->next) {
    printf("\n%sContig %d  %s  %s Contig %d  %d pairs",
	   link->used_in_chain ? " " : 
	   (link->n_pairs > 1 && link->contig2 ? "**" : "*"),
	   link->contig1, link->direction ? "LEFT " : "RIGHT", 
	   link->complement ? "complement" : "", 
	   link->contig2, link->n_pairs);
  }
}

/* greedy algorithm for finding chains of contigs */
/* NEED TO LOOK SYSTEMATICALLY FOR CIRCULAR CHAINS, AND BREAK THEM */
print_contig_chains()
{
  Link *link;
  char *our_alloc();
  int i, j, max_contig, direction, complement, pass;
  Link **link_ptrs;
  int *contig_links[2];
  char *contig_marks;

  if (!n_links) return;
  link_ptrs = (Link **)our_alloc(n_links * sizeof(Link *));
  max_contig = head_link->contig1;

  for (i = 0; i < 2; i++) {
    contig_links[i] = (int *)our_alloc((max_contig + 1) * sizeof(int));
    for (j = 0; j <= max_contig; j++)
      contig_links[i][j] = 0;
  }

  contig_marks = (char *)our_alloc((max_contig + 1) * sizeof(char));
  for (j = 0; j <= max_contig; j++) contig_marks[j] = 0;
  for (link = head_link, i = 0; link; link = link->next, i++) 
    link_ptrs[i] = link;

  qsort(link_ptrs, n_links, sizeof(Link *), compare_link_ptrs);

  for (i = 0; i < n_links; i++) {
    link = link_ptrs[i];
    if (!link->contig1 || !link->contig2) continue;
    direction = link->direction == link->complement;
    if (contig_links[link->direction][link->contig1] || contig_links[direction][link->contig2]) {
      if (contig_links[link->direction][link->contig1] == 
	   (link->complement ? -(link->contig2) : link->contig2)
	  && contig_links[direction][link->contig2] == (link->complement ? -(link->contig1) : link->contig1))
	  link->used_in_chain = 1;
      continue;
    }
    link->used_in_chain = 1;

/*
    printf("\n%d %d %d", link->n_pairs, link->contig1, link->contig2);
*/
    contig_links[link->direction][link->contig1] = 
      link->complement ? -(link->contig2) : link->contig2;
    contig_links[direction][link->contig2] = 
      link->complement ? -(link->contig1) : link->contig1;
  }
  printf("\n\nContig chains:");
  for (pass = 0; pass < 2; pass++) {
    printf("\n Pass: %d", pass + 1);
    for (i = max_contig; i > 0; i--) {
      if (contig_marks[i] || 
	  !contig_links[1][i] && !contig_links[0][i]) 
	continue; /* already printed out, or unlinked */
      if (contig_links[1][i] && !pass) continue;  
        /* only process if nothing to the left,
	   or circular chain on 2d pass*/
      printf("\n  %d  ", i);
      contig_marks[i] = 1;
      direction = complement = 0;
      
      for (j = contig_links[0][i]; j; j = contig_links[direction][j]) {
	if (j < 0) {
	  complement = !complement;
	  direction = !direction;
	  j = -j;
	}
	printf("%s%d  ", complement ? "C " : "", j);
	if (contig_marks[j]) {
	  printf("**REPEATED CONTIG** ");
	  break;
	}
	else contig_marks[j] = 1;
      }
    }
  }
  our_free(contig_links[0]);
  our_free(contig_links[1]);
  our_free(contig_marks);
  our_free(link_ptrs);
}

set_repeat_tags(db)
     Database *db;
{
  int entry1, n;
  char *get_id();
  char c1;
  char *id1;
  Align_info *get_align_entry();
  Align_info *align_entry1;
  Segment *repeat_parse_descrip();
  Segment *seg_list;
  char *get_descrip();

  n = 0;
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    seg_list = repeat_parse_descrip(get_descrip(entry1));
    align_entry1 = get_align_entry(entry1);
    for (; seg_list; seg_list = seg_list->next) {
      n++;
      append_tag(align_entry1, "repeat", seg_list->start, seg_list->end);
      /*      our_free(tag); NEED WAY TO DO THIS */
    }
  }
  
  printf("\n\n%d repeat tags created\n", n);
  if (!n) {
    parameters->repeat_screen &= 2;
    printf("query repeat screening turned off\n");
  }
  printf("\n");
}

