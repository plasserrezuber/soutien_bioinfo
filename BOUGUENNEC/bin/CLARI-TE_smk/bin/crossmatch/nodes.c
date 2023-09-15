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

#define BIGNEG -99999
#define EDIT_WEIGHT 500 /* score assigned to marked high-quality base */
 
extern FILE *fp_log; /* log file */
extern Parameters *parameters; /* exported */
extern double *scaled_err_probs;
int *contig_graph_weights;
 
set_contig_graph_weights()
{
  int i;
  char *our_alloc(); 

  contig_graph_weights = 1 + (int *)our_alloc(257 * sizeof(int));
  for (i = 0; i < 256; i++)
    contig_graph_weights[i] = parameters->contig_graph_weights ? (int)(10000. * scaled_err_probs[i]) : i;
  contig_graph_weights[99] = parameters->contig_graph_weights ? EDIT_WEIGHT + contig_graph_weights[99] : EDIT_WEIGHT;
  contig_graph_weights[-1] = parameters->contig_graph_weights ? contig_graph_weights[0] : -10;
/* value used for N */
}

/* function to search the directed graph corresponding to all matching points in alignments */

copy_node(node1, node2)
     Node *node1, *node2;
{
  node1->score = node2->score;
  node1->best = node2->best;
  node1->link = node2->link;
  node1->entry1 = node2->entry1;
  node1->site = node2->site;
}

static int n_nodes, num_visits;
static Node *node_array;
static int *val, *stack;
static int id, p, max_p, n_components;
static Align_info *align_array;

compare_nodes(node1, node2)
     int *node1, *node2;
{
  int d;

  if (d = align_array[node_array[*node2].entry1].start - align_array[node_array[*node1].entry1].start)
    return d; /* sort so that deepest nodes are at end of list -- to keep recursion depth 
		 for tarjan calls to a minimum */
  else if (d = node_array[*node1].entry1 - node_array[*node2].entry1) return d;
  else if (d = node_array[*node1].site - node_array[*node2].site) return d;
  else return node_array[*node1].link - node_array[*node2].link; 
  /* latter cond'ns to ensure that nodes at beginning of sequence are in
   proper locations */
}
 
node_master(contig, db)
     Contig *contig;
     Database *db;
{
  int i, j, k, i_seq;
  int entry1, i_pass, position, best_length;
  int site, max_score, best_node;
  int l_start, r_start, pos_j;
  Align_info *align_entry;
  Align_info *get_align_entry();
  char *our_alloc();
  int *positions;
  char *get_id(), *get_orig_qual(), *get_adj_qual();
  unsigned char *get_seq(), *get_comp_seq();
  unsigned char *seq;
  char *a_qual, *o_qual;
  int length1, reverse, length, q;
  int last_printed;

  align_array = db->align_array;
  num_visits = 0;
  count_nodes(contig);

  make_nodes(contig);

/* Now check for non-degenerate cycles! */
  n_components = 0;
  fprintf(fp_log,"\n\nContig %d:  %d nodes", contig->index, n_nodes);

  val = (int *)our_alloc(n_nodes * sizeof(int));
  stack = (int *)our_alloc(n_nodes * sizeof(int));
  for (i = 0; i < n_nodes; i++) val[i] = 0;
  id = p = max_p = 0;

  for (i = 0; i < n_nodes; i++) {
    if (!val[i]) nr_tarjan(i); 
  }
  our_free(val);
  our_free(stack);

  fprintf(fp_log,"\n\n %d str. conn. components", n_components);

  /* code to find reads with "disjoined" (i.e. non-overlapping) links
     for (i = 0; i < n_nodes; i++) {
     entry1 = node_array[i].entry1;
     min_break = i;
     max_break = n_nodes;
     if (i == 0 || entry1 != node_array[i - 1].entry1) {
     
     fprintf(fp_log,"\n%-15s  %d", align_array[entry1].db_entry->id, 
     align_array[entry1].start + node_array[i].score + node_array[i].site);
     
     for (j = i; j < n_nodes && node_array[j].entry1 == entry1; j++) {
     entry2 = node_array[j].entry2;
     for (k = i; k < n_nodes && node_array[k].entry1 == entry1; k++) 
     if (node_array[k].entry2 == entry2 && k != j) break;
     if (k < max_break && j < max_break) 
     max_break = k < j ? j : k;
     if (k > min_break && j > min_break) 
     min_break = k < j ? k : j;
     }
     if (min_break > max_break) { 
     fprintf(fp_log,"\n%-15s %d %-15s %d %-15s   %d", 
     align_array[entry1].db_entry->id, node_array[min_break].site, 
     align_array[node_array[min_break].entry2].db_entry->id,
     node_array[max_break].site,
     align_array[node_array[max_break].entry2].db_entry->id, j - i);
     }
     }
     }
     */
  
  max_score = 0;
  best_node = 0;
  for (i = 0; i < n_nodes; i++) {
/*
    if (!strcmp(get_id(node_array[i].entry1), "iv80d05.s1"))
      fprintf(stderr,"iv80d05.s1 %d  %d  ", node_array[i].score, node_array[i].site);
*/
    if (node_array[i].link == BIGNEG && node_array[i].site == 0) {
      if (node_array[i].score > max_score) {
	max_score = node_array[i].score;
	best_node = i;
      }
/*
      if (!strcmp(get_id(node_array[i].entry1), "iv84c12.x1"))
	fprintf(stderr,"\niv84c12.x1 %d", node_array[i].score);
      if (!strcmp(get_id(node_array[i].entry1), "iv80d05.s1"))
	fprintf(stderr,"\niv80d05.s1 %d", node_array[i].score);
*/
    }
  }

  /* extend to left as much as possible -- this is NOT OPTIMAL */
  best_length = 0;
  for (i = best_node; i < n_nodes; i++) {
    if (node_array[i].link == BIGNEG && node_array[i].site == 0 &&
	node_array[i].score == max_score) {
      entry1 = node_array[i].entry1;
      a_qual = get_adj_qual(entry1);
      length = get_seq_length(entry1);
      align_entry = get_align_entry(entry1);
      for (j = 0; j < length; j++) {
	if (q = a_qual[align_entry->reverse ? length - 1 - j : j]) break;
      }
      if (j > best_length) {
	best_length = j;
	best_node = i;
      }
/*
      for (j = i + 1; j < n_nodes && node_array[j].score == max_score
	   && node_array[j].entry1 == entry1; j++)
	if (node_array[j].site > best_length) {
	  best_length = node_array[j].site;
	  best_node = i;
	}
*/
    }
  }

  i_seq = 0;
  entry1 = BIGNEG;
  fprintf(fp_log,"\n\nPath");

  positions = (int *)our_alloc(n_nodes * sizeof(int));
  for (i = 0; i < n_nodes; i++) positions[i] = BIGNEG;

  for (i = best_node; i != BIGNEG; i = node_array[i].best) {
/* should be unnec -- with tarjan
    if (node_array[i].used) {
      fprintf(fp_log,"\nERROR: INFINITE LOOP!");
      fprintf(fp_log,"\n%5d  %-13s %5d  %d  %d", 
	     i_seq, align_array[node_array[i].entry1].db_entry->id, node_array[i].site,
	     i, node_array[i].score);
      break;
    }
    node_array[i].used = 1;
*/
    if (entry1 != node_array[i].entry1) {
      entry1 = node_array[i].entry1; /* equiv. link - so no sequence to
					be appended */
      reverse = get_align_entry(entry1)->reverse;
      k = node_array[i].site;
      append_base_segment(contig, entry1, k + 1, i_seq + 1);
      positions[i] = i_seq;
      if (node_array[i].link > BIGNEG && positions[node_array[i].link] == BIGNEG)
	  positions[node_array[i].link] = positions[i];

      fprintf(fp_log,"\n%5d  %c %-13s %5d ", 
	      i_seq + 1, reverse ? 'C' : ' ', 
	      get_id(entry1), 
	      reverse ? get_seq_length(entry1) - node_array[i].site : node_array[i].site + 1);
/*	     node_array[i].swat_score); */
      
    }
    else {
/*
      fprintf(fp_log,"\n*%5d  %-13s %5d  %d  %d  %d", i_seq, 
      align_array[entry1].db_entry->id, k,
	     node_array[i].site, i, node_array[i].score);
*/
      while (k < node_array[i].site) {
	i_seq++;
	k++;
      }
      contig->base_segment->read_end = k;
      contig->base_segment->contig_end = i_seq;

      for (j = i; j >= 0 && entry1 == node_array[j].entry1 
	   && positions[j] == BIGNEG; j--) {
	positions[j] = i_seq + node_array[j].site - k;
	if (node_array[j].link > BIGNEG) {
	  if (positions[node_array[j].link] == BIGNEG)
	    positions[node_array[j].link] = positions[j];
	  else break;
	}
      }
    }
  }

  fprintf(fp_log,"\nContig length: old %d, new %d", contig->length, i_seq);
  contig->length = i_seq;
  contig->seq = (unsigned char *)our_alloc((i_seq + 1) * sizeof(unsigned char));
  contig->orig_qual = (/* signed */ char *)our_alloc((i_seq + 1) * sizeof(char));
  contig->adj_qual = (/* signed */ char *)our_alloc((i_seq + 1) * sizeof(char));

  i_seq = 0;
  entry1 = BIGNEG;

  for (i = best_node; i != BIGNEG; i = node_array[i].best) {
    if (entry1 != node_array[i].entry1) {
      entry1 = node_array[i].entry1; /* equiv. link - so no sequence to
					be appended */
      reverse = get_align_entry(entry1)->reverse;
      seq = reverse ? get_comp_seq(entry1) : get_seq(entry1);
      a_qual = get_adj_qual(entry1);
      o_qual = get_orig_qual(entry1);
      k = node_array[i].site;
      length1 = get_seq_length(entry1);
    }
    else {
      while (k < node_array[i].site) {
	contig->adj_qual[i_seq] = a_qual[reverse ? length1 - 1 - k : k];
	contig->orig_qual[i_seq] = o_qual[reverse ? length1 - 1 - k : k];
	contig->seq[i_seq++] = seq[k++];
      }
    }
  }
  contig->seq[i_seq] = 0;

  contig->score = max_score;
/* fill in other positions (should ideally fill in nearest entries only!) */
/* make all equiv node go to same position? */
  last_printed = -1;
  for (i_pass = 0; i_pass < 3; i_pass++) {
    entry1 = BIGNEG;
    for (i = 0; i < n_nodes; i++) {
      if (positions[i] != BIGNEG) {
	site = node_array[i].site;
	position = positions[i];
	if (entry1 != node_array[i].entry1) { /* first positioned site for this
						 read */
	  entry1 = node_array[i].entry1;
	  for (j = i - 1; j >= 0 && node_array[j].entry1 == entry1; j--) {
	    positions[j] = position + node_array[j].site - site;
	    if (node_array[j].link > BIGNEG && positions[node_array[j].link] == BIGNEG)
	      positions[node_array[j].link] = positions[j];
	  }
	}
      }
      else if (entry1 == node_array[i].entry1) {
	positions[i] = position + node_array[i].site - site;
	if (node_array[i].link > BIGNEG && positions[node_array[i].link] == BIGNEG)
	  positions[node_array[i].link] = positions[i];
      }
      else if (i_pass == 2) {
	if (node_array[i].entry1 != last_printed) {
	  fprintf(stderr, "\nUNPOSITIONED READ: %s", 
		  get_id(node_array[i].entry1));
	  fprintf(stdout, "\nUNPOSITIONED READ: %s", 
		  get_id(node_array[i].entry1));
	  fprintf(fp_log, "\nUNPOSITIONED READ: %s", 
		  get_id(node_array[i].entry1));
	  last_printed = node_array[i].entry1;
	}
      }
    }
  }
  
  for (i = 0; i < n_nodes; i++) {
    if (node_array[i].link == BIGNEG && node_array[i].site == 0) {
      entry1 = node_array[i].entry1;
      if (positions[i] != BIGNEG) {
	align_entry = align_array + entry1;
/*	align_entry->adjust = positions[i] - align_entry->start; */
	l_start = r_start = positions[i];
	for (j = i + 1; j < n_nodes && node_array[j].entry1 == entry1; j++) {
	  pos_j = positions[j] - node_array[j].site;
	  if (pos_j < l_start) l_start = pos_j;
	  else if (pos_j > r_start) r_start = pos_j;
	}
	align_entry->start = l_start;
	align_entry->end = r_start + get_seq_length(entry1) - 1;
	if (r_start > l_start + 20)
	  fprintf(fp_log,"\nBadly positioned read start: %s, %d - %d",
		  get_id(entry1), l_start, r_start);
      }
      else fprintf(fp_log,"\nUnpositioned read: %s", get_id(entry1));
    }
  }
  our_free(positions);  
  our_free(node_array);  
}

append_base_segment(contig, entry, read_start, contig_start)
     Contig *contig;
     int entry, read_start, contig_start;
{
  char *our_alloc();
  Base_segment *base_segment;

  base_segment = (Base_segment *)our_alloc(sizeof(Base_segment));
  base_segment->next = contig->base_segment;
  base_segment->entry = entry;
  base_segment->read_start = read_start;
  base_segment->contig_start = contig_start;
  contig->base_segment = base_segment;
}


count_nodes(contig)
     Contig *contig;
{
  Align_info *align_entry;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();

  n_nodes = 2 * contig->num_entries; /* + 4 * contig->num_matches;    */
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) { 
    for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
      if (!is_used(pair)) continue;
/*
      if (!is_best(pair)) continue;
*/
      if (pair->entry1 >= pair->entry2) continue; 
      n_nodes += count_pair_nodes(pair);
    }
  }
}

int count_pair_nodes(pair)
  Aligned_pair *pair;
{
  unsigned char *diff;
  int first_disp, last_disp;
  int nn;
  unsigned char *seq1, *seq2;
  unsigned char *get_seq(), *get_comp_seq();
  int n_site1, n_site2, l_site1, l_site2;
  unsigned char *last_mononuc(), *first_mononuc();

  seq1 = get_seq(pair->entry1);
  seq2 = is_reverse(pair) ? get_comp_seq(pair->entry2) : get_seq(pair->entry2);
  nn = 0;
  l_site1 = n_site1 = pair->start1 - 1;
  l_site2 = n_site2 = pair->start2 - 1;
  for (diff = pair->diffs; *diff; diff++) {
    n_site1 += diff_gap1(*diff);
    n_site2 += diff_gap2(*diff);
    if (diff_type(*diff) == 'M') continue;
    if (n_site1 > l_site1 + parameters->node_seg && n_site2 > l_site2 + parameters->node_seg) { 
/* should probably adjust by 1 for type == 'I' or 'D' */
      first_disp = (last_mononuc(seq1 + l_site1) - seq1) - l_site1;
      first_disp -= 2; 
      if (first_disp < parameters->node_space) first_disp = parameters->node_space; 
      last_disp = first_mononuc(seq1 + n_site1 + (diff_type(*diff) == 'I')) - seq1;
      last_disp += 2; 
      if (last_disp > n_site1 - parameters->node_space) last_disp = n_site1 - parameters->node_space; 
      last_disp -= l_site1;
      if (last_disp >= first_disp) {
	nn += 2 * ((last_disp - first_disp) / parameters->node_space + 1);
      }
    }
    l_site1 = n_site1;
    l_site2 = n_site2;
  }
  return nn;
}

make_nodes(contig)
     Contig *contig;
{
  int i, j, k, i_node, length1, length2, entry1;
  int site1, site2;
  int *rev_ptrs;
  char *our_alloc();
  Align_info *align_entry;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  unsigned char *diff;
  int displace, first_disp, last_disp;
  unsigned char *seq1, *seq2;
  unsigned char *get_seq(), *get_comp_seq();
  int *node_ptrs; /* Node **node_ptrs; */
  Node temp_node;
  int n_site1, n_site2, l_site1, l_site2;
  unsigned char *last_mononuc(), *first_mononuc();

  node_array = (Node *)our_alloc(n_nodes * sizeof(Node));
  i_node = 0;
/*
  fprintf(stderr,"\nContig: %d", contig->index);
  notify("");
*/
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
/*
    fprintf(stderr,"%s, ", align_entry->db_entry->id);
    notify("");
*/    
    entry1 = align_entry->seq_entry;
/* nodes at beginning, end of sequence */
    node_array[i_node].entry1 = node_array[i_node + 1].entry1 = entry1;
    node_array[i_node].link = node_array[i_node + 1].link = BIGNEG;
    node_array[i_node].best = node_array[i_node + 1].best = BIGNEG;
    node_array[i_node].score = node_array[i_node + 1].score = BIGNEG;
/*    node_array[i_node].position = node_array[i_node + 1].position = BIGNEG; */
    node_array[i_node].site = 0; 
    node_array[i_node + 1].site = get_seq_length(entry1);
/*    node_array[i_node].swat_score = node_array[i_node + 1].swat_score = 0; */
    
    i_node += 2;
    length1 = get_seq_length(entry1) - 1;
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (!is_used(pair)) continue; /* restrict to pairs that were directly 
				    used in contig construction */
/*
      if (!is_best(pair)) continue;
*/
      if (entry1 >= pair->entry2) continue; 
      length2 = get_seq_length(pair->entry2) - 1;
/*
      set_nodes(i_node, pair);
      set_nodes(i_node + 2, pair);
      if (align_entry->reverse) {
	node_array[i_node].site =  length1 - pair->end1;
	node_array[i_node + 1].site = length2 - pair->end2;
	node_array[i_node + 2].site = length1 - pair->start1;
	node_array[i_node + 3].site = length2 - pair->start2;
      }
      else {
	node_array[i_node].site = pair->start1;
	node_array[i_node + 1].site = pair->start2;
	node_array[i_node + 2].site = pair->end1;
	node_array[i_node + 3].site = pair->end2;
      }
      i_node += 4;
*/
 /*     if (!(diff = pair->diffs)) continue;  or pair->w_diffs */

      seq1 = get_seq(pair->entry1);
      seq2 = is_reverse(pair) ? get_comp_seq(pair->entry2) : get_seq(pair->entry2);
      l_site1 = n_site1 = pair->start1 - 1;
      l_site2 = n_site2 = pair->start2 - 1;
      for (diff = pair->diffs; *diff; diff++) {
	n_site1 += diff_gap1(*diff);
	n_site2 += diff_gap2(*diff);
	if (diff_type(*diff) == 'M') continue;
	if (n_site1 > l_site1 + parameters->node_seg && n_site2 > l_site2 + parameters->node_seg) {
	  /* should probably adjust by 1 for type == 'I' or 'D' */
	  first_disp = (last_mononuc(seq1 + l_site1) - seq1) - l_site1;
	  first_disp -= 2; 
	  if (first_disp < parameters->node_space) first_disp = parameters->node_space; 
	  last_disp = first_mononuc(seq1 + n_site1 + (diff_type(*diff) == 'I')) - seq1;
	  last_disp += 2; 
	  if (last_disp > n_site1 - parameters->node_space) last_disp = n_site1 - parameters->node_space; 
	  last_disp -= l_site1; 
/*
	  if (last_disp < first_disp) 
	    last_disp = first_disp = (diff->site1 - (diff - 1)->site1) / 2;
*/
	  for (displace = first_disp; displace <= last_disp; displace += parameters->node_space) {
	    /*	  displace = (diff->site1 - (diff - 1)->site1) / 2;   / 4 */
	    site1 = l_site1 + displace; 
	    site2 = l_site2 + displace; 
	    set_nodes(i_node, pair);
	    node_array[i_node].site = align_entry->reverse ? length1 - site1 : site1;
	    node_array[i_node + 1].site = align_entry->reverse ? length2 - site2 : site2;
	    i_node += 2;
	  }
	  /*
	    site1 = (diff - 1)->site1 + 3 * displace;
	    site2 = (diff - 1)->site2 + 3 * displace;
	    set_nodes(i_node, pair);
	    if (align_entry->reverse) {
	    node_array[i_node].site =  length1 - site1;
	    node_array[i_node + 1].site = length2 - site2;
	    }
	    else {
	    node_array[i_node].site =  site1;
	    node_array[i_node + 1].site = site2;
	    }
	    i_node += 2;
	    */
	}
	l_site1 = n_site1;
	l_site2 = n_site2;
      }
    }
    /*
      if (num_links > 2) 
      fprintf(fp_log,"\n%-15s  %d links", align_entry->db_entry->id, num_links);
      */
  }
  if (i_node != n_nodes) {
    fprintf(stderr,"\nERROR in no. nodes: %d vs. %d\n", i_node, n_nodes);
    fatalError("");
  }
  n_nodes = i_node;
  node_ptrs = (int *)our_alloc(n_nodes * sizeof(int));
  for (i = 0; i < n_nodes; i++) node_ptrs[i] = i;

  qsort(node_ptrs, n_nodes, sizeof(int), compare_nodes);
 
  rev_ptrs = (int *)our_alloc(n_nodes * sizeof(int));
  for (i = 0; i < n_nodes; i++) rev_ptrs[node_ptrs[i]] = i;
  for (i = 0; i < n_nodes; i++) 
    if (node_array[i].link > BIGNEG)
      node_array[i].link = rev_ptrs[node_array[i].link];
  our_free(rev_ptrs);

  for (i = 0; i < n_nodes; i++) { 
    if (node_ptrs[i] < 0) continue; /* already filled */
    copy_node(&temp_node, &node_array[i]);
    for (j = i; node_ptrs[j] != i; j = k) {
      k = node_ptrs[j];
      copy_node(&node_array[j], &node_array[k]);
      node_ptrs[j] = -1;
    } 
    copy_node(&node_array[j], &temp_node);
    node_ptrs[j] = -1;
  }
  our_free(node_ptrs);
}

set_nodes(i_node, pair)
     int i_node;
     Aligned_pair *pair;
{
  node_array[i_node].entry1 = pair->entry1;
  node_array[i_node].link = i_node + 1;

  node_array[i_node + 1].entry1 = pair->entry2;
  node_array[i_node + 1].link = i_node;

  node_array[i_node].best = node_array[i_node + 1].best = BIGNEG;
  node_array[i_node].score = node_array[i_node + 1].score = BIGNEG;
/*
  node_array[i_node].swat_score = node_array[i_node + 1].swat_score 
    = pair->score;
*/
/*  node_array[i_node].position = node_array[i_node + 1].position = BIGNEG;*/
/*  node_array[i_node].used = node_array[i_node + 2].used = 0; */
}

/* Following is based on Tarjan's algorithm to find strongly connected 
   components of a directed graph -- as implemented in Sedgewick; we use the
   fact that a strongly connected component is pulled off the stack (and
   given scores) before any ancestor of any node in it is */

int tarjan(k)
     int k;
{
  int m, min; 

  num_visits++;
  min = val[k] = ++id;
  stack[p++] = k;
  if (!(p % 1000) && p > max_p) {
    fprintf(stderr, "\nstack depth: %d", p);
    max_p = p;
    fflush(stderr);
  }
  if (k > 0 && node_array[k - 1].entry1 == node_array[k].entry1 
      && node_array[k - 1].site == node_array[k].site) {
    m = val[k - 1] ? val[k - 1] : tarjan(k - 1);
    if (m < min) min = m;
  }
  if (k < n_nodes - 1 && node_array[k + 1].entry1 == node_array[k].entry1) {
    m = val[k + 1] ? val[k + 1] : tarjan(k + 1);
    if (m < min) min = m; /* if (m < min) min = m; */
  }
  if (node_array[k].link > BIGNEG) {
    m = val[node_array[k].link] ? val[node_array[k].link] : tarjan(node_array[k].link);
    if (m < min) min = m; /* if (m < min) min = m; */
  }
  if (min == val[k]) pop_component(k);

  return min;
}


/* version of above that avoids recursive function calls */

typedef struct recursion_state {
  int k, min, calling_point;
  struct recursion_state *next, *prev;
} Recursion_state;

static Recursion_state *curr_state, *head_state;

alloc_recursion_state()
{
  char *our_alloc();
  Recursion_state *old_state;
  int i;

  if (!head_state) {
    head_state = curr_state = (Recursion_state *)our_alloc(sizeof(Recursion_state)); 
    curr_state->next = curr_state->prev = 0;
  }
  old_state = curr_state;
  curr_state = curr_state->next;
  if (!curr_state) {
    curr_state = (Recursion_state *)our_alloc(500 * sizeof(Recursion_state)); 
    old_state->next = curr_state;
    for (i = 0; i < 500; i++) {
      (curr_state + i)->prev = old_state;
      old_state = old_state->next = curr_state + i;
    }
    old_state->next = 0;
  }
}
    
set_recursion_state(k, min, calling_point)
  int k, min, calling_point;
{
  curr_state->k = k;
  curr_state->min = min;
  curr_state->calling_point = calling_point;
}
    
nr_tarjan(k) 
     int k;
{
  int m, min, call_value, calling_point; 

  call_value = k;

 code_start:
  alloc_recursion_state();
  k = call_value;

  num_visits++;
  min = val[k] = ++id;
  stack[p++] = k;
  if (!(p % 1000) && p > max_p) {
    fprintf(stderr, "\nstack depth: %d", p);
    max_p = p;
    fflush(stderr);
  }
  if (k > 0 && node_array[k - 1].entry1 == node_array[k].entry1 
      && node_array[k - 1].site == node_array[k].site) {

    call_value = k - 1;
    calling_point = 0;

    if (val[call_value]) m = val[call_value];
    else {
      set_recursion_state(k, min, calling_point);
      goto code_start;
    re_entry0: /* return here after execution */
      ;
     }

    if (m < min) min = m;
  }
  if (k < n_nodes - 1 && node_array[k + 1].entry1 == node_array[k].entry1) {

    call_value = k + 1;
    calling_point = 1;

    if (val[call_value]) m = val[call_value];
    else {
      set_recursion_state(k, min, calling_point);
      goto code_start;
    re_entry1: /* return here after execution */
      ;
     }

    if (m < min) min = m; 
  }
  if (node_array[k].link > BIGNEG) {

    call_value = node_array[k].link;
    calling_point = 2;

    if (val[call_value]) m = val[call_value];
    else {
      set_recursion_state(k, min, calling_point);
      goto code_start;
    re_entry2: /* return here after execution */
      ;
     }

    if (m < min) min = m; 
  }
  if (min == val[k]) pop_component(k);

  curr_state = curr_state->prev;
  if (curr_state == head_state) {
    return; 
  }

  m = min;

  k = curr_state->k;
  min = curr_state->min;
  calling_point = curr_state->calling_point;

  if (calling_point == 0) goto re_entry0;
  else if (calling_point == 1) goto re_entry1;
  else if (calling_point == 2) goto re_entry2;
}


    
pop_component(k)
     int k;
{
  int entry1, length;
  int p_end, p_start, q, pos;
  char loop_flag, min_flag;
  unsigned char *seq;
  char *get_id();
  unsigned char *get_seq(), *get_comp_seq();
  char *get_adj_qual(), *get_orig_qual();
  char *a_qual, *o_qual;
  int score, best_score, best_i, min_i, max_i, i, j;
  Align_info *align_entry;
  
  p_end = --p;
  while (stack[p] != k) p--;
  p_start = p;
  loop_flag = 0;
  best_score = BIGNEG;
  best_i = BIGNEG;
  for (p = p_start; p <= p_end; p++) {
    min_i = max_i = stack[p];
    val[min_i] = n_nodes + 1;
    if (node_array[min_i].score > BIGNEG) continue; /* previously processed in
						   this loop */
    entry1 = node_array[min_i].entry1;
    for (q = p + 1; q <= p_end; q++) 
      if (node_array[stack[q]].entry1 == entry1) {
	if (stack[q] < min_i) min_i = stack[q];
	if (stack[q] > max_i) max_i = stack[q];
      }
    if (node_array[min_i].site != node_array[max_i].site) {
      if (!loop_flag) {
	loop_flag = 1;
	fprintf(fp_log,"\nLOOP:");
      }
      fprintf(fp_log,"\n  %c %-15s  %4d    %3d", align_array[entry1].reverse ? 'C' : ' ',
	      get_id(entry1), node_array[min_i].site, 
	      node_array[max_i].site - node_array[min_i].site);
    }
    if (max_i < n_nodes - 1 && node_array[max_i + 1].entry1 == entry1) {
      align_entry = align_array + entry1;
      seq = align_entry->reverse ? get_comp_seq(entry1) : get_seq(entry1);
      length = get_seq_length(entry1);
      score = node_array[max_i + 1].score;
      a_qual = get_adj_qual(entry1);
      o_qual = get_orig_qual(entry1);
      for (j = node_array[min_i].site; j < node_array[max_i + 1].site; j++) { 
	/*	  score += align_array[entry1].depth[j]; 
		  if (islower(seq[j]) || seq[j] == 'N' || seq[j] == 'X') ;
		  else score++; 
		  */  
	pos = align_entry->reverse ? length - 1 - j : j;
	q = a_qual[pos];
	score += contig_graph_weights[seq[j] == 'N' && score > 10 ? -1 : 
				      (seq[j] == 'X' ? o_qual[pos] : q)];

/*
	else score += q ? (q == 99 ? EDIT_WEIGHT : q) : o_qual[pos];
*/
	/* + align_entry->orig_qual[pos]*/
      }
      if (score < 0) {
	score = 0;
	node_array[min_i].best = BIGNEG;
      }
      else node_array[min_i].best = max_i + 1;
    }
    else {
      score = 0;
      node_array[min_i].best = BIGNEG;
    }
    node_array[min_i].score = score;
    for (i = min_i + 1; i <= max_i; i++) { 
      node_array[i].score = score;
      node_array[i].best = max_i + 1;
    }
    if (score > best_score) {
      best_score = score;
      best_i = min_i;
    }
  }
  if (!loop_flag) { /* don't allow jumps to another read if a loop exists */
    for (p = p_start; p <= p_end; p++) {
      if (node_array[stack[p]].score < best_score) { 
	node_array[stack[p]].score = best_score;
	node_array[stack[p]].best = best_i;
/*
	if (!strcmp(get_id(entry1), "iv80d05.s1")) {
	  fprintf(stderr, "iv80d05.s1 HERE1 %d ", best_score);
	}
*/
      }
    }
  }
  else { /* allow jumps from first position in group only */
    for (p = p_start; p <= p_end; p++) {
      entry1 = node_array[stack[p]].entry1;
      min_flag = 1;
      for (q = p_start; q <= p_end; q++) 
	if (node_array[stack[q]].entry1 == entry1 && node_array[stack[q]].site < node_array[stack[p]].site) {
	  min_flag = 0;
	  break;
	}
      if (min_flag && node_array[stack[p]].score < best_score) {
/*
	if (!strcmp(get_id(entry1), "iv80d05.s1")) {
	  fprintf(stderr, "iv80d05.s1 HERE2 %d ", best_score);
	}
*/
	node_array[stack[p]].score = best_score;
	node_array[stack[p]].best = best_i;
      }
    }
  }

  p = p_start;
  n_components++;
}
    
typedef struct c_node {
  int entry, entry_pos, contig_pos;
  int contig_jump, entry_jump;
} C_node;

typedef struct contig_place {
  int best_c_node;
  int best_score;
} Contig_place;

compare_c_nodes(c_node1, c_node2)
     C_node *c_node1, *c_node2;
{
  return c_node2->contig_pos - c_node1->contig_pos;
}


static int comp[256];
static int index1[256];
static int **insert_score;

set_vectors(length)
     int length;
{
  int i, j;
  char *our_alloc();

  for (i = 0; i < 256; i++) comp[i] = index1[i] = 0;
  comp['A'] = 'T';
  comp['C'] = 'G';
  comp['G'] = 'C';
  comp['T'] = 'A';
  index1['A'] = 0;
  index1['C'] = 1;
  index1['G'] = 2;
  index1['T'] = 3;

  free_vectors();

  insert_score = (int **)our_alloc(4 * sizeof(int *));
  for (i = 0; i < 4; i++) {
    insert_score[i] = (int *)our_alloc(length * sizeof(int));
    for (j = 0; j < length; j++) insert_score[i][j] = 0;
  }
}

free_vectors()
{
  int i;

  if (insert_score) {
    for (i = 0; i < 4; i++) our_free(insert_score[i]);
    our_free(insert_score);
    insert_score = 0;
  }
}

#define MAX_SCORE_STACK 500
static int n_old_base_segments, n_new_base_segments, old_length, new_length;

/* NB BREAKS IN COVERAGE BY READ ALIGNMENTS WILL CAUSE BREAKS IN CONTIGS. 
   CORRECT THIS BY ALLOWING ORIGINAL CONTIG ITSELF AS A READ?
*/
contig_revise(contig)
     Contig *contig;
{
  unsigned char *seq;
  Align_info *align_entry;
  int length1, length;
  unsigned char *diff;
  int i, j, i_c, i_node, entry1;
  int d_site1, d_site2, p_site1, p_site2, u_site2;
  int d_gap1, d_gap2, type;
  char *our_alloc();
  char *get_id();
  unsigned char *get_comp_seq(), *get_seq();
  char *get_adj_qual(), *get_orig_qual();
  int min_seg, q;
  char *entry_qual, *e_o_qual, *e_a_qual;
  int n_c_nodes, best_score, best_i, score, last_entry;
  C_node *c_nodes;
  Contig_place *contig_places;
  Base_segment *base_segment, *next;
  int o_qual, best_o_qual, shortest_jump, jump, ins_score, ins_score2, score_add, best_avg, inset;
  int use_orig, insert_size, c;
  int displace, max_displace, min_displace, t_displace, t_max_displace, t_min_displace;
  int compress_adjust, j1, k1, reverse;
  int score_stack[MAX_SCORE_STACK];

  set_vectors(contig->length); /* should actually do somewhere else */
  cand_compressions(contig->seq, contig->length);
  compress_adjust = contig_graph_weights[20]; /* MAY BE INAPPROPRIATE IF -contig_graph_weights > 0 */
  use_orig = 0; /* use use_orig = 1 uses orig_qual, rather than adj_qual -- 
		   this should instead be provided externally */
  min_seg = 3; /* minimum size of flanking perfect matching segments */
  inset = 1; /* number of bases to include at start and end of matching segments */
  if (min_seg < 2 * inset - 1) fatalError("Min_seg / inset inconsistency (nodes.c)");
  n_c_nodes = 0;
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
/*     if (align_entry->LLR_score < 0) continue; do not revise using a
						 possibly misplaced read */
    d_site2 = align_entry->m_start - 2;
    d_site1 = align_entry->start + align_entry->m_start - 3;

    for (diff = align_entry->diffs - 1; *(diff + 1); ) {
      get_next_diff(&diff, &d_site1, &d_site2, &d_gap1, &d_gap2, &type);
      if (d_gap1 > d_gap2) d_gap2 = d_gap1;
      if (d_gap2 > min_seg) n_c_nodes += d_gap2 - 1 - 2* (inset - 1);
    }
  }
  c_nodes = (C_node *)our_alloc(n_c_nodes * sizeof(C_node));
  i_node = 0;
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
/*     if (align_entry->LLR_score < 0) continue; do not revise using a
						 possibly misplaced read */
    entry1 = align_entry->seq_entry;
    seq = align_entry->reverse ? get_comp_seq(entry1) : get_seq(entry1);
    length1 = get_seq_length(entry1) - 1;
    p_site2 = d_site2 = align_entry->m_start - 2;
    p_site1 = d_site1 = align_entry->start + align_entry->m_start - 3;
    u_site2 = -1;
    for (diff = align_entry->diffs - 1; *(diff + 1); ) {
      get_next_diff(&diff, &d_site1, &d_site2, &d_gap1, &d_gap2, &type);
      if (d_gap1 > d_gap2) d_gap2 = d_gap1;
      if (d_gap2 > min_seg) {
	if (u_site2 > -1) {
	  c_nodes[i_node - 1].entry_jump = p_site2 + inset + 1 - c_nodes[i_node - 1].entry_pos; 
/* try decreasing by 1 */
	  c_nodes[i_node - 1].contig_jump = p_site1 + inset + 1 - c_nodes[i_node - 1].contig_pos;
	}
	for (i = p_site2 + inset, j = p_site1 + inset; 
	     i <= d_site2 - inset || j <= d_site1 - inset; i++, j++) {
	  c_nodes[i_node].entry = entry1;
	  c_nodes[i_node].entry_pos = i;
	  c_nodes[i_node].contig_pos = j;
	  c_nodes[i_node].entry_jump = 1;
	  c_nodes[i_node].contig_jump = 1;
	  i_node++;
	}
	u_site2 = i;
      }
      p_site2 = d_site2;
      p_site1 = d_site1;
    }
  }
  if (i_node != n_c_nodes) fatalError("Node count");
  qsort(c_nodes, n_c_nodes, sizeof(C_node), compare_c_nodes);
  contig_places = (Contig_place *)our_alloc((contig->length + 1) * sizeof(Contig_place));
  for (i = 0; i <= contig->length; i++) {
    contig_places[i].best_score = contig_places[i].best_c_node = -1;
  }

  for (i = 0; i < n_c_nodes; i++) {
    if (!i || c_nodes[i].contig_pos != c_nodes[i - 1].contig_pos) {
      best_o_qual = 0;
      shortest_jump = 10000;
      best_score = best_avg = -1;
    }
    align_entry = align_array + c_nodes[i].entry;
    entry1 = c_nodes[i].entry;
    reverse = align_entry->reverse;
    entry_qual = use_orig ? get_orig_qual(entry1) : get_adj_qual(entry1);
    e_o_qual = get_orig_qual(entry1);
    seq = reverse ?  get_comp_seq(entry1) : get_seq(entry1);
    length1 = get_seq_length(entry1) - 1;
    jump = c_nodes[i].entry_jump;
    if (jump >= MAX_SCORE_STACK) continue;
    insert_size = jump - c_nodes[i].contig_jump;
    if (insert_size < 0) insert_size = 0;
/*    best_score = contig_places[c_nodes[i].contig_pos].best_score; */

    o_qual = 0;
    for (j1 = 0; j1 < jump; j1++) {
      j = c_nodes[i].entry_pos + j1;
      o_qual += e_o_qual[reverse ? length1 - j : j];
      q = entry_qual[reverse ? length1 - j : j];
      score_add =  contig_graph_weights[seq[j] == 'N' ? -1 : q];
/*
      if (align_entry->LLR_score < 0 && q != EDIT_WEIGHT)
	score_add = score_add > 20 ? score_add - 20 : 0;
*/
/* downweight negative scoring reads -- more likely to be contaminants */
      for (k1 = j1 - 1; k1 >= 0 && score_add > score_stack[k1]; k1--) 
	score_stack[k1 + 1] = score_stack[k1];
      score_stack[k1 + 1] = score_add;
    }
    for (j1 = ins_score = 0; j1 < jump - insert_size; j1++) {
      ins_score += score_stack[j1];
    }
    for (ins_score2 = 0; j1 < jump; j1++) {
      ins_score2 += score_stack[j1];
    }

    j = c_nodes[i].entry_pos + inset;
    if (insert_size == 1 && 
	insert_score[index1[seq[j]]][c_nodes[i].contig_pos + inset - 1] >= 4) {
/*
      j = c_nodes[i].entry_pos - 1;
      j1 = c_nodes[i].entry_pos + jump;
      if (ins_score || 	entry_qual[reverse ? length1 - j : j]
	  || entry_qual[reverse ? length1 - j1 : j1])
*/
      ins_score += ins_score2 + compress_adjust;    
/*      ins_score += compress_adjust;  */
    }
/*     ins_score += ins_score2; readded 5/2/98 -- ???  */  

    score = contig_places[c_nodes[i].contig_pos + c_nodes[i].contig_jump].best_score + ins_score;
    if (score < 0) score = 0;
/*    if (jump <= c_nodes[i].contig_jump || jump < 3 ||  ins_score / (jump - 2) >= 20 ) { */
/*    if (ins_score >= best_avg) */
    if (score > best_score
	|| score == best_score && c_nodes[i].entry_jump < shortest_jump
	|| score == best_score && c_nodes[i].entry_jump == shortest_jump && o_qual >= best_o_qual
	) {
      contig_places[c_nodes[i].contig_pos].best_score = best_score = score;
      contig_places[c_nodes[i].contig_pos].best_c_node = i;
      best_o_qual = o_qual;
      shortest_jump = c_nodes[i].entry_jump;
      best_avg = ins_score;
    }
  }

  best_score = best_i = 0;
  for (i = contig->length; i >= 0; i--) {
    if (contig_places[i].best_score >= best_score) {
      best_score = contig_places[i].best_score;
      best_i = i;
    }
  }

  length = 0;
  for (i = best_i; (i_node = contig_places[i].best_c_node) > -1; 
       i = c_nodes[i_node].contig_pos + c_nodes[i_node].contig_jump) {
    length += c_nodes[i_node].entry_jump;
  }
/*
  fprintf(stderr,"\n\nContig length: old %d, new %d", contig->length, length);
  fprintf(stderr, "\nNew start: %d \n", best_i);
  for (i = 0; i < best_i; i++)
    fprintf(stderr, "%d ", contig->adj_qual[i]);
  for (i = 0; i < best_i; i++)
    fprintf(stderr, "%c", contig->seq[i]);
  fprintf(stderr, "\n\n");
  for (i = best_i; i < best_i + 100; i++)
    fprintf(stderr, "%d ", contig->adj_qual[i]);
  for (i = best_i; i < best_i + 100; i++)
    fprintf(stderr, "%c", contig->seq[i]);
*/

  fprintf(fp_log,"\n\nContig length: old %d, new %d", contig->length, length);
  fprintf(fp_log, "\nNew start: %d ", best_i);
/*
  fprintf(stderr,"\nContig length: old %d, new %d ", contig->length, length); 
  fprintf(stderr," %d  %d ", best_i, best_score); 
*/
  our_free(contig->seq);
  our_free(contig->adj_qual);
  our_free(contig->orig_qual);
  old_length += contig->length;
  new_length += length;
  for (base_segment = contig->base_segment; base_segment; base_segment = next) {
    n_old_base_segments++;
    next = base_segment->next;
    our_free(base_segment);
  }
  contig->base_segment = 0;

  contig->seq = (unsigned char *)our_alloc((length + 1) * sizeof(unsigned char));
  contig->adj_qual = (char *)our_alloc((length + 1) * sizeof(char));
  contig->orig_qual = (char *)our_alloc((length + 1) * sizeof(char));
  contig->length = length;
  i_c = 0;
  last_entry = -1;
  t_min_displace = t_max_displace = t_displace = best_i;
  for (i = best_i; (i_node = contig_places[i].best_c_node) > -1; 
       i = c_nodes[i_node].contig_pos + c_nodes[i_node].contig_jump) {
    if (last_entry != c_nodes[i_node].entry) {
      last_entry = c_nodes[i_node].entry;
      n_new_base_segments++;
      append_base_segment(contig, last_entry, c_nodes[i_node].entry_pos + 1, i_c + 1);
      align_entry = last_entry + align_array;
      reverse = align_entry->reverse;
      length1 = get_seq_length(last_entry) - 1;
      seq = reverse ?  get_comp_seq(last_entry) : get_seq(last_entry);
      e_o_qual = get_orig_qual(last_entry);
      e_a_qual = get_adj_qual(last_entry);
    }
    jump = c_nodes[i_node].entry_jump;
    if (jump > 1/* || jump < c_nodes[i_node].contig_jump*/) {
      if (jump < c_nodes[i_node].contig_jump) {
	c = 'D';
	t_displace += c_nodes[i_node].contig_jump - jump;
	if (t_displace > t_max_displace) t_max_displace = t_displace;
      }
      else if (jump > c_nodes[i_node].contig_jump) {
	c = 'I';
	t_displace += c_nodes[i_node].contig_jump - jump;
	if (t_displace < t_min_displace) t_min_displace = t_displace;
      }
	
      else c = 'S';
      fprintf(fp_log, "\n%c Change at: %d; seg lengths: entry %d, contig %d", 
	      c, i_c + 2, jump - 2 * inset, 
	      c_nodes[i_node].contig_jump - 2 * inset);
    }
    for (j = c_nodes[i_node].entry_pos; 
	 j < c_nodes[i_node].entry_pos + jump; j++, i_c++) {
      contig->seq[i_c] = seq[j];
      contig->adj_qual[i_c] = e_a_qual[reverse ? length1 - j : j];
      contig->orig_qual[i_c] = e_o_qual[reverse ? length1 - j : j];
    }
    contig->base_segment->read_end = j;
    contig->base_segment->contig_end = i_c;
  }
  if (i_c != length) fatalError("Length discrepancy");
  contig->seq[i_c] = 0;
  contig->score = best_score;
/*
  fprintf(stderr, "\n\n");
  for (i = 0; !contig->adj_qual[i] && i < 300; i++)
    fprintf(stderr, "%d ", contig->adj_qual[i]);
  for (i = 0; !contig->adj_qual[i] && i < 300; i++)
    fprintf(stderr, "%c", contig->seq[i]);
  fprintf(stderr, "\n%d\n", i);
  for (j = i; j < 300; j++)
    fprintf(stderr, "%d ", contig->adj_qual[j]);
  for (j = i; j < 300; j++)
    fprintf(stderr, "%c", contig->seq[j]);
*/
  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    displace = max_displace = min_displace = 0;
    for (diff = align_entry->diffs; *(diff + 1); diff++) {
      if (diff_type(*diff) == 'I') {
	displace--;
	if (displace < min_displace) min_displace = displace;
      }
      else if (diff_type(*diff) == 'D') {
	displace++;
	if (displace > max_displace) max_displace = displace;
      }
    }
    displace = align_entry->start;
    align_entry->start = displace + min_displace - t_max_displace;
    align_entry->end = displace + max_displace + get_seq_length(align_entry->seq_entry) - 1 - t_min_displace;
  }

  free_vectors();
  our_free(c_nodes);
  our_free(contig_places);
}

count_base_segments()
{
  fprintf(stderr, " %d old base segments, %.3f Mb\n", 
	  n_old_base_segments, n_old_base_segments * (sizeof(Base_segment) / 1000000.));
  fprintf(stderr, " %d new base segments, %.3f Mb\n", 
	  n_new_base_segments, n_new_base_segments * (sizeof(Base_segment) / 1000000.));
  fprintf(stderr, " total old contig length %d, %.3f Mb\n", 
	 old_length, 3. * old_length / 1000000.);
  fprintf(stderr, " total new contig length %d, %.3f Mb\n", 
	 new_length, 3. * new_length / 1000000.);
}

/*POSSIBLE PROBLEMS WITH FOLLOWING: INSERT_SCORE DOES NOT DISTINGUISH STRAND 
-- AND CODE TO ALLOW FOR DISPLACED POSITION FOR COMPRESSION MAY NOT BE CORRECT.
COMPARE TO CAND_COMPRESSION_MOTIF WHICH MAY BE BETTER */

#define MAX_STEM_LENGTH 50
cand_compressions(seq, length)
     unsigned char *seq;
     int length;
{
  int i, j, k, m, loop_size, loop_score, stem_score, score, id;
  int c, stem_size;
  char stem_seq[MAX_STEM_LENGTH];

  for (i = 8; i < length; i++) {
    for (loop_size = 3; loop_size <= 6; loop_size++) {
      loop_score = 0;
      if (loop_size <= 4 && seq[i - loop_size] == 'G' && seq[i - 1] == 'A'
	  || loop_size == 4 && seq[i - loop_size] == 'C' && seq[i - 1] == 'G')
	loop_score += 2;
      if (loop_size == 3) loop_score++;

      k = i - loop_size - 1;

      stem_seq[0] = comp[seq[k]];
      stem_size = 1;

      for (j = i; j < length - 1 && k > 0; j++, k--) {
 	stem_seq[stem_size++] = comp[seq[k - 1]];
	if (stem_size >= MAX_STEM_LENGTH) break;
	stem_seq[stem_size] = 0;

	stem_score = find_stem_score(stem_seq);
	score = stem_score + loop_score;

	if (seq[k] == comp[seq[j]] && stem_score) {
	  c = comp[seq[k - 1]];
	  id = index1[c];
	  if (insert_score[id][j] < score) insert_score[id][j] = score;
	  for (m = j + 1; m < length && seq[m] == c; m++)
	    if (insert_score[id][m] < score) insert_score[id][m] = score;
	  for (m = j - 1; m >= 0 && seq[m + 1] == c; m--)
	    if (insert_score[id][m] < score) insert_score[id][m] = score;
	}
	if (seq[k - 1] == comp[seq[j]] && stem_score) {
	  c = comp[seq[k]];
	  id = index1[c];
	  if (insert_score[id][j - 1] < score) insert_score[id][j - 1] = score;
	  for (m = j; m < length && seq[m] == c; m++)
	    if (insert_score[id][m] < score) insert_score[id][m] = score;
	  for (m = j - 2; m >= 0 && seq[m + 1] == c; m--)
	    if (insert_score[id][m] < score) insert_score[id][m] = score;
	}
	if (seq[k] != comp[seq[j]]) break;
      }
    }
  }
  for (i = 0; i < length - 8; i++) {
    for (loop_size = 3; loop_size <= 6; loop_size++) {
      loop_score = 0;
      if (loop_size <= 4 && seq[i + loop_size] == 'C' && seq[i + 1] == 'T'
	  || loop_size == 4 && seq[i + loop_size] == 'G' && seq[i + 1] == 'C')
	loop_score += 2;
      if (loop_size == 3) loop_score++;

      k = i + loop_size + 1;

      stem_seq[0] = seq[k];
      stem_size = 1;

      for (j = i; j > 0 && k < length - 1; j--, k++) {
 	stem_seq[stem_size++] = seq[k + 1];
	if (stem_size >= MAX_STEM_LENGTH) break;
	stem_seq[stem_size] = 0;

	stem_score = find_stem_score(stem_seq);
	score = stem_score + loop_score;

	if (seq[k] == comp[seq[j]] && stem_score) {
	  c = comp[seq[k + 1]];
	  id = index1[c];
	  if (insert_score[id][j - 1] < score) insert_score[id][j - 1] = score;
	  for (m = j; m < length && seq[m] == c; m++)
	    if (insert_score[id][m] < score) insert_score[id][m] = score;
	  for (m = j - 2; m >= 0 && seq[m + 1] == c; m--)
	    if (insert_score[id][m] < score) insert_score[id][m] = score;
	}
	if (seq[k + 1] == comp[seq[j]] && stem_score) {
	  c = comp[seq[k]];
	  id = index1[c];
	  if (insert_score[id][j] < score) insert_score[id][j] = score;
	  for (m = j + 1; m < length && seq[m] == c; m++)
	    if (insert_score[id][m] < score) insert_score[id][m] = score;
	  for (m = j - 1; m >= 0 && seq[m + 1] == c; m--)
	    if (insert_score[id][m] < score) insert_score[id][m] = score;
	}
	if (seq[k] != comp[seq[j]]) break;
      }
    }
  }
}

/* following differs from above in that sequence being scanned is assumed to
   be correct -- so looking for exact occurrence of motif; and need to
   find potential locations for deleted bases */

int **cand_compression_motifs(seq, length)
     unsigned char *seq;
     int length;
{
  int i, j, k, m, loop_size, loop_score, stem_score, score, id;
  int c, stem_size;
  char stem_seq[MAX_STEM_LENGTH];
  int **delete_score;

  set_vectors(length);
  delete_score = insert_score;

  for (i = 8; i < length; i++) {
    for (loop_size = 3; loop_size <= 6; loop_size++) {
      loop_score = 0;
      if (loop_size <= 4 && seq[i - loop_size] == 'G' && seq[i - 1] == 'A'
	  || loop_size == 4 && seq[i - loop_size] == 'C' && seq[i - 1] == 'G')
	loop_score += 2;
      if (loop_size == 3) loop_score++;

      k = i - loop_size - 1;

      stem_size = 0;

      for (j = i; j < length - 1 && k > 0 && seq[j] == comp[seq[k]]; 
	   j++, k--) {
 	stem_seq[stem_size++] = seq[j];
	if (stem_size >= MAX_STEM_LENGTH) break;
	stem_seq[stem_size] = 0;

	score = 0;
	if (stem_size >= 2 && (loop_size <= 4 || stem_size >= 3)) {
	  stem_score = find_stem_score(stem_seq);
	  score = stem_score + loop_score;
	}

	if (score) {
	  c = seq[j];
	  if (delete_score[0][j] < score) delete_score[0][j] = score;
	  for (m = j + 1; m < length && seq[m] == c; m++)
	    if (delete_score[0][m] < score) delete_score[0][m] = score;
	  for (m = j - 1; m >= 0 && seq[m] == c; m--)
	    if (delete_score[0][m] < score) delete_score[0][m] = score;
	  c = seq[j - 1];
	  if (delete_score[0][j - 1] < score) delete_score[0][j - 1] = score;
	  for (m = j; m < length && seq[m] == c; m++)
	    if (delete_score[0][m] < score) delete_score[0][m] = score;
	  for (m = j - 2; m >= 0 && seq[m] == c; m--)
	    if (delete_score[0][m] < score) delete_score[0][m] = score;
	}
      }
    }
  }
  for (i = 0; i < length - 8; i++) {
    for (loop_size = 3; loop_size <= 6; loop_size++) {
      loop_score = 0;
      if (loop_size <= 4 && seq[i + loop_size] == 'C' && seq[i + 1] == 'T'
	  || loop_size == 4 && seq[i + loop_size] == 'G' && seq[i + 1] == 'C')
	loop_score += 2;
      if (loop_size == 3) loop_score++;

      k = i + loop_size + 1;

      stem_size = 0;

      for (j = i; j > 0 && k < length - 1 && seq[k] == comp[seq[j]]; 
	   j--, k++) {
 	stem_seq[stem_size++] = seq[k];
	if (stem_size >= MAX_STEM_LENGTH) break;
	stem_seq[stem_size] = 0;

	score = 0;
	if (stem_size >= 2 && (loop_size <= 4 || stem_size >= 3)) {
	  stem_score = find_stem_score(stem_seq);
	  score = stem_score + loop_score;
	}

	if (score) {
	  c = comp[seq[k - 1]];
	  if (delete_score[1][j + 1] < score) delete_score[1][j + 1] = score;
	  for (m = j + 2; m < length && seq[m] == c; m++)
	    if (delete_score[1][m] < score) delete_score[1][m] = score;
	  for (m = j + 1; m >= 0 && seq[m] == c; m--)
	    if (delete_score[1][m] < score) delete_score[1][m] = score;
	  c = comp[seq[k]];
	  if (delete_score[1][j] < score) delete_score[1][j] = score;
	  for (m = j + 1; m < length && seq[m] == c; m++)
	    if (delete_score[1][m] < score) delete_score[1][m] = score;
	  for (m = j - 1; m >= 0 && seq[m] == c; m--)
	    if (delete_score[1][m] < score) delete_score[1][m] = score;
	}
      }
    }
  }
  /* free_vectors(); out per Bonfield suggestions */
  return delete_score;
}

int find_stem_score(stem_seq)
     char *stem_seq;
{
  int stem_size, stem_score, num_gc, i;

  stem_size = strlen(stem_seq);
  if (stem_size < 2) return 0;
  if (stem_seq[stem_size - 1] != 'C' && stem_seq[stem_size - 1] != 'G')
    return 0; /* require terminating G or C (is sometimes violated!) */

  for (i = num_gc = 0; i < stem_size; i++)
    if (stem_seq[i] == 'C' || stem_seq[i] == 'G') num_gc++;
    
  stem_score = 0;
  if (stem_size == 2) {
    if (!strcmp(stem_seq, "AC")) 
      stem_score = 2; /* require strong loop */
    else if (!strcmp(stem_seq, "CC") || !strcmp(stem_seq, "GG")) 
      stem_score = 2; /* 4; */
    else if (!strcmp(stem_seq, "GC")) 
      stem_score = 5; /* allows arb loop */
  }
  else if (stem_size == 3) {
    if (!strcmp(stem_seq, "GCC")) 
      stem_score = 6;
    else if (num_gc == 3)
      stem_score = 4; /* allows arb. 3 base loop */
    else if (num_gc >= 2 && stem_seq[0] != 'T') 
      stem_score = 2; /* so as to require strong loop */
  }
  else if (num_gc > stem_size / 2) {
    stem_score = 5;
  }
  
  return stem_score;
}

/* identify rejected alignments.  */

find_node_rejects(db)
     Database *db;
{
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int entry1, node_rej;

  notify("Finding node rejects ... "); 

  printf("\n\nNo. of node-rejected pairs:");
  node_rej = 0;
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (pair->entry1 >= pair->entry2) continue;
      if (count_pair_nodes(pair)) continue;
      node_rej++;
/*
      printf("\n%s  %d-%d   %s  %d-%d      %d", 
	     align_entry->db_entry->id, pair->start1, pair->end1,
	     align_array[pair->entry2].db_entry->id, pair->start2, pair->end2,
	     pair->score);
*/
      set_reject_node_flag(pair, 1);
      set_reject_node_flag(pair->reversed_pair, 1);
    }
  }
  if (!node_rej) printf(" None.");
  else printf(" %d.", node_rej);
  notify(" Done\n");
}

