/*****************************************************************************
#   Copyright (C) 1994-1998 by Phil Green.                          
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
extern FILE *fp_log; /* log file */


/* In following -- information concerning a contig and its complement
   are stored in the same node, in two lists -- the 0 list refers to
   the uncomplemented contig, and the 1 list to the complemented
   contig.  All "edges" are right edges (i.e. going out of the right side
   of the contig).
   Issues:
    1. Left, right extensions not precisely determined (only to within
      gap_cutoff). As a result there may be too many edges. Make more
      precise, by checking whether implied offset and contig extent really
      imply it sticks out to left or right.
    2. Large strongly connected components may take too long to search.
      A better method than exhaustive search might be to prune list of
      edges, by considering only minimal edges (where e_1 is minimal if
      it cannot be "factored" into e_2 followed by e_3). Note factoring
      requires attention to offsets, will therefore need to store offset
      in the edge structure. Also, try sorting edges by score or # of pairs.
      Need to be careful about def of "minimal" -- e.g. in case of 
      separated direct repeats, no edge may be minimal.
        Do greedy assembly of large components.
    3. Note current method to prevent both a contig and its complement fom
      appearing in the same path (using the in_use flag, which refers to the
      combination of the contig and its complement, and setting for all
      lower strong components) is not optimal, because it could prevent finding
      the optimal path. (Maybe considering both the path and its complement
      will help get around this most of the time). 
      A better general method might therefore be to
      use the idea of 2 for exhaustive search of entire (pruned) set of nodes.
    4. Need code to allow visiting a node more than once, if its depth of
      coverage suggests it is repeated. (Use in_use, but have it take on numerical
      values).
    5. merge_master issues: force joins of selected pairs, by using merge directly
      (rather than calls to pair_merge) -- do test_merge as a check only.
      Make score threshold > 0 (say 10). Also, in reject cutoff, only subtract off
      fraction of pair->score -- this applies in call to test_merge in find_best_paths
      as well.
       Only set is_repeat in first pass (so that only higher-scoring joins have been
      made earlier).
    6. Different path scoring criteria: weight each edge by number of positive
      scoring pairs it represents.
*/

extern Contig *contig_array;
extern int t_num_entries;

typedef struct edge {
  Aligned_pair *pair;
  Contig *contig1, *contig2;
  struct edge *next;
  int reverse, n_pairs;
} Edge;

typedef struct tig_node {
  Contig *contig;
  Edge *edges[2]; /* left, right, edges */
  int val[2]; /* for use with tarjan algorithm */
  int best_score[2]; /* best score for path starting at this tig_node */
  int temp_best_score[2]; /* provisional best score for path starting at this tig_node */
  int component[2]; /* strong component */
  int best_stack_length[2]; /* length of best_stack */
  Edge **best_stack[2]; /* best path within this component, starting at current node,
			 and ending just before jump to next component */
  Edge *best_jump[2]; /* best jump to next component */
  struct tig_node *next;
  int n_allowed_visits; /* no. of times this node can be visited in a single path
			   (depends on its depth of coverage) */
} Tig_node;

Edge **edge_stack, **best_stack;
int stack_index, best_n, best_score, score, n_nodes;
Tig_node *head_node, *base_node;
double n_paths, max_n_paths;
Tig_node **tarjan_stack;
int *tarjan_stack_orient;
int visit_ctr, tarjan_stack_index, n_visited_nodes;

#define WEIGHT_SIZE 1 /* weight for contribution of contig size to score -- 10 */
#define WEIGHT_LLR 0 /* weight for contribution of LLR score to score -- 1 */

alloc_tig_node(contig)
     Contig *contig;
{
  Tig_node *tig_node;
  char *our_alloc();

  tig_node = (Tig_node *)our_alloc(sizeof(Tig_node));
  tig_node->next = head_node;
  tig_node->n_allowed_visits = 1; /* need to make dependent on depth of coverage */
  tig_node->contig = contig;
  tig_node->edges[0] = tig_node->edges[1] = 0;
  tig_node->best_stack[0] = tig_node->best_stack[1] = 0;

  reset_tig_node(tig_node);

  n_nodes++;
  head_node = tig_node;
  if (contig) contig->tig_node = tig_node;
}

reset_tig_node(tig_node)
  Tig_node *tig_node;
{
  tig_node->val[0] = tig_node->val[1] = 0;
  tig_node->component[0] = tig_node->component[1] = 0;
  tig_node->best_score[0] = tig_node->best_score[1] = 0;
  tig_node->best_jump[0] = tig_node->best_jump[1] = 0;
  tig_node->temp_best_score[0] = tig_node->temp_best_score[1] = 0;
  if (tig_node->best_stack[0]) our_free(tig_node->best_stack[0]);
  if (tig_node->best_stack[1]) our_free(tig_node->best_stack[1]);
  tig_node->best_stack[0] = tig_node->best_stack[1] = 0;
  tig_node->best_stack_length[0] = tig_node->best_stack_length[1] = 0;
}

alloc_edge(contig, orientation, pair, reverse)
     Contig *contig;
     int orientation, reverse;
     Aligned_pair *pair;
{
  Tig_node *tig_node;
  Contig *contig1, *contig2;
  Edge *edge;
  char *our_alloc();
  Align_info *get_align_entry();

  if (pair) {
    tig_node = contig->tig_node;
    contig1 = get_align_entry(pair->entry1)->contig;
    contig2 = get_align_entry(pair->entry2)->contig;
    for (edge = tig_node->edges[orientation]; edge; edge = edge->next) 
      if (edge->contig2 == contig2 && edge->reverse == reverse) {
	edge->n_pairs += 1;
	if (pair->LLR_score > edge->pair->LLR_score) 
	  edge->pair = pair;
	return;
      }
  }
  else {
    tig_node = base_node;
    contig1 = 0;
    contig2 = contig;
  }
  edge = (Edge *)our_alloc(sizeof(Edge));
  edge->contig1 = contig1;
  edge->contig2 = contig2;
  edge->pair = pair;
  edge->next = tig_node->edges[orientation];
  edge->reverse = reverse;
  edge->n_pairs = 1;
  tig_node->edges[orientation] = edge;
}

find_best_paths(gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff)
     int gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff;
{
  int i, num_entries, max_entries, tot_entries;
  Contig *contig1, *contig2;
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  char *our_alloc();
  int comp_flag, set_flag, direction, min_score, orientation;
  int left_flag, right_flag, offset, gap, discrep, chim_flag;
  /* allocate nodes */
  Tig_node *tig_node;
  Edge *edge;
  int w_offset;
  Aligned_pair *reject_pair;
  int entry1, entry2;

  alloc_tig_node((Contig *)0);
  base_node = head_node; /* base_node is "dummy", at very end of list */

  for (i = 0; i < t_num_entries; i++) {
    contig1 = contig_array + i;
    if (contig1->num_entries < 2) continue;
    if (contig1->num_entries && !is_anomalous(contig1->first)) {
      for (align_entry = contig1->first; align_entry; align_entry = align_entry->next) {
	for (pair = get_aligned_pairs(align_entry->seq_entry); pair; pair = pair->next) {
	  if (is_used(pair)) continue; 
	  if (pair->LLR_score <= LLR_join_cutoff) continue;
	  if (pair_merge_reject(pair)) continue;
	  align_entry2 = get_align_entry(pair->entry2);
/*
	  chim_flag = is_reject_chimeric(pair);
	  if (is_repeat(pair) || chim_flag) continue; 

	  if (align_entry->segments && align_entry->segments->next 
	      || align_entry2->segments && align_entry2->segments->next)
	    continue;
	  if (align_entry->first_start > align_entry->last_end 
	      || align_entry2->first_start > align_entry2->last_end) 
	    continue;
*/

	  contig2 = align_entry2->contig;
	  if (contig1 == contig2) continue;
	  if (contig2->num_entries < 2) continue;
	  if (is_reverse(pair) != (align_entry->reverse != align_entry2->reverse)) {
	    complement(contig2); /* to ensure contigs are compatible */
	    comp_flag = 1;
	  }
	  else comp_flag = 0;
	  offset = find_pair_offset(pair);
	  discrep = test_merge(contig1, contig2, offset, gap_cutoff, 
			       LLR_reject_cutoff  - pair->LLR_score / 4 , LLR_join_cutoff, 
			       &gap, &min_score, &left_flag, &right_flag, &w_offset, &reject_pair, &entry1, &entry2);
	  set_flag = 0;
	  if (contig2->tig_node) {
	    for (direction = 0; !set_flag && direction < 2; direction++)
	      for (edge = contig2->tig_node->edges[direction]; edge; edge = edge->next) {
		if (edge->pair->reversed_pair == pair && edge->reverse == comp_flag) {
		  set_flag = 1;
		  break;
		}
	      }
	  }

	  if (!discrep && (left_flag == 2 || right_flag == 2)) {
	    if (!contig1->tig_node) {
	      alloc_tig_node(contig1);
	      alloc_edge(contig1, 0, (Aligned_pair *)0, 0); /* edges from base node to this one */
	      alloc_edge(contig1, 1, (Aligned_pair *)0, 0);
	    }
	    if (left_flag == 2) alloc_edge(contig1, 1, pair, comp_flag);
	    if (right_flag == 2) alloc_edge(contig1, 0, pair, comp_flag);
	  }
	  else if (set_flag) {
	    fprintf(stderr, "\nDISCREPANCY: %d/%d  %d %d %d %d", 
		    contig1->num_entries, contig2->num_entries, 
		    min_score, gap, left_flag, right_flag);
	    fprintf(fp_log, "\nDISCREPANCY: %d/%d  %d %d %d %d", 
		    contig1->num_entries, contig2->num_entries, 
		    min_score, gap, left_flag, right_flag);
	    offset = find_pair_offset(pair->reversed_pair);
	    discrep = test_merge(contig2, contig1, offset, gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff, 
				 &gap, &min_score, &left_flag, &right_flag, &w_offset, &reject_pair, &entry1, &entry2);
	    fprintf(stderr, " :: %d %d %d %d", min_score, gap, left_flag, right_flag);  
	    fprintf(fp_log, " :: %d %d %d %d", min_score, gap, left_flag, right_flag);  
	    complement(contig1);
	    complement(contig2);
	    offset = find_pair_offset(pair->reversed_pair);
	    discrep = test_merge(contig2, contig1, offset, gap_cutoff, LLR_reject_cutoff, LLR_join_cutoff, 
				 &gap, &min_score, &left_flag, &right_flag, &w_offset, &reject_pair, &entry1, &entry2);
	    fprintf(stderr, " :: %d %d %d %d", min_score, gap, left_flag, right_flag);  
	    fprintf(fp_log, " :: %d %d %d %d", min_score, gap, left_flag, right_flag);  
	    complement(contig1);
	    complement(contig2);
	    if (!contig1->tig_node) {
	      alloc_tig_node(contig1);
	      alloc_edge(contig1, 0, (Aligned_pair *)0, 0); /* edges from base node to this one */
	      alloc_edge(contig1, 1, (Aligned_pair *)0, 0);
	    }
	  }
	  if (comp_flag) complement(contig2);
	}
      }
    }
  }
  /* allocate edges */
  /* find paths -- need to consider both orientations for each contig; then accept  */

  edge_stack = (Edge **)our_alloc(2 * n_nodes * sizeof(Edge *));
/*
  best_stack = (Edge **)our_alloc(n_nodes * sizeof(Edge *));
*/
  tarjan_stack = (Tig_node **)our_alloc(2 * n_nodes * sizeof(Tig_node *));
  tarjan_stack_orient = (int *)our_alloc(2 * n_nodes * sizeof(int));

  for (i = 0; i < 2 * n_nodes; i++) {
    tarjan_stack[i] = 0;
    tarjan_stack_orient[i] = 0;
  }
  do {
    tarjan_stack_index = visit_ctr = best_score = 0;
    
    tarjan_visit(base_node, 0);
    tarjan_visit(base_node, 1);
    
    orientation = base_node->best_score[1] > base_node->best_score[0]; 
    score = base_node->best_score[orientation];
    if (score) {
      fprintf(stderr, "\n\nScore: %d ", score);
      fprintf(fp_log, "\n\nScore: %d ", score);
      set_in_use(base_node, orientation, -1, 1);
      for (tig_node = head_node; tig_node; tig_node = tig_node->next) {
	reset_tig_node(tig_node);
      }
      fprintf(stderr, "\n");
      fprintf(fp_log, "\n");
    }
  } while (score);

  return;

  fprintf(stderr, "\n%d nodes\n", n_nodes);
  fprintf(fp_log, "\n%d nodes\n", n_nodes);
  max_entries = tot_entries = 0;
  for (tig_node = head_node; tig_node != base_node; tig_node = tig_node->next) {
    fprintf(stderr, "\n%d ", tig_node->contig->num_entries);
    fprintf(fp_log, "\n%d ", tig_node->contig->num_entries);
    if (tig_node->contig->num_entries > max_entries) 
      max_entries = tig_node->contig->num_entries;
    tot_entries += tig_node->contig->num_entries;
    for (edge = tig_node->edges[0], i = best_score = 0; edge; edge = edge->next) {
      i++;
      if (edge->pair && edge->pair->LLR_score > best_score) 
	best_score = edge->pair->LLR_score;
    }  
    fprintf(stderr, "%d [%d] ", i, best_score);
    fprintf(fp_log, "%d [%d] ", i, best_score);
    for (edge = tig_node->edges[1], i = best_score = 0; edge; edge = edge->next) {
      i++;
      if (edge->pair && edge->pair->LLR_score > best_score) 
	best_score = edge->pair->LLR_score;
    }  
    fprintf(stderr, "%d [%d] ", i, best_score);
    fprintf(fp_log, "%d [%d] ", i, best_score);
  }
  fprintf(stderr, "\nMax. no. entries: %d; total no. entries: %d\n", 
	  max_entries, tot_entries);

  fprintf(fp_log, "\nMax. no. entries: %d; total no. entries: %d\n", 
	  max_entries, tot_entries);

  do {
    stack_index = best_score = score = best_n = 0;
    visit_tig_nodes(base_node, 0);
    stack_index = score = 0;
    visit_tig_nodes(base_node, 1);
    if (best_n) {
      fprintf(stderr, "\n%d edges, score %d  ", best_n - 1, best_score);
      fprintf(fp_log, "\n%d edges, score %d  ", best_n - 1, best_score);
      for (i = num_entries = 0; i < best_n; i++) {
	fprintf(stderr, "%d  ", best_stack[i]->contig2->num_entries);
	fprintf(fp_log, "%d  ", best_stack[i]->contig2->num_entries);
	best_stack[i]->contig2->tig_node->n_allowed_visits -= 1; 
	num_entries += best_stack[i]->contig2->num_entries;
	if (i) {
	  set_used_flag(best_stack[i]->pair, 1);
	  set_used_flag(best_stack[i]->pair->reversed_pair, 1);
	  fprintf(stderr, "[%d] - ", best_stack[i]->pair->LLR_score);
	  fprintf(fp_log, "[%d] - ", best_stack[i]->pair->LLR_score);
	}
      }
      fprintf(stderr, "  n_reads %d ", num_entries);
      fprintf(fp_log, "  n_reads %d ", num_entries);
    }
  } while (best_score);

  our_free(edge_stack);
/*
  our_free(best_stack);
*/
}

tarjan_visit(tig_node, orientation)
     Tig_node *tig_node;
     int orientation;
{
  int m, min;
  Tig_node *tig_node2;
  int orient2;
  Edge *edge;

  min = tig_node->val[orientation] = ++visit_ctr;
/*
  fprintf(stderr, "\n%d %d", tig_node->val[orientation], tarjan_stack_index);
*/
  tarjan_stack[tarjan_stack_index] = tig_node;
  tarjan_stack_orient[tarjan_stack_index] = orientation;
  tarjan_stack_index++;
  for (edge = tig_node->edges[orientation]; edge; edge = edge->next) {
    tig_node2 = edge->contig2->tig_node;
    if (tig_node2->n_allowed_visits > 0) {
      orient2 = orientation != edge->reverse;
      m = tig_node2->val[orient2] ? tig_node2->val[orient2] : 
	tarjan_visit(tig_node2, orient2);
      if (m < min) min = m;
    }
  }
/*
  fprintf(stderr, "\n %d %d", min, tig_node->val[orientation]);
*/
  if (min == tig_node->val[orientation])
    pop_strong_component(tig_node, orientation);

  return min;
}

pop_strong_component(tig_node, orientation)
     Tig_node *tig_node;
     int orientation;
{
  int stack_end, s, component, orient2, orient3;
  Tig_node *tig_node2, *tig_node3;
  int i, num_entries, component_size, t_score;
  Contig *contig2;
  Edge *edge, *edge1;
  char *our_alloc();

  component = tig_node->val[orientation];

  stack_end = --tarjan_stack_index;

  while (tarjan_stack[tarjan_stack_index] != tig_node || 
	 tarjan_stack_orient[tarjan_stack_index] != orientation)
    tarjan_stack_index--;
  component_size = stack_end - tarjan_stack_index + 1;
  fprintf(stderr,"\nComponent %d, %d members: ", 
	  component, component_size);
  fprintf(fp_log,"\nComponent %d, %d members: ", 
	  component, component_size);
  for (s = tarjan_stack_index; s <= stack_end; s++) {
    tig_node2 = tarjan_stack[s];
    orient2 = tarjan_stack_orient[s];
    contig2 = tig_node2->contig;
    tig_node2->val[orient2] = n_nodes + n_nodes + 2;
    tig_node2->component[orient2] = component;

    tig_node2->best_stack[orient2] = 
      (Edge **)our_alloc(component_size * sizeof(Edge *));

    tig_node2->best_score[orient2] = WEIGHT_SIZE * (contig2 ? contig2->num_entries : 0);

    for (edge = tig_node2->edges[orient2]; edge; edge = edge->next) {
      tig_node3 = edge->contig2->tig_node;
      orient3 = orient2 != edge->reverse;
      if (!tig_node3->component[orient3] 
	  || tig_node3->component[orient3] == component)
	continue;
      t_score = WEIGHT_LLR * (edge->pair ? edge->pair->LLR_score : 0) 
	+ tig_node3->best_score[orient3];
      if (t_score > tig_node2->temp_best_score[orient2]) {
	edge1 = tig_node2->best_jump[orient2];
	tig_node2->best_jump[orient2] = edge;
	set_in_use(tig_node2, orient2, -1, 0); 
	if (tig_node2->n_allowed_visits) {
	  tig_node2->temp_best_score[orient2] = t_score;
	   tig_node2->best_score[orient2] = 
	     t_score + WEIGHT_SIZE * (contig2 ? contig2->num_entries : 0);
	}
	else {
	  tig_node2->best_jump[orient2] = edge1;
	}
	set_in_use(tig_node2, orient2, 1, 0);
      }
    }
  }

  for (s = tarjan_stack_index; s <= stack_end; s++) {
    tig_node2 = tarjan_stack[s];
    orient2 = tarjan_stack_orient[s];
    contig2 = tig_node2->contig;
    fprintf(stderr, " %d %c ", 
	    contig2 ? contig2->num_entries : 0, 
	    orient2 ? 'C' : ' ');
    fprintf(fp_log, " %d %c ", 
	    contig2 ? contig2->num_entries : 0, 
	    orient2 ? 'C' : ' ');
    if (stack_end > tarjan_stack_index) {
      stack_index = 0;
/*
   edge = tig_node2->best_jump[orient2];
    fprintf(stderr, " orig_best_score %d, best_jump: %d,", best_score, 
	    edge ? edge->contig2->tig_node->component[orient2 != edge->reverse] : 0);
*/
      n_visited_nodes = 0;
      if (tig_node2->n_allowed_visits) {

	score = tig_node2->temp_best_score[orient2];
	score += WEIGHT_SIZE * (contig2 ? contig2->num_entries : 0);
	set_in_use(tig_node2, orient2, -1, 0); /* to keep contig and its complement from appearing
					  in same path; note that this procedure does not
					  guarantee a maximum-scoring path */
	visit_tig_nodes(tig_node2, !orient2, component);
	set_in_use(tig_node2, orient2, 1, 0);
      }

      fprintf(stderr, " n_visited_nodes %d,  ", n_visited_nodes);
      fprintf(fp_log, " n_visited_nodes %d,  ", n_visited_nodes);

    }
  }
  for (s = tarjan_stack_index; s <= stack_end; s++) {
    tig_node2 = tarjan_stack[s];
    orient2 = tarjan_stack_orient[s];
    contig2 = tig_node2->contig;
    if (component_size > 1) fprintf(stderr, "\n");
    if (component_size > 1) fprintf(fp_log, "\n");
    fprintf(stderr, " %d %c ", 
	    contig2 ? contig2->num_entries : 0, 
	    orient2 ? 'C' : ' ');
    fprintf(fp_log, " %d %c ", 
	    contig2 ? contig2->num_entries : 0, 
	    orient2 ? 'C' : ' ');

    edge = tig_node2->best_jump[orient2];
    fprintf(stderr, " %d edges, score %d,  best_jump %d, ", 
	    tig_node2->best_stack_length[orient2], 
 	    tig_node2->best_score[orient2],
	    edge ? edge->contig2->tig_node->component[orient2 != edge->reverse] : 0);
    fprintf(fp_log, " %d edges, score %d,  best_jump %d, ", 
	    tig_node2->best_stack_length[orient2], 
 	    tig_node2->best_score[orient2],
	    edge ? edge->contig2->tig_node->component[orient2 != edge->reverse] : 0);
    if (tig_node2->best_stack_length[orient2]) {
      num_entries = tig_node2->contig->num_entries;
      for (i = 0; i < tig_node2->best_stack_length[orient2]; i++) {
	edge = tig_node2->best_stack[orient2][i];
	fprintf(stderr, "%d  ", edge->contig1->num_entries);
	fprintf(fp_log, "%d  ", edge->contig1->num_entries);
	num_entries += edge->contig1->num_entries;
	if (edge->pair)
	  fprintf(stderr, "[%.1f], ", edge->pair->LLR_score / 10.0);
	  fprintf(fp_log, "[%.1f], ", edge->pair->LLR_score / 10.0);
      }
      fprintf(stderr, "  n_reads %d ", num_entries);
      fprintf(fp_log, "  n_reads %d ", num_entries);
    }
  }
} 

visit_tig_nodes(tig_node, direction, component)
     Tig_node *tig_node;
     int direction, component;
{
  Edge *edge;
  int i, orig_direction, n_edges, t_score, orient2, same_component;
  Tig_node* tig_node2;

  if (!tig_node) fatalError("Undefined tig_node.");
  tig_node->n_allowed_visits -= 1;
  orig_direction = direction;

/*
  for (edge = tig_node->edges[direction], n_edges = 0; edge; edge = edge->next) 
    if (edge->contig2->tig_node->n_allowed_visits > 0) n_edges++;

  if (!n_edges) {
    n_paths++;
    if (! ((int)n_paths % 1000000)) {
      notify(".");
    }
  }
*/  

  for (edge = tig_node->edges[direction]; edge; edge = edge->next) {
    tig_node2 = edge->contig2->tig_node;
    if (tig_node2->n_allowed_visits > 0) {
      direction = orig_direction != edge->reverse;
 /*    t_score = edge->pair ? edge->pair->LLR_score : 0; */
      
      if (tig_node2->component[!direction] != component) /* need ! because going in
							    reverse direction now */
	continue;
      edge_stack[stack_index++] = edge;
      t_score = WEIGHT_SIZE * edge->contig2->num_entries 
	+ WEIGHT_LLR * (edge->pair ? edge->pair->LLR_score : 0);

      score += t_score;
      if (tig_node2->best_score[!direction] < score) {
	for (i = 0; i < stack_index; i++) {
	  tig_node2->best_stack[!direction][i] = edge_stack[stack_index - 1 - i];
	}
	tig_node2->best_stack_length[!direction] = stack_index;
	tig_node2->best_score[!direction] = score;
      }
      n_visited_nodes++;

      if (n_visited_nodes < 1000000)
	visit_tig_nodes(tig_node2, direction, component);

      score -= t_score;
      stack_index--;
    }
  }
  tig_node->n_allowed_visits += 1;
}

set_in_use(tig_node, orientation, value, pair_flag)
     Tig_node *tig_node;
     int orientation, value, pair_flag;
{
  Edge *edge;
  Edge **best_stack;
  int i_stack, stacksize;

  for (edge = tig_node->best_jump[orientation]; edge; 
       edge = tig_node->best_jump[orientation]) {
    tig_node = edge->contig2->tig_node;
    orientation = orientation != edge->reverse;
    tig_node->n_allowed_visits += value;
    best_stack = tig_node->best_stack[orientation];
    stacksize = tig_node->best_stack_length[orientation];
    if (pair_flag && value) {
      fprintf(stderr, "; [%d]  %d/%d", 
	      tig_node->component[orientation], tig_node->contig->num_entries, 
	      tig_node->best_score[orientation]);
      fprintf(fp_log, "; [%d]  %d/%d", 
	      tig_node->component[orientation], tig_node->contig->num_entries, 
	      tig_node->best_score[orientation]);
      if (edge->pair) {
	set_used_flag(edge->pair, 1);
	set_used_flag(edge->pair->reversed_pair, 1);
      }
    }
    for (i_stack = 0; i_stack < stacksize; i_stack++) {
      edge = best_stack[i_stack];
      tig_node = edge->contig1->tig_node;
      orientation = orientation != edge->reverse;
      tig_node->n_allowed_visits += value;
      if (pair_flag && value && edge->pair) {
	set_used_flag(edge->pair, 1);
	set_used_flag(edge->pair->reversed_pair, 1);
	fprintf(stderr, ", %d/%d", tig_node->contig->num_entries, 
	      tig_node->best_score[orientation]);
	fprintf(fp_log, ", %d/%d", tig_node->contig->num_entries, 
	      tig_node->best_score[orientation]);
      }
    }
  }
}
