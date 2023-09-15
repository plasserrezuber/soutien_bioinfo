/*****************************************************************************
#   Copyright (C) 1994-2000, 2006 by Phil Green.                          
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

#define PAIR_BLOCK_SIZE 10000
#define SCORE_BLOCK_SIZE 10000

extern Parameters *parameters;
extern FILE *fp_log; /* log file */
extern int t_num_entries;

#define EXTENT_FUDGE 7 /* adjustment to extent segments to allow for spurious extensions
			  to extents -- also in qual.c */
int num_pairs;
int single_domain;
double mask_frac;

set_domain_vars()
{
  single_domain = parameters->masklevel == 0 || parameters->masklevel == 101;
  mask_frac = parameters->masklevel ? 100.0 / parameters->masklevel : 1.0;
}

Query_domain *update_query_data(entry, score, complement, start, end)
     int entry, score, complement, start, end;
{
  Seq_entry *get_seq_entry(), *seq_entry;
  Query_domain *query_domain, *last_query_domain, *append_query_domain();
  int length, start1, end1, direc, n;
  double intersect;

  seq_entry = get_seq_entry(entry);
  length = seq_entry->seq_length - 1;

  if (complement) {
    n = length - end;
    end = length - start;
    start = n;
  }

  if (!seq_entry->query_domains) 
    return seq_entry->query_domains = append_query_domain(score, start, end); /* check for whether single domain */

  query_domain = seq_entry->query_domains;

  if (!single_domain) {
    length = end - start + 1;

    for (; query_domain; last_query_domain = query_domain, query_domain = query_domain->child[direc]) {
      if (start >= query_domain->start && end <= query_domain->end
	  || start <= query_domain->start && end >= query_domain->end) break; 
      direc = end > query_domain->end; /* can extend to one side only */
      start1 = start > query_domain->start ? start : query_domain->start;
      end1 = direc ? query_domain->end : end;
      intersect = (end1 - start1 + 1) * mask_frac; /* is neg or 0, if no intersection */
      if (intersect >= length || intersect >= query_domain->end - query_domain->start + 1) break;
      /* go to left or right */
    }
    if (!query_domain) 
      return last_query_domain->child[direc] = append_query_domain(score, start, end);
  }

  if (query_domain->best_score < score) { 
    query_domain->next_best = query_domain->best_score;
    query_domain->best_score = score; 
    query_domain->n_best = 1; 
    query_domain->start = start; /* N.B. THIS MAY DISRUPT TREE VALIDITY -- BUT WILL CLEAN UP LATER */
    query_domain->end = end;
  }
  else if (query_domain->best_score == score) { 
    query_domain->n_best += 1; 
  }
  else if (query_domain->next_best < score) { 
    query_domain->next_best = score;
  }

  /*
  if (!single_domain && score >= parameters->minscore) { /* to avoid low-scoring hits affecting domain defs -- altho they
					  can still do it if they are the initial hit */ /*
    if (query_domain->start > start) query_domain->start = start;
    if (query_domain->end < end) query_domain->end = end;
  }            
 */

  return query_domain;
}

static int t_n_pairs, t_edges12, t_n_score_hists, t_n_query_domains, t_n_query_datas;

Query_domain *append_query_domain(score, start, end)
     int score, start, end;
{
  char *our_alloc();
  Query_domain *query_domain;
  static Query_domain *head_query_domain;
  static int n_query_domains;

  if (!n_query_domains) {
    head_query_domain = (Query_domain *)our_alloc(SCORE_BLOCK_SIZE * sizeof(Query_domain));
    t_n_query_domains += SCORE_BLOCK_SIZE;
  }
  query_domain = head_query_domain + n_query_domains;
  n_query_domains = (n_query_domains + 1) % SCORE_BLOCK_SIZE;

  query_domain->start = start;
  query_domain->end = end;
  query_domain->best_score = score;
  query_domain->n_best = 1;
  query_domain->next_best = 0;
  query_domain->score_hist[0] = query_domain->score_hist[1] = 0;
  query_domain->child[0] = query_domain->child[1] = 0;
  query_domain->parent = query_domain;
  query_domain->best_pair = query_domain->recent_pair = 0;
  return query_domain;
}

Score_hist *append_score_hist(query_domain, score, type, count)
     Query_domain *query_domain;
     int score, type, count;
{
  char *our_alloc();
  static Score_hist *head_score;
  static int n_score_hists;
  Score_hist *score_hist, *s_h, *last_s_h;

  /* allocate score_hists in blocks of SCORE_BLOCK_SIZE; implies cannot free individual ones */
  /* CHANGE ORDER -- SO HIGHEST SCORE IS FIRST?? COULD MAKE MORE EFFICIENT */
  for (s_h = query_domain->score_hist[type], last_s_h = 0; s_h && s_h->score < score; last_s_h = s_h, s_h = s_h->next);
  if (s_h && s_h->score == score) {
      s_h->count += count;
      return s_h;
  }
   /* need new one to point to s_h, be pointed to by last_s_h */

  if (!n_score_hists) {
    head_score = (Score_hist *)our_alloc(SCORE_BLOCK_SIZE * sizeof(Score_hist));
    t_n_score_hists += SCORE_BLOCK_SIZE;
  }
  score_hist = head_score + n_score_hists;
  n_score_hists = (n_score_hists + 1) % SCORE_BLOCK_SIZE;

  score_hist->score = score;
  score_hist->count = count;
  score_hist->next = s_h;
  score_hist->best_pair = 0;
  if (last_s_h) {
    last_s_h->next = score_hist;
  }
  else {
    query_domain->score_hist[type] = score_hist;
  }
  return score_hist;
}

Aligned_pair *append_pair(entry1, entry2, reverse_flag)
     int entry1, entry2, reverse_flag;
{
  char *our_alloc();
  Aligned_pair *pair;
  static Aligned_pair *head_pair;
  static int n_pairs;
  Seq_entry *seq_entry;
  Seq_entry *get_seq_entry();
  int n;
  static Query_data *head_query_data;
  static int n_query_datas;

  /* allocate pairs in blocks of PAIR_BLOCK_SIZE; implies cannot free individual ones */
  if (!n_pairs) {
    head_pair = (Aligned_pair *)our_alloc(PAIR_BLOCK_SIZE * sizeof(Aligned_pair));
    t_n_pairs += PAIR_BLOCK_SIZE;
  }
  pair = head_pair + n_pairs;
  n_pairs = (n_pairs + 1) % PAIR_BLOCK_SIZE;

  pair->flags = pair->reject_flags = 0;
  pair->entry1 = entry1;
  pair->entry2 = entry2;
  pair->score = pair->LLR_score = 0;
  seq_entry = get_seq_entry(entry1);
  pair->next = seq_entry->aligned_pairs;
  seq_entry->aligned_pairs = pair;
  pair->reversed_pair = 0;
  pair->diffs = 0;
  pair->start1 = pair->start2 = pair->end1 = pair->end2 = 0;
  pair->spl5 = pair->spl3 = 0;
  if (parameters->keep_query_data) {
    if (!n_query_datas) {
      head_query_data = (Query_data *)our_alloc(SCORE_BLOCK_SIZE * sizeof(Query_data));
      t_n_query_datas += SCORE_BLOCK_SIZE;
    }
    pair->query_data = head_query_data + n_query_datas;
    n_query_datas = (n_query_datas + 1) % SCORE_BLOCK_SIZE;
    if (parameters->splice_edge_length && !reverse_flag) {
      t_edges12 += parameters->splice_edge_alloc;
      pair->query_data->edges12 = (char *)our_alloc(parameters->splice_edge_alloc * sizeof(char));
    }
    else
      pair->query_data->edges12 = 0;
    pair->query_data->splice_sites = 0;
    pair->query_data->query_domain = 0;
  }
  else pair->query_data = 0;
  return pair;
}

make_reversed_pairs()
{
  Aligned_pair *pair;
  int entry1;
  Aligned_pair *get_aligned_pairs();

  notify("Making reversed pairs ...");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) 
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) 
      if (!pair->reversed_pair) {
	make_reversed_pair(pair);
/*
  unsigned char *diff;
  int gap1, gap2;
	for (diff = pair->diffs, gap1 = gap2 = 0; *diff; diff++) {
	  gap1 += diff_gap1(*diff);
	  gap2 += diff_gap2(*diff);
	}
	if (gap1 != pair->end1 - pair->start1 + 2)
	  fprintf(stderr, "\nDifferr1: %d %d", gap1, pair->end1 - pair->start1 + 2);
	if (gap2 != pair->end2 - pair->start2 + 2)
	  fprintf(stderr, "\nDifferr2:  %d %d", gap2, pair->end2 - pair->start2 + 2);
	for (diff = pair->reversed_pair->diffs, gap1 = gap2 = 0; *diff; diff++) {
	  gap1 += diff_gap1(*diff);
	  gap2 += diff_gap2(*diff);
	}
	if (gap1 != pair->reversed_pair->end1 - pair->reversed_pair->start1 + 2)
	  fprintf(stderr, "\nDifferr3:  %d %d", gap1, pair->end1 - pair->start1 + 2);
	if (gap2 != pair->reversed_pair->end2 - pair->reversed_pair->start2 + 2)
	  fprintf(stderr, "\nDifferr4:  %d %d", gap2, pair->end2 - pair->start2 + 2);
*/ 

      }
/* note that, although pairs are being appended to lists as they are being processed,
   they are always added at the beginning (so no infinite loops are possible) */
  notify(" Done\n");
}

make_reversed_pair(pair)
     Aligned_pair *pair;
{
  Aligned_pair *append_pair();
  Aligned_pair *pair2;
  unsigned char *make_reversed_diffs();
  int length1, length2, k, m, s_e_4, s_e_8;
  char *our_alloc();

  pair2 = append_pair(pair->entry2, pair->entry1, 1);
  pair->reversed_pair = pair2;
  pair2->reversed_pair = pair;

  pair2->flags = pair->flags; /* only valid for some flags -- but others are set later */
  pair2->score = pair->score;
  pair2->diffs = make_reversed_diffs(pair->diffs, is_reverse(pair));

  pair2->spl5 = pair->spl5;
  pair2->spl3 = pair->spl3;

  /*
  if (pair->query_data->edges12) {
    s_e_4 = 4 * parameters->splice_edge_length;
    s_e_8 = 2 * s_e_4;
    pair2->query_data->edges12 = (char *)our_alloc((s_e_8 + 1) * sizeof(char));
    pair2->query_data->edges12[s_e_8] = 0;
  }
  */

  if (is_reverse(pair)) {
    length1 = get_seq_length(pair->entry1) - 1;
    length2 = get_seq_length(pair->entry2) - 1;
    pair2->start1 = length2 - pair->end2;
    pair2->start2 = length1 - pair->end1;
    pair2->end1 = length2 - pair->start2;
    pair2->end2 = length1 - pair->start1;
    /*
    if (pair->query_data->edges12) {
      for (k = 0; k < s_e_8; k++) {
	pair2->query_data->edges12[k] = comp_mat[pair->query_data->edges12[s_e_8 - 1 - k]];
      }
    }
    */
  }
  else {
    pair2->start1 = pair->start2;
    pair2->start2 = pair->start1;
    pair2->end1 = pair->end2;
    pair2->end2 = pair->end1;
    /*
    if (pair->query_data->edges12) {
      for (k = 0, m = s_e_4; k < s_e_4; k++, m++) {
	pair2->query_data->edges12[k] = pair->query_data->edges12[m];
	pair2->query_data->edges12[m] = pair->query_data->edges12[k];
      }
    }
    */
  }
}

set_rejectable_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->reject_flags |= 1;
  else pair->reject_flags &= 255 - 1;
}     

set_reverse_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->flags |= 1;
  else pair->flags &= 255 - 1;
}     

set_used_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->flags |= 2;
  else pair->flags &= 255 - 2;
}     

set_best_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->flags |= 4;
  else pair->flags &= 255 - 4;
}     

set_triple_reject_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->flags |= 8;
  else pair->flags &= 255 - 8;
}     

set_repeat_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->flags |= 16;
  else pair->flags &= 255 - 16;
}     

set_reject_chimeric_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->flags |= 32;
  else pair->flags &= 255 - 32;
}     

set_split_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->flags |= 64;
  else pair->flags &= 255 - 64;
}     

set_unaligned_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->flags |= 128;
  else pair->flags &= 255 - 128;
}     

set_reject_self_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->reject_flags |= 2;
  else pair->reject_flags &= 255 - 2;
}     

set_reject_vector_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->reject_flags |= 4;
  else pair->reject_flags &= 255 - 4;
}     

set_reject_qual_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->reject_flags |= 8;
  else pair->reject_flags &= 255 - 8;
}     

set_reject_total_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->reject_flags |= 16;
  else pair->reject_flags &= 255 - 16;
}     

set_left_trunc_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->reject_flags |= 32;
  else pair->reject_flags &= 255 - 32;
}     

set_right_trunc_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->reject_flags |= 64;
  else pair->reject_flags &= 255 - 64;
}     

set_reject_node_flag(pair, value)
     Aligned_pair *pair;
     int value;
{
  if (value) pair->reject_flags |= 128;
  else pair->reject_flags &= 255 - 128;
}     

is_rejectable(pair)
  Aligned_pair *pair;
{
  return pair->reject_flags & 1;
}     

is_reverse(pair)
  Aligned_pair *pair;
{
  return pair->flags & 1;
}     

is_used(pair)
  Aligned_pair *pair;
{
  return pair->flags & 2;
}     

is_best(pair)
  Aligned_pair *pair;
{
  return pair->flags & 4;
}     

is_triple_reject(pair)
  Aligned_pair *pair;
{
  return pair->flags & 8;
}     

is_repeat(pair)
  Aligned_pair *pair;
{
  return pair->flags & 16;
}     

is_reject_chimeric(pair)
  Aligned_pair *pair;
{
  return pair->flags & 32;
}     

is_split(pair)
  Aligned_pair *pair;
{
  return pair->flags & 64;
}     

is_unaligned(pair)
  Aligned_pair *pair;
{
  return pair->flags & 128;
}     

is_reject_self(pair)
  Aligned_pair *pair;
{
  return pair->reject_flags & 2;
}     

is_reject_vector(pair)
  Aligned_pair *pair;
{
  return pair->reject_flags & 4;
}     

is_reject_qual(pair)
  Aligned_pair *pair;
{
  return pair->reject_flags & 8;
}     

is_reject_total(pair)
  Aligned_pair *pair;
{
  return pair->reject_flags & 16;
}     

is_left_trunc(pair)
  Aligned_pair *pair;
{
  return pair->reject_flags & 32;
}     

is_right_trunc(pair)
  Aligned_pair *pair;
{
  return pair->reject_flags & 64;
}     

is_reject_node(pair)
  Aligned_pair *pair;
{
  return pair->reject_flags & 128;
}     

is_reject_any(pair)
  Aligned_pair *pair;
{
  return is_reject_total(pair) || is_reject_vector(pair) || is_reject_self(pair);
}     

extern Aligned_pair *sig_pair;
extern Database *query_db; 

int length1, minscore, near_flag, splice_flag;
Cand_pair *cand_pairs;
Seq_entry *seq_entry;
static Profile *q_profile;
static q_profile_flag;

find_scores(entry1, seq1, rev_store_flag)
     int entry1, rev_store_flag;
     unsigned char *seq1;
{
  int length2, prev_loc, loc, next_loc, found;
  int start, end, prev_offset, offset, new_offset, entry;
  Profile *make_profile_from_seq();
  Segment *segment;
  unsigned char *get_seq(), *get_comp_seq();
  unsigned char *seq2;
  Cand_pair *get_cand_pairs();
  Aligned_pair *pair, *pair2, *last_pair;
  Seq_entry *get_seq_entry();
  Align_info *get_align_entry();

  /* FIX THIS -- MAKING 
  if (is_duplicate(get_align_entry(entry1))) return 0;
  */

  minscore = parameters->minscore;
  near_flag = parameters->minscore != parameters->near_minscore;
  splice_flag = parameters->splice_edge_length;

  seq_entry = get_seq_entry(entry1);
  length1 = get_seq_length(entry1);

  if (!rev_store_flag) {
    if (cand_pairs = get_cand_pairs(entry1)) {
      q_profile = make_profile_from_seq(q_profile, seq1, length1, 0); 
      visit_cand_pairs_scores(cand_pairs);
    }
  }
  else {
    q_profile_flag = 0;
    for (entry = 0; entry < query_db->num_entries; entry++) {
      if (cand_pairs = get_cand_pairs(entry)) {
	if (!q_profile_flag) {
	  q_profile = make_profile_from_seq(q_profile, seq1, length1, 0); 
	  q_profile_flag = 1;
	}
	visit_cand_pairs_scores(cand_pairs); 
	get_seq_entry(entry)->cand_pairs = 0;
      }
    }
  }
  /*
  for (cand_pair = cand_pairs; cand_pair; cand_pair = cand_pair->next) {
    reset_match_pairs();
    seq2 = cand_pair->reverse ? get_comp_seq(cand_pair->entry2) : get_seq(cand_pair->entry2);
    length2 = get_seq_length(cand_pair->entry2);
    prev_loc =   -length2 - parameters->max_intron_length;
    for (segment = cand_pair->band_segments; segment; segment = segment->next) {
      /* if segment is isolated, apply minscore threshold (minscore to 2dary hits); if not, near_minscore */
  /*
      if (near_flag) {
	loc = (segment->start + segment->end) / 2;
	if (sig_pair) prev_loc = sig_pair->start1 - sig_pair->start2;

	minscore = parameters->minscore;
	if (segment != cand_pair->band_segments && parameters->max_intron_length > loc - prev_loc
	    || segment->next && parameters->max_intron_length > (segment->next->start + segment->next->end) / 2 - loc) {
	  minscore = parameters->near_minscore;
	}
      }
      recursive_swat(minscore, parameters->near_minscore, cand_pair, q_profile, seq2, length2, 
		     segment->start, segment->end, 1, length1, 1, length2);
    }
    if (near_flag) {
      /* these are now in reverse order!! */
  /*
      prev_loc = length1 + parameters->max_intron_length;
      last_pair = 0;
      for (pair = seq_entry->aligned_pairs; 
	   pair && pair->entry2 == cand_pair->entry2 && is_reverse(pair) == cand_pair->reverse; 
	   pair = pair->next) {
	if (pair->score >= parameters->minscore) {
	  prev_loc = pair->start1 - pair->start2;
	  last_pair = pair;
	  continue;
	}
	if (prev_loc - (pair->start1 - pair->start2) <= parameters->max_intron_length) {
	  last_pair = pair;
	  continue;
	}
	for (pair2 = pair->next, found = 0; 
	     pair2 && pair2->entry2 == cand_pair->entry2 && is_reverse(pair2) == cand_pair->reverse; 
	     pair2 = pair2->next) {
  /* could also stop when distance exceeds max_intron_length -- but probably takes more time than worth */
  /*
	  if (pair2->score >= parameters->minscore) {
	    found = 1;
	    break;
	  }
	}
	if (found && (pair->start1 - pair->start2) - (pair2->start1 - pair2->start2) <= parameters->max_intron_length) {
	  last_pair = pair;
	  continue;
	}
	if (!last_pair) {
	  seq_entry->aligned_pairs = pair->next;
	}
	else
	  last_pair->next = pair->next;
	our_free(pair->query_data->edges12); /* would be better to free pair!! */ /*
	/* pair->score = 0; */ /*
				 num_pairs--;
				 }
				 for (pair = seq_entry->aligned_pairs; 
				 pair && pair->entry2 == cand_pair->entry2 && is_reverse(pair) == cand_pair->reverse; 
				 pair = pair->next) {
				 find_left_sites(pair, cand_pair, q_profile->seq, seq2);
				 find_right_sites(pair, cand_pair, q_profile->seq, seq2);
				 }
				 for (pair = seq_entry->aligned_pairs; 
				 pair && pair->entry2 == cand_pair->entry2 && is_reverse(pair) == cand_pair->reverse; 
				 pair = pair->next) {
				 prev_offset = offset = pair->start1 - pair->start2;
				 for (pair2 = pair->next; pair2 && pair2->entry2 == cand_pair->entry2 && is_reverse(pair2) == cand_pair->reverse; pair2 = pair2->next) {
				 new_offset = pair2->start1 - pair2->start2;
				 if (offset - new_offset > parameters->max_intron_length) break;
				 if (prev_offset == offset && offset - new_offset > parameters->min_intron_length)
				 prev_offset = new_offset;
				 else if (prev_offset != offset && prev_offset - new_offset > parameters->min_intron_length)
				 break;
				 if (pair->start2 > pair2->end2 + 50) continue; /* large gaps likely mean erroneous data, or intervening pair */ /*
																		   /*  AVOID DOUBLE COUNTING!! */ /*
																						    analyze_splices(pair2, pair, 0, 0, 0, q_profile->seq, seq2); 
																						    }
																						    }
																						    }
																						    }
																						  */
}

visit_cand_pairs_scores(cand_pair)
     Cand_pair *cand_pair;
{
  int length2, prev_loc, loc, next_loc, found;
  int start, end, prev_offset, offset, new_offset, entry2, reverse;
  Segment *segment;
  unsigned char *get_seq(), *get_comp_seq();
  unsigned char *seq2;
  Aligned_pair *pair, *pair2, *last_pair;
  Align_info *get_align_entry();

  if (!cand_pair) return;
  visit_cand_pairs_scores(cand_pair->left);
  entry2 = cand_pair->entry2;
  if (!is_duplicate(get_align_entry(entry2))) {
    reset_match_pairs();
    reverse = cand_pair->reverse;
    seq2 = reverse ? get_comp_seq(entry2) : get_seq(entry2);
    length2 = get_seq_length(entry2);
    prev_loc =  -length2 - parameters->max_intron_length;
    for (segment = cand_pair->band_segments; segment; segment = segment->next) {
      /* if segment is isolated, apply minscore threshold (minscore to 2dary hits); if not, near_minscore */
      if (near_flag) {
	loc = (segment->start + segment->end) / 2;
	if (sig_pair) prev_loc = sig_pair->start1 - sig_pair->start2;

	minscore = parameters->minscore;
	if (segment != cand_pair->band_segments && parameters->max_intron_length > loc - prev_loc
	    || segment->next && parameters->max_intron_length > (segment->next->start + segment->next->end) / 2 - loc) {
	  minscore = parameters->near_minscore;
	}
      }
      recursive_swat(minscore, parameters->near_minscore, cand_pair, q_profile, seq2, length2, 
		     segment->start, segment->end, 1, length1, 1, length2);
    }
    if (near_flag) {
      /* these are now in reverse order!! */
      prev_loc = length1 + parameters->max_intron_length;
      last_pair = 0;
      for (pair = seq_entry->aligned_pairs; 
	   pair && pair->entry2 == entry2 && is_reverse(pair) == reverse; 
	   pair = pair->next) {
	if (pair->score >= parameters->minscore) {
	  prev_loc = pair->start1 - pair->start2;
	  last_pair = pair;
	  continue;
	}
	if (prev_loc - (pair->start1 - pair->start2) <= parameters->max_intron_length) {
	  last_pair = pair;
	  continue;
	}
	for (pair2 = pair->next, found = 0; 
	     pair2 && pair2->entry2 == entry2 && is_reverse(pair2) == reverse; 
	     pair2 = pair2->next) {
	  /* could also stop when distance exceeds max_intron_length -- but probably takes more time than worth */
	  if (pair2->score >= parameters->minscore) {
	    found = 1;
	    break;
	  }
	}
	if (found && (pair->start1 - pair->start2) - (pair2->start1 - pair2->start2) <= parameters->max_intron_length) {
	  last_pair = pair;
	  continue;
	}
	if (!last_pair) {
	  seq_entry->aligned_pairs = pair->next;
	}
	else
	  last_pair->next = pair->next;
	if (parameters->splice_edge_length) {
	  our_free(pair->query_data->edges12); /* would be better to free pair!! */
	  t_edges12 -= parameters->splice_edge_alloc;
	}
	/* pair->score = 0; */
	num_pairs--;
	/*	notify("decrement"); */
      }
    }
    if (splice_flag) {
      for (pair = seq_entry->aligned_pairs; 
	   pair && pair->entry2 == entry2 && is_reverse(pair) == reverse; 
	   pair = pair->next) {
	find_left_sites(pair, cand_pair, q_profile->seq, seq2);
	find_right_sites(pair, cand_pair, q_profile->seq, seq2);
      }
      for (pair = seq_entry->aligned_pairs; 
	   pair && pair->entry2 == entry2 && is_reverse(pair) == reverse; 
	   pair = pair->next) {
	prev_offset = offset = pair->start1 - pair->start2;
	for (pair2 = pair->next; pair2 && pair2->entry2 == entry2 && is_reverse(pair2) == reverse; pair2 = pair2->next) {
	  new_offset = pair2->start1 - pair2->start2;
	  if (offset - new_offset > parameters->max_intron_length) break;
	  if (prev_offset == offset && offset - new_offset > parameters->min_intron_length)
	    prev_offset = new_offset;
	  else if (prev_offset != offset && prev_offset - new_offset > parameters->min_intron_length)
	    break;
	  if (pair->start2 > pair2->end2 + 50) continue; /* large gaps likely mean erroneous data, or intervening pair */
	  /*  AVOID DOUBLE COUNTING!! */
	  analyze_splices(pair2, pair, 0, 0, 0, q_profile->seq, seq2); 
	}
      }
    }
  }
  visit_cand_pairs_scores(cand_pair->right);
}

find_all_scores()
{
  int entry1;
  unsigned char *get_seq();
  
  notify("SWATTING ...");
  
  for (entry1 = 0; entry1 < t_num_entries; entry1++) 
    find_scores(entry1, get_seq(entry1), 0);

  notify(" Done\n");
}

/* CHANGE -- CALL smith_waterman directly */
make_full_pairs(db)
     Database *db;
{
  int entry1, entry2;
  File *file;
  Aligned_pair *get_aligned_pairs();
  unsigned char *seq;
  unsigned char *get_seq();
  Database *sdb;

  notify("Making full pairs ...");

  if (parameters->subject_files) {
    for (file = parameters->subject_files; file; file = file->next) {
      sdb = file->db;
      while (get_next_file_entry(sdb)) {
	entry1 = append_seq_entry(sdb);
	for (entry2 = db->first_entry; entry2 <= db->last_entry; entry2++) {
	  make_new_cand_pair(entry1, entry2, -2, -2, 0, 0);
	  make_new_cand_pair(entry1, entry2, -2, -2, 1, 0);
	}
	find_scores(entry1, sdb->seq_buffer, 0);
	if (!get_aligned_pairs(entry1)) remove_seq_entry(sdb);
	free_cand_pair_blocks();
	free_seg_blocks();
      }
    }
  }
  else {
    for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
      for (entry2 = db->first_entry; entry2 <= entry1; entry2++) {
	make_new_cand_pair(entry1, entry2, -2, -2, 0, 0);
	if (entry2 < entry1) make_new_cand_pair(entry1, entry2, -2, -2, 1, 0);
      }
      find_scores(entry1, get_seq(entry1), 0);
      free_cand_pair_blocks();
      free_seg_blocks();
    }
  }

  notify(" Done\n");
}

print_t_n_pairs()
{
  fprintf(stderr, "Total # pairs: %d, size: %.3f Mbytes; edges12: %d; # score_hists: %d, size: %.3f Mbytes; # query_domains: %d, size: %.3f Mbytes; # query_datas: %d, size: %.3f Mbytes\n", 
	  t_n_pairs, t_n_pairs * (sizeof(Aligned_pair) / 1000000.), t_edges12, 
	  t_n_score_hists, t_n_score_hists * (sizeof(Score_hist) / 1000000.),
	  t_n_query_domains, t_n_query_domains * (sizeof(Query_domain) / 1000000.),
	  t_n_query_datas, t_n_query_datas * (sizeof(Query_data) / 1000000.));
}

/* 
reswat()
{
  int entry1, n_scores;
  Aligned_pair *pair,*pair2;
  int full_smith_waterman();
  int start1, end1, start2, end2, length2;
  Profile *make_profile_from_seq();
  Profile *q_profile;
  int mismatches, insertions, deletions, n_diffs;
  int temp, off1, off2;
  Align_info *align_entry1, *align_entry2;
  unsigned char *make_diffs(), *make_reversed_diffs();
  int orig_score;
  int success, n_successes, n_failures;
  int temp_int;
  unsigned char *temp_diffs;
  Db_entry *db_entry2;

  notify("ReSWATTING ...");

  if (parameters->weak_first) {
FORMAT NOW CHANGED 
    set_score_mat(parameters->matrix, parameters->gap_init, parameters->gap_ext, 0, parameters->penalty);
    reset_n_successes();
  }
  else {
  set_score_mat(0, parameters->weak_penalty - 2, parameters->weak_penalty - 1, 
		0, parameters->weak_penalty);
    n_successes = n_failures = 0;
  }
  n_scores = 0;
  for (entry1 = first_file_only ? num_query_entries : 0;
       entry1 < t_num_entries; entry1++) {
    align_entry1 = align_array + entry1;
    if (!align_entry1->first_pair) continue;
    q_profile = make_profile_from_seq((Profile *)0, align_entry1->db_entry->seq, align_entry1->db_entry->length, 0);
    for (pair = align_entry1->first_pair; pair; pair = pair->next) {
      if (first_file_only || pair->entry2 >= pair->entry1) {
	align_entry2 = align_array + pair->entry2;
	length2 = align_entry2->db_entry->length - 1;
	db_entry2 = align_entry2->db_entry + (is_reverse(pair) ? num_query_entries : 0);
	off1 = pair->start1 - pair->start2;
	off2 = pair->end1 - pair->end2;
	if (off1 > off2) {
	  temp = off1;
	  off1 = off2;
	  off2 = temp;
	}
	if (parameters->weak_first) {
	  recursive_swat(parameters->minscore, nextscore, pair, q_profile, db_entry2, 
			 pair->band_segments->start, pair->band_segments->end,
			 pair->start1 + 1, pair->end1 + 1, pair->start2 + 1, pair->end2 + 1);
	}
	else {
	  pair->w_score = quick_full_smith_waterman(q_profile, db_entry2->seq, db_entry2->length,
						    off1 - parameters->bandwidth - 10, 
						    off2 + parameters->bandwidth + 10, 0,0,0,0, 
						    parameters->minscore, 
						    &orig_score, &success);

	  success ? n_successes++ : n_failures++;
	  get_stats(&start1, &end1, &start2, &end2,
		    &mismatches, &insertions, &deletions);
	  n_scores++;
	  if (pair->w_score >= pair->score && start1 - 1 <= pair->start1
	      && end1 - 1 >= pair->end1) {
	      pair->w_diffs = make_diffs();
	      pair->w_start1 = start1 - 1;
	      pair->w_start2 = start2 - 1;
	      pair->w_end1 = end1 - 1;
	      pair->w_end2 = end2 - 1;
	  }
	  else {
	    fprintf(fp_log, "\nWARNING: w_score, w_start, or w_end inconsistency");
	    fprintf(fp_log, "\n%s  %s  new, old starts: %d  %d; ends:  %d  %d, scores:  %d  %d, reverse %d  ",
		    align_array[pair->entry1].db_entry->id, 
		    align_array[pair->entry2].db_entry->id, 
		    start1 - 1, pair->start1, end1 - 1, pair->end1,
		    pair->w_score, pair->score, is_reverse(pair));
	    pair->w_diffs = pair->diffs; 
	    pair->w_start1 = pair->start1;
	    pair->w_start2 = pair->start2;
	    pair->w_end1 = pair->end1;
	    pair->w_end2 = pair->end2;
	  }
	}
      }
    }
    free_profile(q_profile);
  }
  if (parameters->weak_first) {
    for (entry1 = first_file_only ? num_query_entries : 0;
	 entry1 < t_num_entries; entry1++) {
      align_entry1 = align_array + entry1;
      for (pair = align_entry1->first_pair; pair; pair = pair->next) {
	temp_int = pair->score;
	pair->score = pair->w_score;
	pair->w_score = temp_int;

	temp_int = pair->start1;
	pair->start1 = pair->w_start1;
	pair->w_start1 = temp_int;

	temp_int = pair->end1;
	pair->end1 = pair->w_end1;
	pair->w_end1 = temp_int;

	temp_int = pair->start2;
	pair->start2 = pair->w_start2;
	pair->w_start2 = temp_int;

	temp_int = pair->end2;
	pair->end2 = pair->w_end2;
	pair->w_end2 = temp_int;

	temp_diffs = pair->diffs;
	pair->diffs = pair->w_diffs;
	pair->w_diffs = temp_diffs;
      }
    }
    notify(" Done\n");
    get_n_successes();
  }
  else {
    fprintf(fp_log, "\n%d alignments in reswat.", n_scores);
    notify(" Done\n");
    fprintf(stderr, "Quickalign: %d successes, %d failures\n", n_successes, n_failures);
  }
}
*/

compare_entry1_and_scores(pair_1, pair_2)
     Aligned_pair **pair_1, **pair_2;
{
  Aligned_pair *pair1, *pair2;
  int d;
  
  pair1 = *pair_1;
  pair2 = *pair_2;
  if (d = pair1->entry1 - pair2->entry1) return d;
  /*  if (d = pair1->entry2 - pair2->entry2) return d; */
  return (pair2->score - pair1->score);
}

compare_entry1_and_starts(pair_1, pair_2)
     Aligned_pair **pair_1, **pair_2;
{
  Aligned_pair *pair1, *pair2;
  int d;
  
  pair1 = *pair_1;
  pair2 = *pair_2;
  if (d = pair1->entry1 - pair2->entry1) return d;
  /*  if (d = pair1->entry2 - pair2->entry2) return d; */
  if ((d = !!pair1->score - !!pair2->score) || !pair1->score) return d;
  return pair1->start1 - pair2->start1;
/* (pair1->diffs->site1 - pair2->diffs->site1); */
}

compare_entry12_and_starts(pair_1, pair_2)
     Aligned_pair **pair_1, **pair_2;
{
  Aligned_pair *pair1, *pair2;
  int d;
  
  pair1 = *pair_1;
  pair2 = *pair_2;
  if (d = pair1->entry1 - pair2->entry1) return d;
  /*  if (d = pair1->entry2 - pair2->entry2) return d; */
  if ((d = !!pair1->score - !!pair2->score) || !pair1->score) return d;

  if (d = pair1->entry2 - pair2->entry2) return d;
  if (d = is_reverse(pair1) - is_reverse(pair2)) return d;

  return /* is_reverse(pair1) ? pair2->start2 - pair1->start2 : */ pair1->start2 - pair2->start2;
/* (pair1->diffs->site1 - pair2->diffs->site1); */
}

/* sets pair->used, pair->reject_flags  */
/*
typedef struct score_list {
  int score, count[2];
} Score_list;

Score_list *score_list;
*/

print_matches()
{
  int i_ptr, prev_entry, i, j, n, m, prev_strand;
  Aligned_pair *pair, *prev_pair, *best_pair;
  Aligned_pair **sort_pairs();
  Aligned_pair **pair_pointers;
  int length1, length2;
  int p_c;
  Segment *insert_segment();
  Segment *segments;
  int sort_flag, n_p, n_entries1;
  char *get_id();
  unsigned char *get_seq();
  char *get_adj_qual();
  int splice_edge_len;
  Profile *maximal_profile_segments();
  char *our_alloc();
  /*  Score_list *s_h; */
  int n_s_h, i_s_h, prev_sh_entry, best_score, t_count, type, profile_flag, margin;
  int score_counts[4][1001], start_counts[4][1001], splice_counts[4][100];
  double cum[4], total[4];
  Seq_entry *get_seq_entry(), *seq_entry;
  Score_hist *score_hist;
  Query_domain *query_domain;
  double get_subject_length();

  notify("Printing matches ... ");

  printf("\n\nTotal subject length: %.0f residues, final minscore: %d", get_subject_length(), parameters->minscore);

  printf("\n\nNum. pairs: %d", num_pairs);

  for (i = 0; i < 4; i++)
    for (j = 0; j < 100; j++)
      splice_counts[i][j] = 0;

  for (j = 0; j < 4; j++)
    for (i = 0; i < 1001; i++) 
      score_counts[j][i] = start_counts[j][i] = 0;

  /*
  if (parameters->score_hist) {
    score_list = (Score_list *)our_alloc(1000 * sizeof(Score_list));
  }
  */

  if (parameters->subject_files) {
    sort_flag = 0;
    n_p = num_pairs;
  }
  else {
    sort_flag = 1;
    n_p = 2 * num_pairs;
  }  
  pair_pointers = sort_pairs(compare_entry1_and_scores, sort_flag); /* sort in decreasing score order */
  segments = 0;
  n_entries1 = 0;
  prev_entry = -1;
  printf("\nMaximal single base matches (low complexity regions):");
  if (q_profile) free_profile(q_profile);
  q_profile = 0;

  /* printing rules: minscore (or near_minscore) threshold already imposed (no other matches saved -- ?). Add'l filters:
        -masklevel 101: none (all matches printed)
        -masklevel <= 100:
            -minmargin <= 0:  score >= query_domain->best_score + (int)parameters->minmargin
            -minmargin == 0.5: score == query_domain->best_score AND pair == query_domain->best_pair (must check reversed pair!!)
            -minmargin >= 1:                   "               " AND query_domain->n_best == 1 AND (if minmargin > 1) closest other score is <= minmargin less
  */
  for (i_ptr = 0; i_ptr < n_p; i_ptr++) {
    pair = pair_pointers[i_ptr];

    if (!pair->score) continue;

    if (pair->entry1 != prev_entry) { /* new entry */
      n_entries1++;
      prev_entry = pair->entry1;
      profile_flag = 1;
      /*
      if (parameters->score_flag) {
	query_domain = get_seq_entry(prev_entry)->query_domains;
	test_qd1_merges(query_domain, query_domain); 
      }
      */
      
      /* merge_query_domains(prev_entry); */
      /*
      seq_entry = get_seq_entry(prev_entry);
      best_score = seq_entry->score;
      best_pair = pair;
      if (parameters->score_flag && best_score != pair->score) fatalError("mismatching best scores");
      pair->reject_flags = 2; 
      if (parameters->masklevel < 0) {
	for (type = t_count = 0; type < 2; type++) {
	  for (score_hist = seq_entry->score_hist[type]; score_hist; score_hist = score_hist->next) {
	    if (score_hist->score > best_score + parameters->masklevel) t_count += score_hist->count;
	  }
	}
	if (t_count > 1) pair->reject_flags = 0;
      }

      if (pair->reject_flags) {
	segments = insert_segment((Segment *)0, pair->start1, pair->end1);
	q_profile = maximal_profile_segments(q_profile, get_seq(prev_entry), get_seq_length(prev_entry), 
					     parameters->minscore, get_id(prev_entry), 0);
      }
      */
    }
    /*
    else if (best_pair->reject_flags) {
      if (parameters->minmargin >= 0 && pair->score < best_score - parameters->minmargin) {
	pair->reject_flags = 0;
      }
      else {
	p_c = percent_contained(segments, pair->start1, pair->end1);
	pair->reject_flags = p_c >= parameters->masklevel ? 0 : (p_c > 0 ? 1 : 2);
	if (pair->reject_flags) {
	  segments = insert_segment(segments, pair->start1, pair->end1);
	}
      }
    }
    */

    if (parameters->masklevel >= 101) {
      pair->reject_flags = 2;
    }
    else {
      for (query_domain = pair->reversed_pair->query_data->query_domain;
	   query_domain != query_domain->parent; query_domain = query_domain->parent); 
      best_score = query_domain->best_score;
      pair->reject_flags = pair->score >= best_score + (parameters->minmargin >= 0 ? 0.0 : parameters->minmargin) 
	&& (parameters->minmargin < 1 || query_domain->n_best == 1)
	&& (parameters->minmargin <= 0 || pair->reversed_pair == query_domain->best_pair) ? 2 : 0; 

      if (pair->reject_flags && parameters->minmargin > 1) {
	/*
	for (type = t_count = 0; type < 2; type++) {
	  for (score_hist = query_domain->score_hist[type]; score_hist; score_hist = score_hist->next) {
	    if (score_hist->score > best_score - parameters->minmargin) t_count += score_hist->count;
	  }
	}
	if (t_count > 1) pair->reject_flags = 0;
	*/
	if (best_score - query_domain->next_best < parameters->minmargin)
	  pair->reject_flags = 0;
      }
    }    /* N.B. NO LONGER DISTINGUISHING REJECT_FLAGS == 1 FROM 2 -- SO ASTERISKS NOT PRINTED */
    if (parameters->DNA_flag < 3 && profile_flag && pair->reject_flags) {
	q_profile = maximal_profile_segments(q_profile, get_seq(prev_entry), get_seq_length(prev_entry), 
					     parameters->minscore, get_id(prev_entry), 0);
	profile_flag = 0;
    }
  }
/*   free_segments(segments); */

  our_free(pair_pointers);
  
  splice_edge_len = parameters->splice_edge_length;

  if (splice_edge_len) {
    pair_pointers = sort_pairs(compare_entry12_and_starts, sort_flag);
    /* set_splice_params(); /* some of these should eventually become settable */
  }
  else
    pair_pointers = sort_pairs(compare_entry1_and_starts, sort_flag);

  /* if looking for exons -- also sort by subject & subject strand, and position in subject?? */
  prev_entry = prev_sh_entry = -1;
  prev_pair = 0;

  if (!parameters->discrep_tables) init_diff_hist();

  for (i_ptr = 0; i_ptr < n_p; i_ptr++) {
    pair = pair_pointers[i_ptr];
    if (!pair->score) continue;
    m = !pair->reject_flags + (is_unaligned(pair) ? 2 : 0);
    score_counts[m][pair->score > 1000 ? 1000 : pair->score] += 1;
    /*    if (pair->score == 16) */
    start_counts[m][pair->start1 > 1000 ? 1000 : pair->start1] += 1;
   /* if new pair for an entry, print; -- should coordinate with prev_entry instead to eliminate this variable */
    if (parameters->score_hist) {
      if (pair->entry1 != prev_sh_entry) {
	if (prev_sh_entry >= 0) {
	  if (prev_sh_entry != prev_entry) printf("\n");
	  print_score_list(prev_sh_entry, get_seq_entry(prev_sh_entry)->query_domains);
	}
	/*	score_list->score = score_list->count[0] = score_list->count[1] = 0; 
	n_s_h = 0;
	*/
	prev_sh_entry = pair->entry1;
      }
      /*
      for (s_h = score_list; pair->score < s_h->score; s_h++);
      if (pair->score == s_h->score) s_h->count[!!is_unaligned(pair)] += 1;
      else {
	n_s_h++;
	for (i_s_h = n_s_h; i_s_h > s_h - score_list; i_s_h--) {
	  score_list[i_s_h].score = score_list[i_s_h - 1].score;
	  score_list[i_s_h].count[0] = score_list[i_s_h - 1].count[0];
	  score_list[i_s_h].count[1] = score_list[i_s_h - 1].count[1];
	}
	s_h->score = pair->score;
	s_h->count[!!is_unaligned(pair)] = 1;
      }
      */
    }

    if (is_unaligned(pair)) {
      splice_counts[!pair->reject_flags][63 & pair->spl5] += 1;
      splice_counts[2 + !pair->reject_flags][63 & pair->spl3] += 1;
    }
    if (!pair->reject_flags) continue; 
    if (pair->entry1 != prev_entry) {
      prev_entry = pair->entry1;
      if (splice_edge_len) {
	filter_matches(i_ptr, n_p, pair_pointers, 1);
	analyze_multiple(i_ptr, n_p, pair_pointers, 1);
      }
      else {
	printf("\n");
      }
/* following inserted 10/25/98 */      
      free_int_qual();
      length1 = get_seq_length(pair->entry1);
      set_int_qual(get_adj_qual(pair->entry1), length1, 0);
      if (!pair->score) continue;
      /*      if (!pair->reject_flags) continue; /* because may have been set to 0 */
    }
    if (splice_edge_len && (pair->reject_flags & 16) && !(pair->reject_flags & 4) && !(pair->reject_flags & 8)) {
      /* analyze_splices(prev_pair, pair, 1, 1, 1, (char *)0, (char *)0); */

      prev_pair = pair;
    }
    printf("\n");

    if (parameters->qual_scores) {
      revise_score(pair);
      for (query_domain = pair->reversed_pair->query_data->query_domain;
	   query_domain != query_domain->parent; query_domain = query_domain->parent); 
      margin = query_domain->n_best > 1 ? 0 : query_domain->best_score - query_domain->next_best;
      if (pair->score > margin) pair->score = margin;
      print_pair(pair, stdout, 0);
    }
    else 
      print_pair(pair, stdout, 1);

    length2 = get_seq_length(pair->entry2);
    if (parameters->alignments) {
      print_alignment_from_diffs(pair->reversed_pair);
    }
    /*
      set_int_qual(get_adj_qual(pair->entry1), length1, 0);
    */

    if (parameters->discrep_lists) print_diffs(pair, length2); 

    if (parameters->discrep_tables) {
      init_diff_hist();
      incr_diff_hist(get_seq(pair->entry1),
		     pair->diffs, length1, length2, pair->start1, 
		     pair->start2, pair->end1, pair->end2, 0, 0 /* is_unaligned(pair) */);
      print_diff_hist("Discrepancy histogram");
      free_diff_hist();
    }
    else if (parameters->qual_flag) {

      incr_diff_hist(get_seq(pair->entry1),
		     pair->diffs, length1, length2, pair->start1, 
		     pair->start2, pair->end1, pair->end2, 0, 0 /* is_unaligned(pair) */);
    }  

    /*
      free_int_qual();
    */  
/*     if (pair->score > length1 / 10 || pair->score > length2 / 10) {
          only print diffs and histogram for "long" matches -- to avoid repeats etc. 
*/
  }
  if (prev_sh_entry >= 0) {
    print_score_list(prev_sh_entry, get_seq_entry(prev_sh_entry)->query_domains);
  }
 
  /*
    for (entry1 = first_file_only ? num_query_entries : 0;
    entry1 < t_num_entries; entry1++) {
    if (!align_array[entry1].first_pair) continue;
    for (pair = align_array[entry1].first_pair; pair;
    pair = pair->next) {
    printf("\n");
    print_pair(pair, stdout);
    }
    }
    */

  if (!parameters->discrep_tables) {      
    printf("\n\n%d matching entries (first file).", n_entries1);
    print_diff_hist("Discrepancy summary:");
    free_diff_hist();
  }

  if (splice_edge_len) {
    print_splice_hists();
  }

  printf("\n\nScore histogram: displayed unspliced matches/cum/frac, other unspliced, displayed spliced, other spliced:");
  for (j = 0; j < 4; j++) {
    total[j] = cum[j] = 0;
    for (i = 1000; i >= 0; i--) 
      total[j] += score_counts[j][i];
  }
  for (i = 1000; i >= 0; i--) {
    if (score_counts[0][i] + score_counts[1][i]) {
      printf("\n%s%5d", i == 1000 ? ">" : "", i);
      for (j = 0; j < 4; j++) {
	cum[j] += score_counts[j][i];
	printf(" %4d %5.0f %.3f  ", score_counts[j][i], cum[j], cum[j] ? cum[j] / total[j] : 0.0);
      }
    }
  }

  printf("\n\nStart histogram: displayed unspliced matches/cum/frac, other unspliced, displayed spliced, other spliced:");
  for (j = 0; j < 4; j++) {
    total[j] = cum[j] = 0;
    for (i = 1000; i >= 0; i--) 
      total[j] += start_counts[j][i];
  }
  for (i = 1000; i >= 0; i--) {
    if (start_counts[0][i] + start_counts[1][i]) {
      printf("\n%s%5d", i == 1000 ? ">" : "", i);
      for (j = 0; j < 4; j++) {
	cum[j] += start_counts[j][i];
	printf(" %4d %5.0f %.3f  ", start_counts[j][i], cum[j], cum[j] ? cum[j] / total[j] : 0.0);
      }
    }
  }

  printf("\n");

  notify(" Done\n");
  our_free(pair_pointers);
  if (parameters->output_bcdsites) {
    notify("Writing BCDsites...");
    open_diffseg_file();
    write_diffsegnodes();
    notify(" Done\n");
  }

  printf("\n\nSplice histogram : displayed 5' splice, other 5' splice: displayed 3' splice, other 3' splice:");
  for (j = 0; j < 4; j++) {
    total[j] = 0;
    for (i = 0; i < 100; i++) 
      total[j] += splice_counts[j][i];
  }

  for (i = 0; i < 100; i++) {
    if (splice_counts[0][i] + splice_counts[1][i] + splice_counts[2][i] + splice_counts[3][i]) {
      printf("\n%3d %c%c%c ", i, "ACGTN"[i / 16], "ACGTN"[(i / 4) % 4], "ACGTN"[i % 4]);
      for (j = 0; j < 4; j++) {
	printf(" %4d %.3f  ", splice_counts[j][i], splice_counts[j][i] ? splice_counts[j][i] / total[j] : 0.0);
      }
    }
  }
}

int trim_flag;

merge_qds(db, t_f)
     Database *db;
     int t_f;
{
  Query_domain *query_domain;
  int entry, keep_flag, best_score;
  Seq_entry *get_seq_entry();
  Aligned_pair *get_aligned_pairs(), *pair;

  if (parameters->masklevel == 101) return;

  trim_flag = t_f;
  if (parameters->masklevel > 0) {
    for (entry = db->first_entry; entry < db->first_entry + db->num_entries; entry++) {
      query_domain = get_seq_entry(entry)->query_domains;
      if (trim_flag)
	init_qd_trims(query_domain);
      test_qd1_merges(query_domain, query_domain);
    }
  }
  
  for (entry = 0; entry < t_num_entries; entry++) { /* note also considering entries without any pairs!! */
    for (pair = get_aligned_pairs(entry); pair; pair = pair->next) {
      if (!pair->score) continue;
      for (query_domain = pair->query_data->query_domain; query_domain != query_domain->parent; 
	   query_domain = query_domain->parent); 
      best_score = query_domain->best_score;
      keep_flag = pair->score >= best_score + (parameters->minmargin >= 0 ? 0.0 : parameters->minmargin) 
	&& (parameters->minmargin < 1 || query_domain->n_best == 1)
	&& (parameters->minmargin <= 0 || pair == query_domain->best_pair);

      if (keep_flag && parameters->minmargin > 1) {
	if (best_score - query_domain->next_best < parameters->minmargin)
	  keep_flag = 0;
      }
      if (!keep_flag) pair->score = 0;
    }
  }
}

init_qd_trims(query_domain)
     Query_domain *query_domain;
{
  if (!query_domain) return;
  
  if (query_domain == query_domain->parent) {
    query_domain->trim_start = query_domain->start;
    query_domain->trim_end = query_domain->end;
  }
  init_qd_trims(query_domain->child[0]);
  init_qd_trims(query_domain->child[1]);
}

test_qd1_merges(query_domain1, query_domain2)
     Query_domain *query_domain1, *query_domain2;
{
  if (!query_domain1) return;

  test_qd2_merges(query_domain1, query_domain2);
  test_qd1_merges(query_domain1->child[0], query_domain2);
  test_qd1_merges(query_domain1->child[1], query_domain2);
}

test_qd2_merges(query_domain1, query_domain2)
     Query_domain *query_domain1, *query_domain2;
{
  if (!query_domain2) return;

  /* avoid testing pairs twice */

  if (query_domain1 < query_domain2
      && query_domain1 == query_domain1->parent
      && query_domain2 == query_domain2->parent) {
    compare_query_domains(query_domain1, query_domain2);
    /*
    fprintf(stderr, " (%d,%d)", 
	    query_domain1->best_score, query_domain2->best_score);
    */
  }

  test_qd2_merges(query_domain1, query_domain2->child[0]);
  test_qd2_merges(query_domain1, query_domain2->child[1]);
}

/* merge 1st into 2d (assumed higher scoring; both are equiv class reps) */
merge_query_domains(query_domain1, query_domain2)
     Query_domain *query_domain1, *query_domain2;
{
  int type;
  Score_hist *append_score_hist(), *s_h;

  query_domain1->parent = query_domain2;
  if (query_domain1->best_score == query_domain2->best_score) {
    query_domain2->n_best += query_domain1->n_best;
    if (query_domain1->next_best > query_domain2->next_best)
      query_domain2->next_best = query_domain1->next_best;
  }
  else if (query_domain1->best_score > query_domain2->next_best)
      query_domain2->next_best = query_domain1->best_score;

  /* keep best_pair, start, end from higher scoring domain */

  /* now merge score_hists */
  for (type = 0; type < 2; type++) {
    for (s_h = query_domain1->score_hist[type]; s_h; s_h = s_h->next) {
      append_score_hist(query_domain2, s_h->score, type, s_h->count);
    }
  }
}

compare_query_domains(query_domain1, query_domain2)
  Query_domain *query_domain1, *query_domain2;
{
  int start, end, intersect, i;
  Query_domain *qd[2];

  /*
  for (parent1 = query_domain1->parent; parent1 != parent1->parent; parent1 = parent1->parent);
  query_domain1->parent = parent1;
  for (parent2 = query_domain2->parent; parent2 != parent2->parent; parent2 = parent2->parent);
  query_domain2->parent = parent2;
  if (parent1 == parent2) return; /* already merged */

  start = query_domain1->start > query_domain2->start ? query_domain1->start : query_domain2->start;
  end = query_domain1->end > query_domain2->end ? query_domain2->end : query_domain1->end;
  intersect = (end - start + 1);
  if (intersect <= 0) return; /* is neg or 0, if no intersection */
  intersect *= mask_frac; 

  i = query_domain1->best_score < query_domain2->best_score;
  qd[!i] = query_domain1;
  qd[i] = query_domain2;
    
  if (intersect >= query_domain1->end - query_domain1->start  + 1 || intersect >= query_domain2->end - query_domain2->start + 1) {
    /* sig overlap -- so in same equiv class */
    merge_query_domains(qd[0], qd[1]);
  }
  else if (trim_flag && !is_unaligned(qd[0]->best_pair) && qd[0]->best_score < qd[1]->best_score) { /* require score strictly less */
    if (qd[0]->trim_start >= qd[1]->start && qd[0]->trim_start <= qd[1]->end) /* 2d cond'n superfluous? */
      qd[0]->trim_start = qd[1]->end + 1;
    if (qd[0]->trim_end <= qd[1]->end && qd[0]->trim_end >= qd[1]->start)
      qd[0]->trim_end = qd[1]->start - 1;
    if (qd[0]->trim_end < qd[0]->trim_start)
      merge_query_domains(qd[0], qd[1]);
  }
}

print_score_list(entry, query_domain)
     int entry;
     Query_domain *query_domain;
{
  char *get_id();
  /*  Score_list *s_h; */
  Score_hist *s_h2[2];
  int j, score, splice_flag, t_count, last_score, margin, margin2;

  if (!query_domain) return;

  print_score_list(entry, query_domain->child[0]);

  if (query_domain == query_domain->parent) {

    if (single_domain) {
      query_domain->start = 0;
      query_domain->end = get_seq_length(entry) - 1;
    }

    splice_flag = parameters->spliced_word_gapsize || parameters->spliced_word_gapsize2;
    printf("\n%s%s %d..%d all scores(counts): ", parameters->tags ? "SCORE_HISTOGRAM  " : "", get_id(entry), query_domain->start + 1, query_domain->end + 1);
    /*
      for (s_h = score_list; s_h->score; s_h++) {
      if (splice_flag)
      printf("%d(%d;%d) ", s_h->score, s_h->count[0], s_h->count[1]);
      else
      printf("%d(%d) ", s_h->score, s_h->count[0]);
      }

      printf("  :: ");
    */
    last_score = parameters->min_record_score;
    for (s_h2[0] = query_domain->score_hist[0], s_h2[1] = query_domain->score_hist[1]; s_h2[0] || s_h2[1]; ) {
      if (s_h2[0]) {
	score = s_h2[0]->score;
	if (splice_flag && s_h2[1] && score > s_h2[1]->score) score = s_h2[1]->score;
      }
      else score = s_h2[1]->score;
      printf("%d(", score);
      t_count = 0;
      for (j = 0; j < (splice_flag ? 2 : 1); j++) {
	if (j) printf(";");
	if (!s_h2[j] || score < s_h2[j]->score) printf("0");
	else {
	  printf("%d", s_h2[j]->count);
	  t_count += s_h2[j]->count;
	  s_h2[j] = s_h2[j]->next;
	}
      }
      printf(") ");
      margin = t_count == 1 ? score - last_score : 0;
      last_score = score > parameters->min_record_score ? score : parameters->min_record_score;
      if (!t_count) {
	fprintf(stderr, "\n%d  %d %d", score, s_h2[0] ? s_h2[0]->score : 0, s_h2[1] ? s_h2[1]->score : 0);
	fatalError("0 counts");
      }
    }
    /* printf(" %d %d", query_domain->best_score, query_domain->n_best); */
    /*
    margin2 = query_domain->n_best > 1 ? 0 : query_domain->best_score - query_domain->next_best;
    if (margin2 != margin) printf(" MARGIN %d", margin2);
    */
  }
  print_score_list(entry, query_domain->child[1]);

}

count_diffs(pair, sub_counts, indel_counts, pass, add_inform_subs, add_inform_indels, add_left_sub, add_right_sub) /* based on print_diffs() in qual.c  pair = pair_pointers[i_ptr]; */
     Aligned_pair *pair;
     int *sub_counts, *indel_counts;
     int pass;
     int *add_inform_subs, *add_inform_indels, *add_left_sub, *add_right_sub;
{
  unsigned char *diff, *diff1;
  int i, j, mult_flag, mid;
  int site1, site2, type;
  int d_site1, d_site2;
  unsigned char *get_seq();

  *add_inform_subs = *add_inform_indels = 0;
  *add_left_sub = pair->start1 - 1;
  *add_right_sub = pair->end1 + 1;
  site1 = pair->start1 - 1;
  site2 = pair->start2 - 1;
  mid = (pair->end1 + pair->start1) / 2;
  
  for (diff = pair->diffs; *(diff + 1); diff++) {
    site1 += diff_gap1(*diff);
    site2 += diff_gap2(*diff);
    if ((type = diff_type(*diff)) == 'M') continue;
    /* note reversal of type : so that refers to change in 1st sequence with respect to the 2d */
    mult_flag = 0;
    /*    q = qual[site1]; */
    diff1 = diff;
    if (type == 'D') {
      for (diff1 = diff + 1; diff_type(*diff1) == 'D' && !diff_gap2(*diff1); diff1++);
      i = diff1 - diff; 
      diff1--;
      if (i > 1) {
	mult_flag = 1;
	/* I: insertion, i = sizeof, site1 + 1) = printed location; site1 up to site1 + i - 1 are inserted bases in seq */;
	for (j = site1; j < site1 + i; j++) {
	  if (!pass) indel_counts[j] += 1;
	  else if (indel_counts[j] != indel_counts[j - 1]) *add_inform_indels += 1; /* condition implies this site is informative */
	}
     }
    }  
    else if (type == 'I') {
      for (diff1 = diff + 1; diff_type(*diff1) == 'I' && !diff_gap1(*diff1); diff1++);
      i = diff1 - diff; 
      diff1--;
      if (i > 1) {
	mult_flag = 1;
	/* deletion, i = sizeof, site1 + 1 = printed location */;
	if (!pass) 
	  indel_counts[site1] += 1;
	else if (indel_counts[site1] != indel_counts[site1 - 1]) 
	  *add_inform_indels += 1;
      }	
    }  
    if (!mult_flag) {	
      if (type == 'S') {
	if (!pass) 
	  sub_counts[site1] += 1;
	else if (sub_counts[site1] != sub_counts[site1 - 1]) {
	  *add_inform_subs += 1;
	  if (site1 < mid) *add_left_sub = site1;
	  else if (*add_right_sub > site1) *add_right_sub = site1;
	}
      }
      else {
	if (!pass) 
	  indel_counts[site1] += 1;
	else if (indel_counts[site1] != indel_counts[site1 - 1]) 
	  *add_inform_indels += 1;
      }
    }

      ; /* print -- reverse D and I -- site1 + 1 = printed location, seq[site1] = base */

   /* printed position in subject: is_reverse(pair) ? length2 - site2 : site2 + 1); where length2 = get_seq_length(pair->entry2); */
	   
    d_site1 = site1;
    d_site2 = site2;
    if (diff < diff1)
      for (diff++; diff <= diff1; diff++) {
	d_site1 += diff_gap1(*diff);
	d_site2 += diff_gap2(*diff);
      }
    /* printed bases
    for (i = site1 > 6 ? site1 - 6 : 0; i < length && i < d_site1 + 7; i++) 
      printf("%c", i < site1 || i > d_site1 ? tolower(seq[i]) : seq[i]);
    */

    diff = diff1;
    site1 = d_site1;
    site2 = d_site2;

  }
}

print_pair(pair, fp, no_zero)
     Aligned_pair *pair;
     FILE *fp;
     int no_zero;
{
  int length1;
  Aligned_pair *pair2;
  int mismatches, insertions, deletions;
  char *get_id();
  char *id1, *id2;

  id1 = get_id(pair->entry1);
  id2 = get_id(pair->entry2);
  if (no_zero && !pair->score) {
    fprintf(fp, "0  %s  %s", id1, id2);
    return;
  }

  length1 = pair->end1 - pair->start1 + 1;

  get_discreps(pair->diffs, &mismatches, &insertions, &deletions);

  fprintf(fp, "%s%s%4d %5.2f %4.2f %4.2f  %s    %5d %5d (%d)  %c %s  ", 
	  is_unaligned(pair) ? "SPLICED " : "",
	  parameters->tags ? "ALIGNMENT  " : "", /* indicate intron location AND which strand ss is on*/
	  pair->score, 
	  100.0 * (float)mismatches / length1,
	  100.0 * (float)insertions / length1,
	  100.0 * (float)deletions / length1,
	  id1, 1 + pair->start1, 1 + pair->end1,
	  get_seq_length(pair->entry1) - 1 - pair->end1,
	  is_reverse(pair) ? 'C' : ' ',
	  id2);

  if (is_reverse(pair)) {
    pair2 = pair->reversed_pair;
    fprintf(fp, " (%d) %5d %5d",
	   get_seq_length(pair->entry2) - 1 - pair2->end1,
	    1 + pair2->end1, 1 + pair2->start1);
  }
  else fprintf(fp, "  %5d %5d (%d)",
	       1 + pair->start2, 1 + pair->end2,
	       get_seq_length(pair->entry2) - 1 - pair->end2);
  fprintf(fp, " %c", pair->reject_flags < 2 ? '*' : ' ');
}

check_pairs()
{
  int entry1, length;
  Aligned_pair *pair, *pair_rev;
  Aligned_pair *get_aligned_pairs();

  notify("Checking pairs ...");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (!pair->score) continue;
      length = get_seq_length(entry1);
      if (!parameters->qual_scores && !parameters->matrix && (pair->end1 - pair->start1 + 1 < parameters->near_minscore ||
	  pair->end1 >= length))
	fatalError("Pair list corruption1"); /* assumes penalty = -1 */
      pair_rev = pair->reversed_pair;
      if (pair != pair_rev->reversed_pair)
	fatalError("Pair list corruption2");
      length = get_seq_length(pair_rev->entry1);
      if (!parameters->qual_scores && !parameters->matrix && (pair_rev->end1 - pair_rev->start1 + 1 < parameters->near_minscore ||
	  pair_rev->end1 >= length)) {
	fprintf(stderr,"\n%d %d %d ", pair_rev->start1, pair_rev->end1, length);
	fatalError("Pair list corruption3");
      }
    }
  }
  notify(" Done\n");
}

screen_seqs(db)
     Database *db;
{
  int entry1, i;
  Aligned_pair *pair;
  FILE *fq;
  FILE *fopenWrite();
  char file_name[100];
  Aligned_pair *get_aligned_pairs();
  unsigned char *seq;
  char *get_id(), *get_descrip();
  unsigned char *get_seq();

  strcpy(file_name, db->file->name);
  strcat(file_name,".screen");
  fq = fopenWrite(file_name);

  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    seq = get_seq(entry1);
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (pair->score < parameters->minscore) continue;
      for (i = pair->start1; i <= pair->end1; i++) seq[i] = 'X';
    }
    write_entry(fq, get_id(entry1), get_descrip(entry1), seq);
  }

  printf("\nScreened sequences written to  %s\n", file_name);
  fclose(fq);
}

/* find best pairs among all which overlap */

find_best()
{
  int entry1, better_flag, i, offset_diff;
  Aligned_pair *pair, *pair1, *pairs;
  Aligned_pair *get_aligned_pairs();
  int offset_diff_hist[100];

/*  fprintf(fp_log,"\n\nMultiple alignments:"); */

  for (i = 0; i < 100; i++) offset_diff_hist[i] = 0;

  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    pairs = get_aligned_pairs(entry1);
    for (pair = pairs; pair; pair = pair->next) {
      if (pair->entry1 >= pair->entry2) continue; 
      better_flag = 0;      
      for (pair1 = pairs; pair1; pair1 = pair1->next) {
	if (pair1->entry2 == pair->entry2 && pair1 != pair 
	    && is_reverse(pair) == is_reverse(pair1)) {
	  offset_diff = abs(pair->offset - pair1->offset);
	  offset_diff_hist[offset_diff > 99 ? 99 : offset_diff] += 1;
	  if (pair1->LLR_score >= pair->LLR_score) {
	    if (pair1->start1 < pair->end1 && pair1->end1 > pair->start1 && 
	      (pair1->LLR_score > pair->LLR_score || is_best(pair1) || 
	       pair1->LLR_score == pair->LLR_score && pair1->score > pair->score))
		better_flag= 1;
	  }
/*
	  if (pair1->w_start1 != pair->w_start1 && pair1->w_end1 != pair->w_end1) 
	    fprintf(fp_log, "\n%s %s  %d-%d %d-%d   %d-%d %d-%d    %d-%d %d-%d   %d-%d %d-%d  %d  %d ", 
		    align_entry->db_entry->id, 
		    align_array[pair->entry2].db_entry->id,
		    pair->w_start1, pair->w_end1, pair1->w_start1, pair1->w_end1, 
		    pair->w_start2, pair->w_end2, pair1->w_start2, pair1->w_end2,
		    pair->start1, pair->end1, pair1->start1, pair1->end1, 
		    pair->start2, pair->end2, pair1->start2, pair1->end2, pair->score, pair1->score);
*/	  
	}
      }
      set_best_flag(pair->reversed_pair, !better_flag);
      set_best_flag(pair, !better_flag);
    }
  }
  fprintf(fp_log,"\n\nHistogram of relative pair offsets, for read pairs with multiple alignments:");
  for (i = 0; i < 100; i++)
    if (offset_diff_hist[i])
	fprintf(fp_log,"\n%2d  %4d", i, offset_diff_hist[i]);
} 

/* an already used pair has overlapping extent -- indicates have tandem repeat;
	      so don't use this pair. Thus only the highest-scoring non-rejected overlap 
	      involving the tandem repeat will be used. */
/*
int test_tandem(pair)
     Aligned_pair *pair;
{
  Aligned_pair *pair1;

  for (pair1 = align_array[pair->entry1].first_pair; pair1; pair1 = pair1->next) {
    if (pair1->used && pair1->entry2 == pair->entry2) {
      if (pair1->start1 < pair->end1 && pair1->end1 > pair->start1)
	break;     
    }
  }
  return (pair1 ? pair1->score : 0);

}
*/

find_triple_rejects()
{
  int entry1;
  Aligned_pair *pair, *pair1, *pair2, *pairs;
  Aligned_pair *get_aligned_pairs();

  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    pairs = get_aligned_pairs(entry1);
    for (pair = pairs; pair; pair = pair->next) {
      if (pair->LLR_score < 0) continue;
      if (!is_best(pair)) continue;
      if (entry1 == pair->entry2) continue; 
      for (pair1 = pairs; pair1; pair1 = pair1->next) {
	if (pair1->LLR_score >= -10) continue;
	if (!is_best(pair1)) continue;
	if (pair1->entry1 == pair1->entry2 || pair1->entry2 == pair->entry2) continue; 
	for (pair2 = get_aligned_pairs(pair1->entry2); pair2; pair2 = pair2->next) 
	  if (pair2->entry2 == pair->entry2 && is_best(pair2) && pair2->LLR_score >= 0 
	      && is_reverse(pair2) == (is_reverse(pair) != is_reverse(pair1))) {
	    set_triple_reject_flag(pair, 1);
	    set_triple_reject_flag(pair1, 1);
	    set_triple_reject_flag(pair2, 1);
	    set_triple_reject_flag(pair->reversed_pair, 1);
	    set_triple_reject_flag(pair1->reversed_pair, 1);
	    set_triple_reject_flag(pair2->reversed_pair, 1);
/* SHOULD ALSO CHECK POSITIONING -- TO AVOID UNNECESSARY REJECTS ???? */
	  }
      }
    }
  }
}

/* find average offset for pair */
find_mean_offsets()
{
  int entry1;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int offset, d, match_size, g1;
  double t_offset;
  unsigned char *diff;

  notify("Finding mean offsets ... ");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      offset = pair->start1 - pair->start2;
      t_offset = match_size = 0;
      for (diff = pair->diffs; *diff; diff++) {
	g1 = diff_gap1(*diff);
	if ((d = diff_type(*diff)) == 'D') {
	  offset++;
	  g1--;
	}
	else if (d == 'I') offset--;
	match_size += g1;
	t_offset += g1 * offset;
      }
      if (offset != pair->end1 - pair->end2) fatalError("Offset discrepancy");
      pair->offset = t_offset / match_size;
    }
  }
/*
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    align_entry = align_array + entry1;
    for (pair = align_entry->first_pair; pair; pair = pair->next) {
      if ((d = abs(pair->offset + pair->reversed_pair->offset)) > 1)
	fprintf(stderr,"\n %d %d", d, is_reverse(pair));
    }
  }
*/
  notify(" Done\n");
}


/* sets align_entry->first_start, ->last_end,  ->segments, ->fwd_conf_segs, and ->rev_conf_segs;
   depends upon pair->hq_mismatch (plus database, etc.)
*/

/* note first_start, last_end get changed later (in find_extents); may not be 
necessary. Note also set orig_qual to 0 in  initial part of read (before first_start
since it is likely to be vector */
 
find_segments(db)
     Database *db;
{
  Align_info *align_entry;
  Align_info *get_align_entry();
  int entry1;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Segment *insert_segment();
  Segment *segment;

  notify("Finding segments ... ");
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    align_entry->segments = 0;
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (same_subclone(pair->entry1, pair->entry2) 
/* FOLLOWING TAKEN OUT 7/20/98 */
	  /* && !different_chemistry(pair->entry1, pair->entry2) */) continue; 
      align_entry->segments =
	insert_segment(align_entry->segments, pair->start1 + EXTENT_FUDGE, pair->end1 - EXTENT_FUDGE);
    }
    for (segment = align_entry->segments; segment; segment = segment->next){
      segment->start -= EXTENT_FUDGE;
      segment->end += EXTENT_FUDGE;
    }
    /* add and subtract EXTENT_FUDGE to compensate for spurious matching nucs at ends of
       true matches */
  }
  notify(" Done\n");
}

find_starts_ends(db, chim_test)
     Database *db;
     int chim_test;
{
  Align_info *align_entry;
  Align_info *get_align_entry();
  int entry1;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Segment *insert_segment();

  notify("Finding starts/ends ... ");
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    align_entry->last_end = align_entry->rev_last_end = 0;
    align_entry->first_start = align_entry->rev_first_start = get_seq_length(entry1);
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (same_subclone(pair->entry1, pair->entry2) 
/* FOLLOWING TAKEN OUT 7/20/98 */
	  /* && !different_chemistry(pair->entry1, pair->entry2) */) continue; 
      if (chim_test && is_reject_chimeric(pair)) continue;
      if (pair->start1 < align_entry->first_start)
	align_entry->first_start = pair->start1;
      if (pair->end1 > align_entry->last_end)
	align_entry->last_end = pair->end1;
      if (is_reverse(pair)) {
	if (pair->start1 < align_entry->rev_first_start)
	  align_entry->rev_first_start = pair->start1;
	if (pair->end1 > align_entry->rev_last_end)
	  align_entry->rev_last_end = pair->end1;
      }
    }
  }
  notify(" Done\n");
}

old_find_extents(db, pass_level)
     Database *db;
     int pass_level;
{
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();
  int entry1, length1;
  unsigned char *get_seq(), *get_comp_seq();
  char *get_orig_qual(), *get_adj_qual(), *get_id(), *get_descrip();
  char *oqual1, *oqual2;
  static char *qual_buffer;
  static int qual_buffer_length;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int i, j, q1, q2, q, x, length2, last_disc, last_disc2, site1;
  unsigned char *diff;
  unsigned char *seq, *seq2;
  int start0, end0, g1, g2, last_start, strong_score, max_score, start, end, start1, end1, start2, end2;
  int left_first_start, left_last_end;
  int right_first_start, right_last_end;
  int first_start, rev_first_start, last_end, rev_last_end;
  int str_confirm;
  int diff_site, diff_site2, d;
  unsigned char *last_mononuc(), *first_mononuc();
  unsigned char *weak_last_mononuc(), *weak_first_mononuc();
  char *adj_qual;
  static int *best_same_match, *best_opp_match, *best_same_discrep, *best_opp_discrep, *test_qual;
  static int *n_same_matches, *n_opp_matches, *n_w_same_matches, *n_w_opp_matches;
  int *best_discrep, *best_match, *n_matches, *n_w_matches;
  static int vec_length;
  char *our_alloc();
  int w_last_disc, next_w_last_disc, w_last_disc2, w_start, next_w_start, w_end, j_start, j_end; 
  int confirm_flag, w_site1;
  int nuc, k_start, k_end, n_m;
  int create_test_file;
  char file_name[50], buf[10];
  FILE *fopenWrite();
  FILE *fp, *fq;

  notify("Finding extents ... ");

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage");
  if (pass_level == 2 && parameters->create_test_file) {
    create_test_file = parameters->create_test_file;
    sprintf(buf, ".%d", create_test_file);
    strcpy(file_name, "test_qual.fasta");
    strcat(file_name, buf);
    fp = fopenWrite(file_name);
    strcat(file_name, ".qual");
    fq = fopenWrite(file_name);
  }
  else create_test_file = 0;

  mark_save_block(); /* temporary segments will be created; want to eliminate later */
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    seq = get_seq(entry1);
    length1 = get_seq_length(entry1);
    oqual1 = get_orig_qual(entry1);
    if (vec_length < length1 + 3) {
      if (vec_length) {
	our_free(best_same_match - 1);
	our_free(best_same_discrep - 1);
	our_free(best_opp_match - 1);
	our_free(best_opp_discrep - 1);
	our_free(test_qual - 1);
	our_free(n_same_matches - 1);
	our_free(n_w_same_matches - 1);
	our_free(n_opp_matches - 1);
	our_free(n_w_opp_matches - 1);
      } 
      vec_length = length1 + 3;
      best_same_match = 1 + (int *)our_alloc(vec_length * sizeof(int));
      best_same_discrep = 1 + (int *)our_alloc(vec_length * sizeof(int));
      best_opp_match = 1 + (int *)our_alloc(vec_length * sizeof(int));
      best_opp_discrep = 1 + (int *)our_alloc(vec_length * sizeof(int));
      test_qual = 1 + (int *)our_alloc(vec_length * sizeof(int));
      n_same_matches = 1 + (int *)our_alloc(vec_length * sizeof(int));
      n_w_same_matches = 1 + (int *)our_alloc(vec_length * sizeof(int));
      n_opp_matches = 1 + (int *)our_alloc(vec_length * sizeof(int));
      n_w_opp_matches = 1 + (int *)our_alloc(vec_length * sizeof(int));
    }
    for (i = 0; i < length1; i++) {
      best_opp_match[i] = best_same_match[i] =  0;
      best_same_discrep[i] = best_opp_discrep[i] = 0;
      n_same_matches[i] = n_w_same_matches[i] = n_opp_matches[i] = n_w_opp_matches[i] = 0;
    }
    left_first_start = right_first_start
      = first_start = rev_first_start = length1;
    left_last_end = right_last_end = last_end = rev_last_end = 0;

    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      str_confirm = is_reverse(pair) || different_chemistry(pair->entry1, pair->entry2);
      if (same_subclone(pair->entry1, pair->entry2) && !str_confirm) continue; 
/* NOTE THAT DYE-TERMINATOR READ MAY SPURIOUSLY CONFIRM LINK FROM CHIMERIC SUBCLONE -- ALTHO NOT
   LIKELY */
/* RESTORE?
*/ 
      if (pass_level && pair->LLR_score < 0 && !is_used(pair)) continue; 
      if (pass_level == 2 && !is_used(pair)) continue; 
      /* check whether this restriction is appropriate -- e.g. for HQ mismatch detection ??? */
      if (!strcmp(get_id(entry1), "l177p_337.s1"))
	  fprintf(stderr, "\n%d  %s  %s", pass_level, get_id(entry1), get_id(pair->entry2));
/* Is following necessary??
      if (is_reject_any(pair)) {
	if (reject_left(pair)) {
	  if (pair->start1 < left_first_start) left_first_start = pair->start1;
	  if (pair->end1 > left_last_end) left_last_end = pair->end1;
	}
	if (reject_right(pair)) {
	  if (pair->start1 < right_first_start) right_first_start = pair->start1;
	  if (pair->end1 > right_last_end) right_last_end = pair->end1;
	}
	continue;
      }
*/
      best_discrep = str_confirm ? best_opp_discrep : best_same_discrep;
      best_match = str_confirm ? best_opp_match : best_same_match;
      n_matches = str_confirm ? n_opp_matches : n_same_matches;
      n_w_matches = str_confirm ? n_w_opp_matches : n_w_same_matches;

      if (create_test_file) {
	if ((create_test_file % 2) == str_confirm) 
	  for (i = 0; i < length1; i++) test_qual[i] = 0;
	else continue;
      }
	
      align_entry2 = get_align_entry(pair->entry2);
      seq2 = is_reverse(pair) ? get_comp_seq(pair->entry2) : get_seq(pair->entry2);
      length2 = get_seq_length(pair->entry2);
      oqual2 = get_orig_qual(pair->entry2);
      if (is_reverse(pair)) {
	if (length2 + 3 > qual_buffer_length) {
	  if (qual_buffer_length) our_free(qual_buffer - 1);
	  qual_buffer_length = length2 + 3;
	  qual_buffer = 1 + (char *)our_alloc(qual_buffer_length * sizeof(char));
	}
	for (i = 0; i < length2; i++) 
	  qual_buffer[i] = oqual2[length2 - 1 - i];
	oqual2 = qual_buffer;  
      }

      diff_site = pair->start1 - 1;
      last_disc = last_mononuc(seq + diff_site) - seq;
      w_last_disc = weak_last_mononuc(seq + diff_site) - seq; /* diff_site; */
      w_start = w_last_disc + 1;
      diff_site2 = pair->start2 - 1;
      last_disc2 = diff_site2 + (last_disc - diff_site);
      w_last_disc2 = diff_site2 + (w_last_disc - diff_site);

      last_start = diff_site + 1;
      strong_score = max_score = 0;
      start0 = end0 = -1;
/*
      start0 = align_entry->db_entry->length;      
      end0 = 0;
*/
      for (diff = pair->diffs; *diff; diff++) {
	g1 = diff_gap1(*diff);
	g2 = diff_gap2(*diff);
	d = diff_type(*diff);
	strong_score += (g1 > g2 ? g1 : g2) - (d != 'M'); 
	diff_site += g1;
	diff_site2 += g2;
	if (d == 'M') continue;

	if (strong_score >= max_score) {
	  start0 = last_start;
	  end0 = d == 'I' ? diff_site : diff_site - 1;
	  max_score = strong_score;
	}
	if (strong_score < -parameters->confirm_penalty) {
	  strong_score = 0;
	  last_start = diff_site + 1;
	}
	else strong_score += parameters->confirm_penalty;

	site1 = first_mononuc(seq + diff_site + (d == 'I')) - seq;
	w_site1 = weak_first_mononuc(seq + diff_site + (d == 'I')) - seq; /* diff_site + (d == 'I') */
	w_end = w_site1 - 1; /* (d == 'I' || d == 'D' || !*(diff + 1)) ? w_site1 - 1 : diff_site - 1; */
	next_w_last_disc = weak_last_mononuc(seq + diff_site) - seq; /* diff_site; */
	next_w_start = next_w_last_disc + 1; /* (d == 'I' || d == 'D') ? next_w_last_disc + 1 : diff_site + 1; */
      
	if (diff_site2 < length2) {
	  q2 = oqual2[diff_site2];
	  if (q2 >= 98) q2 = 0; /* edited bases get quality 0, for this purpose */
	  if (d == 'D') {
	    nuc = seq[diff_site];
	    for (j_start = diff_site - 1; seq[j_start] == nuc; j_start--);
	    j_start++;
	    for (j_end = diff_site + 1; seq[j_end] == nuc; j_end++);
	    j_end--;
	    k_start = diff_site2;
	    do {
	      q = oqual2[k_start];
	      if (q < q2) q2 = q;
	    } while (seq2[k_start--] == nuc);

	    k_end = diff_site2 + 1;
	    do {
	      q = oqual2[k_end];
	      if (q < q2) q2 = q;
	    } while (seq2[k_end++] == nuc);

/* change -- to only affect 1st, last bases? */
	  }
	  else if (d == 'I') {
	    nuc = seq2[diff_site2];
	    for (j_start = diff_site; seq[j_start] == nuc; j_start--);
	    for (j_end = diff_site + 1; seq[j_end] == nuc; j_end++);

	    for (k_start = diff_site2 - 1; seq2[k_start] == nuc; k_start--) {
	      q = oqual2[k_start];
	      if (q < q2) q2 = q;
	    }

	    for (k_end = diff_site2 + 1; seq2[k_end] == nuc; k_end++) {
	      q = oqual2[k_end];
	      if (q < q2) q2 = q;
	    }

	/*     j_start =  w_site1; diff_site; */
	/*     j_end = next_w_last_disc; diff_site + (d == 'I'); */
	  }
	  else {
	    j_start = diff_site;
	    j_end = diff_site;
	  }
	  if (j_start < 0) j_start = 0;
	  if (j_end >= length1) j_end = length1 - 1;

	  for (i = j_start; i <= j_end; i++) 
	    if (best_discrep[i] < q2) best_discrep[i] = q2;
	}

	confirm_flag = site1 > last_disc + 5; /* parameters->confirm_length; */
	if (w_site1 > w_last_disc + 3) { 
	  start = last_disc + parameters->confirm_trim + 1; 
	  start2 = last_disc2 + parameters->confirm_trim + 1; 
	  end = site1 - parameters->confirm_trim - 1;
	  end2 = start2 + (end - start);
	  w_start = w_last_disc + 1;
/*
	  start = (diff - 1)->site1 + parameters->confirm_trim + 1;
	  if (last_disc + 1 > start) start = last_disc + 1;
	  end = diff->site1 - parameters->confirm_trim - 1;
	  if (site1 - 1 < end) end = site1 - 1;
*/
/*
	  if (start < start0) start0 = start; 
	  end0 = end;
*/
/* a second reverse confirming read for a base counts as a forward read (to
   symmetrize the 2 + 1 rule) */

	  if (confirm_flag) {
	    for (i = w_start; i < start; i++) n_w_matches[i] += 1;
	    for (i = start, j = start2; i <= end; i++, j++) {
	      q2 = oqual2[j];
	      if (q2 < 98 && best_match[i] < q2) best_match[i] = q2;
	      n_matches[i] += 1;
	    }
	    for (i = end + 1; i <= w_end; i++) n_w_matches[i] += 1;
	    if (str_confirm && (create_test_file == 1 || create_test_file == 3)) {
	      for (i = start, j = start2; i <= end; i++, j++) {
		q2 = oqual2[j];
		if (q2 >= 98) continue;
		q1 = oqual1[i];
		if (create_test_file == 3 && q1 && q2) continue;
		if (create_test_file == 1 && (!q1 || !q2)) continue;
		if (!q1) q1 = 1;
		if (!q2) q2 = 1;
		test_qual[i] = q1 + q2 < 90 ? q1 + q2 : 90;
	      }
	    }
	    if (!str_confirm && (create_test_file == 2 || create_test_file == 4)) {
	      for (i = start, j = start2; i <= end; i++, j++) {
		q2 = oqual2[j];
		if (q2 >= 98) continue;
		q1 = oqual1[i];
		if (create_test_file == 4 && q1 && q2) continue;
		if (create_test_file == 2 && (!q1 || !q2)) continue;
		if (!q1) q1 = 1;
		if (!q2) q2 = 1;
		test_qual[i] = q1 < q2 ? q2 : q1;
	      }
	    }
	  }
	  else {
	    for (i = w_start; i <= w_end; i++) n_w_matches[i] += 1;
	  }
	  if (str_confirm && create_test_file == 5) {
	    for (i = w_start; i <= w_end; i++) {
	      if (!confirm_flag || i < start || i > end) {
		q1 = oqual1[i];
		test_qual[i] = q1 + 15;
	      }
	    }
	  }
	  if (!str_confirm && create_test_file == 6) {
	    for (i = w_start; i <= w_end; i++) {
	      if (!confirm_flag || i < start || i > end) {
		q1 = oqual1[i];
		if (!q1) q1 = 1;
		test_qual[i] = q1;
	      }
	    }
	  }
	}
	last_disc = last_mononuc(seq + diff_site) - seq;
	last_disc2 = diff_site2 + (last_disc - diff_site);
	w_last_disc = next_w_last_disc; /* weak_last_mononuc(seq + diff_site) - seq;  diff_site; */
	w_last_disc2 = diff_site2 + (w_last_disc - diff_site);
	w_start = next_w_start;
      }
/* note that confirm_score not being used as criterion in the above */
/*
      if (strong_score >= parameters->confirm_score && !is_reject_chimeric(pair)) {
*/
      if (!is_reject_chimeric(pair)) {
	if (str_confirm && start0 < rev_first_start) rev_first_start = start0;
	if (str_confirm && end0 > rev_last_end) rev_last_end = end0;
	if (start0 < first_start) first_start = start0;
	if (end0 > last_end) last_end = end0;
      }
      if (create_test_file) {
	for (j = 0; j < length1 && !test_qual[j]; j++);
	if (j < length1) {
	  fprintf(fp,">%s", get_id(entry1));
	  for (j = 0; j < length1; j++) {
	    if (!(j % 50)) fprintf(fp,"\n");
	    fprintf(fp,"%c", seq[j]);
	  }
	  fprintf(fp,"\n");
	  fprintf(fq, ">%s", get_id(entry1));
	  for (j = 0; j < length1; j++) {
	    if (!(j % 50)) fprintf(fq,"\n");
	    fprintf(fq, "%d ", test_qual[j]);
	  }
	  fprintf(fq,"\n");
	}
      }
    }

    adj_qual = get_adj_qual(entry1);
    for (i = 0; i <= align_entry->last_vec; i++) adj_qual[i] = 0;
      
    for (; i < length1; i++) {
      q1 = oqual1[i];
      if (q1 == 99) adj_qual[i] = 99;
      else if (seq[i] == 'N' || seq[i] == 'X') adj_qual[i] = 0;
      else {
/*
*/
	if (n_same_matches[i] + n_w_same_matches[i] && best_same_match[i] < 15) 
	  best_same_match[i] = 15; 

	if (n_opp_matches[i] + n_w_opp_matches[i] && best_opp_match[i] < 15) 
	  best_opp_match[i] = 15; 
/*
	if (best_same_discrep[i] < 15) best_same_discrep[i] = 0;
	if (best_opp_discrep[i] < 15) best_opp_discrep[i] = 0;
*/
	if (q1 == 98) q1 = 0; /* edited low quality base */

	x = best_same_discrep[i] + best_opp_discrep[i]; /* + 5; */
	n_m = n_same_matches[i] + n_w_same_matches[i] + n_opp_matches[i] + n_w_opp_matches[i];
	if (i >= first_start && i <= last_end || n_m) {
/*
*/
	  if (q1 < 15) q1 = 15;
	  q1 = q1 > best_same_match[i] ? q1 : best_same_match[i];
/*
	  q1 += n_m > 2 ? 15 : n_m * 7;
*/
	  q1 += best_opp_match[i];
	}
	else {
	  if (x < 25) x = 25;
	}
/*
	x = q1 - x;
	q1 += x > 0 ? 0 : x < -20 ? -20 : x;
*/
	if (x > q1) x = q1;
	if (/* pass_level < 2 && */ x > 20) x = 20;
	q1 -= x;
	adj_qual[i] = q1 < 2 ? 2 : q1 > 90 ? 90 : q1;
	/*  if (pass_level > 1) q1 -= best_same_discrep[i] + best_opp_discrep[i]; */
      }
    }

/*
    for (segment = align_entry->segments, i_bit = 1; segment; 
       segment = segment->next, i_bit *= 2) 
      if (i_bit & align_entry->chimera_bits) 
	for (i = segment->start; i <= segment->end; i++)
	  adj_qual[i] = 0;
*/

/* RESTORE?
*/
    align_entry->first_start = first_start;
    align_entry->last_end = last_end;

    align_entry->rev_first_start = rev_first_start;
    align_entry->rev_last_end = rev_last_end;

/*
    if (left_first_start < weak_first_start - parameters->maxgap)
      printf("\n%s   %d-%d    confirmed only by left-incomplete matches", 
	     align_entry->db_entry->id, left_first_start + 1, weak_first_start);
    if (left_last_end > weak_last_end + parameters->maxgap)
      printf("\n%s   %d-%d    confirmed only by left-incomplete matches", 
	     align_entry->db_entry->id, weak_last_end + 2, left_last_end + 1);
    if (right_first_start < weak_first_start - parameters->maxgap)
      printf("\n%s   %d-%d    confirmed only by right-incomplete matches", 
	     align_entry->db_entry->id, right_first_start + 1, weak_first_start);
    if (right_last_end > weak_last_end + parameters->maxgap)
      printf("\n%s   %d-%d    confirmed only by right-incomplete matches", 
	     align_entry->db_entry->id, weak_last_end + 2, right_last_end + 1);
*/
  }
  free_seg_blocks();
  notify(" Done\n");
  if (create_test_file) {
    fclose(fp);
    fclose(fq);
    exit(1);
  }
}


/* co-ordinates of best_discrep and best_match:
   1st coord indicates whether inside (0) or outside (1) the quality trimmed region
       of the matching read;
   2d coord indicates whether is highest (0) or 2d highest (1) quality;
   3d coord indicates whether
             for best_match: whether aligned base is solid (0) or weak (1)
             for best_discrep: whether discrepancy is deletion, insertion or substitution
   4th coord gives quality value.
*/

find_extents(db, pass_level)
     Database *db;
     int pass_level;
{
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();
  int entry1, length1;
  unsigned char *get_seq(), *get_comp_seq();
  char *get_orig_qual(), *get_adj_qual(), *get_id(), *get_descrip();
  char *oqual1, *oqual2;
  static char *qual_buffer;
  static int qual_buffer_length;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int i, j, k, m, mm, n, q1, q2, q, x, length2, length2m1, last_disc, last_disc2, site1;
  unsigned char *diff;
  unsigned char *seq, *seq2;
  int start0, end0, g1, g2, last_start, strong_score, max_score, start, end, start1, end1, start2, end2;
  int left_first_start, left_last_end;
  int right_first_start, right_last_end;
  int first_start, rev_first_start, last_end, rev_last_end;
  int str_confirm;
  int diff_site, diff_site2, d;
  unsigned char *last_mononuc(), *first_mononuc();
  unsigned char *weak_last_mononuc(), *weak_first_mononuc();
  char *adj_qual;
  static int *test_qual;
  static int ****best_same_match, ****best_opp_match, 
      *****best_same_discrep, *****best_opp_discrep;
  static int *n_same_matches, *n_opp_matches, *n_w_same_matches, *n_w_opp_matches;
  int *n_matches, *n_w_matches;
  int *****best_discrep, ****best_match;
  static int vec_length;
  char *our_alloc();
  int w_last_disc, next_w_last_disc, w_last_disc2, w_start, next_w_start, w_end, j_start, j_end; 
  int confirm_flag, w_site1;
  int nuc, k_start, k_end, n_m;
  int create_test_file;
  int w_start2, y, max_x, max_x2, max_y, max_y2, k2, q_flag;
  int q_start2, q_end2;
  char file_name[50], buf[10];
  FILE *fopenWrite();
  FILE *fp, *fq;
  int nuc_lookup[256];
  int chemistry;
  int seg_start, seg_end;
  Segment *segment;

  notify("Finding extents ... ");

  if (parameters->DNA_flag >= 3) fatalError("illicit compressed quality usage");

  for (i = 0; i < 256; i++) nuc_lookup[i] = 0;
  nuc_lookup['A'] = 0;
  nuc_lookup['C'] = 1;
  nuc_lookup['G'] = 2;
  nuc_lookup['T'] = 3;
  
  if (!vec_length) {
    best_same_match = (int ****)our_alloc(2 * sizeof(int ***));
    best_opp_match = (int ****)our_alloc(2 * sizeof(int ***));
    best_same_discrep = (int *****)our_alloc(2 * sizeof(int ****));
    best_opp_discrep = (int *****)our_alloc(2 * sizeof(int ****));
    
    for (k = 0; k < 2; k++) {
      best_same_match[k] = (int ***)our_alloc(2 * sizeof(int **));
      best_opp_match[k] = (int ***)our_alloc(2 * sizeof(int **));
      best_same_discrep[k] = (int ****)our_alloc(2 * sizeof(int ***));
      best_opp_discrep[k] = (int ****)our_alloc(2 * sizeof(int ***));
      for (i = 0; i < 2; i++) {
	best_same_match[k][i] = (int **)our_alloc(2 * sizeof(int *));
	best_opp_match[k][i] = (int **)our_alloc(2 * sizeof(int *));
	best_same_discrep[k][i] = (int ***)our_alloc(4 * sizeof(int **)); /* corrected from 3 per James Bonfield's comment */
	best_opp_discrep[k][i] = (int ***)our_alloc(4 * sizeof(int **)); 
	for (m = 0; m < 4; m++) {
	  best_same_discrep[k][i][m] = (int **)our_alloc(4 * sizeof(int *));
	  best_opp_discrep[k][i][m] = (int **)our_alloc(4 * sizeof(int *));
	}
      }
    }
  }

  if (pass_level == 2 && parameters->create_test_file) {
    create_test_file = parameters->create_test_file;
    sprintf(buf, ".%d", create_test_file);
    strcpy(file_name, "test_qual.fasta");
    strcat(file_name, buf);
    fp = fopenWrite(file_name);
    strcat(file_name, ".qual");
    fq = fopenWrite(file_name);
  }
  else create_test_file = 0;

  mark_save_block(); /* temporary segments will be created; want to eliminate later */
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    chemistry = align_entry->chemistry;
    seq = get_seq(entry1);
    length1 = get_seq_length(entry1);
    oqual1 = get_orig_qual(entry1);
    if (vec_length < length1 + 3) {
      if (vec_length) {
	for (k = 0; k < 2; k++) {
	  for (i = 0; i < 2; i++) {
	    for (j = 0; j < 2; j++) {
	      our_free(best_same_match[k][i][j] - 1);
	      our_free(best_opp_match[k][i][j] - 1);
	    }
	    for (j = 0; j < 3; j++) {
	      for (m = 0; m < 4; m++) {
		our_free(best_same_discrep[k][i][j][m] - 1);
		our_free(best_opp_discrep[k][i][j][m] - 1);
	      }
	    }
	  }
	}
	our_free(test_qual - 1);
	our_free(n_same_matches - 1);
	our_free(n_w_same_matches - 1);
	our_free(n_opp_matches - 1);
	our_free(n_w_opp_matches - 1);
      } 
      vec_length = length1 + 3;
      for (k = 0; k < 2; k++) {
	for (i = 0; i < 2; i++) {
	  for (j = 0; j < 2; j++) {
	    best_same_match[k][i][j] = 1 + (int *)our_alloc(vec_length * sizeof(int));
	    best_opp_match[k][i][j] = 1 + (int *)our_alloc(vec_length * sizeof(int));
	  }
	  for (j = 0; j < 3; j++) {
	    for (m = 0; m < 4; m++) {
	      best_same_discrep[k][i][j][m] = 1 + (int *)our_alloc(vec_length * sizeof(int));
	      best_opp_discrep[k][i][j][m] = 1 + (int *)our_alloc(vec_length * sizeof(int));
	    }
	  }
	}
      }
      test_qual = 1 + (int *)our_alloc(vec_length * sizeof(int));
      n_same_matches = 1 + (int *)our_alloc(vec_length * sizeof(int));
      n_w_same_matches = 1 + (int *)our_alloc(vec_length * sizeof(int));
      n_opp_matches = 1 + (int *)our_alloc(vec_length * sizeof(int));
      n_w_opp_matches = 1 + (int *)our_alloc(vec_length * sizeof(int));
    }
    for (k = 0; k < 2; k++) {
      for (i = 0; i < 2; i++) {
	for (j = 0; j < 2; j++) {
	  for (m = 0; m < length1; m++) {
	    best_same_match[k][i][j][m] = best_opp_match[k][i][j][m] = 0;
	  }
	}
	for (j = 0; j < 3; j++) {
	  for (m = 0; m < 4; m++) 
	    for (n = 0; n < length1; n++) 
	      best_same_discrep[k][i][j][m][n] = best_opp_discrep[k][i][j][m][n] = 0;
	}
      }
    }
    for (m = 0; m < length1; m++) 
      if (oqual1[m] < 98) {
/*
	k = m < align_entry->qual_start || m > align_entry->qual_end;
	best_same_match[k][0][0][m] = oqual1[m];
*/
	best_same_match[0][0][0][m] = oqual1[m];
      }

    for (i = 0; i < length1; i++) 
      n_same_matches[i] = n_w_same_matches[i] = n_opp_matches[i] = n_w_opp_matches[i] = 0;
    left_first_start = right_first_start
      = first_start = rev_first_start = length1;
    left_last_end = right_last_end = last_end = rev_last_end = 0;

    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (pass_level && pair->LLR_score < 0 && !is_used(pair)) continue; 
      if (pass_level == 2 && !is_used(pair)) continue; 
      /* check whether this restriction is appropriate -- e.g. for HQ mismatch detection ??? */
      str_confirm = is_reverse(pair) || different_chemistry(pair->entry1, pair->entry2);
      if (same_subclone(pair->entry1, pair->entry2) && !str_confirm) continue; 
/* NOTE THAT DYE-TERMINATOR READ MAY SPURIOUSLY CONFIRM LINK FROM CHIMERIC SUBCLONE -- ALTHO NOT
   LIKELY */

/* Is following necessary??
      if (is_reject_any(pair)) {
	if (reject_left(pair)) {
	  if (pair->start1 < left_first_start) left_first_start = pair->start1;
	  if (pair->end1 > left_last_end) left_last_end = pair->end1;
	}
	if (reject_right(pair)) {
	  if (pair->start1 < right_first_start) right_first_start = pair->start1;
	  if (pair->end1 > right_last_end) right_last_end = pair->end1;
	}
	continue;
      }
*/
      best_discrep = str_confirm ? best_opp_discrep : best_same_discrep;
      best_match = str_confirm ? best_opp_match : best_same_match;
      n_matches = str_confirm ? n_opp_matches : n_same_matches;
      n_w_matches = str_confirm ? n_w_opp_matches : n_w_same_matches;

      if (create_test_file) {
	if ((create_test_file % 2) == str_confirm) 
	  for (i = 0; i < length1; i++) test_qual[i] = 0;
	else continue;
      }
	
      align_entry2 = get_align_entry(pair->entry2);
      seq2 = is_reverse(pair) ? get_comp_seq(pair->entry2) : get_seq(pair->entry2);
      length2 = get_seq_length(pair->entry2);
      length2m1 = length2 - 1;
      oqual2 = get_orig_qual(pair->entry2);
      q_start2 = is_reverse(pair) ? length2m1 - align_entry2->qual_end : align_entry2->qual_start;
      q_end2 = is_reverse(pair) ? length2m1 - align_entry2->qual_start : align_entry2->qual_end;
      if (is_reverse(pair)) {
	if (length2 + 3 > qual_buffer_length) {
	  if (qual_buffer_length) our_free(qual_buffer - 1);
	  qual_buffer_length = length2 + 3;
	  qual_buffer = 1 + (char *)our_alloc(qual_buffer_length * sizeof(char));
	}
	for (i = 0; i < length2; i++) 
	  qual_buffer[i] = oqual2[length2m1 - i];
	oqual2 = qual_buffer;  
      }

      diff_site = pair->start1 - 1;
      last_disc = last_mononuc(seq + diff_site) - seq;
      w_last_disc = weak_last_mononuc(seq + diff_site) - seq; /* diff_site; */
      w_start = w_last_disc + 1;
      diff_site2 = pair->start2 - 1;
      last_disc2 = diff_site2 + (last_disc - diff_site);
      w_last_disc2 = diff_site2 + (w_last_disc - diff_site);

      last_start = diff_site + 1;
      strong_score = max_score = 0;
      start0 = end0 = -1;
/*
      start0 = align_entry->db_entry->length;      
      end0 = 0;
*/
      for (diff = pair->diffs; *diff; diff++) {
	g1 = diff_gap1(*diff);
	g2 = diff_gap2(*diff);
	d = diff_type(*diff);
	strong_score += (g1 > g2 ? g1 : g2) - (d != 'M'); 
	diff_site += g1;
	diff_site2 += g2;
	if (d == 'M') continue;

	if (strong_score >= max_score) {
	  start0 = last_start;
	  end0 = d == 'I' ? diff_site : diff_site - 1;
	  max_score = strong_score;
	}
	if (strong_score < -parameters->confirm_penalty) {
	  strong_score = 0;
	  last_start = diff_site + 1;
	}
	else strong_score += parameters->confirm_penalty;

	site1 = first_mononuc(seq + diff_site + (d == 'I')) - seq;
	w_site1 = weak_first_mononuc(seq + diff_site + (d == 'I')) - seq; /* diff_site + (d == 'I') */
	w_end = w_site1 - 1; /* (d == 'I' || d == 'D' || !*(diff + 1)) ? w_site1 - 1 : diff_site - 1; */
	next_w_last_disc = weak_last_mononuc(seq + diff_site) - seq; /* diff_site; */
	next_w_start = next_w_last_disc + 1; /* (d == 'I' || d == 'D') ? next_w_last_disc + 1 : diff_site + 1; */
      
	if (diff_site2 < length2) {
	  q2 = oqual2[diff_site2];
	  if (q2 >= 98) q2 = 0; /* edited bases get quality 0, for this purpose */
	  k2 = diff_site2 < q_start2 || diff_site2 > q_end2;
	  if (d == 'D') {
	    nuc = seq[diff_site];
	    for (j_start = diff_site - 1; seq[j_start] == nuc; j_start--);
	    j_start++;
	    for (j_end = diff_site + 1; seq[j_end] == nuc; j_end++);
	    j_end--;
	    k_start = diff_site2;
	    do {
	      q = oqual2[k_start];
	      if (q < q2) q2 = q;
	    } while (seq2[k_start--] == nuc);

	    k_end = diff_site2 + 1;
	    do {
	      q = oqual2[k_end];
	      if (q < q2) q2 = q;
	    } while (seq2[k_end++] == nuc);

/* change -- to only affect 1st, last bases? */
	  }
	  else if (d == 'I') {
	    nuc = seq2[diff_site2];
	    for (j_start = diff_site; seq[j_start] == nuc; j_start--);
	    for (j_end = diff_site + 1; seq[j_end] == nuc; j_end++);

	    for (k_start = diff_site2 - 1; seq2[k_start] == nuc; k_start--) {
	      q = oqual2[k_start];
	      if (q < q2) q2 = q;
	    }

	    for (k_end = diff_site2 + 1; seq2[k_end] == nuc; k_end++) {
	      q = oqual2[k_end];
	      if (q < q2) q2 = q;
	    }

	/*     j_start =  w_site1; diff_site; */
	/*     j_end = next_w_last_disc; diff_site + (d == 'I'); */
	  }
	  else {
	    j_start = diff_site;
	    j_end = diff_site;
	  }
	  if (j_start < 0) j_start = 0;
	  if (j_end >= length1) j_end = length1 - 1;

	  mm = nuc_lookup[seq2[diff_site2]];
	  k = d == 'D' ? 0 : (d == 'I' ? 1 : 2);
	  for (i = j_start; i <= j_end; i++) {
	    if (best_discrep[k2][0][k][mm][i] < q2) {
	      best_discrep[k2][1][k][mm][i] = best_discrep[k2][0][k][mm][i];
	      best_discrep[k2][0][k][mm][i] = q2;
	    }
	    else if (best_discrep[k2][1][k][mm][i] < q2) {
	      best_discrep[k2][1][k][mm][i] = q2;
	    }
	  }
	}

	confirm_flag = site1 > last_disc + 5; /* parameters->confirm_length; */
	if (w_site1 > w_last_disc + 3) { 
	  start = last_disc + parameters->confirm_trim + 1; 
	  start2 = last_disc2 + parameters->confirm_trim + 1; 
	  end = site1 - parameters->confirm_trim - 1;
	  end2 = start2 + (end - start);
	  w_start = w_last_disc + 1;
	  w_start2 = w_last_disc2 + 1;
	  
/*
	  start = (diff - 1)->site1 + parameters->confirm_trim + 1;
	  if (last_disc + 1 > start) start = last_disc + 1;
	  end = diff->site1 - parameters->confirm_trim - 1;
	  if (site1 - 1 < end) end = site1 - 1;
*/
/*
	  if (start < start0) start0 = start; 
	  end0 = end;
*/
/* a second reverse confirming read for a base counts as a forward read (to
   symmetrize the 2 + 1 rule) */
/*
	  if (confirm_flag) {
*/
	  for (i = w_start, j = w_start2; i <= w_end; i++, j++) {
	    k = i < start || i > end;
	    q2 = oqual2[j];
	    k2 = j < q_start2 || j > q_end2;
	    if (q2 < 98) {
	      if (q2 < 15) q2 = 15; /* have different (smaller) reward 
				       if !confirm_flag ? */
	      if (best_match[k2][0][k][i] < q2) {
		best_match[k2][1][k][i] = best_match[k2][0][k][i];
		best_match[k2][0][k][i] = q2;
	      }
	      else if (best_match[k2][1][k][i] < q2) 
		best_match[k2][1][k][i] = q2;
	    }
	    if (k) n_w_matches[i] += 1;
	    else n_matches[i] += 1;
	  }
	  if (create_test_file) {
	    if (str_confirm && (create_test_file == 1 || create_test_file == 3)) {
	      for (i = start, j = start2; i <= end; i++, j++) {
		q2 = oqual2[j];
		if (q2 >= 98) continue;
		q1 = oqual1[i];
		if (create_test_file == 3 && q1 && q2) continue;
		if (create_test_file == 1 && (!q1 || !q2)) continue;
		if (!q1) q1 = 1;
		if (!q2) q2 = 1;
		test_qual[i] = q1 + q2 < 90 ? q1 + q2 : 90;
	      }
	    }
	    if (!str_confirm && (create_test_file == 2 || create_test_file == 4)) {
	      for (i = start, j = start2; i <= end; i++, j++) {
		q2 = oqual2[j];
		if (q2 >= 98) continue;
		q1 = oqual1[i];
		if (create_test_file == 4 && q1 && q2) continue;
		if (create_test_file == 2 && (!q1 || !q2)) continue;
		if (!q1) q1 = 1;
		if (!q2) q2 = 1;
		test_qual[i] = q1 < q2 ? q2 : q1;
	      }
	    }
	    /*
	       }
	       else {
		 for (i = w_start; i <= w_end; i++) n_w_matches[i] += 1;
	       }
	       */
	    if (str_confirm && create_test_file == 5) {
	      for (i = w_start; i <= w_end; i++) {
		if (!confirm_flag || i < start || i > end) {
		  q1 = oqual1[i];
		  test_qual[i] = q1 + 15;
		}
	      }
	    }
	    if (!str_confirm && create_test_file == 6) {
	      for (i = w_start; i <= w_end; i++) {
		if (!confirm_flag || i < start || i > end) {
		  q1 = oqual1[i];
		  if (!q1) q1 = 1;
		  test_qual[i] = q1;
		}
	      }
	    }
	  } /* if (create_test_file) */
	} /* if (w_site1 > w_last_disc + 3) */
	last_disc = last_mononuc(seq + diff_site) - seq;
	last_disc2 = diff_site2 + (last_disc - diff_site);
	w_last_disc = next_w_last_disc; /* weak_last_mononuc(seq + diff_site) - seq;  diff_site; */
	w_last_disc2 = diff_site2 + (w_last_disc - diff_site);
	w_start = next_w_start;
      }
/* note that confirm_score not being used as criterion in the above */
/*
      if (strong_score >= parameters->confirm_score && !is_reject_chimeric(pair)) {
*/
      if (!is_reject_chimeric(pair)) {
	if (str_confirm && start0 < rev_first_start) rev_first_start = start0;
	if (str_confirm && end0 > rev_last_end) rev_last_end = end0;
	if (start0 < first_start) first_start = start0;
	if (end0 > last_end) last_end = end0;
      }
      if (create_test_file) {
	for (j = 0; j < length1 && !test_qual[j]; j++);
	if (j < length1) {
	  fprintf(fp,">%s", get_id(entry1));
	  for (j = 0; j < length1; j++) {
	    if (!(j % 50)) fprintf(fp,"\n");
	    fprintf(fp,"%c", seq[j]);
	  }
	  fprintf(fp,"\n");
	  fprintf(fq, ">%s", get_id(entry1));
	  for (j = 0; j < length1; j++) {
	    if (!(j % 50)) fprintf(fq,"\n");
	    fprintf(fq, "%d ", test_qual[j]);
	  }
	  fprintf(fq,"\n");
	}
      }
    }

    adj_qual = get_adj_qual(entry1);
    for (i = 0; i <= align_entry->last_vec; i++) adj_qual[i] = 0;
    if (pass_level < 2) {
      if (align_entry->qual_start < length1 - 1)
	for (; i < align_entry->qual_start; i++) 
	  if (oqual1[i] != 99) adj_qual[i] = 0;
      if (segment = align_entry->segments) {
	seg_start = segment->start;
	for (; segment; segment = segment->next)
	  if (!segment->next) seg_end = segment->end;
	for (; i < seg_start; i++) 
	  if (oqual1[i] != 99) adj_qual[i] = 0;
      }
    }

    for (; i < length1; i++) {
      q1 = oqual1[i];
      if (q1 == 99) adj_qual[i] = 99;
      else if (seq[i] == 'N' || seq[i] == 'X' || q1 == 98) adj_qual[i] = 0;
      else {
	q_flag = 1;

	max_x = max_y = 0;
	max_x2 = max_y2 = 0;
	for (k = 0; k < 1 /* 2 */; k++) {
	  for (m = 0; m < 2; m++) {
	    for (n = 0; n < 2; n++) {
	      x = best_same_match[k][m][n][i];
	      y = best_opp_match[k][m][n][i];

	      if (y < 15) y = 0;

	      if (k || m) {
		x = x < 15 ? 0 : 10;
		if (y > 10) y = 10;
	      }
	      else if (n) {
		x = x < 15 ? 0 : 15;
		if (y > 15) y = 15;
	      }
	      if (max_x < x && !k) {
		max_x2 = max_x;
		max_x = x;
	      }
	      else if (max_x2 < x) {
		max_x2 = x;
		if (k) q_flag = 0;
	      }

	      if (max_y < y && !k) {
		max_y2 = max_y;
		max_y = y;
	      }
	      else if (max_y2 < y && (q_flag || !k)) {
		max_y2 = y;
	      }
	    }
	  }
	}
/*
	q1 = best_same_match[0][0][0][i];
	if (q1 < best_same_match[0][0][1][i])
	  q1 = best_same_match[0][0][1][i];

	x = best_same_match[0][1][0][i];
	if (x > 10) x = 10;
	else if (x < best_same_match[1][0][0][i]) {
	  x = best_same_match[1][0][0][i];
	  q_flag = 0;
	}
	q1 += x;

	y = best_opp_match[0][0][0][i];
	if (y < best_opp_match[0][0][1][i])
	  y = best_opp_match[0][0][1][i];

	x = best_opp_match[0][1][0][i];
	if (x > 10) x = 10;
	else if (x < best_opp_match[1][0][0][i] && q_flag) {
	  x = best_opp_match[1][0][0][i];
	  q_flag = 0;
	}

	y += x;
	q1 += y;
*/
/*
	fprintf(stderr, "%d ", max_x2);
*/
	q1 = max_x + max_y;

	if (chemistry != 1) { /* don't further increment dye terminator reads */
	  max_x2 += max_y2;
	  if (max_x2 > 10) max_x2 = 10;
	  q1 += max_x2;
	}

/* Following helps low end distribution -- but hard to justify theoretically
	if (i >= first_start && i <= last_end && q1 < 15) q1 = 15;
*/
	max_x = 0;
	n_m = n_opp_matches[i] + n_w_opp_matches[i] + n_same_matches[i] + n_w_same_matches[i];

	if (!n_m && pass_level < 2) { /* omit for pass_level == 0 ? */
	  for (k = 0; k < 3; k++) {
	    for (mm = 0; mm < 4; mm++) {
	      x = best_same_discrep[0][0][k][mm][i] + best_opp_discrep[0][0][k][mm][i]; 
	      x = x >= 15 ? x : 0;
	      y = best_same_discrep[0][1][k][mm][i] + best_opp_discrep[0][1][k][mm][i]; 
	      y = y >= 15 ? 10 : 0;

	      x += y;
	      if (max_x < x) max_x = x;
	    }
	  }
	  if (max_x > 15) max_x = 15;
	}
	x = max_x;

/* following discounts bases that are not confirmed by an opposite strand
   or opposite chemistry read -- provided they are in the range that is
   confirmed by some such read. In particular this will minimize effects
   on LLR-score of compressions in dye-primer reads, or G-dropouts in
   (old) dye-terminator reads. However it will also tend to favor incorporation
   of "contaminants" (including reads from regions of assembly where there
   is not yet a strongly confirming read). This could affect
   consensus, or conceivably assembly accuracy.
      In pass 0, discounting is complete, so that all matches utilized
   in pass 1. In pass 1 discounting is partial (but strong).
      In final pass (prior to construction of consensus) no
   discounting is performed), to avoid affecting consensus generation.

   Modified version: specific to chemistry. If (old) dye terminator chemistry,
     and there is no confirming opposite strand or opposite chemistry read at
     that site, and there is a discrepant opposite strand or opposite 
     chemistry read which substitutes a 'G' at that position, then the
     quality is discounted to 0 (in first two passes).
   Similarly if dye primer chemistry, and no confirming opp str or chem read,
     and there is a discrepant opp strand or chem read which inserts a 'G' or
     'C', then quality is discounted to 0 (in first two passes).
*/

	
	n_m = n_opp_matches[i] + n_w_opp_matches[i];

	if (!n_m) {
	  if (chemistry == 1) {
	    if (best_opp_discrep[0][0][2][2][i] >= 15
		|| pass_level < 2 && 
		(best_opp_discrep[0][0][2][0][i] >= 15
		|| best_opp_discrep[0][0][2][1][i] >= 15
		|| best_opp_discrep[0][0][2][3][i] >= 15))
	       x = q1;
	  }
	  else if (pass_level < 2 && chemistry == 0 
		   && (best_opp_discrep[0][0][1][1][i] >= 15
		       || best_opp_discrep[0][0][1][2][i] >= 15))
	    x = q1; 
	}

/* similar to above -- except discount any base (outside high-quality
   region) that is not confirmed by any other read
*/

	if (pass_level < 2) {
	  n_m += n_same_matches[i] + n_w_same_matches[i];

	  if (!n_m && 
	      (i < align_entry->qual_start || i > align_entry->qual_end)) {
	    if (!pass_level) x = q1;
	    else if (x < 20) x = 20;
	  }
	}

	/* CHECK VALUES ABOVE */

/* issues wrt old version: promotion of quality to 15 if 
   i >= first_start && i <= last_end or there is any match (n_m > 0);
*/
/*
	fprintf(stderr, "%d %d| ", q1, x);
*/
	q1 -= x;
	adj_qual[i] = q1 < 2 ? 2 : q1 > 90 ? 90 : q1;
      }
    }
    if (pass_level < 2 && align_entry->segments) 
      for (i = seg_end + 1; i < length1; i++) 
	if (adj_qual[i] != 99) adj_qual[i] = 0;
 
    if (chemistry == 3)
/*AVOID ALL OF ABOVE IN THIS CASE!!!! */
      for (i = 0; i < length1; i++) adj_qual[i] = oqual1[i];
/*
    fprintf(stderr, "\n%d %d", align_entry->qual_start, align_entry->qual_end);

    fprintf(stderr,"\n%s %d %d\n", get_id(entry1), pass_level, chemistry);
    for (i = 0; i < length1; i++) fprintf(stderr, "%d ",adj_qual[i]);
*/
/*
    for (segment = align_entry->segments, i_bit = 1; segment; 
       segment = segment->next, i_bit *= 2) 
      if (i_bit & align_entry->chimera_bits) 
	for (i = segment->start; i <= segment->end; i++)
	  adj_qual[i] = 0;
*/

/* RESTORE?
*/
    align_entry->first_start = first_start;
    align_entry->last_end = last_end;

    align_entry->rev_first_start = rev_first_start;
    align_entry->rev_last_end = rev_last_end;

/*
    if (left_first_start < weak_first_start - parameters->maxgap)
      printf("\n%s   %d-%d    confirmed only by left-incomplete matches", 
	     align_entry->db_entry->id, left_first_start + 1, weak_first_start);
    if (left_last_end > weak_last_end + parameters->maxgap)
      printf("\n%s   %d-%d    confirmed only by left-incomplete matches", 
	     align_entry->db_entry->id, weak_last_end + 2, left_last_end + 1);
    if (right_first_start < weak_first_start - parameters->maxgap)
      printf("\n%s   %d-%d    confirmed only by right-incomplete matches", 
	     align_entry->db_entry->id, right_first_start + 1, weak_first_start);
    if (right_last_end > weak_last_end + parameters->maxgap)
      printf("\n%s   %d-%d    confirmed only by right-incomplete matches", 
	     align_entry->db_entry->id, weak_last_end + 2, right_last_end + 1);
*/
  }
  free_seg_blocks();
  notify(" Done\n");
  if (create_test_file) {
    fclose(fp);
    fclose(fq);
    exit(1);
  }
}
 
/* return index of 1st occurrence of 3d nuc following input position */
unsigned char *last_mononuc(seq1)
     unsigned char *seq1;
{
  unsigned char nuc1, nuc2;

  seq1++;
  nuc1 = *seq1;
  for (seq1++; *seq1 == nuc1; seq1++);
  nuc2 = *seq1;
  if (nuc2) for (seq1++; *seq1 == nuc1 || *seq1 == nuc2; seq1++);
  return seq1;
}

/* return index of 1st occurrence of new nuc following input position */
unsigned char *weak_last_mononuc(seq1)
     unsigned char *seq1;
{
  unsigned char nuc1, nuc2;

  seq1++;
  nuc1 = *seq1;
  for (seq1++; *seq1 == nuc1; seq1++);
  return seq1;
}

/*
int last_mononuc(diff, seq1, seq2)
     Diff *diff;
     unsigned char *seq1, *seq2;
{
  int last_disc, j, k;
  unsigned char nuc1, nuc2;

  nuc1 = nuc2 = 0;
  last_disc = diff->site1;
  if (diff->type == 'D') { 
    nuc1 = seq1[last_disc];
    for (j = last_disc + 1; seq1[j] == nuc1; j++);
    last_disc = j - 1;
  }
  else if (diff->type == 'I') {
    nuc1 = seq2[diff->site2];
    for (j = last_disc + 1; seq1[j] == nuc1; j++);
    last_disc = j - 1;
  }
  else if (diff->type == 'B' || diff->type == 'S') {
    for (j = last_disc + 1; seq1[j] == seq2[diff->site2]; j++);
    for (k = last_disc + 1; seq1[k] == seq1[last_disc]; k++);
    last_disc = j > k ? j - 1 : k - 1;
  }
  return last_disc;
}
*/
/* return index of 1st occurrence of 3d nuc preceding input position  */
unsigned char *first_mononuc(seq1)
     unsigned char *seq1;
{
  unsigned char nuc1, nuc2;

  seq1--;
  nuc1 = *seq1;
  for (seq1--; *seq1 == nuc1; seq1--);
  nuc2 = *seq1;
  if (nuc2) for (seq1--; *seq1 == nuc1 || *seq1 == nuc2; seq1--);
  return seq1;
}	    

/* return index of 1st occurrence of new nuc preceding input position  */
unsigned char *weak_first_mononuc(seq1)
     unsigned char *seq1;
{
  unsigned char nuc1, nuc2;

  seq1--;
  nuc1 = *seq1;
  for (seq1--; *seq1 == nuc1; seq1--);
  return seq1;
}	    

/*
int first_mononuc(diff, seq1, seq2)
     Diff *diff;
     unsigned char *seq1, *seq2;
{
  int site1, j, k;

  site1 = diff->site1;
  if (diff->type == 'D') {
    for (j = site1 - 1; seq1[j] == seq1[site1]; j--);
    site1 = j + 1;
  }
  else if (diff->type == 'I') {
    for (j = site1; seq1[j] == seq2[diff->site2]; j--);
    site1 = j + 1;
  }
  else if (diff->type == 'E' || diff->type == 'S') {
    for (j = site1 - 1; seq1[j] == seq2[diff->site2]; j--);
    for (k = site1 - 1; seq1[k] == seq1[site1]; k--);
    site1 = j < k ? j + 1 : k + 1;
  }
  return site1;
}	    
*/

/* recheck extents against overlapping reads in contig -- to see if any have changed */  
reset_extents(db)
     Database *db;
{
  Align_info *align_entry;
  Align_info *get_align_entry();
  int entry1, temp_start, temp_end;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  
  printf("\n\nRevised extents   start (old)    end  (old) ");
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    temp_end = 0;
    temp_start = get_seq_length(entry1);
    align_entry = get_align_entry(entry1);
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (!pair->score) continue;
      if (!is_used(pair)) continue;
      /*
	align_entry2 = align_array + pair->entry2; 
	if (align_entry2->contig != align_entry->contig 
	|| is_reverse(pair) != (align_entry->reverse != align_entry2->reverse)) continue;
	if (align_entry2->start > align_entry->end
	|| align_entry2->end < align_entry->start) continue;
	*/
      if (pair->start1 < temp_start) temp_start = pair->start1;
      if (pair->end1 > temp_end) temp_end = pair->end1;
    }
    if (temp_start > align_entry->first_start
	|| temp_end < align_entry->last_end) {
      /*
	printf("\n %-15s    %5d   (%d)    %5d  (%d)",
	align_array[entry1].db_entry->id, temp_start, align_entry->first_start, temp_end, align_entry->last_end);
	*/
      align_entry->first_start = temp_start;
      align_entry->last_end = temp_end;
    }
  }
}

find_unaligned_pairs()
{
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int entry1;
  notify("Finding unaligned pairs ... ");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) 
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) 
      set_unaligned_flag(pair, outside_match(pair));
  
  notify(" Done\n");
}

extern Database *query_db; 

output_nonmatching_queries()
{
  int entry1, length, i, q;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  FILE *fopenWrite(), *fp, *fq;
  char file_name[500], *qual;
  char *get_id(), *get_descrip(), *get_orig_qual(), *id, *descrip;
  unsigned char *get_seq(), *seq;

  sprintf(file_name, "%s.nonmatching", parameters->query_files->name);
  if (parameters->compact_qual) {
    strcat(file_name, ".calf");
    fp = fopenWrite(file_name);
    fprintf(fp, "%c%c", 0, 0); /* Include an ASCII header?? */
  }
  else {
    fp = fopenWrite(file_name);
    strcat(file_name, ".qual");
    fq = fopenWrite(file_name);
  }

  for (entry1 = query_db->first_entry; entry1 <= query_db->last_entry; entry1++) {
    /*  for (entry1 = 0; entry1 < t_num_entries; entry1++) { incorrect!!! */
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (pair->score) break;
    }
    if (!pair) {
      id = get_id(entry1);
      descrip = get_descrip(entry1);
      seq = get_seq(entry1);
      if (parameters->compact_qual) {
	length = get_seq_length(entry1);
	if (strlen(id) || strlen(descrip)) {
	  fprintf(fp, "%c%c%s", 62, 0, id);
	  if (strlen(descrip)) fprintf(fp, " %s", descrip);
	  fprintf(fp, "%c%c", 0, 62);
	}
	for (i = 0; i <= length; i++) fputc(seq[i], fp); /* includes terminating 0 */
      }
      else {
	fprintf(fp, ">%s %s\n%s\n", id, descrip, seq);
	fprintf(fq, ">%s\n", id);
	qual = get_orig_qual(entry1);
	for (i = 0; i < length; i++) fprintf(fq, " %d", qual[i]);
	fprintf(fq, "\n");
      }
    }
  }
  fclose(fp);
  if (!parameters->compact_qual) 
    fclose(fq);
}

Aligned_pair **sort_pairs(sort_function, both)
     int (*sort_function)();
     int both; /* 0 if only one copy of each pair is to be used, 1 if both are */
{
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  Aligned_pair **pair_pointers;
  int entry1;
  int i_ptr, np;
  char *our_alloc();

  notify("Sorting pairs ...");
  i_ptr = 0;
  np = both ? 2 * num_pairs : num_pairs; /* this currently too large! */
  pair_pointers = (Aligned_pair **)our_alloc(np * sizeof(Aligned_pair *));
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (both || pair->entry2 > pair->entry1 
	  || pair->entry2 == pair->entry1 && pair->start1 < pair->start2) {
	/* N.B. This is inadequate for case pair->entry2 == pair->entry1 and is_reverse(pair) */
	pair_pointers[i_ptr++] = pair;
      }
/*
      if (parameters->weak_first && pair->entry2 == pair->entry1 && pair->w_start1 == pair->w_start2) {
	fprintf(stderr, "\nERROR: SAME START  -- %d, orientation %d\n", 
		pair->w_start1, is_reverse(pair));
	exit(1);
      }
*/
    }
  }
  if (i_ptr != np) {
    fprintf(stderr, "\nERROR: DISCREPANCY in # pairs -- %d, %d, %d\n", i_ptr, np, num_pairs);
    exit(1);
  }
  qsort(pair_pointers, np, sizeof(Aligned_pair *), sort_function);
  notify(" Done");
  return pair_pointers; /* REMEMBER TO FREE LATER */
}

/* THE FOLLOWING NEEDS WORK -- SINCE CHANGES TO first_start AND last_end */
print_coverage(db) 
     Database *db;
{
  Align_info *align_entry, *align_entry2;
  Align_info *get_align_entry();
  int i, j, n_pairs, length1, length2;
  int entry1;
  double t_length, t_uc_length;
  double t_confirmed_length, t_str_confirmed_length, n_confirmed_reads;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int histogram[100];
  int prev_start, prev_end, rel_start, rel_end, confirmed_bases;
  char *our_alloc();
  int *depths;
  int max_depth, depth;
  int h_size;
  char *get_adj_qual();
  char *adj_qual;

  notify("Printing coverage ...");
  h_size = 26;
  for (i = 0; i < h_size; i++) 
    histogram[i] = 0;
  t_length = t_uc_length = t_confirmed_length = t_str_confirmed_length = n_confirmed_reads = 0;
  confirmed_bases = 0;
  
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    align_entry = get_align_entry(entry1);
    if (is_anomalous(align_entry)) continue;
    adj_qual = get_adj_qual(entry1);
    prev_start = length1 = get_seq_length(entry1);
    prev_end = 0;
    depths = (int *)our_alloc((length1 + 1) * sizeof(int));
    for (i = 0; i <= length1; i++) depths[i] = 0;
    n_pairs = max_depth = 0;
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      if (pair->LLR_score <= 0) continue;
      if (!is_best(pair)) continue;
      if (same_subclone(pair->entry1, pair->entry2)) continue; /* ???? */
      align_entry2 = get_align_entry(pair->entry2);
      length2 = get_seq_length(pair->entry2) - 1;
      n_pairs++; 
 
      depths[pair->start1] += 1;
      depths[pair->end1 + 1] -= 1;
      
      if (pair->entry2 < entry1 
	  && align_entry2->first_start < align_entry2->last_end) { 
	/* find previously covered part */
	rel_start = pair->start1 - pair->start2;
	rel_end = pair->end1 - pair->end2;
	if (is_reverse(pair)) {
	  rel_start += length2 - align_entry2->last_end;
	  rel_end += length2 - align_entry2->first_start;
	}
	else {
	  rel_start += align_entry2->first_start;
	  rel_end += align_entry2->last_end;
	}
	if (rel_start < prev_start) prev_start = rel_start;
	if (rel_end > prev_end) prev_end = rel_end;
      }
    }
    if (n_pairs) {
      t_length += length1;
      max_depth = depth = 0;
      for (i = 0; i < length1; i++) {
	depth += depths[i];
	if (depth > max_depth) max_depth = depth;
      }
      n_confirmed_reads += 1;
      if (align_entry->first_start < align_entry->last_end)
	t_confirmed_length += 
	  align_entry->last_end - align_entry->first_start + 1;
      for (i = 0; !adj_qual[i] && i < length1; i++);
      for (j = length1 - 1; !adj_qual[j] && j > i; j--);
      t_uc_length += j - i + 1;
      if (align_entry->first_start < align_entry->last_end) {
	if (prev_start > prev_end) 
	  confirmed_bases += align_entry->last_end - align_entry->first_start + 1;
	else {
	  if (prev_start > align_entry->first_start) 
	    confirmed_bases += prev_start - align_entry->first_start;
	  if (prev_end < align_entry->last_end) 
	    confirmed_bases += align_entry->last_end - prev_end;
	}
      }
    }
    our_free(depths);
    histogram[max_depth >= h_size ? h_size - 1 : max_depth] += 1;
    if (align_entry->rev_last_end > align_entry->rev_first_start)
      t_str_confirmed_length += align_entry->rev_last_end - align_entry->rev_first_start;
  }
  
  printf("\n\nNo. confirmed reads: %.0f", n_confirmed_reads);
  if (n_confirmed_reads) 
    printf("\nAvg. length: %.1f, confirmed: %.1f, str. confirmed: %.1f, trimmed: %.1f", 
	   t_length / n_confirmed_reads, 
	   t_confirmed_length / n_confirmed_reads,
	   t_str_confirmed_length / n_confirmed_reads,
	   t_uc_length / n_confirmed_reads);
  if (!parameters->subject_files)
    printf("\nPreliminary clone size estimate: %d bp, depth of coverage: %.1f", 
	   confirmed_bases, confirmed_bases ? t_confirmed_length / confirmed_bases : 0.0);
  printf("\n\nDepth histogram (max_depth, #reads, cum #reads):\n");
  for (i = h_size - 1, max_depth = 0; i >= 0; i--)
    if (histogram[i]) {
      max_depth += histogram[i];
      printf("\n%2d   %3d     %3d", i, histogram[i], max_depth);
    }
  notify(" Done\n");
}

/*
slide_indels()
{
  Align_info *align_entry, *align_entry2;
  Aligned_pair *pair;
  Diff *diff;
  unsigned char *seq1, *seq2;
  int entry1;

  for (entry1 = first_file_only ? num_query_entries : 0; 
       entry1 < t_num_entries; entry1++) {
    align_entry = align_array + entry1;
    seq1 = align_entry->db_entry->seq;
    for (pair = align_entry->first_pair; pair; pair = pair->next) {
      seq2 = (align_array[pair->entry2].db_entry + 
				  (is_reverse(pair) ? num_query_entries : 0))->seq;
      
      for (diff = pair->diffs + 1; diff->type != 'E'; diff++) {
	if (diff->type == 'S') continue;
	if (diff->type == 'I') {
	  if (diff->site1 < (diff - 1)->site1) {
	     previous diff has already been adjusted -- so make this one compatible 
	    diff->site1 = (diff - 1)->site1;
	    diff->site2 = (diff - 1)->site2 + 1;
	  }
	  else {
	    while (seq1[diff->site1 + 1] == seq2[diff->site2]
		   && seq1[diff->site1 + 1] != 'N' && seq1[diff->site1 + 1] != 'X') {
	     condition needed to exclude possibility that the while condition holds
	       but there is another diff (which must be at an 'N', since otherwise
	       rearrangement would increase score) at diff->site1 + 1 --
	       dealing with this would necessitate rearranging order of diffs 
	      diff->site1 += 1;
	      diff->site2 += 1;
	    }
	  }
	}
	else if (diff->type == 'D') {
	  if (diff->site2 < (diff - 1)->site2) {
	     previous diff has already been adjusted -- so make this one compatible 
	    diff->site1 = (diff - 1)->site1 + 1;
	    diff->site2 = (diff - 1)->site2;
	  }
	  else {
	    while (seq2[diff->site2 + 1] == seq1[diff->site1]
		   && seq2[diff->site2 + 1] != 'N' && seq2[diff->site2 + 1] != 'X') {
	      diff->site1 += 1;
	      diff->site2 += 1;
	    }
	  }
	}
      }
    }
  }
}
*/
analyze_discreps(db)
     Database *db;
{
  int i, j, k, length2;
  int entry1;
  Aligned_pair *pair;
  Aligned_pair *get_aligned_pairs();
  int spacing_hist[2][100];
  int n_subs[2], n_dels[2], n_inserts[2];
  double n_bases[2];
  char alphabet[15];
  unsigned char *diff;
  unsigned char *get_seq(), *get_comp_seq();
  unsigned char *seq1, *seq2;
  char *get_adj_qual();
  char *adj_qual1, *adj_qual2;
  int total, last_s, last_d, last_i, spacing, h_size;
  char c1, c2;
  /* signed */ char q1, q2;
  int nuc_counts[2][20][20];
  int qual_counts[2][20][20]; /* needs to be adjusted for higher quality values! */
  int site1, site2, diff_t;
  int counts[20];

  strcpy(alphabet, "ACGTNXZ");
  h_size = 26;
  for (i = 0; i < h_size; i++) 
    spacing_hist[0][i] = spacing_hist[1][i] = 0;
  
  for (k = 0; k < 2; k++) 
    for (i = 0; i < 20; i++) 
      for (j = 0; j < 20; j++) 
	nuc_counts[k][i][j] = qual_counts[k][i][j] = 0;
  
  for (i = 0; i < 2; i++)
    n_subs[i] = n_dels[i] = n_inserts[i] = n_bases[i] = 0;
  
  for (entry1 = db->first_entry; entry1 <= db->last_entry; entry1++) {
    seq1 = get_seq(entry1); 
    adj_qual1 = get_adj_qual(entry1);
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
/*      if (!is_best(pair)) continue; */
      if (pair->LLR_score <= 0) continue;
      if (same_subclone(pair->entry1, pair->entry2)) continue; /* ???? */
      length2 = get_seq_length(pair->entry2) - 1;
      last_s = last_d = last_i = -1;
      seq2 = is_reverse(pair) ? get_comp_seq(pair->entry2) : get_seq(pair->entry2);
      adj_qual2 = get_adj_qual(pair->entry2);
      site1 = pair->start1 - 1;
      site2 = pair->start2 - 1;
      for (diff = pair->diffs; *(diff + 1); diff++) {
	site1 += diff_gap1(*diff);
	site2 += diff_gap2(*diff);
	if ((diff_t = diff_type(*diff)) == 'M') continue;
	if (diff_t == 'S') {
	  c1 = seq1[site1];
	  c2 = seq2[site2];
	  q1 = adj_qual1[site1] / 10;
	  q2 = adj_qual2[is_reverse(pair) ? 
				  length2 - site2 : site2] / 10;
	  qual_counts[is_reverse(pair)][q1][q2] += 1;
	  /*
	    if (q1 >= 2 && q2 >= 2) {
	    fprintf(stderr, "\n\n");
	    for (i = site1 - 5; i <= site1 + 5; i++)
	    fprintf(stderr, "%c", seq1[i]);
	    fprintf(stderr, "\n");
	    for (i = site2 - 5; i <= site2 + 5; i++)
	    fprintf(stderr, "%c", seq2[i]);
	    }
	    */
	  n_subs[is_reverse(pair)] += 1;
	  i = strchr(alphabet, c1) ? strchr(alphabet, c1) - alphabet : strlen(alphabet) - 1;
	  j = strchr(alphabet, c2) ? strchr(alphabet, c2) - alphabet : strlen(alphabet) - 1; 	
	  nuc_counts[is_reverse(pair)][i][j] += 1;
	}
	else if (diff_t == 'D') {
	  n_dels[is_reverse(pair)] += 1;
	  if (last_i > -1) {
	    spacing = site1 - last_i;
	    if (spacing >= h_size) spacing = h_size - 1;
	    spacing_hist[is_reverse(pair)][spacing] += 1;
	    last_i = -1;
	  }
	}
	else if (diff_t == 'I') {
	  n_inserts[is_reverse(pair)] += 1;
	  if (last_d > -1) {
	    spacing = site1 - last_d + 1;
	    if (spacing >= h_size) spacing = h_size - 1;
	    spacing_hist[is_reverse(pair)][spacing] += 1;
	    last_d = -1;
	  }
	}
	if (diff_t == 'I') last_i = site1;
	else if (diff_t == 'D') last_d = site1;
      }
      n_bases[is_reverse(pair)] += pair->end1 - pair->start1 + 1;
    }
  }

  for (j = 0; j < 20; j++) {
    counts[j] = 0;
    for (i = 0; i < 2; i++)
      for (k = 0; k < 20; k++)
	counts[j] += qual_counts[i][k][j] + qual_counts[i][j][k];
  }
  
  for (i = 0; i < 2; i++) {
    printf("\n\n%s confirmed bases: %.0f", i ? "Reverse" : "Forward", n_bases[i]);
    if (n_bases[i])
      printf("\n Subs: %d (%.2f%%), dels: %d (%.2f%%), inserts: %d (%.2f%%)", 
	     n_subs[i], 100.0 * n_subs[i] / n_bases[i], 
	     n_dels[i], 100.0 * n_dels[i] / n_bases[i], 
	     n_inserts[i], 100.0 * n_inserts[i] / n_bases[i]);
    printf("\n\nSubstitutions by nucleotide:\n ");
    for (j = 0; j < strlen(alphabet); j++) 
      printf("      %c", alphabet[j]);
    printf("    Total");
    for (k = 0; k < strlen(alphabet); k++) {
      total = 0;
      printf("\n%c", alphabet[k]);
      for (j = 0; j < strlen(alphabet); j++) {
	total += nuc_counts[i][k][j];
	printf(" %6d", nuc_counts[i][k][j]);
      }
      printf("   %6d", total);
    }
    printf("\n\nSubstitutions by quality: \n   ");
    for (j = 0; j < 20; j++) 
      if (counts[j]) printf(" %6d", j); /* ASSUMES MAX QUAL <= 19 */
    printf("    Total");
    for (k = 0; k < 20; k++) {
      if (!counts[k]) continue;
      total = 0;
      printf("\n%d  ", k);
      for (j = 0; j < 20; j++) {
	if (counts[j]) {
	  total += qual_counts[i][k][j];
	  printf(" %6d", qual_counts[i][k][j]);
	}
      }
      printf("   %6d", total);
    }
    
    printf("\n\nHistogram of spacings between adjacent indel pairs:\n");
    for (j = 0; j < h_size; j++)
      if (spacing_hist[i][j]) {
	printf("\n%s%2d   %3d   ", 
	       j == h_size - 1 ? ">=" : "  ", j, spacing_hist[i][j]);
      }
  }
}

get_discreps(diffs, mismatches, insertions, deletions)
     unsigned char *diffs;
     int *mismatches, *insertions, *deletions;
{
  char type, p_type;
  int n;

  *mismatches = *insertions = *deletions = 0;
  for (; *(diffs + 1); diffs++) {
    type = diff_type(*diffs);
    if (type == 'S') *mismatches += 1;
    else if (type == 'D') *deletions += 1;
    else if (type == 'I') *insertions += 1;
    else if (type == 'm') {
      n = parameters->word_intron_margin;
      p_type = diff_type(*(diffs - 1));
      if (p_type == 'D') *deletions -= n;
      else if (p_type == 'I') *insertions -= n;
      for (diffs++; diff_type(*diffs) == 'M'; diffs++);
      diffs += n;
    }
  }
}

/* following no longer used; is a previous attempt to define quality

find_depths(contig)
     Contig *contig;
{
  char *our_alloc();
  int i, entry1, length1, c_depth;
  int *depth;
  Aligned_pair *pair;
  Align_info *align_entry;
  Diff *diff;

  for (align_entry = contig->first; align_entry; align_entry = align_entry->next) {
    entry1 = align_entry - align_array;
    length1 = align_entry->db_entry->length;
    depth = align_entry->depth = (int *)our_alloc((length1 + 1) * sizeof(int));
    for (i = 0; i < length1; i++) depth[i] = 0; 
    for (pair = align_entry->first_pair; pair; pair = pair->next) {
      if (!pair->used) continue;
      if (!pair->score) continue;
      if (align_entry->reverse) {
	depth[length1 - 1 - pair->end1] += 1;
	depth[length1 - pair->start1] -= 1;
      }
      else {
	depth[pair->start1] += 1;
	depth[pair->end1 + 1] -= 1;
      }
    }
    for (i = 0, c_depth = 1; i < length1; i++) { 
      c_depth += depth[i];
      depth[i] = c_depth;
    }
    for (pair = align_entry->first_pair; pair; pair = pair->next) {
      if (!pair->used) continue;
      for (diff = pair->diffs + 1; diff->type != 'E'; diff++) {
	if (align_entry->reverse) 
	  depth[length1 - 1 - diff->site1] -= 1;
	else depth[diff->site1] -= 1;
      }
    }
    for (i = 0, c_depth = 0; i < length1; i++) { 
      c_depth += depth[i];
    }
    if (c_depth/length1 > 20) {
      fprintf(stderr, "\n%s Avg depth: %.1f", 
	      align_array[entry1].db_entry->id, c_depth/(float)length1);
    }
    
  }
}
    
*/

outside_match(pair)
     Aligned_pair *pair;     
{
  Align_info *align_entry;
  Align_info *get_align_entry();
  int excess;

  align_entry = get_align_entry(pair->entry1);
  excess =  
    (align_entry->reverse ? get_seq_length(pair->entry1) - align_entry->m_end : 
     align_entry->m_start - 1) - pair->start1;
  if (excess > .5 * (pair->end1 - pair->start1 + 1)) return 1;
  excess =  pair->end1 -
    (align_entry->reverse ? get_seq_length(pair->entry1) - align_entry->m_start : 
     align_entry->m_end - 1);
  if (excess > .5 * (pair->end1 - pair->start1 + 1)) return 1;
  return 0;
}

/* fuse pairs with overlapping alignment: only relevant with gap1
   (as opposed to banded) search
*/
compare_entry2(pair_1, pair_2)
     Aligned_pair **pair_1, **pair_2;
{
  Aligned_pair *pair1, *pair2;
  int d;
  
  pair1 = *pair_1;
  pair2 = *pair_2;
  if (d = pair1->entry2 - pair2->entry2) return d;
  if (d = is_reverse(pair1) - is_reverse(pair2)) return d;
  if (d = pair1->start1 - pair2->start1) return d;
  return pair1->end1 - pair2->end1;
}

fuse_all_pairs()
{
  int entry1, i, n_p_1, i_p_1, j_p_1;
  Aligned_pair *get_aligned_pairs(), *pair, *pair1, *pair2, **pair_ptrs;
  char *our_alloc();

  notify("fusing pairs..."); 
  for (entry1 = 0; entry1 < t_num_entries; entry1++) { /* note also considering entries without any pairs!! */
    if (!(pair = get_aligned_pairs(entry1))) continue;
    for (n_p_1 = 0; pair; pair = pair->next) n_p_1++;
    pair_ptrs = (Aligned_pair **)our_alloc(n_p_1 * sizeof(Aligned_pair *));
    for (pair = get_aligned_pairs(entry1), i_p_1 = 0; pair; pair = pair->next) pair_ptrs[i_p_1++] = pair;
    qsort(pair_ptrs, n_p_1, sizeof(Aligned_pair *), compare_entry2);
    for (i_p_1 = 0; i_p_1 < n_p_1; i_p_1++) {
      pair1 = pair_ptrs[i_p_1];
      for (j_p_1 = i_p_1 + 1; j_p_1 < n_p_1; j_p_1++) {
	pair2 = pair_ptrs[j_p_1];
	if (pair2->entry2 < pair1->entry2) fatalError("pair sort");
	if (pair2->entry2 != pair1->entry2 || is_reverse(pair1) != is_reverse(pair2)) break;
	if (pair2->start1 < pair1->start1) fatalError("pair sort");
	fuse_pairs(pair1, pair2);
      }
    }
    our_free(pair_ptrs);
  }

  notify("done\n");
}

fuse_pairs(pair1, pair2)
     Aligned_pair *pair1, *pair2;
{
  unsigned char *diff1, *diff2, *diff, *diffs, *set_diff_block();
  int mid, pos1, pos2, gap, gap1, gap2, orig_pos2, n_diffs, i_diff, large_flag, length, score;
  Query_domain *query_domain;

  if (
      !pair1->score || !pair2->score /* already merged */
      /* || pair1->entry1 != pair2->entry1 -- tested outside */
      /* || pair1->entry2 != pair2->entry2  -- tested outside */
      /* || is_reverse(pair1) != is_reverse(pair2) */
      /* || pair1->start1 >= pair2->start1 */ || pair1->end1 < pair2->start1 /* || pair1->end1 >= pair2->end1 */ /* improper overlap */
      || pair1->start2 > pair2->end2 || pair1->end2 < pair2->start2
      || is_unaligned(pair1) && is_unaligned(pair2) && (pair1->spl3 & 64) != (pair2->spl3 & 64) /* spliced on opp strands */
      /* || abs((pair1->end1 - pair1->end2) - (pair2->start1 - pair2->start2)) > 2 /* offsets different */
      ) return 0;

  if (pair1->start1 == pair2->start1 && pair1->start2 == pair2->start2
      || pair1->end1 == pair2->end1 && pair1->end2 == pair2->end2) { /* same positions aligned; choose higher scoring
									or smaller; maybe incorrect?
								     */
    if (pair2->score < pair1->score 
	|| pair2->score == pair1->score && pair2->end1 - pair2->start1 >= pair1->end1 - pair1->start1) {
      pair2->score = 0;
    }
    else pair1->score = 0;
    return 1;
  }

  mid = (pair2->start1 + pair1->end1) / 2;

  for (diff1 = pair1->diffs; *diff1; diff1++); /* avoid by inputting length?? */

  /* fprintf(stderr," %d %d %d %d %d %d %d ", mid, pair1->start1, pair1->end1, pair1->start2, pair1->end2, pair2->start2, pair2->end2); */
  for (diff1--, pos1 = pair1->end1 + 1, pos2 = pair1->end2 + 1; diff1 >= pair1->diffs; diff1--){
    /* fprintf(stderr," %c %d %d %d ", diff_type(*diff1), diff_gap1(*diff1), diff_gap2(*diff1), pos2); */
    if (diff_type(*diff1) == 'm') {
      notify("mdiff ");
      /*
      fprintf(stderr," %d %d %d %d %d %d %d ", mid, pair1->start1, pair1->end1, pair1->start2, pair1->end2, pair2->start2, pair2->end2); 
      */
      return 0;
    }
    pos1 -= diff_gap1(*diff1); /* pos at next diff */
    pos2 -= diff_gap2(*diff1);
    if (pos1 == mid) { /* mid is at discrepancy; change it, until get one that isn't */
      mid--;
    }
    else if (pos1 < mid) {
      gap1 = mid - pos1;
      pos2 += gap1;
      pos1 = mid;
      break;
    }
  }
  if (pos1 != mid) {
    notify("WARNING: midpoint not found");
    return 0;
  }
  orig_pos2 = pos2;

  for (diff2 = pair2->diffs, pos1 = pair2->start1 - 1, pos2 = pair2->start2 - 1; *diff2; diff2++) {
    if (diff_type(*diff2) == 'm') {
      notify("m diff");
      return 0;
    }
    gap = diff_gap1(*diff2);
    if (pos1 + gap > mid) { /* possibly not correct ? */
      gap2 = pos1 + gap - mid;
      pos2 += mid - pos1;
      pos1 = mid;
      break;
    }
    else {
      pos1 += gap; /* pos at current diff */
      pos2 += diff_gap2(*diff2);
    }
  }
  if (pos1 != mid) {
    notify("WARNING: midpoint not found");
    return 0;
  }
  if (pos2 != orig_pos2) {
    fprintf(stderr, "MIDPOINT DIFF: %d %d", pos2, orig_pos2);
    return 0;
  }
  /* diffs: from pair1->diffs to diff1 - 1 (possibly empty) 
            'M' at mid
            from diff2 to end
  */
  /*
  for (diff = pair2->diffs, pos1 = pair2->start1 - 1, pos2 = pair2->start2 - 1, large_flag = 0; *diff; diff++) {
        if (diff_type(*diff) == 'm') 
      large_flag = !large_flag;
    pos1 += diff_gap1(*diff);
    if (!large_flag)
      pos2 += diff_gap2(*diff);
  }
  fprintf(stderr, " [%d,%d]", pos1 - (pair2->end1 + 1), pos2 - (pair2->end2 + 1));
  */

  for (diff = diff2; *diff; diff++);

  n_diffs = (diff1 - pair1->diffs) + (diff - diff2); /* don't count last diff from pair2 */

  diffs = set_diff_block(n_diffs);

  for (diff = pair1->diffs, i_diff = 0; diff < diff1; diff++)
    diffs[i_diff++] = *diff;

  diffs[i_diff++] = set_diff(gap1, 'M');
 
  *diff2 = set_diff(gap2, diff_type(*diff2));

  for (diff = diff2; *diff; diff++)
    diffs[i_diff++] = *diff;

  if (i_diff != n_diffs + 1) fatalError("diff count");

  score = pair1->score < pair2->score ? pair1->score : pair2->score;
  length = pair1->end1 - pair2->start1 + 1;
  if (length > score) length = score;

  pair1->score = pair1->score + pair2->score - length; /* APPROXIMATE!! & FAILS WITH NON-STANDARD SCORING */
  pair2->score = 0;


  pair1->end1 = pair2->end1;
  pair1->end2 = pair2->end2;
  
  if (parameters->score_flag) {
    for (query_domain = pair1->query_data->query_domain; query_domain != query_domain->parent; 
	 query_domain = query_domain->parent);
    /*    printf("***HERE2: %d %d", pair1->score, query_domain->best_score); */
    if (pair1->score > query_domain->best_score) { 
      query_domain->next_best = query_domain->best_score;
      query_domain->best_score = pair1->score;
      query_domain->n_best = 1;
      query_domain->best_pair = pair1;
      if (is_reverse(pair1)) {
	length = get_seq_length(pair1->entry2) - 1;
	query_domain->start = length - pair1->end2;
	query_domain->end = length - pair1->start2;
      }
      else {
	query_domain->start = pair1->start2;
	query_domain->end = pair1->end2;
      }
    }
    else if (pair1->score == query_domain->best_score) {
      query_domain->n_best += 1;

    }
  }

  if (is_unaligned(pair2)) {
    set_unaligned_flag(pair1, 1); /* flag as spliced */
    pair1->spl3 = pair2->spl3; /* would be better to retain this all! */
    pair1->spl5 = pair2->spl5;
  }
  /*
  fprintf(stderr, "*");
  for (diff = pair1->diffs; *diff; diff++)
    fprintf(stderr, " %c %d", diff_type(*diff), diff_gap1(*diff));
  fprintf(stderr, "*");
  for (diff = pair2->diffs; *diff; diff++)
    fprintf(stderr, " %c %d", diff_type(*diff), diff_gap1(*diff));
  fprintf(stderr, "*");
  for (diff = diffs; *diff; diff++)
    fprintf(stderr, " %c %d", diff_type(*diff), diff_gap1(*diff));
  */

  pair1->diffs = diffs;

  /*
  pair1->query_domain = 
  score_hist = 
  score_hist->best_pair
  notify("S");
  */
  return 1;
}

static int revise_minscore;

revise_scores(query_db)
     Database *query_db; 
{
  Aligned_pair *pair;
  int entry1;
  Aligned_pair *get_aligned_pairs();
  Query_domain *query_domain;
  Seq_entry *get_seq_entry();
  char *get_id();

  notify("Revising scores ...");
  revise_minscore = parameters->minscore;

  for (entry1 = 0; entry1 < t_num_entries; entry1++) 
    for (pair = get_aligned_pairs(entry1); pair; pair = pair->next) {
      pair->score = pair->score < revise_minscore ? 0 : pair->score - revise_minscore;
    }
  for (entry1 = query_db->first_entry; entry1 < query_db->first_entry + query_db->num_entries; entry1++) {
    query_domain = get_seq_entry(entry1)->query_domains;
    revise_qd_scores(query_domain);
  }

  parameters->min_record_score = parameters->minscore = 0;
  notify("done");
}

revise_qd_scores(query_domain)
     Query_domain *query_domain;
{
  Score_hist *s_h;

  if (!query_domain) return;
  query_domain->best_score -= revise_minscore;
  query_domain->next_best -= revise_minscore;
  if (query_domain->best_score < 0) query_domain->best_score = 0;
  if (query_domain->next_best < 0) query_domain->next_best = 0;
  for (s_h = query_domain->score_hist[0]; s_h; s_h = s_h->next) 
     s_h->score -= revise_minscore; 
  for (s_h = query_domain->score_hist[1]; s_h; s_h = s_h->next) 
     s_h->score -= revise_minscore; 

  revise_qd_scores(query_domain->child[0]);
  revise_qd_scores(query_domain->child[1]);
}
