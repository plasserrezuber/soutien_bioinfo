/*****************************************************************************
#   Copyright (C) 2006 by Phil Green.                          
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

int splice_edge_len;

/* 
  strategy for cDNA matching to genome :
    find matches -- both high-scoring (minscore) anchoring matches, and lower scoring 'marginal' ones (near_minscore) in reasonably close proximity
     (and same orientation; at present, no other filters are applied at this point, but some of the ones used below could
      be). The marginal matches are statistically
     significant (given location constraints) 
    find_left_sites(), find_right_sites(). find candidate sites near boundaries of matches, and scan for 'extensions' (matches to adjacent segments of reads);
       -- the latter need to be appended to the match list
    filter_matches() filter (but don't remove) matches: when there are multiple matches involving the same segment of the cdna
       identify informative sites (positions where some of the matches to the genome support the read base, and others don't)
       and tag matches that have excessive numbers of discrepancies (some subtleties in 'excessive')
    group_matches() need to allow chimeric reads. Identify parts of reads (seg_equiv classes) (and corresponding matches) that
       should be grouped together, because they correspond to same underlying gene. 
    analyze_multiple()  find consistent 'parses' -- sets of matches which are consistent with a gene match, simply on the
          basis of consistent order and orientation of matches between cdna and genome, and plausible intron size
          conditions on parse order of successive segments in read is same as that in genome; minimal extension length
            for read segments; maximal allowed overlap of read segments. implied intron size either close
            to zero (implying essentially contiguous in genome) or have intron length compatible with real gene.


          do parses both allowing, and disallowing, the filtered out matches.
          
        use dynamic programming to find best parses

        begin by eliminating marginal matches which don't fit into a parse.

        issues regarding parses:  
          problematic cases: 
           multiple parses of same seg_equiv class: 
              ambiguous links (to right)
	      > 1 start (ambiguous links to left)
	      multiple independent parses
           parts of cdna missing from all parses -- considering both included and deleted 'marginal' matches
           single cluster includes multiple seg_equiv classes (suggests parsing error)
           distinct seg_equiv classes include overlapping parts of reads
  
       check plausibility of splice sites; also see if this can be used to discriminate among multiple parses.
       
       evaluate splice junctions, identify ambiguous cases using probabilistic model.

  CURRENTLY SCANS CAN GO PAST END OF GENOME SEQUENCE -- NEED TO FIX THIS -- check for 0 character

improvements: 
  HANDLE MIXED READS!! subtract highest peak.
  revisit setting of reject_flags for matches -- may not be best way to handle
  control sequences for evaluating performance: random genomic segments + simulated errors.

  handle more than 2 ambiguous links

  hardwired constants

  maximize seg_equiv length over requirement that only one match per subcluster

  check necessity of cand_splice, splice_site structure members; reduce aligned_pair storage (eliminate
    stored alignment when possible).
  make default splice_site->offset clearly out of range.

  some functions (e.g. score_adjust) assume default scoring -- need to account for this.

  explicit test for duplicated read segments (all match same genomic segment, and only that one)
   identify primer occurrences.

   Most apparent 'missing exon' cases are low qual sequence where sequence is not adequately scanned
  penalty in bits now!! 
   improve memory usage: release sequence from pairs (currently needed only to print out alignment, once splice sites obtained
        (in find_left_sites and find_right_sites)  --
       not really needed); reduce storage per pair

 parameters to explore: bandwidth; max_intron_length; splice_edge_len; also minmatch; try using smaller minmatch
       but with option to extend it when frequency is high (won't work -- because frequency is on query strand!!)

    reduce memory requirements for pairs: put all splicing related stuff in single structure; reduce sequence storage (only end of alignment -- not flanking genomic).
    check for being near end of genomic sequence -- current analyze_splices ignores that.

  is reject_flags actually necessary in set_matches?

 have a -cdna_mode parameter.

 look for adjacent matches (as opposed to ones separated by intron) -- may be getting ignored.

 choose among ambig links by spliceability

 'slop' in choosing best parse should not be allowed to ignore short exons!!

 marking of accepted matches:also mark next best!

 increase penalty for large neg overhangs (or not!!)

 in chains, need to check for consistency of strand

 print outs: gene structures for each EST; when more than one is compatible, label separately with OR.
  combine data from multiple reads. Eliminate alternatives that require data errors (i.e. neg overhangs) in same locations, as unlikely (altho could be polymorphisms!).

 some parses are implausible -- when have duplicated exons; remove these!!

*/

typedef struct cand_splice {
  struct cand_splice *next;
  Splice_site *left_site, *right_site;
  int penalty, strand, missing_exon, left_align_overhang, right_align_overhang, left_cdna_overhang, right_cdna_overhang;
} Cand_splice;

typedef struct seg_equiv { /* segment equivalence classes */
  struct seg_equiv *parent;
  Aligned_pair *pair;
  Segment *segment;
  int start, end, length, id, n_chains, n_segments, best_length[2], best_start[2], start_count[2];
  char usable;
} Seg_equiv;

Cand_splice *head_cand_splice, *best_cand_splice, *next_best_cand_splice;
int *score_adjust, *cdna_indel_adjust, *genome_indel_adjust;
int min_penalty_gap;
int splice_penalty_hist[41], missing_exon_hist[56], length_discrep_hist[21];
char *exon_seq;
double base_llr, intron_coeff;

set_splice_params()
{
  int i, j;
  char *our_alloc();
  double q, mean_intron_size;

  min_penalty_gap = 6;
  mean_intron_size = 500;
  
  splice_edge_len = parameters->splice_edge_length;

  score_adjust = (int *)our_alloc(3 * splice_edge_len * sizeof(int));
  cdna_indel_adjust = score_adjust + splice_edge_len;
  genome_indel_adjust = score_adjust + 2 * splice_edge_len;

  for (i = 0; i < 41; i++) splice_penalty_hist[i] = 0;
  for (i = 0; i < 56; i++) missing_exon_hist[i] = 0;
  exon_seq = (char *)our_alloc(100 * sizeof(char));

  for (i = 0; i < 21; i++) length_discrep_hist[i] = 0;
  
  q = 1 / mean_intron_size;
  base_llr = log(q / .75) / log(2.0) + 3;
  intron_coeff = log(1 - q) / log(2.0);

  printf("\nBase_llr %.2f, intron_coeff: %.8f", base_llr, intron_coeff); 
}

/* want to penalize: 
   negative overhangs (reflect sequence errors) -- BUT should still be some detectable similarity,
   positive overhangs (long ones unlikely to occur by chance -- barring some biological phenomenon!!)
   missing exons: very short (+/-1) reflect indels; 
                  small ones rare;

   roughly speaking, each penalty point is 1 bit (so -- each extra overhang base is 2 points)		 

   examine gap characters -- to see if anything missed
   cases of large missing exons: want to make sure this part of read doesn't match somewhere else, instead
*/

append_cand_splice(strand, left_splice_site, right_splice_site, overlap, intron_length, left_pair, right_pair)
     int strand, overlap, intron_length;
     Splice_site *left_splice_site, *right_splice_site;
     Aligned_pair *left_pair, *right_pair;
{
  char *our_alloc();
  Cand_splice *cand_splice;
  int missing_exon, penalty;
  int left_cdna_overhang, right_cdna_overhang, n, check_right_overhang, check_left_overhang;
  int left_align_overhang, right_align_overhang;
  
  if (left_splice_site->type == 5 && right_splice_site->type != 2 || right_splice_site->type == 5 && left_splice_site->type != 2) return; 
  /* U12 constraint: AC must pair with ATATCC */
 
  penalty = left_splice_site->penalty + right_splice_site->penalty;  /* 0; */

  left_cdna_overhang = left_align_overhang = -left_splice_site->loc;
  if (left_align_overhang > 0) 
    left_cdna_overhang -= left_splice_site->cdna_overhang_adjust; 

  right_cdna_overhang = right_align_overhang = right_splice_site->loc;
  if (right_align_overhang > 0) 
    right_cdna_overhang -= right_splice_site->cdna_overhang_adjust; 

  missing_exon = right_cdna_overhang + left_cdna_overhang - overlap;

  check_right_overhang = check_left_overhang = 1;

  if (missing_exon < parameters->min_exon_length || intron_length / 2 < parameters->min_intron_length /* intron too small to fit!! */ ) 
    penalty += abs(missing_exon) * 5; /* 5 bits per inserted or deleted base */
  else penalty += 8; /* 1/256 prob. of missed exon -- too low? */

  if (missing_exon == 1 || missing_exon == -1) {
    n = (missing_exon + 1) / 2; /* 0 for deletion in read, 1 for insertion */
    check_right_overhang = right_splice_site->overhang_check[n];
    check_left_overhang = left_splice_site->overhang_check[n];
  }

  if (check_right_overhang) {
    if (right_align_overhang < 0) penalty += 5; /* error prob at particular location */
    if (right_align_overhang < -5) penalty += 5; /* truncated earlier than explainable by a single error in read */
  }

  if (check_left_overhang) {
    if (left_align_overhang < 0) penalty += 5; 
    if (left_align_overhang < -5) penalty += 5; 
  }

  cand_splice = (Cand_splice *)our_alloc(sizeof(Cand_splice));
  cand_splice->left_site = left_splice_site;
  cand_splice->right_site = right_splice_site;
  cand_splice->penalty = penalty;
  cand_splice->strand = strand;
  cand_splice->missing_exon = missing_exon;
  cand_splice->left_align_overhang = left_align_overhang; 
  cand_splice->right_align_overhang = right_align_overhang;
  cand_splice->left_cdna_overhang = left_cdna_overhang; 
  cand_splice->right_cdna_overhang = right_cdna_overhang;
  cand_splice->next = head_cand_splice;
  head_cand_splice = cand_splice;
  if (!best_cand_splice || penalty < best_cand_splice->penalty) {
    next_best_cand_splice = best_cand_splice;
    best_cand_splice = cand_splice;
  }
  else if (!next_best_cand_splice || penalty < next_best_cand_splice->penalty) {
    next_best_cand_splice = cand_splice;
  }
}

print_edges(edge_len, left_pair, right_pair)
     int edge_len;
     Aligned_pair *left_pair, *right_pair;
     
{
  int i;

  printf("\n");
  for (i = 6 * edge_len; i < 8 * edge_len; i++) {
    if (i == 7 * edge_len) printf(" ");
    printf("%c", left_pair->query_data->edges12[i]);
  }
  printf("..");
  for (i = 4 * edge_len; i < 6 * edge_len; i++) {
    if (i == 5 * edge_len) printf(" ");
    printf("%c", right_pair->query_data->edges12[i]);
  }
  printf("\n");
  for (i = 2 * edge_len; i < 4 * edge_len; i++) {
    if (i == 3 * edge_len) printf(" ");
    printf("%c", left_pair->query_data->edges12[i]);
  }
  printf("..");
  for (i = 0; i < 2 * edge_len; i++) {
    if (i == edge_len) printf(" ");
    printf("%c", right_pair->query_data->edges12[i]);
  }
}


print_cand_splice(cand_splice)
     Cand_splice *cand_splice;
{
  printf(" (%d: %d, %d, %d, %d)", cand_splice->penalty,
	 cand_splice->strand, cand_splice->missing_exon, 
	 cand_splice->left_align_overhang, cand_splice->right_align_overhang);
}     

free_cand_splices()
{
  Cand_splice *cand_splice;

  for (cand_splice = head_cand_splice; cand_splice; cand_splice = cand_splice->next)
    our_free(cand_splice);
  next_best_cand_splice = best_cand_splice = head_cand_splice = 0;
}

analyze_splices(prev_pair, pair, print_flag, incr_flag, orient, seq1, seq2)
     Aligned_pair *pair, *prev_pair;
     int print_flag, incr_flag, orient;
     char *seq1, *seq2;
{
  int i, j, k, m, n, c, d, cum, found, start, end, penalty, score, max_score;
  int left_cdna_overhang, right_cdna_overhang, missing_exon, strand, best_start, n_best, best_m, first_cand_splice;
  Aligned_pair *left_pair, *right_pair;
  Cand_splice *cand_splice;
  int overlap, intron_length, fail_flag, left_extend, right_extend;
  static int intron_num, prev_strand;
  Splice_site *left_site, *right_site, *best_left_site, *best_right_site;
  char *get_id();
  Splice_site *get_best_site();

  left_cdna_overhang = right_cdna_overhang = missing_exon = -1000;

  if (orient) {
    if (prev_pair) {
      overlap = prev_pair->end1 - pair->start1 + 1;
      left_extend = pair->start1 - prev_pair->start1;
      right_extend = pair->end1 - prev_pair->end1;
      intron_length = pair->start2 - prev_pair->end2 - 1 + overlap;
    }

    if (!prev_pair || prev_pair->entry1 != pair->entry1 || prev_pair->entry2 != pair->entry2 
	|| is_reverse(prev_pair) != is_reverse(pair) || intron_length > parameters->max_intron_length) {
      /* in new gene */
      intron_num = 0;
      prev_strand = -1;
      return -1;
    }
    left_pair = (is_reverse(pair) ? pair : prev_pair)->reversed_pair;
    right_pair = (is_reverse(pair) ? prev_pair : pair)->reversed_pair;
  }
  else {
    left_pair = prev_pair;
    right_pair = pair;
    overlap = prev_pair->end2 - pair->start2 + 1;
    left_extend = pair->start2 - prev_pair->start2;
    right_extend = pair->end2 - prev_pair->end2;
    intron_length = pair->start1 - prev_pair->end1 - 1 + overlap;
  }

  if (intron_length < parameters->min_intron_length || overlap > parameters->max_overlap
      || left_extend < parameters->min_exon_length || right_extend < parameters->min_exon_length) 
    return -1;

  intron_num++;
  if (print_flag)
    printf("\nEdges overlap: %d, intron_len: %d, intron #: %d", 
	   overlap, intron_length, intron_num);
  

  /* find most likely locations: near splice boundary, within (if no errors) or slightly without 
     for both locations to be outside, must have two errors -- implies inaccurate region of sequence
     so ignore!

     need to allow for gap chars within genomic part of alignment!! */

  /* do need to worry about edges of seq!! 
  if (!orient) {
    for (i = 0; i < splice_edge_len; i++) {
      if (left_string[i] != seq1[left_pair->end1 + 1 +  i])
	fatalError("seq out of register");
      if (right_string[-1 - i] != seq1[right_pair->start1 - 1 - i])
	fatalError("seq out of register2");
    }
  }
  */

  for (i = 0; i < 2; i++) {
    for (left_site = left_pair->query_data->splice_sites; left_site; left_site = left_site->next) {
      if (left_site->side || left_site->strand != i) continue;
      for (right_site = right_pair->query_data->splice_sites; right_site; right_site = right_site->next) {
	if (!right_site->side || right_site->strand != i) continue;
 	append_cand_splice(i, left_site, right_site, overlap, intron_length, left_pair, right_pair);
      }
    }
  }
  fail_flag = 0;
  if (best_cand_splice) {
    if (best_cand_splice->penalty < 0) fatalError("negative splice penalty");

    if (print_flag) {
      best_left_site = get_best_site(left_pair, 0);
      best_right_site = get_best_site(right_pair, 1);
      print_cand_splice(best_cand_splice);
      if (best_left_site != best_cand_splice->left_site && abs(best_left_site->intron_size - intron_length) <= 1) {
	printf(" site discrepancy: %d %d vs. %d %d ", 
	       best_cand_splice->left_site->strand, best_cand_splice->left_site->loc, best_left_site->strand, best_left_site->loc);
	fail_flag = 1;
      }
      if (best_left_site->llr > 0 && (best_left_site != best_cand_splice->left_site || abs(best_left_site->intron_size - intron_length) > 1)) {
	fail_flag = 1;
	printf("L:%.1f (%.1f) (%d)  ",  best_left_site->llr, best_cand_splice->left_site->llr, best_left_site->intron_size);
      }
      if (best_right_site->llr > 0 && (best_right_site != best_cand_splice->right_site || abs(best_right_site->intron_size - intron_length) > 1)) {
	fail_flag = 1;
	printf("R:%.1f (%.1f) (%d)  ",  best_right_site->llr, best_cand_splice->right_site->llr, best_right_site->intron_size);
      }
      if (next_best_cand_splice) 
	  printf(" penalty gap: %d ",  next_best_cand_splice->penalty - best_cand_splice->penalty);
    }
    if (incr_flag) {
      splice_penalty_hist[best_cand_splice->penalty > 40 ? 40 : best_cand_splice->penalty] += 1;
      i = best_cand_splice->missing_exon + 5;
      missing_exon_hist[i < 0 ? 0 : i > 55 ? 55 : i] += 1;
    }
  }
  if (print_flag && best_cand_splice && best_cand_splice->penalty <= 20 && best_cand_splice->missing_exon <= 30 && best_cand_splice->missing_exon >= parameters->min_exon_length)
    printf(" MISSING EXON (%d: %s): %d ", left_pair->entry2, get_id(left_pair->entry2), best_cand_splice->missing_exon);

  if (!best_cand_splice || best_cand_splice->penalty > 20 || next_best_cand_splice && best_cand_splice->penalty + min_penalty_gap > next_best_cand_splice->penalty) {
    if (print_flag) {
      printf(" failure ");
      if (best_cand_splice && best_cand_splice->penalty <= 20) {
	for (cand_splice = head_cand_splice; cand_splice; cand_splice = cand_splice->next) {
	  if (cand_splice == best_cand_splice) continue;
	  if (cand_splice->penalty < best_cand_splice->penalty + min_penalty_gap) {
	    print_cand_splice(cand_splice);
	    if (cand_splice->penalty == best_cand_splice->penalty) printf("* ");
	  }
	}
      }
      print_edges(splice_edge_len, left_pair, right_pair);
    }
  }
  else if (fail_flag) {
      print_edges(splice_edge_len, left_pair, right_pair);
  }
  else {
    if (prev_strand > -1 && best_cand_splice->strand != prev_strand)
      if (print_flag) printf(" STRAND CHANGE");

    prev_strand = best_cand_splice->strand; 
  }
  if (print_flag) {
    printf("\n");
  }
  penalty = !best_cand_splice ? 2000 : best_cand_splice->penalty;

  if (!orient && best_cand_splice && best_cand_splice->penalty <= 20) {
    first_cand_splice = 1;
    for (cand_splice = head_cand_splice; cand_splice; cand_splice = cand_splice->next) {
      if (cand_splice->penalty >= best_cand_splice->penalty + min_penalty_gap) continue;
      left_cdna_overhang = cand_splice->left_cdna_overhang;
      right_cdna_overhang = cand_splice->right_cdna_overhang;
      missing_exon = cand_splice->missing_exon;
      strand = cand_splice->strand;

      if (/* left_cdna_overhang >= 0 && right_cdna_overhang >= 0 && */ missing_exon >= parameters->min_exon_length
	  && missing_exon <= 30) {
	printf("\nCand missing exon (%d %s): ", pair->entry2, get_id(pair->entry2));
	print_cand_splice(cand_splice);
	print_cand_splice(best_cand_splice);
	if (first_cand_splice) {
	  printf(" [intron: %d, overlap: %d %d..%d %d..%d] ", 
		 intron_length, overlap, prev_pair->end1, pair->start1, prev_pair->end2, pair->start2); 
	  first_cand_splice = 0;
	}
	strncpy(exon_seq, strand ? "AC" : "AG", 2);
	strncpy(exon_seq + 2, seq2 + left_pair->end2 + 1 - left_cdna_overhang, missing_exon);
	strncpy(exon_seq + missing_exon + 2, strand ? "CT" : "GT", 2);
	exon_seq[missing_exon + 4] = 0;
	start = left_pair->end1 + 1 - left_cdna_overhang;
	end = right_pair->start1 - 1 + right_cdna_overhang - (missing_exon + 3);
	for (i = start, max_score = 0; i <= end; i++) {
	  /*      printf("%c", seq1[i]); */
	  if (exon_seq[0] == seq1[i] && exon_seq[1] == seq1[i + 1]
	      || exon_seq[missing_exon + 2] == seq1[i + missing_exon + 2] 
	      && exon_seq[missing_exon + 3] == seq1[i + missing_exon + 3]) {
	    for (j = n = m = score = 0; j < missing_exon + 4; j++) {
	      if (exon_seq[j] == 'N') n++;
	      else if (exon_seq[j] == seq1[i + j]) {
		score++;
		m++;
	      }
	      else score -= 2;
	      if (score < 0) score = 0;
	      else if (max_score < score) {
		max_score = score;
		best_start = i;
		n_best = 1;
		best_m = m;
	      }
	      else if (max_score == score)
		n_best++;
	    }
	  }
	}
	if (max_score >= 7) {
	  printf(" ");
	  for (j = 0; j < missing_exon + 4; j++) {
	    c = exon_seq[j]; 
	    printf("%c", c == seq1[best_start + j] ? c : tolower(c));
	  }
	  printf(" %d %d ***(%d,%d)  (%d hits)", best_start - start + 2, end - best_start + 2, best_m, max_score, n_best);

	  /* CHECK potential lead/trail pieces (more likely than small internal!!) */
	}
      }
    }
  }

  free_cand_splices();

  return penalty;
}

check_best_site(pair, side, cand_pair)
  Aligned_pair *pair;
  int side;
  Cand_pair *cand_pair;
{
  Splice_site *site;
  Segment *segment;
  char *get_id();
  int cdna_overhang;
  Splice_site *get_best_site();

  site = get_best_site(pair, side);
  if (!site || site->llr <= 0) return;
 
  for (segment = cand_pair->band_segments; segment; segment = segment->next) { /* NOT RELIABLE -- BECAUSE SCORE THRESHOLD MAY HAVE BEEN MISSED */
    if (segment->start < site->offset && segment->end > site->offset) break;
  }
  if (segment) return;

  cdna_overhang = site->side ? site->loc : -site->loc;

  if (cdna_overhang > 0) cdna_overhang -= site->cdna_overhang_adjust;
  printf("\nFound %d %2d %2d %2d %4d   %8d %6.2f  %s %d..%d %d..%d", 
	 site->side,  cdna_overhang, site->penalty, site->m_len, site->intron_size, 
	 site->offset, site->llr, get_id(pair->entry2), pair->start2, pair->end2, pair->start1, pair->end1); /* , segment ? "" : "***" */ /* not quite right!! */
}

Splice_site *get_best_site(pair, side)
  Aligned_pair *pair;
  int side;
{
  float best_llr;
  Splice_site *best_site, *site;

  best_llr = -10000;
  best_site = 0;
  for (site = pair->query_data->splice_sites; site; site = site->next) {
    if (site->side != side) continue;
    if (best_llr < site->llr) {
      best_llr = site->llr;
      best_site = site;
    }
  }
  return best_site;
}

int extend[3];

find_left_sites(pair, cand_pair, seq1, seq2)
  Aligned_pair *pair;
  Cand_pair *cand_pair;
  char *seq1, *seq2;
{
  int i, j, k, c, d, cum, max_len, loc;
  char ss_seq[20];
  char *string;
  Splice_site *site;

  for (site = pair->query_data->splice_sites; site; site = site->next) 
    if (!site->side) return; /* previously done! */

  for (i = cum = 0; i < splice_edge_len; i++) {
    cdna_indel_adjust[i] = cum += pair->query_data->edges12[7 * splice_edge_len - i - 1] == '-';
  }
  for (i = cum = 0; i < splice_edge_len; i++) {
    genome_indel_adjust[i] = cum += pair->query_data->edges12[3 * splice_edge_len - i - 1] == '-';
  }
  for (i = cum = 0; i < splice_edge_len; i++) {
    c = pair->query_data->edges12[7 * splice_edge_len - i - 1];
    d = pair->query_data->edges12[3 * splice_edge_len - i - 1];
    score_adjust[i] = cum += c == d ? 1 : c == '-' || d == '-' ? -4 : c == 'N' || d == 'N' ? 0 : -2;

    /* assumes default scoring */

    if (cum < 0) {
      printf("\nLeft %d %d", i, cum);
      fatalError("negative alignment score");
    }
  }  

  string = pair->query_data->edges12 + 3 * splice_edge_len; 

  /* check for single mismatch/inserted/deleted base */
  for (i = 0; i < splice_edge_len; i++) {
    c = string[4 * splice_edge_len + i];
    d = string[i + 1];
    if (c != d) break;
  }
  extend[0] = i; /* extension past deleted base in read */

  for (i = 0; i < splice_edge_len; i++) {
    c = string[4 * splice_edge_len + i + 1];
    d = string[i];
    if (c != d) break;
  }
  extend[1] = i; /* extension past inserted base in read */

  for (i = 0; i < splice_edge_len; i++) {
    c = string[4 * splice_edge_len + i + 1];
    d = string[i + 1];
    if (c != d) break;
  }
  extend[2] = i; /* extension past mismatch read */

  for (i = -splice_edge_len; i < splice_edge_len - 1; i++) {
    if (i < 0 && string[i] == '-') continue;

    max_len = 9;
    strncpy(ss_seq, i < 0 ? string + i : seq1 + pair->end1 + 1 + i, max_len);
    ss_seq[max_len] = 0;

    for (j = k = 0; j <= max_len; j++) {
      if (ss_seq[j] != '-') {
	ss_seq[k] = ss_seq[j];
	if (!ss_seq[k]) break;
	k++;
      }
    }

    if (i >= 0 || string[i + 1] != '-') {
      loc = i - (i > -splice_edge_len && i <= 0 && string[i - 1] == '-');
      try_splice_site(seq1, seq2, pair, 0, loc, ss_seq, k);
    }
    else {
      loc = i;
      try_splice_site(seq1, seq2, pair, 0, loc, ss_seq, k);
    }
  }
  max_len = 9;
  ss_seq[max_len] = 0;
  for (i = splice_edge_len - 1; i < 3 * splice_edge_len; i++) {
    strncpy(ss_seq, seq1 + pair->end1 + 1 + i, max_len);
    try_splice_site(seq1, seq2, pair, 0, i, ss_seq, k);
  }

  check_best_site(pair, 0, cand_pair);

}

find_right_sites(pair, cand_pair, seq1, seq2)
  Aligned_pair *pair;
  Cand_pair *cand_pair;
  char *seq1, *seq2;
{
  int i, j, k, c, d, cum, end, max_len, loc;
  char ss_seq[20];
  char *string;
  Splice_site *site;

  for (site = pair->query_data->splice_sites; site; site = site->next) 
    if (site->side) return; /* previously done! */

  for (i = cum = 0; i < splice_edge_len; i++) {
    cdna_indel_adjust[i] = cum += pair->query_data->edges12[5 * splice_edge_len + i] == '-';
  }
  for (i = cum = 0; i < splice_edge_len; i++) {
    genome_indel_adjust[i] = cum += pair->query_data->edges12[splice_edge_len + i] == '-';
  }

  for (i = cum = 0; i < splice_edge_len; i++) {
    c = pair->query_data->edges12[5 * splice_edge_len + i];
    d = pair->query_data->edges12[splice_edge_len + i];
    score_adjust[i] = cum += c == d ? 1 : c == '-' || d == '-' ? -4 : c == 'N' || d == 'N' ? 0 : -2;
    if (cum < 0) {
      printf("\nRight %d %d", i, cum);
      fatalError("negative alignment score");
    }
  }

  string = pair->query_data->edges12 + splice_edge_len;

  /* check for single mismatch/inserted/deleted base */
  for (i = 0; i < splice_edge_len; i++) {
    c = string[4 * splice_edge_len - i - 1];
    d = string[-i - 2];
    if (c != d) break;
  }
  extend[0] = i; /* extension past deleted base in read */

  for (i = 0; i < splice_edge_len; i++) {
    c = string[4 * splice_edge_len - i - 2];
    d = string[-i - 1];
    if (c != d) break;
  }
  extend[1] = i; /* extension past inserted base in read */

  for (i = 0; i < splice_edge_len; i++) {
    c = string[4 * splice_edge_len - i - 2];
    d = string[-i - 2];
    if (c != d) break;
  }
  extend[2] = i; /* extension past mismatch read */

  for (i = -splice_edge_len; i < splice_edge_len - 1; i++) {

    if (string[i] == '-') continue;
    end = string[i + 1] == '-' ? i + 2 : i + 1;
    max_len = end + 1 + splice_edge_len < 9 ?  end + 1 + splice_edge_len : 9;
    strncpy(ss_seq, string + end + 1 - max_len, max_len);
    ss_seq[max_len] = 0;
    for (j = k = 0; j <= max_len; j++) {
      if (ss_seq[j] != '-') {
	ss_seq[k] = ss_seq[j];
	if (!ss_seq[k]) break;
	k++;
      }
    }
    if (k >= 2) {
      if (string[i + 1] != '-') {
	loc = i + 2 + (i < splice_edge_len - 2 && string[i + 2] == '-');
	try_splice_site(seq1, seq2, pair, 1, loc, ss_seq, k);
      }
      else if (i < splice_edge_len - 2) {
	loc = i + 3;
	try_splice_site(seq1, seq2, pair, 1, loc, ss_seq, k);
      }
    }
  }
  max_len = 9;
  ss_seq[max_len] = 0;
  for (i = -3 * splice_edge_len; i < -splice_edge_len; i++) {
    strncpy(ss_seq, seq1 + pair->start1 + i + 2 - max_len, max_len);
    loc = i + 2;
    try_splice_site(seq1, seq2, pair, 1, loc, ss_seq, max_len);
  }

  check_best_site(pair, 1, cand_pair);
}

try_splice_site(seq1, seq2, pair, side, loc, seq, len)
     Aligned_pair *pair;
     int side, loc, len;
     char *seq1, *seq2, *seq;
{
  Splice_site *splice_site;

  /* from Levine & Durbin, for U12 sites, ..AC only occurred with AT;  the AT is actually ATATCC;
   categories found (decreasing order) were GTAG, ATAC, ATAG, GTAT, ATAT, GTGG, ATAA, GTAA 
  */

  if (side == 0) {
    if (!strncmp(seq, "GT", 2)) {
      increment_splice_site(seq1, seq2, pair, 0, side, loc, seq, 0);
    }
    if (!strncmp(seq, "GC", 2)) {
      increment_splice_site(seq1, seq2, pair, 0, side, loc, seq, 1);
    }
    if (!strncmp(seq, "ATATCC", 6)) {
      increment_splice_site(seq1, seq2, pair, 0, side, loc, seq, 2);
    }
    if (!strncmp(seq, "CT", 2)) {
      increment_splice_site(seq1, seq2, pair, 1, side, loc, seq, !strncmp(seq, "CTC", 3) ? 4 : 3);
    }
    if (!strncmp(seq, "GT", 2)) {
      increment_splice_site(seq1, seq2, pair, 1, side, loc, seq, 5);
    }
  }

  else {
    if (!strncmp(seq + len - 2, "AG", 2)) {
      increment_splice_site(seq1, seq2, pair, 0, side, loc, seq, !strncmp(seq + len - 3, "GAG", 3) ? 4 : 3);
    }
    if (!strncmp(seq + len - 2, "AC", 2)) {
      increment_splice_site(seq1, seq2, pair, 0, side, loc, seq, 5); 
    }
    if (!strncmp(seq + len - 2, "AC", 2)) {
      increment_splice_site(seq1, seq2, pair, 1, side, loc, seq, 0);
    }
    if (!strncmp(seq + len - 2, "GC", 2)) {
      increment_splice_site(seq1, seq2, pair, 1, side, loc, seq, 1);
    }
    if (!strncmp(seq + len - 6, "GGATAT", 6)) {
      increment_splice_site(seq1, seq2, pair, 1, side, loc, seq, 2); 
    }
  }
}

/*
  likelihood ratio scoring of candidate matches:
  prob under intron model: 
     prob(data at alignment end|intron) * prob(intron_size|intron) * prob(matchlength|intron) = a * b * c

  under non-intron model: 

     prob(data at alignment end|no intron) * prob(matchlength|no intron) = d * e

  a: log_2 is -penalty (bits) (+ -1 bit for strand choice). This is clearly an approximation!! In particular is assuming
           that prob of extension of score r past splice site is 4^-(max(r - 1,0)) which cannot be right -- is lower than
           this (and should take length, not just score, into account). Also the canonical sites (GT, HAG) are treated as prob 1.0 whereas are slightly less
	   (altho prob not much; GC and GAG are given probs of 1/16, prob too high).
  b: assume geometric, is q * (1-q)^n  where q = 1/mean intron size, n = intron size. Again -- not a very good approx!!
  c: depends on read error rate; but will take to be 1.0
  d: prob of sequence (= .25 * .25) (only looks at 2 bases -- so again approximate; should correct for GAG in particular).
  e: .75 (.25)^m

  SHOULD TABULATE VALUES OCCURRING -- TO GET MORE REASONABLE PROb ESTIMATES


then log(abc/de) = -penalty - 1 + log(q/.75) + n log(1 - q) - (m + 2) log(.25) = -penalty + log(q/.75) + n log(1 - q) + 2 * m + 3
This can be considered a test of the whole data (since match probs at non-sites cancel).

This doesn't really address multiple testing issue; 
*/

increment_splice_site(seq1, seq2, pair, strand, side, loc, seq, type)
     Aligned_pair *pair;
     int strand, side, loc, type;
     unsigned char *seq1, *seq2, *seq;
{
  Splice_site *splice_site, *splice_site2;
  char *our_alloc();
  int i, j, n, c_len, test_len, align_overhang, penalty, cdna_overhang, genome_overhang, max_i;
  int agree_flag, cdna_displace;
  double best_llr, max_llr;

  c_len = 15;
  test_len = c_len + 3;

  splice_site = (Splice_site *)our_alloc(sizeof(Splice_site)); /* do more efficiently!! */ 

  splice_site->next = pair->query_data->splice_sites;
  pair->query_data->splice_sites = splice_site;

  splice_site->loc = loc;
  splice_site->type = type;
  splice_site->strand = strand;
  splice_site->side = side;  
  splice_site->m_len = 0;
  splice_site->llr = -1000;
  splice_site->intron_size = 0;
  splice_site->cdna_overhang_adjust = splice_site->genome_overhang_adjust = 0;
  splice_site->j_offset = 0;
  penalty = 0;

  if (type == 1 || type == 4) penalty += 6; /* 1/64 for GC or GAG -- probably too high, at least for C. elegans!! */

  if (type == 2) penalty += 10; /* U12 penalty */

  genome_overhang = cdna_overhang = align_overhang = side ? loc : -loc;
  if (align_overhang > 0) {
    n = score_adjust[align_overhang - 1];
    /* 2 bits per extending base (allow one for free -- due to splice site prefs) */
    if (n > 0)
      penalty += 2 * n - 1;
    cdna_overhang -= splice_site->cdna_overhang_adjust = cdna_indel_adjust[align_overhang - 1];
    if (cdna_overhang < 0) fatalError("read overhang adjustment");
    genome_overhang -= splice_site->genome_overhang_adjust = genome_indel_adjust[align_overhang - 1];
    if (genome_overhang < 0) fatalError("read overhang adjustment");
  }
  /* strcpy(splice_site->seq, seq); */

  splice_site->overhang_check[0] = -align_overhang > extend[0] + 1;
  splice_site->overhang_check[1] = -align_overhang > extend[1];
  splice_site->overhang_check[2] = -align_overhang > extend[2] + 1;

  splice_site->penalty = penalty;
  
  splice_site->offset = -100000; /* MAKE THIS CLEARLY OUT OF RANGE!! */
  
  if (cdna_overhang < -5 || type == 5 || type == 2
      /* || (splice_site->overhang_check[0] && splice_site->overhang_check[1] && splice_site->overhang_check[2])*/ ) 
    return; /* don't allow signif truncated alignment (implying multiple seq error), U12 introns, or truncation due to more than one data error
	       for this scan */

  if (cdna_overhang < 0) penalty += 5; /* could also exclude these altogether (if want to prevent seq errors!!) */

  exon_seq[test_len] = 0;
  exon_seq[1] = strand ? 'C' : 'G';
  exon_seq[0] = side ? 'T' : 'A';

  /* NEED FOLLOWING: CHECK SEQUENCE LENGTHS; BETTER COORDINATE DETERMINATION; CHECKING AGAINST EXONS ALREADY FOUND 
   SPEED UP !! */

  splice_site->intron_size = parameters->max_intron_length;
  max_llr = base_llr - penalty + 2 * (test_len - 1);

  cdna_displace = -cdna_overhang;
  /*
  if (cdna_displace > 0 && splice_site->overhang_check[2]) { 
    if (!splice_site->overhang_check[0]) 
      cdna_displace -= 1;
    else if (!splice_site->overhang_check[1]) 
      cdna_displace += 1; 
  }
  */
  if (!side) {
    /* check if out of range */
    for (j = 2, i = pair->end2 + 1 + cdna_displace - 1; j < test_len; )
      exon_seq[j++] = seq2[i++];
    /* maximum i which could allow a positive LLR: max_llr + (max_i - pair->end1 + 1 + align_overhang) * intron_coeff >= 0
       => max_i <= -max_llr / intron_coeff + pair->end1 - 1 - align_overhang */
    max_i = -max_llr / intron_coeff + pair->end1 - 1;
    if (max_i > pair->end1 + 10000) max_i = pair->end1 + 10000;
  }
  else { /* in reverse order now!! */
 
   for (j = 2, i = pair->start2 - 1 - cdna_displace + 1; j < test_len; )
      exon_seq[j++] = seq2[i--];

    /* maximum i which could allow a positive LLR: max_llr + (pair->start1 + align_overhang - max_i + 1) * intron_coeff >= 0
       => max_i >= max_llr / intron_coeff + pair->start1 + align_overhang + 1 */
    max_i = max_llr / intron_coeff + pair->start1 + 1;
    if (max_i < pair->start1 - 10000) max_i = pair->start1 - 10000;
  }

  scan(splice_site, pair, seq1, genome_overhang, cdna_overhang, max_i, test_len);
}

/* input: splice_site, c_len, pair, seq1, exon_seq, genome_overhang, cdna_overhang, base_llr, intron_coeff, max_i  */

scan(splice_site, pair, seq1, genome_overhang, cdna_overhang, max_i, test_len)
  Splice_site *splice_site;
  Aligned_pair *pair;
  unsigned char *seq1;
  int genome_overhang, cdna_overhang, max_i, test_len;
{
  int i, i1, j, max_j, score, max_score, offset_displace, intron_adjust, intron_size, j_offset;
  double best_llr, base, llr;

  base = base_llr - (splice_site->penalty + (cdna_overhang < 0 ? 5 : 0));

  splice_site->llr = best_llr = -20.0 + base;
 
  offset_displace = splice_site->side ? -2 - (pair->start2 + cdna_overhang - 1) : 2 - (pair->end2 + 1 - cdna_overhang);
  intron_adjust = splice_site->side ? pair->start1 + genome_overhang + 1 : -pair->end1 + 1 + genome_overhang;

  if (!splice_site->side) {
    for (i = pair->end1 + 1; i < max_i; i++) {
      if (exon_seq[0] != seq1[i] || exon_seq[1] != seq1[i + 1]) continue; /* could also eliminate GAG sites at this point */
      /* try all three exon_seqs here!! */
      for (j_offset = -1; j_offset < 2; j_offset++) {
	for (j = 3 + j_offset, max_j = score = max_score = 2, i1 = i + 2; j < test_len; j++, i1++) 
	  if (exon_seq[j] != 'N') {
	    score += (exon_seq[j] == seq1[i1] ? 1 : -2);
	    /*	  if (score < 0) break; */
	    if (max_score < score) {
	      max_score = score;
	      max_j = j - 1 - j_offset;
	    }
	  }
	if (max_score < 4) continue;

	intron_size = intron_adjust + i;
	llr = intron_size * intron_coeff + 2 * max_score + base;
	if (best_llr < llr && intron_size >= parameters->min_intron_length) {
	  splice_site->offset = i + offset_displace;
	  splice_site->m_len = max_j; 
	  splice_site->intron_size = intron_size;
	  splice_site->llr = best_llr = llr;
	  splice_site->j_offset = j_offset;
	  if (best_llr > 20.0) return;
	}
      }
    }
  }
  else {
    for (i = pair->start1 - 1; i > max_i; i--) {
      if (exon_seq[0] != seq1[i] || exon_seq[1] != seq1[i - 1]) continue;
      for (j_offset = -1; j_offset < 2; j_offset++) {
	for (j = 3 + j_offset, max_j = score = max_score = 2, i1 = i - 2; j < test_len; j++, i1--)
	  if (exon_seq[j] != 'N') {
	    score += (exon_seq[j] == seq1[i1] ? 1 : -2);
	    if (max_score < score) {
	      max_score = score;
	      max_j = j - 1 - j_offset;
	    }
	  }
	if (max_score < 4) continue;

	intron_size = intron_adjust - i; 
	llr = intron_size * intron_coeff + 2 * max_score + base;
	if (best_llr < llr && intron_size >= parameters->min_intron_length) {
	  splice_site->offset = i + offset_displace;
	  splice_site->m_len = max_j;
	  splice_site->intron_size = intron_size;
	  splice_site->llr = best_llr = llr;
	  splice_site->j_offset = j_offset;
	  if (best_llr > 20.0) return;
	}
      }
    }
  }
}

print_splice_hists()
{
  int i;

  printf("\n\nSplice penalty histogram:");
  for (i = 0; i <= 40; i++)
    if (splice_penalty_hist[i])
      printf("\n%2d %5d", i, splice_penalty_hist[i]);
  printf("\n\nMissing exon size histogram:");
  for (i = 0; i < 56; i++)
    if (missing_exon_hist[i])
      printf("\n%2d %5d", i - 5, missing_exon_hist[i]);
  printf("\n\nIntron length discrepancy histogram:");
  for (i = 0; i < 21; i++)
    if (length_discrep_hist[i])
      printf("\n%2d  %5d", i, length_discrep_hist[i]);

}

/* analyze multiple matches for a given query -- to discriminate among them:
 idea: only interested in multiply covered regions in query (& ignore ends of matches?) 
 in 1st pass thru matches, identify positions within those which are discrepant in some matches, not in others -- these are informative.
 in 2d pass thru matches, identify how many agreements and disagreements at informative sites
         eliminate matches which have too many disagreements; don't eliminate matches which are only informative ones at a site 
 issues: can be more than one diff at same position in alignment: this affects interaction with coverage. (now handle indels & subs differently)
  discreps near ends: maybe should use these to trim; then want to check if anything remains
  indels can show up in different positions, depending on orientation
  basecalling errors

*/


int max_length, max_n_matches, n_matches, n_seg_equiv;
int *coverage, *new_coverage, *sub_counts, *os_clusters, *os_sub_clusters, *indel_counts, *subs, *used_marks, *best_chain[2], *best_chain2[2], *best_link[2], *best_link2[2]; 
                        /* subscripts are 0: quality restricted (based on informative subs); 1 unrestricted */
Aligned_pair **matches;
Seg_equiv *seg_equivs, **pair_seg_equivs;

set_matches(i_ptr, n_p, pair_pointers)
     int i_ptr, n_p;
     Aligned_pair **pair_pointers;
{
  char *our_alloc();
  int i, entry1, i_p, i_match;
  Aligned_pair *pair;

  entry1 = pair_pointers[i_ptr]->entry1;

  for (i_p = i_ptr, n_matches = 0; i_p < n_p; i_p++) { /* loop through relevant pairs */
    pair = pair_pointers[i_p];
    if (pair->entry1 != entry1) break;
    if (!pair->score) continue;
    if (!pair->reject_flags || pair->reject_flags & 4) continue; 
    n_matches++;
  }

  if (max_n_matches < n_matches) {

    if (max_n_matches) {
      our_free(subs);
      our_free(matches);
      our_free(seg_equivs);
      our_free(pair_seg_equivs);
      our_free(os_clusters);
      our_free(os_sub_clusters);
      our_free(used_marks);
      for (i = 0; i < 2; i++) {
	our_free(best_chain[i]);
	our_free(best_chain2[i]);
	our_free(best_link[i]);
	our_free(best_link2[i]);
      }
    }

    subs = (int *)our_alloc(n_matches * sizeof(int));
    matches = (Aligned_pair **)our_alloc(n_matches * sizeof(Aligned_pair *));
    seg_equivs = (Seg_equiv *)our_alloc(n_matches * sizeof(Seg_equiv));
    pair_seg_equivs = (Seg_equiv **)our_alloc(n_matches * sizeof(Seg_equiv *));
    os_clusters = (int *)our_alloc(n_matches * sizeof(int));
    os_sub_clusters = (int *)our_alloc(n_matches * sizeof(int));
    used_marks = (int *)our_alloc(n_matches * sizeof(int));
    for (i = 0; i < 2; i++) {
      best_chain[i] = (int *)our_alloc(n_matches * sizeof(int));
      best_chain2[i] = (int *)our_alloc(n_matches * sizeof(int));
      best_link[i] = (int *)our_alloc(n_matches * sizeof(int));
      best_link2[i] = (int *)our_alloc(n_matches * sizeof(int));
    }
    max_n_matches = n_matches;
  }

  for (i_p = i_ptr, i_match = 0; i_p < n_p; i_p++) { /* loop through relevant pairs */
    pair = pair_pointers[i_p];
    if (pair->entry1 != entry1) break;
    if (!pair->score) continue;
    if (!pair->reject_flags || pair->reject_flags & 4) continue; 
    matches[i_match++] = pair;
  } 
  if (i_match != n_matches) fatalError("match count");

  return n_matches;
}

check_max_length(length)
     int length;
{
  char *our_alloc();

  if (max_length >= length) return;
  notify("ALLOCATING SPLICE INFO");
  if (max_length) {
    our_free(coverage);
    our_free(new_coverage);
    our_free(sub_counts);
    our_free(indel_counts);
  }
  coverage = (int *)our_alloc((length + 1) * sizeof(int));
  new_coverage = (int *)our_alloc((length + 1) * sizeof(int));
  sub_counts = (int *)our_alloc((length + 1) * sizeof(int));
  indel_counts = (int *)our_alloc((length + 1) * sizeof(int));
  max_length = length;
}

filter_matches(i_ptr, n_p, pair_pointers, print_flag)
     int i_ptr, n_p, print_flag;
     Aligned_pair **pair_pointers;
{
  int i, c, entry1, length, n_sub_inform, n_indel_inform, n_poss, n_found, n_0, n_2, inform_subs, inform_indels, n_poss_subs, n_poss_indels, left_sub, right_sub, max_n_subs, n_delete_match, found_sig;
  int i_match, j_match, hist[21], best_subs, cum, intron_length, overlap, max_2, max_depth, old_n_matches;
  Aligned_pair *pair, *prev_pair;

  max_n_subs = 10; /* maximum no. of informative substitutions */

  if (!set_matches(i_ptr, n_p, pair_pointers)) return;

  entry1 = pair_pointers[i_ptr]->entry1;
  length = get_seq_length(entry1);

  check_max_length(length);

  /* improve following: multiple matches within read to same genomic segment are causing apparent truncation of parse */


  for (i_match = 0; i_match < n_matches; i_match++) subs[i_match] = 0;

  max_2 = n_sub_inform = n_indel_inform = 0; 

  if (print_flag) printf("\n\n%d matches;", n_matches);

  if (n_matches > 1) {
    for (i = 0; i <= length; i++) {
      coverage[i] = sub_counts[i] = indel_counts[i] = 0;
      new_coverage[i] = -1;
    }

    for (i_match = 0; i_match < n_matches; i_match++) {
      pair = matches[i_match];
      if (pair->entry1 != entry1) fatalError("match error");
      coverage[pair->start1] += 1;
      coverage[pair->end1 + 1] -= 1; 
    }

    for (i = 0; i < 21; i++) hist[i] = 0;
    for (i = n_2 = max_2 = 0; i <= length; i++) {
      if (i) coverage[i] += coverage[i - 1];
      c = coverage[i];
      hist[c < 20 ? c : 20] += 1;
      if (c >= 2) {
	n_2++;
	if (max_2 < n_2) max_2 = n_2;
      }
      else 
	n_2 = 0;
    }

    if (coverage[length]) fatalError("coverage calculation");

    if (print_flag) {
      for (c = 20, cum = 0; c >= 0; c--) {
	cum += hist[c];
	if (cum > splice_edge_len) break; /* not appropriate algorithm!! -- */
      }
      
      printf(" maximum coverage depth: %d; ", c);
    }

    if (max_2 > splice_edge_len) { /* maximum doubly covered stretch is > splice_edge_len */

      /* initial pass -- to find diffs */
      for (i_match = 0; i_match < n_matches; i_match++) { 
	pair = matches[i_match];
	count_diffs(pair, sub_counts, indel_counts, 0, &inform_subs, &inform_indels, &left_sub, &right_sub);
      } 

      /* note that following doesn't take into account fact that different matches can have independent types of errors at a position -- so may overcount informative sites */

      for (i = n_sub_inform = n_indel_inform = 0; i < length; i++) {
	if (sub_counts[i]) {
	  if (sub_counts[i] < coverage[i]) { /* so coverage must be at least 2 */
	    n_sub_inform++;
	    sub_counts[i] = 1;
	  }
	  else sub_counts[i] = 0;
	}
	if (indel_counts[i]) {
	  if (indel_counts[i] < coverage[i]) { 
	    n_indel_inform++;
	    indel_counts[i] = 1;
	  }
	  else indel_counts[i] = 0;
	}
      }

      if (print_flag) {
	printf(" informative sites: %d subs, %d indels;", n_sub_inform, n_indel_inform);
      }
    }
  }

  if (n_sub_inform) {

    /* what about discrepancies near ends of alignment -- those less relevant (due to pseudo-extensions) */
    /* make use of fact that start - 1 and end + 1 must be discrepant?? -- but probably won't want to use those! */
    for (i = 1; i <= length; i++) {
      sub_counts[i] += sub_counts[i - 1]; /* now gives cumulative no. ; */
      indel_counts[i] += indel_counts[i - 1]; 
    }

    for (i_match = n_delete_match = 0; i_match < n_matches; i_match++) { 
      pair = matches[i_match];

      /* first check whether there are any!! */
      n_poss_subs = sub_counts[pair->end1] - (pair->start1 ? sub_counts[pair->start1 - 1] : 0);
      n_poss_indels = indel_counts[pair->end1] - (pair->start1 ? indel_counts[pair->start1 - 1] : 0);
      inform_subs = inform_indels = 0;
      if (n_poss_subs + n_poss_indels)
	count_diffs(pair, sub_counts, indel_counts, 1, &inform_subs, &inform_indels, &left_sub, &right_sub);
      if (print_flag > 1)
	printf(" (%d of %d, %d of %d -- span %d)", inform_subs, n_poss_subs, inform_indels, n_poss_indels, right_sub - left_sub - 1);
      if (n_poss_subs < inform_subs || n_poss_indels < inform_indels) fatalError("sub or indel counts");
      subs[i_match] = inform_subs; 
      if (inform_subs > max_n_subs) {
	pair->reject_flags |= 4;
	n_delete_match++;
	continue;
      }
      for (i = pair->start1; i <= pair->end1; i++) {
	if (new_coverage[i] == -1 || new_coverage[i] > inform_subs) 
	  new_coverage[i] = inform_subs;
      }
    }

    if (n_delete_match) {
      if (print_flag)
	printf(" %d matches deleted;", n_delete_match);
    }

    for (i_match = 0; i_match < n_matches; i_match++) { 
      pair = matches[i_match];
      if (pair->reject_flags & 4) continue;
      for (i = pair->start1, n_0 = 0; i <= pair->end1; i++) {
	if (subs[i_match] <= new_coverage[i] + 1) n_0++; /* allow one extra informative mismatch than the best one in this region */
      }

      if (n_0 < 10) {
	pair->reject_flags |= 8;
	if (print_flag) printf(" 0");
      }
      else {
	if (print_flag) printf(" 1");
      }
    }
  }
}

sort_seg_equivs(seg_equiv1, seg_equiv2)
     Seg_equiv *seg_equiv1, *seg_equiv2;
{
  return seg_equiv1->start - seg_equiv2->start;
}

analyze_multiple(i_ptr, n_p, pair_pointers, print_flag)
     int i_ptr, n_p, print_flag;
     Aligned_pair **pair_pointers;
{
  int i_cluster, i_sub_cluster, pair_length, intron_length, overlap, i_match, j_match, next_j_match, n_0, flag, new_chain, short_flag, n, i_seg_equiv, j_seg_equiv, left_extend, right_extend, max_gap, start, end, excess, min_excess;
  Seg_equiv *seg_equiv, *seg_equiv2;
  Aligned_pair *pair, *prev_pair, *next_pair;
  Splice_site *get_best_site();
  Splice_site *site, *next_site;
  Segment *head_all_segment, *head_chain0_segment, *segment, *segment2, *insert_segment();
  Segment *find_gap_segment_list();

  if (!set_matches(i_ptr, n_p, pair_pointers)) return;

  for (i_match = i_cluster = i_sub_cluster = 0; i_match < n_matches; i_match++) { 
    pair = matches[i_match];
    if (i_match) {
      prev_pair = matches[i_match - 1];
      overlap = prev_pair->end1 - pair->start1 + 1;
      intron_length = pair->start2 - prev_pair->end2 - 1 + overlap;
    }
    if (i_match && prev_pair->entry1 == pair->entry1 && prev_pair->entry2 == pair->entry2 && is_reverse(prev_pair) == is_reverse(pair)
	&& intron_length <= parameters->max_intron_length) { /* making use of sorted order here!! */
      if (pair->start2 > prev_pair->end2) i_sub_cluster++;
    }
    else {
      i_cluster++;
      i_sub_cluster++;
    }

    os_clusters[i_match] = i_cluster; /* orientation-spacing clusters */
    os_sub_clusters[i_match] = i_sub_cluster; /* orientation-spacing subclusters (correspond to non-overlapping genomic segments involved in matches) 
					       potential splices must relate distinct subclusters within same cluster, and also preserve order within query 
					       quality of parse reflects 
					       1) whether omitted part of read is present in any other parse; 
					       2) quality of matches; 
					       3) quality of splice 
					      */
  }

  /* NEED TO DEAL WITH DUPLICATED SEGMENTS IN READ!! (SOME OF WHICH MIGHT BE ALT SPLICING) */
  /* equivalence classes in read: parts that consistently match within some genomic region (but what about duplicated segments?)
       1. group overlapping query segments into 'segments'
       2. group segments into equiv classes (via spliceability) 
       3. find parses (of non-rejected matches) that span the entire set of segments in an equiv class.
   */
  short_flag = group_matches();

  mark_save_block();

  head_all_segment = head_chain0_segment = 0;

  for (i_match = n_matches - 1; i_match >= 0; i_match--) { /* reverse order -- so chains in proper order */
    if (!pair_seg_equivs[i_match]) continue;
    pair = matches[i_match];
    head_all_segment = insert_segment(head_all_segment, pair->start1, pair->end1);
    pair_length = pair->end1 - pair->start1 + 1;
    best_chain[0][i_match] = (pair->reject_flags & 8) ? 0 : pair_length;
    best_chain[1][i_match] = pair_length;
    best_chain2[0][i_match] = best_chain2[1][i_match] = 0;
    best_link2[0][i_match] = best_link[0][i_match] = -1;
    best_link2[1][i_match] = best_link[1][i_match] = -1;
    for (j_match = i_match + 1; j_match < n_matches && os_clusters[i_match] == os_clusters[j_match]; j_match++) {
      if (!pair_seg_equivs[j_match]) continue;
      if (os_sub_clusters[i_match] == os_sub_clusters[j_match]) continue;
      next_pair = matches[j_match];
      overlap = pair->end1 - next_pair->start1 + 1;
      left_extend = next_pair->start1 - pair->start1;
      right_extend = next_pair->end1 - pair->end1;
      if (overlap > parameters->max_overlap || left_extend < parameters->min_exon_length
	  || right_extend < parameters->min_exon_length) continue;
      intron_length = next_pair->start2 - pair->end2 - 1 + overlap;
      if (intron_length < parameters->min_intron_length && abs(intron_length) > 5) continue; /* exclude matches that are neither adjacent,
											   nor imply plausible introns */
      /*
      n = analyze_splices(pair, next_pair, 0, 0, 1));
      if (n < 0 || n  > 10) {
	continue;
      }
      */
      if (!(pair->reject_flags & 8)) {
	new_chain = best_chain[0][j_match] + pair_length;
	if (best_chain[0][i_match] < new_chain) { /* probably not best statistic to rank on!! */
	  best_link2[0][i_match] = best_link[0][i_match];
	  best_chain2[0][i_match] = best_chain[0][i_match];
	  best_link[0][i_match] = j_match;
	  best_chain[0][i_match] = new_chain;
	}
	else if (best_chain2[0][i_match] < new_chain) {
	  best_link2[0][i_match] = j_match;
	  best_chain2[0][i_match] = new_chain;
	}
      }
      new_chain = best_chain[1][j_match] + pair_length;
      if (best_chain[1][i_match] < new_chain) { /* probably not best statistic to rank on!! */
	best_link2[1][i_match] = best_link[1][i_match];
	best_chain2[1][i_match] = best_chain[1][i_match];
	best_link[1][i_match] = j_match;
	best_chain[1][i_match] = new_chain;
      }
      else if (best_chain2[1][i_match] < new_chain) {
	best_link2[1][i_match] = j_match;
	best_chain2[1][i_match] = new_chain;
      }

      /*
      if (print_flag)
	printf(" (%d, %d, %d)%s", j_match, i_match, intron_length, (pair->reject_flags & 8) != (next_pair->reject_flags & 8) ? "*" : "");
      */
    }
  }

  for (i_match = n_matches - 1; i_match >= 0; i_match--) { /* reverse order -- to compress parses */
    if (!pair_seg_equivs[i_match]) continue;
    pair = matches[i_match];
    seg_equiv = pair_seg_equivs[i_match]->parent;
    if (seg_equiv->best_length[1] <= best_chain[1][i_match]) {
      seg_equiv->start_count[1] = seg_equiv->best_length[1] < best_chain[1][i_match] ? 1 : seg_equiv->start_count[1] + 1;
      seg_equiv->best_length[1] = best_chain[1][i_match];
      seg_equiv->best_start[1] = i_match;
    }
    if (pair->reject_flags & 8) continue;
    if (seg_equiv->best_length[0] <= best_chain[0][i_match]) {
      seg_equiv->start_count[0] = seg_equiv->best_length[0] < best_chain[0][i_match] ? 1 : seg_equiv->start_count[0] + 1;
      seg_equiv->best_length[0] = best_chain[0][i_match];
      seg_equiv->best_start[0] = i_match;
    }
  }


  for (i_match = n_matches - 1; i_match >= 0; i_match--) {
    if (!pair_seg_equivs[i_match]) continue;
    pair = matches[i_match];
    if (pair->reject_flags & 8) continue;
    seg_equiv = pair_seg_equivs[i_match]->parent;
    if (best_chain[0][i_match] > seg_equiv->best_length[0] - parameters->min_exon_length) { 
      /* doesn't allow for 'slop'; alternative is seg_equiv->length - parameters->min_exon_length) { */
      seg_equiv->n_chains += 1;
      if (print_flag) {
	printf(" %schain(%d)", 
	       best_chain[0][i_match] > seg_equiv->best_length[1] - parameters->min_exon_length ? "" : "*", 
	       seg_equiv->id);
      }
      intron_length = 0;
      next_site = 0;
      for (j_match = i_match; j_match > -1; j_match = next_j_match) {
	pair_seg_equivs[j_match]->parent->segment = insert_segment(pair_seg_equivs[j_match]->parent->segment,
								   matches[j_match]->start1, matches[j_match]->end1);
	head_chain0_segment = insert_segment(head_chain0_segment, matches[j_match]->start1, matches[j_match]->end1);
	next_j_match = best_link[0][j_match];
	if (print_flag) {
	  printf(":");
	  site = next_site ? next_site : get_best_site(matches[j_match]->reversed_pair, !is_reverse(matches[j_match]));
	  if (site) {
	    n = abs(site->intron_size - intron_length);
	    if (intron_length) 
	      length_discrep_hist[n > 20 ? 20 : n] += 1;
	    if (n > 2 && site->llr > 0)
	      printf("[%.1f%s,%d vs %d]", 
		     site->llr, 
		     site->llr >= 7.0 ? (intron_length ? "***" : "**") : "", 
		     site->intron_size, intron_length);
	  }
	  printf("%d", j_match);
	  if (next_j_match > -1) {
	    next_pair = matches[next_j_match];
	    pair = matches[j_match];
	    intron_length = next_pair->start2 - pair->end2 + pair->end1 - next_pair->start1;
	    next_site = get_best_site(matches[next_j_match]->reversed_pair, !is_reverse(matches[next_j_match]));
	    if (next_site)
	      intron_length += (next_site->cdna_overhang_adjust + next_site->j_offset) - next_site->genome_overhang_adjust; 
	  }
	  else intron_length = 0;
	  site = get_best_site(matches[j_match]->reversed_pair, is_reverse(matches[j_match]));
	  if (site) {
	    if (intron_length) 
	      intron_length += (site->cdna_overhang_adjust + site->j_offset) - site->genome_overhang_adjust; 

	    n = abs(site->intron_size - intron_length);
	    if (intron_length) 
	      length_discrep_hist[n > 20 ? 20 : n] += 1;

	    if (n > 2 && site->llr > 0)
	      printf("[%.1f%s,%d vs %d]", site->llr, site->llr >= 7.0 ? (intron_length ? "***" : "**") : "", site->intron_size, intron_length);
	  }
	}
	matches[j_match]->reject_flags |= 16;
	if (next_j_match > -1 && best_chain2[0][j_match] > best_chain[0][j_match] - 10) 
	  if (print_flag)
	    printf("[AMBIGUOUS LINK:%d:%d or %d:%d]", 
		   next_j_match, best_link[0][next_j_match], best_link2[0][j_match], best_link[0][best_link2[0][j_match]]);
      }
    }
    /*
    if (seg_equiv->best_length[1] < best_chain[0][i_match]) {
      seg_equiv->best_length[1] = best_chain[0][i_match];
      seg_equiv->best_start[1] = i_match;
    }
    */
  }

  if (print_flag) printf(" segments (equiv_class): ");
  for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    if (!seg_equiv->usable) continue; 
    if (print_flag)
      printf(" %d..%d (%d)", seg_equiv->start, seg_equiv->end, seg_equiv->parent->id);
  }

  for (i_seg_equiv = flag = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    if (!seg_equiv->usable) continue; 
    if (seg_equiv == seg_equiv->parent) {
      if (seg_equiv->n_chains != 1) {
	if (!flag) {
	  if (print_flag)
	    printf("; %s segment equiv classes: ", seg_equiv->n_chains > 1 ? "MULTIPLE": "UNSUCCESSFUL");
	  flag = 1;
	}
	if (print_flag)
	  printf(" class %d, %d chains", seg_equiv->id, seg_equiv->n_chains);
	if (!seg_equiv->n_chains) {
	  if (print_flag) {
	    printf(" best length: %d;",  seg_equiv->best_length[1]);
	    printf(" chain(%d)", seg_equiv->id);
	    if (seg_equiv->start_count[1] > 1) printf("MULTIPLE STARTS"); 
	  }
	  for (j_match = seg_equiv->best_start[1]; j_match > -1; j_match = best_link[0][j_match]) {
	    if (print_flag)
	      printf(":%d", j_match);
	    matches[j_match]->reject_flags |= 16;
	    if (best_link[0][j_match] > -1 && best_chain2[0][j_match] > best_chain[0][j_match] - 10) /* HARDWIRED -- CHANGE THIS */
	      if (print_flag)
		printf("[AMBIGUOUS LINK:%d or %d]", 
		       best_link[0][j_match], best_link2[0][j_match]);
	  }

	  next_pair = 0;
	  analyze_splices((Aligned_pair *)0, (Aligned_pair *)0, 1, 1, 1, (char *)0, (char *)0); /* to reset intron number */
	  for (j_match = seg_equiv->best_start[1]; j_match > -1; j_match = next_j_match) {
	    used_marks[j_match] = 1;
	    next_j_match = best_link[0][j_match];
	    if (next_j_match <= -1) break;
	    pair = next_pair ? next_pair : matches[j_match];
	    next_pair = matches[next_j_match];
	    if ((pair->reject_flags & 16) && !(pair->reject_flags & 8)
		&& (next_pair->reject_flags & 16) && !(next_pair->reject_flags & 8))
	      analyze_splices(pair, next_pair, 1, 1, 1, (char *)0, (char *)0);
	  }
	}
      }
    }
  }
  if (short_flag) printf(" SEGMENT LENGTHENED");

  for (i_match = 0; i_match < n_matches; i_match++) 
    used_marks[i_match] = 0;
  
  for (i_match = n_matches - 1; i_match >= 0; i_match--) {
    if (!pair_seg_equivs[i_match]) continue;
    pair = matches[i_match];
    if (pair->reject_flags & 8) continue;
    if (used_marks[i_match]) continue; /* SHOULD APPLY THIS ABOVE, ALSO -- TO PREVENT REUSE!! */
    seg_equiv = pair_seg_equivs[i_match]->parent;

    if (best_chain[0][i_match] > seg_equiv->best_length[0] - parameters->min_exon_length) { 
      next_pair = 0;
      analyze_splices((Aligned_pair *)0, (Aligned_pair *)0, 1, 1, 1, (char *)0, (char *)0); /* to reset intron number */
      for (j_match = i_match; j_match > -1; j_match = next_j_match) {
	used_marks[j_match] = 1;
	next_j_match = best_link[0][j_match];
	if (next_j_match <= -1) break;
	pair = next_pair ? next_pair : matches[j_match];
	next_pair = matches[next_j_match];
	if ((pair->reject_flags & 16) && !(pair->reject_flags & 8)
	    && (next_pair->reject_flags & 16) && !(next_pair->reject_flags & 8))
	  analyze_splices(pair, next_pair, 1, 1, 1, (char *)0, (char *)0);
      }
    }
  }

  max_gap = find_max_gap_list(head_chain0_segment, head_all_segment);
  if (max_gap >= parameters->min_exon_length) {
    printf("\nWARNING: parts of cdna missing from parses:");
    segment = find_gap_segment_list(head_chain0_segment, head_all_segment);
    for (; segment; segment = segment->next) printf(" %d..%d", segment->start, segment->end);
  }
  /* the following check probably not useful */
  for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    if (!seg_equiv->usable) continue; 
    if (seg_equiv != seg_equiv->parent) continue;
    if (abs(seg_equiv->best_length[0] - seg_equiv->best_length[1]) > 10 || abs(seg_equiv->best_length[0] - seg_equiv->length) > 10)
      printf("\nWARNING: seg_equiv lengths differ: %d: %d %d %d", 
	     seg_equiv->id, seg_equiv->best_length[0], seg_equiv->best_length[1], seg_equiv->length);
  }

  /* check whether any eliminated match had unique portion of read */
  for (i_match = 0; i_match < n_matches; i_match++) { 
    pair = matches[i_match];
    if (pair->score >= parameters->minscore) continue;
    if (pair_seg_equivs[i_match]) continue;
    min_excess = pair->end1 - pair->start1 + 1;
    for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
      seg_equiv = seg_equivs + i_seg_equiv;
      if (!seg_equiv->usable || !seg_equiv->parent->usable) continue;
      if (pair->start1 <= seg_equiv->end && pair->end1 >= seg_equiv->start) {
	start = pair->start1 > seg_equiv->start ? pair->start1 : seg_equiv->start;
	end = pair->end1 < seg_equiv->end ? pair->end1 : seg_equiv->end;
	excess =  (pair->end1 - pair->start1) - (end - start);
	if (min_excess > excess) min_excess = excess;
      }
    }
    if (min_excess >= parameters->min_exon_length) 
      printf("\nWARNING: DELETED UNIQUE MATCH: %d..%d  %d", pair->start1, pair->end1, min_excess);
  }    

  for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    if (seg_equiv != seg_equiv->parent) continue;
    if (!seg_equiv->usable) continue; 
    for (j_seg_equiv = i_seg_equiv + 1; j_seg_equiv < n_seg_equiv; j_seg_equiv++) {
      seg_equiv2 = seg_equivs + j_seg_equiv;
      if (seg_equiv2 != seg_equiv2->parent) continue;
      if (!seg_equiv2->usable) continue; 
      for (segment = seg_equiv->segment; segment; segment = segment->next)
	for (segment2 = seg_equiv2->segment; segment2; segment2 = segment2->next) 
	  if (segment->start < segment2->end - 5 && segment->end > segment2->start + 5)
	    printf("\nWARNING: overlapping segments from distinct equiv classes: %d..%d %d..%d",
		   segment->start, segment->end,segment2->start, segment2->end);
    }
  }

  for (i_match = 0; i_match < n_matches; i_match++) { 
    seg_equiv = pair_seg_equivs[i_match];
    if (!seg_equiv) continue;
    seg_equiv = seg_equiv->parent;
    for (j_match = i_match + 1; j_match < n_matches && os_clusters[i_match] == os_clusters[j_match]; j_match++) { 
      if (!pair_seg_equivs[j_match]) continue;
      if (seg_equiv != pair_seg_equivs[j_match]->parent) 
	printf("\nWARNING: clustered matches in different seg_equiv classes: %d %d", i_match, j_match);
      break; /* so only print out changes within a cluster */
    }
  }
  free_seg_blocks();
}

/* input:matches[], n_matches, os_clusters, os_sub_clusters
   modified: seg_equivs, pair_seg_equivs
   output: short_flag, n_seg_equiv, 
*/

/* 

  NEED: to make sure that same part of read not being used for two different genes; this not fully
        tested yet. AND all of read should be accounted for.

*/

group_matches()
{
  Segment *head_segment, *segment, *insert_segment();
  int short_flag, i_match, j_match, start, end, n_segments, i_seg_equiv, n, overlap, intron_length, left_extend, right_extend, excess, min_excess;
  Aligned_pair *pair, *prev_pair, *next_pair;
  Seg_equiv *seg_equiv, *seg_equiv2;
  char *get_id();

  mark_save_block();

  head_segment = 0;
  short_flag = 0;

  /* clear paralog matches are grouped into same classes */

  for (i_match = 0; i_match < n_matches; i_match++) { 
    pair = matches[i_match];
    if (pair->score < parameters->minscore) continue;
    start = pair->start1 + 10;
    end = pair->end1 - 10;
    if (start > end) {
      short_flag = 1;
      start = end = (start + end) / 2;
    }
    head_segment = insert_segment(head_segment, start, end);
  }

  for (segment = head_segment, n_segments = 0; segment; segment = segment->next, n_segments++) {
    seg_equiv = seg_equivs + n_segments;
    seg_equiv->start = segment->start - 10; 
    seg_equiv->end = segment->end + 10;
    seg_equiv->length = seg_equiv->end - seg_equiv->start + 1;
    seg_equiv->pair = 0;
  }
  i_seg_equiv = n_segments;

  /* marginal matches become separate classes (because don't want them to force joins) */

  for (i_match = 0; i_match < n_matches; i_match++) {
    pair = matches[i_match];
    if (pair->score >= parameters->minscore) continue;
    seg_equiv = seg_equivs + i_seg_equiv;
    seg_equiv->start = pair->start1;
    seg_equiv->end = pair->end1;
    seg_equiv->length = seg_equiv->end - seg_equiv->start + 1;
    seg_equiv->pair = pair;
    i_seg_equiv++;
  }
  n_seg_equiv = i_seg_equiv;
  qsort(seg_equivs, n_seg_equiv, sizeof(Seg_equiv), sort_seg_equivs);

  for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    seg_equiv->parent = seg_equiv;
    seg_equiv->id = seg_equiv->n_chains = seg_equiv->best_length[0] = seg_equiv->best_length[1] = 0;
    seg_equiv->best_start[0] = seg_equiv->best_start[1] = -1;
    seg_equiv->start_count[0] = seg_equiv->start_count[1] = 0;
    seg_equiv->usable = 1;
    seg_equiv->n_segments = 1;
    seg_equiv->segment = 0;
  }

  for (i_match = 0; i_match < n_matches; i_match++) {
    pair = matches[i_match];
    for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
      seg_equiv = seg_equivs + i_seg_equiv;
      if (pair->start1 >= seg_equiv->start && pair->end1 <= seg_equiv->end) {
	if (pair->score >= parameters->minscore && !seg_equiv->pair || seg_equiv->pair == pair) {
	  pair_seg_equivs[i_match] = seg_equiv;
	  break;
	}
      }
    }
    if (i_seg_equiv == n_seg_equiv) fatalError("missing segment");
  }

    /* GREEDY PARSE DOES NOT NEC. GIVE BEST RESULTS!! -- KEEP TRACK OF MULTIPLE POSSIBILITIES!! */
  /* following unnecessary if there is only one segment!! */
  for (i_match = n_matches - 1; i_match >= 0; i_match--) { /* reverse order -- so chains in proper order */
    pair = matches[i_match];
    for (j_match = i_match + 1; j_match < n_matches && os_clusters[i_match] == os_clusters[j_match]; j_match++) {
      if (os_sub_clusters[i_match] == os_sub_clusters[j_match]) continue;
      next_pair = matches[j_match];
      overlap = pair->end1 - next_pair->start1 + 1;
      left_extend = next_pair->start1 - pair->start1;
      right_extend = next_pair->end1 - pair->end1;
      /*
      if (!strcmp(get_id(pair->entry1), "B06_1672146_SL2_3_041229.f"))
	printf("\nHere: %d %d %d %d %d", i_match, j_match, overlap, left_extend, right_extend);
      */
      if (overlap > parameters->max_overlap || left_extend < parameters->min_exon_length
	  || right_extend < parameters->min_exon_length) continue;
                /* can cause small exon to be missed -- if spurious extension of other match */
      intron_length = next_pair->start2 - pair->end2 - 1 + overlap;
      /*
      if (!strcmp(get_id(pair->entry1), "B06_1672146_SL2_3_041229.f"))
	printf("\nHere: %d %d %d", i_match, j_match, intron_length);
      */
      if (intron_length < parameters->min_intron_length && abs(intron_length) > 5) continue; /* exclude matches that are neither adjacent,

																						   nor imply plausible introns */
      /* more stringent criterion for marginal matches? following is probably too stringent, because these matches
	 are likely to be 'genuine' if they don't overlap a better match!!
      if ((pair->score < parameters->minscore || next_pair->score < parameters->minscore)
	  && abs(overlap) > 10 && intron_length > 500) 
							
	continue;
      */
      /*
      if (!strcmp(get_id(pair->entry1), "B06_1672146_SL2_3_041229.f"))
	printf("\nHere3: %d %d", i_match, j_match);
      */
      for (seg_equiv = pair_seg_equivs[i_match]; seg_equiv != seg_equiv->parent; seg_equiv = seg_equiv->parent);
      for (seg_equiv2 = pair_seg_equivs[j_match]; seg_equiv2 != seg_equiv2->parent; seg_equiv2 = seg_equiv2->parent);
      if (seg_equiv != seg_equiv2) {
	if (!seg_equiv2->pair) { /* favor merging with non-marginal match */
	  seg_equiv->parent = seg_equiv2;
	  seg_equiv2->n_segments += seg_equiv->n_segments;
	}
	else {
	  seg_equiv2->parent = seg_equiv;
	  seg_equiv->n_segments += seg_equiv2->n_segments;
	}
      }
    }
  }

  /* if a marginally significant match is in one equiv class, but signif overlaps another, it is probably spurious,
     so delete */
  for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    for (seg_equiv2 = seg_equiv; seg_equiv2 != seg_equiv2->parent; seg_equiv2 = seg_equiv2->parent);
    seg_equiv->parent = seg_equiv2;
  }

  /*
  for (i_match = 0; i_match < n_matches; i_match++) { 
    pair = matches[i_match];
    pair_seg_equivs[i_match] = pair_seg_equivs[i_match]->parent;
  }
  */

  for (i_match = 0; i_match < n_matches; i_match++) { 
    pair = matches[i_match];
    if (pair->score >= parameters->minscore) continue;
    seg_equiv2 = pair_seg_equivs[i_match]->parent; 
    if (seg_equiv2->n_segments == 1) {
      seg_equiv2->usable = 0;
      pair_seg_equivs[i_match] = 0;
      continue; /* unlinked marginal match is ignored */
    }
    /* IS FOLLOWING NECESSARY?? */
    for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
      seg_equiv = seg_equivs + i_seg_equiv;
      if (seg_equiv->pair) continue;
      if (pair->start1 < seg_equiv->end - 10 && pair->end1 > seg_equiv->start + 10) {
	if (seg_equiv->parent != seg_equiv2) {
	  pair_seg_equivs[i_match]->usable = 0;
	  pair_seg_equivs[i_match] = 0;
	  break;
	}
      }
    }
  }    

  for (i_seg_equiv = n = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    if (!seg_equiv->usable) continue; 
    seg_equiv->n_segments = 1; /* reset to 1 */
    if (seg_equiv == seg_equiv->parent) {
      seg_equiv->id = n;
      n++;
    }
  }

  for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    if (!seg_equiv->usable) continue; 
    /* should also exclude lengths of filtered pairs */
    if (seg_equiv != seg_equiv->parent) {
      seg_equiv->parent->n_segments += 1;
      seg_equiv->parent->length += seg_equiv->length;
    }
  }

  for (i_seg_equiv = 0; i_seg_equiv < n_seg_equiv; i_seg_equiv++) {
    seg_equiv = seg_equivs + i_seg_equiv;
    if (!seg_equiv->usable) continue; 
    if (!seg_equiv->parent->usable) {
      printf("\nWARNING: deleted parent --  restored ");
      seg_equiv->parent = seg_equiv;
    }
  }

  free_seg_blocks();

  return short_flag;
}
