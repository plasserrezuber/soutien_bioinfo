/*****************************************************************************
#   Copyright (C) 2008 by Phil Green.                          
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

#include "calf.h"

File_stat *head_file_stat;

extern void *malloc(), *realloc();

append_file_stat(file_name)
     char *file_name;
{
  FILE *fopen();
  File_stat *file_stat;
  static int n_files;

  n_files++;
  if (n_files >= FOPEN_MAX) { /* assuming one output file */
    fprintf(stderr, "FATAL ERROR: at most %d total files (input + output) can be processed",
	    FOPEN_MAX);
    exit(1);
  }
  
  file_stat = (File_stat *)malloc(sizeof(File_stat));

  file_stat->fp = fopen(file_name, "r");
  if (!file_stat->fp) {
    fprintf(stderr, "\n ERROR: FILE %s NOT FOUND \n", file_name);
    exit(1);
  }
  strncpy(file_stat->name, file_name, 100);

  file_stat->ref_loc = file_stat->ref_num = file_stat->type = file_stat->prev_type = file_stat->ref_nuc = file_stat->next_ref_nuc = file_stat->half = file_stat->n_active = file_stat->t_n_marked = file_stat->ascii_only = 0;
  file_stat->read_stat = 0;
  file_stat->head_block = file_stat->tail_block = 0;
  file_stat->next = head_file_stat;
  head_file_stat = file_stat;
/* order of files in linked list is reversed from command line */
  for (file_stat = file_stat->next; file_stat; file_stat = file_stat->next)
    if (!strcmp(file_stat->name, file_name)) fatal_error("duplicate input file names");
}

/* get new record header byte; should extend, to allow return code */
get_record_header(file_stat)
     File_stat *file_stat;
{
  int c, type;

  c = fgetc(file_stat->fp);

  if (c == EOF) fatal_error("unexpected end of file");

  file_stat->type = type = c & 3;
  file_stat->prev_type = (c >> 2) & 3;
  file_stat->ref_nuc = c >> 4;

  if (c && !type) fatal_error("type 0 record header");

  if (file_stat->ref_nuc && type != 1) fatal_error("non-type 1 record header w/ non-zero nuc");

  if (type != 1 && file_stat->n_active) fatal_error("discordant start & end markers");

  if (type == 3) {
    file_stat->half = 0;
    read_type_3(file_stat, 1);
  }
  else if (type != 1 || file_stat->ref_nuc) 
    file_stat->ref_loc += 1;

  if (!file_stat->prev_type) {
    file_stat->ref_loc += 1; /* needed when starting new ref -- to distinguish gap cases */
    file_stat->ref_num += 1; 
  }
}

read_type_2(file_stat)
     File_stat *file_stat;
{
  int i;
  long int size;

 /* uncovered size-only segment -- read 4 bytes big-endian */

  for (i = size = 0; i < 4; i++) size = (size << 8) + fgetc(file_stat->fp);
  if (fgetc(file_stat->fp)) fatal_error("type 2 record terminator");
  if (!size) fatal_error("0 size type 2 record");
  file_stat->ref_loc += size - 1; /* last base -- since was already at first base */

  get_record_header(file_stat);
}

read_type_3(file_stat, new_rec_flag)
     File_stat *file_stat;
     int new_rec_flag;
{
  int c;

  if (!file_stat->half) {
    c = fgetc(file_stat->fp);
    file_stat->ref_nuc = c >> 4;
    file_stat->next_ref_nuc = c & 15;
  }
  else {
    file_stat->ref_nuc = file_stat->next_ref_nuc;
  }

  if (!file_stat->ref_nuc) {
    if (file_stat->half) c = fgetc(file_stat->fp);
    if (c) fatal_error("corrupt type 3 record"); /* expect record terminator */
    if (new_rec_flag) fatal_error("empty type 3 record"); 

    get_record_header(file_stat);
    return;
  }
  file_stat->ref_loc += 1;

  file_stat->half = !file_stat->half;
}

/* 
   use fact that ftell reports origin 0 co-ord of next byte to be written or read
      fseek(fp_out, n, SEEK_SET) positions to write or read n-th byte 
      note reading or writing advances position by 1

   if fp_out == 0, then is read only
*/

read_write_type_1(file_stat, fp_out, ref_offset_flag, header_flag)
     File_stat *file_stat;
     FILE *fp_out;
     int ref_offset_flag, header_flag;
{
  FILE *fp;
  int c, n, new_n, new_q, ref_n, q, needs_ref;
  int start, first, sign, ref_sign, bases_since_start;
  long int byte_offset, new_byte_offset, ref_offset, new_ref_offset, input_start_loc, output_start_loc, output_ptr_loc, fp_loc;
  Entry *lookup(), *entry;
  long int read_offset();
  long int max_2, max_4, max_off;
  Read_stat *read_stat, *new_read_stat, *last_read_stat, *prev_read_stat;
  
  fp = file_stat->fp;
  read_stat = file_stat->read_stat; /* next read_stat to be utilized */
  last_read_stat = 0; /* previous active one */
  bases_since_start = 0; /* # bases since last start marker; needed to identify 'tags' (0 length
			    reads) by checking end marker byte */
  while (c = fgetc(fp)) {
    if (fp_out) fputc(c, fp_out); 
    if (c == 192) { /* '*' byte, indicating unaligned region; find terminator */
      do {
	c = fgetc(fp);
	if (fp_out) fputc(c, fp_out);
      } while (c != 192);
      continue;
    }
    q = c & 63;
    if (q == 63 && !n) { /* end marker byte */
      file_stat->n_active -= 1;
      if (file_stat->n_active < 0) fatal_error("discordant start & end markers");

      if (!bases_since_start) { /* this is a tag (0 length read) so remove current read_stat rather
				 than prev one */
	last_read_stat = read_stat;
	read_stat = read_stat->next;
      }
      /* remove last_read_stat from list */
      if (!last_read_stat) fatal_error("0 read_stat");
      prev_read_stat = last_read_stat->prev;
      if (prev_read_stat) prev_read_stat->next = read_stat;
      else file_stat->read_stat = read_stat;
      if (read_stat) read_stat->prev = prev_read_stat;
      realloc(last_read_stat, 0);
      last_read_stat = prev_read_stat;
    }
    else if (q > 61) { /* start marker byte */
      new_n = n = c >> 6;
      new_q = n && ref_offset_flag == 2 ? 63 : !ref_offset_flag ? 62 : q;
      file_stat->n_active += 1;
      start = c;
      input_start_loc = ftell(fp) - 1; 
      if (fp_out) {
	output_start_loc = ftell(fp_out) - 1;
      }
      c = fgetc(fp);
      if (!c) { /* 0 byte starting ascii header */
	do { 
	  if (header_flag) fputc(c, fp_out);  
	  c = fgetc(fp);
	} while (c); /* find next 0 byte, which ends header */
	if (header_flag) fputc(c, fp_out);  
	c = fgetc(fp);
      }  /* now at strand/mapqual byte */
      if (fp_out) fputc(c, fp_out);  

      /* insert a new read_stat into the list */
      new_read_stat = (Read_stat *)malloc(sizeof(Read_stat));
      new_read_stat->next = read_stat;
      new_read_stat->prev = last_read_stat;
      if (read_stat) read_stat->prev = new_read_stat;
      if (last_read_stat) last_read_stat->next = new_read_stat;
      else file_stat->read_stat = new_read_stat;

      read_stat = new_read_stat;

      read_stat->strand = c >> 7;
      read_stat->map_qual = c & 127;
      read_stat->contin = -1; /* default value: no continuation */
      read_stat->ref_loc = file_stat->ref_loc;

      bases_since_start = 0; /* reset at new start marker */

      if (n > 0) { /* there is a byte (and possibly refseq) offset */
	first = fgetc(fp);
	sign = (first >> 3) & 1; /* sign bit for offset */
	read_stat->contin = first >> 3; /* retain sign bit, plus others */

	new_byte_offset = byte_offset = read_offset(fp, n, (long int)(first & 7));

	ref_offset = 0;

	if (q == 63) { /* get refseq offset as well */
	  c = fgetc(fp);
	  ref_sign = c >> 7;
	  ref_offset = read_offset(fp, n, (long int)(c & 127));
	  if (ref_sign != sign && ref_offset) fatal_error("refseq offset sign"); /* signs should agree, unless refseq offset is 0 */
	}
	new_ref_offset = ref_offset;
	if (fgetc(fp) != start) fatal_error("start marker copy");
      
	if (fp_out) {
	  if (byte_offset) {
	    output_ptr_loc = ftell(fp_out); /* 1st ptr byte */
	    if (!sign) { /* points downstream */
	      new_n = 3; /* increase allocation to 6 bytes */
	      append_entry(file_stat, input_start_loc, output_start_loc, output_ptr_loc);
	    }
	    else { /* points upstream */
	      /* lookup location of start in index */
	      entry = lookup(file_stat, input_start_loc - byte_offset);
	      new_byte_offset = output_start_loc - entry->new_loc;
	      new_ref_offset = file_stat->ref_loc - entry->ref_loc;
	      if (q == 63 && ref_offset != -new_ref_offset) /* assuming ref offset symmetry */
		fatal_error("reference offset discrepancy");
	      needs_ref = q != new_q && new_q == 63; /* assuming ref offset symmetry */
	      if (new_byte_offset != byte_offset || needs_ref) {
		/* alter pointer bytes at prev start -- change low order 3 bits + subsequent 5 bytes */
		fseek(fp_out, entry->new_ptr_loc - output_ptr_loc, SEEK_END); 
		c = fgetc(fp_out);  
		fseek(fp_out, (long)(-1), SEEK_CUR); /* back up, to rewrite */
		write_offset(fp_out, 3, c & 248, new_byte_offset); /* keep 1st five bits of 1st byte */ 
		if (needs_ref)
		  write_offset(fp_out, 3, 0, new_ref_offset); 
		fseek(fp_out, (long)0, SEEK_END);
	      }
	      /* find min # bytes needed for current pointer */

	      max_2 = (8 << 8) - 1; /* max size of pointer storable in 2 bytes (less sign &
				       other special bits */
	      max_4 = (8 << 24) - 1;

	      max_off = new_byte_offset > new_ref_offset >> 4 ? new_byte_offset : new_ref_offset >> 4; /* no special bits with ref_offset; it can exceed new_byte_offset (if there are type 2 or type 3 records */ 

	      new_n = max_off <= max_2 ? 1 : max_off <= max_4 ? 2 : 3;
	    }
	  } /* if byte_offset */

	  if (new_n > 0) {
	    if (new_n != n || new_q != q) { /* must revise start marker */
	      fseek(fp_out, output_start_loc, SEEK_SET);
	      start = (new_n << 6) + new_q;
	      fputc(start, fp_out); 
	      fseek(fp_out, (long)0, SEEK_END);
	    }
	    /* write offset pointers */
	    write_offset(fp_out, new_n, first & 248, new_byte_offset); /* keep 1st five bits of 1st byte (which include sign) */
	    if (new_q == 63) 
	      write_offset(fp_out, new_n, ref_sign << 7, new_ref_offset); 
	  }
	  fputc(start, fp_out); 
	} /* if fp_out */
      } /* pointer case */
    } /* start marker */
    else { /* aligned base or gap char */
      if (!read_stat) fatal_error("0 read stat");
      bases_since_start++;
      last_read_stat = read_stat;
      /* any code to consider aligned base goes here */
      read_stat = read_stat->next; /* next one */
    }
  }
  /* end of record reached; get next header */
  get_record_header(file_stat);
}

/* NEED TWO'S COMPLEMENT ADJUSTMENT */

long int read_offset(fp, n, offset)
     FILE *fp;
     int n; 
     long int offset;
{
  int i;

  for (i = 1; i < 2 * n; i++)  
    offset = (offset << 8) + fgetc(fp);

  return offset;
}

/* NEED TWO'S COMPLEMENT ADJUSTMENT */

write_offset(fp, n, add, offset)
     FILE *fp;
     int n, add; /* add is added to first byte only */
     long int offset;
{
  int shift;

  for (shift = 16 * n - 8; shift >= 0; shift -= 8, add = 0)
    fputc(add + ((offset >> shift) & 255), fp);
}

/* strategy: 
   keep linked list of blocks of entries, append new entries to end; 
   lookup: find block by linear search starting at most recent, 
     then do binary search on it, mark when found, count # of found marks
   periodically delete all found
   always lookup last block first
*/

#define BLOCKSIZE 10000

append_entry(file_stat, old_loc, new_loc, new_ptr_loc)
     File_stat *file_stat;
     long int old_loc, new_loc, new_ptr_loc;
{
  Entry_block *block, *new_block;
  Entry *entry;

  block = file_stat->tail_block;
  if (!block || block->last_entry >= block->first_entry + BLOCKSIZE) {
    new_block = (Entry_block *)malloc(sizeof(Entry_block));
    new_block->last_entry = new_block->first_entry = (Entry *)malloc(BLOCKSIZE * sizeof(Entry));
    new_block->next = 0;
    new_block->n_marked = 0;

    if (block) block->next = new_block;
    else file_stat->head_block = new_block;

    new_block->prev = block;
    block = file_stat->tail_block = new_block;
  }

  entry = block->last_entry;
  entry->old_loc = old_loc;
  entry->new_loc = new_loc;
  entry->new_ptr_loc = new_ptr_loc;
  entry->ref_loc = file_stat->ref_loc;
  entry->marked = 0;
  block->last_entry += 1;
}

Entry *lookup(file_stat, old_loc)
     File_stat *file_stat;
     long int old_loc;
{
  Entry_block *block;
  Entry *mid_entry;
  int min, max, mid;

  for (block = file_stat->tail_block; block && block->first_entry->old_loc > old_loc; block = block->prev);
  if (!block) fatal_error("entry not found");
  min = 0;
  max = (block->last_entry - block->first_entry) - 1;
  do {
    mid = (min + max) / 2;
    mid_entry = block->first_entry + mid;
    if (mid_entry->old_loc >= old_loc) max = mid;
    else min = mid + 1;
  } while (min < max);
  if (mid_entry->old_loc != old_loc) fatal_error("entry not found");
  if (mid_entry->marked) fatal_error("entry found twice");
  mid_entry->marked = 1;
  block->n_marked += 1;
  file_stat->t_n_marked += 1;
  if (file_stat->t_n_marked >= 2 * BLOCKSIZE) cleanup(file_stat);
  return mid_entry;
}

cleanup(file_stat)
     File_stat *file_stat;
{
  Entry_block *old_block, *new_block;
  Entry *old_entry, *new_entry;

  new_block = file_stat->head_block;
  new_entry = new_block->first_entry;
  file_stat->t_n_marked = 0; /* anticipating */
  /* could make following faster by skipping blocks for which n_marked = 0 */
  for (old_block = file_stat->head_block; old_block; old_block = old_block->next) {
    for (old_entry = old_block->first_entry; old_entry < old_block->last_entry; old_entry++) {
      if (old_entry->marked) continue;
      if (new_entry >= new_block->first_entry + BLOCKSIZE) {
	new_block->last_entry = new_block->first_entry + BLOCKSIZE;
	new_block->n_marked = 0;
	new_block = new_block->next;
	new_entry = new_block->first_entry;
      }
      if (new_entry != old_entry) {
	new_entry->old_loc = old_entry->old_loc;
	new_entry->new_loc = old_entry->new_loc;
	new_entry->ref_loc = old_entry->ref_loc;
	new_entry->marked = 0;
      }
      new_entry++;
    }
  }
  new_block->last_entry = new_entry; 
  new_block->n_marked = 0;
  for (old_block = new_block->next; old_block; old_block = old_block->next) {
    realloc(old_block->first_entry, 0);
    realloc(old_block, 0);
  }
  new_block->next = 0;
  file_stat->tail_block = new_block;
}

/* at end should be no unmarked entries!! -- check for this */
 
fatal_error(message)
     char *message;
{
  fprintf(stderr, "\nFATAL ERROR: %s\n", message);
  exit(1);
}

