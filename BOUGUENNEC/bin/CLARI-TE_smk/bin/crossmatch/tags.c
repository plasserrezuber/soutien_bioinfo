/*****************************************************************************
#   Copyright (C) 1997-2007 by Phil Green.                          
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
static int total_n_tags;

/* checks whether segment from start to end is contained (to within tag_slop)
   a tag of the stated type
*/

int contained_in_tag(entry, type, start, end, complement, tag_slop)
     int entry, start, end, complement, tag_slop;
     char *type;
{
  Tag *tag;
  int length, temp;
  Align_info *align_entry;
  Align_info *get_align_entry();

  align_entry = get_align_entry(entry);
  if (!align_entry->tags) return 0;

  if (complement) {
    length = get_seq_length(entry);
    temp = start;
    start = length + 1 - end;
    end = length + 1 - temp;
  }
  start += tag_slop;
  end -= tag_slop;

  for (tag = align_entry->tags; tag; tag = tag->next) {
    if (strcmp(tag->type, type)) continue;
    if (start >= tag->start && end <= tag->end) return 1;
  }
  return 0;
}

int contained_in_tag_list(tag, start, end, tag_slop)
     Tag *tag;
     int start, end, tag_slop;
{

  if (!tag) return 0;

  start += tag_slop;
  end -= tag_slop;

  for (; tag; tag = tag->next) {
    if (start >= tag->start && end <= tag->end) return 1;
  }
  return 0;
}

append_tag(align_entry, type, start, end)
     Align_info *align_entry;
     int start, end;
     char *type;
{
  Tag *tag;
  char *our_alloc();

  tag = (Tag *)our_alloc(sizeof(Tag));
  total_n_tags++;
  tag->next = align_entry->tags;
  align_entry->tags = tag;
  tag->start = start;
  tag->end = end;
  strcpy(tag->type, type);
}

/* N.B. Following is for tags on reads only!! */
/* prints padded read positions -- computed in write_contigs from the original ones */

write_tags(fp, align_entry)
     FILE *fp;
     Align_info *align_entry;
{
  Tag *tag;
  char *get_id();

  for (tag = align_entry->tags; tag; tag = tag->next) {
    if (tag->start < -1) 
      fprintf(fp, "\n\nRT{\n%s %s phrap %d %d %s\n}", 
	      get_id(align_entry->seq_entry), tag->type, -1 - tag->start, -1 - tag->end,
	      parameters->date);
    
    else
      fprintf(fp, "\n\nWR{\n%s %s phrap %s\n}", 
	      get_id(align_entry->seq_entry), tag->type, parameters->date);
    
  }
}

count_tags(align_entry, whole_read_infos_ptr, read_tags_ptr)
     Align_info *align_entry;
     int *whole_read_infos_ptr;
     int *read_tags_ptr;
{
  Tag *tag;

  *whole_read_infos_ptr = 0;
  *read_tags_ptr = 0;

  for (tag = align_entry->tags; tag; tag = tag->next) {
    if (tag->start >= 0)
      ++(*read_tags_ptr);
    else
      ++(*whole_read_infos_ptr);
    
  }
}

print_t_n_tags()
{
  fprintf(stderr, "Total # tags: %d, size: %.3f Mbytes\n",
	  total_n_tags, total_n_tags * sizeof(Tag) / 1000000.);
}
