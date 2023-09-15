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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct entry {
  long int old_loc, new_loc, new_ptr_loc, ref_loc;
  char marked;
} Entry;

typedef struct entry_block {
  struct entry_block *next, *prev;
  Entry *first_entry, *last_entry; /* actually last + 1 */
  int n_marked;
} Entry_block;

typedef struct read_stat {
  struct read_stat *next, *prev;
  long int header; /* header location in file */
  long int ref_loc; /* ref seq location of read start */
  int strand, map_qual, contin;
} Read_stat;

typedef struct file_stat {
  char name[101];
  struct file_stat *next;
  Read_stat *read_stat;
  FILE *fp;
  Entry_block *head_block, *tail_block; /* 1st, last block of entries in index */
  int t_n_marked; /* # marked entries in inbox */
  long int ref_loc; /* ref seq location */
  int ref_num; /* number of ref seq currently in */
  int type, prev_type, ref_nuc, next_ref_nuc, half; /* record type, ref_nuc, next_ref_nuc 
						       half = 0 for upper 4 bits, 1 for lower 4 bits (in type 2 records) */
  int n_active; /* # active reads */
  int ascii_only;
} File_stat;

