/*****************************************************************************
#   Copyright (C) 1997 by Phil Green.                          
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
extern int t_num_entries;

main(argc,argv)
     int argc;
     char *argv[];
{
  int entry1;
  char file_name[100];
  char *get_id(), *get_descrip();
  unsigned char *get_seq();
  FILE *fq;
  FILE *fopenWrite();
  unsigned char *seq;
  Database *query_db;
  Database *append_db();
  Profile *maximal_profile_segments(), *q_profile;

  get_parameters(argc, argv, "loco");
  query_db = append_db(parameters->query_files, 1, 0, 1); 
  printf("\nLow complexity regions (maximal single base matches):");

  if (parameters->screen) {
    strcpy(file_name, query_db->file->name);
    strcat(file_name, ".loco");
    fq = fopenWrite(file_name);
  }

  q_profile = 0;
  for (entry1 = query_db->first_entry; entry1 <= query_db->last_entry; entry1++) {
    seq = get_seq(entry1);
    q_profile = maximal_profile_segments(q_profile, seq, get_seq_length(entry1), parameters->minscore, get_id(entry1), parameters->screen);
    if (parameters->screen)
      write_entry(fq, get_id(entry1), get_descrip(entry1), seq);
  }
}
