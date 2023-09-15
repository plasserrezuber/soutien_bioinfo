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

extern Parameters *parameters; /* from get_parameters() */
extern Database *query_db;

main(argc,argv)
     int argc;
     char *argv[];
{

  /*
  fprintf(stderr," %d %d %d ", 15 >> 3 << 3, (15 >> 3) << 3, 15 >> (3 << 3));
  exit(1);
  */
  /* 
  int i, j, k;
  double x, y, z;

  for (i = 1; i < 12; i++) {
    x = pow(2.0, -i / 12.0);
    printf("\n%2d %.4f  ", i, x);
    z = 3.0;
    for (j = 2; j < 21; j++)
      for (k = j + 1; k < 21; k++) {
	y = k / (double)j; 
	if (z > 261 * fabs(y * x - 1)) { 
	  z = 261 * fabs(y * x - 1);
	  printf(" %d/%d: %.4f", k, j, z);
	}
      }
  }
  printf("\n");
  exit(1);
  /* */

#if defined(GCPHRAP)
  get_parameters(argc, argv, "gccross_match");
#else
  get_parameters(argc, argv, "cross_match");
#endif

/* read in sequence files, initialize array of aligned read info, and log file of ancillary information,
 and find pairwise matches */
  if (parameters->score_flag)
    set_domain_vars();

  readin_and_match();

  print_matches();
/*  print_coverage(query_db); */
/*  analyze_discreps(query_db); Needs revision*/

  if (parameters->screen) screen_seqs(query_db);

  check_pairs();
  check_alloc_list();
  print_t_n_pairs();
  print_t_n_segments();
  print_t_n_diffs();
  notify("\n");
  printf("\n");
  set_time("end");
  print_times();
  return 0;
}

