/*****************************************************************************
#   Copyright (C) 1999 by Phil Green.                          
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

call_write_exp_files(n_sing, n_contigs, contig_ptrs)
     int n_sing, n_contigs;
     Contig **contig_ptrs;
{
#if defined(GCPHRAP)
  write_exp_files(n_sing, n_contigs, contig_ptrs);
  return;
#else
  return;
#endif
}

Database *call_append_db(file, in_memory, store_comps, store_aligns)
     File *file;
     int in_memory, store_comps, store_aligns;
{
  Database *append_db(), *append_exp_db();

#ifdef GCPHRAP
  if (parameters->exp_input)
    return append_exp_db(file, in_memory, store_comps, store_aligns);
#endif

  return append_db(file, in_memory, store_comps, store_aligns);
}
