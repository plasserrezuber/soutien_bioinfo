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
extern Database *query_db;
extern int t_num_entries;

typedef struct cluster {
  int entry, counts;
} Cluster;

main(argc,argv)
     int argc;
     char *argv[];
{
  int i, n, entry1, i_cluster, low, high, mid, final;
  int num_clusters;
  char *our_alloc();
  int class, class1, entry, final_class;
  Cluster *clusters;
  int histogram[50001];
  char *get_id(), *get_descrip();
  unsigned char *get_seq();
  FILE *fp;
  FILE *fopenWrite();
  int max_size, max_entry;

  get_parameters(argc, argv, "cluster");
  readin_and_match();

  num_clusters = 0;
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    class1 = get_equiv_class(entry1);
    for (class = class1, entry = entry1; class != entry;
	 entry = class, class = get_equiv_class(entry));
    if (entry1 == class) num_clusters++;
    else {
      final_class = class;
      for (class = class1, entry = entry1; class != final_class; 
	   entry = class, class = get_equiv_class(entry))
	set_equiv_class(entry, final_class);
    }
  }

  printf("\nNo. clusters: %d\n", num_clusters);
  clusters = (Cluster *)our_alloc(num_clusters * sizeof(Cluster));
  for (entry1 = i_cluster = 0; entry1 < t_num_entries; entry1++) {
    if ((class = get_equiv_class(entry1)) == entry1) {
      clusters[i_cluster].entry = entry1;
      clusters[i_cluster].counts = 1;
      i_cluster++;
    }
    else {
      low = 0;
      high = i_cluster - 1;
      if (clusters[low].entry == class) {
	final = low;
      }
      else if (clusters[high].entry == class) {
	final = high;
      }
      else {
	final = -1;
	do {
	  mid = (low + high) / 2;
	  if (class < clusters[mid].entry) high = mid;
	  else if (class > clusters[mid].entry) low = mid;
	  else {
	    final = mid;
	  }
	} while (final < 0);
      }
      clusters[final].counts += 1;
    }
  }
  for (i = 0; i < 50001; i++) histogram[i] = 0;
  max_size = 0;
  for (i_cluster = 0; i_cluster < num_clusters; i_cluster++) {
    i = clusters[i_cluster].counts;
    if (max_size < i) {
      max_size = i;
      max_entry = clusters[i_cluster].entry;
    }
    histogram[i < 50000 ? i : 50000] += 1;
    if (i > 100)
      printf("\n%5d  %s", i, get_id(clusters[i_cluster].entry));
  }
  printf("\n\nCluster size  #clusters\n");

  for (i = n = 0; i < 50001; i++)
    if (histogram[i]) {
      printf("\n%4d     %5d", i, histogram[i]);
      n += i * histogram[i];
    }

  fp = fopenWrite("maxcluster");
  for (entry1 = 0; entry1 < t_num_entries; entry1++) {
    if (max_entry == get_equiv_class(entry1)) {
      write_entry(fp, get_id(entry1), get_descrip(entry1), get_seq(entry1));
    }
  }

  if (n != t_num_entries) fatalError("histogram counts");
  printf("\n");
  notify("\n");
  return 0;
}
