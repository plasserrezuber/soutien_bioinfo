#include <ctype.h>
#include <math.h>
#include <stdio.h>

main(argc,argv)
     int argc;
     char *argv[];
{
  int i, j, n_entries, length, fix_len;
  FILE *fp;
  char alph[50], alphsize;

  alphsize = 4;
  strcpy(alph,"ACGT");
  /*  strcpy(alph,"ACDEFGHIKLMNPQRSTVWY"); */
  if (argc < 2) {
    fprintf(stderr, "Specify no. of sequences\n");
    exit(1);
  }
  sscanf(argv[1], "%d", &n_entries);
  if (argc > 2) fix_len = 1;
  else fix_len = 0;
  fp = stdout; /*fopen("random.sdb", "w"); */
  for (i = 0; i < n_entries; i++) {
    fprintf(fp,">RandSeq%d",i);
    length = fix_len ? 300 : 100 + random() % 400; 
    for (j = 0; j < length; j++) {
      if (!(j % 50)) fprintf(fp,"\n");
      fprintf(fp, "%c", alph[random() % alphsize]);
    }
    fprintf(fp,"\n");  
  }
}
