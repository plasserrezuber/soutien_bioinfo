/*****************************************************************************
#   Copyright (C) 1993-1996 by Phil Green.                          
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

#include <stdio.h>
/*   reads file in .seq format, and creates a new file in FASTA format */

main(argc, argv)
     int argc;
     char *argv[];
{
  FILE *fp, *fq;
  FILE *fopenRead(), *fopenWrite();
  int i_file, i, i_start, i_middle, i_end, i_successes;
  char c;
  int seq_len, start_trim, keep_len;
  char read_buffer[1000], seq_buffer[10000];

  i_successes = 0;
  fq = fopenWrite("temp.fasta");
  for (i_file = 1; i_file < argc; i_file++) {
    fp = fopenRead(argv[i_file]);
    fgets(read_buffer, 1000, fp);
    sscanf(read_buffer, "; %d %d %dABI", &seq_len, &start_trim, &keep_len);
    i_start = 0;
    i_middle = 0; /* start_trim; */
    i_end = keep_len; /* start_trim + keep_len; */
    while (fgets(read_buffer, 1000, fp)) {
      if (read_buffer[1] == '<') {
	;
/* ignore initial trimmed sequence -- is sequencing vector in general 
	for (i = 2; c = read_buffer[i]; i++) {
	  if (isalpha(c)) seq_buffer[i_start++] = tolower(c);
	  else if (c == '-') seq_buffer[i_start++] = 'N';
	}
*/
      }
      else if (read_buffer[1] == '>') {
	for (i = 2; c = read_buffer[i]; i++) {
	  if (isalpha(c)) seq_buffer[i_end++] = tolower(c);
	  else if (c == '-') seq_buffer[i_end++] = 'N';
	}
      }
      else for (i = 0; c = read_buffer[i]; i++) {
	  if (isalpha(c)) seq_buffer[i_middle++] = c;
	  else if (c == '-') seq_buffer[i_middle++] = 'N';
	}
     }
/*
    if (i_start != start_trim || i_middle != start_trim + keep_len
	|| i_end != seq_len)
*/
    if (i_middle != keep_len || i_end != seq_len - start_trim)
      printf("\nERROR: %d vs. %d, %d vs. %d, %d vs. %d", i_start, start_trim, 
	     i_middle, start_trim + keep_len,
	     i_end, seq_len);
    else i_successes++;
    fprintf(fq,">%s",argv[i_file]);
    for (i = 0; i < i_end; i++) {
      if (!(i%50)) fprintf(fq,"\n");
      fprintf(fq,"%c",seq_buffer[i]);
    }
    fprintf(fq,"\n");

    if (!(i_successes % 50)) printf("%d\n", i_successes);
    fclose(fp);
  }
  printf("%d files out of %d successfully read", i_successes, i_file - 1);
  fclose(fq);
}


