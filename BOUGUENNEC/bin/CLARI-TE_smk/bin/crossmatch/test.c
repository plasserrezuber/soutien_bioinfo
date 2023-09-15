#include <stdio.h>
main()
{
  double x, y;
  int i;

  printf("%d %d\n",sizeof(int),sizeof(int *));
  x = 0;
  fprintf(stdout,"line2line1");
  for (i = 0; i < 8000; i++) printf("l");
/*  fflush(stdout); */
  y = 1/x;
}
