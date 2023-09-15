main(argc,argv)
     int argc;
     char *argv[];
{
  int i, c;
  char *signed_c;
  unsigned char *unsigned_c;
  char array[256];
  unsigned char unsigned_array[256];

  for (i = 0; i < 256; i++) {
    array[i] = i;
    unsigned_array[i] = i;
    printf("\n%d %d %d", i, array[i], unsigned_array[i]); 
  }
}
