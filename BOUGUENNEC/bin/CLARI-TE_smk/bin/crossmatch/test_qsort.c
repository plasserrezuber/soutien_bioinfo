compare_doubles(i,j)
     double *i, *j;
{
  return *i - *j;
}

main()
{
  char *our_alloc();
  int i, n;
  double *x;
  double y, z;
    
  for (i = 0; i < 10; i++) {
    n = random();
    y = n % 100000;
    z = y / 10.0;
    printf("%d %f %f\n ", n, y, z);
  }
  n = 1000000;
  x = (double *)our_alloc(n * sizeof(double));
  for (i = 0; i < n; i++) {
    x[i] = (random() % (10 * n)) / 10.0;
    if (i < 10) printf("%f\n ", x[i]);
  }
  qsort(x, n, sizeof(double), compare_doubles);
  for (i = 0; i < n; i += n / 10)
    printf("%6d  %6.3f\n", i, x[i]);
}
