#include <math.h>

main()
{
  int i;
  double n;
  double x, y, a_n, b_n, prev_b_n, log2pi, f_n;

  log2pi = log(2 * 3.14159265);
  for (i = 1, n = 10; i <= 5; i++, n *= 10) {
    b_n = 1.0;
    do {
      prev_b_n = b_n;
      b_n = sqrt(2 * log(n / prev_b_n) - log2pi);
    } while (fabs(prev_b_n - b_n) > .0001);
    a_n = 1.0 / b_n;
    for (x = -3.0; x <= 5.0; x += .5) {
      y = a_n * x + b_n;
      f_n = pow(.5 + .5 * erf(y / sqrt(2.0)), n);
      printf("\n%.0f %.1f %.6f %.6f", n, x, f_n, f_n - exp(-exp(-x)) );
    }
  }
}

normal_cdf_residual(x)
     double x;
{
  return .5 * erfc(x / sqrt(2.0));
}
