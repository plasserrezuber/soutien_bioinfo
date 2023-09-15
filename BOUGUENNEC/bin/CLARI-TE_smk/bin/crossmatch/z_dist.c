#define TOL 0.00000001
#include <math.h>

main(argc, argv)
     int argc;
     char *argv[];
{
  float z_max, z_incr;
  double e_last, z, e, e_diff, cum_p, mean, var;

  sscanf(argv[1],"%f",&z_max);
  sscanf(argv[2],"%f",&z_incr);

  e_last = 0;
  cum_p = mean = var = 0.0;
  for (z = -z_max; z <= z_max; z += z_incr) {
    
    e = exp(-exp(-1.282 * (z + z_incr / 2.0) - .577));
    e_diff = e - e_last;
    cum_p += e_diff;
    mean += e_diff * z;
/*     if (fabs(mean) < TOL) mean = 0.0; */
    var += e_diff * z * z;
    e_last = e;
  }
  printf("\n%.6f %.6f %.6f\n", cum_p, mean, var - mean * mean);
}
 /* above choice of parameters for Gumbel give mean of .000168, variance of
1.000857 */
