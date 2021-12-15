
#ifdef __cplusplus
extern "C" {
#endif
double brent(double ax, double bx, double cx, double (*f)(double), double tol,  double *xmin);

#ifdef __cplusplus
}
#endif
