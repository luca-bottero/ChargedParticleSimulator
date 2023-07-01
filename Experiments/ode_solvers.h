void EulerStep  (double, double, double, double *, int, void (*F)(double, double, double *, double *));
void RK2Step    (double, double, double, double *, int, void (*F)(double, double, double *, double *));
void RK4Step    (double, double *, double, double *, int, void (*F)(double, double, double *, double *));
void Boris2ndOrderStep  (double, double, double, double *, int, void (*F)(double, double *, double *), void (*G)(double, double *, double *));
void cross_prod (double *, double *, double *);
