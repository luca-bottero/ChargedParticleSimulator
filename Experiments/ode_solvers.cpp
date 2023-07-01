# include "ode_solvers.h"
# include <math.h>

# define C 299792458.
# define e 1.602176634e-19
# define m_e 9.1093837015e-31
# define m_p 1.67262192369e-27

void RK4Step(double t, double alpha, double dt, double * Y, int Neq, void (*RHS_Func)(double, double, double *, double *)){
    double Y1[Neq], k1[Neq], k2[Neq], k3[Neq], k4[Neq];

    RHS_Func(t, alpha,Y,k1);

    for (int i=0; i < Neq; i++) Y1[i] = Y[i] + 0.5 * dt * k1[i];

    RHS_Func(t + dt*0.5, alpha, Y1, k2);

    for (int i=0; i < Neq; i++) Y1[i] = Y[i] + 0.5 * dt * k2[i];

    RHS_Func(t + dt*0.5, alpha, Y1, k3);

    for (int i=0; i < Neq; i++) Y1[i] = Y[i] + dt * k3[i];

    RHS_Func(t + dt, alpha, Y1, k4);

    for (int i = 0; i < Neq; i++) Y[i] += dt * (k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]) / 6.;
}

void Boris2ndOrderStep(double t, double dt, double alpha, double * Y, int Neq, 
                        void (*E_Field)(double, double *, double *), void (*B_Field)(double, double *, double *)){
    // Calculates one iteration of the Boris-push algorithm of the 2nd order for a signle particle

    double x_half[3], u_minus[3], u_plus[3];
    double E[3], b[3];          // electrical field at n+1/2; scaled magnetic field
    
    double gamma = 1./sqrt(1. - (Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])/(C*C));     // Lorentz factor

    // drift
    for(int i=0; i<3; i++) x_half[i] = Y[i] + 0.5*dt*Y[i+3];   

    // Get field value at x_half
    double pos[3] = {x_half[0], x_half[1], x_half[2]};
    E_Field(t + 0.5*dt, pos, E);
    B_Field(t + 0.5*dt, pos, b);

    // kick
    for(int i=0; i<3; i++) u_minus[i] = Y[i+3]*gamma + 0.5*alpha*dt*E[i]; 

    double gamma_minus = sqrt(1. + (u_minus[0]*u_minus[0] + u_minus[1]*u_minus[1] + u_minus[2]*u_minus[2])/(C*C));

    // Get the b field which is defined as h*B/(2*gamma)
    for(int i=0; i<3; i++) b[i] *= alpha*dt*0.5/gamma_minus;

    // rotate
    double u_minus_cross_b[3], fraction[3], rotation_cross[3];
    double correction_factor = 2. / (1. + b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    cross_prod(u_minus, b, u_minus_cross_b);
    for(int i=0; i<3; i++) fraction[i] = (u_minus[i] + u_minus_cross_b[i]) * correction_factor;
    cross_prod(fraction, b, rotation_cross);
    for(int i=0; i<3; i++) u_plus[i] = u_minus[i] + rotation_cross[i];

    //gamma = 1./sqrt(1. - (Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])/(C*C)); 

    // kick
    for(int i=0; i<3; i++) Y[i+3] = u_plus[i] + 0.5*alpha*dt*E[i];

    double gamma_plus = sqrt(1. + (Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])/(C*C)); 

    // Transform the relativistic corrected velocity back to velocity
    for(int i=0; i<3; i++) Y[i+3] /= gamma_plus;
    // drift
    for(int i=0; i<3; i++) Y[i] = x_half[i] + 0.5*dt*Y[i+3];   
}

void cross_prod(double * a, double * b, double * res){
    // computes res = a x b

    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
}










