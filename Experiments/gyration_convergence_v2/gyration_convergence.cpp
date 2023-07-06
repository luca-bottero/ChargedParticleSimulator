# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include "../ode_solvers.h"

# define C 299792458.
# define e 1.602176634e-19
# define m_e 9.1093837015e-31
# define m_p 1.67262192369e-27

# define N_MAX 1e5
# define N_ROUND_MIN 1
# define N_ROUND_MAX 100

using namespace std;

void Set_Initial_Values (double *);
void RHS_Func   (double, double *, double *, double *);
void E_Field (double, double *, double *);
void B_Field (double, double *, double *);

int main(){
    
    double t;
    double t_in = 0., t_fin;
    int n;

    double dt;
    double T;

    int N_PARTICLES = 1;
    int Neq = 6 * N_PARTICLES;

    double Y[Neq];

    double alpha = -e/(m_e);

    // Result correctness

    ofstream fdata;
    fdata.precision(20);
    fdata.open("radius_error.dat");

    // gamma = 10
    Set_Initial_Values(Y);
    Y[4] = sqrt(1. - 0.01) * C;
    double gamma = 1./sqrt(1. - Y[4]*Y[4]/(C*C));
    Y[0] = gamma*m_e*Y[4]/e;            // set x pos to radius value   

    double x_start = Y[0];     

    cout << Y[4] << " " << gamma << endl;

    t_fin = 2*M_PI*m_e/(e)*gamma;

    for (int i=1; i < 41; i++){
        dt = (t_fin - t_in)/(double) 100;

        t = t_in;

        for (int j = 0; j < 100*i; j++){
            Boris2ndOrderStep(t, dt, alpha, Y, Neq, E_Field, B_Field);
            // cout << j << endl;
            t += dt;
        }
        cout << i << " " << dt << " " << Y[0] << " " << Y[1] << " " << abs(sqrt(Y[0]*Y[0] + Y[1]*Y[1]) - x_start);
        cout << " " << atan2(Y[1],Y[0]) << endl;
        cout << Y[0] << " " << Y[1] << " " << t/t_fin << endl;

        // i  dt   r-r_0   phi   dist
        fdata << i << " " << dt << " " << abs(sqrt(Y[0]*Y[0] + Y[1]*Y[1]) - x_start);
        fdata << " " << abs(atan2(Y[1], Y[0])) << " " << sqrt(Y[1]*Y[1] + (Y[0] - x_start)*(Y[0] - x_start)) << endl;        
    }
    fdata << endl << endl;
    cout << endl;

    // gamma = 10^4
    Set_Initial_Values(Y);
    Y[4] = sqrt(1. - 1./1e8) * C;

    gamma = 1./sqrt(1. - Y[4]*Y[4]/(C*C));
    Y[0] = gamma*m_e*Y[4]/e;            // set x pos to radius value   

    x_start = Y[0];     

    t_fin = 2*M_PI*m_e/(e)*gamma;

    for (int i=1; i < 41; i++){
        dt = (t_fin - t_in)/(double) 100;

        t = t_in;

        for (int j = 0; j < 100*i; j++){
            Boris2ndOrderStep(t, dt, alpha, Y, Neq, E_Field, B_Field);
            t += dt;
        }
        cout << i << " " << dt << " " << Y[0] << " " << Y[1] << " " << abs(sqrt(Y[0]*Y[0] + Y[1]*Y[1]) - x_start);
        cout << " " << atan2(Y[1],Y[0]) << endl;

        fdata << i << " " << dt << " " << abs(sqrt(Y[0]*Y[0] + Y[1]*Y[1]) - x_start);
        fdata << " " << abs(atan2(Y[1], Y[0])) << " " << sqrt(Y[1]*Y[1] + (Y[0] - x_start)*(Y[0] - x_start)) << endl;
    }
}


void Set_Initial_Values(double * Y){
    Y[0] = 0. *C;  // x
    Y[1] = 0. *C;  // y
    Y[2] = 0. *C;  // z
    Y[3] = 0. *C;  // v_x
    Y[4] = 0. *C;  // v_y
    Y[5] = 0. *C;  // v_z
}

void E_Field(double t, double * pos, double * E){
    E[0] = 0.;
    E[1] = 0.;
    E[2] = 0.;
}

void B_Field(double t, double * pos, double * B){
    B[0] = 0.;
    B[1] = 0.;
    B[2] = 1.;
}

void RHS_Func(double t, double * alpha, double * Y, double * rhs){
    // This is the RHS of the equation for a single particle with charge-mass ratio = alpha

    double E[3], B[3];
    double pos[3] = {Y[0], Y[1], Y[2]};
    E_Field(t, pos, E);
    B_Field(t, pos, B);

    double vec_prod[3];
    double vel[3] = {Y[3], Y[4], Y[5]};
    cross_prod(vel, B, vec_prod);

    for(int i=0; i<3; i++) rhs[i] = Y[i+3];          
    for(int i=0; i<3; i++) Y[i+3] = alpha[0] * (C * E[i]  + vec_prod[i]);    // velocity
    
    double gamma = 1. / sqrt(1. - (Y[3]*Y[3] + Y[4]*Y[4] + Y[5]*Y[5])/(C*C));

    for(int i=0; i<3; i++) rhs[i+3] = Y[i+3]/gamma;
}
