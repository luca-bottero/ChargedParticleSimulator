# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include "../ode_solvers.h"

# define NATURAL 1
# define SI 2

# define UNIT_SYSTEM SI

# if UNIT_SYSTEM == NATURAL
    // Defined by c = h_slash = 1
    # define C 1.
    # define e 0.303//8.5424546e-2//sqrt(4.* M_PI * 0.0072973525628)       // the last term is the fine-structure constant
    # define m_e 0.51099895e6                        // about 0.5 MeV
    # define m_p 938.27208816e6
    # define TESLA 1.444027e3
    # define METER 1.97327e-7
# elif UNIT_SYSTEM == SI
    # define C 299792458.
    # define e 1.602176634e-19
    # define m_e 9.1093837015e-31
    # define m_p 1.67262192369e-27
    # define TESLA 1.
    # define METER 1.
#endif

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

    int N_PARTICLES = 1;
    int Neq = 6 * N_PARTICLES;

    double Y[Neq];

    double alpha = -e/(m_e);

    // Result correctness

    ofstream fdata;
    fdata.precision(20);
    fdata.open("output.dat");

    for(int i=1; i<16; i++){
        Set_Initial_Values(Y);
        Y[4] = (1. - 1./(pow(10,i))) * C;
        double gamma = 1./sqrt(1. - Y[4]*Y[4]/(C*C));
        Y[0] = gamma*m_e*Y[4]/e;            // set x pos to radius value   

        double x_start = Y[0];     

        t_fin = 2*M_PI*m_e/(e)*gamma;
        n = 360;

        dt = (t_fin - t_in)/(double) n;
        t = t_in;

        for (int i = 0; i < n; i++){
            Boris2ndOrderStep(t, dt, alpha, Y, Neq, E_Field, B_Field);
            t += dt;
        }
        cout << i << " " << sqrt(Y[0]*Y[0] + Y[1]*Y[1])/x_start << endl;

        fdata << i << " " << Y[0] - x_start << " " << Y[1] << " ";
        fdata << sqrt((Y[0] - x_start)*(Y[0] - x_start) + Y[1]*Y[1]) << " ";
        fdata << sqrt(Y[0]*Y[0] + Y[1]*Y[1])/x_start << endl;



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
    B[0] = 0. *TESLA;
    B[1] = 0. *TESLA;
    B[2] = 1. *TESLA ;
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
