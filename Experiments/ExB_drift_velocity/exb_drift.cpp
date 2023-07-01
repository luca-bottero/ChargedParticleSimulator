# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include "../ode_solvers.h"

using namespace std;

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

double E_0 = 1.;
# define B_0 -1.

void Set_Initial_Values (double *);
void RHS_Func   (double, double *, double *, double *);
void E_Field (double, double *, double *);
void B_Field (double, double *, double *);

int main(){
    cout << setiosflags(ios::scientific) << setprecision(10);

    double t;
    double t_in = 0., t_fin = 1e-10;
    int n = 1000000;

    int N_PARTICLES = 1;
    int Neq = 6 * N_PARTICLES;

    double Y[Neq];

    double alpha = -e/(m_e);

    ofstream fdata;
    fdata.precision(20);
    fdata.open("output.dat");

    // Weak field regime
    for(int i = -9; i < 0; i++){
        E_0 = 1.*pow(10,i)*C;

        double gamma = 1./sqrt(1. - E_0*E_0/(B_0*B_0*C*C));
        t_fin = 1.1 * 2. * M_PI * m_e / (abs(B_0) * e);

        double dt = (t_fin - t_in)/(double) n;    
        
        t = t_in;
        Set_Initial_Values(Y);

        bool reverse = false;

        for (int i = 0; i < n; i++){
            Boris2ndOrderStep(t, dt, alpha, Y, Neq, E_Field, B_Field);
            t += dt;
            if (Y[3] > 0.) reverse = true;
            if (Y[3] < 0. && reverse) break;
        }

        fdata << E_0 << " " << 1. - abs((Y[1]/t)/(E_0/B_0)) << " " << t/dt << endl;
    }
    // Strong field regime
    for(int i = 1; i < 6; i++){
        E_0 = C*(1. - pow(10,-i));

        double gamma = 1./sqrt(1. - E_0*E_0/(B_0*B_0*C*C));
        t_fin = 2. * M_PI *m_e / (abs(B_0) * e)*1e3;

        double dt = (t_fin - t_in)/(double) n;    
        
        t = t_in;
        Set_Initial_Values(Y);

        bool reverse = false;

        for (int i = 0; i < n; i++){
            Boris2ndOrderStep(t, dt, alpha, Y, Neq, E_Field, B_Field);
            t += dt;
            if (Y[3] > 0.) reverse = true;
            if (Y[3] < 0. && reverse) break;
        }

        fdata << E_0 << " " << 1. - abs((Y[1]/t)/(E_0/B_0)) << " " << t/dt << endl;
    }

    fdata.close();
}

void Set_Initial_Values(double * Y){
    Y[0] = 0. ;  // x
    Y[1] = 0. ;  // y
    Y[2] = 0. ;  // z
    Y[3] = 0.;  // v_x
    Y[4] = 0.;  // v_y
    Y[5] = 0.;  // v_z
}

void E_Field(double t, double * pos, double * E){
    E[0] = E_0;
    E[1] = 0.;
    E[2] = 0.;
}

void B_Field(double t, double * pos, double * B){
    B[0] = 0.;
    B[1] = 0.;
    B[2] = B_0;
}



