# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include "../ode_solvers.h"
# include <omp.h>

using namespace std;

# define N_ITERS 1000
# define N_PARTICLES 1000
# define L 10.

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

void Set_Initial_Values (double *);
void RHS_Func   (double, double *, double *, double *);
void E_Field (double, double *, double *);
void B_Field (double, double *, double *);
void Save_B_Field_Strength(double, double, double, double);

int main(){
    cout << setiosflags(ios::scientific) << setprecision(10);

    double t;
    double t_in = 0., t_fin = 2e-4;

    int Neq = 6;
    double Y[Neq];
    double alpha = e/m_p;

    double dt = (t_fin - t_in)/(double) N_ITERS;

    ofstream fdata;
    fdata.precision(20);
    fdata.open("output.dat");

    srand((unsigned) time(0));

    Save_B_Field_Strength(-10.*L, 10.*L, -10.*L, 10.*L);
    
    for(int j=0; j<N_PARTICLES; j++){
        t = t_in;
        Set_Initial_Values(Y);

        fdata << t << " " << Y[0] << " " << Y[1] << " " << Y[2] <<  endl;
        for (int i = 0; i < N_ITERS; i++){
            Boris2ndOrderStep(t, dt, alpha, Y, Neq, E_Field, B_Field);
            t += dt;
            fdata << t << " " << Y[0] << " " << Y[1] << " " << Y[2] <<  endl;
        }
        fdata << endl << endl;
    }

    fdata.close();
}

void Save_B_Field_Strength(double x_min, double x_max, double y_min, double y_max){
    ofstream fdata;
    fdata.precision(20);
    fdata.open("b_field.dat");

    double B[3];
    double pos[3];
    double dx = (x_max - x_min)/20.;
    double dy = (y_max - y_min)/20.;

    for(double x=x_min; x < x_max; x += dx){
        for(double y=y_min; y < y_max; y += dy){
            pos[0] = x;
            pos[1] = y;
            pos[2] = 0.;

            B_Field(0., pos, B);
            fdata << pos[0] << " " << pos[1] << " " << B[0] << " " << B[1] << " " << B[2] << " ";
            fdata << sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) << endl;
        }
    }
}

void Set_Initial_Values(double * Y){
    double theta = (double)rand()/(double) RAND_MAX*2*M_PI;

    Y[0] = ((double) rand()/(double) RAND_MAX * 2. - 1.)*L ;  // x
    Y[1] = ((double) rand()/(double) RAND_MAX * 2. - 1.)*L ;  // y
    Y[2] = 0. ;  // z
    Y[3] = 0.1 * cos(theta)*C;  // v_x
    Y[4] = 0.1 * sin(theta)*C;  // v_y
    Y[5] = 0.;  // v_z
}

void E_Field(double t, double * pos, double * E){
    E[0] = 0.;
    E[1] = 0.;
    E[2] = 0.5*C;
}

void B_Field(double t, double * pos, double * B){
    B[0] = pos[1]/L;
    B[1] = pos[0]/L;
    B[2] = 0.;
}



