# include <iostream>
# include <math.h>
# include <iomanip>
# include <fstream>
# include "../ode_solvers.h"
# include <omp.h>

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
    # define mu_0 4. * M_PI * 1e-7
    # define e 1.602176634e-19
    # define m_e 9.1093837015e-31
    # define m_p 1.67262192369e-27
    # define TESLA 1.
    # define METER 1.
#endif

# define MAG_MOMENT 8e22    // Approximate Earth's magnetic dipole
# define RADIUS_EARTH 6371e3
# define N_ITERS 1e5
# define N_PARTICLES 10

void Set_Initial_Values (double *);
void RHS_Func   (double, double *, double *, double *);
void E_Field (double, double *, double *);
void B_Field (double, double *, double *);
void Save_B_Field_Strength(double, double, double, double);

int main(){
    cout << setiosflags(ios::scientific) << setprecision(10);

    double t;
    double t_in = 0., t_fin = 1.5;

    int Neq = 6;
    double Y[Neq];
    double alpha = e/m_p;

    double dt = (t_fin - t_in)/(double) N_ITERS;

    ofstream fdata;
    fdata.precision(20);
    fdata.open("output.dat");

    srand((unsigned) time(0));
    
    for(int j=0; j<N_PARTICLES; j++){
        t = t_in;
        Set_Initial_Values(Y);

        fdata << t << " " << Y[0]/RADIUS_EARTH << " " << Y[1]/RADIUS_EARTH << " " << Y[2]/RADIUS_EARTH <<  endl;
        for (int i = 0; i < N_ITERS; i++){
            Boris2ndOrderStep(t, dt, alpha, Y, Neq, E_Field, B_Field);
            t += dt;
            fdata << t << " " << Y[0]/RADIUS_EARTH << " " << Y[1]/RADIUS_EARTH << " " << Y[2]/RADIUS_EARTH <<  endl;
        }
        fdata << endl << endl;
    }

    fdata.close();
}

void Set_Initial_Values(double * Y){
    double K = 100e6;           // 100 MeV
    double M_eV = 938.272e6;    // Proton mass in eV

    double theta_pos = (double)rand()/(double) RAND_MAX*2*M_PI;
    double r = RADIUS_EARTH*(1.2 + 1.8*(double)rand()/(double)RAND_MAX); // Inner van Allen Belt

    double theta_vel = (2. - 2.*(double)rand()/(double) RAND_MAX)*M_PI;
    double phi_vel = (double)rand()/(double) RAND_MAX*2*M_PI;
    K *= (0.8 + 0.7*(double)rand()/(double) RAND_MAX);      // Randomize kinetic energy of particle
    double speed = C*sqrt(1. - 1./((1 + K/(M_eV))*(1 + K/(M_eV))));

    Y[0] = r * cos(theta_pos);  // x
    Y[1] = r * sin(theta_pos);  // y
    Y[2] = 0. ;  // z
    Y[3] = speed * sin(theta_vel)*cos(phi_vel);  // v_x
    Y[4] = speed * sin(theta_vel)*cos(theta_vel);  // v_y
    Y[5] = speed * cos(theta_vel);  // v_z

    double B[3];
    double pos[3] = {Y[0], Y[1], Y[2]};
    B_Field(0., pos, B);
    r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
}

void E_Field(double t, double * pos, double * E){
    E[0] = 0.;
    E[1] = 0.;
    E[2] = 0.;
}

void B_Field(double t, double * pos, double * B){
    double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    double coupling_const = MAG_MOMENT/(r*r*r*r*r)*mu_0*0.25/M_PI;

    B[0] = 3. * pos[0] * pos[2] * coupling_const;
    B[1] = 3. * pos[1] * pos[2] * coupling_const;
    B[2] = 3. * (3.*pos[2]*pos[2] - r*r) * coupling_const;
}



