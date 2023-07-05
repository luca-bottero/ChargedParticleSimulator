# include "field_sources.h"
# include <math.h>

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

# define A_R1 8. * M_SQRT2 / M_PI
# define A_Z1 8. * M_SQRT2 / M_PI
# define B_Z1 2. * M_SQRT2 / M_PI

# define MAX_3N(i,j,k) i > j? (i > k? i: k): (j > k? j: k)

class Loop {
    public:
        double l_pos[3];
        double l_dir[3];    // if current in positive, + direction is given by right hand rule
        double radius;
        double B_cylindrical[3];
        double I;
        void CalculateCurrent(double t){
            this->I = 1.;
        }

        void PosLocalToAbsolute(double *pos){
            // Transform a position from the local reference frame to the absolute reference frame            
        }

        void CalculateBField(double t, double *pos){
            CalculateCurrent(t);

            double z;
            double h;
            double W;
            double alpha_z, beta_z;

            W = 2.*h*radius/(h*h + z*z + radius*radius);
            W *= W;

            // B_r

            B_cylindrical[0] = mu_0*0.25*I*radius*radius*z*h;
            B_cylindrical[0] /= pow(h*h + radius*radius + z*z, 2.5);
            B_cylindrical[0] *= 3 + A_R1*W/(1.-W);

            // B_z
            alpha_z = 3. + A_Z1*W/(1.-W);
            beta_z = 1. + B_Z1*W/(1.-W);
            B_cylindrical[1] = mu_0*I*0.25/pow(h*h + radius*radius + z*z, .5);
            B_cylindrical[1] *= 2.*radius*beta_z/(h*h + radius*radius + z*z) - 0.25*W*alpha_z;

            // trnasoform to absolute reference frame
        } 
};



