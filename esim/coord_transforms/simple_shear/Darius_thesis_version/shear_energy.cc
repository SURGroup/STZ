#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include "common.hh"
#include "shear_sim.hh"
#include <iostream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000; // STZ formation energy, eZ/kB [K]
const double kB=1.3806503e-23; // Boltzman constant, [J/K]
const double kB_metal = 8.617330350e-5; // Boltzmann constant, [eV/K]
const double eV_per_J = 6.241509e18; // Converts eV to Joules
const double Ez = TZ * kB * eV_per_J; // STZ formation Energy, [eV]

void syntax_message() {
    fputs("Syntax: ./shear_energy <run_type> <input_file> <beta> <u0>\n\n"
          "<run_type> is either:\n"
          "     \"direct\" for direct simulation, or\n"
          "     \"qs\" for quasi-static simulation. \n\n"
          "<input_file> expects the relative path for a text file containing potential energy data. \n\n"
          "<beta> is the initial energy scaling parameter, [1/eV]. \n\n"
          "<u0> is an energy offset value, [eV].\n",stderr);
    exit(1);
}

int main(int argc,char **argv) {

    // Check command-line arguments
    if (argc < 5) syntax_message();
    bool qs = true;
    if (strcmp(argv[1], "direct") == 0) qs = false;
    else if (strcmp(argv[1], "qs") != 0) syntax_message();
    double beta, u0;
    beta = strtod(argv[3],NULL);
    u0 = strtod(argv[4],NULL);

    // Set default value for yield stress
    double rho_ = 6125; // density (kg/m^3)
    double sy_ = 0.85; // yield stress (GPa)
    double Delta_ = 8000; // activation temperature (K)
    double Omega_ = 300; // activation volume (A**3)
    double l_chi = 4.01;
    double c0_ = 0.3; // plastic work fraction
    double eps0_ = 0.03; // Local strain
    double mu_ = 20; // shear modulus (GPa)
    double K_ = 143.791; // bulk modulus (GPa)
    double theta_ = 300; // bath temperature (K)
    double chiinf_ = 2730; // upper limiting effective temperature (K)
    bool chiinf_bool = true; //
    double kappa_ = 0; // Thermal diffusion length scale
    double Omegaeps0_ = 11.647; // Activation volume and local strain product (A^3)

    string endref = "pe.Cu.MD.100.txt";
    ifstream endreffile; // Final reference MD PE file

    if (argc > 5)
    {
        printf("\nThe following values were set at runtime:\n\n");
        int i_arg = 5;
        while (i_arg < argc) 
        {
            if (strcmp(argv[i_arg], "chi_inf") == 0)
            {
                i_arg++;
                chiinf_ = strtod(argv[i_arg], NULL);
                chiinf_bool = true;
                printf("Steady-state effective temperature chi_inf = %g K.\n", chiinf_);
            }  
            if (strcmp(argv[i_arg], "theta") == 0)
            {
                i_arg++;
                theta_ =strtod(argv[i_arg], NULL);
                printf("Bath temperature theta = %g K.\n", theta_);
            }
            if (strcmp(argv[i_arg], "sy") == 0)
            {
                i_arg++;
                sy_ = strtod(argv[i_arg],NULL);
                printf("Yield stress sy = %g GPa.\n", sy_);
            }
            if (strcmp(argv[i_arg], "endref") == 0)
            {
                i_arg++;
                string endref = argv[i_arg];
                cout << "Reading final state from file " << endref << "." << endl;
            } 
            if (strcmp(argv[i_arg], "Delta") == 0)
            {
                i_arg++;
                Delta_ = strtod(argv[i_arg], NULL);
                printf("STZ activation barrier Delta = %g K.\n", Delta_);
            }
            if (strcmp(argv[i_arg], "Omega") == 0)
            {
                i_arg++;
                Omega_ = strtod(argv[i_arg], NULL);
                printf("STZ activation volume Omega = %g A^3.\n", Omega_);
            }
            if (strcmp(argv[i_arg], "chi_len") == 0)
            {
                i_arg++;
                l_chi = strtod(argv[i_arg], NULL);
                printf("Chi diffusion length scale l_chi =  %g Anstroms.\n", l_chi);
            }
            if (strcmp(argv[i_arg], "c0") == 0)
            {
                i_arg++;
                c0_ = strtod(argv[i_arg], NULL);
                printf("Plastic work fraction c0 = %g.\n", c0_);
            }
            if (strcmp(argv[i_arg], "mu") == 0)
            {
                i_arg++;
                mu_ = strtod(argv[i_arg], NULL);
                printf("Shear modulus mu = %g GPa.\n", mu_);
            }
            if (strcmp(argv[i_arg], "K") == 0)
            {
                i_arg++;
                K_ = strtod(argv[i_arg], NULL);
                printf("Bulk modulus K = %g GPa.\n", K_);
            }
            if (strcmp(argv[i_arg], "rho") == 0)
            {
                i_arg++;
                rho_ = strtod(argv[i_arg], NULL);
                printf("Density rho = %g kg/m^3.\n", rho_);
            }
            if (strcmp(argv[i_arg], "eps0") == 0)
            {
                i_arg++;
                eps0_ = strtod(argv[i_arg], NULL);
                printf("Local strain eps0 = %g.\n", eps0_);
            }
            if (strcmp(argv[i_arg], "Omegaeps0") == 0)
            {
                i_arg++;
                Omegaeps0_ = strtod(argv[i_arg], NULL);
                printf("STZ activation volume Omega times epsilon0 = %g A^3.\n", Omegaeps0_);
            } 
            if (strcmp(argv[i_arg], "kappa") == 0)
            {
                i_arg++;
                kappa_ = strtod(argv[i_arg], NULL);
                printf("kappa = %g.\n", kappa_);
            }
            i_arg++;
        }
    }

    if (chiinf_bool)
    {   
        chiinf_ /= TZ;
        printf("\nThe upper-limiting effective temperature was set using a user-defined value.\n");
        printf("The effective temperature is %g K.\n", chiinf_*TZ);
    }
    else 
    {
        printf("\nSetting the upper-limiting effective temperature T_inf using");
        printf("\nthe maximum pe value obtained in the final MD snapshot.\n\n");
        // Open final MD PE reference file and determine maximum PE
        char cstr[endref.size() + 1];
        vector<double> refPEfin;
        strcpy(cstr, &endref[0]);

        endreffile.open(cstr);

        if (!endreffile.is_open()) {
            std::cerr << "There was an issue opening the file." << endl;
            exit(1);
        }
    
        double num = 0.0;
        endreffile.ignore(80, '\n'); // skip the single line header
    
        while (endreffile >> num) {
            refPEfin.push_back(num);
        }
    
        endreffile.close();
        double maxPEf = *max_element(refPEfin.begin(), refPEfin.end());
        chiinf_ = 0.95*beta*(maxPEf - u0)*Ez/kB_metal/TZ;
        printf("The maximum PE from the final MD snapshot is %g eV.\n", maxPEf);
        printf("For beta = %g and E0 = %g eV, T_inf = %g K.\n\n", beta, u0, chiinf_*TZ);
    }     
    
    // Elasticity Parameters
    const double s_y = sy_*1e9;                 // Yield stress (Pa)
    
    // STZ model parameters (based on Vitreloy 1 BMG)
    const double c0 = c0_;
    const double tau0 = 1e-13;                  // Vibration timescale, (s)
    const double kappa = kappa_;
    const double Delta = Delta_/TZ;
    const double Omega = Omega_*22/17;        // Fix typo in 2015 JCP
    const double eps0 = eps0_;
    const double theta = theta_/TZ; // Bath Temperature, [--]
    const double Omegaeps0 = Omegaeps0_*1e-30*s_y/(TZ*kB);
    const double chi_inf = chiinf_;

    // Output filenames
    const char dfn[] = "sct_d.out", qfn[] = "sct_q.out";

    // Output fields and miscellaneous flags. 1-u, 2-v, 4-p, 8-q, 16-s, 32-tau,
    // 64-chi, 128-tem, 256-dev, 512-X, 1024-Y, 2048-(total strain components),
    // 4096-(total strain invariants), 8192-(Lagrangian tracer output).
    // const unsigned int fflags=32|128|512|1024|2048|4096|8192;
    const unsigned int fflags=32|128|512|1024|2048|4096;

    // Elasticity related parameters (based on Vitreloy 1 BMG)
    //const double s_y = sy_*1e9;                 // Yield stress (GPa)
    const double K = K_*1e9/s_y;            // Rescaled bulk modulus (--)
    //const double mu = mu_*1e9/s_y;                 // Rescaled Shear modulus (--)    
    const double mu = mu_;      // Rescaled Shear Modulus (GPa)

    // Other parameters. Note that parameters labeled (*) are far larger
    // than realistic values, but allow for a direct--quasistatic
    // comparison.
    //const double le = 0.08;                 // Length scale (m)
    const double rho0 = rho_ ;               // Density (kg/m^3)
    const double sca = 2e4;
    //const double t_scale = le/sqrt(mu*s_y/rho0)*sca; // Plasticity timescale (*)
    const double visc = 0.02;               // Viscous damping
    //const double chi_len=l_chi*1e-10;              // Dpl-mediated diffusion (m)
    const double chi_len=l_chi;              // Dpl-mediated diffusion (m)
    //const double adapt_fac = 2e-3; 
    const double adapt_fac = 1e-4; 

    // MD Simulation-Based Continuum Parameters
    const int x_grid = 32;              // grid-points in x-direction
    const int y_grid = 32;              // grid-points in y-direction
    const double x_beg = 0.0;           // x-origin (m)
    //const double x_end = 4e-8;         // x-terminus (m)
    //const double x_end = 400e-10;         // x-terminus (m)
    const double x_end = 400;         // Thesis x-terminus (m)
    const double y_beg = 0.0;           // y-origin (m)
    //const double y_end = 4e-8;         // y-terminus (m)
    //const double y_end = 400e-10;         // y-terminus (m)
    const double y_end = 400;         // Thesis  y-terminus (m)
    const double gam_max = 0.5;         // MD max strain
    
    //double dx = (x_end-x_beg)/x_grid;
    //double dx = 1.1e-2;
    double dx = 4e-8; // MD length (m)
    //const double t_scale = dx/sqrt(mu*s_y/rho0)*sca;
    //const double t_scale = dx/sqrt(mu*s_y/rho0);

    //const double t_scale = 0.01/sqrt(mu*s_y/rho0); // Thesis
    const double t_scale = 0.011/sqrt(mu*s_y/rho0)*sca; // Thesis

    // Set parameters for either a periodic or non-periodic test
    const bool y_prd=true;
    double u_bdry, lamb, tf;
    if(y_prd) {

        // Periodic test
        u_bdry=0.;
        //lamb=1e8; // strain rate (--/1)
        lamb=1e-4; // Thesis
        tf=gam_max/lamb; // final time (s)
    } else {

        // Non-periodic test
        u_bdry=1e-7*sca;
        lamb=0;tf=3e6/sca;
    }

    // Make the output directory if it doesn't already exist
    mkdir(qs? qfn : dfn, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Initialize an STZ plasticity model
    //stz_dynamics_linear_athermal stz(TZ, c0, tau0, kappa, Delta, Omegaeps0, chi_inf, theta, 0);
    stz_dynamics_linear stz(TZ, c0, tau0, kappa, Delta, Omegaeps0, chi_inf, theta, 0);
    // Initialize the simulation
    shear_sim sim(x_grid, y_grid, x_beg, x_end, y_beg, y_end, mu, K, visc, chi_len, t_scale, adapt_fac, u_bdry, lamb, &stz, y_prd, fflags, qs? qfn : dfn);
    //shear_sim sim(x_grid, y_grid, 0, 1, 0, 1, mu, K, visc, chi_len, t_scale, adapt_fac, u_bdry, lamb, &stz, y_prd, fflags, qs? qfn : dfn);
    sim.init_fields(0, 580, 220);

    // Open the input file and read in the grid dimensions and scaling
    FILE *fp=safe_fopen(argv[2],"r");
    int mm,nn;
    double chi_base,chi_scale;
    if(fscanf(fp,"%d %d %lg %lg",&mm,&nn,&chi_base,&chi_scale)!=4) {
        fputs("Error reading header information\n",stderr);
        return 1;
    }

    // Check that the dimensions make sense, and allocate memory
    if(mm<=0 || nn<=0) {
        fputs("Grid dimensions are invalid\n",stderr);
        return 1;
    }
    double *f=new double[mm*nn];

    // Read in the energy values from the file
    for(int j=0;j<nn;j++) {
        for(int i=0;i<mm;i++) {
            if(fscanf(fp,"%lg",f+i+mm*j)!=1) {
                fputs("Error reading energy information\n",stderr);
                return 1;
            }
        }
    }

    // Map energy values onto effective temperature, chi
    for(int j=0; j<nn; j++) {
        for(int i=0; i<mm; i++) {
            f[i+mm*j] = beta*(f[i+mm*j]-u0)*Ez/kB_metal;
        }
    }

    printf("# chi grid dimensions : %d by %d\n"
           "# chi base scale      : %g K\n"
           "# chi range scale     : %g K\n",mm,nn,chi_base,chi_scale);
    sim.initialize_chi_bicubic(mm,nn,f,chi_base/TZ,chi_scale/TZ);

    // Free the dynamically allocated memory
    delete [] f;

    // Carry out the simulation using the selected simulation method
    // qs?sim.solve_quasistatic(tf,n_frames,5):sim.solve(tf,160);
    int n_frames = 100;
    int steps = 100;
    double app_t = tf/n_frames/steps;
    // Output to console
    printf("\nThe plasticity timestep is : %g s.\n", t_scale);
    printf("The continuum timestep is : %g s.\n", app_t);
    printf("The final time is : %g s.\n\n", tf);
    qs?sim.solve_quasistatic(tf,n_frames,steps):sim.solve(tf,160);
    
    exit(0);
}
