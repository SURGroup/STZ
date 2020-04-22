#include <cstring>
#include "common.hh"
#include "shear_sim.hh"
#include "mat.hh"

using std::isnan;
using std::abs;

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for the fields, and calls the
 * routine to initialize the fields.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *              vertical directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] mu_ the elastic shear modulus.
 * \param[in] K_ the elastic bulk modulus.
 * \param[in] viscosity_ the viscosity.
 * \param[in] t_scale_ a time scale associated with the plastic deformation.
 * \param[in] u_bdry_ the boundary velocity.
 * \param[in] adapt_fac_ a factor used in the adaptive timestepping of
 *              plasticity.
 * \param[in] filename_ the filename of the output directory. */
shear_sim::shear_sim(const int m_, const int n_, const double ax_, const double bx_,
        const double ay_, const double by_, const double mu_,
        const double K_, const double viscosity_, const double chi_len_,
        const double t_scale_, const double adapt_fac_, const double u_bdry_,
        const double lamb_, stz_dynamics *stz_, bool y_prd_,
        const unsigned int fflags_, const char *filename_)
    : m(m_), gm(m+4), n(n_), lsn(y_prd_?n:n+1), gn(n+4), mn(m*n), gmn(gm*gn),
    fflags(fflags_), ax(ax_), ay(ay_), bx(bx_), by(by_), lx(bx-ax), ly(by-ay),
    dx(lx/m), dy(ly/n), xsp(1./dx), ysp(1./dy), mu(mu_), mu_inv(1./mu_), K(K_),
    viscosity(viscosity_), chi_len(chi_len_), t_scale(t_scale_), adapt_fac(adapt_fac_),
    u_bdry(u_bdry_), lamb(lamb_), stz(stz_), y_prd(y_prd_),
    filename(filename_), fbase(new c_field[gmn]), fm(fbase+2*gm+2), tr(NULL),
    time(0.), f_num(0), qsm(*this), buf(new float[m>=63?m+6:64]) {

    // Set STZ related parameters
    TZ=stz->TZ;
    chi_inf=stz->chi_inf;
}

/** The class destructor frees the dynamically allocated memory. */
shear_sim::~shear_sim() {
    delete [] buf;
    if(tr!=NULL) {
        delete [] tr;
        for(int j=0;j<4;j++) delete bi[j];
    }
    delete [] fbase;
}

/** Initializes the simulation fields.
 * \param[in] chi_case the type of setup to use for chi.
 * \param[in] tem_base the base value of effective temperature (in Kelvin) to
 *              use.
 * \param[in] tem_delta the delta value of effective temperature (in Kelvin) to
 *              use. */
void shear_sim::init_fields(int chi_case, double tem_base, double tem_delta) {
    double chi0 = tem_base/TZ, chi1 = tem_delta/TZ;

    // Initialize table for sine wave case
    double sx[16], sy[16];
    if(chi_case == 2) {
        for(int i = 0; i < 16; i++) {
            sx[i] = ax + lx*(i + 0.5)*0.0625;
            sy[i] = 0.2*sin(M_PI*0.0625*(2*i+1));
        }
    }

    // Loop over the rows of the grid
#pragma omp parallel for
    for(int j = 0; j < ( y_prd? n : n+1 ); j++) {

        // Calculate the cell-cornered and cell-centered y-reference maps
        double Y = ay + j*dy, Ys = Y + 0.5*dy;

        // Pointer to the first element in row j
        c_field *fp = fm + gm*j;

        // Loop over row j
        for(int i = 0; i < m; i++, fp++) {

            // Reference to the current grid location
            c_field &f = *fp;

            // Calculate the cell-cornered and cell-centered x-reference maps
            double X = ax + i*dx, Xs = X + 0.5*dx;

            // Set stresses and velocities so that, if a velocity is input,
            // it matches up with the corresponding lambda transformation.
            f.u = Y*u_bdry; f.v = 0;
            f.p = f.q = f.s = f.tau = 0;

            // Zeroth case: a small Gaussian blip in chi at the
            // origin
            if (chi_case == 0) f.chi = chi0 + chi1*exp(-200*(Xs*Xs+Ys*Ys));
            else if (chi_case == 1) {

                // First case: a rotated line of higher chi
                const double csa = cos(30*M_PI/180.), sia = sin(30*M_PI/180.);
                double Xr = csa*Xs + sia*Ys, Yr = -sia*Xs + csa*Ys;
                Xr = (Xr < -1)? Xr+1 : ((Xr > 1)? Xr-1 : 0);
                f.chi = chi0 + chi1*exp(-200*(Xr*Xr + Yr*Yr));
            } else {

                // Second case: a sequence of Gaussian blips in
                // a sine wave pattern
                double R=1e100, Rt, Xr, Yr;
                for(int k = 0; k < 16; k++) {
                    Xr = sx[k] - Xs;
                    if (Xr<-lx*0.5) Xr+=lx;
                    if (Xr>lx*0.5) Xr-=lx;
                    Yr = sy[k] - Ys;
                    Rt = Xr*Xr + Yr*Yr;
                    if (Rt<R) R = Rt;
                }
                f.chi = chi0 + chi1*exp(-200*R);
            }

            // Set reference map field to match initial coordinate
            f.X = X; f.Y = Y;
        }
    }
    set_boundaries();
}

/** Carries out the simulation for a specified time interval using the direct
 * simulation method, periodically saving the output.
 * \param[in] t_start the simulation time to start from.
 * \param[in] t_end the simulation time to end at.
 * \param[in] frames the number of frames to save. */
void shear_sim::solve(double duration, int frames) {
    double time_interval = duration/frames, t0, t1, t2;

    // Choose a timestep to be lower than the restrictions due to elastic waves
    // and viscosity
    const double cfl       = 0.25*dx,
                 visc_cond = 0.25*dx*dx / viscosity,
                 dt = (cfl < visc_cond)? cfl : visc_cond;
    printf("CFL: %g\n", dt/dx);
    printf("Viscosity Condition: %g\n", viscosity*dt*(1./(xsp*xsp + ysp*ysp)));
    const int l=static_cast<int>(time_interval/dt)+1;
    const double adt=time_interval/l;

    // Output the initial fields and record initial time
    write_files(0);
    puts("# Output frame 0");
    t0 = wtime();

    for(int k = 1; k <= frames; k++) {

        // Perform the explicit timestep updates, and then check the grid
        // for any strange behavior
        for(int j=0;j<l;j++) step_forward(adt);
        if(fflags&4096) check_grid_wtf();

        // Output the fields
        t1 = wtime();
        write_files(k);

        // Print diagnostic information
        t2 = wtime();
        printf("# Output frame %d [%d, %.8g s, %.8g s]\n", k, l, t1 - t0, t2 - t1);
        t0 = t2;
    }
}

/** Carries out the simulation for a specified time interval using the
 * quasi-static simulation method, periodically saving the output.
 * \param[in] t_start the simulation time to start from.
 * \param[in] t_end the simulation time to end at.
 * \param[in] frames the number of frames to save.
 * \param[in] steps the number of quasi-static steps to take per frame. */
void shear_sim::solve_quasistatic(double duration, int frames, int steps) {
    double time_interval = duration/frames, t0, t1, t2, isteps = 1.0/steps;

    // Output the initial fields
    if(f_num==0) {write_files(0);puts("# Output frame 0");}
    t0 = wtime();

    for(int k = 1; k <= frames; k++) {

        // Perform the quasistatic timestep updates, and then check the grid
        // for any strange behavior
        for(int j = 0; j < steps; j++) step_forward_quasistatic(time_interval*isteps);
        if(fflags&4096) check_grid_wtf();
        t1 = wtime();

        // Output the fields
        write_files(k);

        // Print diagnostic information
        t2 = wtime();
        printf("# Output frame %d [%d, %.8g s, %.8g s] {%.2f}\n", k, steps, t1-t0, t2-t1, qsm.tp.avg_iters());
        t0 = t2;
    }
}

/** Checks the grid for any strange behavior, looking for large field values or
 * NaNs. Note that this currently only scans one layer of ghost points, since
 * for the non-periodic case some points in the second ghost layer aren't used.
 */
void shear_sim::check_grid_wtf() {
    c_field curr_field;
    for (int jj = -1; jj < n+1; jj++)
        for (int ii = -1; ii < m+1; ii++){
            curr_field = *(fm + jj*gm + ii);
            if (curr_field.weird()){
                printf("wtf on grid point: (%d, %d)\n", ii, jj);
                exit(1);
            }
        }
}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void shear_sim::step_forward(double dt) {

    // Declare some useful constants. L is the first Lamé parameter.
    const double third = 1./3, L = K - 2*third*mu, lt = lamb*time, ltsq = lt*lt;
    double hx = 0.5*xsp*dt, hy = 0.5*ysp*dt,
           hxxv = 2*xsp*hx*viscosity, hyyv = 2*ysp*hy*viscosity,
           hcv  = 2*(hxxv + hyyv);

#pragma omp parallel for
    for(int j = 0; j < n; j++) for(int i = 0; i < m; i++) {
        st_field f1, f2;

        // Create references to the fields in the neighboring gridpoints,
        c_field *fp = fm + (gm*j + i), *flp = fp - 1, *frp = fp + 1;
        c_field &f = *fp, &fl = *flp, &fr = *frp, &fd = fp[-gm], &fu = fp[gm];

        // uc, vc hold the current (central?) velocities
        // ux, vx will hold the upwinded x-derivatives of the velocities
        // uy, vy will hold the upwinded y-derivatives of the velocities
        // uvisc, vvisc will hold the viscous smoothing terms
        // and dtmp will hold the result of the adaptive_plastic_term
        // algorithm, which is (2*mu*Dpl*dt/sbar)
        double uc, vc, ux, vx, uy, vy, uvisc, vvisc, dtmp;

        // Unless this is the zeroth row, or it's a y-periodic domain,
        // calculate the update to the reference map and velocity
        if(y_prd || j > 0) {
            double Xx, Xy, Yx, Yy, fx, fy;
            uc = f.u; vc = f.v;

            // Calculate ENO derivatives of regular fields
            uc>0?rmv_eno2(ux, vx, Xx, Yx, hx, fp, -1):rmv_eno2(ux, vx, Xx, Yx, -hx, fp, 1);
            vc>0?rmv_eno2(uy, vy, Xy, Yy, hy, fp, -gm):rmv_eno2(uy, vy, Xy, Yy, -hy, fp, gm);

            // Calculate net force due to stress imbalance
            // Fills fx = -dp'_dx - dq'_dx + ds'_dx + dtau_dy
            //       fy = -dp'_dy - dq'_dy - ds'_dy + dtau_dx
            // Here we use hxa and hxb
            net_force(fx, fy, -hx, -hy, fp);

            // Calculate viscosity. The second and third line in both
            // calculations is the extra term for computing the Laplacian in
            // the untransformed coordinates.
            uvisc = hxxv*(fr.u + fl.u) + hyyv*(fu.u + fd.u) - hcv*uc
                + dt*ltsq*xsp*xsp*(fr.u - 2*f.u + fl.u)
                - .5*lt*dt*xsp*ysp*(frp[gm].u - frp[-gm].u - flp[gm].u + flp[-gm].u);

            vvisc = hxxv*(fr.v + fl.v) + hyyv*(fu.v + fd.v) - hcv*vc
                + dt*ltsq*xsp*xsp*(fr.v - 2*f.v + fl.v)
                - .5*lt*dt*xsp*ysp*(frp[gm].v - frp[-gm].v - flp[gm].v + flp[-gm].v);

            // Apply the update to velocity. Recall that in dimensionless units
            // mu = rho.
            f.cu = -uc*ux - vc*uy + mu_inv*(fx + uvisc) - 2*dt*lamb*vc;
            f.cv = -uc*vx - vc*vy + mu_inv*(fy + vvisc);

            // Calculate update to reference map fields. Simple advective
            // equation dX/dt = dY/dt = 0.
            f.cX = -uc*Xx - vc*Xy;
            f.cY = -uc*Yx - vc*Yy;
        }

        // Calculate deformation rate, staggered velocity, and spin
        uc = 0.25*(f.u + fr.u + fu.u + frp[gm].u);
        vc = 0.25*(f.v + fr.v + fu.v + frp[gm].v);
        ux =  hx*(-f.u + fr.u - fu.u + frp[gm].u);
        uy =  hy*(-f.u - fr.u + fu.u + frp[gm].u);
        vx =  hx*(-f.v + fr.v - fu.v + frp[gm].v);
        vy =  hy*(-f.v - fr.v + fu.v + frp[gm].v);
        double old_p = calc_p(f);

        // Calculate ENO derivatives of staggered fields
        f1 = uc>0?st_eno2(hx, fp, -1):st_eno2(-hx, fp, 1);
        f2 = vc>0?st_eno2(hy, fp, -gm):st_eno2(-hy, fp, gm);

        // Calculate the adaptive effective temperature update
        f.cchi  = f.chi;
        dtmp    = adaptive_plastic_term(f.dev(lt), f.cchi, f.ddev, dt);
        f.cchi -= uc*f1.chi + vc*f2.chi;

        // Calculate updates to stress
/*
 *        f.cp   = -uc*f1.p - vc*f2.p + 2*lamb*lt*third*mu*dt
 *            - third*( dtmp*(2*lt*f.tau - ltsq*(f.p + f.q + f.s + old_p))
 *            - 2*(f.tau - lt*mu)*(uy + vx) - (K + f.p - 2*f.q)*(ux + vy)
 *            + 2*f.s*two_omega - ltsq*( (K - 2*third*mu)*vy + (K + 4*third*mu)*ux ) );
 *
 *        f.cq   = -uc*f1.q - vc*f2.q - dtmp*(f.q - .5*third*ltsq*old_p)
 *            + third*(dt*ltsq*third*mu - (f.tau - lt*mu)*(uy + vx) + (f.p - 2*f.q - mu)*(ux + vy)
 *            - .5*ltsq*((K - 2*third*mu)*vy + (K + 16*third*mu)*ux)
 *            + f.s*(ux + vy) );
 *
 *        f.cs   = -uc*f1.s - vc*f2.s - lamb*lt*mu*dt + .5*ltsq*((K - 2*third*mu)*vy + (K + 4*third)*ux)
 *            - two_omega*(f.tau - lt*mu) + (mu - f.p - f.q)*(ux - vy) - dtmp*(f.s + .5*ltsq*old_p);
 *
 *        f.ctau = -uc*f1.tau - vc*f2.tau + dt*lamb*mu + ltsq*mu*vx - lt*(K + third*mu)*(ux + vy)
 *            + two_omega*f.s + (mu - f.p - f.q)*(vx + uy) - dtmp*(f.tau - lt*old_p);
 */

        f.cp   = -uc*f1.p - vc*f2.p + 2*lamb*lt*third*mu*dt
            + third*( -dtmp*(2*lt*f.tau - ltsq*(f.p + f.q + f.s + old_p))
            - 2*(f.tau - lt*mu)*(uy + vx) - (3*L + 2*mu + f.p - 2*f.q)*(ux + vy)
            + 2*f.s*(vy - ux) - ltsq*(L*(vy + ux) + 2*mu*ux));

        f.cq   = -uc*f1.q - vc*f2.q - dtmp*(f.q - .5*third*ltsq*old_p)
            + third*(dt*lt*lamb*mu - (f.tau - lt*mu)*(uy + vx) + (f.p - 2*f.q - mu - .5*ltsq*L)*(ux + vy)
            - ltsq*mu*ux + f.s*(vy - ux) );

        f.cs   = -uc*f1.s - vc*f2.s - lamb*lt*mu*dt - dtmp*(f.s + .5*ltsq*old_p)
            - (f.tau - lt*mu)*(vx - uy) + (f.p + f.q - mu)*(vy - ux)
            + .5*ltsq*L*(ux + vy) + ltsq*mu*ux;

        f.ctau = -uc*f1.tau - vc*f2.tau + dt*lamb*mu + ltsq*mu*vx - lt*(mu + L)*(ux + vy)
            + (vx - uy)*f.s + (mu - f.p - f.q)*(vx + uy) - dtmp*(f.tau - lt*old_p);

    }

    // Add in the chi-diffusion
    chi_diffusion(dt);

    // Apply the updates to the fields
#pragma omp parallel for
    for(int j = 0; j < n; j++) {
        if (j>0) {for(int i = 0; i < m; i++) fm[i + gm*j].update();}
        else for(int i = 0; i < m; i++) fm[i + gm*j].update_staggered();
    }

    // Update the simulation time, and set the boundaries
    time += dt;
    set_boundaries();
}

/** Steps the simulation fields forward using the quasistatic update procedure.
 * \param[in] dt the time step to use. */
void shear_sim::step_forward_quasistatic(double dt) {

    // Carry out the intermediate step that considers the advective and
    // plasticity terms
    advection_step(dt);

    // Apply the projection step to ensure quasistaticity
    projection_step(dt);
}

/** Carries out an intermediate step forward, by evaluating the advective and
 * plasticity terms in the equations only.
 * \param[in] dt the time step to use. */
void shear_sim::advection_step(double dt) {
    double hx = 0.5*xsp*dt, hy = 0.5*ysp*dt;
    double third = 1./3., lt = lamb*time, ltsq = lt*lt;

#pragma omp parallel for
    for(int j = 0; j < n; j++) for(int i = 0; i < m; i++) {
        double uc, vc, ux, uy, vx, vy, old_p, dtmp;
        st_field f1, f2;

        // Create references to the fields in the neighboring gridpoints,
        c_field *fp = fm + (i + gm*j), &f=*fp, &fr=fp[1], &fu=fp[gm], &fru=fp[gm+1];

        // Calculate update to reference map fields
        if (j>0||y_prd) {
            double Xx, Xy, Yx, Yy;
            uc = f.u; vc = f.v;

            // Calculate ENO derivatives of the reference map
            uc>0?rm_eno2(Xx, Yx, hx, fp, -1):rm_eno2(Xx, Yx, -hx, fp, 1);
            vc>0?rm_eno2(Xy, Yy, hy, fp, -gm):rm_eno2(Xy, Yy, -hy, fp, gm);

            // Calculate update to reference map fields
            f.cX = -uc*Xx - vc*Xy;
            f.cY = -uc*Yx - vc*Yy;
        }

        // Calculate deformation rate, staggered velocity, and spin
        uc = 0.25*(f.u + fr.u + fu.u + fru.u);
        vc = 0.25*(f.v + fr.v + fu.v + fru.v);
        ux =  hx*(-f.u + fr.u - fu.u + fru.u);
        uy =  hy*(-f.u - fr.u + fu.u + fru.u);
        vx =  hx*(-f.v + fr.v - fu.v + fru.v);
        vy =  hy*(-f.v - fr.v + fu.v + fru.v);
        old_p = calc_p(f);

        // Calculate ENO derivatives of staggered fields
        // f1 will hold all staggered ENO x-derivatives
        // f2 will hold all staggered ENO y-derivatives
        f1=uc>0?st_eno2(hx, fp, -1):st_eno2(-hx, fp, 1);
        f2=vc>0?st_eno2(hy, fp, -gm):st_eno2(-hy, fp, gm);

        // Calculate the adaptive effective temperature update
        f.cchi  = f.chi;
        dtmp    = adaptive_plastic_term(f.dev(lt), f.cchi, f.ddev, dt);
        f.cchi -= uc*f1.chi + vc*f2.chi;

        // Calculate updates to stresses
        f.cp   = -uc*f1.p - vc*f2.p - third*dtmp*(2*lt*f.tau - ltsq*(f.p + f.q + f.s) - ltsq*old_p)
               - 2*third*f.tau*(uy + vx) - third*(f.p - 2*f.q)*(ux + vy) + 2*third*f.s*(vy - ux);
        f.cq   = -uc*f1.q - vc*f2.q - dtmp*(f.q - ltsq/6*old_p) - third*f.tau*(uy + vx)
               + third*(f.p - 2*f.q)*(ux + vy) + third*f.s*(vy - ux);
        f.cs   = -uc*f1.s - vc*f2.s - dtmp*(f.s + ltsq/2*old_p) - f.tau*(vx - uy) + (f.p + f.q)*(vy - ux);
        f.ctau = -uc*f1.tau - vc*f2.tau - dtmp*(f.tau - lt*old_p) - (f.p + f.q)*(uy + vx) + f.s*(vx - uy);
    }

    // Add in chi-diffusion. This must be done separately, since it
    // requires knowledge of the previously computed and stored D_pl.
    chi_diffusion(dt);

    // Apply the updates to the fields
#pragma omp parallel for
    for(int j = 0; j < n; j++) {
        if (j>0||y_prd) {for(int i = 0; i < m; i++) fm[i + gm*j].update_others();}
        else for(int i = 0; i < m; i++) fm[i + gm*j].update_staggered();
    }

    // Update the simulation time, and set the boundaries
    time+=dt;
    set_boundaries();
}

/** Adds in contributions to chi diffusion to change in the chi field. */
void shear_sim::chi_diffusion(double dt) {
    double dfac=0.25*chi_len*chi_len*mu_inv,dfacx=dfac*xsp*xsp,dfacy=dfac*ysp*ysp;
    double lt=lamb*time,lfac=dy*xsp*lt*0.25;

#pragma omp parallel for
    for(int j=0;j<n;j++) for(int i=0;i<m;i++) {
        c_field *fp=fm+(i+gm*j),&f=*fp,&fl=fp[-1],&fr=fp[1],&fd=fp[-gm],&fu=fp[gm];

        // Compute D_chi terms
        double Dd=f.ddev+fd.ddev,Dl=f.ddev+fl.ddev,
               Dr=f.ddev+fr.ddev,Du=f.ddev+fu.ddev,

        // Compute transformed chi gradients
               Td=Dd*(f.chi-fd.chi-lfac*(fr.chi+fp[1-gm].chi-fl.chi-fp[-1-gm].chi)),
               Tl=Dl*(f.chi-fl.chi),Tr=Dr*(fr.chi-f.chi),
               Tu=Du*(fu.chi-f.chi-lfac*(fp[1+gm].chi+fr.chi-fp[-1+gm].chi-fl.chi)),

        // Compute pre-transformed diffusivities
               diffx=dfacx*(Tr-Tl),diffy=dfacy*(Tu-Td);

        // Add diffusion term
        f.cchi+=(1-lt)*diffx+diffy;
    }
}

/** Carries out a projection step of the stresses using a multigrid solve,
 * to enforce the quasistaticity constraint.
 * \param[in] dt the timestep to use. */
void shear_sim::projection_step(double dt) {
    const double third = 1./3., L = K-2*third*mu;
    double hx   = 0.5*xsp*dt, hy   = 0.5*ysp*dt,
           xfac = 0.5/dt*xsp, yfac = 0.5/dt*ysp,
           lt = lamb*time, ltsq = lt*lt;
    vec *src=qsm.b,*vel=qsm.z;
    int j;

    // Initialize constants in the algebraic system. Set up the preconditioning
    // matrix.
    qsm.init(viscosity/dt,lt);
    mat pre(qsm.regci*mat(1.,lt,0,1.));

    if (y_prd) {

        // Calculate the source term for the projection step, with no boundary
        // conditions - we turn on periodicity in the multigrid solve.
#pragma omp parallel for
        for(j = 0; j < n; j++) for(int i = 0; i < m; i++) {
            int ij = i + m*j;
            c_field *fp = fm + (i + gm*j);
            net_force(src[ij].x, src[ij].y, xfac, yfac, fp);
            src[ij] = pre*src[ij];
        }
    } else {

        // Set lower boundary condition
        vec xbc(bdry_vel(),0.);
        for(j = 0; j < m; j++) src[j] = -xbc;

        // Calculate the source term for the projection step
#pragma omp parallel for
        for(j = 1; j < n; j++) for(int i = 0; i < m; i++) {
            int ij = i + m*j;
            c_field *fp = fm + (i + gm*j);
            net_force(src[ij].x, src[ij].y, xfac, yfac, fp);
            src[ij] = pre*src[ij];
        }

        // Set upper boundary condition
        for(j = mn; j < mn + m; j++) src[j] = xbc;
    }

    // Set up the algebraic systems in the grid hierarchy and carry out the
    // multigrid solve
    qsm.solve();

    // Apply the correction to the stress tensor, based upon the computed
    // velocity. In addition, set up a representation of the velocity on
    // the primary grid.
#pragma omp parallel for
    for(j = 0; j < n; j++) {
        for(int i = 0; i < m; i++) {
            c_field &f = fm[i + gm*j];
            int d=i==m-1?1-m:1,e=j==lsn-1?m*(1-lsn):m;
            vec *vp = vel + (i + m*j);
            f.u = vp->x;
            f.v = vp->y;
            double ux = hx*(-vp->x + vp[d].x - vp[e].x + vp[d+e].x),
                   uy = hy*(-vp->x - vp[d].x + vp[e].x + vp[d+e].x),
                   vx = hx*(-vp->y + vp[d].y - vp[e].y + vp[d+e].y),
                   vy = hy*(-vp->y - vp[d].y + vp[e].y + vp[d+e].y);

            // Apply the changes to the fields
            f.p   += 2*third*lt*lamb*mu*dt + 2*third*lt*mu*uy
                   - (L + third*L*ltsq + 2*third*mu)*(vy + ux)
                   - 2*third*mu*ltsq*ux + 2*third*mu*lt*vx;;
            f.q   += third*lt*lamb*mu*dt + third*lt*mu*uy - (L*ltsq/6 + third*mu)*(vy + ux)
                   - third*ltsq*mu*ux + third*lt*mu*vx;
            f.s   += -lamb*lt*mu*dt + lt*mu*(vx - uy) + (L*ltsq/2 - mu)*vy + (L*ltsq/2 + mu + ltsq*mu)*ux;
            f.tau += lamb*mu*dt + mu*uy - lt*(L + mu)*(vy + ux) + mu*(1 + ltsq)*vx;
        }
    }
    set_boundaries();
}

/** Calculates the net force at a cell corner using the stress tensors on the
 * four neighboring cell centers.
 * \param[out] (fx,fy) the components of force.
 * \param[in] (xf,yf) prefactors to apply to the x and y derivatives that go
 *                    into the computation of force.
 * \param[in] fp a pointer to the grid point to consider. */
inline void shear_sim::net_force(double &fx, double &fy, double xf, double yf, c_field *fp) {
    c_field &f0 = fp[-gm-1], &f1 = fp[-gm], &f2 = fp[-1], &f3 = *fp;

    // Calculate the net force in the x direction
    fx = (-f0.p + f1.p - f2.p + f3.p - f0.q + f1.q - f2.q + f3.q + f0.s - f1.s + f2.s - f3.s)*xf
      + (f0.tau + f1.tau - f2.tau - f3.tau)*yf;

    // Calculate the net force in the y direction
    fy = (f0.tau - f1.tau + f2.tau - f3.tau)*xf
      + (-f0.p - f1.p + f2.p + f3.p - f0.q - f1.q + f2.q + f3.q - f0.s - f1.s + f2.s + f3.s)*yf;
}

/** Evaluates the squared L2 norm between the simulation fields and the
 * simulation fields in another instance of the shear_sim class.
 * \param[in] ss another shear_sim class.
 * \param[in] l2 a four-entry array in which to store the output. */
void shear_sim::l2_comparison(shear_sim &ss, double *l2) {
    *l2 = l2[1] = l2[2] = l2[3] = 0;

#pragma omp parallel for
    for (int j = 0; j < n; j++) {
        c_field *fp = fm + gm*j, *fo = ss.fm + gm*j, *fe = fp + m;
        double c0 = 0, c1 = 0, c2 = 0, c3 = 0;

        while (fp < fe) {

            // Compute the differences between the fields in the two simulations
            double du = fp->u - fo->u, dv = fp->v - fo->v, dp = fp->p - fo->p,
                   dq = fp->q - fo->q, ds = fp->s - fo->s, dtau = fp->tau - fo->tau,
                   dchi = fp->chi - fo->chi, dX = fp->X - fo->X, dY = fp->Y - fo->Y;

            // Compute the contributions to the L2 norms
            c0 += du*du + dv*dv;
            c1 += dtau*dtau + ds*ds + 3*dp*dp + 6*dq*dq;
            c2 += dchi*dchi;
            c3 += dX*dX + dY*dY;
            fp++; fo++;
        }

        // Scale the zeroth line of the regular fields according to the
        // trapezoidal rule
        if (j==0) {c0 *= 0.5; c3 *= 0.5;}

        // Add the totals for this line to the accumulators
#pragma omp critical
        {*l2 += c0; l2[1] += c1; l2[2] += c2; l2[3] += c3;}
    }

    // Deal with the extra line for the regular fields
    c_field *fp = fm + gmn,
            *fo = ss.fm + gmn,
            *fe = fp + m;
    while (fp < fe) {
        double du = fp->u - fo->u, dv = fp->v - fo->v,
               dX = fp->X - fo->X, dY = fp->Y - fo->Y;
        *l2   += 0.5*(du*du + dv*dv);
        l2[3] += 0.5*(dX*dX + dY*dY);
        fp++; fo++;
    }

    // Normalize by the grid size
    double fac = dx*dy;
    *l2   *= fac; l2[1] *= fac;
    l2[2] *= fac; l2[3] *= fac;
}

/** Evaluates the squared L2 norm between the simulation fields and the
 * simulation fields in another instance of the shear_sim class. It assumes
 * that THIS simulation is in the transformed frame, while the other simulation
 * is in the physical frame.
 * \param[in] ss another shear_sim class.
 * \param[in] l2 a four-entry array in which to store the output. */
void shear_sim::l2_comparison_transform(shear_sim &ss, double *l2) {
    *l2 = l2[1] = l2[2] = l2[3] = 0;

#pragma omp parallel for
    for (int j = 0; j < n; j++) {

        // Compute the physical grid point corresponding to the current grid
        // point. In this case, the y dimension is untransformed, and so j is
        // simply invariant.
        double j_phys = j;

        // Starting and ending grid points of the computation.
        c_field *fp = fm + gm*j,
                *fo = ss.fm + gm*j,
                *fe = fp + m;

        // Accumulators, and interpolation storage
        double c0 = 0, c1 = 0, c2 = 0, c3 = 0, interp_vals[7];

        // Keep track of the current X grid point.
        int i = 0;
        while (fp < fe) {

            // Convert the current X grid point to the equivalent
            // x grid point
            double i_phys = i + lamb*time*xsp*(ay + j*dy);

            // Calculate the interpolation for the physical system.
            ss.bilin_interp(i_phys, j_phys, interp_vals);

            // Unpack the interpolation values
            double phys_u   = interp_vals[0], phys_v   = interp_vals[1],
                   phys_p   = interp_vals[2], phys_q   = interp_vals[3],
                   phys_s   = interp_vals[4], phys_tau = interp_vals[5],
                   phys_chi = interp_vals[6],
                   du = fp->u - phys_u, dv = fp->v - phys_v, dp = fp->p - phys_p,
                   dq = fp->q - phys_q, ds = fp->s - phys_s, dtau = fp->tau - phys_tau,
                   dchi = fp->chi - phys_chi, dX = fp->X - fo->X, dY = fp->Y - fo->Y;

            c0 += du*du + dv*dv;
            c1 += dtau*dtau + ds*ds + 3*dp*dp + 6*dq*dq;
            c2 += dchi*dchi;
            c3 += dX*dX + dY*dY;
            fp++; fo++; i++;
        }

        // Scale the zeroth line of the regular fields according to the
        // trapezoidal rule
        if (j==0) {c0 *= 0.5; c3 *= 0.5;}

        // Add the totals for this line to the accumulators
#pragma omp critical
        {*l2 += c0; l2[1] += c1; l2[2] += c2; l2[3] += c3;}
    }

    // Deal with the extra line for the regular fields
    c_field *fp = fm + gmn, *fo = ss.fm + gmn, *fe = fp + m;

    // Physical j coordinate is on the top boundary, just like for
    // the transformed frame.
    double j_phys = n;

    // Will hold the results of the bilinear interpolation.
    double interp_vals[7];

    int i = 0;
    while (fp < fe) {
        // Convert the current X grid point to the equivalent
        // x grid point.
        double i_phys = i + lamb*time*xsp*(ay + n*dy);

        // Calculate the interpolation for the physical system.
        ss.bilin_interp(i_phys, j_phys, interp_vals);

        // Unpack the interpolation values.
        double phys_u   = interp_vals[0], phys_v   = interp_vals[1],
               du = fp->u - phys_u, dv = fp->v - phys_v,
               dX = fp->X - fo->X, dY = fp->Y - fo->Y;

        *l2   += 0.5*(du*du + dv*dv);
        l2[3] += 0.5*(dX*dX + dY*dY);
        fp++; fo++; i++;
    }

    // Normalize by the grid size
    double fac = dx*dy;
    *l2   *= fac; l2[1] *= fac;
    l2[2] *= fac; l2[3] *= fac;
}

/* Bilinearly interpolates the field values to the (potentially fractional,
 * nonexistent) grid point (i, j).
 * \param[in] (ii,jj) the coordinate position to interpolate at.
 * \param[in] results an array in which to store the results. */
void shear_sim::bilin_interp(double ii, double jj, double (&results)[7]) {

    // Determine the bounding square of grid points
    int left = static_cast<int>(ii), down = static_cast<int>(jj),
        right = left + 1, up = down + 1;

    // Get a handle on the field values
    c_field *fdl = (fm + left  + down*gm), *fdr = (fdl + 1),
            *ful = (fdl + gm), *fur = (fdr + gm);

    // Interpolate each field
    double f=right-ii, g=ii-left;
    for(int k=0; k<7; k++)
        results[k] = (up-jj)*(f*fdl->fval(k) + g*fdr->fval(k))
                   + (jj-down)*(f*ful->fval(k) + g*fur->fval(k));
}

/** Adaptively calculates the effect of plasticity on the shear stress
 * components and the effective temperature, by dividing up the timestep
 * into smaller intervals to ensure that the change in shear stress remains
 * small.
 * \param[in] sbar the magnitude of the shear stress at a grid point.
 * \param[in,out] chiv the effective temperature at the grid point.
 * \param[out] ddev a scaled measure of Dpl, needed for chi diffusion.
 * \param[in] dt the timetep.
 * \return The factor to scale the shear stresses by. */
double shear_sim::adaptive_plastic_term(double sbar,double &chiv,double &ddev,double dt) {
    const double small_dev_cutoff=1e-12;
    const double theta=0;
    const bool no_chidot=false;
    bool adapt=true;
    double dtl=dt, osbar=sbar, dplas, dchi1, dchi2, adt;

    // To avoid division-by-zero errors, check for the case when the
    // deviatoric stress is small, and just return with no changes
    if(sbar<small_dev_cutoff) return 0;
    do {

        // Calculate the plastic deformation
        dplas=2*t_scale*mu*stz->Dplastic(sbar,chiv,dchi1,dchi2);

        // If the plastic deformation is large, then adaptively
        // choose the timestep so as not to make a large alteration in
        // the deviatoric stress magnitude
        if(fabs(dtl*dplas)>adapt_fac) {
            adt=adapt_fac/fabs(dplas);
            dtl-=adt;
        } else {
            adapt=false;
            adt=dtl;
        }
        sbar-=adt*dplas;

        // Calculate the change in chi if needed
        if(!no_chidot) chiv+=(dchi1*(chi_inf-chiv)+dchi2*(theta-chiv))*t_scale*adt;
    } while(adapt);

    // Calculate the quantity dt*2*mu*Dpl/sbar, which will be used in the
    // finite-difference update
    ddev=osbar-sbar;
    return ddev/osbar;
}

/** Sets the fields in the ghost regions according to the boundary conditions.
 */
void shear_sim::set_boundaries() {

    // Set the top and bottom boundaries
    const int g=n*gm;
    if (y_prd) {
        for (c_field *fp=fm; fp<fm+m; fp++) {
            fp[-2*gm].prd_bc(fp[g-2*gm],0,-ly);
            fp[-gm].prd_bc(fp[g-gm],0,-ly);
            fp[g].prd_bc(*fp,0,ly);
            fp[g+gm].prd_bc(fp[gm],0,ly);
        }
    } else {
        c_field *fp=fm;
        double U=bdry_vel(),bp=bdry_pos();
        for(int i=0; i<m; i++, fp++) {
            double xl=ax+bp+i*dx,xu=ax-bp+i*dx;
            fp->set_corner(-U,0,xl,ay);fp[-gm].extrap(*fp,fp[gm]);
            fp[-gm].set_corner(-U,0,xl,ay-dy);fp[-2*gm].extrap(fp[-gm],*fp);

            fp[g].set_corner(U,0,xu,by);fp[g].extrap(fp[g-gm],fp[g-2*gm]);
            fp[g+gm].set_corner(U,0,xu,by+dy);fp[g+gm].extrap(fp[g],fp[g-gm]);
        }
    }

    // Loop over the "ghost columns" and set them equal to the corresponding
    // points from the other end of the row.
    for(c_field *fp=fbase+2; fp<fbase+gmn; fp+=gm) {
        fp[-2].prd_bc(fp[m-2],-lx,0);
        fp[-1].prd_bc(fp[m-1],-lx,0);
        fp[m].prd_bc(*fp,lx,0);
        fp[m+1].prd_bc(fp[1],lx,0);
    }
}

/** Sets the fields in the ghost regions according to the boundary conditions.
 */
void shear_sim::set_strain_boundaries() {

    // Set the top and bottom boundaries
    const int g=n*gm;
    if (y_prd) {
        for (c_field *fp=fm; fp<fm+m; fp++) {
            fp[-2*gm].strain_prd_bc(fp[g-2*gm]);
            fp[-gm].strain_prd_bc(fp[g-gm]);
            fp[g].strain_prd_bc(*fp);
            fp[g+gm].strain_prd_bc(fp[gm]);
        }
    } else {
        c_field *fp=fm;
        for(int i=0; i<m; i++, fp++) {
            fp[-gm].strain_extrap(*fp,fp[gm]);
            fp[-2*gm].strain_extrap(fp[-gm],*fp);
            fp[g].strain_extrap(fp[g-gm],fp[g-2*gm]);
            fp[g+gm].strain_extrap(fp[g],fp[g-gm]);
        }
    }

    // Loop over the "ghost columns" and set them equal to the corresponding
    // points from the other end of the row.
    for(c_field *fp=fbase+2; fp<fbase+gmn; fp+=gm) {
        fp[-2].strain_prd_bc(fp[m-2]);
        fp[-1].strain_prd_bc(fp[m-1]);
        fp[m].strain_prd_bc(*fp);
        fp[m+1].strain_prd_bc(fp[1]);
    }
}

/** Calculates one-sided derivatives of the cell-cornered fields using the
 * second-order ENO2 scheme.
 * \param[out] (ud,vd) the computed velocity derivatives.
 * \param[out] (Xd,Yd) the computed reference map derivatives.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] fp a pointer to the field at which to evaluate the derivatives.
 * \param[in] d the step between the grid points to use. */
inline void shear_sim::rmv_eno2(double &ud, double &vd, double &Xd, double &Yd, double hs, c_field *fp, int d) {
    ud = hs*eno2(fp[-d].u, fp->u, fp[d].u, fp[2*d].u);
    vd = hs*eno2(fp[-d].v, fp->v, fp[d].v, fp[2*d].v);
    rm_eno2(Xd, Yd, hs, fp, d);
}

/** Calculates one-sided derivatives of the reference map fields using the
 * second-order ENO2 scheme.
 * \param[out] (ud,vd) the computed velocity derivatives.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] fp a pointer to the field at which to evaluate the derivatives.
 * \param[in] d the step between the grid points to use. */
inline void shear_sim::rm_eno2(double &Xd, double &Yd, double hs, c_field *fp, int d) {
    Xd = hs*eno2(fp[-d].X, fp->X, fp[d].X, fp[2*d].X);
    Yd = hs*eno2(fp[-d].Y, fp->Y, fp[d].Y, fp[2*d].Y);
}

/** Calculates one-sided derivatives of the cell-centered fields using the
 * second-order ENO2 scheme.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] fp a pointer to the field at which to evaluate the derivatives.
 * \param[in] d the step between the grid points to use.
 * \return The computed derivatives. */
inline st_field shear_sim::st_eno2(double hs,c_field *fp,int d) {
    c_field &f0=fp[-d],&f1=*fp,&f2=fp[d],&f3=fp[2*d];
    return st_field(hs*eno2(f0.p,f1.p,f2.p,f3.p),hs*eno2(f0.q,f1.q,f2.q,f3.q),
        hs*eno2(f0.s,f1.s,f2.s,f3.s),hs*eno2(f0.tau,f1.tau,f2.tau,f3.tau),
        hs*eno2(f0.chi,f1.chi,f2.chi,f3.chi));
}

/** Calculates the ENO derivative using a sequence of values at
 * four gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double shear_sim::eno2(double p0,double p1,double p2,double p3) {
    return fabs(p0-2*p1+p2)>fabs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void shear_sim::write_files(int k) {

    // Output the main simulation fields
    if(fflags&1) output("u",0,k);
    if(fflags&2) output("v",1,k);
    if(fflags&4) output("p",2,k);
    if(fflags&8) output("q",3,k);
    if(fflags&16) output("s",4,k);
    if(fflags&32) output("tau",5,k);
    if(fflags&64) output("chi",6,k);
    if(fflags&128) output("tem",7,k);
    if(fflags&256||true) output("dev",8,k);
    if(fflags&512) output("X",9,k);
    if(fflags&1024) output("Y",10,k);

    // If any other type of output is requested, it will need strain, so
    // calculate that now
    if(fflags&14336) compute_strain();

    // Output the components of strain
    if(fflags&2048) {
        output("Exx",11,k);
        output("Exy",12,k);
        output("Eyy",13,k);
    }

    // Output the invariants of strain
    if(fflags&4096) {
        output("E_I",14,k);
        output("E_II",15,k);
        output("Edev",16,k);
    }

    // Output the tracers and chi field pulled back to the reference state
    if((fflags&8192)&&tr!=NULL) {
        if(k>0) update_tracers();
        output_tracers(k);
        output_tracers_matrix("teml",4,k);
    }
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename. */
void shear_sim::output(const char *prefix,const int mode,const int sn) {

    // Determine whether to output a staggered field or not
    bool st = (mode >= 2) && (mode <= 8);

    // Number of grid points in the x direction.
    int l = st? m : y_prd? m : m + 1;

    // Assemble the output filename and open the output file.
    char *bufc = reinterpret_cast<char*>(buf);
    sprintf(bufc, "%s/%s.%d", filename, prefix, sn);
    FILE *outf = safe_fopen(bufc, "wb");

    // Output the first line of the file.
    // See gnuplot.sourceforge.net/docs_4.2/node330.html
    int i, j;
    float *bp = buf + 1, *be = buf + m + 1;

    // Set the first element of the first row to be the number of x points. The
    // rest of the first row is all of the x coordinates.
    *buf = l;
    for(i = 0; i < l; i++) *(bp++) = ax + (st? (i + 0.5)*dx : i*dx);
    fwrite(buf, sizeof(float), l+1, outf);

    // Output the field values to the file.
    c_field *fp = fm;

    // In the non y-periodic case, there is one more nonstaggered field than staggered
    // field, so that the top boundary can be included. In the y-periodic case, there
    // is the same number of staggered and nonstaggered fields.
    int max_y = st? n : y_prd? n : n + 1;
    for(j = 0; j < max_y; j++) {

        // First element in every row for binary matrix format is the y coordinate
        *buf = ay + (st? (j + 0.5)*dy : j*dy); bp = buf + 1;
        switch(mode) {
            case  0: while(bp < be) *(bp++) = (fp++)->u;      break;
            case  1: while(bp < be) *(bp++) = (fp++)->v;      break;
            case  2: while(bp < be) *(bp++) = (fp++)->p;      break;
            case  3: while(bp < be) *(bp++) = (fp++)->q;      break;
            case  4: while(bp < be) *(bp++) = (fp++)->s;      break;
            case  5: while(bp < be) *(bp++) = (fp++)->tau;    break;
            case  6: while(bp < be) *(bp++) = (fp++)->chi;    break;
            case  7: while(bp < be) *(bp++) = (fp++)->chi*TZ; break;
            case  8: while(bp < be) *(bp++) = (fp++)->dev(lamb*time); break;
            case  9: while(bp < be) *(bp++) = (fp++)->X;      break;
            case 10: while(bp < be) *(bp++) = (fp++)->Y;      break;
            case 11: while(bp < be) *(bp++) = (fp++)->cu;     break;
            case 12: while(bp < be) *(bp++) = (fp++)->cv;     break;
            case 13: while(bp < be) *(bp++) = (fp++)->cp;     break;
            case 14: while(bp < be) *(bp++) = (fp++)->cq;     break;
            case 15: while(bp < be) *(bp++) = (fp++)->cs;     break;
            case 16: while(bp < be) *(bp++) = (fp++)->ctau;
        }

        // Pass over the ghost points
        fp += 4;

        // Something with the reference map
        if (!st) *bp = (mode == 9)? buf[1] + bx - ax : buf[1];

        // Write it out
        fwrite(buf, sizeof(float), l+1, outf);
    }

    // Close the file
    fclose(outf);
}

/** Calculates the Green-Lagrange strain and stores its components and
 * invariants in temporary values with the main data structure. */
void shear_sim::compute_strain() {

    // Create transformation matrix
    mat T=mat(1,lamb*time,0,1);

#pragma omp parallel
    {
        double J,Xx,Xy,Yx,Yy,tfac;
        mat F;sym_mat E;

        // Loop over the main grid
#pragma omp for
        for(int j=0;j<n;j++) {
            for(c_field *fp=fm+j*gm,*fe=fp+m;fp<fe;fp++) {

                // Calculate the Jacobian of the reference map
                Xy=0.5*ysp*(fp[gm].X-fp[-gm].X);
                Yy=0.5*ysp*(fp[gm].Y-fp[-gm].Y);
                Xx=0.5*xsp*(fp[1].X-fp[-1].X);
                Yx=0.5*xsp*(fp[1].Y-fp[-1].Y);

                // Compute the deformation gradient tensor. Pre-multiply by the
                // transformation matrix to obtain the physical deformation
                // gradient.
                J=1/(Xx*Yy-Xy*Yx);
                F=T*mat(Yy*J,-Xy*J,-Yx*J,Xx*J);

                // Compute the Green-Lagrange strain tensor E
                E=F.ATA();
                E.a-=1.;E.d-=1.;E*=0.5;

                // Store components of E
                fp->cu=E.a;fp->cv=E.b;fp->cp=E.d;

                // Store invariants of E (interpreted as a 3x3 tensor):
                // 1. I_A=tr(E),
                // 2. II_A=0.5*(tr(E)^2-tr(E^2)),
                // 3. ||E-1/3.*tr(E)*I||_F.           (Frobenius norm)
                fp->cq=E.trace();
                fp->cs=E.invariant2();
                tfac=1/3.*fp->cq;
                E.a-=tfac;E.d-=tfac;
                fp->ctau=sqrt(E.mod_sq()+tfac*tfac);
            }
        }
    }

    // Set the strain values in the ghost regions. This is needes to ensure the
    // boundary values are correct on output, and to ensure that bicubic
    // interpolation of strain works correctly.
    set_strain_boundaries();
}

/** Initializes the effective temperature field as a bicubic interpolation of a
 * given array, scaling the array values.
 * \param[in] (mm,nn) the dimensions of the given array.
 * \param[in] arr a pointer to the array.
 * \param[in] lo an additive constant to scale by.
 * \param[in] sca an multiplicative constant to scale by. */
void shear_sim::initialize_chi_bicubic(int mm,int nn,double *arr,double lo,double sca) {
    int mmnn=mm*nn;
    double *xd=new double[3*mmnn],*yd=xd+mmnn,*xyd=yd+mmnn;

    // Calculate x derivative of the given array
#pragma omp parallel for
    for(int j=0;j<nn;j++) {
        int ij=j*mm;
        for(int i=0;i<mm;i++,ij++) xd[ij]=0.5*(arr[i==mm-1?ij-mm+1:ij+1]-arr[i==0?ij+mm-1:ij-1]);
    }

    // Calculate y and xy derivatives of the given array
#pragma omp parallel for
    for(int j=0;j<nn;j++) {
        int ij=mm*j,ijd,iju;
        double fac;

        // Set rows to consider and the prefactor in the derivative
        ijd=j==0?(y_prd?ij+(mmnn-mm):ij):ij-mm;
        iju=j==nn-1?(y_prd?ij+(mm-mmnn):ij):ij+mm;
        fac=!y_prd&&j==0&&j==mm-1?1.:0.5;

        // Compute the derivatives for this row
        for(int ije=ij+mm;ij<ije;ij++,ijd++,iju++) {
            yd[ij]=fac*(arr[iju]-arr[ijd]);
            xyd[ij]=fac*(xd[iju]-xd[ijd]);
        }
    }

    // Assemble the effective temperature using bicubic interpolation
    double xsca=dx*mm/lx,ysca=dy*nn/ly;
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double Y=(j+0.5)*ysca-0.5,mY;
        int jj,jjhi,d;
        if(y_prd) {
            if(Y<0) Y+=nn;
            jjhi=nn-1;
        } else jjhi=nn-2;
        jj=static_cast<int>(Y);
        if(jj<0) jj=0;else if(jj>jjhi) jj=jjhi;
        d=jj==nn-1?mm-mmnn:mm;
        Y-=jj;mY=1.-Y;
        double b0=Y*mY*mY,b1=mY*mY*(1+2*Y),b2=Y*Y*(3-2*Y),b3=-Y*Y*mY;
        for(int i=0;i<m;i++) {
            double X=(i+0.5)*xsca-0.5,mX;
            if(X<0) X+=mm;
            int ii=static_cast<int>(X),ij,ijp;
            if(ii>=mm) ii-=mm;
            ij=ii+mm*jj;ijp=(ii==mm-1?ii-mm:ii)+1+mm*jj;
            X-=ii;mX=1-X;
            double c0=X*mX*mX,c1=mX*mX*(1+2*X),c2=X*X*(3-2*X),c3=-X*X*mX;
            fm[i+gm*j].chi=lo+sca*(b0*(xyd[ij]*c0+yd[ij]*c1+yd[ijp]*c2+xyd[ijp]*c3)
                +b1*(xd[ij]*c0+arr[ij]*c1+arr[ijp]*c2+xd[ijp]*c3)
                +b2*(xd[ij+d]*c0+arr[ij+d]*c1+arr[ijp+d]*c2+xd[ijp+d]*c3)
                +b3*(xyd[ij+d]*c0+yd[ij+d]*c1+yd[ijp+d]*c2+xyd[ijp+d]*c3));
        }
    }

    set_boundaries();

    // Delete the temporary space
    delete [] xd;
}

/** Calculates the normalizing factors for the gaussian blurring.
 * \param[in] l the length scale of the gaussian, measured in gridpoints.
 * \param[in] cut the number of gridpoints to cut the PDF off at.
 * \return The normalization factor for the variance. */
double shear_sim::gaussian_normal_factor(double llinv,int cut) {
    double vfac=1,q,q2;
    for(int j=1;j<cut;j++) {
        q=exp(-j*j*llinv);
        for(int k=0;k<=j;k++) {
            q2=q*exp(-k*k*llinv);
            vfac+=(k==0||k==j?4:8)*q2*q2;
        }
    }
    printf("%g\n",vfac);
    return sqrt(vfac);
}

/** Initializes the simulation fields, setting the effective temperature field
 * to be a smoothed
 * \param[in] l the smoothing length scale, measured in gridpoints.
 * \param[in] tem_avg the average effective temperature.
 * \param[in] tem_stdev the standard deviation of the effective temperature. */
void shear_sim::initialize_random(double l,double tem_avg,double tem_stdev) {

    // Calculate cutoff and normalizing Gaussian factors
    int cut=int(l*gaussian_fac),nx=y_prd?n:n+(cut<<1);
    double llinv=0.5/(l*l),chi0=tem_avg/TZ,chi1=tem_stdev/(gaussian_normal_factor(llinv,cut)*TZ);

    // Create grid of Gaussian noise
    double *no=new double[m*nx],*np=no,*ne=no+m*nx;
    for(;np<ne-1;np+=2) box_muller(*np,np[1]);
    if(np==ne-1) {
        double r1;
        box_muller(*np,r1);
    }

    // Initialize the simulation fields
    np=y_prd?no:no+m*cut;
#pragma omp parallel for
    for(int j=0;j<=n;j++) {
        double Y=ay+j*dy;
        c_field *fp=fm+gm*j;
        for(int i=0;i<m;i++,fp++) {
            c_field &f=*fp;
            double gau=0,*npy;

            // Set stresses and velocities to zero
            f.u=f.v=0;

            // Set reference map field to match initial coordinate
            f.X=ax+i*dx;f.Y=Y;

            // Calculate chi as a Gaussian smoothing of the random field
            if(j<n) {
                for(int cj=j-cut;cj<=j+cut;cj++) {
                    npy=np+(y_prd?step_mod(cj,n):cj)*m;
                    for(int ci=i-cut;ci<=i+cut;ci++)
                        gau+=npy[step_mod(ci,m)]*exp(-llinv*((cj-j)*(cj-j)+(ci-i)*(ci-i)));
                }
                f.chi=chi0+gau*chi1;
                f.p=f.q=f.s=f.tau=0;
            }
        }
    }

    // Remove Gaussian noise grid
    delete [] no;
}

/** Calculates two normally distributed random numbers.
 * \param[out] (n0,n1) two normally distributed numbers. */
void shear_sim::box_muller(double &r0,double &r1) {
    double x,y,r;
    do {
        x=2.0/RAND_MAX*rand()-1;
        y=2.0/RAND_MAX*rand()-1;
        r=x*x+y*y;
    } while(r==0.0||r>1.0);
    double d=sqrt(-2.0*log(r)/r);
    r0=x*d;
    r1=y*d;
}

/** Calculates the tractions on the top and bottom surfaces.
 * \param[out] (bx,by) the components of the traction on the bottom surface.
 * \param[out] (tx,ty) the components of the traction on the top surface. */
void shear_sim::tractions(double &bx,double &by,double &tx,double &ty) {

    bx=by=tx=ty=0;

    // Calculate bottom traction, extrapolating the cell-centered stresses
    // to the boundary
    c_field *fp=fm,*fe=fp+m;
    while(fp<fe) {
        bx+=3*fp->tau-fp[m].tau;
        by+=3*fp->syy()-fp[m].syy();
        fp++;
    }

    // Calculate top traction
    fp=fm+(n-2)*m,fe=fp+m;
    while(fp<fe) {
        tx-=3*fp[m].tau-fp->tau;
        ty-=3*fp[m].syy()-fp->syy();
        fp++;
    }

    // Normalize by grid spacing, also taking into account factor of two in
    // extrapolation formula
    bx*=0.5*dx;by*=0.5*dx;
    tx*=0.5*dx;ty*=0.5*dx;
}

/** Initializes a grid of Lagrangian tracers.
 * \param[in] (mm,nn) the. */
void shear_sim::initialize_tracers(int mm,int nn) {

    // Remove any previous tracers, and allocate memory. Set up bicubic
    // interpolation for X, Y, pressure, dev(), chi if needed.
    if(tr==NULL) {
        *bi=new bicubic_interp(*this,7);
        bi[1]=new bicubic_interp(*this,8);
        bi[2]=new bicubic_interp(*this,2);
        bi[3]=new bicubic_interp(*this,6);
        bi[4]=new bicubic_interp(*this,9);
        bi[5]=new bicubic_interp(*this,10);
        bi[6]=new bicubic_interp(*this,11);
        bi[7]=new bicubic_interp(*this,12);
        bi[8]=new bicubic_interp(*this,13);
        bi[9]=new bicubic_interp(*this,14);
    } else delete [] tr;
    ntrace=mm*nn;
    tr=new double[ntrace*6];

    // Initialize tracers on a regular grid
    double yy,dxx=(bx-ax)/mm,dyy=(by-ay)/nn,*trp=tr;
    for(int j=0;j<nn;j++) {
        yy=ay+dyy*(j+0.5);
        for(int i=0;i<mm;i++) {

            // Reference point of tracer
            *(trp++)=ax+dxx*(i+0.5);*(trp++)=yy;

            // Current position of tracer
            *(trp++)=ax+dxx*(i+0.5);*(trp++)=yy;

            // Velocity estimate of tracer
            *(trp++)=0;*(trp++)=0;
        }
    }
    tr_time=time;tr_m=mm;tr_n=nn;
}

/** Update tracer positions using the reference map. */
void shear_sim::update_tracers() {
    double dt=time-tr_time,dtinv=1/dt;
    double xref,yref,xx,yy,*trp=tr,*tre=tr+6*ntrace;
    while(trp<tre) {
        xref=*trp;yref=trp[1];

        // Estimate the tracer's position using the velocity estimate
        xx=trp[2]+dt*trp[4];
        yy=trp[3]+dt*trp[5];

        // Perform Newton step to correct the tracer position
        if(correct_tracer(xref,yref,xx,yy)) {

            // If the Newton step was successful, then update the
            // tracer position the velocity estimate
            trp[4]=(xx-trp[2])*dtinv;
            trp[5]=(yy-trp[3])*dtinv;
            trp[2]=xx;
            trp[3]=yy;
            trp+=6;
        } else {

            // If the Newton step failed (which will probably be
            // very rare) then delete the tracer from memory
            ntrace--;tre-=6;
            *trp=*tre;trp[1]=tre[1];trp[2]=tre[2];
            trp[3]=tre[3];trp[4]=tre[4];trp[5]=tre[5];
            fprintf(stderr,"# Lost tracer at %g %g\n",trp[2]+dt*trp[4],trp[3]+dt*trp[5]);
        }
    }
}

/** Corrects tracer position to be anchored to the reference map field.
 * \param[in] (xref,yref) the reference map position to anchor to.
 * \param[in,out] (xx,yy) the tracer position to be corrected.
 * \return True if the anchoring was successful, false otherwise. */
bool shear_sim::correct_tracer(double xref,double yref,double &xx,double &yy) {
    const double tol=std::numeric_limits<double>::epsilon(),tolnewt=tol*tol*1e6;
    const double ilx=1/lx,ily=1/ly,bigstep=lx*lx+ly*ly;
    double fx,fy,fxx,fxy,fyx,fyy,det,delx,dely,xs,ys;

    int count=0,iter=0;
    while(count<4) {

        // Check for too many iterations
        if(iter>100) return false;

        // Compute function values and Jacobian
        xs=floor((xx-ax)*ilx)*lx;
        ys=floor((yy-ay)*ily)*ly;
        fx=bi[0]->f_grad_f(xx-xs,yy-ys,fxx,fxy)-xref+xs;
        fy=bi[1]->f_grad_f(xx-xs,yy-ys,fyx,fyy)-yref+ys;

        // Check for convergence
        if(fx*fx+fy*fy<tolnewt) count++;

        // Bail out if the determinant is within the machine epsilon
        det=fxx*fyy-fxy*fyx;
        if(det<tol&&det>-tol) return false;

        // Compute update and bail out if there's a huge step
        delx=( fyy*fx-fxy*fy)/det;
        dely=(-fyx*fx+fxx*fy)/det;
        if(delx*delx+dely*dely>bigstep) return false;

        // Apply updates
        xx-=delx;yy-=dely;iter++;
    }
    return true;
}

/** Outputs the tracer positions, and bicubic interpolations of several
 * different fields.
 * \param[in] sn the suffix of the tracer file. */
void shear_sim::output_tracers(int sn) {
    const double ilx=1/lx,ily=1/ly;

    // Assemble the output filename and open the output file
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/tra.%d",filename,sn);
    FILE *outf=safe_fopen(bufc,"w");

    // Loop over the tracers and output the interpolated fields at their
    // positions
    double xx,yy,xw,yw,chi;
    for(double *trp=tr,*tre=tr+6*ntrace;trp<tre;trp+=6) {
        xx=trp[2];yy=trp[3];
        xw=xx-floor((xx-ax)*ilx)*lx;
        yw=yy-floor((yy-ay)*ily)*ly;
        chi=bi[3]->f(xw,yw);
        fprintf(outf,"%.12g %.12g %.12g %.12g %.12g %.12g %.12g %.10g %.10g %.10g %.10g %.10g %.10g\n",xx,yy,
            *trp,trp[1],bi[2]->f(xw,yw),chi,chi*TZ,bi[4]->f(xw,yw),bi[5]->f(xw,yw),bi[6]->f(xw,yw),
            bi[7]->f(xw,yw),bi[8]->f(xw,yw),bi[9]->f(xw,yw));
    }

    // Close the file
    fclose(outf);
}

/** Outputs a 2D field in original Lagrangian frame to file.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename. */
void shear_sim::output_tracers_matrix(const char *prefix,const int mode,const int sn) {
    const double ilx=1/lx,ily=1/ly;

    // Check for possible errors
    if(tr_m>m+1) fatal_error("Output buffer not large enough",1);
    if(ntrace!=tr_m*tr_n) fatal_error("Can't do matrix output after tracers are lost",1);

    // Assemble the output filename and open the output file
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the first line of the file
    int i,j;
    float *bp=buf+1,*be=buf+tr_m+1;
    double delx=lx/tr_m,dely=(by-ay)/tr_n,*trp=tr;
    *buf=tr_m;
    for(i=0;i<tr_m;i++) *(bp++)=ax+(i+0.5)*delx;
    fwrite(buf,sizeof(float),tr_m+1,outf);

    // Output the field values to the file. Cases 0 and 1: the current x
    // and y position of this Lagrangian tracer. Cases (2,3,4,5) correspond
    // to (p,chi,tem), respectively.
    for(j=0;j<tr_n;j++) {
        *buf=ay+(j+0.5)*dely;bp=buf+1;
        switch(mode) {
            case 0: while(bp<be) {*(bp++)=trp[2];trp+=6;} break;
            case 1: while(bp<be) {*(bp++)=trp[3];trp+=6;} break;
            default: while(bp<be) {
                double xx=trp[2],xw=xx-floor((xx-ax)*ilx)*lx,
                       yy=trp[3],yw=yy-floor((yy-ay)*ily)*ly;
                *(bp++)=mode<4?bi[mode]->f(xw,yw):bi[3]->f(xw,yw)*TZ;
                trp+=6;
            }
        }
        fwrite(buf,sizeof(float),tr_m+1,outf);
    }

    // Close the file
    fclose(outf);
}
