#include <cstring>
#include "common.hh"
#include "shear_sim.hh"
#include "mat3.hh"

using std::isnan;
using std::abs;
using std::string;

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
              const double ay_, const double by_, const double mu_, const double K_,
              const double visc_, const double t_scale_, const double adapt_fac_, const double v0_,
              stz_dynamics *stz_, const unsigned int fflags_, const char *filename_)
    : m(m_), gm(m+4), n(n_), lsn(n+1), gn(n+4), mn(m*n), gmn(gm*gn),
    fflags(fflags_), ax(ax_), ay(ay_), bx(bx_), by(by_), lx(bx-ax), ly(by-ay),
    dx(lx/m), dy(ly/n), xsp(1./dx), ysp(1./dy), mu(mu_), mu_inv(1./mu_), K(K_),
    viscosity(visc_), t_scale(t_scale_), adapt_fac(adapt_fac_), v0(v0_),
    stz(stz_), filename(filename_), fbase(new c_field[gmn]), fm(fbase+2*gm+2), time(0.),
    f_num(0), qsm(*this), bdrys_g(new bdry_field[gm]), bdrys(bdrys_g + 2), buf(new float[m>=63?m+6:64]) {

    // Set STZ related parameters
    TZ=stz->TZ;
    chi_inf=stz->chi_inf;
}

/** The class destructor frees the dynamically allocated memory. */
shear_sim::~shear_sim() {
    delete [] buf;
    delete [] fbase;
    delete [] bdrys_g;
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
    for(int j = 0; j < n+1; j++) {

        // Calculate the cell-cornered and cell-centered y-reference maps
        double Y  = ay + j*dy, 
               Ys = Y + 0.5*dy;

        // Pointer to the first element in row j
        c_field *fp = fm + gm*j;

        // Loop over row j
        for(int i = 0; i < m+1; i++, fp++) {

            // Reference to the current grid location
            c_field &f = *fp;

            // Calculate the cell-cornered and cell-centered x-reference maps
            double X  = ax + i*dx, 
                   Xs = X + 0.5*dx;

            // Transformed velocities and stresses start at simply zero.
            f.u   = 0; f.v = 0;
            f.s11 = f.s12 = f.s22 = f.s33 = 0;

            // Zeroth case: a small Gaussian blip in chi at the
            // origin
            if (chi_case == 0) {
                f.chi = chi0 + chi1*exp(-200*(Xs*Xs+Ys*Ys));
            }
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

    // Set the boundary profiles.
    bdry_field *bdp = bdrys_g;
    for (int i = 0; i < gm; i++, bdp++) { bdp->tp = by; bdp->bp = ay; }

    // Handle the boundary conditions.
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
                 visc_cond = 0.25*dx*dx/viscosity,
                 dt        = (cfl < visc_cond)? cfl : visc_cond;

    printf("CFL: %g\n", dt/dx);
    printf("Viscosity Condition: %g\n", viscosity*dt*(1./(xsp*xsp + ysp*ysp)));
    const int l = static_cast<int>(time_interval/dt) + 1;
    const double adt = time_interval/l;

    // Output the initial fields and record initial time
    write_files(0);
    puts("# Output frame 0");
    t0 = wtime();

    for(int k = 1; k <= frames; k++) {

        // Perform the explicit timestep updates, and then check the grid
        // for any strange behavior
        for(int j = 0; j < l; j++) step_forward(adt);
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
    for (int jj = -2; jj < n+2; jj++)
        for (int ii = -2; ii < m+2; ii++){
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
    double hx   = 0.5*xsp*dt, hy = 0.5*ysp*dt,
           hxxv = 2*xsp*hx*viscosity, hyyv = 2*ysp*hy*viscosity,
           hcv  = 2*(hxxv + hyyv),
           L_val   = calc_L(),
           Ld_val  = calc_Ldot(),
           Ldd_val = calc_Lddot();

    update_bdry_profiles(dt, hx);

#pragma omp parallel for
    for(int j = 0; j < n; j++) {
        for(int i = 0; i < m; i++) {
            st_field f1, f2;
            double x_val = ax + dx*i; // Current X position (cell corner).
            double y_val = ay + dy*j; // Current Y position (cell corner).

            // Create references to the fields in the neighboring gridpoints,
            c_field *fp  = fm + (gm*j + i), 
                    *flp = fp - 1, 
                    *frp = fp + 1;

            c_field &f  = *fp, 
                    &fl = *flp, 
                    &fr = *frp, 
                    &fd = fp[-gm], 
                    &fu = fp[gm];

            bdry_field *bdry_p  = bdrys  + i, // Current point on the boundary.
                       *bdry_pl = bdry_p - 1, // One to the left.
                       *bdry_pr = bdry_p + 1; // One to the right.

            // uc, vc hold the central velocities
            // ux, vx will hold the upwinded x-derivatives of the velocities
            // uy, vy will hold the upwinded y-derivatives of the velocities
            // uvisc, vvisc will hold the viscous smoothing terms
            // and dtmp will hold the result of the adaptive_plastic_term
            // algorithm, which is (2*mu*Dpl*dt/sbar)
            double uc, vc, ux, vx, uy, vy, uvisc, vvisc, dtmp;

            // Wx = dW_dX, Cx = dC_dX, likewise for Wxx, Cxx
            // W_advx, C_advx, = x derivatives of Wadv and Cadv
            // W/Cadvadv = second advective derivatives of W and C.
            // Wr, W, Wl = W values to the right, center, and left.
            // sxx, sxy are the XX and XY components of Sigma at cell corner.
            double Wx, Cx, sxx, sxy, Wxx, Cxx,
                   W_advx, C_advx, W_advadv, C_advadv,
                   Wr, W, Wl, Cr, C, Cl;
            double Xx, Xy, Yx, Yy, fx, fy;

            /**** Velocity and reference map updates. ****/
            uc = f.u; 
            vc = f.v;

            // Calculate ENO derivatives of regular fields
            (uc > 0)? rmv_eno2(ux, vx, Xx, Yx, hx, fp,  -1) : rmv_eno2(ux, vx, Xx, Yx, -hx, fp,  1);
            (vc > 0)? rmv_eno2(uy, vy, Xy, Yy, hy, fp, -gm) : rmv_eno2(uy, vy, Xy, Yy, -hy, fp, gm);

            // Calculate ENO derivatives of the advective derivatives of the boundary terms.
            (uc > 0)? bdry_adv_eno2(W_advx, C_advx, hx, bdry_p, -1) : bdry_adv_eno2(W_advx, C_advx, -hx, bdry_p, 1);

            // Compute the advective derivatives of W and C (no y-dependence).
            // TODO: Note this is missing the \p/\p t term! How do we even compute that?
            W_advadv = uc*W_advx;
            C_advadv = uc*C_advx;

            // Calculate net force due to stress imbalance.
            // fx = d_sxx_dx + d_sxy_dy
            // fy = d_sxy_dx + d_syy_dy
            net_force(fx, fy, hx, hy, fp);

            //printf("Net force: %d %d %d %g %g\n", curr_step, i, j, fx, fy);

            // Calculate viscous smoothing term.
            // TODO: Transform the Laplacian to impose smoothing in the physical domain?
            // TODO: Check with Chris if we really need to do this.
            uvisc = hxxv*(fr.u + fl.u) + hyyv*(fu.u + fd.u) - hcv*uc;
            vvisc = hxxv*(fr.v + fl.v) + hyyv*(fu.v + fd.v) - hcv*vc;
            
            //printf("Viscous terms: %d %d %d %g %g\n", curr_step, i, j, uvisc, vvisc);

            // Compute the x derivative of W.
            W   = bdry_p->comp_W();  C  = bdry_p->comp_C();
            Wr  = bdry_pr->comp_W(); Cr = bdry_pr->comp_C();
            Wl  = bdry_pl->comp_W(); Cl = bdry_pl->comp_C();
            Wx  = xsp*(Wr - W);
            Wxx = xsp*xsp*(Wr - 2*W + Wl);
            Cx  = xsp*(Cr - C);
            Cxx = xsp*xsp*(Cr - 2*C + Cl);

            //printf("cs: %d i: %d j: %d W: %g Wr: %g Wl: %g Wx: %g Wxx: %g Cx: %g Cxx:%g\n", curr_step, i, j, W, Wr, Wl, Wx, Wxx, Cx, Cxx);

            // Compute the required stress components at the current corner.
            sxx = .25*(f.s11 + fl.s11 + fd.s11 + flp[-gm].s11);
            sxy = .25*(f.s12 + fl.s12 + fd.s12 + flp[-gm].s12);

            // Apply the update to velocity. Recall that in dimensionless units mu = rho.
            f.cu = -uc*ux - vc*uy + mu_inv*(fx + uvisc + dt*Wx/W*sxx) - dt/L_val*(Ldd_val*x_val + 2*Ld_val*uc);
            f.cv = -uc*vx - vc*vy + mu_inv*vvisc/W + (dt*mu_inv*(3*sxy*Wx + sxx*Wx*(Cx + y_val*Wx)/W + sxx*(Cxx + y_val*Wxx))
                    + mu_inv*(W*fy + (Cx + y_val*Wx)*fx) - (W_advadv*y_val + 2*dt*bdry_p->comp_Wadv()*vc + C_advadv) )/W;

            //printf("Velocity updates: %d %d %d %g %g\n", curr_step, i, j, f.cu, f.cv);

            // Calculate update to reference map fields. Simple advective
            // equation dX/dt = dY/dt = 0.
            f.cX = -uc*Xx - vc*Xy;
            f.cY = -uc*Yx - vc*Yy;

            // Calculate velocity and its derivatives at cell center.
            uc = 0.25*(f.u + fr.u + fu.u + frp[gm].u);
            vc = 0.25*(f.v + fr.v + fu.v + frp[gm].v);
            ux = .5*xsp*(-f.u + fr.u - fu.u + frp[gm].u);

            /*** Stress update. ***/
            // Calculate matrices and everything we need for them.
            mat3 Sigma        = f.Sigma();
            
            // These values are computed at cell centers, but because there's no y-dependence, they look like cell corners.
            double dCdot_dX   = xsp*(bdry_pr->comp_Cadv() - bdry_p->comp_Cadv());
            double dC_dX_dot  = dCdot_dX - ux*Cx;
            double dWdot_dX   = xsp*(bdry_pr->comp_Wadv() - bdry_p->comp_Wadv());
            double dW_dX_dot  = dWdot_dX - ux*Wx;

            // Wdot and W at cell centers.
            double Wdot       = .5*(bdry_p->comp_Wadv() + bdry_pr->comp_Wadv());
            W                 = .5*(W + Wr);

            // Note that Ts, Ts_dot, Ts_inv_trans compute these matrices at cell-centers, so we need to modify y_val.
            mat3 Ts           = calc_Ts(L_val, Cx, y_val+.5*dy, Wx, W);
            mat3 Ts_dot       = calc_Ts_dot(Ld_val, dC_dX_dot, y_val+.5*dy, dW_dX_dot, vc, Wx, Wdot);
            mat3 Ts_inv_trans = calc_Ts_inv_trans(L_val, Cx, y_val+.5*dy, Wx, W);

            // Note that calc_l_trans uses x_val and y_val to compute the velocities, so x_val and y_val correspond to corners.
            mat3 l            = calc_l_trans(f, fr, fu, fp[gm+1], *bdry_p, *bdry_pr, L_val, Ld_val, x_val, y_val, Ts_inv_trans, i, j).transpose();
            mat3 sig          = calc_sig(Ts, Sigma);
            mat3 sig0         = calc_sig0(sig);
            double sbar       = calc_sbar(sig0);
            f.dev = sbar; // Save sbar as dev, because we've already computed it, and it requires derivatives and such.
            mat3 CD           = dt*(mat3((K - 2./3.*mu)*l.trace()) + mu*(l + l.transpose()));

            // Calculate ENO derivatives of transformed stress and effective temperature.
            f1 = (uc > 0)? st_eno2(hx, fp, -1)  : st_eno2(-hx, fp, 1);  // f1 holds the x derivatives.
            f2 = (vc > 0)? st_eno2(hy, fp, -gm) : st_eno2(-hy, fp, gm); // f2 holds the y derivatives.

            // Calculate the adaptive effective temperature update
            f.cchi  = f.chi;
            dtmp    = adaptive_plastic_term(sbar, f.cchi, f.ddev, dt); // 2*mu*dt*Dpl/sbar
            f.cchi -= uc*f1.chi + vc*f2.chi;

            // C:D^pl = 2*mu*sig0/sbar*Dpl, and dtmp = 2*mu*dt*Dpl/sbar (note dt included!)
            mat3 CDpl        = dtmp*sig0;

            // Compute terms that have transposes in the equation to save computation.
            mat3 lsig        = dt*l*sig;
            mat3 Tsdot_term  = dt*Ts_dot*Sigma*(Ts.transpose());

            // Combine all the matrices for one grand updates.
            mat3 all_updates = Ts_inv_trans.transpose()*(CD - CDpl + lsig + lsig.transpose() + dt*l.trace()*sig - Tsdot_term - Tsdot_term.transpose())*Ts_inv_trans;
            //mat3 all_updates = Ts_inv_trans.transpose()*(CD - CDpl - Tsdot_term - Tsdot_term.transpose())*Ts_inv_trans;

            // Calculate updates to stress by unpacking the update and including the advective terms.
            f.cs11 = all_updates.a11 - uc*f1.s11 - vc*f2.s11;
            f.cs12 = all_updates.a12 - uc*f1.s12 - vc*f2.s12;
            f.cs22 = all_updates.a22 - uc*f1.s22 - vc*f2.s22;
            f.cs33 = all_updates.a33 - uc*f1.s33 - vc*f2.s33;


            if(isnan(f.cu)   || fabs(f.cu)   > 1e10) { printf("u is weird! (%d %d) \n", i, j);}
            if(isnan(f.cv)   || fabs(f.cv)   > 1e10) { printf("v is weird! (%d %d) \n", i, j);}
            if(isnan(f.cs11) || fabs(f.cs11) > 1e10) { printf("s11 is weird! (%d %d) \n", i, j);}
            if(isnan(f.cs12) || fabs(f.cs12) > 1e10) { printf("s12 is weird! (%d %d) \n", i, j);}
            if(isnan(f.cs22) || fabs(f.cs22) > 1e10) { printf("s22 is weird! (%d %d) \n", i, j);}
            if(isnan(f.cs33) || fabs(f.cs33) > 1e10) { printf("s33 is weird! (%d %d) \n", i, j);}

            //printf("CD term: %d %d %d %g %g %g %g\n", curr_step, i, j, CD.a11, CD.a12, CD.a22, CD.a33);
            //printf("CDpl term: %d %d %d %g %g %g %g\n", curr_step, i, j, CDpl.a11, CDpl.a12, CDpl.a22, CDpl.a33);
            //printf("lsig term: %d %d %d %g %g %g %g\n", curr_step, i, j, lsig.a11, lsig.a12, lsig.a22, lsig.a33);
            //printf("Ts_inv_trans: %d %d %d %g %g %g %g\n", curr_step, i, j, Ts_inv_trans.a11, Ts_inv_trans.a12, Ts_inv_trans.a22, Ts_inv_trans.a33);
            //printf("Ts_dot: %d %d %d %g %g %g %g\n", curr_step, i, j, Tsdot_term.a11, Tsdot_term.a12, Tsdot_term.a22, Tsdot_term.a33);
            //printf("Stress updates: %d %d %d %g %g %g %g\n", curr_step, i, j, f.cs11, f.cs12, f.cs22, f.cs33);
            
            /*
             *if (CD.mod() > 1)         { printf("CD term: %d %d %d %g %g\n", curr_step, i, j, CD.mod(), l.mod()); }
             *if (CDpl.mod() > 1)       { printf("CDpl term: %d %d %d %g\n", curr_step, i, j, CDpl.mod()); }
             *if (lsig.mod() > 1)       { printf("lsig term: %d %d %d %g\n", curr_step, i, j, lsig.mod()); }
             *if (Tsdot_term.mod() > 1) { printf("Ts_dot: %d %d %d %g\n", curr_step, i, j, Tsdot_term.mod()); }
             */
        }
    }

    // Apply the updates to the fields
#pragma omp parallel for
    for(int j = 0; j < n; j++) {
        if (j > 0)    { for(int i = 0; i < m; i++) fm[i + gm*j].update(); }
        else          { for(int i = 0; i < m; i++) fm[i + gm*j].update_staggered(); }
    }

    for(int i = 0; i < m; i++) { bdrys[i].update(); }

    // Update the simulation time, and set the boundaries
    time += dt;
    set_boundaries();
    curr_step += 1;
}

/* Compute the updates for T and B. */
void shear_sim::update_bdry_profiles(double dt, double hx) {
    // Field pointers for looping over top and bottom boundaries.
    c_field    *fbp, *ftp;
    bdry_field *fp;

    // Upwinded x derivatives of current point on top and bottom boundaries.
    double Tx, T_advx, Bx, B_advx;

    // Velocity on the top and bottom, for advective derivatives.
    double uct, ucb;
    double dW_dx, dtp_dx, dbp_dx, ddtp_ddx, ddbp_ddx, ds22_dy, ds12_dy, sxx, ds11_dx, W;

    // Loop over the top and bottom boundaries simultaneously and compute the updates.
#pragma omp parallel for
    for(int i = 0; i < m; i++) { 
        // Set up the pointers.
        fbp = fm + i;
        ftp = fm + n*gm + i;
        fp  = bdrys + i;

        // Store the velocities.
        uct = ftp->u;
        ucb = fbp->u;

        // Calculate ENO derivatives of the boundary terms.
        (uct > 0)? top_bdry_eno2(Tx, T_advx, hx, fp, -1) : top_bdry_eno2(Tx, T_advx, -hx, fp, 1);
        (ucb > 0)? bot_bdry_eno2(Bx, B_advx, hx, fp, -1) : bot_bdry_eno2(Bx, B_advx, -hx, fp, 1);

        // Compute the update for the top and bottom profiles.
        fp->ctp = dt*fp->tp_adv - uct*Tx;
        fp->cbp = dt*fp->bp_adv - ucb*Bx;

        // Compute the needed derivatives for the top boundary.
        sxx         = .25*(ftp->s11 + ftp[-1].s11 + ftp[-gm].s11 + ftp[-gm-1].s11);
        W           = fp->comp_W();

        /* Should these be centered differences? */
        dW_dx       = xsp*(fp[1].comp_W() - W);
        dtp_dx      = xsp*(fp[1].tp - fp->tp); 
        ddtp_ddx    = xsp*xsp*(fp[1].tp - 2*fp->tp + fp[-1].tp);
        ds22_dy     = .5*ysp*(ftp->s22 - ftp[-gm].s22 + ftp[-1].s22 - ftp[-gm-1].s22);
        ds12_dy     = .5*ysp*(ftp->s12 - ftp[-gm].s12 + ftp[-1].s12 - ftp[-gm-1].s12);
        ds11_dx     = .5*xsp*(ftp->s11 - ftp[-1].s11  + ftp[-gm].s11 - ftp[-gm-1].s11);
        fp->ctp_adv = dt*(sxx*dW_dx*dtp_dx/W + sxx*ddtp_ddx + W*ds22_dy + dtp_dx*(ds12_dy + ds11_dx)) - uct*T_advx;

        // And needed derivatives for the bottom boundary.
        sxx         = .25*(fbp->s11 + fbp[-1].s11 + fbp[-gm].s11 + fbp[-gm-1].s11);
        dbp_dx      = xsp*(fp[1].bp - fp->bp); 
        ddbp_ddx    = xsp*xsp*(fp[1].bp - 2*fp->bp + fp[-1].bp);
        ds22_dy     = .5*ysp*(fbp->s22 - fbp[-gm].s22 + fbp[-1].s22  - fbp[-gm-1].s22);
        ds12_dy     = .5*ysp*(fbp->s12 - fbp[-gm].s12 + fbp[-1].s12  - fbp[-gm-1].s12);
        ds11_dx     = .5*xsp*(fbp->s11 - fbp[-1].s11  + fbp[-gm].s11 - fbp[-gm-1].s11);
        fp->cbp_adv = dt*(sxx*dbp_dx*dW_dx/W + sxx*ddbp_ddx + W*ds22_dy + dbp_dx*(ds12_dy + ds11_dx)) - ucb*B_advx;

        //printf("profile updates: %d %d %g %g %g %g\n", curr_step, i, fp->ctp, fp->cbp, fp->ctp_adv, fp->cbp_adv);
    }
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
void shear_sim::advection_step(double dt) { ; }

/** Carries out a projection step of the stresses using a multigrid solve,
 * to enforce the quasistaticity constraint.
 * \param[in] dt the timestep to use. */
void shear_sim::projection_step(double dt) { ; }

/** Calculates the net force at a cell corner using the stress tensors at the
 *  four neighboring cell centers.
 * \param[out] (fx, fy) the components of force.
 * \param[in]  (xf, yf) prefactors to apply to the x and y derivatives that go
 *                     into the computation of force.
 * \param[in]  fp a pointer to the grid point to consider. */
inline void shear_sim::net_force(double &fx, double &fy, double xf, double yf, c_field *fp) {
    c_field &f0 = fp[-gm-1], 
            &f1 = fp[-gm], 
            &f2 = fp[-1], 
            &f3 = *fp;

    // Calculate the net force in the x direction
    fx = (f3.s11 - f2.s11 + f1.s11 - f0.s11)*xf + (f3.s12 - f1.s12 + f2.s12 - f0.s12)*yf;

    // Calculate the net force in the y direction
    fy = (f3.s12 - f2.s12 + f1.s12 - f0.s12)*xf + (f3.s22 - f1.s22 + f2.s22 - f0.s22)*yf;
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
double shear_sim::adaptive_plastic_term(double sbar, double &chiv, double &ddev, double dt) {
    const double small_dev_cutoff = 1e-12;
    const double theta            = 0;
    const bool no_chidot          = false;
    bool adapt   = true;
    double dtl   = dt, 
           osbar = sbar, 
           dplas, 
           dchi1, 
           dchi2, 
           adt;

    // To avoid division-by-zero errors, check for the case when the
    // deviatoric stress is small, and just return with no changes
    if (sbar < small_dev_cutoff) return 0;
    do {
        // Calculate the plastic deformation
        dplas = 2*t_scale*mu*stz->Dplastic(sbar, chiv, dchi1, dchi2);

        // If the plastic deformation is large, then adaptively
        // choose the timestep so as not to make a large alteration in
        // the deviatoric stress magnitude
        if(fabs(dtl*dplas) > adapt_fac) {
            adt = adapt_fac/fabs(dplas);
            dtl -= adt;
            //printf("wtf. sbar: %g, dplas: %g, dt: %g, dtl: %g, dtl*dplas: %g, adt: %g, adt*dplas: %g\n", sbar, dplas, dt, dtl, dtl*dplas, adt, adt*dplas);
        } else {
            adapt = false;
            adt   = dtl;
        }
        sbar -= adt*dplas;
        //printf("sbar: %g\n", sbar);

        // Calculate the change in chi if needed
        if(!no_chidot) chiv += (dchi1*(chi_inf - chiv) + dchi2*(theta - chiv))*t_scale*adt;
    } while(adapt);


    // Calculate the quantity dt*2*mu*Dpl/sbar, which will be used in the
    // finite-difference update
    ddev = osbar - sbar;
    return ddev/osbar;
}

/** Sets the fields in the ghost regions according to the boundary conditions. */
void shear_sim::set_boundaries() {
    // Set the top and bottom boundaries
    const int g = n*gm;

    c_field *fp = fm;
    for(int i = 0; i <= m; i++, fp++) {
        fp->v   = 0;
        fp[g].v = 0;
        fp->u   = 2*fp[gm].u - fp[2*gm].u; // Interpolate U to the bottom boundary.
        fp[g].u = 2*fp[g-gm].u - fp[g-2*gm].u; // Interpolate U to the top boundary.
    }

    fp = fm;
    for(int j = 0; j <= n; j++, fp += gm) {
        fp->v   = 0;
        fp->u   = 0;
        fp[m].v = 0;
        fp[m].u = 0;
    }

    fp = fm - 2;
    // Handle the top and bottom boundary conditions in U (free) and V (zero), and also extrapolate over the boundary.
    for(int i = -2; i < m + 2; i++, fp++) {
        fp[-gm].extrap(*fp, fp[gm]);
        fp[-gm].s12   = -fp->s12; // s12, s22 = 0 on the bottom boundary.
        fp[-gm].s22   = -fp->s22;
        fp[-2*gm].extrap(fp[-gm], *fp); // Extrapolate the rest of them.

        
        fp[g].s12 = -fp[g-gm].s12; // s12, s22 = 0 on the top boundary.
        fp[g].s22 = -fp[g-gm].s22;
        fp[g].s11 = 2*fp[g-gm].s11 - fp[g-2*gm].s11; // Extrapolate the rest.
        fp[g].s33 = 2*fp[g-gm].s33 - fp[g-2*gm].s33;
        fp[g].chi = 2*fp[g-gm].chi - fp[g-2*gm].chi;
        fp[g+gm].extrap(fp[g], fp[g-gm]);
    }

    // Handle the right and left boundary conditions in U and V (both zero), and extrapolate over the boundaries.
    fp = fm - 2*gm;
    for (int j = -2; j < n + 2; j++, fp += gm) {
        // Extrapolate over the left boundary, and set the boundary velocities.
        fp[-1].extrap(*fp, fp[1]);
        fp[-2].extrap(fp[-1], *fp);

        // Extrapolate over the right boundary, and set the boundary velocities.
        fp[m].s11 = 2*fp[m-1].s11 - fp[m-2].s11;
        fp[m].s12 = 2*fp[m-1].s12 - fp[m-2].s12;
        fp[m].s22 = 2*fp[m-1].s22 - fp[m-2].s22;
        fp[m].s33 = 2*fp[m-1].s33 - fp[m-2].s33;
        fp[m].chi = 2*fp[m-1].chi - fp[m-2].chi;
        fp[m+1].extrap(fp[m], fp[m-1]);
    }

    // Extrapolate the boundary profiles over the left and right boundaries.
    bdry_field *bp = bdrys;
    bdrys[-1].extrap(*bp, bp[1]);
    bdrys[-2].extrap(bp[-1], *bp);
    bdrys[m].extrap(bp[m-1], bp[m-2]);
    bdrys[m+1].extrap(bp[m], bp[m-1]);
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

/** Calculates one-sided derivatives of cell-cornered fields using the second-order ENO2 scheme,
 * specifically for the advective derivative of the boundary terms, for use in computing the second
 * advective derivative terms. */
inline void shear_sim::bdry_adv_eno2(double &Wadv_d, double &Cadv_d, double hs, bdry_field *fp, int d) {
    Wadv_d = hs*eno2(fp[-d].comp_Wadv(), fp->comp_Wadv(), fp[d].comp_Wadv(), fp[2*d].comp_Wadv());
    Cadv_d = hs*eno2(fp[-d].comp_Cadv(), fp->comp_Cadv(), fp[d].comp_Cadv(), fp[2*d].comp_Cadv());
}

/* Calculates one-sided derivatives of cell-cornered fields using the second order ENO2 scheme,
 * specifically for the top boundary. */
inline void shear_sim::top_bdry_eno2(double &p_d, double &p_dd, double hs, bdry_field *fp, int d) {
    p_d  = hs*eno2(fp[-d].tp    , fp->tp    , fp[d].tp    , fp[2*d].tp);
    p_dd = hs*eno2(fp[-d].tp_adv, fp->tp_adv, fp[d].tp_adv, fp[2*d].tp_adv);
}

/* Calculates one-sided derivatives of cell-cornered fields using the second order ENO2 scheme,
 * specifically for the bottom boundary. */
inline void shear_sim::bot_bdry_eno2(double &p_d, double &p_dd, double hs, bdry_field *fp, int d) {
    p_d  = hs*eno2(fp[-d].bp    , fp->bp    , fp[d].bp    , fp[2*d].bp);
    p_dd = hs*eno2(fp[-d].bp_adv, fp->bp_adv, fp[d].bp_adv, fp[2*d].bp_adv);
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
inline st_field shear_sim::st_eno2(double hs, c_field *fp, int d) {
    c_field &f0 = fp[-d],
            &f1 = *fp,
            &f2 = fp[d],
            &f3 = fp[2*d];

    return st_field(hs*eno2(f0.s11, f1.s11, f2.s11, f3.s11),
                    hs*eno2(f0.s12, f1.s12, f2.s12, f3.s12),
                    hs*eno2(f0.s22, f1.s22, f2.s22, f3.s22),
                    hs*eno2(f0.s33, f1.s33, f2.s33, f3.s33),
                    hs*eno2(f0.chi, f1.chi, f2.chi, f3.chi));
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
    if(fflags&1)     output("u",   0,  k);
    if(fflags&2)     output("v",   1,  k);
    if(fflags&4)     output("s11", 2,  k);
    if(fflags&8)     output("s12", 3,  k);
    if(fflags&16)    output("s22", 4,  k);
    if(fflags&32)    output("s33", 5,  k);
    if(fflags&64)    output("chi", 6,  k);
    if(fflags&128)   output("tem", 7,  k);
    if(fflags&256)   output("dev", 8,  k);
    if(fflags&512)   output("X",   9,  k);
    if(fflags&1024)  output("Y",   10, k);
    if(fflags&2048)  { output("top", 11, k); output("bot", 12, k); }
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename. */
void shear_sim::output(const char *prefix, const int mode, const int sn) {
    // Determine whether to output a staggered field or not
    bool st = ((mode >= 2) && (mode <= 8));

    // Number of grid points in the x direction.
    int l  = st? m : m + 1;

    // Assemble the output filename and open the output file.
    char *bufc = reinterpret_cast<char*>(buf);
    sprintf(bufc, "%s/%s.%d", filename, prefix, sn);
    FILE *outf = safe_fopen(bufc, "wb");

    // Output the first line of the file.
    // See gnuplot.sourceforge.net/docs_4.2/node330.html
    int i, j;
    float *bp = buf + 1, 
          *be = buf + m + 1;

    // Set the first element of the first row to be the number of x points. The
    // rest of the first row is all of the x coordinates.
    *buf = l;
    for(i = 0; i < l; i++) *(bp++) = ax + (st? (i + 0.5)*dx : i*dx);
    fwrite(buf, sizeof(float), l+1, outf);

    // Output the field values to the file.
    if (mode < 11) {
        c_field *fp = fm;
        int max_y   = st? n : n + 1;
        for (j = 0; j < max_y; j++) {
            int i = 0;

            // First element in every row for binary matrix format is the y coordinate
            *buf = ay + (st? (j + 0.5)*dy : j*dy); bp = buf + 1;
            switch(mode) {
                case  0: while(bp < be) { *(bp++) = (fp++)->u; i++; }      break;
                case  1: while(bp < be) { *(bp++) = (fp++)->v; i++; }      break;
                case  2: while(bp < be) { *(bp++) = (fp++)->s11; i++; }    break;
                case  3: while(bp < be) { *(bp++) = (fp++)->s12; i++; }    break;
                case  4: while(bp < be) { *(bp++) = (fp++)->s22; i++; }    break;
                case  5: while(bp < be) { *(bp++) = (fp++)->s33; i++; }    break;
                case  6: while(bp < be) { *(bp++) = (fp++)->chi; i++; }    break;
                case  7: while(bp < be) { *(bp++) = (fp++)->chi*TZ; i++; } break;
                case  8: while(bp < be) { *(bp++) = (fp++)->dev; i++; }    break;
                case  9: while(bp < be) { *(bp++) = (fp++)->X; i++; }      break;
                case 10: while(bp < be) { *(bp++) = (fp++)->Y; i++; }      break;
                case 11: while(bp < be) { *(bp++) = (fp++)->cu; i++; }     break;
                case 12: while(bp < be) { *(bp++) = (fp++)->cv; i++; }     break;
                case 13: while(bp < be) { *(bp++) = (fp++)->cs11; i++; }   break;
                case 14: while(bp < be) { *(bp++) = (fp++)->cs12; i++; }   break;
                case 15: while(bp < be) { *(bp++) = (fp++)->cs22; i++; }   break;
                case 16: while(bp < be) { *(bp++) = (fp++)->cs33; i++; }
            }

            // Handle the last point.
            // Why do we have to do it this way?!?!!?
            if (!st) { 
                *bp = (mode == 9)? buf[1] + bx - ax : (*fp).fval(mode); 
                if (mode == 0 && (*fp).fval(mode) != 0) printf("u val: %g %d %d\n", (*fp).fval(mode), i-1, j);
                if (mode == 1 && (*fp).fval(mode) != 0) printf("v val: %g %d %d\n", (*fp).fval(mode), i-1, j);
            }

            // Pass over the ghost points
            fp += 4;

            // Write it out
            fwrite(buf, sizeof(float), l+1, outf);
        }
    }
    // Output the boundary profile values to the file.
    else {
        bdry_field *bdp = bdrys;
        if (mode == 11) { *buf = by; bp = buf+1; while (bp < be+1) *(bp++) = bdp->tp; }
        else            { *buf = ay; bp = buf+1; while (bp < be+1) *(bp++) = bdp->bp; }

        // Write it out
        fwrite(buf, sizeof(float), l+1, outf);
    }

    // Close the file
    fclose(outf);
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
    int cut=int(l*gaussian_fac),nx=n+(cut<<1);
    double llinv=0.5/(l*l),chi0=tem_avg/TZ,chi1=tem_stdev/(gaussian_normal_factor(llinv,cut)*TZ);

    // Create grid of Gaussian noise
    double *no=new double[m*nx],*np=no,*ne=no+m*nx;
    for(;np<ne-1;np+=2) box_muller(*np,np[1]);
    if(np==ne-1) {
        double r1;
        box_muller(*np,r1);
    }

    // Initialize the simulation fields
    np=no+m*cut;
#pragma omp parallel for
    for(int j=0;j<=n;j++) {
        double Y=ay+j*dy;
        c_field *fp=fm+m*j;
        for(int i=0;i<m;i++,fp++) {
            c_field &f=*fp;
            double gau=0;

            // Set stresses and velocities to zero
            f.u=f.v=0;

            // Set reference map field to match initial coordinate
            f.X=ax+i*dx;f.Y=Y;

            // Calculate chi as a gaussian smoothing of the random field
            if(j<n) {
                for(int cj=j-cut;cj<=j+cut;cj++) for(int ci=i-cut;ci<=i+cut;ci++)
                    gau+=np[cj*m+step_mod(ci,m)]*exp(-llinv*((cj-j)*(cj-j)+(ci-i)*(ci-i)));
                f.chi=chi0+gau*chi1;
                f.s11=f.s12=f.s22=f.s33=0;
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

/* Compute the spatial component of the transformation gradient. */
mat3 shear_sim::calc_Ts(double L, double dC_dX, double Y, double dW_dX, double W) {
    return mat3(L             , 0, 0,
                dC_dX + Y*dW_dX, W, 0,
                0             , 0, 1);
}

/* Compute the inverse transpose of the spatial component of the transformation gradient. */
mat3 shear_sim::calc_Ts_inv_trans(double L, double dC_dX, double Y, double dW_dX, double W) {
    double Linv = 1./L;
    double Winv = 1./W;
    return mat3(Linv, -Linv*Winv*(dC_dX + Y*dW_dX), 0,
                   0,                         Winv, 0,
                   0,                            0, 1);
}

/* Compute the advective derivative of the spatial component of the transformation gradient. */
mat3 shear_sim::calc_Ts_dot(double Ldot, double dC_dX_dot, double Y, double dW_dX_dot, double VY, double dW_dX, double Wdot) {
    return mat3(                              Ldot,    0, 0,
                dC_dX_dot + Y*dW_dX_dot + VY*dW_dX, Wdot, 0,
                                                 0,    0, 0);
}

/* Computes the transpose of l. Computes the transpose because we use the convention l_{ij} = \partial_j u_i.
 * However, the nice transformation formula T_S^{-T}\nabla_X u is only valid for l_{ij} = \partial_i u_j. */
mat3 shear_sim::calc_l_trans(c_field &f, c_field &fr, c_field &fu, c_field &fur, 
                  bdry_field &bf, bdry_field &bfr, double L, double Ldot, double X, double Y,
                  mat3 &Ts_inv_trans, int i, int j){

   // Compute the boundary profile terms at the current and adjacent grid point.
   // Note we are using these to compute the velocities, so they are located at cell corners.
   double W     = bf.comp_W();
   double Wdot  = bf.comp_Wadv();
   double Wr    = bfr.comp_W();
   double Wdotr = bfr.comp_Wadv();
   double Cdot  = bf.comp_Cadv();
   double Cdotr = bfr.comp_Cadv();
   
   // No periodicity ever, so no wrapping is necessary.
   double Xr = X + dx;
   double Yu = Y + dy;

   // Compute the physical velocities at the neighboring grid points.
   double uphys    = Ldot*X + L*f.u;
   double vphys    = Wdot*Y + W*f.v + Cdot;

   double uphys_r  = Ldot*Xr + L*fr.u;
   double vphys_r  = Wdotr*Y + Wr*fr.v + Cdotr;

   double uphys_u  = Ldot*X + L*fu.u;
   double vphys_u  = Wdot*Yu + W*fu.v + Cdot;

   double uphys_ur = Ldot*Xr + L*fur.u;
   double vphys_ur = Wdotr*Yu + Wr*fur.v + Cdotr;

   // Compute the physical velocity derivatives with respect to the transformed coordinates at cell centers.
   double du_dX = .5*xsp*(uphys_r - uphys + uphys_ur - uphys_u);
   double du_dY = .5*ysp*(uphys_u - uphys + uphys_ur - uphys_r);
   double dv_dX = .5*xsp*(vphys_r - vphys + vphys_ur - vphys_u);
   double dv_dY = .5*ysp*(vphys_u - vphys + vphys_ur - vphys_r);

   //printf("Velocities: %d %d %g %g %g %g\n", i, j, vphys, vphys_r, vphys_u, vphys_ur);
   //printf("Velocity derivatives: %d %d %g %g %g %g\n\n", i, j, du_dX, du_dY, dv_dX, dv_dY);

   // Transform the transformed derivatives to physical derivatives.
   return Ts_inv_trans*mat3(du_dX, dv_dX, 0,
                            du_dY, dv_dY, 0,
                                0,     0, 0);
}

/* Calculate sbar, given as the Frobenius norm of sig0 with a factor of .5. */
double shear_sim::calc_sbar(mat3 &sig0) {
    return sqrt(.5*(sig0.a11*sig0.a11 + 2*sig0.a12*sig0.a12 + sig0.a22*sig0.a22 + sig0.a33*sig0.a33));
}
