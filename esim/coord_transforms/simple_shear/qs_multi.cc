#include "qs_multi.hh"
#include "shear_sim.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent shear_sim
 * class.
 * \param[in] ss a reference to the parent shear_sim class. */
qs_multi::qs_multi(shear_sim &ss) : m(ss.m), n(ss.lsn), mn(m*n),
    y_prd(ss.y_prd), K(ss.K), mu(ss.mu), L(K-2./3.*mu), dx(ss.dx), dy(ss.dy),
    xsp(ss.xsp), ysp(ss.ysp), xsp2(xsp*xsp), ysp2(ysp*ysp), xyfac(xsp2 + ysp2),
    acc(tgmg_accuracy(2.,1e4)), regz(mat(0)), rego(mat(1)), b(new vec[mn]),
    z(new vec[mn]), mg(*this,b,z), tp() {}

/** The class destructor frees the dynamically allocated memory. */
qs_multi::~qs_multi() {
    delete [] z;
    delete [] b;
}

/** Sets constants that are used during the multigrid solve.
 * \param[in] viscdt the current value of the viscosity divided by the
 *                   timestep.
 * \param[in] lt the constant lambda*time, which features in many terms in the
 *               linear system. */
void qs_multi::init(double viscdt, double lt) {
    double ltsq = lt*lt, L = K - 2./3.*mu;

    // Set 2x2 matrices that are used within the multigrid solve
    regh   = mat((L + 2*mu + ltsq*mu)*xsp2, (lt*mu + lt*ltsq*mu)*xsp2,
                         -lt*(L + mu)*xsp2,        mu*(1 + ltsq)*xsp2);
    regv   = mat(mu*ysp2,      lt*mu*ysp2,
                       0, (L + 2*mu)*ysp2);
    regc   = -2.*(regh+regv);
    regci  = regc.inverse();
    regd  = mat(-0.5*xsp*ysp*lt*mu, 0.25*xsp*ysp*(L + mu - 2*ltsq*mu),
                0.25*xsp*ysp*(L + mu), -0.25*xsp*ysp*lt*(L + 3*mu));

    // Multiply through the matrices by the inverse central element to
    // precondition the system
    regh=regci*regh;
    regv=regci*regv;
    regd=regci*regd;

    // Set up the multigrid hierarchy
    mg.setup();
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<qs_multi,vec,mat>;
template class tgmg_base<tgmg_level<vec,mat>,vec,mat>;
