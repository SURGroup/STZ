#include "qs_multi.hh"
#include "shear_sim.hh"
#include "dmg_sim.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent shear_sim
 * class.
 * \param[in] ss a reference to the parent shear_sim class. */
qs_multi::qs_multi(shear_sim &ss) : m(ss.m), n(ss.n+1), mn(m*n),
	x_prd(true), y_prd(false), z(ss.vel), K(ss.K+ss.mu/3.0),
	qK(.25*K*ss.xsp*ss.ysp), Kmu(K+ss.mu), mu(ss.mu), dx(ss.dx), dy(ss.dy),
	xsp(ss.xsp), ysp(ss.ysp), xsp2(ss.xsp*ss.xsp), ysp2(ss.ysp*ss.ysp),
	xyfac(xsp2+ysp2), tmufac1(Kmu*xsp2+mu*ysp2), tmufac2(mu*xsp2+Kmu*ysp2),
	bcs(ss.bcs), viscdt_prev(-1), regz(mat(0)), acc(tgmg_accuracy(2*(tmufac1+tmufac2),1e4)),
	mg(*this,ss.src,z), tp() {}

/** The constructor sets many internal constants from the parent dmg_sim
 * class.
 * \param[in] ss a reference to the parent dmg_sim class. */
qs_multi::qs_multi(dmg_sim &ss) : m(ss.m), n(ss.n+1), mn(m*n),
	x_prd(true), y_prd(false), z(ss.vel), K(ss.K+ss.mu/3.0),
	qK(.25*K*ss.xsp*ss.ysp), Kmu(K+ss.mu), mu(ss.mu), dx(ss.dx), dy(ss.dy),
	xsp(ss.xsp), ysp(ss.ysp), xsp2(ss.xsp*ss.xsp), ysp2(ss.ysp*ss.ysp),
	xyfac(xsp2+ysp2), tmufac1(Kmu*xsp2+mu*ysp2), tmufac2(mu*xsp2+Kmu*ysp2),
	bcs(ss.bcs), viscdt_prev(-1), regz(mat(0)), acc(tgmg_accuracy(2*(tmufac1+tmufac2),1e4)),
	mg(*this,ss.src,z), tp() {}

/** Sets constants that are used during the multigrid solve.
 * \param[in] viscdt the current value of the viscosity divided by the
 *		     timestep. */
void qs_multi::init(double viscdt) {

	if(viscdt!=viscdt_prev) {
		viscdt_prev=viscdt;

		// Set elasticity-related constants
		Kmx=(Kmu+viscdt)*xsp2;
		Kmy=(Kmu+viscdt)*ysp2;
		mx=(mu+viscdt)*xsp2;
		my=(mu+viscdt)*ysp2;
		cxx=-2*(tmufac1+viscdt*xyfac);
		cyy=-2*(tmufac2+viscdt*xyfac);

		// Set 2x2 matrices that are used within the multigrid solve
		regh=mat(Kmx,0,0,my);
		regv=mat(mx,0,0,Kmy);
		regc=mat(cxx,0,0,cyy);
		regci=mat(1/cxx,0,0,1/cyy);
		regcb=mat(bcs*cxx,0,0,bcs*cyy);
		regcbi=mat(1/(bcs*cxx),0,0,1/(bcs*cyy));

		mg.setup();
	}
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<qs_multi,vec,mat>;
template class tgmg_base<tgmg_level<vec,mat>,vec,mat>;
