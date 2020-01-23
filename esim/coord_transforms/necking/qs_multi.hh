#ifndef QS_MULTI_HH
#define QS_MULTI_HH

#include <cstdio>
#include <cstdlib>

#include "vec.hh"
#include "mat.hh"
#include "tgmg.hh"

class shear_sim;

struct qs_multi {
    const int m, n, mn;
    static const bool x_prd=true;
    const bool y_prd;
    static const char gs_mode=0;
    const double K, mu, L;
    const double dx, dy, xsp, ysp, xsp2, ysp2, xyfac;
    const double acc;
    const mat regz;
    const mat rego;
    mat regh,regv,regd,regc,regci;
    vec* const b;
    vec* const z;
    qs_multi(shear_sim &ss);
    ~qs_multi();
    void init(double viscdt, double lt);
    inline bool bdry(int i, int ij) {
        return !y_prd && (top_bdry(ij) || bot_bdry(ij));
    }
    inline bool top_bdry(int ij) { return ij >= mn-m; }
    inline bool bot_bdry(int ij) { return ij < m;     }
    inline mat a_dl(int i, int ij)     { return bdry(i,ij) ? regz   : regd;  }
    inline mat a_dr(int i, int ij)     { return bdry(i,ij) ? regz   : -regd; }
    inline mat a_ul(int i, int ij)     { return bdry(i,ij) ? regz   : -regd; }
    inline mat a_ur(int i, int ij)     { return bdry(i,ij) ? regz   : regd;  }
    inline mat a_cl(int i, int ij)     { return bdry(i,ij) ? regz   : regh;  }
    inline mat a_cr(int i, int ij)     { return bdry(i,ij) ? regz   : regh;  }
    inline mat a_dc(int i, int ij)     { return bdry(i,ij) ? regz   : regv;  }
    inline mat a_uc(int i, int ij)     { return bdry(i,ij) ? regz   : regv;  }
    inline mat a_cc(int i, int ij)     { return rego;  }
    inline vec inv_cc(int i, int ij, vec v) { return v; }
    inline vec mul_a(int i, int ij) {
        if(y_prd) {

            // Set pointers to handle the x periodicity
            vec *zp = z + ij, *zd = zp + (i==0? m-1 : -1),
            *zu = zp + (i == m-1? 1-m : 1);

            // Set displacements to handle the y periodicity, and compute the
            // result
            int us = top_bdry(ij)? -m*(n-1) :  m;
            int ds = bot_bdry(ij)?  m*(n-1) : -m;
            return regd*(zd[ds]-zu[ds]-zd[us]+zu[us])
                  +regv*(zp[ds]+zp[us])+regh*(*zd+*zu);
        } else {

            // Deal with the Dirichlet boundary conditions for the top and
            // bottom boundaries
            if(bot_bdry(ij)||top_bdry(ij)) return vec(0.);

            // Set pointers to handle the x periodicity, and compute the result
            vec *zp = z + ij, *zd = zp + (i==0? m-1 : -1),
            *zu = zp + (i == m-1? 1-m : 1);
            return regd*(zd[-m]-zu[-m]-zd[m]+zu[m])
                  +regv*(zp[-m]+zp[m])+regh*(*zd+*zu);
        }
    }
    inline void solve() {
        if(!mg.solve_v_cycle(tp)) exit(1);
    }
    /** The multigrid solver. */
    tgmg<qs_multi,vec,mat> mg;
    tgmg_predict tp;
};

#endif
