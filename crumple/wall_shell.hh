#ifndef WALL_SHELL_HH
#define WALL_SHELL_HH

#include "voro++.hh"
using namespace voro;

#include <cmath>

struct wall_shell : public wall {
    public:
        wall_shell(double xc_,double yc_,double zc_,double rc,double sc,int w_id_=-99)
            : w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), lc(rc-sc), uc(rc+sc) {}
        bool point_inside(double x,double y,double z) {
            double rsq=(x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc);
            return rsq>lc*lc&&rsq<uc*uc;
        }
        template<class v_cell>
        bool cut_cell_base(v_cell &c,double x,double y,double z) {
            double xd=x-xc,yd=y-yc,zd=z-zc,dq=xd*xd+yd*yd+zd*zd,dq2;
            if (dq>1e-5) {
                dq2=2*(sqrt(dq)*lc-dq);
                dq=2*(sqrt(dq)*uc-dq);
                return c.nplane(xd,yd,zd,dq,w_id)&&c.nplane(-xd,-yd,-zd,-dq2,w_id);
            }
            return true;
        }
        bool cut_cell(voronoicell &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
        bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
    private:
        const int w_id;
        const double xc,yc,zc,lc,uc;
};

#endif
