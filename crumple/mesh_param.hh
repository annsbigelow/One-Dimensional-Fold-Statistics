#ifndef MESH_PARAM_HH
#define MESH_PARAM_HH

struct mesh_param {
    /** The spring constant. */
    const double K;
    /** The drag constant. */
    const double drag;
    /** The dashpot damping constant. */
    const double B;
    /** The edge bending force. */
    const double kappa;
    /** The effective diameter of repulsive potential. */
    double sigma;
    /** Whether to include dashpots. */
    const bool dashpot;
    /** Whether to use the bendable sheet model, which also has connections
     * between triangle edges. */
    const bool bsheet_model;
    /** Whether to add pairwise repulsion for contact avoidance. */
    const bool repulsion;
    /** Whether to fix the mesh boundaries. */
    const bool fix_boundary;
    mesh_param(double K_,double drag_,bool fix_boundary_) : K(K_), drag(drag_),
        B(0.),  kappa(0.), sigma(0.), dashpot(false),
	bsheet_model(false), repulsion(false), fix_boundary(fix_boundary_) {}
    mesh_param(double K_,double drag_,double B_,bool fix_boundary_) : K(K_), drag(drag_),
	B(B_), kappa(0.), sigma(0.), dashpot(true),
	bsheet_model(false), repulsion(false), fix_boundary(fix_boundary_) {}
    mesh_param(double K_,double drag_,double B_,double kappa_,bool fix_boundary_) : K(K_), drag(drag_),
	B(B_), kappa(kappa_), sigma(0.), dashpot(true),
	bsheet_model(true), repulsion(false), fix_boundary(fix_boundary_) {}
    mesh_param(double K_,double drag_,double kappa_,bool dashpot_,bool fix_boundary_) : K(K_), drag(drag_),
	B(0.), kappa(kappa_), sigma(0.), dashpot(dashpot_),
	bsheet_model(false), repulsion(false), fix_boundary(fix_boundary_) {}
    mesh_param(double K_,double drag_,double B_,double kappa_) : K(K_), drag(drag_),
        B(B_), kappa(kappa_), sigma(sqrt(3)), dashpot(true),
        bsheet_model(true), repulsion(true), fix_boundary(false) {}
};

#endif
