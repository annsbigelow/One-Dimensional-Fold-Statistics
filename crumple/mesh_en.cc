#include "mesh.hh"

#include <cstdio>
#include <cmath>
#include <cstring>

/** Global function used to link the GSL optimizer back to the class routine
 * for evaluating the energy. */
double mesh_f(const gsl_vector *v,void *params) {
    return reinterpret_cast<mesh*>(params)->fun(v->data);
}

/** Global function used to link the GSL optimizer back to the class routine
 * for evaluating the derivative of the energy. */
void mesh_df(const gsl_vector *v,void *params,gsl_vector *df) {
    reinterpret_cast<mesh*>(params)->dfun(v->data,df->data);
}

/** Global function used to link the GSL optimizer back to the class routine
 * for evaluating the energy and its derivative. */
void mesh_fdf(const gsl_vector *v,void *params,double *f,gsl_vector *df) {
    mesh *mp=reinterpret_cast<mesh*>(params);
    mp->dfun(v->data,df->data);
    *f=mp->fun(v->data);
}

/** Checks that the finite-difference derivative of the energy gives the
 * acceleration. */
void mesh::check_deriv(double t_) {
    fzt=t_;
    int i;
    const double h=1e-5,hfac=0.5/h;
    double e,*fder=new double[6*n],*fp,l2=0.;

    // Perform centered finite differences of the energy function
    for(i=0;i<3*n;i++) {
        pts[i]+=h;e=fun(pts);
        pts[i]-=2*h;fder[i]=hfac*(e-fun(pts));
        pts[i]+=h;
    }

    // Compute the acceleration and print out a comparison of the two
    dfun(pts,fder+3*n);
    puts("# Index, Position, Fdiff energy, -Acceleration");
    for(fp=fder,i=0;i<n;i++,fp+=3) {
        e=*fp-fp[3*n];l2+=e*e;
        e=fp[1]+fp[3*n+1];l2+=e*e;
        e=fp[2]+fp[3*n+2];l2+=e*e;
        printf("%d  %g %g %g  %g %g %g  %g %g %g\n",i,pts[3*i],pts[3*i+1],
               pts[3*i+2],*fp,fp[1],fp[2],fp[3*n],fp[3*n+1],fp[3*n+2]);
    }
    printf("# RMS error: %g\n",sqrt(l2/static_cast<double>(3*n)));

    delete [] fder;
}

void mesh::minimize_energy(double t_) {
    fzt=t_;
    int iter=0,status,piter=0;
    double t0=0;
    int pdof=3*n;

    // Inititalize the function to be minimized
    gsl_multimin_function_fdf my_func;
    my_func.n=pdof;
    my_func.f=mesh_f;
    my_func.df=mesh_df;
    my_func.fdf=mesh_fdf;
    my_func.params=this;

    // Set the starting point for the minimization
    gsl_vector *x=gsl_vector_alloc(3*n);
    memcpy(x->data,pts,pdof*sizeof(double));

    // Initialize the BFGS minimizer
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    T=gsl_multimin_fdfminimizer_conjugate_fr;
    s=gsl_multimin_fdfminimizer_alloc(T,pdof);
    gsl_multimin_fdfminimizer_set(s,&my_func,x,0.01,1e-3);

    // Do the BFGS iterations
    piter=0;
    do {
        if(iter%100==0) min_message(piter,iter,s,t0);
        iter++;
        if(gsl_multimin_fdfminimizer_iterate(s)) break;
        status=gsl_multimin_test_gradient(s->gradient,1e-10);
    } while(status==GSL_CONTINUE&&iter<10000);

    // Print final information
    min_message(piter,iter,s,t0);
    puts(status==GSL_SUCCESS?"Minimum found":"Minimization failed");

    // Copy the vector contents into the coefficient array within the class
    memcpy(pts,s->x->data,pdof*sizeof(double));

    // Free the dynamically allocated memory
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
}

/** Prints a message about the status of the energy minimization.*/
void mesh::min_message(int &piter,int iter,gsl_multimin_fdfminimizer *s,double &t0) {
    if(piter==iter) {
        printf("Iteration %d, residual %.10f\n",iter,s->f);
        t0=wtime();
    } else {
        double t1=wtime();
        printf("Iteration %d, residual %.10f [%.6g s/iter]\n",iter,s->f,(t1-t0)/(iter-piter));
        t0=t1;
        piter=iter;
    }
}
