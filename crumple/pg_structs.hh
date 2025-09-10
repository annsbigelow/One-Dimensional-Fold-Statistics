#ifndef PG_STRUCTS_HH
#define PG_STRUCTS_HH

struct point_info {
    /** The numerical ID of this point. */
    int id;
    /** The x position of this point. */
    double x;
    /** The y position of this point. */
    double y;
    /** The z position of this point. */
    double z;
    inline void set(int id_,double *p) {
        id=id_;x=*p;y=p[1];z=p[2];
    }
};

struct oflow_info {
    /** The block index that this point belongs in. */
    int ijk;
    /** The memory slot within this block that the point belongs in. */
    int s;
    /** The numerical ID of this point. */
    int id;
    /** The x position of this point. */
    double x;
    /** The y position of this point. */
    double y;
    /** The z position of this point. */
    double z;
    inline void set(double ijk_,int s_,int id_,double *p) {
        ijk=ijk_;s=s_;id=id_;
        x=*p;y=p[1];z=p[2];
    }
    inline void put(point_info *p) {
        p->id=id;p->x=x;p->y=y;p->z=z;
    }
};

#endif
