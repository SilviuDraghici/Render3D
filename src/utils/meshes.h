#ifndef MESHES_H
#define MESHES_H

#include "utils.h"
#include "objects.h"

extern double num_intersection_tests;
extern double num_intersect_calls;

class Intersectable {
    public:
    virtual Intersectable *intersect(struct ray *r, double *lambda,
                           point *bary_coords) = 0;
};

class TriangleFace : public Intersectable {
   protected:
    point p1, p2, p3;

   public:
    friend std::ostream &operator<<(std::ostream &strm, const TriangleFace &a);
    TriangleFace();
    TriangleFace(point p1, point p2, point p3);
    TriangleFace(double x1, double y1, double z1, double x2, double y2,
                 double z2, double x3, double y3, double z3);
    Intersectable *intersect(struct ray *r, double *lambda, point *bary_coords);
    double min_x();
    double min_y();
    double min_z();
    double max_x();
    double max_y();
    double max_z();
    virtual point normal(point *bary_coords);
};

class TriangleFace_N : public TriangleFace {
    point n1, n2, n3;

   public:
    friend std::ostream &operator<<(std::ostream &strm,
                                    const TriangleFace_N &a);
    using TriangleFace::TriangleFace;
    void setNormals(point n1, point n2, point n3);
    point normal(point *bary_coords);
};

class BoundingBox : public Intersectable {
    int depth;

   protected:
    Intersectable *c1, *c2;
    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;

   public:
    void setChildren(TriangleFace_N *faces, int start, int end);
    void setBounds(double min_x, double max_x, double min_y, double max_y,
                   double min_z, double max_z);
    Intersectable *intersect(struct ray *r, double *lambda, point *bary_coords);
};

class VisibleBoundingBox : public BoundingBox {
    double width = 0.02;

   public:
    color col;
    VisibleBoundingBox();
    color getCol();
    Intersectable *intersect(struct ray *r, double *lambda, point *bary_coords);
};

class Mesh : public Object {
    // used for reading in Mesh
    int num_vertices = 0;
    point *vertices = NULL;

    int num_normals;
    point *normals = NULL;

    int num_faces = 0;
    TriangleFace_N *faces = NULL;

    point bary_coords;

    BoundingBox *box = NULL;

   public:
    using Object::Object;
    void setMesh(const char *filename, point *cam_gaze);
    void intersect(struct ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);
};

int comp_face_max_x(const void *a, const void *b);
int comp_face_max_y(const void *a, const void *b);
int comp_face_max_z(const void *a, const void *b);

#endif