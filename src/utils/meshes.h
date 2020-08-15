#ifndef MESHES_H
#define MESHES_H

#include "utils.h"
#include "objects.h"
#include "BVH.h"

class TriangleFace : public Primitive {
   protected:
    point p1, p2, p3;

   public:
    friend std::ostream &operator<<(std::ostream &strm, const TriangleFace &a);
    TriangleFace();
    TriangleFace(point p1, point p2, point p3);
    TriangleFace(double x1, double y1, double z1, double x2, double y2,
                 double z2, double x3, double y3, double z3);
    double intersect(struct ray *r, double lambda);
    void intersect(struct ray *r, double *lambda, point *bary_coords);
    bool isprim();
    double min_x() const;
    double min_y() const;
    double min_z() const;
    double max_x() const;
    double max_y() const;
    double max_z() const;
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

class Mesh : public Object {
    // used for reading in Mesh
    int num_vertices = 0;
    point *vertices = NULL;

    int num_normals;
    point *normals = NULL;

    int num_faces = 0;
    PrimitiveData *faces = NULL;

    point bary_coords;

    BVH bvh;

   public:
    using Object::Object;
    void setMesh(const char *filename);
    void set_canonical_bounds();
    void intersect(struct ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);
};

#endif