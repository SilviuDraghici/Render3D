#ifndef MESHES_H
#define MESHES_H

#include "utils.h"
#include "objects.h"

extern double num_intersection_tests;
extern double num_intersect_calls;

class BVH_Node {
    public:
    virtual BVH_Node *intersect(struct ray *r, double *lambda,
                           point *bary_coords) = 0;
    virtual bool isFace() = 0;
    virtual double min_x() = 0;
    virtual double min_y() = 0;
    virtual double min_z() = 0;
    virtual double max_x() = 0;
    virtual double max_y() = 0;
    virtual double max_z() = 0;
};

class TriangleFace : public BVH_Node {
   protected:
    point p1, p2, p3;

   public:
    friend std::ostream &operator<<(std::ostream &strm, const TriangleFace &a);
    TriangleFace();
    TriangleFace(point p1, point p2, point p3);
    TriangleFace(double x1, double y1, double z1, double x2, double y2,
                 double z2, double x3, double y3, double z3);
    BVH_Node *intersect(struct ray *r, double *lambda, point *bary_coords);
    bool isFace();
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

class BoundingBox : public BVH_Node {
    int depth;

   public:
    BVH_Node *c1, *c2;
    double b_min_x, b_max_x;
    double b_min_y, b_max_y;
    double b_min_z, b_max_z;
    void setChildren(TriangleFace *faces[], int start, int end);
    void setBounds(double min_x, double max_x, double min_y, double max_y,
                   double min_z, double max_z);
    double min_x();
    double min_y();
    double min_z();
    double max_x();
    double max_y();
    double max_z();
    BVH_Node *intersect(struct ray *r, double *lambda, point *bary_coords);
    bool isFace();
};

class BoundingBox_Visible : public BoundingBox {
    double width = 0.02;

   public:
    color col;
    BoundingBox_Visible();
    color getCol();
    BVH_Node *intersect(struct ray *r, double *lambda, point *bary_coords);
};

class Mesh : public Object {
    // used for reading in Mesh
    int num_vertices = 0;
    point *vertices = NULL;

    int num_normals;
    point *normals = NULL;

    int num_faces = 0;
    TriangleFace **faces = NULL;

    point bary_coords;

    BoundingBox *box = NULL;

   public:
    using Object::Object;
    void setMesh(const char *filename);
    void intersect(struct ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);
};

class PQ_Node {
   public:
    double lambda;
    BVH_Node *node;
    PQ_Node(double l, BVH_Node *n) {
        lambda = l;
        node = n;
    }
    PQ_Node(double l) {
        lambda = l;
        node = NULL;
    }

    friend bool operator<(const PQ_Node &lhs, const PQ_Node &rhs) {
        return lhs.lambda < rhs.lambda;
    }
    friend bool operator>(const PQ_Node &lhs, const PQ_Node &rhs) {
        return lhs.lambda > rhs.lambda;
    }

    friend std::ostream &operator<<(std::ostream &strm, const PQ_Node &a) {
        return strm << a.lambda;
    }
};

int comp_face_max_x(const void *a, const void *b);
int comp_face_max_y(const void *a, const void *b);
int comp_face_max_z(const void *a, const void *b);

#endif