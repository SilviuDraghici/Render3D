#ifndef MESHES_H
#define MESHES_H

#include <list>

#include "utils.h"
#include "objects.h"
#include "BVH.h"

void adjust_indexes(int& i1, int& i2, int& i3, int array_size);

class TriangleFace : public Primitive {
   protected:
    point p1, p2, p3;

   public:
    friend std::ostream &operator<<(std::ostream &strm, const TriangleFace &a);
    TriangleFace();
    TriangleFace(point p1, point p2, point p3);
    TriangleFace(double x1, double y1, double z1, double x2, double y2,
                 double z2, double x3, double y3, double z3);
    double intersect(struct Ray *r, double lambda);
    void intersect(struct Ray *r, double *lambda, point *bary_coords);
    bool isprim();
    double min_x() const;
    double min_y() const;
    double min_z() const;
    double max_x() const;
    double max_y() const;
    double max_z() const;
    virtual point normal(point *bary_coords);
    virtual void texture_coordinates(double *a, double *b, point *bary_coords);
};

class TriangleFace_N : virtual public TriangleFace {
    point n1, n2, n3;

   public:
    friend std::ostream &operator<<(std::ostream &strm,
                                    const TriangleFace_N &a);
    using TriangleFace::TriangleFace;
    void setNormals(point n1, point n2, point n3);
    point normal(point *bary_coords);
};

class TriangleFace_T : virtual public TriangleFace {
    double p1_u, p1_v, p2_u, p2_v, p3_u, p3_v;

   public:
    friend std::ostream &operator<<(std::ostream &strm,
                                    const TriangleFace_T &a);
    using TriangleFace::TriangleFace;
    void set_texture_coordinates(double p1_u, double p1_v, double p2_u, double p2_v, double p3_u, double p3_v);
    void texture_coordinates(double *a, double *b, point *bary_coords);
};

class TriangleFace_T_N : public TriangleFace_T, public TriangleFace_N {
   public:
    friend std::ostream &operator<<(std::ostream &strm,
                                    const TriangleFace_T_N &a);
    using TriangleFace_T::TriangleFace_T;
};

class Mesh : public Object {
    // used for reading in Mesh
    int num_vertices = 0;
    point *vertices = NULL;

    int num_texture_coords = 0;
    double *texture_coords = NULL;

    int num_normals;
    point *normals = NULL;

    int num_faces = 0;
    PrimitiveData *faces = NULL;

    point bary_coords;

    

    TriangleFace* buildFace(std::string& line);

    std::string dir;
    std::list<std::string> mtl_file_list;
    std::list<std::string> mtl_list;

   public:
    BVH bvh;
    using Object::Object;
    void setMesh(const std::string& filename);
    void setMaterial(const std::string& mtllib_line);
    void set_canonical_bounds();
    void intersect(struct Ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);
};

//I will probably delete these
class normalFunctor{
  public:
    virtual point operator()(point &bary_coords) = 0;
};

class FaceNormal : public normalFunctor {
    point operator()(point &bary_coords);
};

class InterpolateNormal : public normalFunctor {
    point operator()(point &bary_coords);
};
#endif