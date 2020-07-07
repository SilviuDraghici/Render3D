#ifndef BVH_H
#define BVH_H

#include "objects.h"

#define BB BoundingBox

//debug
extern double num_intersection_tests;
extern double num_intersect_calls;

enum class BuildMethod { SAH,
                         MidSplit };

enum class SearchMethod { DFS,
                          BFS };

class BVH_Node {
   public:
    virtual BVH_Node *intersect(struct ray *r, double *lambda,
                                point *bary_coords) = 0;
    virtual bool isFace() = 0;
    virtual double min_x() const = 0;
    virtual double min_y() const = 0;
    virtual double min_z() const = 0;
    virtual double max_x() const = 0;
    virtual double max_y() const = 0;
    virtual double max_z() const = 0;
};

struct PrimitiveData {
    BVH_Node *primitive;
    point min, max, center;
    PrimitiveData &operator=(const BVH_Node *prim) {
        primitive = (BVH_Node *)prim;
        min.x = prim->min_x();
        min.y = prim->min_y();
        min.z = prim->min_z();

        max.x = prim->max_x();
        max.y = prim->max_y();
        max.z = prim->max_z();

        center = 0.5 * min + 0.5 * max;
        return *this;
    }
};

class BVH {
    BVH_Node *root = NULL;

    //build methods
    void SAH_build(PrimitiveData prims[], int num_prims);
    void midsplit_build(PrimitiveData prims[], int num_prims);
    
    typedef void (BVH::*BVH_build)(PrimitiveData prims[], int num_prims);
    BVH_build build_ptr = &BVH::SAH_build;

    //search methods
    BVH_Node *bfs_pqueue_search(struct ray *ray);
    BVH_Node *dfs_search(struct ray *ray);

    typedef BVH_Node *(BVH::*BVH_search)(struct ray *ray);
    BVH_search search_ptr = &BVH::bfs_pqueue_search;

   public:

    void set_build_method(BuildMethod split_method);
    inline void build(PrimitiveData prims[], int num_prims) {
        // build the BVH using the selected split method
        (this->*build_ptr)(prims, num_prims);
    }

    void set_search_method(SearchMethod search_method);
    inline BVH_Node *search(struct ray *ray) {
        // search for the first hit primitve in the BVH using the
        // selected search method
        return (this->*search_ptr)(ray);
    }
};

class BoundingBox : public BVH_Node {
    int depth;

   public:
    BVH_Node *c1, *c2;
    double b_min_x, b_max_x;
    double b_min_y, b_max_y;
    double b_min_z, b_max_z;
    void setBounds(double min_x, double max_x,
                   double min_y, double max_y,
                   double min_z, double max_z);

    double min_x() const;
    double min_y() const;
    double min_z() const;
    double max_x() const;
    double max_y() const;
    double max_z() const;
    BVH_Node *intersect(struct ray *r, double *lambda, point *bary_coords);
    bool isFace();
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

#endif  // ifdef BVH_H