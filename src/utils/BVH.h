#ifndef BVH_H
#define BVH_H

#include "objects.h"

//debug
extern double num_intersection_tests;
extern double num_bvh_searches;

enum class BuildMethod { SAH,
                         MidSplit };

enum class SearchMethod { DFS,
                          BFS };

struct PrimitiveData {
    Primitive *mprimitive;
    Bounds bound;
    point center;
    PrimitiveData &operator=(const Primitive *prim) {
        mprimitive = (Primitive *)prim;
        bound.min.x = prim->min_x();
        bound.min.y = prim->min_y();
        bound.min.z = prim->min_z();

        bound.max.x = prim->max_x();
        bound.max.y = prim->max_y();
        bound.max.z = prim->max_z();

        center = 0.5 * bound.min + 0.5 * bound.max;
        return *this;
    }
};

class BVH {

    //build methods
    void SAH_build(PrimitiveData prims[], int num_prims);
    void midsplit_build(PrimitiveData prims[], int num_prims);

    typedef void (BVH::*BVH_build)(PrimitiveData prims[], int num_prims);
    BVH_build build_ptr = &BVH::SAH_build;

    //search methods
    Primitive *bfs(struct Ray *ray);
    Primitive *dfs(struct Ray *ray);

    typedef Primitive *(BVH::*BVH_search)(struct Ray *ray);
    BVH_search search_ptr = &BVH::bfs;

   public:
    Primitive *root = NULL;
    void set_build_method(BuildMethod split_method);
    void print();
    inline void build(PrimitiveData prims[], int num_prims) {
        // build the BVH using the selected split method
        (this->*build_ptr)(prims, num_prims);
    }

    void set_search_method(SearchMethod search_method);
    inline Primitive *search(struct Ray *ray) {
        // search for the first hit primitve in the BVH using the
        // selected search method
        return (this->*search_ptr)(ray);
    }
};

class BoundingBox : public Primitive {
    int depth;

   public:
    Primitive *c1, *c2;
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
    double intersect(struct Ray *r, double lambda);
    bool isprim();
};

class Search_Node {
   public:
    double lambda;
    Primitive *node;
    Search_Node(double l, Primitive *n) {
        lambda = l;
        node = n;
    }
    Search_Node(double l) {
        lambda = l;
        node = NULL;
    }

    friend bool operator<(const Search_Node &lhs, const Search_Node &rhs) {
        return lhs.lambda < rhs.lambda;
    }
    friend bool operator>(const Search_Node &lhs, const Search_Node &rhs) {
        return lhs.lambda > rhs.lambda;
    }

    friend std::ostream &operator<<(std::ostream &strm, const Search_Node &a) {
        return strm << a.lambda;
    }
};

int comp_face_max_x(const void *a, const void *b);
int comp_face_max_y(const void *a, const void *b);
int comp_face_max_z(const void *a, const void *b);

#endif  // ifdef BVH_H