#include "BVH.h"

#include <algorithm>
#include <iostream>
#include <queue>

#include "ray.h"

//debug
double num_intersection_tests;
double num_intersect_calls;

int num_buckets = 20;
struct bucket {
    double count = 0;
    Bounds bound;
};

Axis BVH_Node::longestAxis() const {
    double x = max_x() - min_x();
    double y = max_y() - min_y();
    double z = max_z() - min_z();
    if (x > y && x > z) {
        return Axis::X;
    } else if (y > z) {
        return Axis::Y;
    } else {
        return Axis::Z;
    }
}

inline void union_bounds(BVH_Node *a, BVH_Node *b, BoundingBox *union_box) {
    //assigns the bounds of union_box to fit a and b
    union_box->b_min_x = MIN(a->min_x(), b->min_x());
    union_box->b_min_y = MIN(a->min_y(), b->min_y());
    union_box->b_min_z = MIN(a->min_z(), b->min_z());
    union_box->b_max_x = MAX(a->max_x(), b->max_x());
    union_box->b_max_y = MAX(a->max_y(), b->max_y());
    union_box->b_max_z = MAX(a->max_z(), b->max_z());
}

void union_bounds(Bounds &a, point &b, Bounds &union_box) {
    union_box.min.x = MIN(a.min.x, b.x);
    union_box.min.y = MIN(a.min.y, b.y);
    union_box.min.z = MIN(a.min.z, b.z);
    union_box.max.x = MAX(a.max.x, b.x);
    union_box.max.y = MAX(a.max.y, b.y);
    union_box.max.z = MAX(a.max.z, b.z);
}

void union_bounds(Bounds &a, Bounds &b, Bounds &union_box) {
    //assigns the bounds of union_box to fit a and b
    union_box.min.x = MIN(a.min.x, b.min.x);
    union_box.min.y = MIN(a.min.y, b.min.y);
    union_box.min.z = MIN(a.min.z, b.min.z);
    union_box.max.x = MAX(a.max.x, b.max.x);
    union_box.max.y = MAX(a.max.y, b.max.y);
    union_box.max.z = MAX(a.max.z, b.max.z);
}

Axis Bounds::longestAxis() const {
    double x = max.x - min.x;
    double y = max.y - min.y;
    double z = max.z - min.z;
    if (x > y && x > z) {
        return Axis::X;
    } else if (y > z) {
        return Axis::Y;
    } else {
        return Axis::Z;
    }
}

double Bounds::offset(const point &p, Axis dim) const {
    // returns how far along the point p is along the bounding box in the given dimension
    // i.e 0 if p[dim] == min[dim] and 1 if p[dim] == max[dim]
    double pdim = p[dim] - min[dim];
    return pdim / (max[dim] - min[dim]);
}

double Bounds::surface_area() const {
    double x = max.x - min.x;
    double y = max.y - min.y;
    double z = max.z - min.z;
    return 2 * ((x * y) + (x * z) + (y * z));
}

double BoundingBox::min_x() const { return b_min_x; }
double BoundingBox::min_y() const { return b_min_y; }
double BoundingBox::min_z() const { return b_min_z; }
double BoundingBox::max_x() const { return b_max_x; }
double BoundingBox::max_y() const { return b_max_y; }
double BoundingBox::max_z() const { return b_max_z; }

BVH_Node *BoundingBox::intersect(struct ray *ray, double *lambda,
                                 point *bary_coords) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical box.

    // debug :
    num_intersection_tests++;

    if ((b_min_x < ray->p0.x && ray->p0.x < b_max_x) &&
        (b_min_y < ray->p0.y && ray->p0.y < b_max_y) &&
        (b_min_z < ray->p0.z && ray->p0.z < b_max_z)) {
        *lambda = THR;
        return NULL;
    }

    point p;

    // current intersection lambda
    double b_lambda;

    // y-z plane box face at min_x
    b_lambda = (b_min_x - ray->p0.x) / ray->d.x;
    if (THR < b_lambda) {
        rayPosition(ray, b_lambda, &p);
        if ((b_min_y < p.y && p.y < b_max_y) &&
            (b_min_z < p.z && p.z < b_max_z)) {
            *lambda = b_lambda;
        }
    }

    // y-z plane box face at max_x
    b_lambda = (b_max_x - ray->p0.x) / ray->d.x;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(ray, b_lambda, &p);
        if ((b_min_y < p.y && p.y < b_max_y) &&
            (b_min_z < p.z && p.z < b_max_z)) {
            *lambda = b_lambda;
        }
    }

    // x-z plane box face at min_y
    b_lambda = (b_min_y - ray->p0.y) / ray->d.y;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(ray, b_lambda, &p);
        if ((b_min_x < p.x && p.x < b_max_x) &&
            (b_min_z < p.z && p.z < b_max_z)) {
            *lambda = b_lambda;
        }
    }

    // x-z plane box face at max_y
    b_lambda = (b_max_y - ray->p0.y) / ray->d.y;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(ray, b_lambda, &p);
        if ((b_min_x < p.x && p.x < b_max_x) &&
            (b_min_z < p.z && p.z < b_max_z)) {
            *lambda = b_lambda;
        }
    }

    // x-y plane box face at min_z
    b_lambda = (b_min_z - ray->p0.z) / ray->d.z;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(ray, b_lambda, &p);
        if ((b_min_x < p.x && p.x < b_max_x) &&
            (b_min_y < p.y && p.y < b_max_y)) {
            *lambda = b_lambda;
        }
    }

    // x-y plane box face at max_z
    b_lambda = (b_max_z - ray->p0.z) / ray->d.z;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(ray, b_lambda, &p);
        if ((b_min_x < p.x && p.x < b_max_x) &&
            (b_min_y < p.y && p.y < b_max_y)) {
            *lambda = b_lambda;
        }
    }

    return NULL;
}

bool BoundingBox::isFace() { return false; }

BVH_Node *midsplit(PrimitiveData prims[], int start, int end) {
    int num_prims = end - start;
    if (num_prims == 1) {
        return prims[start].primitive;
    }

    BoundingBox *b_box = new BoundingBox;

    Bounds centroid_bounds;
    for (int i = start; i < end; i++) {
        union_bounds(centroid_bounds, prims[i].center, centroid_bounds);
    }
    Axis max_axis = centroid_bounds.longestAxis();

    double mid = (centroid_bounds.min[max_axis] + centroid_bounds.max[max_axis]) / 2;

    PrimitiveData *mid_ptr = std::partition(&prims[start], &prims[end],
                                            [max_axis, mid](const PrimitiveData &p) {
                                                return (p.center[max_axis] < mid);
                                            });

    //find index of first primitive greater than midpoint
    int split = mid_ptr - &prims[0];
    if (split == start || split == end) {
        split = (start + end) / 2;
    }
    //std::cout << start << " " << split << " " << end << "\n";

    b_box->c1 = midsplit(prims, start, split);
    b_box->c2 = midsplit(prims, split, end);

    //set bounds for bounding box
    union_bounds(b_box->c1, b_box->c2, b_box);

    return b_box;
}

BVH_Node *SAH_split(PrimitiveData prims[], int start, int end) {
    int num_prims = end - start;
    if (num_prims <= 4) {
        return midsplit(prims, start, end);
    }

    BoundingBox *b_box = new BoundingBox;

    Bounds bound;
    for (int i = start; i < end; i++) {
        union_bounds(bound, prims[i].bound, bound);
    }

    Bounds centroid_bounds;
    for (int i = start; i < end; i++) {
        union_bounds(centroid_bounds, prims[i].center, centroid_bounds);
    }
    Axis max_axis = centroid_bounds.longestAxis();

    struct bucket buckets[num_buckets];

    int b;
    for (int i = start; i < end; i++) {
        b = MIN(num_buckets - 1, num_buckets * centroid_bounds.offset(prims[i].center, max_axis));
        buckets[b].count++;
        union_bounds(buckets[b].bound, prims[i].bound, buckets[b].bound);
    }

    double split_cost[num_buckets - 1];
    for (int i = 0; i < num_buckets - 1; i++) {
        Bounds b0, b1;
        int count0 = 0, count1 = 0;

        for (int j = 0; j <= i; j++) {
            union_bounds(b0, buckets[j].bound, b0);
            count0 += buckets[j].count;
        }

        for (int j = i + 1; j < num_buckets; j++) {
            union_bounds(b1, buckets[j].bound, b1);
            count1 += buckets[j].count;
        }
        split_cost[i] = 0.125 + (b0.surface_area() * count0 + b1.surface_area() * count1) / bound.surface_area();
    }

    double min_cost = split_cost[0];
    int min_cost_split = 0;
    for (int i = 1; i < num_buckets - 1; i++) {
        if (split_cost[i] < min_cost) {
            min_cost = split_cost[i];
            min_cost_split = i;
        }
    }

    PrimitiveData *mid_ptr = std::partition(&prims[start], &prims[end],
                                            [=](const PrimitiveData &p) {
                                                int b = MIN(num_buckets - 1, num_buckets * centroid_bounds.offset(p.center, max_axis));
                                                return b <= min_cost_split;
                                            });

    //find index of first primitive greater than midpoint
    int split = mid_ptr - &prims[0];
    if (split == start || split == end) {
        split = (start + end) / 2;
    }

    b_box->c1 = SAH_split(prims, start, split);
    b_box->c2 = SAH_split(prims, split, end);

    //set bounds for bounding box
    union_bounds(b_box->c1, b_box->c2, b_box);

    return b_box;
}

void BVH::set_build_method(BuildMethod split_method) {
    if (split_method == BuildMethod::SAH) {
        build_ptr = &BVH::SAH_build;
    } else if (split_method == BuildMethod::MidSplit) {
        build_ptr = &BVH::midsplit_build;
    }
}

void BVH::SAH_build(PrimitiveData prims[], int num_prims) {
    root = SAH_split(prims, 0, num_prims);
}

void BVH::midsplit_build(PrimitiveData prims[], int num_prims) {
    root = midsplit(prims, 0, num_prims);
}

void BVH::set_search_method(SearchMethod search_method) {
    if (search_method == SearchMethod::BFS) {
        search_ptr = &BVH::bfs_pqueue_search;
    } else if (search_method == SearchMethod::DFS) {
        search_ptr = &BVH::dfs_search;
    }
}

BVH_Node *BVH::bfs_pqueue_search(struct ray *ray) {
    //std::cout << "pqueue search\n";
    double lambda = INFINITY;
    point bary_coords;
    BVH_Node *closest_prim = NULL;

    std::priority_queue<PQ_Node, std::vector<PQ_Node>, std::greater<PQ_Node>> bvh_queue;

    root->intersect(ray, &lambda, &bary_coords);
    if (lambda < INFINITY) {
        bvh_queue.emplace(lambda, root);

        // debug :
        num_intersect_calls++;
    }

    BVH_Node *currentNode;
    BoundingBox *currentBox;
    while (!bvh_queue.empty()) {
        currentNode = bvh_queue.top().node;
        bvh_queue.pop();

        if (currentNode->isFace()) {
            closest_prim = currentNode;
            break;
        }

        currentBox = ((BoundingBox *)currentNode);

        lambda = INFINITY;
        currentBox->c1->intersect(ray, &lambda, &bary_coords);
        if (lambda < INFINITY) {
            bvh_queue.emplace(lambda, currentBox->c1);
        }

        lambda = INFINITY;
        currentBox->c2->intersect(ray, &lambda, &bary_coords);
        if (lambda < INFINITY) {
            bvh_queue.emplace(lambda, currentBox->c2);
        }
    }
    return closest_prim;
}

BVH_Node *BVH::dfs_search(struct ray *ray) {
    //todo
    std::cerr << "DFS search not implemented!!!\n";
    return NULL;
}

int comp_face_max_x(const void *a, const void *b) {
    PrimitiveData *pa = (PrimitiveData *)a, *pb = (PrimitiveData *)b;
    double result = pa->bound.max.x - pb->bound.max.x;
    return result > 0;
}

int comp_face_max_y(const void *a, const void *b) {
    PrimitiveData *pa = (PrimitiveData *)a, *pb = (PrimitiveData *)b;
    double result = pa->bound.max.y - pb->bound.max.y;
    return result > 0;
}

int comp_face_max_z(const void *a, const void *b) {
    PrimitiveData *pa = (PrimitiveData *)a, *pb = (PrimitiveData *)b;
    double result = pa->bound.max.z - pb->bound.max.z;
    return result > 0;
}