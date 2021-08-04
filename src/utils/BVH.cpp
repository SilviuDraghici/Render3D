#include "BVH.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <stack>

#include "ray.h"

//debug
double num_intersection_tests;
double num_bvh_searches;

int num_buckets = 20;
struct bucket {
    double count = 0;
    Bounds bound;
};

inline void union_bounds(Primitive *a, Primitive *b, BoundingBox *union_box) {
    //assigns the bounds of union_box to fit a and b
    union_box->b_min_x = MIN(a->min_x(), b->min_x());
    union_box->b_min_y = MIN(a->min_y(), b->min_y());
    union_box->b_min_z = MIN(a->min_z(), b->min_z());
    union_box->b_max_x = MAX(a->max_x(), b->max_x());
    union_box->b_max_y = MAX(a->max_y(), b->max_y());
    union_box->b_max_z = MAX(a->max_z(), b->max_z());
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

double BoundingBox::intersect(struct Ray *ray, double lambda) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical box.

    double tmin = -INFINITY;
    double tmax = INFINITY;

    double _tmin;
    double _tmax;
    
    // x axis slab
    _tmin = (min_x() - ray->p0.x) / ray->d.x;
    _tmax = (max_x() - ray->p0.x) / ray->d.x;
    if (_tmin > _tmax) std::swap(_tmin, _tmax);
    if (_tmin > tmin) tmin = _tmin;
    if (_tmax < tmax) tmax = _tmax;

    // y axis slab
    _tmin = (min_y() - ray->p0.y) / ray->d.y;
    _tmax = (max_y() - ray->p0.y) / ray->d.y;
    if (_tmin > _tmax) std::swap(_tmin, _tmax);
    if (tmin > _tmax || _tmin > tmax) return INFINITY;
    if (_tmin > tmin) tmin = _tmin;
    if (_tmax < tmax) tmax = _tmax;

    // z axis slab
    _tmin = (min_z() - ray->p0.z) / ray->d.z;
    _tmax = (max_z() - ray->p0.z) / ray->d.z;
    if (_tmin > _tmax) std::swap(_tmin, _tmax);
    if (tmin > _tmax || _tmin > tmax) return INFINITY;
    if (_tmin > tmin) tmin = _tmin;
    if (_tmax < tmax) tmax = _tmax;

    if (tmin > lambda) return INFINITY;
    if (tmin < THR && tmax < THR) return INFINITY;
    if (tmin < THR) tmin = THR;
    
    return tmin;
}

bool BoundingBox::isprim() { return false; }

Primitive *midsplit(PrimitiveData prims[], int start, int end) {
    int num_prims = end - start;
    if (num_prims == 1) {
        return prims[start].mprimitive;
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

Primitive *SAH_split(PrimitiveData prims[], int start, int end) {
    int num_prims = end - start;
    if (num_prims <= 4) {
        return midsplit(prims, start, end);
    }

    BoundingBox *b_box = new BoundingBox;

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
        split_cost[i] = b0.surface_area() * count0 + b1.surface_area() * count1;
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
        search_ptr = &BVH::bfs;
    } else if (search_method == SearchMethod::DFS) {
        search_ptr = &BVH::dfs;
    }
}

Primitive *BVH::bfs(struct Ray *ray) {
    double lambda = INFINITY;
    Primitive *closest_prim = NULL;

    std::priority_queue<Search_Node, std::vector<Search_Node>, std::greater<Search_Node>> bvh_queue;

    lambda = root->intersect(ray, lambda);
    if (lambda < INFINITY) {
        bvh_queue.emplace(lambda, root);

        //std::cout << "root:" << lambda << "\n";

        // debug :
        num_bvh_searches++;
        num_intersection_tests++;
    }

    Primitive *currentNode;
    BoundingBox *currentBox;
    while (!bvh_queue.empty()) {
        currentNode = bvh_queue.top().node;
        bvh_queue.pop();

        if (currentNode->isprim()) {
            closest_prim = currentNode;
            break;
        }

        // debug :
        num_intersection_tests += 2;

        currentBox = ((BoundingBox *)currentNode);

        lambda = INFINITY;
        lambda = currentBox->c1->intersect(ray, lambda);
        if (lambda < INFINITY) {
            if (lambda < THR) std::cout << ((Object *)currentBox->c1)->name << " " << lambda << "\n";
            bvh_queue.emplace(lambda, currentBox->c1);
        }

        lambda = INFINITY;
        lambda = currentBox->c2->intersect(ray, lambda);
        if (lambda < INFINITY) {
            if (lambda < THR) std::cout << ((Object *)currentBox->c2)->name << " " << lambda << "\n";
            bvh_queue.emplace(lambda, currentBox->c2);
        }
    }
    return closest_prim;
}

Primitive *BVH::dfs(struct Ray *ray) {
    double lambda = INFINITY, l1, l2;
    Primitive *closest_prim = NULL;

    std::stack<Search_Node> bvh_stack;

    l1 = root->intersect(ray, l1);
    if (l1 < INFINITY) {
        bvh_stack.emplace(l1, root);

        // debug :
        num_bvh_searches++;
        num_intersection_tests++;
    }

    Primitive *currentNode;
    BoundingBox *currentBox;
    while (!bvh_stack.empty()) {
        currentNode = bvh_stack.top().node;
        l1 = bvh_stack.top().lambda;
        bvh_stack.pop();

        if (l1 < lambda) {
            if (currentNode->isprim()) {
                closest_prim = currentNode;
                lambda = l1;
            } else {
                // debug :
                num_intersection_tests += 2;

                currentBox = ((BoundingBox *)currentNode);

                l1 = INFINITY;
                l1 = currentBox->c1->intersect(ray, l1);

                l2 = INFINITY;
                l2 = currentBox->c2->intersect(ray, l2);

                if (l1 < l2) {
                    if (l2 < lambda) bvh_stack.emplace(l2, currentBox->c2);
                    if (l1 < lambda) bvh_stack.emplace(l1, currentBox->c1);
                } else {
                    if (l1 < lambda) bvh_stack.emplace(l1, currentBox->c1);
                    if (l2 < lambda) bvh_stack.emplace(l2, currentBox->c2);
                }
            }
        }
    }

    return closest_prim;
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

void printNode(Primitive *node, int space) {
    if (node->isprim()) {
        for (int i = 0; i < space; i++) {
            std::cout << "               ";
        }
        point min = point(node->min_x(), node->min_y(), node->min_z());
        point max = point(node->max_x(), node->max_y(), node->max_z());
        std::cout << ((Object *)node)->name << " " << min << " " << max << "\n";
        return;
    }

    printNode(((BoundingBox *)node)->c1, space + 1);

    for (int i = 0; i < space; i++) {
        std::cout << "                   ";
    }
    point min = point(node->min_x(), node->min_y(), node->min_z());
    point max = point(node->max_x(), node->max_y(), node->max_z());
    std::cout << min << " " << max << "\n";

    printNode(((BoundingBox *)node)->c2, space + 1);
}

void BVH::print() {
    if (root != NULL)
        printNode(root, 0);
}