#include "BVH.h"

#include <iostream>
#include <queue>

#include "ray.h"

//debug
double num_intersection_tests;
double num_intersect_calls;

double max_split_checks = 100;

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

    //set bounds for bounding box
    double x, y, z;
    b_box->b_min_x = b_box->b_min_y = b_box->b_min_z = INFINITY;
    b_box->b_max_x = b_box->b_max_y = b_box->b_max_z = -INFINITY;
    for (int i = start; i < end; i++) {
        x = prims[i].min.x;
        if (x < b_box->b_min_x) {
            b_box->b_min_x = x;
        }
        x = prims[i].max.x;
        if (b_box->b_max_x < x) {
            b_box->b_max_x = x;
        }

        y = prims[i].min.y;
        if (y < b_box->b_min_y) {
            b_box->b_min_y = y;
        }
        y = prims[i].max.y;
        if (b_box->b_max_y < y) {
            b_box->b_max_y = y;
        }

        z = prims[i].min.z;
        if (z < b_box->b_min_z) {
            b_box->b_min_z = z;
        }
        z = prims[i].max.z;
        if (b_box->b_max_z < z) {
            b_box->b_max_z = z;
        }
    }

    int split = (start + end) / 2;
    b_box->c1 = midsplit(prims, start, split);
    b_box->c2 = midsplit(prims, split, end);

    return b_box;
}

void BVH::set_build_method(BuildMethod split_method){
    if(split_method == BuildMethod::SAH){
        build_ptr = &BVH::SAH_build;
    } else if(split_method == BuildMethod::MidSplit){
        build_ptr = &BVH::midsplit_build;
    }
}

void BVH::SAH_build(PrimitiveData prims[], int num_prims){

}
void BVH::midsplit_build(PrimitiveData prims[], int num_prims) {
    int (*compare)(const void *, const void *) = comp_face_max_x;
    /*if (abs(max_y - min_y)) {
        compare = comp_face_max_y;
    } else if (scale == abs(max_z - min_z)) {
        compare = comp_face_max_z;
    }*/
    qsort(prims, num_prims, sizeof(PrimitiveData), compare);
    root = midsplit(prims, 0, num_prims);
}

void BVH::set_search_method(SearchMethod search_method){
    if(search_method == SearchMethod::BFS){
        search_ptr = &BVH::bfs_pqueue_search;
    } else if(search_method == SearchMethod::DFS){
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

BVH_Node *BVH::dfs_search(struct ray *ray){
    //todo
    std::cerr << "DFS search not implemented!!!\n";
    return NULL;
}

int comp_face_max_x(const void *a, const void *b) {
    PrimitiveData *pa = (PrimitiveData *)a, *pb = (PrimitiveData *)b;
    double result = pa->max.x - pb->max.x;
    return result > 0;
}

int comp_face_max_y(const void *a, const void *b) {
    PrimitiveData *pa = (PrimitiveData *)a, *pb = (PrimitiveData *)b;
    double result = pa->max.y - pb->max.y;
    return result > 0;
}

int comp_face_max_z(const void *a, const void *b) {
    PrimitiveData *pa = (PrimitiveData *)a, *pb = (PrimitiveData *)b;
    double result = pa->max.z - pb->max.z;
    return result > 0;
}