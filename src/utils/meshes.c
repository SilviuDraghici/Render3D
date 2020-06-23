#include "meshes.h"

#include <fstream>
#include <iostream>
#include <queue>
#include <limits>

#include "ray.h"

#define BB BoundingBox
//#define BB BoundingBox_Visible

//#define DEBUG

double num_intersection_tests;
double num_intersect_calls;

TriangleFace::TriangleFace() {}

TriangleFace::TriangleFace(point p1, point p2, point p3) {
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
}

TriangleFace::TriangleFace(double x1, double y1, double z1, double x2,
                           double y2, double z2, double x3, double y3,
                           double z3) {
    p1.x = x1, p1.y = y1, p1.z = z1;
    p2.x = x2, p2.y = y2, p2.z = z2;
    p3.x = x3, p3.y = y3, p3.z = z3;
}

double TriangleFace::min_x() { return MIN(p1.x, MIN(p2.x, p3.x)); }
double TriangleFace::min_y() { return MIN(p1.y, MIN(p2.y, p3.y)); }
double TriangleFace::min_z() { return MIN(p1.z, MIN(p2.z, p3.z)); }
double TriangleFace::max_x() { return MAX(p1.x, MAX(p2.x, p3.x)); }
double TriangleFace::max_y() { return MAX(p1.y, MAX(p2.y, p3.y)); }
double TriangleFace::max_z() { return MAX(p1.z, MAX(p2.z, p3.z)); }

BVH_Node *TriangleFace::intersect(struct ray *ray, double *lambda,
                                  point *bary_coords) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified triangle.

    // debug :
    num_intersection_tests++;

    point e12 = p2 - p1;
    point e23 = p3 - p2;

    point normal = cross(&e12, &e23);
    double denom = dot(&normal, &normal);

    // ray plane intersection calculation
    point pi, r = p1 - ray->p0;
    double d_dot_n = dot(&ray->d, &normal);
    double r_dot_n = dot(&r, &normal);
    double t_lambda = r_dot_n / d_dot_n;

    if (t_lambda < THR)  //|| *lambda < t_lambda)
        return NULL;     // intersection with plane is negative

    double u, v, w;

    rayPosition(ray, t_lambda, &pi);

    // e12 has already been calculated
    point v1i = pi - p1;
    point crossp = cross(&e12, &v1i);
    u = dot(&crossp, &normal);
    if (u < 0) return NULL;  // outside first edge

    // e23 has already been calculated
    point v2i = pi - p2;
    crossp = cross(&e23, &v2i);
    v = dot(&crossp, &normal);
    if (v < 0) return NULL;  // outside second edge

    point e31 = p1 - p3;
    point v3i = pi - p3;
    crossp = cross(&e31, &v3i);
    w = dot(&crossp, &normal);
    if (w < 0) return NULL;  // outside third edge

    // intersection is within triangle
    bary_coords->z = u / denom;
    bary_coords->x = v / denom;
    bary_coords->y = w / denom;
    *lambda = t_lambda;
    return this;
}

bool TriangleFace::isFace() { return true; }

point TriangleFace::normal(point *bary_coords) {
    point e12 = p2 - p1;
    point e23 = p3 - p2;

    return cross(&e12, &e23);
}

std::ostream &operator<<(std::ostream &strm, const TriangleFace &a) {
    return strm << "p1: " << a.p1 << " p2: " << a.p2 << " p3: " << a.p3;
}

void TriangleFace_N::setNormals(point n1, point n2, point n3) {
    this->n1 = n1, this->n2 = n2, this->n3 = n3;
}

point TriangleFace_N::normal(point *bary_coords) {
    // std::cout << "n1: " << n1 << " n2: " << n2 << " n3: " << n3 << "\n";
    return n1 * bary_coords->x + n2 * bary_coords->y + n3 * bary_coords->z;
}

std::ostream &operator<<(std::ostream &strm, const TriangleFace_N &a) {
    return strm << "p1: " << a.p1 << " p2: " << a.p2 << " p3: " << a.p3
                << "n1: " << a.n1 << " n2: " << a.n2 << " n3: " << a.n3;
}

void BoundingBox::setBounds(double min_x, double max_x, double min_y,
                            double max_y, double min_z, double max_z) {
    this->b_min_x = min_x, this->b_max_x = max_x;
    this->b_min_y = min_y, this->b_max_y = max_y;
    this->b_min_z = min_z, this->b_max_z = max_z;
}

void BoundingBox::setChildren(TriangleFace_N *faces, int start, int end) {
    // std::cout << "start: " << start << " end: " << end << "\n";
    double x, y, z;
    b_min_x = b_min_y = b_min_z = INFINITY;
    b_max_x = b_max_y = b_max_z = -INFINITY;
    for (int i = start; i < end; i++) {
        x = faces[i].min_x();
        if (x < b_min_x) {
            b_min_x = x;
        }
        x = faces[i].max_x();
        if (b_max_x < x) {
            b_max_x = x;
        }

        y = faces[i].min_y();
        if (y < b_min_y) {
            b_min_y = y;
        }
        y = faces[i].max_y();
        if (b_max_y < y) {
            b_max_y = y;
        }

        z = faces[i].min_z();
        if (z < b_min_z) {
            b_min_z = z;
        }
        z = faces[i].max_z();
        if (b_max_z < z) {
            b_max_z = z;
        }
    }

    if (end - start == 2) {
        c1 = &faces[start];
        c2 = &faces[start + 1];
    } else {
        int mid = (start + end) / 2;
        if (mid - start == 1) {
            c1 = &faces[start];
        } else {
            c1 = new BB;
            ((BB *)c1)->setChildren(faces, start, mid);
        }

        if (end - mid == 1) {
            c1 = &faces[mid];
        } else {
            c2 = new BB;
            ((BB *)c2)->setChildren(faces, mid, end);
        }
    }
}

double BoundingBox::min_x() { return b_min_x; }
double BoundingBox::min_y() { return b_min_y; }
double BoundingBox::min_z() { return b_min_z; }
double BoundingBox::max_x() { return b_max_x; }
double BoundingBox::max_y() { return b_max_y; }
double BoundingBox::max_z() { return b_max_z; }

BVH_Node *BoundingBox::intersect(struct ray *ray, double *lambda,
                                 point *bary_coords) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical box.

    // debug :
    num_intersection_tests++;

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

BoundingBox_Visible::BoundingBox_Visible() {
    col = color(drand48(), drand48(), drand48());
}

color BoundingBox_Visible::getCol() { return col; }

BVH_Node *BoundingBox_Visible::intersect(struct ray *ray, double *lambda,
                                         point *bary_coords) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical box.

    // debug :
    num_intersection_tests++;

    point p;
    color col = getCol();
    bary_coords->x = col.R;
    bary_coords->y = col.G;
    bary_coords->z = col.B;
    // current intersection lambda
    double b_lambda;
    double a, b;
    // y-z plane box face at min_x
    b_lambda = (b_min_x - ray->p0.x) / ray->d.x;
    rayPosition(ray, b_lambda, &p);
    if (THR < b_lambda && (b_min_y < p.y && p.y < b_max_y) &&
        (b_min_z < p.z && p.z < b_max_z)) {
        a = (p.y - b_min_y) / (b_max_y - b_min_y);
        b = (p.z - b_min_z) / (b_max_z - b_min_z);

        if ((0 <= a && a <= 0 + width || 1 - width <= a && a <= 1) ||
            (0 <= b && b <= 0 + width || 1 - width <= b && b <= 1)) {
            *lambda = b_lambda;
            bary_coords->w = -1;
#if 0
                printf("hit box %f\n", *lambda);
#endif
            return NULL;
        }
        *lambda = b_lambda;
    }

    // y-z plane box face at max_x
    b_lambda = (b_max_x - ray->p0.x) / ray->d.x;
    rayPosition(ray, b_lambda, &p);
    if (THR < b_lambda && (b_min_y < p.y && p.y < b_max_y) &&
        (b_min_z < p.z && p.z < b_max_z)) {
        a = (p.y - b_min_y) / (b_max_y - b_min_y);
        b = (p.z - b_min_z) / (b_max_z - b_min_z);
        if ((0 <= a && a <= 0 + width || 1 - width <= a && a <= 1) ||
            (0 <= b && b <= 0 + width || 1 - width <= b && b <= 1)) {
            *lambda = b_lambda;
            bary_coords->w = -1;
#if 0
                printf("hit box %f\n", *lambda);
#endif
            return NULL;
        }
        if (b_lambda < *lambda) {
            *lambda = b_lambda;
        }
    }

    // x-z plane box face at min_y
    b_lambda = (b_min_y - ray->p0.y) / ray->d.y;
    rayPosition(ray, b_lambda, &p);
    if (THR < b_lambda && (b_min_x < p.x && p.x < b_max_x) &&
        (b_min_z < p.z && p.z < b_max_z)) {
        a = (p.x - b_min_x) / (b_max_x - b_min_x);
        b = (p.z - b_min_z) / (b_max_z - b_min_z);
        if ((0 <= a && a <= 0 + width || 1 - width <= a && a <= 1) ||
            (0 <= b && b <= 0 + width || 1 - width <= b && b <= 1)) {
            *lambda = b_lambda;
            bary_coords->w = -1;
#if 0
                printf("hit box %f\n", *lambda);
#endif
            return NULL;
        }
        if (b_lambda < *lambda) {
            *lambda = b_lambda;
        }
    }

    // x-z plane box face at max_y
    b_lambda = (b_max_y - ray->p0.y) / ray->d.y;
    rayPosition(ray, b_lambda, &p);
    if (THR < b_lambda && (b_min_x < p.x && p.x < b_max_x) &&
        (b_min_z < p.z && p.z < b_max_z)) {
        a = (p.x - b_min_x) / (b_max_x - b_min_x);
        b = (p.z - b_min_z) / (b_max_z - b_min_z);
        if ((0 <= a && a <= 0 + width || 1 - width <= a && a <= 1) ||
            (0 <= b && b <= 0 + width || 1 - width <= b && b <= 1)) {
            *lambda = b_lambda;
            bary_coords->w = -1;
#if 0
               printf("hit box %f\n", *lambda);
#endif
            return NULL;
        }
        if (b_lambda < *lambda) {
            *lambda = b_lambda;
        }
    }

    // x-y plane box face at min_z
    b_lambda = (b_min_z - ray->p0.z) / ray->d.z;
    rayPosition(ray, b_lambda, &p);
    if (THR < b_lambda && (b_min_x < p.x && p.x < b_max_x) &&
        (b_min_y < p.y && p.y < b_max_y)) {
        a = (p.x - b_min_x) / (b_max_x - b_min_x);
        b = (p.y - b_min_y) / (b_max_y - b_min_y);
        if ((0 <= a && a <= 0 + width || 1 - width <= a && a <= 1) ||
            (0 <= b && b <= 0 + width || 1 - width <= b && b <= 1)) {
            *lambda = b_lambda;
            bary_coords->w = -1;
#if 0
                printf("hit box %f\n", *lambda);
#endif
            return NULL;
        }
        if (b_lambda < *lambda) {
            *lambda = b_lambda;
        }
    }

    // x-y plane box face at max_z
    b_lambda = (b_max_z - ray->p0.z) / ray->d.z;
    rayPosition(ray, b_lambda, &p);
    if (THR < b_lambda && (b_min_x < p.x && p.x < b_max_x) &&
        (b_min_y < p.y && p.y < b_max_y)) {
        a = (p.x - b_min_x) / (b_max_x - b_min_x);
        b = (p.y - b_min_y) / (b_max_y - b_min_y);
        if ((0 <= a && a <= 0 + width || 1 - width <= a && a <= 1) ||
            (0 <= b && b <= 0 + width || 1 - width <= b && b <= 1)) {
            *lambda = b_lambda;
            bary_coords->w = -1;
#if 0
                printf("hit box %f\n", *lambda);
#endif
            return NULL;
        }
        if (b_lambda < *lambda) {
            *lambda = b_lambda;
        }
    }

    if (*lambda < INFINITY) {
        *lambda = INFINITY;
        BVH_Node *c1_result = c1->intersect(ray, lambda, bary_coords);

        if (c2 == NULL || bary_coords->w == -1) {
            return c1_result;
        } else {
            double c_lambda = INFINITY;
            point c_bary_coords;
            BVH_Node *c2_result = c2->intersect(ray, &c_lambda, &c_bary_coords);
#if 0 
                printf("hit box %f %f \n", *lambda, c_lambda);
#endif
            if (c_lambda < *lambda || c_bary_coords.w == -1) {
                *lambda = c_lambda;
                *bary_coords = c_bary_coords;
                return c2_result;
            } else {
#if 0
                printf("else %f %f\n", *lambda, bary_coords->w);
#endif
                return c1_result;
            }
        }
    }

    return NULL;
}

void Mesh::setMesh(const char *filename) {
    frontAndBack = 0;

    box = new BB;

    std::string line;
    std::ifstream mesh_obj(filename);

    std::streampos first_vertex, first_normal, first_face;

    num_vertices = 0, num_faces = 0;

    double x, y, z, scale;
    double min_x, max_x, avg_x, min_y, max_y, avg_y, min_z, max_z, avg_z;
    min_x = min_y = min_z = INFINITY;
    max_x = max_y = max_z = -INFINITY;
    avg_x = avg_y = avg_z = 0;

    if (!mesh_obj.is_open()) {
        std::cout << "Unable to open mesh file: " << filename << "\n";
        return;
    }

    // scrub through comments
    while (getline(mesh_obj, line) && line.rfind("#", 0) == 0) {
        // sscanf(line.c_str(), "# Vertices: %d", &num_vertices);
        sscanf(line.c_str(), "# Faces: %d", &num_faces);
    }

    first_vertex = (int)mesh_obj.tellg() - (line.size() + 1);
    mesh_obj.seekg(first_vertex);

    // count the number of vertices
    while (getline(mesh_obj, line) && line.rfind("v ", 0) == 0) {
        num_vertices++;
        sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
        avg_x += x, avg_y += y, avg_z += z;
        if (x < min_x) {
            min_x = x;
        }
        if (max_x < x) {
            max_x = x;
        }

        if (y < min_y) {
            min_y = y;
        }
        if (max_y < y) {
            max_y = y;
        }

        if (z < min_z) {
            min_z = z;
        }
        if (max_z < z) {
            max_z = z;
        }
    }

    // calculate average diatance from 0 so that mesh can be centered
    avg_x /= num_vertices, avg_y /= num_vertices, avg_z /= num_vertices;
    // calculate max dimension so mesh can be scaled down to 1x1x1 ish
    scale =
        MAX(abs(max_x - min_x), MAX(abs(max_y - min_y), abs(max_z - min_z)));

    num_vertices += 1;  // << overcount by 1
    vertices = (point *)malloc(num_vertices * sizeof(point));

    // read vertices into array
    mesh_obj.seekg(first_vertex);

    // start filling at 1 so face's lookup doesn't need to subtract 1
    for (int i = 1; i < num_vertices; i++) {
        getline(mesh_obj, line);
        sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
        x = (x - avg_x) / scale, y = (y - avg_y) / scale,
        z = (z - avg_z) / scale;
        vertices[i] = point(x, y, z);
    }

    // save location of first normal
    first_normal = mesh_obj.tellg();

    // count number over vertex normals
    num_normals = 1;  // overcount
    while (getline(mesh_obj, line) && line.rfind("vn ", 0) == 0) {
        // std::cout << line << "\n";
        num_normals++;
    }
    if (num_normals > 1) {
        normals = (point *)malloc(num_normals * sizeof(point));
        mesh_obj.seekg(first_normal);
        for (int i = 1; i < num_normals; i++) {
            getline(mesh_obj, line);
            // std::cout << line << "\n";
            sscanf(line.c_str(), "vn %lf %lf %lf", &x, &y, &z);
            normals[i] = point(x, y, z);
        }
    }

    // save location of first face
    first_face = (int)mesh_obj.tellg();

    if (num_faces == 0) {
        // count number of faces
        while (getline(mesh_obj, line) && line.rfind("f ", 0) == 0) {
            num_faces++;
        }
    }

    if (num_normals > 1) {
        faces = new TriangleFace_N[num_faces];
    } else {
        // faces = new TriangleFace[num_faces];
    }

    mesh_obj.clear();
    mesh_obj.seekg(first_face);

    int p1, p2, p3, n1, n2, n3;
    getline(mesh_obj, line);
    if (sscanf(line.c_str(), "f %d %d %d", &p1, &p2, &p3) == 3) {
        mesh_obj.seekg(first_face);
        for (int i = 0; i < num_faces; i++) {
            getline(mesh_obj, line);
            // std::cout << line << "\n";
            sscanf(line.c_str(), "f %d %d %d", &p1, &p2, &p3);
            faces[i] = TriangleFace_N(vertices[p1], vertices[p2], vertices[p3]);
            if (num_normals > 1) {
                ((TriangleFace_N *)faces)[i].setNormals(
                    normals[p1], normals[p2], normals[p3]);
            }
        }
    } else {
        mesh_obj.seekg(first_face);
        for (int i = 0; i < num_faces; i++) {
            getline(mesh_obj, line);
            // std::cout << line << "\n";
            sscanf(line.c_str(), "f %d//%d %d//%d %d//%d", &p1, &n1, &p2, &n2,
                   &p3, &n3);
            faces[i] = TriangleFace_N(vertices[p1], vertices[p2], vertices[p3]);
            ((TriangleFace_N *)faces)[i].setNormals(normals[n1], normals[n2],
                                                    normals[n3]);
        }
    }

    // the faces are sorted by max vertex along one of the axis.
    // which axis is sorted is chosen based on which axis
    // is closest to being perdendicular to the
    // camera's gaze direction
    int (*compare)(const void *, const void *) = comp_face_max_x;
    if (scale == abs(max_y - min_y)) {
        compare = comp_face_max_y;
    } else if (scale == abs(max_z - min_z)) {
        compare = comp_face_max_z;
    }

    qsort(faces, num_faces, sizeof(TriangleFace_N), compare);

    box->setChildren(faces, 0, num_faces);

    mesh_obj.close();
    free(vertices);
    free(normals);
}
double EPSILON = std::numeric_limits<double>::epsilon();
class B {
   public:
    double lambda;
    BVH_Node *node;
    B(double l, BVH_Node *n) {
        lambda = l;
        node = n;
    }
    B(double l) {
        lambda = l;
        node = NULL;
    }

    friend bool operator<(const B &lhs, const B &rhs) {
        return lhs.lambda < rhs.lambda;// < -EPSILON;
    }
    friend bool operator>(const B &lhs, const B &rhs) {
        return lhs.lambda > rhs.lambda;// > EPSILON;
    }

    friend std::ostream &operator<<(std::ostream &strm, const B &a) {
        return strm << a.lambda;
    }
};

void Mesh::intersect(struct ray *ray, double *lambda, struct point *p,
                     struct point *n, double *a, double *b) {
    *lambda = INFINITY;

    // debug :
    num_intersect_calls++;

    /*
    std::priority_queue<B, std::vector<B>, std::greater<B>> q;
    for (double n : {1., 8.8, 5., 6., 3., 4., 0., 9., 7., 2.}) {
        q.emplace(n);
    }
    while (!q.empty()) {
        std::cout << q.top() << " ";
        q.pop();
    }
    std::cout << '\n';
    */

#ifdef DEBUG
    std::cout << " Mesh call\n";
#endif

    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    double f_lambda = INFINITY;
    point bary_coords;
    TriangleFace *closest_face = NULL;

    std::priority_queue<B, std::vector<B>, std::greater<B>> bvh_queue;

    box->intersect(&ray_transformed, &f_lambda, &bary_coords);
    if (f_lambda < INFINITY) {
        bvh_queue.emplace(f_lambda, box);
    }

    double curr_lamb = INFINITY;
    point f_bary_coords;
    BVH_Node *currentNode = NULL;
    BoundingBox *currentBox = NULL;
    while (!bvh_queue.empty()) {
        curr_lamb = bvh_queue.top().lambda;
        currentNode = bvh_queue.top().node;
        bvh_queue.pop();

        if (curr_lamb < *lambda) {
            if (currentNode->isFace()) {
                    *lambda = curr_lamb;
                    closest_face = (TriangleFace *)currentNode;
            } else {
                currentBox = ((BoundingBox *)currentNode);

                f_lambda = INFINITY;
                currentBox->c1->intersect(&ray_transformed, &f_lambda, &f_bary_coords);
                if (f_lambda < INFINITY) {
                    bvh_queue.emplace(f_lambda, currentBox->c1);
                }

                if (currentBox->c2 != NULL) {
                    f_lambda = INFINITY;
                    currentBox->c2->intersect(&ray_transformed, &f_lambda, &f_bary_coords);
                    if (f_lambda < INFINITY) {
                        bvh_queue.emplace(f_lambda, currentBox->c2);
                    }
                }
            }
        }
    }

    if (bary_coords.w == -1) {
        *lambda = f_lambda;
        *a = -1;
        *p = bary_coords;
        return;
    }

    if (closest_face != NULL) {
        closest_face->intersect(&ray_transformed, &f_lambda, &bary_coords);
        //*lambda = f_lambda;
        *n = closest_face->normal(&bary_coords);
        // std::cout << typeid(closest_face).name() << "\n";
        if (num_normals > 1) {
            // std::cout << "closest face: " << *closest_face << "\n";
            *n = ((TriangleFace_N *)closest_face)->normal(&bary_coords);
        }
        normalTransform(n, n, this);
        normalize(n);
        rayPosition(ray, *lambda, p);
    } else {
        *lambda = -1;
    }
}

int comp_face_max_x(const void *a, const void *b) {
    TriangleFace_N *ta = (TriangleFace_N *)a, *tb = (TriangleFace_N *)b;
    double result = ta->max_x() - tb->max_x();
    return result > 0;
}

int comp_face_max_y(const void *a, const void *b) {
    TriangleFace_N *ta = (TriangleFace_N *)a, *tb = (TriangleFace_N *)b;
    double result = ta->max_y() - tb->max_y();
    return result > 0;
}

int comp_face_max_z(const void *a, const void *b) {
    TriangleFace_N *ta = (TriangleFace_N *)a, *tb = (TriangleFace_N *)b;
    double result = ta->max_z() - tb->max_z();
    return result > 0;
}