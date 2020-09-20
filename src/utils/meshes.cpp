#include "meshes.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <queue>

#include "ray.h"

//#define DEBUG

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

double TriangleFace::min_x() const { return MIN(p1.x, MIN(p2.x, p3.x)); }
double TriangleFace::min_y() const { return MIN(p1.y, MIN(p2.y, p3.y)); }
double TriangleFace::min_z() const { return MIN(p1.z, MIN(p2.z, p3.z)); }
double TriangleFace::max_x() const { return MAX(p1.x, MAX(p2.x, p3.x)); }
double TriangleFace::max_y() const { return MAX(p1.y, MAX(p2.y, p3.y)); }
double TriangleFace::max_z() const { return MAX(p1.z, MAX(p2.z, p3.z)); }

void TriangleFace::intersect(struct Ray *ray, double *lambda,
                             point *bary_coords) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified triangle.

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
        return;          // intersection with plane is negative

    double u, v, w;

    rayPosition(ray, t_lambda, &pi);

    // e12 has already been calculated
    point v1i = pi - p1;
    point crossp = cross(&e12, &v1i);
    u = dot(&crossp, &normal);
    if (u < 0) return;  // outside first edge

    // e23 has already been calculated
    point v2i = pi - p2;
    crossp = cross(&e23, &v2i);
    v = dot(&crossp, &normal);
    if (v < 0) return;  // outside second edge

    point e31 = p1 - p3;
    point v3i = pi - p3;
    crossp = cross(&e31, &v3i);
    w = dot(&crossp, &normal);
    if (w < 0) return;  // outside third edge

    // intersection is within triangle
    bary_coords->z = u / denom;
    bary_coords->x = v / denom;
    bary_coords->y = w / denom;
    *lambda = t_lambda;
    return;
}

double TriangleFace::intersect(struct Ray *ray, double lambda) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified triangle.

    point e12 = p2 - p1;
    point e23 = p3 - p2;

    point normal = cross(&e12, &e23);

    // ray plane intersection calculation
    point pi, r = p1 - ray->p0;
    double d_dot_n = dot(&ray->d, &normal);
    double r_dot_n = dot(&r, &normal);
    double t_lambda = r_dot_n / d_dot_n;

    if (t_lambda < THR || lambda < t_lambda)
        return INFINITY;  // intersection with plane is negative
                          // or larger than current lambda

    double u, v, w;

    rayPosition(ray, t_lambda, &pi);

    // e12 has already been calculated
    point v1i = pi - p1;
    point crossp = cross(&e12, &v1i);
    u = dot(&crossp, &normal);
    if (u < 0) return INFINITY;  // outside first edge

    // e23 has already been calculated
    point v2i = pi - p2;
    crossp = cross(&e23, &v2i);
    v = dot(&crossp, &normal);
    if (v < 0) return INFINITY;  // outside second edge

    point e31 = p1 - p3;
    point v3i = pi - p3;
    crossp = cross(&e31, &v3i);
    w = dot(&crossp, &normal);
    if (w < 0) return INFINITY;  // outside third edge

    return t_lambda;
}

bool TriangleFace::isprim() { return true; }

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

void Mesh::setMesh(const char *filename) {
    frontAndBack = 0;

    std::string line;
    std::ifstream mesh_obj(filename);

    std::streampos first_vertex, first_normal, first_face;

    num_vertices = 0, num_faces = 0, num_normals = 0;

    int v = 1, vn = 1, f = 0;
    double x, y, z, scale;
    int p1, p2, p3, n1, n2, n3;
    double min_x, max_x, avg_x, min_y, max_y, avg_y, min_z, max_z, avg_z;
    min_x = min_y = min_z = INFINITY;
    max_x = max_y = max_z = -INFINITY;
    avg_x = avg_y = avg_z = 0;

    if (!mesh_obj.is_open()) {
        std::cout << "Unable to open mesh file: " << filename << "\n";
        return;
    }

    while (getline(mesh_obj, line)) {
        if (line.rfind("v ", 0) == 0) {
            //count vertices and update average location
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
        } else if (line.rfind("vn ", 0) == 0) {
            //count vertex normals
            //std::cout << line << "\n";
            num_normals++;
        } else if (line.rfind("f ", 0) == 0) {
            //count faces
            num_faces++;
        }
    }

    //go to begining of file
    mesh_obj.clear();
    mesh_obj.seekg(0);

    // calculate average distance from 0 so that mesh can be centered
    avg_x /= num_vertices, avg_y /= num_vertices, avg_z /= num_vertices;
    // calculate max dimension so mesh can be scaled down to 1x1x1 ish
    scale = MAX(abs(max_x - min_x), MAX(abs(max_y - min_y), abs(max_z - min_z)));

    num_vertices += 1;  // << overcount by 1
    vertices = (point *)malloc(num_vertices * sizeof(point));

    num_normals += 1;
    if (num_normals > 1) {
        normals = (point *)malloc(num_normals * sizeof(point));
    }

    faces = (PrimitiveData *)malloc(num_faces * sizeof(PrimitiveData));
    while (getline(mesh_obj, line)) {
        if (line.rfind("v ", 0) == 0) {
            sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
            x = (x - avg_x) / scale, y = (y - avg_y) / scale,
            z = -(z - avg_z) / scale;
            vertices[v] = point(x, y, z);
            v++;
        } else if(line.rfind("vn ", 0) == 0){
            sscanf(line.c_str(), "vn %lf %lf %lf", &x, &y, &z);
            normals[vn] = point(x, y, -z);
            vn++;
        } else if (line.rfind("f ", 0) == 0){
            if (sscanf(line.c_str(), "f %d %d %d", &p1, &p2, &p3) == 3){
                faces[f] = new TriangleFace(vertices[p3], vertices[p2], vertices[p1]);
                f++;
            } else if(sscanf(line.c_str(), "f %d//%d %d//%d %d//%d", &p1, &n1, &p2, &n2, &p3, &n3) == 6){
                faces[f] = new TriangleFace_N(vertices[p1], vertices[p2], vertices[p3]);
                ((TriangleFace_N *)faces[f].mprimitive)->setNormals(normals[n1], normals[n2], normals[n3]);
                f++;
            }
        }
    }
    bvh.set_build_method(BuildMethod::MidSplit);
    bvh.set_search_method(SearchMethod::BFS);
    bvh.build(faces, num_faces);

    mesh_obj.close();

    free(vertices);
    free(normals);
    free(faces);
}

void Mesh::set_canonical_bounds() {
    w_bound.min.x = bvh.root->min_x();
    w_bound.min.y = bvh.root->min_y();
    w_bound.min.z = bvh.root->min_z();
    w_bound.max.x = bvh.root->max_x();
    w_bound.max.y = bvh.root->max_y();
    w_bound.max.z = bvh.root->max_z();
}

void Mesh::intersect(struct Ray *ray, double *lambda, struct point *p,
                     struct point *n, double *a, double *b) {
    *lambda = INFINITY;

    struct Ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    TriangleFace *closest_face;
    //((bvh).*(bvh.search))(ray);
    closest_face = (TriangleFace *)bvh.search(&ray_transformed);

    point bary_coords;
    if (closest_face != NULL) {
        closest_face->intersect(&ray_transformed, lambda, &bary_coords);
        *n = closest_face->normal(&bary_coords);
        normalTransform(n, n, this);
        normalize(n);
        rayPosition(ray, *lambda, p);
    }
}
