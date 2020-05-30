#include "objects.h"

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <string>

#include "mappings.h"
#include "ray.h"
#include "utils.h"

Object::Object(double r = 1, double g = 1, double b = 1) {
    col.R = r;
    col.G = g;
    col.B = b;
    rt.alpha = 1;
    rt.shinyness = 2;
    pt.LSweight = 1;
    r_index = 1;
    refl_sig = 0;
    texImg = NULL;
    normalMap = NULL;
    alphaMap = NULL;
    frontAndBack = 0;
    isLightSource = 0;
    T = I();
    next = NULL;
}

void Object::set_rayTrace_properties(double ambient, double diffuse,
                                     double specular, double global,
                                     double alpha, double shiny) {
    rt.ambient = ambient;
    rt.diffuse = diffuse;
    rt.specular = specular;
    rt.global = global;
    rt.alpha = alpha;
    rt.shinyness = shiny;

    // create a similar set of params incase the object is drawn in pathtrace
    // mode
    pt.diffuse = alpha * (ambient + diffuse);
    pt.reflect = alpha * (1 - diffuse);
    pt.refract = 1 - alpha;

    // ensure sum to 1
    double sum = pt.diffuse + pt.reflect + pt.refract;
    pt.diffuse /= sum;
    pt.reflect /= sum;
    pt.refract /= sum;
}

void Object::set_pathTrace_properties(double diffuse, double reflect,
                                      double refract) {
    pt.diffuse = diffuse;
    pt.reflect = reflect;
    pt.refract = refract;

    // ensure sum to 1
    double sum = pt.diffuse + pt.reflect + pt.refract;
    pt.diffuse /= sum;
    pt.reflect /= sum;
    pt.refract /= sum;

    // create a similar set of parameters for ray tracing
    rt.ambient = 0.1 * diffuse;
    rt.diffuse = 0.9 * diffuse;
    rt.specular = reflect;
    rt.global = reflect;
    rt.alpha = 1 - refract;
    rt.shinyness = reflect * 10;
}

void Object::surfaceCoordinates(double a, double b, double *x, double *y,
                                double *z) {
    fprintf(stderr, "Object::surfaceCoordinates\n");
}

void Object::randomPoint(double *x, double *y, double *z) {
    point r = 0;
    r = T * r;
    *x = r.x;
    *y = r.y;
    *z = r.z;
    fprintf(stderr, "Object::randomPoint\n");
}

///////////////////////////////Plane//////////////////////////////////////

Plane::Plane(double r = 1, double g = 1, double b = 1) : Object(r, g, b) {
    frontAndBack = 1;
}

void Plane::intersect(struct ray *ray, double *lambda, struct point *p,
                      struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical plane.

    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);
    *lambda = -1;
    struct point norm;
    // normal of canonical plane
    norm.x = 0;
    norm.y = 0;
    norm.z = 1;
    struct point p1;

    double l;
    double d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0) {
        p1.x = -ray_transformed.p0.x;
        p1.y = -ray_transformed.p0.y;
        p1.z = -ray_transformed.p0.z;

        l = dot(&(p1), &norm) / d_dot_n;
        // Check if the intersection point is inside the plane
        rayPosition(&ray_transformed, l, p);
        double x = p->x;
        double y = p->y;
        if (fabs(p->z) < THR && fabs(p->x) <= 1 + THR &&
            fabs(p->y) <= 1 + THR) {
            *lambda = l;
            rayPosition(ray, l, p);
            // printf("n: (%f, %f, %f)\n", n->x, n->y, n->z);
            normalTransform(&norm, n, this);
            // printf("nt: (%f, %f, %f)\n", n->x, n->y, n->z);
        }

        *a = (x + 1) / 2;
        *b = (-y + 1) / 2;
    }
}

void Plane::surfaceCoordinates(double a, double b, double *x, double *y,
                               double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2
    // parameters a,b in [0,1]. 'a' controls displacement from the left side of
    // the plane, 'b' controls displacement from the bottom of the plane.

    struct point p;
    p.x = -2 * a + 1;
    p.y = -2 * b + 1;
    p.z = 0;
    p.w = 1;

    p = T * p;
    // matVecMult(plane->T, &p);

    *x = p.x;
    *y = p.y;
    *z = p.z;
}

void Plane::randomPoint(double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the
    // plane Sapling should be uniform, meaning there should be an equal chance
    // of getting any spot on the plane

    double a = drand48();
    double b = drand48();
    surfaceCoordinates(a, b, x, y, z);
}

void Sphere::intersect(struct ray *ray, double *lambda, struct point *p,
                       struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical sphere.

    struct ray ray_transformed;
    double l1 = -1, l2 = -1;
    rayTransform(ray, &ray_transformed, this);
    solveQuadratic(&ray_transformed, &l1, &l2);

    *lambda = -1;
    if (MAX(l1, l2) >= THR) {
        if (MIN(l1, l2) <= THR) {
            *lambda = MAX(l1, l2);
        } else {
            *lambda = MIN(l1, l2);
        }
    }

    double x, y, z;
    if (*lambda != -1) {
        rayPosition(&ray_transformed, *lambda, p);
        n->x = p->x;
        n->y = p->y;
        n->z = p->z;
        x = p->x;
        y = p->y;
        z = p->z;
        n->w = 1;

        normalize(n);

        normalTransform(n, n, this);

        // use the lambda to set the position along the untransformed ray
        rayPosition(ray, *lambda, p);
        *a = 0.5 + (atan2(z, x)) / (2 * PI);
        *b = 0.5 - (asin(y)) / (PI);
    }
}

void Sphere::surfaceCoordinates(double a, double b, double *x, double *y,
                                double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2
    // parameters a,b in [0,1]. 'a' controls displacement from the left side of
    // the plane, 'b' controls displacement from the bottom of the plane.

    struct point p;
    p.x = cos(a) * sin(b);
    p.y = sin(a) * sin(b);
    p.z = cos(b);
    p.w = 1;

    // matVecMult(sphere->T, &p);
    p = T * p;

    *x = p.x;
    *y = p.y;
    *z = p.z;
}

void Sphere::randomPoint(double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the
    // plane Sapling should be uniform, meaning there should be an equal change
    // of gedtting any spot on the plane

    double a = drand48() * 2 * PI;
    double b = drand48() * 2 - 1;
    b = acos(b);
    surfaceCoordinates(a, b, x, y, z);
}

void Box::intersect(struct ray *ray, double *lambda, struct point *p,
                    struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical box.

    // default large number to compare intersection lamdas to
    double MAX = 1000000000;

    *lambda = MAX;

    // current intersection lambda
    double b_lambda;

    // use the objects inverse transform to get a ray in the cannonical space
    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    // y-z plane box face at x = -0.5
    b_lambda = (-0.5 - ray_transformed.p0.x) / ray_transformed.d.x;
    if (THR < b_lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->y && p->y < 0.5) && (-0.5 < p->z && p->z < 0.5)) {
            // hit the face at x = -0.5
            *lambda = b_lambda;
            *n = point(-1, 0, 0);
        }
    }

    // y-z plane box face at x = 0.5
    b_lambda = (0.5 - ray_transformed.p0.x) / ray_transformed.d.x;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->y && p->y < 0.5) && (-0.5 < p->z && p->z < 0.5)) {
            // hit the face at x = -0.5
            *lambda = b_lambda;
            *n = point(1, 0, 0);
        }
    }

    // x-z plane box face at y = -0.5
    b_lambda = (-0.5 - ray_transformed.p0.y) / ray_transformed.d.y;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->x && p->x < 0.5) && (-0.5 < p->z && p->z < 0.5)) {
            // hit the face at y = -0.5
            *lambda = b_lambda;
            *n = point(0, -1, 0);
        }
    }

    // x-z plane box face at y = 0.5
    b_lambda = (0.5 - ray_transformed.p0.y) / ray_transformed.d.y;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->x && p->x < 0.5) && (-0.5 < p->z && p->z < 0.5)) {
            // hit the face at y = -0.5
            *lambda = b_lambda;
            *n = point(0, 1, 0);
        }
    }

    // x-y plane box face at z = -0.5
    b_lambda = (-0.5 - ray_transformed.p0.z) / ray_transformed.d.z;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->x && p->x < 0.5) && (-0.5 < p->y && p->y < 0.5)) {
            // hit the face at z = -0.5
            *lambda = b_lambda;
            *n = point(0, 0, -1);
        }
    }

    // x-y plane box face at z = 0.5
    b_lambda = (0.5 - ray_transformed.p0.z) / ray_transformed.d.z;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->x && p->x < 0.5) && (-0.5 < p->y && p->y < 0.5)) {
            // hit the face at z = -0.5
            *lambda = b_lambda;
            *n = point(0, 0, 1);
        }
    }

    // printf("intersect point: (%f, %f, %f)\n", p->x, p->y, p->z);

    // if there is an intersection, update the normal and intersection point
    if (*lambda != MAX) {
        normalTransform(n, n, this);
        rayPosition(ray, *lambda, p);
    } else {  // return no intersection flag
        *lambda = -1;
    }

    /*
    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);
    *lambda = -1;
    struct point norm;
    struct point temp_p;
    struct point p1;
    struct point e;
    double x;
    double y;

    double l = -1;
    double d_dot_n;
    e.x = ray_transformed.p0.x;
    e.y = ray_transformed.p0.y;
    e.z = ray_transformed.p0.z;

    norm.x = 1;
    norm.y = 0;
    norm.z = 0;
    norm.w = 1;
    d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0)
    {
      p1.x = 0.5 - e.x;
      p1.y = -e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l > THR)
      {
        // Check if the intersection point is inside the plane
        rayPosition(&ray_transformed, l, &temp_p);

        if (fabs(temp_p.x - 0.5) < THR && fabs(temp_p.z) < 0.5 + THR &&
    fabs(temp_p.y) < 0.5 + THR && (*lambda == -1 || l < *lambda))
        {
          rayPosition(ray, l, p);
          *lambda = l;
          normalTransform(&norm, n, this);
        }
      }
    }

    l = -1;
    norm.x = -1;
    norm.y = 0;
    norm.z = 0;
    norm.w = 1;
    d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0)
    {
      p1.x = -0.5 - e.x;
      p1.y = -e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l > THR)
      {
        // Check if the intersection point is inside the plane

        rayPosition(&ray_transformed, l, &temp_p);

        if (fabs(temp_p.x + 0.5) < THR && fabs(temp_p.z) < 0.5 + THR &&
    fabs(temp_p.y) < 0.5 + THR && (*lambda == -1 || l < *lambda))
        {
          rayPosition(ray, l, p);
          *lambda = l;
          normalTransform(&norm, n, this);
        }
      }
    }

    l = -1;
    norm.x = 0;
    norm.y = 1;
    norm.z = 0;
    norm.w = 1;

    d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0)
    {
      p1.x = -e.x;
      p1.y = 0.5 - e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l >= 0)
      {
        // Check if the intersection point is inside the plane
        rayPosition(&ray_transformed, l, &temp_p);

        if (fabs(temp_p.z) < 0.5 + THR && fabs(temp_p.x) < 0.5 + THR &&
    fabs(temp_p.y - 0.5) < THR && (*lambda == -1 || l < *lambda))
        {
          rayPosition(ray, l, p);

          *lambda = l;
          normalTransform(&norm, n, this);
          x = p->x;
          y = p->y;
        }
      }
    }

    l = -1;
    norm.x = 0;
    norm.y = -1;
    norm.z = 0;
    norm.w = 1;

    d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0)
    {
      p1.x = -e.x;
      p1.y = -0.5 - e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
        // Check if the intersection point is inside the plane
        rayPosition(&ray_transformed, l, &temp_p);

        if (fabs(temp_p.z) < 0.5 + THR && fabs(temp_p.x) < 0.5 + THR &&
    fabs(temp_p.y + 0.5) < THR && (*lambda == -1 || l < *lambda))
        {
          rayPosition(ray, l, p);
          *lambda = l;
          normalTransform(&norm, n, this);
        }
      }
    }

    l = -1;
    norm.x = 0;
    norm.y = 0;
    norm.z = 1;
    norm.w = 1;

    d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0)
    {
      p1.x = -e.x;
      p1.y = -e.y;
      p1.z = 0.5 - e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
        // Check if the intersection point is inside the plane
        rayPosition(&ray_transformed, l, &temp_p);

        if (fabs(temp_p.z - 0.5) < THR && fabs(temp_p.x) < 0.5 + THR &&
    fabs(temp_p.y) < 0.5 + THR && (*lambda == -1 || l < *lambda))
        {
          rayPosition(ray, l, p);
          *lambda = l;
          normalTransform(&norm, n, this);
        }
      }
    }

    l = -1;
    norm.x = 0;
    norm.y = 0;
    norm.z = -1;
    norm.w = 1;

    d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0)
    {
      p1.x = -e.x;
      p1.y = -e.y;
      p1.z = -0.5 - e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
        // Check if the intersection point is inside the plane
        rayPosition(&ray_transformed, l, &temp_p);

        if (fabs(temp_p.z + 0.5) < THR && fabs(temp_p.x) < 0.5 + THR &&
    fabs(temp_p.y) < 0.5 + THR && (*lambda == -1 || l < *lambda))
        {
          rayPosition(ray, l, p);
          *lambda = l;

          normalTransform(&norm, n, this);
        }
      }
    }

    if (*lambda > -1)
    {
      x = p->x;
      y = p->y;
      if (this->texImg != NULL)
      {
        *a = (x + 1) / 2;
        *b = (-y + 1) / 2;
      }
    }
    */
}

Triangle::Triangle(double r = 1, double g = 1, double b = 1) : Object(r, g, b) {
    frontAndBack = 1;
    // equilateral triangle on x-y plane unit circle centered at 0
    p1.x = -0.866, p1.y = -0.5;
    p2.x = 0.866, p2.y = -0.5;
    p3.x = 0.0, p3.y = 1.0;
}

void Triangle::intersect(struct ray *ray, double *lambda, struct point *p,
                         struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified triangle.

    // default no intersection
    *lambda = -1;

    // use the objects inverse transform to get a ray in the cannonical space
    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    point e12 = p2 - p1;
    point e23 = p3 - p2;

    point normal = cross(&e12, &e23);
    // std::cout << "normal: " << normal << std::endl;

    // ray plane intersection calculation
    point r = p1 - ray_transformed.p0;
    double d_dot_n = dot(&ray_transformed.d, &normal);
    double r_dot_n = dot(&r, &normal);
    double t_lambda = r_dot_n / d_dot_n;

    if (THR < t_lambda) {  // check that intersection with plane is positive
        rayPosition(&ray_transformed, t_lambda, p);
        // e12 has already been calculated
        point v1i = *p - p1;

        point crossp = cross(&e12, &v1i);
        if (dot(&crossp, &normal) >= 0) {  // check within first edge
            // e23 has already been calculated
            point v2i = *p - p2;
            crossp = cross(&e23, &v2i);
            if (dot(&crossp, &normal) >= 0) {  // check within second edge
                point e31 = p1 - p3;
                point v3i = *p - p3;
                crossp = cross(&e31, &v3i);
                if (dot(&crossp, &normal) >= 0) {  // check within third edge
                    // intersection is within triangle
                    *lambda = t_lambda;
                    normalTransform(&normal, n, this);
                    normalize(n);
                    rayPosition(ray, *lambda, p);
                }
            }
        }
    }
}

void Triangle::setPoints(double x1, double y1, double z1, double x2, double y2,
                         double z2, double x3, double y3, double z3) {
    p1.x = x1, p1.y = y1, p1.z = z1;
    p2.x = x2, p2.y = y2, p2.z = z2;
    p3.x = x3, p3.y = y3, p3.z = z3;
}

void Triangle::setPoints(point p1, point p2, point p3) {
    this->p1 = p1, this->p2 = p2, this->p3 = p3;
}

Polygon::Polygon(double r = 1, double g = 1, double b = 1) : Object(r, g, b) {
    frontAndBack = 1;
    p = (point *)malloc(3 * sizeof(point));
    e = (point *)malloc(3 * sizeof(point));
    // equilateral triangle on x-y plane unit circle centered at 0
    p[0].x = -0.866, p[0].y = -0.5;
    p[1].x = 0.866, p[1].y = -0.5;
    p[2].x = 0.0, p[2].y = 1.0;
    calculate_edge_vectors();
}

void Polygon::intersect(struct ray *ray, double *lambda, struct point *pi,
                        struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified polygon.

    // default no intersection
    *lambda = -1;

    // use the objects inverse transform to get a ray in the cannonical space
    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    point e12 = p[1] - p[0];
    point e23 = p[2] - p[1];

    point normal = cross(&e[0], &e[1]);
    // std::cout << "normal: " << normal << std::endl;

    // ray plane intersection calculation
    point r = p[0] - ray_transformed.p0;
    double d_dot_n = dot(&ray_transformed.d, &normal);
    double r_dot_n = dot(&r, &normal);
    double p_lambda = r_dot_n / d_dot_n;

    if (THR < p_lambda) {  // check that intersection with plane is positive
        rayPosition(&ray_transformed, p_lambda, pi);
        point v, crossp;
        int in_edges = 0;
        for (int i = 0; i < numPoints; i++) {
            v = *pi - p[i];
            crossp = cross(&e[i], &v);
            if (dot(&crossp, &normal) < 0) break;
            in_edges++;
        }
        if (in_edges == numPoints) {
            // intersection is within polygon
            *lambda = p_lambda;
            normalTransform(&normal, n, this);
            normalize(n);
            rayPosition(ray, *lambda, pi);
        }
    }
}

void Polygon::setNumPoints(int num) {
    free(p);
    free(e);
    numPoints = num;
    currPoint = 0;
    p = (point *)malloc(numPoints * sizeof(point));
    e = (point *)malloc(numPoints * sizeof(point));
}

void Polygon::addPoint(point point) {
    p[currPoint] = point;
    currPoint++;
    if (currPoint == numPoints) {
        calculate_edge_vectors();
    }
}

void Polygon::addPoint(double x, double y, double z) {
    p[currPoint] = point(x, y, z);
    currPoint++;
    if (currPoint == numPoints) {
        calculate_edge_vectors();
    }
}

void Polygon::calculate_edge_vectors() {
    // calculate vectors for all edges
    for (int i = 0; i < numPoints - 1; i++) {
        e[i] = p[i + 1] - p[i];
    }
    // last edge goes from last point to first point
    e[numPoints - 1] = p[0] - p[numPoints - 1];
}

void Mesh::setMesh(const char *filename) {
    frontAndBack = 0;

    std::string line;
    std::ifstream mesh_obj(filename);

    std::streampos first_normal, first_face;

    double x, y, z, scale;
    double min_x, max_x, avg_x, min_y, max_y, avg_y, min_z, max_z, avg_z;
    min_x = min_y = min_z = INFINITY;
    max_x = max_y = max_z = -INFINITY;
    avg_x = avg_y = avg_z = 0;

    if (mesh_obj.is_open()) {
        
        // count the number of vertices
        num_vertices = 0;
        while (getline(mesh_obj, line) && line.rfind("v ", 0) == 0) {
            num_vertices++;
            sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
            avg_x += x, avg_y += y, avg_z +=z;
            if(x < min_x) {
                min_x = x;
            } if (max_x < x){
                max_x = x;
            }

            if(y < min_y) {
                min_y = y;
            } if (max_y < y){
                max_y = y;
            }

            if(z < min_z) {
                min_z = z;
            } if (max_z < z){
                max_z = z;
            }
        }

        //calculate average diatance from 0 so that mesh can be centered
        avg_x /= num_vertices, avg_y /= num_vertices, avg_z /= num_vertices;
        //calculate max dimension so mesh can be scaled down to 1x1x1 ish
        scale = MAX(abs(max_x - min_x), MAX(abs(max_y - min_y),abs(max_z - min_z)));
        //std::cout << scale << "\n";

        num_vertices += 1;  // << overcount by 1
        vertices = (point *)malloc(num_vertices * sizeof(point));

        // read vertices into array
        mesh_obj.seekg(0);
        
        // keeps first vertex empty so face's lookup doesn't need to subtract 1
        for (int i = 1; i < num_vertices; i++) {
            getline(mesh_obj, line);
            sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
            vertices[i] = point((x - avg_x)/scale, (y - avg_y)/scale, (z - avg_z)/scale);
        }

        //std::cout << "min: "<< min_x << " " << min_y << " " << min_z << "\n";
        //std::cout << "max: "<< max_x << " " << max_y << " " << max_z << "\n";
        //std::cout << "avg: "<< avg_x << " " << avg_y << " " << avg_z << "\n";

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
        first_face = (int)mesh_obj.tellg() ;

        // count number of faces
        num_faces = 0;  // first face has already been scanned
        while (getline(mesh_obj, line) && line.rfind("f ", 0) == 0) {
            num_faces++;
        }

        if (num_normals > 1) {
            faces = new TriangleFace_N[num_faces];
        } else {
            // faces = new TriangleFace[num_faces];
        }

        mesh_obj.clear();
        mesh_obj.seekg(first_face);
        int p1, p2, p3;
        //std::cout << num_faces;
        for (int i = 0; i < num_faces; i++) {
            getline(mesh_obj, line);
            //std::cout << line << "\n";
            sscanf(line.c_str(), "f %d %d %d", &p1, &p2, &p3);
            faces[i] = TriangleFace_N(vertices[p1], vertices[p2], vertices[p3]);
            if (num_normals > 1) {
                ((TriangleFace_N *)faces)[i].setNormals(
                    normals[p1], normals[p2], normals[p3]);
            }
        }

        num_vertices--;
        mesh_obj.close();
        free(vertices);
        free(normals);
    } else {
        std::cout << "Unable to open mesh file";
    }
}

std::ostream &operator<<(std::ostream &strm, const TriangleFace &a) {
    return strm << "p1: " << a.p1 << " p2: " << a.p2 << " p3: " << a.p3;
}

std::ostream &operator<<(std::ostream &strm, const TriangleFace_N &a) {
    return strm << "p1: " << a.p1 << " p2: " << a.p2 << " p3: " << a.p3
                << "n1: " << a.n1 << " n2: " << a.n2 << " n3: " << a.n3;
}

void Mesh::intersect(struct ray *ray, double *lambda, struct point *p,
                     struct point *n, double *a, double *b) {
    *lambda = INFINITY;

    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    double f_lambda;
    TriangleFace *closest_face = NULL;
    point f_coords;
    for (int i = 0; i < num_faces; i++) {
        f_lambda = INFINITY;
        faces[i].intersect(&ray_transformed, &f_lambda, &f_coords);
        if (f_lambda < *lambda) {
            *lambda = f_lambda;
            closest_face = &faces[i];
            bary_coords = f_coords;
        }
    }

    if (*lambda < INFINITY) {
        point pi;
        rayPosition(&ray_transformed, *lambda, &pi);
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

void TriangleFace::intersect(struct ray *ray, double *lambda,
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

    if (t_lambda < THR) return;  // intersection with plane is negative

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
}

point TriangleFace::normal(point *bary_coords) {
    point e12 = p2 - p1;
    point e23 = p3 - p2;

    return cross(&e12, &e23);
}

void TriangleFace_N::setNormals(point n1, point n2, point n3) {
    this->n1 = n1, this->n2 = n2, this->n3 = n3;
}

point TriangleFace_N::normal(point *bary_coords) {
    // std::cout << "n1: " << n1 << " n2: " << n2 << " n3: " << n3 << "\n";
    return n1 * bary_coords->x + n2 * bary_coords->y + n3 * bary_coords->z;
}

void insertObject(Object *o, Object **list) {
    if (o == NULL) return;
    // Inserts an object into the object list.
    if (*(list) == NULL) {
        *(list) = o;
        (*(list))->next = NULL;
    } else {
        o->next = (*(list))->next;
        (*(list))->next = o;
    }
}

inline void normalTransform(struct point *n, struct point *n_transformed,
                            Object *obj) {
    // Computes the normal at an affinely transformed point given the original
    // normal and the object's inverse transformation. From the notes:
    // n_transformed=A^-T*n normalized.
    double x = obj->Tinv.T[0][0] * n->x + obj->Tinv.T[1][0] * n->y +
               obj->Tinv.T[2][0] * n->z;
    double y = obj->Tinv.T[0][1] * n->x + obj->Tinv.T[1][1] * n->y +
               obj->Tinv.T[2][1] * n->z;
    double z = obj->Tinv.T[0][2] * n->x + obj->Tinv.T[1][2] * n->y +
               obj->Tinv.T[2][2] * n->z;
    n_transformed->x = x;
    n_transformed->y = y;
    n_transformed->z = z;
    n_transformed->w = 1;
    normalize(n_transformed);
}

struct pointLS *newPLS(struct point *p0, double r, double g, double b) {
    // Allocate a new point light sourse structure. Initialize the light
    // source to the specified RGB colour
    // Note that this is a point light source in that it is a single point
    // in space, if you also want a uniform direction for light over the
    // scene (a so-called directional light) you need to place the
    // light source really far away.

    struct pointLS *ls = (struct pointLS *)calloc(1, sizeof(struct pointLS));
    if (!ls)
        fprintf(stderr, "Out of memory allocating light source!\n");
    else {
        memcpy(&ls->p0, p0,
               sizeof(struct point));  // Copy light source location

        ls->col.R = r;  // Store light source colour and
        ls->col.G = g;  // intensity
        ls->col.B = b;
    }
    return (ls);
}

void insertPLS(struct pointLS *l, struct pointLS **list) {
    if (l == NULL) return;
    // Inserts a light source into the list of light sources
    if (*(list) == NULL) {
        *(list) = l;
        (*(list))->next = NULL;
    } else {
        l->next = (*(list))->next;
        (*(list))->next = l;
    }
}
