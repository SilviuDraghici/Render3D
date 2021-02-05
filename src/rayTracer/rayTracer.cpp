#include "rayTracer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../utils/affineTransforms.h"
#include "../utils/buildscene.h"
#include "../utils/camera.h"
#include "../utils/color.h"
#include "../utils/imageProcessor.h"
#include "../utils/mappings.h"
#include "../utils/objects.h"
#include "../utils/ray.h"
#include "../utils/timer.h"
#include "../utils/utils.h"

//#define DEBUG

static Scene *scene;
void rayTraceMain(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "RayTracer: Can not parse input parameters\n");
        fprintf(stderr, "USAGE: Render3D 0 size antialias output_name\n");
        fprintf(stderr, "   size = Image size (both along x and y)\n");
        fprintf(stderr,
                "   antialias = A single digit, 0 disables antialiasing. "
                "Anything else enables antialiasing\n");
        fprintf(
            stderr,
            "   output_name = Name of the output file, e.g. MyRender.ppm\n");
        exit(0);
    }

    Scene sc;
    scene = &sc;

    scene->sx = atoi(argv[2]);
    strcpy(&output_name[0], argv[4]);

    if (6 <= argc) {
        sc.frame = atoi(argv[5]) - 1;
    }

    unsigned char *rgbIm;
    // Allocate memory for the new image
    outImage = newImage(scene->sx, scene->sx, sizeof(unsigned char));
    if (!outImage) {
        fprintf(stderr, "Unable to allocate memory for raytraced image\n");
        exit(0);
    } else {
        rgbIm = (unsigned char *)outImage->rgbdata;
    }

    // Camera center is at (0,0,-1)
    scene->cam_pos = point(0, 0, -1);

    // To define the gaze vector, we choose a point 'pc' in the scene that
    // the camera is looking at, and do the vector subtraction pc-e.
    // Here we set up the camera to be looking at the origin.

    scene->cam_gaze_point = point(0, 0, 0);
    scene->cam_gaze = scene->cam_gaze_point - scene->cam_pos;

    scene->cam_up = point(0, 1, 0);

    scene->cam_focal = -1;

    Timer buildscene_timer("Buildscene");
    buildscene_timer.start();
    buildScene(scene);
    buildscene_timer.end();
    buildscene_timer.print_elapsed_time(std::cerr);

    view *cam;  // Camera and view for this scene
    cam = setupView(&(scene->cam_pos), &(scene->cam_gaze), &(scene->cam_up),
                    scene->cam_focal, -2, 2, 4);

    if (cam == NULL) {
        fprintf(stderr,
                "Unable to set up the view and camera parameters. Our of "
                "memory!\n");
        // cleanup(object_list, light_list, texture_list);
        deleteImage(outImage);
        exit(0);
    }

    setPixelStep(scene, cam, scene->sx, scene->sx);

    Ray ray;
    double depth;
    color col;
    color background = 0;

    std::cout << "Rendering...\n";

    //debug :
    num_intersection_tests = 0;
    num_bvh_searches = 0;
    Timer raytracing_timer("Raytracing");
    raytracing_timer.start();

    double start_j = 0, end_j = scene->sx, start_i = 0, end_i = scene->sx;
#ifdef DEBUG
    start_i = 831, start_j = 815;
    end_i = start_i + 1, end_j = start_j + 1;
#endif

    for (int j = start_j; j < end_j; j++) {  // For each of the pixels in the image
        for (int i = start_i; i < end_i; i++) {
#ifdef DEBUG
            printf("-------pixel: (%d %d)-------\n", i, j);
#endif
            getRayFromPixel(scene, &ray, cam, i, j);
            depth = 1;
            col = 0;
            rayTrace(&ray, depth, &col, NULL);
            if (col.R == -1) {
                col = background;
            }
            rgbIm[3 * (j * outImage->sx + i) + 0] =
                (unsigned char)(255 * col.R);
            rgbIm[3 * (j * outImage->sx + i) + 1] =
                (unsigned char)(255 * col.G);
            rgbIm[3 * (j * outImage->sx + i) + 2] =
                (unsigned char)(255 * col.B);
        }  // end for i
    }      // end for j
    raytracing_timer.end();
    // Output rendered image
    imageOutput(outImage, output_name);

    std::cout << "Ray Tracing Done!\n";

    free(cam);
    
    //debug :
    raytracing_timer.print_elapsed_time(std::cerr);
    std::cerr << "Average number of intersection tests: " << num_intersection_tests / num_bvh_searches << "\n";
}

void rayTrace(Ray *ray, int depth, color *col, Object *Os) {
    // Trace one ray through the scene.
    //
    // Parameters:
    //   *ray   -  A pointer to the ray being traced
    //   depth  -  Current recursion depth for recursive raytracing
    //   *col   - Pointer to an RGB colour structure so you can return the
    //   object colour
    //            at the intersection point of this ray with the closest scene
    //            object.
    //   *Os    - 'Object source' is a pointer to the object from which the ray
    //            originates so you can discard self-intersections due to
    //            numerical errors. NULL for rays originating from the center of
    //            projection.

    double lambda;       // Lambda at intersection
    double a = 0, b;     // Texture coordinates
    Object *obj = NULL;  // Pointer to object at intersection
    point p;      // Intersection point
    point n;      // Normal at intersection
    color I;      // Colour returned by shading function

    if (depth > scene->rt_max_depth)  // Max recursion depth reached. Return
                                      // invalid colour.
    {
        *col = -1;
        return;
    }
    findFirstHit(scene, ray, &lambda, Os, &obj, &p, &n, &a, &b);

    if (lambda == -1 || obj == NULL) {
        *col = -1;
        return;
    }
    rtShade(obj, &p, &n, ray, depth, a, b, col);
}

void rtShade(Object *obj, point *p, point *n, Ray *ray,
             int depth, double a, double b, color *col) {
    // This function implements the shading model as described in lecture. It
    // takes
    // - A pointer to the first object intersected by the ray (to get the colour
    // properties)
    // - The coordinates of the intersection point (in world coordinates)
    // - The normal at the point
    // - The ray (needed to determine the reflection direction to use for the
    // global component, as well as for
    //   the Phong specular component)
    // - The current recursion depth
    // - The (a,b) texture coordinates (meaningless unless texture is enabled)
    //
    // Returns:
    // - The colour for this ray (using the col pointer)
    //

    color tmp_col;  // Accumulator for colour components
    color objcol;   // Colour for the object in R G and B
    double alpha;

    // This will hold the colour as we process all the components of
    // the Phong illumination model
    tmp_col = 0;
    *col = 0;

    // set object color
    textureMap(obj, a, b, &objcol);

    // check for normal map
    normalMap(obj, a, b, n);

    // check for alpha map
    alphaMap(obj, a, b, &alpha, obj->rt.alpha);

    // vector from intersection point to camera
    point c;
    c = -ray->d;
    normalize(&c);

    double ra = obj->rt.ambient, rd = obj->rt.diffuse, rs = obj->rt.specular,
           rg = obj->rt.global;

    color ambient;
    ambient = objcol * ra * alpha;

    color diffuse;

    color specular;

    color refract = 0;

    color global;

    PointLS *light = scene->rt_point_light_list;

    while (light != NULL) {
        // ray from intersection point to light source
        Ray pToLight;
        // pToLight.p0 = *p - n * THR;
        pToLight.p0 = *p;

        pToLight.d = light->p0 - *p;
        double light_lambda;
        // we don't care about these details of the hit
        Object *dummyobj = NULL;
        point dummyp;
        double dummya;

        findFirstHit(scene, &pToLight, &light_lambda, obj, &dummyobj, &dummyp,
                     &dummyp, &dummya, &dummya);

        if (!(THR < light_lambda && light_lambda < 1)) {
            // vector from intersection point to light
            point s;
            s = light->p0 - *p;
            normalize(&s);

            double n_dot_s = MAX(0, dot(n, &s));
            if (obj->frontAndBack) {
                n_dot_s = MAX(0, abs(dot(n, &s)));
            }
            diffuse += objcol * rd * light->col * n_dot_s * alpha;

            // perfect light ray reflection
            Ray light_reflection;
            normalize(&pToLight.d);
            pToLight.d = -pToLight.d;
            rayReflect(&pToLight, p, n, &light_reflection);
            double c_dot_m = dot(&c, &(light_reflection.d));

            specular += light->col * rs *
                        pow(MAX(0, c_dot_m), obj->rt.shinyness) * alpha;
        }

        light = light->next;
    }

    if (alpha < 1) {
        // dvec = r * bvec + (r*c − srt(1 - r*r*(1-c*c))) * nvec
        point norm = *n;
        Ray refractRay;
        double s, R_Shlick;
        rayRefract(ray, obj, p, &norm, &refractRay, &s, &R_Shlick);
        if (s > 0) {
            rayTrace(&refractRay, depth + 1, &refract, obj);
            // printf("refract r: %f  g: %f  b: %f\n", refract.R, refract.G,
            // refract.B);
            refract = refract * (1 - alpha) * obj->col;
        } else {
            ambient += ambient * (1 - alpha);
            diffuse += diffuse * (1 - alpha);
            specular += specular * (1 - alpha);
        }
    }

    // perfect camera ray reflection
    Ray cam_reflection;
    rayReflect(ray, p, n, &cam_reflection);

    if (obj->refl_sig > 0) {
        rt_brandished_trace(&cam_reflection, obj, &global, depth);
    } else {
        rayTrace(&cam_reflection, depth + 1, &global, obj);
    }
    if (global.R == -1) {
        global = 0;
    } else {
        global *= rg;
    }

    col->R = MIN(1, ambient.R + diffuse.R + specular.R + refract.R + global.R);
    col->G = MIN(1, ambient.G + diffuse.G + specular.G + refract.G + global.G);
    col->B = MIN(1, ambient.B + diffuse.B + specular.B + refract.B + global.B);
    return;
}

void rt_brandished_trace(Ray *ray, Object *obj, color *col,
                         int depth) {
    Ray brandished_ray;
    memcpy(&brandished_ray, ray, sizeof(Ray));
    color brandished_color = 0;

    double num_samples = 10;
    for (int i = 0; i < num_samples; i++) {
        brandished_ray.d.x = rand_normal_dist(ray->d.x, obj->refl_sig);
        brandished_ray.d.y = rand_normal_dist(ray->d.y, obj->refl_sig);
        brandished_ray.d.z = rand_normal_dist(ray->d.z, obj->refl_sig);
        normalize(&brandished_ray.d);
        rayTrace(&brandished_ray, depth + 1, &brandished_color, obj);
        if (brandished_color.R == -1) {
            brandished_color = 0;
        }
        *col += brandished_color;
    }
    *col /= num_samples;
}
