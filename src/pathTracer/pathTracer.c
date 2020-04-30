#include "pathTracer.h"

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../utils/affinetransforms.h"
#include "../utils/buildscene.h"
#include "../utils/camera.h"
#include "../utils/color.h"
#include "../utils/imageProcessor.h"
#include "../utils/mappings.h"
#include "../utils/objects.h"
#include "../utils/ray.h"
#include "../utils/utils.h"

int samples_per_update = 25;

static int MAX_DEPTH;
unsigned long int NUM_RAYS;

double total_weight;

//array of light sources
struct object **light_listt;
int num_lights;
int curr_light;

void pathTraceMain(int argc, char *argv[]) {
    struct color col;  // Return colour for pixels

    if (argc < 6) {
        fprintf(stderr, "PathTracer: Can not parse input parameters\n");
        fprintf(stderr, "USAGE: Render3D 1 size rec_depth num_samples output_name\n");
        fprintf(stderr, "   size = Image size (both along x and y)\n");
        fprintf(stderr, "   rec_depth = Recursion depth\n");
        fprintf(stderr, "   num_samples = Number of samples per pixel\n");
        fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
        exit(0);
    }

    int sx = atoi(argv[2]);
    MAX_DEPTH = atoi(argv[3]);
    int num_samples = atoi(argv[4]);
    strcpy(&output_name[0], argv[5]);

    double *rgbIm;
    // Allocate memory for the new image
    outImage = newImage(sx, sx, sizeof(double));

    if (!outImage) {
        fprintf(stderr, "Unable to allocate memory for pathtraced image\n");
        exit(0);
    } else {
        rgbIm = (double *)outImage->rgbdata;
    }

    //array of weights per pixel used to scale image colors
    double *wght;  // Holds weights for each pixel - to provide log response
    double wt;
    wght = (double *)calloc(sx * sx, sizeof(double));
    if (!wght) {
        fprintf(stderr, "Unable to allocate memory for pathTracw]e weights\n");
        exit(0);
    }
    for (int i = 0; i < sx * sx; i++) {
        *(wght + i) = 1.0;
    }

    // Camera center is at (0,0,-1)
    cam_pos = point(0, 0, -1);

    // To define the gaze vector, we choose a point 'pc' in the scene that
    // the camera is looking at, and do the vector subtraction pc-e.
    // Here we set up the camera to be looking at the origin.

    cam_gaze_point = point(0, 0, 0);
    cam_gaze = cam_gaze_point - cam_pos;

    cam_up = point(0, 1, 0);

    cam_focal = -1;

    buildScene();

    //count number of lights
    num_lights = 0;
    struct object *curr_obj = object_list;
    while (curr_obj != NULL) {
        if (curr_obj->isLightSource) {
            num_lights++;
        }
        curr_obj = curr_obj->next;
    }

    //create array of light pointers
    light_listt = (struct object **)malloc(num_lights * sizeof(struct object *));
    num_lights = 0;
    curr_obj = object_list;
    while (curr_obj != NULL) {
        if (curr_obj->isLightSource) {
            light_listt[num_lights] = curr_obj;
            num_lights++;
        }
        curr_obj = curr_obj->next;
    }
    curr_light = 0;

    struct view *cam;  // Camera and view for this scene
    cam = setupView(&cam_pos, &cam_gaze, &cam_up, cam_focal, -2, 2, 4);

    if (cam == NULL) {
        fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
        //cleanup(object_list, light_listt, texture_list);
        deleteImage(outImage);
        exit(0);
    }

    setPixelStep(cam, sx, sx);

    normalizeLightWeights(object_list);

    NUM_RAYS = 0;

    time_t t1, t2;
    t1 = time(NULL);

    fprintf(stderr, "Rendering...\n");
    struct ray ray;
    int k, j, i;
    double is, js;

    t1 = time(NULL);

    for (k = 1; k <= num_samples; k++) {
        fprintf(stderr, "\r%d/%d", k, num_samples);
        //fflush(stderr);
#pragma omp parallel for schedule(dynamic, 1) private(i, j, wt, ray, col, is, js)
        for (j = 0; j < sx; j++) {  // For each of the pixels in the image
            for (i = 0; i < sx; i++) {
                col = 0;

                // Random sample within the pixel's area
                is = i + drand48() - 0.5;
                js = j + drand48() - 0.5;

                getRayFromPixel(&ray, cam, is, js);

                ray.pt.ray_col = 1;
                ray.pt.expl_col = 0;
                ray.pt.isLightRay = 0;

                wt = *(wght + i + (j * sx));

                PathTrace(&ray, 1, &col, NULL, NULL);
                rgbIm[3 * (j * outImage->sx + i) + 0] += col.R * pow(2, -log(wt));
                rgbIm[3 * (j * outImage->sx + i) + 1] += col.G * pow(2, -log(wt));
                rgbIm[3 * (j * outImage->sx + i) + 2] += col.B * pow(2, -log(wt));

                wt += col.R;
                wt += col.G;
                wt += col.B;
                *(wght + i + (j * sx)) = wt;
            }  // end for i
        }      // end for j

        if (k % samples_per_update == 0) {  // update output image
            dataOutput(rgbIm, sx, output_name);
        }

    }  // End for k

    t2 = time(NULL);

    // Output rendered image
    dataOutput(rgbIm, sx, output_name);

    fprintf(stderr, "\nPath Tracing Done!\n");

    fprintf(stderr, "Total number of rays created: %ld\n", NUM_RAYS);
    fprintf(stderr, "Rays per second: %.0f\n", (double)NUM_RAYS / (double)difftime(t2, t1));
}

void PathTrace(struct ray *ray, int depth, struct color *col, struct object *Os, struct object *explicit_l) {
    // Trace one light path through the scene.
    //
    // Parameters:
    //   *ray   -  A pointer to the ray being traced
    //   depth  -  Current recursion depth for recursive raytracing
    //   *col   - Pointer to an RGB colour structure so you can return the object colour
    //            at the intersection point of this ray with the closest scene object.
    //   *Os    - 'Object source' is a pointer to the object from which the ray
    //            originates so you can discard self-intersections due to numerical
    //            errors. NULL for rays originating from the center of projection.
    NUM_RAYS++;
    double lambda;              // Lambda at intersection
    double a, b;                // Texture coordinates
    struct object *obj = NULL;  // Pointer to object at intersection
    struct point p;             // Intersection point
    struct point n;             // Normal at intersection
    struct point d;
    struct color objcol;  // Colour for the object in R G and B
    double diffuse, reflect, refract;
    double R_Shlick = 1;
    double dice;  // Handy to keep a random value
    double max_col;

    struct object *explt = explicit_l;

    if (depth > MAX_DEPTH)  // Max recursion depth reached. Return black (no light coming into pixel from this path).
    {
        *col = ray->pt.expl_col;  // These are accumulators, initialized at 0. Whenever we find a source of light these
                                  // get incremented accordingly. At the end of the recursion, we return whatever light
                                  // we accumulated into these three values.
        return;
    }

    findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
    if (obj == NULL || lambda < THR) {
        *col = ray->pt.expl_col;
        return;
    }

    // set object color
    textureMap(obj, a, b, &objcol);

    // check for normal map
    normalMap(obj, a, b, &n);

    // check for alpha map
    if (obj->alphaMap == NULL) {
        diffuse = obj->pt.diffuse;
        reflect = obj->pt.reflect;
        refract = obj->pt.refract;
    } else {
        alphaMap(obj, a, b, &diffuse, obj->pt.diffuse);
        refract = 1 - diffuse;
        reflect = obj->pt.reflect;

        //scale down so total is still 1
        diffuse *= (1 - reflect);
        refract *= (1 - reflect);
    }

    //all the other terms use this
    ray->pt.ray_col *= objcol;

    // if hit light source
    if (obj->isLightSource) {
        *col = ray->pt.ray_col + ray->pt.expl_col;
        if (ray->pt.isLightRay == 0) {
            if (explt == obj) {  // the ray cast is the same as explicit
                *col = ray->pt.expl_col;
            }
        }
        return;
    }

    //make sure normal is on side of incoming ray
    if (obj->frontAndBack && dot(&ray->d, &n) > 0) {
        //printf("flip\n");
        n.x *= -1;
        n.y *= -1;
        n.z *= -1;
    }

    memcpy(&ray->p0, &p, sizeof(struct point));

    dice = drand48();
    if (dice <= diffuse) {  //diffuse
        // ******************** Importance sampling ************************
        cosWeightedSample(&n, &d);
        ray->d = d;
        // ******************** Random sampling ************************
        //hemiSphereCoordinates(&n, &d);

        //double d_dot_n = dot(&d, &n);
        //ray->d = d * d_dot_n;

        // explicit light sampling
        if (ray->pt.isLightRay == 0) {
            // ray from intersection point to light source
            struct ray pToLight;
            pToLight.p0 = p;

            double prob = 0;
            dice = drand48();
            for (int i = 0; i < num_lights; i++) {
                prob += light_listt[i]->pt.LSweight;
                if (dice < prob) {
                    curr_light = i;
                    break;
                }
            }

            double x, y, z;
            light_listt[curr_light]->randomPoint(light_listt[curr_light], &x, &y, &z);
            pToLight.d.x = x - p.x;  //0 - p.px;
            pToLight.d.y = y - p.y;  //9.95 - p.py;
            pToLight.d.z = z - p.z;  //5 - p.pz;
            pToLight.d.w = 1;

            double light_lambda;
            struct object *obstruction = NULL;
            struct point lightp;
            struct point nls;
            double La, Lb;

            findFirstHit(&pToLight, &light_lambda, obj, &obstruction, &lightp, &nls,
                         &La, &Lb);

            //printf("source: %s, obstruction: %s\n", obj->label, obstruction->label);
            explt = light_listt[curr_light];
            if (obstruction == light_listt[curr_light]) {
                //printf("total weight: %f\n", total_weight);
                double A = total_weight * light_listt[curr_light]->pt.LSweight;
                double dxd = pToLight.d.x * pToLight.d.x + pToLight.d.y * pToLight.d.y + pToLight.d.z * pToLight.d.z;
                normalize(&pToLight.d);
                double n_dot_l = fabs(dot(&n, &pToLight.d));
                double nls_dot_l = fabs(dot(&nls, &pToLight.d));
                double w = MIN(1, (A * n_dot_l * nls_dot_l) / (dxd));

                struct color light_col;
                // set light color
                textureMap(light_listt[curr_light], La, Lb, &light_col);

                ray->pt.expl_col += ray->pt.ray_col * light_col * w;
            }
        }
        // end of explicit light sampling
    } else if (dice <= diffuse + reflect) {  //reflective
        explt = NULL;
        struct ray ray_reflected;
        rayReflect(ray, &p, &n, &ray_reflected);
        //burnished reflection
        ray->d.x = rand_normal_dist(ray_reflected.d.x, obj->refl_sig);
        ray->d.y = rand_normal_dist(ray_reflected.d.y, obj->refl_sig);
        ray->d.z = rand_normal_dist(ray_reflected.d.z, obj->refl_sig);
        normalize(&ray->d);
    } else {  //refractive
        explt = NULL;
        struct ray refractRay;
        double s, R_Shlick;
        rayRefract(ray, obj, &p, &n, &refractRay, &s, &R_Shlick);
        if (s > 0 && drand48() < (1 - R_Shlick)) { /*keep refraction*/
        } else {
            // refract the ray because of total internal reflection (s <= 0)
            // or dice role resulted in accumulating relected portion (drand48() > (1 - R_Shlick))
            rayReflect(ray, &p, &n, &refractRay);
        }
        ray->d = refractRay.d;
    }

    dice = drand48();
    max_col = MAX(MAX(ray->pt.ray_col.R, ray->pt.ray_col.G), MAX(ray->pt.ray_col.R, ray->pt.ray_col.B));
    if (sqrt(dice) < max_col) {
        return PathTrace(ray, depth + 1, col, obj, explt);
    } else {
        *col = ray->pt.expl_col;
        return;
    }
}

void normalizeLightWeights(struct object *object_list) {
    // Update light source weights - will give you weights for each light source that add up to 1
    struct object *obj = object_list;
    total_weight = 0;
    while (obj != NULL) {
        if (obj->isLightSource)
            total_weight += obj->pt.LSweight;
        obj = obj->next;
    }
    obj = object_list;
    while (obj != NULL) {
        if (obj->isLightSource) {
            obj->pt.LSweight /= total_weight;
        }
        obj = obj->next;
    }
}