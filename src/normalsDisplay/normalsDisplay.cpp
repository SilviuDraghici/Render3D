#include "normalsDisplay.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>
#include <sstream>

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

static Scene *scene;
void normalsDisplayMain(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Normals Display: Can not parse input parameters\n");
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

    std::string dims = argv[2];
    if (std::find(dims.begin(), dims.end(), 'x') != dims.end()) {
        std::stringstream sdims(dims);
        std::string dim;
        getline(sdims,dim, 'x');
        scene->sx = atoi(dim.c_str());
        getline(sdims, dim, 'x');
        scene->sy = atoi(dim.c_str());
    } else {
        scene->sx = atoi(argv[2]);
        scene->sy = atoi(argv[2]);
    }
    //std::cout << "Image size: " << scene->sx << "x" << scene->sy << std::endl;

    
    strcpy(&output_name[0], argv[4]);

    if (6 <= argc) {
        sc.frame = atoi(argv[5]) - 1;
    }

    double* rgbIm;
    // Allocate memory for the new image
    outImage = newImage(scene->sx, scene->sy, sizeof(double));
    if (!outImage) {
        fprintf(stderr, "Unable to allocate memory for raytraced image\n");
        exit(0);
    } else {
        rgbIm = (double* )outImage->rgbdata;
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
    double wleft,wtop,wwidth,wheight;
    scaleCameraPlainByImageDimesions(wleft, wtop, wwidth, wheight, scene->sx, scene->sy);
    cam = setupView(&(scene->cam_pos), &(scene->cam_gaze), &(scene->cam_up),
                    scene->cam_focal, wleft,wtop,wwidth,wheight);

    if (cam == NULL) {
        fprintf(stderr,
                "Unable to set up the view and camera parameters. Our of "
                "memory!\n");
        // cleanup(object_list, light_list, texture_list);
        deleteImage(outImage);
        exit(0);
    }

    setPixelStep(scene, cam, scene->sx, scene->sy);

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

    double start_j = 0, end_j = scene->sy, start_i = 0, end_i = scene->sx;

    for (int j = start_j; j < end_j; j++) {  // For each of the pixels in the image
        for (int i = start_i; i < end_i; i++) {
            getRayFromPixel(scene, &ray, cam, i, j);
            depth = 1;
            col = 0;
            rayTrace(&ray, &col, NULL);
            if (col.R == -1) {
                col = background;
            }
            rgbIm[3 * (j * outImage->sx + i) + 0] = col.R;
            rgbIm[3 * (j * outImage->sx + i) + 1] = col.G;
            rgbIm[3 * (j * outImage->sx + i) + 2] = col.B;
        }  // end for i
    }      // end for j
    raytracing_timer.end();
    // Output rendered image
    PNGImageOutput(outImage, output_name);

    std::cout << "Ray Tracing Done!\n";

    free(cam);
    
    //debug :
    raytracing_timer.print_elapsed_time(std::cerr);
    std::cerr << "Average number of intersection tests: " << num_intersection_tests / num_bvh_searches << "\n";
}

void rayTrace(Ray *ray, color *col, Object *Os) {
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

    findFirstHit(scene, ray, &lambda, Os, &obj, &p, &n, &a, &b);

    if (lambda == -1 || obj == NULL) {
        *col = -1;
        return;
    }
    col->R = (n.x + 1) / 2;
    col->G = (n.y + 1) / 2;
    col->B = (n.z + 1) / 2;
}
