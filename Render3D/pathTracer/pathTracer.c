#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pathTracer.h"

#include "../utils/affinetransforms.h"
#include "../utils/buildscene.h"
#include "../utils/camera.h"
#include "../utils/color.h"
#include "../utils/imageProcessor.h"
#include "../utils/mappings.h"
#include "../utils/ray.h"
#include "../utils/utils.h"
#include "../utils/objects.h"

static int MAX_DEPTH;

void pathTraceMain(int argc, char *argv[]){
     double *wght;         // Holds weights for each pixel - to provide log response


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

   
   unsigned char *rgbIm;
   // Allocate memory for the new image
   outImage = newImage(sx, sx);
   if (!outImage) {
      fprintf(stderr, "Unable to allocate memory for raytraced image\n");
      exit(0);
   } else {
      rgbIm = (unsigned char *)outImage->rgbdata;
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

   struct view *cam;  // Camera and view for this scene
   cam = setupView(&cam_pos, &cam_gaze, &cam_up, cam_focal, -2, 2, 4);

    if (cam == NULL) {
      fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
      //cleanup(object_list, light_list, texture_list);
      deleteImage(outImage);
      exit(0);
   }

   setPixelStep(cam, sx, sx);

   fprintf(stderr, "Path Tracing Done!\n");

   

   // Output rendered image
   imageOutput(outImage, output_name);
}