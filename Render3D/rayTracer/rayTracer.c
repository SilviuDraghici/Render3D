#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rayTracer.h"

#include "../utils/imageProcessor.h"
#include "../utils/affinetransforms.h"

void rayTrace(int argc, char *argv[]) {
   if (argc < 5) {
      fprintf(stderr, "RayTracer: Can not parse input parameters\n");
      fprintf(stderr, "USAGE: Render3D 0 size rec_depth antialias output_name\n");
      fprintf(stderr, "   size = Image size (both along x and y)\n");
      fprintf(stderr, "   rec_depth = Recursion depth\n");
      fprintf(stderr, "   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
      fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
      exit(0);
   }


   int sx = atoi(argv[2]);
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

   for (int j = 0; j < sx; j++) {  // For each of the pixels in the image
      for (int i = 0; i < sx; i++) {
         rgbIm[3 * (j * outImage->sx + i) + 0] = (unsigned char)(255 * 1);
         rgbIm[3 * (j * outImage->sx + i) + 1] = (unsigned char)(255 * 0.25);
         rgbIm[3 * (j * outImage->sx + i) + 2] = (unsigned char)(255 * 0.25);
      }  // end for i
   }     // end for j

   fprintf(stderr, "Done!\n");

   // Output rendered image
   imageOutput(outImage, output_name);
}