#include <stdio.h>
#include <stdlib.h>

#include "rayTracer/rayTracer.h"

void modeError();

int main(int argc, char *argv[]) {
   if (argc < 2) {
      modeError();
   } else {
      int mode = atof(argv[1]);
      if (mode == 0) {
         rayTraceMain(argc, argv);

      } else if (mode == 1) {
         if (argc < 5) {
            fprintf(stderr, "PathTracer: Can not parse input parameters\n");
            fprintf(stderr, "USAGE: Render3D 1 size rec_depth num_samples output_name\n");
            fprintf(stderr, "   size = Image size (both along x and y)\n");
            fprintf(stderr, "   rec_depth = Recursion depth\n");
            fprintf(stderr, "   num_samples = Number of samples per pixel\n");
            fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
            exit(0);
         }
      } else {
         modeError();
      }
   }
}

void modeError() {
   fprintf(stderr, "Render3D: Can not parse input parameters\n");
   fprintf(stderr, "USAGE: Render3D mode ...\n");
   fprintf(stderr, "   mode : 0 for Ray tracing\n");
   fprintf(stderr, "          1 for Path tracing\n");
   fprintf(stderr, "Specify a mode to receive further instructions\n");
   exit(0);
}