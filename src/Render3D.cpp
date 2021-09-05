#include <stdio.h>
#include <stdlib.h>

#include "pathTracer/pathTracer.h"
#include "rayTracer/rayTracer.h"
#include "normalsDisplay/normalsDisplay.h"
#include "colorsDisplay/colorsDisplay.h"

void modeError();

int main(int argc, char *argv[]) {
    if (argc < 2) {
        modeError();
    } else {
        int mode = atof(argv[1]);
        if (mode == 0) {
            rayTraceMain(argc, argv);
        } else if (mode == 1) {
            pathTraceMain(argc, argv);
        } else if (mode == 2) {
            normalsDisplayMain(argc, argv);
        } else if (mode == 3) {
            colorsDisplayMain(argc, argv);
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