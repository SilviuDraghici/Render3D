#include "objects.h"
#include "affineTransforms.h"
#include "utils.h"
#include "color.h"

#ifndef BUILDSCENE_H
#define BUILDSCENE_H

void buildScene(void) {
#include "buildscene.c"  // <-- Import the scene definition!
}

#endif