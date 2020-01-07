#include "objects.h"
#include "affineTransforms.h"
#include "utils.h"
#include "color.h"
#include "mappings.h"

#ifndef BUILDSCENE_H
#define BUILDSCENE_H

void buildScene(void) {
#include "buildscenem.c"  // <-- Import the scene definition!
}

#endif