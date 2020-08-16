#ifndef BUILDSCENE_H
#define BUILDSCENE_H

#include "objects.h"
#include "scene.h"
#include "meshes.h"

//int min_for_bvh = 10;

void buildScene(Scene *scene);
void insertObject(Object *o, Scene *scene);

#endif