#ifndef L_SYSTEMS_H
#define L_SYSTEMS_H

#include "scene.h"

void new_Flower(Scene *scene, matrix &hierarchyMat, color *petalCol);
void new_Branch(Scene *scene, matrix &hierarchyMat, color *col, double numBranches, double maxRotation, double curr_depth, double depth);
void new_FlTree(Scene *scene, matrix &hierarchyMat, color *col, double distFromC, double numBranches, double maxRotation, double depth);
#endif