#include "L-Systems.h"
#include "random.h"

#define b_to_f 0.666
double yscale = 2;
color stemcol = color(0.2, 0.9, 0.25);
matrix stemMat = Sc(0.25, yscale, 0.25) * RotX(-PI / 2);

void new_Flower(Scene *scene, matrix &hierarchyMat, color *petalCol) {
    color f;

    if (petalCol == NULL) {
        //printf("making color\n");
        f.R = xor128();
        f.G = xor128();
        f.B = xor128();
        petalCol = &f;
    }
    //printf("flower.R = %f, flower.G = %f, flower.B = %f;\n", petalCol->R,petalCol->G, petalCol->B);
    Object *o;
    //stem

    double fscale = 3;  // the initial flower was too small
    double stemscale = 0.5;
    //o = newCyl(.05, .95, 0, 0, stemR, stemG, stemB, 1, 1, 6);
    o = new Cylinder(stemcol);
    o->set_rayTrace_properties(.05, .95, 0, 0, 1, 6);
    o->T *= RotX(-PI / 2);
    o->T *= Sc(fscale * 0.1, fscale * stemscale, fscale * 0.1);
    o->T *= hierarchyMat;
    o->invert_and_bound();
    scene->insertObject(o);

    //newSphere(ra,  rd, rs, rg,  r, g, b, alpha, r_index, shiny)
    //o = newSphere(.05, .95, .95, .75, 1 - petalCol->R, 1 - petalCol->G, 1 - petalCol->B, 1, 1, 6);
    color cinv = petalCol->inverse();
    o = new Sphere(cinv);
    o->set_rayTrace_properties(.05, .95, .95, .75, 1, 6);
    o->T *= Sc(fscale * 0.3, fscale * 0.1, fscale * 0.3);
    o->T *= Tr(0, fscale * stemscale, 0);
    o->T *= hierarchyMat;
    o->invert_and_bound();
    scene->insertObject(o);

    for (int i = 0; i < 5; i++) {
        //o = newSphere(.05, .95, 0, 0, petalCol->R, petalCol->G, petalCol->B, 1, 1, 6);
        o = new Sphere(*petalCol);
        o->set_rayTrace_properties(.05, .95, 0, 0, 1, 6);
        o->T *= Tr(0.5, 0, 0);
        o->T *= Sc(fscale * 0.5, fscale * 0.005, fscale * 0.3);
        o->T *= RotZ(0.25);
        o->T *= Tr(0.1, fscale * stemscale - 0.04, 0);
        o->T *= RotY(i * 2 * PI / 5);
        o->T *= hierarchyMat;
        o->invert_and_bound();
        scene->insertObject(o);
    }
}

void new_Branch(Scene *scene, matrix &hierarchyMat, color *col, double numBranches, double maxRotation, double curr_depth, double depth) {
    Object *o;
    //this is what the default shape of the cylinder should be
    //not part of the hierarchy
    o = new Cylinder(stemcol);
    o->set_rayTrace_properties(.05, .95, 0, 0, 1, 6);
    o->T *= stemMat;
    o->T *= hierarchyMat;
    o->invert_and_bound();
    scene->insertObject(o);

    if (curr_depth + 1 >= depth) {
        return;
    }

    matrix newTransforms;
    double angle = xor128() * 2 * maxRotation;
    double dice;
    color f;

    for (int i = 0; i < numBranches; i++) {
        newTransforms = I();
        newTransforms *= Sc(0.8);
        newTransforms *= Tr(0, yscale / 2 + 0.5, 0);

        newTransforms *= RotX(PI / 4);
        newTransforms *= RotY(angle + i * 2 * PI / numBranches);
        newTransforms *= Tr(0, yscale / 2 + 0.5, 0);
        newTransforms *= hierarchyMat;

        dice = xor128();
        if (dice <= b_to_f) {  //new flower
            new_Flower(scene, newTransforms, col);
        } else {
            new_Branch(scene, newTransforms, col, numBranches, maxRotation, curr_depth + 1, depth);
        }
    }
}

void new_FlTree(Scene *scene, matrix &hierarchyMat, color *col, double distFromC, double numBranches, double maxRotation, double depth) {
    Object *o;

    //this is what the default shape of the cylinder should be
    //not part of the hierarchy
    o = new Cylinder(stemcol);
    o->set_rayTrace_properties(.05, .95, 0, 0, 1, 6);
    o->T *= stemMat;
    o->T *= hierarchyMat;
    o->invert_and_bound();
    scene->insertObject(o);

    if (depth <= 1) {
        return;
    }

    matrix newTransforms;
    double angle = xor128() * 2 * maxRotation;
    for (int i = 0; i < numBranches; i++) {
        newTransforms = I();
        newTransforms *= Sc(0.8);
        newTransforms *= Tr(0, yscale / 2 + 0.5, 0);
        newTransforms *= RotX(PI / 4);
        newTransforms *= RotY(angle + i * 2 * PI / numBranches);
        newTransforms *= Tr(0, yscale / 2 + 0.5, 0);
        newTransforms *= hierarchyMat;
        new_Branch(scene, newTransforms, col, numBranches, maxRotation, 1, depth);
    }
}