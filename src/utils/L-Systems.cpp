#include "L-Systems.h"

color stemcol = color(0.2, 0.9, 0.25);

void new_flower(Scene *scene, color *petalCol, matrix &hierarchyMat) {
    std::cout << "new Flower!\n";

    color f;

    if (petalCol == NULL) {
        //printf("making color\n");
        f.R = drand48();
        f.G = drand48();
        f.B = drand48();
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