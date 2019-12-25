void buildScene(void)
{
    // Sets up all objets in the scene. This involves creating each object,
    // defining the transformations needed to shape and position it as
    // desired, specifying the reflectance properties (albedos and colours)
    // and setting up textures where needed.
    //
    // NOTE: Light sources are now EXCLUSIVELY area light sources. They are
    //       defined by regular objects whose 'isLightSource' flag is set
    //       to 1. Therefore, you can create light sources with any shape
    //       and colour using the same object primitives and transforms
    //       you're using to set up the scene.
    //
    // To create hierarchical objects:
    //    You must keep track of transformations carried out by parent objects
    //    as you move through the hierarchy. Declare and manipulate your own
    //    transformation matrices (use the provided functions in utils.c to
    //    compound transformations on these matrices). When declaring a new
    //    object within the hierarchy
    //    - Initialize the object
    //    - Apply any object-level transforms to shape/rotate/resize/move
    //      the object using regular object transformation functions
    //    - Apply the transformations passed on from the parent object
    //      by pre-multiplying the matrix containing the parent's transforms
    //      with the object's own transformation matrix.
    //    - Compute and store the object's inverse transform as usual.
    //
    // NOTE: After setting up the transformations for each object, don't
    //       forget to set up the inverse transform matrix!

    struct object3D *o;
    struct point3D p;

    // Cornell box
    o = newSphere(1.0, 0.0, 0.0, .75, .25, .25, .05, 1.4); // Left
    o->label[0] = 'L';
    o->label[1] = 'E';
    o->label[2] = 'F';
    o->label[3] = 'T';
    o->label[4] = '\0';
    Scale(o, 500, 500, 500);
    Translate(o, -510, 0, 5);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    insertObject(o, &object_list);

    o = newSphere(1.0, 0.0, 0.0, .25, .25, .75, .05, 1.4); // Right
    o->label[0] = 'R';
    o->label[1] = 'I';
    o->label[2] = 'G';
    o->label[3] = 'H';
    o->label[4] = 'T';
    o->label[5] = '\0';
    Scale(o, 500, 500, 500);
    Translate(o, 510, 0, 5);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    insertObject(o, &object_list);

    o = newSphere(1.0, 0.0, 0.0, .75, .75, .75, .05, 1.4); // Back
    o->label[0] = 'B';
    o->label[1] = 'A';
    o->label[2] = 'C';
    o->label[3] = 'K';
    o->label[4] = '\0';
    Scale(o, 500, 500, 500);
    Translate(o, 0, 0, 515);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    insertObject(o, &object_list);

    o = newSphere(1.0, 0.0, 0.0, .75, .75, .75, .02, 1.4); // Bottom
    o->label[0] = 'B';
    o->label[1] = 'O';
    o->label[2] = 'T';
    o->label[3] = 'T';
    o->label[4] = 'O';
    o->label[5] = 'M';
    o->label[6] = '\0';
    Scale(o, 500, 500, 500);
    Translate(o, 0, -510, 5);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    insertObject(o, &object_list);

    o = newSphere(1.0, 0.0, 0.0, .75, .75, .75, .05, 1.4); // Top
    o->label[0] = 'T';
    o->label[1] = 'O';
    o->label[2] = 'P';
    o->label[3] = '\0';
    Scale(o, 500, 500, 500);
    Translate(o, 0, 510, 5);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    insertObject(o, &object_list);

    // Planar light source at top
    o = newPlane(1.00, 0.00, 0.0, 1.0, 1.0, 1.0, 0.0, 1.54);
    o->label[0] = 'L';
    o->label[1] = '1';
    o->label[2] = '\0';
    Scale(o, .5, 2.5, 1);
    RotateX(o, PI / 2);
    Translate(o, 0, 9.995, 5);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    o->isLightSource = 1;
    insertObject(o, &object_list);


    char mesh_file[] = "teapot.obj";
    struct object3D *m;
    m = newMesh(mesh_file, 1.00, 0.00, 0.0, 1.0, 1.0, 1.0, 0.0, 1.54);
    //Scale(m,5, 1,1);
    //Translate(o,0,0,5);
    Translate(m, 4, -3.75, 6.5);
    invert(&m->T[0][0], &m->Tinv[0][0]);
    insertObject(m, &object_list);

    double xscale = m->max->px - m->min->px;
    double yscale = m->max->py - m->min->py;
    double zscale = m->max->pz - m->min->pz;

    double xtrans = (m->max->px - m->min->px)/2;
    double ytrans = (m->max->py - m->min->py)/2;
    double ztrans = (m->max->pz - m->min->pz)/2;

    printf("max x %f y %f z %f\n", m->max->px, m->max->py, m->max->pz);
    printf("min x %f y %f z %f\n", m->min->px, m->min->py, m->min->pz);

    printf("scale x %f y %f z %f\n", xscale, yscale, zscale);
    o=newBox(1.00, 0.00, 0.0, 1.0, 1.0, 1.0, 0.0, 1.54);
    //Translate(o, xtrans, ytrans, ztrans);
    memcpy(&o->T[0][0], &m->T[0][0], 16 * sizeof(double));

    struct point3D translate;
    translate.px = m->T[0][3];
    translate.py = m->T[1][3];
    translate.pz = m->T[2][3];

    o->T[0][3] = 0;
    o->T[1][3] = 0;
    o->T[2][3] = 0;
    Scale(o, xscale, yscale, zscale);

    o->T[0][3] = translate.px;
    o->T[1][3] = translate.py;
    o->T[2][3] = translate.pz;
    Translate(o, 0, ytrans, 0);
    //Translate(o, translate.px, translate.py, translate.pz);

    invert(&o->T[0][0], &o->Tinv[0][0]);

    printmatrix(m->T);
     //Scale(o, xscale, yscale, zscale);

    // Scale(o, 5, 1,1);

    // Translate(o, 4, -3.75, 6.5);
    //invert(&o->T[0][0], &o->Tinv[0][0]);
//   memcpy(&o->Tinv[0][0], &m->Tinv[0][0], 16 * sizeof(double));

    //insertObject(o, &object_list);

}
