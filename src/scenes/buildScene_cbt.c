
Object *o;
struct point p;
struct pointLS *l;

struct textureNode *t_list = NULL;

const char *checkers = "scenes/checkers.ppm";

//const char *mesh_file = "scenes/triangular_prism.obj";
const char *mesh_file = "scenes/tea_pot.obj";
//const char *mesh_file = "scenes/bunny.obj";

drand48();
drand48();
drand48();

// Cornell box
scene->cam_pos = point(0, 0, -15);
//scene->cam_gaze_point = point(0, 0, 0);
//scene->cam_gaze = cam_gaze_point - cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;

o = new Mesh(1.0, 1.0, 0.0);
((Mesh *)o)->setMesh(mesh_file);
strcpy(o->label, "Mesh test");
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
o->T *= RotY(scene->frame * PI / 120);
o->T *= RotX(scene->frame * PI / 120);
o->T *= Sc(15.0);
o->T *= Tr(0, 0, 5.5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

int draw_box = 0;
if (draw_box){
    // Left
    o = new Plane(.75, .25, .25);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Left Wall");
    o->T *= RotY(PI / 2);
    o->T *= Sc(25);
    o->T *= Tr(-10, 0, 5);
    invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
    insertObject(o, &(scene->object_list));

    // Right
    o = new Plane(.25, .25, .75);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Right Wall");
    o->T *= RotY(PI / 2);
    o->T *= Sc(25);
    o->T *= Tr(10, 0, 5);
    invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
    insertObject(o, &(scene->object_list));

    // Back
    o = new Plane(.75, .75, .75);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Back Wall");
    loadTexture(o, checkers, 1, &t_list);
    //o->T *= RotateZ(o, PI/4);
    o->T *= Sc(10);
    o->T *= Tr(0, 0, 15);
    invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
    insertObject(o, &(scene->object_list));

    // Front
    o = new Plane(.3, 1, .3);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Front Wall");
    //o->T *= RotateZ(o, PI/4);
    o->T *= Sc(10);
    o->T *= Tr(0, 0, -15.1);
    invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
    insertObject(o, &(scene->object_list));

    // Bottom
    o = new Plane(.75, .75, .75);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Bottom Wall");
    o->T *= RotX(PI / 2);
    o->T *= Sc(25);
    o->T *= Tr(0, -10, 5);
    invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
    insertObject(o, &(scene->object_list));

    // Top
    o = new Plane(.75, .75, .75);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Top Wall");
    o->T *= RotX(-PI / 2);
    o->T *= Sc(25);
    o->T *= Tr(0, 10, 5);
    invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
    insertObject(o, &(scene->object_list));
}

// Planar light source at top
o = new Plane(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
strcpy(o->label, "Top Light");
o->T *= Sc(0.5, 2.5, 1);
o->T *= RotX(PI / 2);
o->T *= Tr(0, 9.995, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
o->isLightSource = 1;
o->pt.LSweight *= 4;  // <- scale weight by scale
insertObject(o, &(scene->object_list));

p.x = 0;
p.y = 9.9;
p.z = 5;
p.w = 1;
l = newPLS(&p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));

p.x = 0;
p.y = 0;
p.z = -15;
p.w = 1;
l = newPLS(&p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));
