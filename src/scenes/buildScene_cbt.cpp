
Object *o;
struct point p;
struct pointLS *l;

struct textureNode *t_list = NULL;

const char *checkers = "scenes/checkers.ppm";

//const char *mesh_file = "scenes/triangular_prism_8.obj";
//const char *mesh_file = "scenes/tea_pot_992.obj";
//const char *mesh_file = "scenes/tea_pot_N_992.obj";
//const char *mesh_file = "scenes/wineglass_13k.obj";
//const char *mesh_file = "scenes/bunny_69k.obj";
//const char *mesh_file = "scenes/hand_134k.obj";
const char *mesh_file = "scenes/lucy_1M.obj";


//dragon meshes
//const char *mesh_file = "scenes/dragon_030k.obj";
//const char *mesh_file = "scenes/dragon_108k.obj";
//const char *mesh_file = "scenes/dragon_435k.obj";
//const char *mesh_file = "scenes/dragon_871k.obj";

drand48();
drand48();
drand48();

// Cornell box
scene->cam_pos = point(0, 0, -15);
//scene->cam_gaze_point = point(0, 0, 0);
//scene->cam_gaze = cam_gaze_point - cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;

scene->rt_max_depth = 1;
scene->pt_max_depth = 20;

o = new Mesh(74/255.0, 255/255.0, 249/255.0);
((Mesh *)o)->setMesh(mesh_file);
strcpy(o->label, "Mesh test");
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
o->T *= RotX(scene->frame * PI / 120);
o->T *= RotZ(scene->frame * PI / 120);
//o->T *= RotX(-PI / 2);
o->T *= Sc(10.0);
o->T *= Tr(0, 0, 5.5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

bool draw_box = scene->path_tracing_mode;
if (draw_box){

    o->set_color(1,1,1);
    o->set_pathTrace_properties(0.0, 0.0, 1.0);

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