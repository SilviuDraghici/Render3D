
//table top
Object *o;
struct point p;
PointLS *l;

const char *checkers = "scenes/checkers.ppm";

//const char *mesh_file = "scenes/triangular_prism_8.obj";
//const char *mesh_file = "scenes/tea_pot_992.obj";
//const char *mesh_file = "scenes/tea_pot_N_992.obj";
const char *mesh_file = "scenes/wineglass_13k.obj";
//const char *mesh_file = "scenes/bunny_69k.obj";
//const char *mesh_file = "scenes/hand_134k.obj";
//const char *mesh_file = "scenes/lucy_1M.obj";


//dragon meshes
//const char *mesh_file = "scenes/dragon_030k.obj";
//const char *mesh_file = "scenes/dragon_108k.obj";
//const char *mesh_file = "scenes/dragon_435k.obj";
//const char *mesh_file = "scenes/dragon_871k.obj";

drand48();
drand48();
drand48();

// Cornell box
scene->cam_pos = point(0, 5, -15);
scene->cam_gaze_point = point(0, -2, 0);
scene->cam_gaze = scene->cam_gaze_point - scene->cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;

scene->rt_max_depth = 5;
scene->pt_max_depth = 100;

o = new Mesh(1.0, 1.0, 1.0);
((Mesh *)o)->setMesh(mesh_file);
strcpy(o->label, "Mesh test");
//o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->set_pathTrace_properties(0.0, 0.0, 1.0);
o->r_index = 1.54;
o->T *= RotY(scene->frame * PI / 120);
o->T *= RotX(scene->frame * PI / 120);
o->T *= Sc(5.0);
o->T *= Tr(0, 0, 0);
o->invert_and_bound();
scene->insertObject(o);

o = new Plane(.75, .75, .75);
    o->set_pathTrace_properties(0.7, 0.3, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "tabletop");
    o->T *= RotX(PI / 2);
    o->T *= Sc(10,10,10);
    o->T *= Tr(0, -2.25, 0);
    o->invert_and_bound();
    scene->insertObject(o);

bool draw_box = 1;//scene->path_tracing_mode;
if (draw_box){

    // Left
    o = new Plane(.75, .25, .25);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Left Wall");
    o->T *= RotY(PI / 2);
    o->T *= Sc(10);
    o->T *= Tr(-10, 0, 0);
    o->invert_and_bound();
    scene->insertObject(o);

    // Right
    o = new Plane(.25, .25, .75);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Right Wall");
    o->T *= RotY(PI / 2);
    o->T *= Sc(25);
    o->T *= Tr(10, 0, 0);
    o->invert_and_bound();
    scene->insertObject(o);

    // Back
    o = new Plane(.75, .75, .75);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Back Wall");
    loadTexture(o, checkers, 1, scene);
    //o->T *= RotateZ(o, PI/4);
    o->T *= Sc(10);
    o->T *= Tr(0, 0, 10);
    o->invert_and_bound();
    scene->insertObject(o);

    // Front
    o = new Plane(.3, 1, .3);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Front Wall");
    //o->T *= RotateZ(o, PI/4);
    o->T *= Sc(10);
    o->T *= Tr(0, 0, -15.1);
    o->invert_and_bound();
    scene->insertObject(o);

    // Bottom
    o = new Plane(.75, .75, .75);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Bottom Wall");
    o->T *= RotX(PI / 2);
    o->T *= Sc(25);
    o->T *= Tr(0, -10, 0);
    o->invert_and_bound();
    scene->insertObject(o);

    // Top
    o = new Plane(.75, .75, .75);
    o->set_pathTrace_properties(1.0, 0.0, 0.0);
    o->r_index = 1.4;
    strcpy(o->label, "Top Wall");
    o->T *= RotX(-PI / 2);
    o->T *= Sc(25);
    o->T *= Tr(0, 10, 0);
    o->invert_and_bound();
    scene->insertObject(o);
}

p.x = 8;
p.y = 8;
p.z = 8;

// Planar light source at top
o = new Plane(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
strcpy(o->label, "Top Light");
//o->T *= Sc(0.5, 2.5, 1);
o->T *= RotX(PI / 4);
o->T *= RotY(-PI * 3 / 4 + scene->frame * PI / 120);
o->T *= Tr(p.x*1.01, p.y*1.01, p.z*1.01);
o->isLightSource = 1;
o->pt.LSweight *= 4;  // <- scale weight by scale
o->invert_and_bound();
scene->insertObject(o);

p.w = 1;
l = new PointLS(p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));

p.x = 0;
p.y = 0;
p.z = -15;
p.w = 1;
l = new PointLS(p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));
