
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

// Cornell box
scene->cam_pos = point(0, 0, -15);
//scene->cam_gaze_point = point(0, 0, 0);
//scene->cam_gaze = cam_gaze_point - cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;

scene->rt_max_depth = 1;
scene->pt_max_depth = 20;

o = new Mesh(74 / 255.0, 255 / 255.0, 249 / 255.0);
((Mesh *)o)->setMesh(mesh_file);
strcpy(o->label, "mesh test");
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
o->T *= Sc(10.0);
o->T *= RotX(scene->frame * PI / 120);
o->T *= RotZ(scene->frame * PI / 120);
//o->T *= RotX(-PI / 2);
//o->T *= Tr(4, 0, 0);
o->invert_and_bound();
scene->insertObject(o);

/*
o = new Sphere(74 / 255.0, 255 / 255.0, 249 / 255.0);
//((Mesh *)o)->setMesh(mesh_file);
strcpy(o->label, "left");
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
o->T *= RotX(scene->frame * PI / 120);
o->T *= RotZ(scene->frame * PI / 120);
//o->T *= RotX(-PI / 2);
//o->T *= Sc(10.0);
o->T *= Tr(-4, 0, 0);
o->invert_and_bound();
insertObject(o, scene);
*/

matrix m = I();
m *= Sc(0.4);
m *= RotX(- PI/4);
m *= Tr(-1.8, 1, -1.2);
color pc = color(0.019842, 0.378377, 0.678877);
new_flower(scene, &pc, m);

int num_spheres = 100;
for (int i = 0; i < num_spheres; i++) {
    o = new Cylinder(1, 1, 0);
    strcpy(o->label, ("yellow " + std::to_string(i)).c_str());
    o->set_rayTrace_properties(.05, .95, .35, .35, 1, 6);
    o->T *= RotX(i * 2 * PI / num_spheres);
    o->T *= Sc(0.1);
    o->T *= Tr(0, 0, 5.5);
    o->T *= RotY(i * 2 * PI / num_spheres);
    o->invert_and_bound();
    scene->insertObject(o);
}

for (int i = 0; i < num_spheres; i++) {
    o = new Sphere(1, 0, 1);
    strcpy(o->label, ("purple " + std::to_string(i)).c_str());
    o->set_rayTrace_properties(.05, .95, .35, .35, 1, 6);
    o->T *= Sc(0.1);
    o->T *= Tr(0, 0, 6.5);
    o->T *= RotY(i * 2 * PI / num_spheres);
    o->T *= RotX(PI / 4);
    o->invert_and_bound();
    scene->insertObject(o);
}

for (int i = 0; i < num_spheres; i++) {
    o = new Sphere(0, 1, 1);
    strcpy(o->label, ("cyan " + std::to_string(i)).c_str());
    o->set_rayTrace_properties(.05, .95, .35, .35, 1, 6);
    o->T *= Sc(0.1);
    o->T *= Tr(0, 0, 7.5);
    o->T *= RotY(i * 2 * PI / num_spheres);
    o->T *= RotX(-PI / 4);
    o->invert_and_bound();
    scene->insertObject(o);
}

bool draw_box = scene->path_tracing_mode;
if (draw_box || true) {
    //o->set_color(1,1,1);
    //o->set_pathTrace_properties(0.0, 0.0, 1.0);

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
    loadTexture(o, checkers, 1, &t_list);
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

// Planar light source at top
o = new Plane(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
strcpy(o->label, "Top Light");
o->T *= Sc(0.5, 2.5, 1);
o->T *= RotX(PI / 2);
o->T *= Tr(0, 9.995, 0);
o->invert_and_bound();
o->isLightSource = 1;
o->pt.LSweight *= 4;  // <- scale weight by scale
scene->insertObject(o);

p.x = 0;
p.y = 9.9;
p.z = 0;
p.w = 1;
l = newPLS(&p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));

l = newPLS(&scene->cam_pos, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));
