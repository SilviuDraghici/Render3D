#ifndef BUILDSCENE_H
#define BUILDSCENE_H

#include "objects.h"
#include "meshes.h"

typedef struct scene_struct {
    //general settings:
    int sx = 1024, sy = 1024;

    //max depth has different effect in rt vs. pt so it can
    //be set seperately
    int rt_max_depth = 4;
    int pt_max_depth = 10;

    int pt_num_samples = 1000;

    bool path_tracing_mode = 0;

    //set 1 to use antialiasing in the raytrcer
    //path tracer implicitly has antialiasing
    int rt_antialiasing = 0;

    //variables for the camera capturing the scene
    struct point cam_pos;
    struct point cam_up;
    struct point cam_gaze;
    struct point cam_gaze_point;
    double cam_focal;  // should be negative

    double du, dv;


    //variables for the scene itself
    Object *object_list = NULL;

    //this is for the ray tracer
    struct pointLS *rt_point_light_list = NULL;

    //this can be used for creating animations
    int frame = 0;
} Scene;

void buildScene(Scene *scene);

#endif