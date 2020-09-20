#ifndef SCENE_H
#define SCENE_H

#include "BVH.h"
#include "utils.h"

struct textureNode {
    char name[1024];
    int type;
    image *im;
    textureNode *next;
    ~textureNode(){
        free(im->rgbdata);
        free(im);
        delete(next);
    }
};

class Scene {
    public:
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
    point cam_pos;
    point cam_up;
    point cam_gaze;
    point cam_gaze_point;
    double cam_focal;  // should be negative

    double du, dv;

    //variables for the scene itself
    int num_objects = 0;
    Object *object_list = NULL;

    BVH *bvh = NULL;

    //this is for the ray tracer
    PointLS *rt_point_light_list = NULL;

    textureNode *texture_list = NULL;

    //this can be used for creating animations
    int frame = 0;
    void insertObject(Object *o) {
        if (o == NULL) return;
        num_objects++;
        
        // Inserts an object into the object list.
        if (object_list == NULL) {
            object_list = o;
            object_list->next = NULL;
        } else {
            o->next = object_list->next;
            object_list->next = o;
        }
    }

    ~Scene(){
        delete(rt_point_light_list);
        delete(bvh);
        delete(texture_list);
    }
};

#endif