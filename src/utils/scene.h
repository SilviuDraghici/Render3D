#ifndef SCENE_H
#define SCENE_H

#include <list>

#include "BVH.h"
#include "utils.h"
#include "textureNode.h"
#include "MeshFactory.h"

class Scene {
    public:

    Scene():
    meshFactory(object_list, texture_list)
    {}

    //general settings:
    int sx = 1024, sy = 1024;

    //max depth has different effect in rt vs. pt so it can
    //be set seperately
    int rt_max_depth = 4;
    int pt_max_depth = 10;

    int pt_num_samples = 1000;

    bool path_tracing_mode = 0;

    double exposure = 1.0;

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
    std::list<Object *> object_list;

    BVH *bvh = NULL;

    //this is for the ray tracer
    PointLS *rt_point_light_list = NULL;

    std::list<textureNode *> texture_list;

    MeshFactory meshFactory;

    //this can be used for creating animations
    int frame = 0;

    void insertObject(Object *o) {
        if (o == NULL) return;
        object_list.push_back(o);
    }

    ~Scene(){
        //std::cout << "Objects:\n";
        //for(Object* obj: object_list){
        //    std::cout << obj->name << "\n";
        //}
        delete(rt_point_light_list);
        delete(bvh);
        for (textureNode * t_node : texture_list){
            delete(t_node);
        }
        texture_list.clear();
    }
};

#endif