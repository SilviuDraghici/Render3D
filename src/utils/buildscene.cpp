#include "buildscene.h"

#include "affineTransforms.h"
#include "camera.h"
#include "color.h"
#include "mappings.h"
#include "objects.h"
#include "L-Systems.h"
#include "utils.h"

void buildScene(Scene *scene) {
#include "../scenes/buildscene_cbt.cpp"
    
    PrimitiveData *prims = (PrimitiveData *)malloc(scene->num_objects * sizeof(PrimitiveData));
    scene->bvh = new BVH;
    int i = 0;
    Object *prim = scene->object_list;
    while(prim != NULL){
        prims[i] = prim;
        i++;
        prim = prim->next;
    }
    scene->bvh->set_build_method(BuildMethod::MidSplit);
    scene->bvh->set_search_method(SearchMethod::BFS);
    scene->bvh->build(prims, scene->num_objects);
    //scene->bvh->print();
    std::cout << scene->num_objects << " objects in scene\n";
    //point min = point(scene->bvh->root->min_x(), scene->bvh->root->min_y(), scene->bvh->root->min_z());
    //point max = point(scene->bvh->root->max_x(), scene->bvh->root->max_y(), scene->bvh->root->max_z());
    //std::cout << "b bound:"<< min << " " << max << "\n";
}

void insertObject(Object *o, Scene *scene) {
    Object **list = &scene->object_list;
    if (o == NULL) return;
    scene->num_objects++;
    // Inserts an object into the object list.
    if (*(list) == NULL) {
        *(list) = o;
        (*(list))->next = NULL;
    } else {
        o->next = (*(list))->next;
        (*(list))->next = o;
    }
}