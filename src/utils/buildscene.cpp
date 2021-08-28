#include "buildscene.h"

#include "affineTransforms.h"
#include "camera.h"
#include "color.h"
#include "mappings.h"
#include "objects.h"
#include "L-Systems.h"
#include "utils.h"

void buildScene(Scene *scene) {
#include "../../scenes/buildScene_sdb.cpp"
    
    PrimitiveData *prims = (PrimitiveData *)malloc(scene->object_list.size() * sizeof(PrimitiveData));
    scene->bvh = new BVH;
    int i = 0;
    
    for(Object* const& prim : scene->object_list){
        prims[i] = prim;
        i++;
    }

    scene->bvh->set_build_method(BuildMethod::SAH);
    scene->bvh->set_search_method(SearchMethod::BFS);
    scene->bvh->build(prims, scene->object_list.size());
    //scene->bvh->print();
    std::cout << scene->object_list.size() << " objects in scene\n";
    //point min = point(scene->bvh->root->min_x(), scene->bvh->root->min_y(), scene->bvh->root->min_z());
    //point max = point(scene->bvh->root->max_x(), scene->bvh->root->max_y(), scene->bvh->root->max_z());
    //std::cout << "b bound:"<< min << " " << max << "\n";
    free(prims);
}
