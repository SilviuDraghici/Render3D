#ifndef MESHFACTORY_H
#define MESHFACTORY_H

#include <list>

#include "utils.h"
#include "objects.h"
#include "BVH.h"
#include "meshes.h"
#include "color.h"
#include "textureNode.h"
#include "material.h"

#define MAX_REFL_SIG 0.3

class MeshFactory{
  public:
    MeshFactory(std::list<Object*>& o_l, std::list<textureNode*>& t_l);
    ~MeshFactory();
    void setDefaultColor(color& c);
    void setDefaultColor(double r, double g, double b);
    void setTransform(matrix m);
    void loadMeshFile(const std::string& filename);

  private:
    Mesh* mesh = NULL;
    std::list<Object*>& object_list;
    std::list<textureNode*>& texture_list;

    int num_faces_in_object, first_face_in_object;
    std::string mtl_name;
    std::string object_name;

    TriangleFace* buildFace(std::string& line);
    void loadMaterialFile(const std::string& mtllib_line);
    void buildMesh();
    
    int num_vertices;
    point *vertices = NULL;

    int num_texture_coords;
    double *texture_coords = NULL;

    int num_normals;
    point *normals = NULL;

    int num_faces;
    PrimitiveData *faces = NULL;

    std::string dir;
    std::string material_file_name;
    std::list<std::string> mtl_file_list;
    std::list<material> mtl_list;
    
    color default_color = color(74 / 255.0, 255 / 255.0, 249 / 255.0);
    matrix transformation;
};
#endif