#include "MeshFactory.h"

#include <fstream>
#include <iostream>
#include <algorithm>

inline int count_vertices(std::string &face){
    int n = 0;
    for(int i = 1; i < face.size(); i++){
        if(face[i] != ' ')
            n++;
        while( i < face.size() && face[i] != ' ')
          i++;
    }
    return n;
}

MeshFactory::MeshFactory(std::list<Object*>& o_l, std::list<textureNode*>& t_l):
object_list(o_l), texture_list(t_l){}

MeshFactory::~MeshFactory(){

}

void MeshFactory::setDefaultColor(color& c){
    default_color = c;
}

void MeshFactory::setDefaultColor(double r, double g, double b){
    default_color = color(r, g, b);
}

void MeshFactory::setTransform(matrix& m){
    transformation = m;
}

void MeshFactory::setTransform(matrix m){
    transformation = m;
}

void MeshFactory::buildMesh(const std::string& filename){
    mesh = new Mesh(default_color);
    mesh->T = transformation;
    mesh->set_pathTrace_properties(1.0, 0.0, 0.0);
    mesh->r_index = 1.54;

    int last_dir_sep = filename.find_last_of("/");
    if(last_dir_sep != std::string::npos){
        dir = filename.substr(0, last_dir_sep + 1);
    } else if ((last_dir_sep = filename.find_last_of("\\")) != std::string::npos){
        dir = filename.substr(0, last_dir_sep + 1);
    } else {
        dir = "";
    }
    //std::cout << "dir: " << dir << "\n";
    
    std::string line;
    std::ifstream mesh_obj(filename);

    num_vertices = 0, num_faces = 0, num_normals = 0;

    int v = 1, t = 1, vn = 1, f = 0;
    double x, y, z, scale;
    double min_x, max_x, avg_x, min_y, max_y, avg_y, min_z, max_z, avg_z;
    min_x = min_y = min_z = INFINITY;
    max_x = max_y = max_z = -INFINITY;
    avg_x = avg_y = avg_z = 0;

    if (!mesh_obj.is_open()) {
        std::cout << "Unable to open mesh file: " << filename << "\n";
        return;
    }

    while (getline(mesh_obj, line)) {
        if(line.find('\r') != std::string::npos){
            //remove '\r' from windows format objs
            line.erase(line.find('\r'));
        }

        if(line.rfind("mtllib ", 0) == 0){
            setMaterial(line);
        } else if (line.rfind("v ", 0) == 0) {
            //count vertices and update average location
            num_vertices++;
            sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
            avg_x += x, avg_y += y, avg_z += z;
            if (x < min_x) {
                min_x = x;
            }
            if (max_x < x) {
                max_x = x;
            }

            if (y < min_y) {
                min_y = y;
            }
            if (max_y < y) {
                max_y = y;
            }

            if (z < min_z) {
                min_z = z;
            }
            if (max_z < z) {
                max_z = z;
            }
        } else if (line.rfind("vn ", 0) == 0) {
            //count vertex normals
            //std::cout << line << "\n";
            num_normals++;
        } else if (line.rfind("vt ", 0) == 0){
            num_texture_coords++;
        }else if (line.rfind("f ", 0) == 0) {
            //count faces
            num_faces += count_vertices(line) - 2;
        }
    }

    //std::cout << "num vertices: " << num_vertices << "\n";
    //std::cout << "num normals: " << num_normals << "\n";
    //std::cout << "num faces: " << num_faces << "\n";

    //go to begining of file
    mesh_obj.clear();
    mesh_obj.seekg(0);

    // calculate average distance from 0 so that mesh can be centered
    avg_x /= num_vertices, avg_y /= num_vertices, avg_z /= num_vertices;
    // calculate max dimension so mesh can be scaled down to 1x1x1 ish
    scale = MAX(abs(max_x - min_x), MAX(abs(max_y - min_y), abs(max_z - min_z)));

    num_vertices += 1;  // << overcount by 1
    vertices = (point *)malloc(num_vertices * sizeof(point));

    num_texture_coords +=1;
    if (num_texture_coords > 1) {
        //note one texture coord is 2 doubles
        texture_coords = (double *)malloc(num_texture_coords * 2 * sizeof(double));
    }

    num_normals += 1;
    if (num_normals > 1) {
        normals = (point *)malloc(num_normals * sizeof(point));
    }

    faces = (PrimitiveData *)malloc(num_faces * sizeof(PrimitiveData));
    int num_vertices;
    int start_index = 2, end_index;
    std::string face_string, first_vertex;
    while (getline(mesh_obj, line)) {
        if(line.find('\r') != std::string::npos){
            //remove '\r' from windows format objs
            line.erase(line.find('\r'));
        }
        if (line.rfind("v ", 0) == 0) {
            sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
            x = (x - avg_x) / scale, y = (y - avg_y) / scale,
            z = -(z - avg_z) / scale;
            vertices[v] = point(x, y, z);
            v++;
        } else if (line.rfind("vt ", 0) == 0){
            sscanf(line.c_str(), "vt %lf %lf", &x, &y);
            texture_coords[t *2] = x;
            texture_coords[(t * 2) + 1] = y;
            t++;
        } else if (line.rfind("vn ", 0) == 0){
            sscanf(line.c_str(), "vn %lf %lf %lf", &x, &y, &z);
            normals[vn] = point(x, y, -z);
            vn++;
        } else if (line.rfind("f ", 0) == 0){
            num_vertices = count_vertices(line);
            //the following code if for splitting polygons into triangles 
            start_index = line.find_first_not_of(" ", 1);//begining of first vertex
            end_index = line.find(" ", start_index);//end of first vertex
            first_vertex.assign(line, 0, end_index);// all triangles will share the first vertex
            //printf("\nline: [%s]\n", line.c_str());
            //printf("fver: [%s]\n", first_vertex.c_str());
            
            for(int i = 0; i < num_vertices - 2; i++){
                //if(f > num_faces) std::cout << "num faces Exceded!\n";
                face_string = first_vertex;
                for(int j=0; j<2; j++){
                    start_index = line.find_first_not_of(" ", end_index);
                    end_index = line.find(" ", start_index + 1);
                    if(end_index == std::string::npos){
                        end_index=line.length();
                    }
                    face_string += " " + line.substr(start_index, end_index-start_index);
                }
                end_index = start_index;//go back so the last vertex is repeated
                //printf("flin: [%s]\n", face_string.c_str());
                faces[f] = buildFace(face_string);
                f++;
            }
        }
        //if( f > 5) break;
    }

    //print the materials in this object
    std::cout << mtl_list.size() << " [";
    for (auto const &i: mtl_list) {
        std::cout << i << ", ";
    } std::cout << "]\n";


    mesh->bvh.set_build_method(BuildMethod::MidSplit);
    mesh->bvh.set_search_method(SearchMethod::BFS);
    mesh->bvh.build(faces, num_faces);

    mesh_obj.close();

    free(vertices);
    free(texture_coords);
    free(normals);
    free(faces);

    mesh->invert_and_bound();
    object_list.push_front(mesh);
}

void MeshFactory::setMaterial(const std::string &mtllib_line){
    size_t start_index = mtllib_line.find_first_of(" ") + 1;
    std::string file_name = dir + mtllib_line.substr(start_index);
    
    if(std::find(mtl_file_list.begin(), mtl_file_list.end(), file_name) == mtl_file_list.end()){
        mtl_file_list.push_front(file_name);
        std::cout << "\n" << file_name << "\n";
        std::cout << "###########################################\n";

        std::string line;
        std::ifstream mesh_obj(file_name);
        if (!mesh_obj.is_open()) {
            std::cout << "Unable to open material file: " << file_name << "\n";
            return;
        }

        while (getline(mesh_obj, line)) {
            if(line.find('\r') != std::string::npos){
                //remove '\r' from windows format objs
                line.erase(line.find('\r'));
            }
            //remove leading whitespacce
            start_index = line.find_first_not_of(" ");
            line = (start_index == std::string::npos) ? "" : line.substr(start_index);

            //std::cout << line << "\n";
            std::string mtl_name;
            std::string texture_name;
            if(line.rfind("newmtl ", 0) == 0){
                start_index = line.find_first_of(" ") + 1;
                mtl_name = file_name + "::" + line.substr(start_index);
                if(std::find(mtl_list.begin(), mtl_list.end(), mtl_name) == mtl_list.end()){
                    mtl_list.push_front(mtl_name);
                }
                //std::cout << "mtl: " << mtl_name << "\n";
            } else if (line.rfind("map_Kd ", 0) == 0){
                start_index = line.find_first_of(" ") + 1;
                texture_name = dir + line.substr(start_index);
                for (std::string::size_type i = 0; i < texture_name.size(); i++) {
                    texture_name[i] = (texture_name[i] == '\\') ? '/' : texture_name[i];
                }
                std::cout << "texture: " << texture_name << std::endl;
                //loadTexture(this, texture_name, 1, scene);
            }
        }
    }
}

TriangleFace* MeshFactory::buildFace(std::string& line){
    //std::cout << "["<< line << "]\n";
    int p1, p2, p3, t1, t2, t3, n1, n2, n3;
    if (sscanf(line.c_str(), "f %d %d %d", &p1, &p2, &p3) == 3){
        adjust_indexes(p1, p2, p3, num_vertices);
        return new TriangleFace(vertices[p3], vertices[p2], vertices[p1]);
    } else if(sscanf(line.c_str(), "f %d/%d %d/%d %d/%d", &p1, &t1, &p2, &t2, &p3, &t3) == 6){
        adjust_indexes(p1, p2, p3, num_vertices);
        adjust_indexes(t1, t2, t3, num_texture_coords);
        TriangleFace_T *face = new TriangleFace_T(vertices[p1], vertices[p2], vertices[p3]);
        face->set_texture_coordinates(texture_coords[t1 * 2],texture_coords[t1 * 2 + 1],
                                      texture_coords[t2 * 2],texture_coords[t2 * 2 + 1],
                                      texture_coords[t3 * 2],texture_coords[t3 * 2 + 1]);
        return face;
    } else if(sscanf(line.c_str(), "f %d//%d %d//%d %d//%d", &p1, &n1, &p2, &n2, &p3, &n3) == 6){
        adjust_indexes(p1, p2, p3, num_vertices);
        adjust_indexes(n1, n2, n3, num_normals);
        TriangleFace_N * face = new TriangleFace_N(vertices[p1], vertices[p2], vertices[p3]);
        face->setNormals(normals[n1], normals[n2], normals[n3]);
        return face;
    } else if(sscanf(line.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d", &p1, &t1, &n1, &p2, &t2, &n2, &p3, &t3, &n3) == 9){
        adjust_indexes(p1, p2, p3, num_vertices);
        adjust_indexes(t1, t2, t3, num_texture_coords);
        adjust_indexes(n1, n2, n3, num_normals);
        TriangleFace_T_N *face = new TriangleFace_T_N(vertices[p1], vertices[p2], vertices[p3]);
        face->setNormals(normals[n1], normals[n2], normals[n3]);
        face->set_texture_coordinates(texture_coords[t1 * 2],texture_coords[t1 * 2 + 1],
                                      texture_coords[t2 * 2],texture_coords[t2 * 2 + 1],
                                      texture_coords[t3 * 2],texture_coords[t3 * 2 + 1]);
        return face;
    }
    return NULL;
}