#include "MeshFactory.h"

#include <fstream>
#include <iostream>
#include <algorithm>

#include "mappings.h"

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

void MeshFactory::setTransform(matrix m){
    transformation = m;
}

void MeshFactory::loadMeshFile(const std::string& filename){
    //mesh = new Mesh(default_color);
    //mesh->name = filename;
    //mesh->T = transformation;
    //mesh->set_pathTrace_properties(1.0, 0.0, 0.0);
    //mesh->r_index = 1.54;

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

    num_vertices = num_texture_coords = num_normals = num_faces = 0;

    int num_objects=0;
    //int* object_face_counts = NULL;
    material* current_material = NULL;

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
            loadMaterialFile(line);
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
        } else if (line.rfind("f ", 0) == 0) {
            //count faces
            num_faces += count_vertices(line) - 2;
        } else if (line.rfind("usemtl ", 0) == 0){
            num_objects++;
        }
    }

    //std::cout << "num vertices: " << num_vertices << "\n";
    //std::cout << "num normals: " << num_normals << "\n";
    //std::cout << "num faces: " << num_faces << "\n";
    //std::cout << "num objects: " << num_objects << "\n";

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

    num_faces_in_object = 0, first_face_in_object = 0;
    int start_index = 2, end_index;
    std::string face_string, first_vertex;
    while (getline(mesh_obj, line)) {
        if(line.find('\r') != std::string::npos){
            //remove '\r' from windows format objs
            line.erase(line.find('\r'));
        }
        
        if (line.rfind("usemtl ", 0) == 0){
            //build the previous object if it isnt null and has more than 0 faces.
            
            if(!object_name.empty() && num_faces_in_object > 0){
                buildMesh();
                object_list.push_front(mesh);
            }
            start_index = line.find_first_not_of(" ", 6);
            object_name = material_file_name + "::" + line.substr(start_index);
            mtl_name = material_file_name + "::" + line.substr(start_index);
            num_faces_in_object = 0;
        } else if (line.rfind("v ", 0) == 0) {
            sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z);
            x = (x - avg_x) / scale, y = (y - avg_y) / scale,
            z = -(z - avg_z) / scale;
            vertices[v] = point(x, y, z);
            v++;
        } else if (line.rfind("vt ", 0) == 0){
            sscanf(line.c_str(), "vt %lf %lf", &x, &y);
            texture_coords[t *2] = x;
            // obj define texture coords in y axis to start from the bottom
            // instead of the top so 1 - y makes textures right side up
            texture_coords[(t * 2) + 1] = 1 - y;
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
                
                // determine material of certain face
                //if (f == 1232600){
                //    std::cout << "mtl: " << mtl_name << "\n";
                //}
                
                faces[f] = buildFace(face_string);
                f++;
                num_faces_in_object++;
            }
        }
        //if( f > 5) break;
    }

    if(num_faces_in_object > 0){
        buildMesh();
        object_list.push_front(mesh);
    }

    //print the materials in this object
    //std::cout << mtl_list.size() << "[\n";
    //for (auto const &i: mtl_list) {
    //    std::cout << i << "\n";
    //} std::cout << "]\n";

    //mesh->bvh.set_build_method(BuildMethod::MidSplit);
    //mesh->bvh.set_search_method(SearchMethod::BFS);
    //mesh->bvh.build(faces, num_faces);

    mesh_obj.close();

    free(vertices);
    free(texture_coords);
    free(normals);
    free(faces);

    //mesh->invert_and_bound();
    //object_list.push_front(mesh);
}

void MeshFactory::loadMaterialFile(const std::string &mtllib_line){
    size_t start_index = mtllib_line.find_first_of(" ") + 1;
    material_file_name = dir + mtllib_line.substr(start_index);
    
    if(std::find(mtl_file_list.begin(), mtl_file_list.end(), material_file_name) == mtl_file_list.end()){
        mtl_file_list.push_front(material_file_name);
        //std::cout << "\n" << material_file_name << "\n";
        //std::cout << "###########################################\n";

        std::string line;
        std::ifstream mesh_obj(material_file_name);
        if (!mesh_obj.is_open()) {
            std::cout << "Unable to open material file: " << material_file_name << "\n";
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
                mtl_name = material_file_name + "::" + line.substr(start_index);
                if(std::find(mtl_list.begin(), mtl_list.end(), mtl_name) == mtl_list.end()){
                    mtl_list.push_front(material(mtl_name));
                }
                //std::cout << "mtl: " << mtl_name << "\n";
            } else if(line.rfind("Ke ", 0) == 0){
                //Ke is a light source in Physically based rendering
                double r,g,b;
                sscanf(line.c_str(), "Ke %lf %lf %lf", &r, &g, &b);
                mtl_list.front().col_diffuse = {r,g,b};
                mtl_list.front().is_light_source = 1;
            } else if(line.rfind("Ka ", 0) == 0){
                double r,g,b;
                sscanf(line.c_str(), "Ka %lf %lf %lf", &r, &g, &b);
                if(!mtl_list.front().is_light_source){
                    mtl_list.front().col_ambient = {r,g,b};
                }
            } else if(line.rfind("Kd ", 0) == 0){
                double r,g,b;
                sscanf(line.c_str(), "Kd %lf %lf %lf", &r, &g, &b);
                if(!mtl_list.front().is_light_source){
                    mtl_list.front().col_diffuse = {r,g,b};
                }
            } else if(line.rfind("Ks ", 0) == 0){
                double r,g,b;
                sscanf(line.c_str(), "Ks %lf %lf %lf", &r, &g, &b);
                if(!mtl_list.front().is_light_source){
                    mtl_list.front().col_specular = {r,g,b};
                }
            } else if (line.rfind("Ns ", 0) == 0){
                sscanf(line.c_str(), "Ns %i", &mtl_list.front().Ns);
            } else if (line.rfind("d ", 0) == 0){
                //d is how opaque an object is. 1 means fully opaque
                sscanf(line.c_str(), "d %lf", &mtl_list.front().alpha);
            } else if (line.rfind("Ti ", 0) == 0){
                //Td is how transparent an object is. 1 means fully transparent
                sscanf(line.c_str(), "Ti %lf", &mtl_list.front().alpha);
                //alpha 1 means fully opaque
                mtl_list.front().alpha = 1 - mtl_list.front().alpha;
            } else if (line.rfind("Ni ", 0) == 0){
                sscanf(line.c_str(), "Ni %lf", &mtl_list.front().index_of_refraction);
            } else if (line.rfind("map_Kd ", 0) == 0){
                start_index = line.find_first_of(" ") + 1;
                texture_name = dir + line.substr(start_index);
                for (std::string::size_type i = 0; i < texture_name.size(); i++) {
                    texture_name[i] = (texture_name[i] == '\\') ? '/' : texture_name[i];
                }
                //std::cout << "texture: [" << texture_name << "]\n";
                //std::cout << "materials: [";
                //for (material m: mtl_list){
                //    std::cout << m << ",";
                //}std::cout << "]\n";
                textureNode *texture = loadTexture(texture_name, 1, texture_list);
                if(texture){
                    mtl_list.front().im = texture->im;
                }
                //mesh->texImg = mtl_list.front().im; 
                
                //set_texture(mesh, texture_name, 1, texture_list);
                //std::cout << "mtl: "<< mtl_list.front().name << "\n";
            }
        }
    }
}

void MeshFactory::buildMesh(){
    // std::cout << "mtl name: " << mtl_name << "\n";
    material mtl = *std::find(mtl_list.begin(), mtl_list.end(), mtl_name);
    color col = mtl.col_ambient + mtl.col_diffuse + mtl.col_specular;
    double diffuse, reflect, refract, refl_sig;
    diffuse = sqrt(mtl.col_ambient.R * mtl.col_ambient.R +
                   mtl.col_ambient.G * mtl.col_ambient.G +
                   mtl.col_ambient.B * mtl.col_ambient.B) +
              sqrt(mtl.col_diffuse.R * mtl.col_diffuse.R +
                   mtl.col_diffuse.G * mtl.col_diffuse.G +
                   mtl.col_diffuse.B * mtl.col_diffuse.B);
    double color_length = sqrt(col.R * col.R + col.G * col.G + col.B * col.B);
    
    diffuse = diffuse / color_length;
    reflect = 1 - diffuse;
    
    double specular_length = sqrt(mtl.col_specular.R * mtl.col_specular.R +
                                  mtl.col_specular.G * mtl.col_specular.G +
                                  mtl.col_specular.B * mtl.col_specular.B);
    refl_sig = MAX_REFL_SIG;
    if (specular_length > 0) {
        if( 0 <= mtl.Ns && mtl.Ns < 1000){
            // convert 0 - 1000 range to 0 to 1 where 1000 maps to 0
            refl_sig -= (MAX_REFL_SIG * (mtl.Ns/1000.0));
        } else if (1000 <= mtl.Ns){
            refl_sig = 0;
        }
    }
    
    diffuse = mtl.alpha * diffuse;
    reflect = mtl.alpha * reflect;
    refract = 1 - mtl.alpha;

    double sum = diffuse + reflect + refract;
    diffuse /= sum;
    reflect /= sum;
    refract /= sum;

    double max_col = max(col.R, max(col.G, col.B));
    if (max_col > 1 && !mtl.is_light_source) {
        col.R /= max_col;
        col.G /= max_col;
        col.B /= max_col;
    }

    if (mtl.is_light_source){
        //std::cout << "col: " << col << "\n";
        MeshLight* ml = new MeshLight(col); 
        ml->isLightSource = mtl.is_light_source;
        ml->pt.LSweight = 10;
        ml->bvh.set_build_method(BuildMethod::MidSplit);
        ml->bvh.set_search_method(SearchMethod::DFS);
        ml->bvh.build(faces + first_face_in_object, num_faces_in_object);
        ml->buildLightFaceList();
        mesh = ml;
    } else {
        mesh = new Mesh(col);
        mesh->bvh.set_build_method(BuildMethod::MidSplit);
        mesh->bvh.set_search_method(SearchMethod::DFS);
        mesh->bvh.build(faces + first_face_in_object, num_faces_in_object);
    }
    first_face_in_object += num_faces_in_object;
    
    //mesh->set_rayTrace_properties(0.1, double diffuse,
    //                                 double specular, double global,
    //                                 double alpha, double shiny)
    mesh->set_pathTrace_properties(diffuse, reflect, refract);
    mesh->refl_sig = refl_sig;
    mesh->name = object_name;
    mesh->T = transformation;
    mesh->r_index = mtl.index_of_refraction;

    /*
    std::cout << "\n------------------------\n";
    std::cout << "mtl: " << mtl.name << "\n";
    std::cout << "color: " << mesh->col << "\n";
    std::cout << "diffuse: " << mesh->pt.diffuse << "\n";
    std::cout << "reflect: " << mesh->pt.reflect << "\n";
    std::cout << "refract: " << mesh->pt.refract << "\n";
    std::cout << "refl_sig: " << mesh->refl_sig << "\n";
    std::cout << "r_index: "<< mesh->r_index << "\n";
    */
    
    mesh->texImg = mtl.im;
    mesh->invert_and_bound();
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