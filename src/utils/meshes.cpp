#include "meshes.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <algorithm>

#include "ray.h"
#include "random.h"

int count_vertices(std::string &face){
    int n = 0;
    for(int i = 1; i < face.size(); i++){
        if(face[i] != ' ')
            n++;
        while( i < face.size() && face[i] != ' ')
          i++;
    }
    return n;
}

TriangleFace::TriangleFace() {}

TriangleFace::TriangleFace(point p1, point p2, point p3) {
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
}

TriangleFace::TriangleFace(double x1, double y1, double z1, double x2,
                           double y2, double z2, double x3, double y3,
                           double z3) {
    p1.x = x1, p1.y = y1, p1.z = z1;
    p2.x = x2, p2.y = y2, p2.z = z2;
    p3.x = x3, p3.y = y3, p3.z = z3;
}

double TriangleFace::min_x() const { return MIN(p1.x, MIN(p2.x, p3.x)); }
double TriangleFace::min_y() const { return MIN(p1.y, MIN(p2.y, p3.y)); }
double TriangleFace::min_z() const { return MIN(p1.z, MIN(p2.z, p3.z)); }
double TriangleFace::max_x() const { return MAX(p1.x, MAX(p2.x, p3.x)); }
double TriangleFace::max_y() const { return MAX(p1.y, MAX(p2.y, p3.y)); }
double TriangleFace::max_z() const { return MAX(p1.z, MAX(p2.z, p3.z)); }

void TriangleFace::intersect(struct Ray *ray, double *lambda,
                             point *bary_coords) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified triangle.

    point e12 = p2 - p1;
    point e23 = p3 - p2;

    point normal = cross(&e12, &e23);
    double denom = dot(&normal, &normal);

    // ray plane intersection calculation
    point pi, r = p1 - ray->p0;
    double d_dot_n = dot(&ray->d, &normal);
    double r_dot_n = dot(&r, &normal);
    double t_lambda = r_dot_n / d_dot_n;

    if (t_lambda < THR)  //|| *lambda < t_lambda)
        return;          // intersection with plane is negative

    double u, v, w;

    rayPosition(ray, t_lambda, &pi);

    // e12 has already been calculated
    point v1i = pi - p1;
    point crossp = cross(&e12, &v1i);
    u = dot(&crossp, &normal);
    if (u < 0) return;  // outside first edge

    // e23 has already been calculated
    point v2i = pi - p2;
    crossp = cross(&e23, &v2i);
    v = dot(&crossp, &normal);
    if (v < 0) return;  // outside second edge

    point e31 = p1 - p3;
    point v3i = pi - p3;
    crossp = cross(&e31, &v3i);
    w = dot(&crossp, &normal);
    if (w < 0) return;  // outside third edge

    // intersection is within triangle
    bary_coords->z = u / denom;
    bary_coords->x = v / denom;
    bary_coords->y = w / denom;
    *lambda = t_lambda;
    return;
}

double TriangleFace::intersect(struct Ray *ray, double lambda) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified triangle.

    point e12 = p2 - p1;
    point e23 = p3 - p2;

    point normal = cross(&e12, &e23);

    // ray plane intersection calculation
    point pi, r = p1 - ray->p0;
    double d_dot_n = dot(&ray->d, &normal);
    double r_dot_n = dot(&r, &normal);
    double t_lambda = r_dot_n / d_dot_n;

    if (t_lambda < THR || lambda < t_lambda)
        return INFINITY;  // intersection with plane is negative
                          // or larger than current lambda

    double u, v, w;

    rayPosition(ray, t_lambda, &pi);

    // e12 has already been calculated
    point v1i = pi - p1;
    point crossp = cross(&e12, &v1i);
    u = dot(&crossp, &normal);
    if (u < 0) return INFINITY;  // outside first edge

    // e23 has already been calculated
    point v2i = pi - p2;
    crossp = cross(&e23, &v2i);
    v = dot(&crossp, &normal);
    if (v < 0) return INFINITY;  // outside second edge

    point e31 = p1 - p3;
    point v3i = pi - p3;
    crossp = cross(&e31, &v3i);
    w = dot(&crossp, &normal);
    if (w < 0) return INFINITY;  // outside third edge

    return t_lambda;
}

bool TriangleFace::isprim() { return true; }

point TriangleFace::normal(point *bary_coords) {
    point e12 = p2 - p1;
    point e23 = p3 - p2;

    return cross(&e12, &e23);
}

void TriangleFace::texture_coordinates(double *a, double *b, point *bary_coords) {return;}

double TriangleFace::surfaceArea(){
    point e12 = p2 - p1;
    point e23 = p3 - p2;

    point normal = cross(&e12, &e23);
    return normal.length()/2.0;
}

void TriangleFace::surfaceCoordinates(double a, double b, point& p){
    //https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle-in-3d
    double sqrt_a = sqrt(a);
    p = (1 - sqrt_a) * p1 + (sqrt_a * (1 - b)) * p2 + (b * sqrt_a) * p3;
}

std::ostream &operator<<(std::ostream &strm, const TriangleFace &a) {
    return strm << "p1: " << a.p1 << " p2: " << a.p2 << " p3: " << a.p3;
}

void TriangleFace_N::setNormals(point n1, point n2, point n3) {
    this->n1 = n1, this->n2 = n2, this->n3 = n3;
}

point TriangleFace_N::normal(point *bary_coords) {
    return n1 * bary_coords->x + n2 * bary_coords->y + n3 * bary_coords->z;
}

void TriangleFace_T::set_texture_coordinates(double p1_u, double p1_v, double p2_u, double p2_v, double p3_u, double p3_v){
    this->p1_u=p1_u, this->p1_v=p1_v,this->p2_u=p2_u, this->p2_v=p2_v, this->p3_u=p3_u, this->p3_v=p3_v;
}

void TriangleFace_T::texture_coordinates(double *a, double *b, point *bary_coords){
    *a = p1_u * bary_coords->x + p2_u * bary_coords->y + p3_u * bary_coords->z;
    *b = p1_v * bary_coords->x + p2_v * bary_coords->y + p3_v * bary_coords->z;
}

std::ostream &operator<<(std::ostream &strm, const TriangleFace_N &a) {
    return strm << "p1: " << a.p1 << " p2: " << a.p2 << " p3: " << a.p3
                << "n1: " << a.n1 << " n2: " << a.n2 << " n3: " << a.n3;
}

void BoundingBox::setBounds(double min_x, double max_x, double min_y,
                            double max_y, double min_z, double max_z) {
    this->b_min_x = min_x, this->b_max_x = max_x;
    this->b_min_y = min_y, this->b_max_y = max_y;
    this->b_min_z = min_z, this->b_max_z = max_z;
}

void Mesh::setMesh(const std::string& filename) {
    //std::cout << "\nBuilding Mesh: " << this->label << "\n";
    frontAndBack = 0;

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
    //std::cout << mtl_list.size() << " [";
    //for (auto const &i: mtl_list) {
    //    std::cout << i << ", ";
    //} std::cout << "]\n";


    bvh.set_build_method(BuildMethod::MidSplit);
    bvh.set_search_method(SearchMethod::BFS);
    bvh.build(faces, num_faces);

    mesh_obj.close();

    free(vertices);
    free(texture_coords);
    free(normals);
    free(faces);
    //std::cout<< "Build complete!\n";
}

void Mesh::setMaterial(const std::string &mtllib_line){
    size_t start_index = mtllib_line.find_first_of(" ") + 1;
    std::string file_name = dir + mtllib_line.substr(start_index);
    
    if(std::find(mtl_file_list.begin(), mtl_file_list.end(), file_name) == mtl_file_list.end()){
        mtl_file_list.push_front(file_name);
        //std::cout << "\n" << file_name << "\n";
        //std::cout << "###########################################\n";

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
                //std::cout << "texture: " << texture_name << std::endl;
                //loadTexture(this, texture_name, 1, scene);
            }
        }
    }
}

TriangleFace* Mesh::buildFace(std::string& line){
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

void Mesh::set_canonical_bounds() {
    w_bound.min.x = bvh.root->min_x();
    w_bound.min.y = bvh.root->min_y();
    w_bound.min.z = bvh.root->min_z();
    w_bound.max.x = bvh.root->max_x();
    w_bound.max.y = bvh.root->max_y();
    w_bound.max.z = bvh.root->max_z();
}

void Mesh::intersect(struct Ray *ray, double *lambda, struct point *p,
                     struct point *n, double *a, double *b) {
    *lambda = INFINITY;

    struct Ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    TriangleFace *closest_face;
    //((bvh).*(bvh.search))(ray);
    closest_face = (TriangleFace *)bvh.search(&ray_transformed);

    point bary_coords;
    if (closest_face != NULL) {
        closest_face->intersect(&ray_transformed, lambda, &bary_coords);
        *n = closest_face->normal(&bary_coords);
        normalTransform(n, n, this);
        normalize(n);
        closest_face->texture_coordinates(a, b, &bary_coords);
        rayPosition(ray, *lambda, p);
    }
}

void MeshLight::randomPoint(double *x, double *y, double *z){
    point p;

    TriangleFace* face = NULL;
    
    double a = xor128();
    double b = xor128();
    
    double prob = 0;
    for (MeshLight::LightFace& light : lightFaceList) {
        prob += light.faceWeight;
        if (a < prob) {
            face = light.face;
            break;
        }
    }
    face->surfaceCoordinates(a,b, p);
    p = T * p;
    *x = p.x;
    *y = p.y;
    *z = p.z;
}

void MeshLight::buildLightFaceList(){
    traverseBVH(bvh.root, lightFaceList);
    std::cout << "lightFaceList size: " << lightFaceList.size() << "\n";
    double total_area = 0;
    for(auto it = lightFaceList.begin(); it != lightFaceList.end(); it++){
        total_area += it->faceWeight;
    }

    for(auto it = lightFaceList.begin(); it != lightFaceList.end(); it++){
        it->faceWeight /= total_area;
        std::cout << "lightFace Surface Area: " << it->faceWeight << "\n";
    }
}

void MeshLight::traverseBVH(Primitive* bvhNode, std::list<LightFace>& lightFaceList){
    if(bvhNode == NULL) return;
    if(bvhNode->isprim()){
        TriangleFace* face = (TriangleFace* )bvhNode;
        lightFaceList.push_back({face->surfaceArea(), face});
        return;
    }
    BoundingBox *currentBox = (BoundingBox* )bvhNode;
    traverseBVH(currentBox->c1, lightFaceList);
    traverseBVH(currentBox->c2, lightFaceList);
}

void adjust_indexes(int& i1, int& i2, int& i3, int array_size){
    if (i1 < 0) i1 = array_size + i1;
    if (i2 < 0) i2 = array_size + i2;
    if (i3 < 0) i3 = array_size + i3;
}
