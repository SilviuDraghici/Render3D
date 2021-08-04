#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <iostream>

#include "imageProcessor.h"

struct material {
    std::string name;
    color col;
    image *im = NULL;
    
    material(std::string& name):
    name(name)
    {}

    ~material(){
        //std::cout << "deleteing material " << name << "\n";
    }

    bool operator==(const std::string& name){
        return this->name == name;
    }
};

inline std::ostream &operator<<(std::ostream &strm, const material &mat) {
    return strm << mat.name;
}

inline bool operator==(const material* lhs, const std::string& name){
    return lhs->name == name;
}

#endif