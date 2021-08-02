#ifndef TEXTURENODE_H
#define TEXTURENODE_H

#include <string>
#include <iostream>

#include "imageProcessor.h"

struct textureNode {
    std::string name;
    image *im;
    
    bool operator==(const std::string& name){
        return this->name == name;
    }
    
    ~textureNode(){
        std::cout << "deleteing texture " << name << "\n";
        free(im->rgbdata);
        free(im);
    }
};

inline bool operator==(const textureNode* lhs, const std::string& name){
    return lhs->name == name;
}

#endif