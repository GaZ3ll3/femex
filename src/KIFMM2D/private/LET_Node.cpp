//
// Created by lurker on 8/11/16.
//

#include "LET_Node.h"

Point Point::operator+(const Point &rhs) const {
    Point result;
    result.x = this->x + rhs.x;
    result.y = this->y + rhs.y;
    return result;
}

Point Point::operator*(scalar_t rhs) {
    Point result;
    result.x = this->x * rhs;
    result.y = this->y * rhs;
    return result;
}

LET_Node::LET_Node(index_t nLevel, index_t nodeIndex) {
    parent = nullptr;
    for (index_t k = 0; k < 4; ++k) {child[k] = nullptr;}
    for (index_t k = 0; k < 8; ++k) {neighbor[k] = nullptr;}
    for (index_t k = 0; k < 27; ++k) {interaction[k] = nullptr;}

    nNeighbor = 0;nInteraction = 0;
    isLeaf = false; isEmpty = false;
    this->nLevel = nLevel;
    this->nodeIndex = nodeIndex;
    this->chargeComputed = false;
}

LET_Node::~LET_Node() {
    for (index_t k = 0; k < 4; ++k) {
        if (child[k] != nullptr) {
            delete child[k];
            child[k] = nullptr;
        }
    }
}