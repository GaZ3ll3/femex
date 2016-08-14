//
// Created by lurker on 8/11/16.
//

#ifndef KIFMM2D_LET_NODE_H
#define KIFMM2D_LET_NODE_H


#include "utils.h"

using namespace Eigen;
using std::vector;

class Point {
public:
    scalar_t x;
    scalar_t y;
    Point():x(0.), y(0.){}
    Point(scalar_t x_, scalar_t y_) : x(x_), y(y_) {}
    Point operator+(const Point& rhs) const ;
    Point operator*(scalar_t rhs);
};


class LET_Node {
public:
    LET_Node* parent;
    LET_Node* child[4];
    LET_Node* neighbor[8];
    LET_Node* interaction[27];

    index_t nLevel;
    index_t nNeighbor;
    index_t nInteraction;
    index_t nodeIndex;

    Point center;
    Point radius;

    ulong N;

    VectorXi index;
    vector<Point> location;

    vector<Point> shiftedUpEquivalentSurface;
    vector<Point> shiftedDownEquivalentSurface;
    vector<Point> shiftedUpCheckSurface;
    vector<Point> shiftedDownCheckSurface;

    VectorXd charge;
    VectorXd potential;

    VectorXd upwardEquivalent;
    VectorXd downwardEquivalent;

    bool isLeaf;
    bool isEmpty;
    bool chargeComputed;

    LET_Node(index_t nLevel, index_t nodeIndex);
    ~LET_Node();

};


#endif //KIFMM2D_LET_NODE_H
