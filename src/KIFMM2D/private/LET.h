//
// Created by lurker on 8/11/16.
//

#ifndef KIFMM2D_LET_H
#define KIFMM2D_LET_H

#include "LET_Node.h"

class LET {
public:
    LET();

    virtual ~LET();


    LET_Node* root;
    ulong N;
    index_t m;
    index_t nSurface;

    vector<Point> standardUpEquivalentSurface;
    vector<Point> standardUpCheckSurface;
    vector<Point> standardDownEquivalentSurface;
    vector<Point> standardDownCheckSurface;



    index_t rank;

    index_t maxLevel;
    MatrixXd chargeTree;
    vector<Point> locationTree;

    Point center;
    Point radius;

    /*
     * utilities variables for timing
     */
    struct timespec start, finish;
    double elapsed = 0;
    /*
     * caching matrix for GMRES part
     */
    Option option;
    size_t cacheIndex;
    vector<MatrixXd> cache;

    void tic() {clock_gettime(CLOCK_MONOTONIC, &start);}
    void toc() {clock_gettime(CLOCK_MONOTONIC, &finish);elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;}

    void initialize(const index_t nSurface, const vector<Point>& location, double* const charge, const ulong N, const index_t rank, Option option);
    void getCenterRadius(const vector<Point>& location, Point& center, Point& radius);
    void maxAndMinCoordinates(const vector<Point>& vec, scalar_t& maxX, scalar_t& maxY, scalar_t& minX, scalar_t& minY);

    void getStandardSurface(const index_t nSurface, scalar_t radius, vector<Point>&  surface);
    void getLocalSurface(const Identity id, Point& center, Point& radius, vector<Point>& surface);

    void getUpwardEquivalent_acc(index_t nSurface, vector<Point> &location, vector<Point> &upwardCheckSurface,  VectorXd &charge, VectorXd &density);
    void getUpwardEquivalent_acc_cache(index_t nSurface, vector<Point> &location, vector<Point> &upwardCheckSurface,  VectorXd &charge, VectorXd &density);
    void getUpwardEquivalent_acc_fast(index_t nSurface, vector<Point> &location, vector<Point> &upwardCheckSurface,  VectorXd &charge, VectorXd &density);



    void getUpwardEquivalent_inv(index_t nSurface,vector<Point> &upwardEquivalentSurface,
                                 vector<Point> &upwardCheckSurface, VectorXd &rhs, VectorXd &density);
    void getUpwardEquivalent_inv_cache(index_t nSurface,vector<Point> &upwardEquivalentSurface,
                                     vector<Point> &upwardCheckSurface, VectorXd &rhs, VectorXd &density);
    void getUpwardEquivalent_inv_fast(index_t nSurface,vector<Point> &upwardEquivalentSurface,
                                     vector<Point> &upwardCheckSurface, VectorXd &rhs, VectorXd &density);


    void getDownwardEquivalent_acc(index_t nSurface, vector<Point>& location, vector<Point>& downCheckSurface,
                               VectorXd& charge, VectorXd& density);
    void getDownwardEquivalent_acc_cache(index_t nSurface, vector<Point>& location, vector<Point>& downCheckSurface,
                               VectorXd& charge, VectorXd& density);
    void getDownwardEquivalent_acc_fast(index_t nSurface, vector<Point>& location, vector<Point>& downCheckSurface,
                               VectorXd& charge, VectorXd& density);

    void getDownwardEquivalent_inv(index_t nSurface, vector<Point> &downwardEquivalentSurface, vector<Point> &downCheckSurface,
                                   VectorXd &rhs, VectorXd &density);
    void getDownwardEquivalent_inv_cache(index_t nSurface, vector<Point> &downwardEquivalentSurface, vector<Point> &downCheckSurface,
                                   VectorXd &rhs, VectorXd &density);
    void getDownwardEquivalent_inv_fast(index_t nSurface, vector<Point> &downwardEquivalentSurface, vector<Point> &downCheckSurface,
                                   VectorXd &rhs, VectorXd &density);


    void interactionList(LET_Node*&node);
    void parentToChild(LET_Node*& node);

    void resetPotential(LET_Node*& node);

    void transferPotential(LET_Node*& node, VectorXd& potential);

    void downPass(LET_Node*& node, VectorXd& potential);

    void getCharge(LET_Node*& node);

    void kernelEvaluate(vector<Point>& targetLocation, vector<Point>& sourceLocation, MatrixXd& K);

    void assignChildren(LET_Node*& node);
    void assignSiblings(LET_Node*& node);
    void assignCousin(LET_Node*& node, index_t neighborIndex);
    void buildTree(LET_Node*& node);


    void display(LET_Node*& node);

    /*
     * pseudo inverse using SVD. lhs = K\rhs.
     */
    void pseudoInverse( MatrixXd& K, VectorXd& rhs, VectorXd& lhs);

    virtual scalar_t kernelFunction(Point& a, Point& b) = 0;


};


#endif //KIFMM2D_LET_H
