//
// Created by lurker on 8/11/16.
//

#include "LET.h"
#include <ctime>


using namespace Eigen;
using namespace std;


LET::LET() {
	root = nullptr;
	rank = 0;
	m = 0;
	nSurface = 0;
	N = 0;
	maxLevel = 0;
	cacheIndex = 0;
	option = Option::Cache;
}

void LET::initialize(const index_t nSurface, const vector<Point> &location, double *const charge, const ulong N, const index_t rank, Option option) {
    this->rank = rank;
    this->N = N;
    this->m = 1;
    this->maxLevel = 0;
    this->nSurface = nSurface;
    this->chargeTree = Map<MatrixXd>(charge, N, size_t(1));
    this->locationTree = location;
    this->option = option;
    /*
     * important for acceleration
     */
    this->cacheIndex = 0;


    standardUpCheckSurface.resize(nSurface);
    standardUpEquivalentSurface.resize(nSurface);
    standardDownCheckSurface.resize(nSurface);
    standardDownEquivalentSurface.resize(nSurface);

    getStandardSurface(nSurface, 3.0, standardUpCheckSurface);
    getStandardSurface(nSurface, 1.0, standardUpEquivalentSurface);
    getStandardSurface(nSurface, 1.0 , standardDownCheckSurface);
    getStandardSurface(nSurface, 3.0, standardDownEquivalentSurface);

    getCenterRadius(location, center, radius);

    maxLevel = 0;

    root = new LET_Node(0, 0);
    root->nNeighbor = 0; root->nInteraction = 0;
    root->N = N;
    root->center = center;
    root->radius = {0.5, 0.5};
    root->index.setLinSpaced(N, 0, N-1);

    assignChildren(root);

    buildTree(root);

}

void LET::assignChildren(LET_Node *&node) {
    if (node->N == 0) {
        node->isLeaf = true;
        node->isEmpty = true;
    }
    else {
        node->potential = VectorXd::Zero(node->N);
        node->upwardEquivalent = VectorXd::Zero(nSurface);
        node->downwardEquivalent = VectorXd::Zero(nSurface);



        getLocalSurface(UE, node->center, node->radius, node->shiftedUpEquivalentSurface);
        getLocalSurface(UC, node->center, node->radius, node->shiftedUpCheckSurface);
        getLocalSurface(DE, node->center, node->radius, node->shiftedDownEquivalentSurface);
        getLocalSurface(DC, node->center, node->radius, node->shiftedDownCheckSurface);


        for (size_t k = 0; k < node->N; ++k) {
            node->location.push_back(locationTree[(size_t)node->index(k)]);
        }



        /*
         * when particles inside current box is beyond 4 * rank, split it into 4 children.
         * else it is a leaf.
         */
        if (node->N < (ulong) 4 * rank) {
            node->isLeaf = true;
            /*
             * 1. at leaf node, calculate upward equivalent density at each upwardEquivalentSurface.
             */
            getCharge(node);

            VectorXd contribution; contribution = VectorXd::Zero(nSurface);

            if (option == Option::Cache) {
				getUpwardEquivalent_acc_cache(nSurface, node->location,
									node->shiftedUpCheckSurface, node->charge, contribution);

				getUpwardEquivalent_inv_cache(nSurface, node->shiftedUpEquivalentSurface,
						node->shiftedUpCheckSurface, contribution, node->upwardEquivalent);
            }
            else {
				getUpwardEquivalent_acc_fast(nSurface, node->location,
									node->shiftedUpCheckSurface, node->charge, contribution);

				getUpwardEquivalent_inv_fast(nSurface, node->shiftedUpEquivalentSurface,
						node->shiftedUpCheckSurface, contribution, node->upwardEquivalent);
            }


            if (this->maxLevel < node->nLevel) {this->maxLevel = node->nLevel;}
        }
        else {
            for (index_t k = 0; k < 4; ++k) {
                node->child[k] = new LET_Node(index_t(node->nLevel + 1), k);
                node->child[k]->parent = node;
                node->child[k]->center.x = node->center.x + ((k & 1) - 0.5) * node->radius.x;
                node->child[k]->center.y = node->center.y + ((k >> 1) - 0.5) * node->radius.y;
                node->child[k]->radius.x = node->radius.x * 0.5;
                node->child[k]->radius.y = node->radius.y * 0.5;
                node->child[k]->N = 0;
            }

            for(size_t k=0; k<node->N; ++k){
                if(locationTree[(size_t)node->index(k)].x<node->center.x){
                    if(locationTree[(size_t)node->index(k)].y<node->center.y){
                        node->child[0]->index.conservativeResize(node->child[0]->N+1);
                        node->child[0]->index(node->child[0]->N)	=	node->index(k);
                        ++node->child[0]->N;
                    }
                    else{
                        node->child[2]->index.conservativeResize(node->child[2]->N+1);
                        node->child[2]->index(node->child[2]->N)	=	node->index(k);
                        ++node->child[2]->N;
                    }
                }
                else{
                    if(locationTree[(size_t)node->index(k)].y<node->center.y){
                        node->child[1]->index.conservativeResize(node->child[1]->N+1);
                        node->child[1]->index(node->child[1]->N)	=	node->index(k);
                        ++node->child[1]->N;
                    }
                    else{
                        node->child[3]->index.conservativeResize(node->child[3]->N+1);
                        node->child[3]->index(node->child[3]->N)	=	node->index(k);
                        ++node->child[3]->N;
                    }
                }
            }


            VectorXd contribution;contribution = VectorXd::Zero(nSurface);
            for (index_t k = 0; k < 4; ++k) {
                assignChildren(node->child[k]);
                // bottom up

                if (!node->child[k]->isEmpty) {
                    /*
                     * 2. each child's upward equivalent density will contribute to parent upward equivalent density.
                     */
                	if (option == Option::Cache) {
                		getUpwardEquivalent_acc_cache(nSurface, node->child[k]->shiftedUpEquivalentSurface,
                                        node->shiftedUpCheckSurface, node->child[k]->upwardEquivalent, contribution);
                	}
                	else {
                		getUpwardEquivalent_acc_fast(nSurface, node->child[k]->shiftedUpEquivalentSurface,
                		                                        node->shiftedUpCheckSurface, node->child[k]->upwardEquivalent, contribution);
                	}

                }
            }
            if (option == Option::Cache) {
            	getUpwardEquivalent_inv_cache(nSurface,  node->shiftedUpEquivalentSurface,
            			node->shiftedUpCheckSurface, contribution, node->upwardEquivalent);
            }
            else {
            	getUpwardEquivalent_inv_fast(nSurface,  node->shiftedUpEquivalentSurface,
            	            			node->shiftedUpCheckSurface, contribution, node->upwardEquivalent);
            }
        }
    }
}


void LET::buildTree(LET_Node *&node) {
    if(!node->isEmpty){
        if(!node->isLeaf){
            assignSiblings(node);
            for(unsigned short k=0;k<8;++k){
                if(node->neighbor[k]!=nullptr){
                    if(!node->neighbor[k]->isLeaf  && !node->neighbor[k]->isEmpty){
                        assignCousin(node,k);
                    }
                }
            }
            for(unsigned short k=0;k<4;++k) {
                buildTree(node->child[k]);
            }
        }
    }
}

void LET::resetPotential(LET_Node *&node) {
    if (node != nullptr) {
        node->potential = VectorXd::Zero(node->potential.rows());
        node->downwardEquivalent = VectorXd::Zero(node->downwardEquivalent.rows());

        for (index_t k = 0; k < 4; ++k) {
            resetPotential(node->child[k]);
        }
    }
}

void LET:: assignSiblings(LET_Node*& node){
//	Assign siblings to child[0]
    node->child[0]->neighbor[3]	=	node->child[1];
    node->child[0]->neighbor[5]	=	node->child[2];
    node->child[0]->neighbor[4]	=	node->child[3];

//	Assign siblings to child[1]
    node->child[1]->neighbor[7]	=	node->child[0];
    node->child[1]->neighbor[6]	=	node->child[2];
    node->child[1]->neighbor[5]	=	node->child[3];

//	Assign siblings to child[2]
    node->child[2]->neighbor[1]	=	node->child[0];
    node->child[2]->neighbor[2]	=	node->child[1];
    node->child[2]->neighbor[3]	=	node->child[3];

//	Assign siblings to child[3]
    node->child[3]->neighbor[0]	=	node->child[0];
    node->child[3]->neighbor[1]	=	node->child[1];
    node->child[3]->neighbor[7]	=	node->child[2];

    for(unsigned short k=0;k<4;++k){
        node->child[k]->nNeighbor+=3;
    }
}



void LET:: assignCousin(LET_Node*& node, index_t neighorIndex){
//	Assigning children of neighbor 0
    if(neighorIndex==0 && node->neighbor[0] != nullptr ){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[0]->child[0];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[0]->child[1];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[0]->child[2];
        node->child[0]->neighbor[0]					=	node->neighbor[0]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[0]->child[0];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[0]->child[1];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[0]->child[2];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[0]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[0]->child[0];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[0]->child[1];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[0]->child[2];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[0]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[0]->child[0];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[0]->child[1];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[0]->child[2];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[0]->child[3];

//	Update neighbor count.

        node->child[0]->nNeighbor					=	node->child[0]->nNeighbor+1;
    }
//	Assigning children of neighbor 1
    else if(neighorIndex==1 && node->neighbor[1] != nullptr ){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[1]->child[0];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[1]->child[1];
        node->child[0]->neighbor[1]					=	node->neighbor[1]->child[2];
        node->child[0]->neighbor[2]					=	node->neighbor[1]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[1]->child[0];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[1]->child[1];
        node->child[1]->neighbor[0]					=	node->neighbor[1]->child[2];
        node->child[1]->neighbor[1]					=	node->neighbor[1]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[1]->child[0];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[1]->child[1];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[1]->child[2];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[1]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[1]->child[0];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[1]->child[1];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[1]->child[2];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[1]->child[3];

//	Update neighbor count.

        node->child[0]->nNeighbor					=	node->child[0]->nNeighbor+2;
        node->child[1]->nNeighbor					=	node->child[1]->nNeighbor+2;
    }
//	Assigning children of neighbor 2
    else if (neighorIndex==2 && node->neighbor[2] != nullptr ){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[2]->child[0];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[2]->child[1];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[2]->child[2];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[2]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[2]->child[0];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[2]->child[1];
        node->child[1]->neighbor[2]					=	node->neighbor[2]->child[2];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[2]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[2]->child[0];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[2]->child[1];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[2]->child[2];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[2]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[2]->child[0];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[2]->child[1];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[2]->child[2];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[2]->child[3];

//	Update neighbor count.

        node->child[1]->nNeighbor					=	node->child[1]->nNeighbor+1;
    }
//	Assigning children of neighbor 3
    else if(neighorIndex==3 && node->neighbor[3] != nullptr ){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[3]->child[0];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[3]->child[1];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[3]->child[2];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[3]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

        node->child[1]->neighbor[3]					=	node->neighbor[3]->child[0];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[3]->child[1];
        node->child[1]->neighbor[4]					=	node->neighbor[3]->child[2];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[3]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[3]->child[0];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[3]->child[1];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[3]->child[2];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[3]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

        node->child[3]->neighbor[2]					=	node->neighbor[3]->child[0];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[3]->child[1];
        node->child[3]->neighbor[3]					=	node->neighbor[3]->child[2];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[3]->child[3];

//	Update neighbor count.

        node->child[1]->nNeighbor					=	node->child[1]->nNeighbor+2;
        node->child[3]->nNeighbor					=	node->child[3]->nNeighbor+2;
    }
//	Assigning children of neighbor 4
    else if(neighorIndex==4 && node->neighbor[4] != nullptr ){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[4]->child[0];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[4]->child[1];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[4]->child[2];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[4]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[4]->child[0];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[4]->child[1];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[4]->child[2];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[4]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[4]->child[0];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[4]->child[1];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[4]->child[2];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[4]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

        node->child[3]->neighbor[4]					=	node->neighbor[4]->child[0];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[4]->child[1];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[4]->child[2];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[4]->child[3];

//	Update neighbor count.

        node->child[3]->nNeighbor					=	node->child[3]->nNeighbor+1;
    }
//	Assigning children of neighbor 5
    else if(neighorIndex==5 && node->neighbor[5] != nullptr ){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[5]->child[0];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[5]->child[1];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[5]->child[2];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[5]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[5]->child[0];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[5]->child[1];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[5]->child[2];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[5]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

        node->child[2]->neighbor[5]					=	node->neighbor[5]->child[0];
        node->child[2]->neighbor[4]					=	node->neighbor[5]->child[1];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[5]->child[2];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[5]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

        node->child[3]->neighbor[6]					=	node->neighbor[5]->child[0];
        node->child[3]->neighbor[5]					=	node->neighbor[5]->child[1];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[5]->child[2];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[5]->child[3];

//	Update neighbor count.

        node->child[2]->nNeighbor					=	node->child[2]->nNeighbor+2;
        node->child[3]->nNeighbor					=	node->child[3]->nNeighbor+2;
    }
//	Assigning children of neighbor 6
    else if (neighorIndex==6 && node->neighbor[6] != nullptr ){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[6]->child[0];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[6]->child[1];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[6]->child[2];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[6]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[6]->child[0];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[6]->child[1];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[6]->child[2];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[6]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[6]->child[0];
        node->child[2]->neighbor[6]					=	node->neighbor[6]->child[1];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[6]->child[2];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[6]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[6]->child[0];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[6]->child[1];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[6]->child[2];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[6]->child[3];

//	Update neighbor count.

        node->child[2]->nNeighbor					=	node->child[2]->nNeighbor+1;
    }
//	Assigning children of neighbor 7
    else if (neighorIndex==7 && node->neighbor[7] != nullptr ){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[7]->child[0];
        node->child[0]->neighbor[7]					=	node->neighbor[7]->child[1];
        node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[7]->child[2];
        node->child[0]->neighbor[6]					=	node->neighbor[7]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[7]->child[0];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[7]->child[1];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[7]->child[2];
        node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[7]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[7]->child[0];
        node->child[2]->neighbor[0]					=	node->neighbor[7]->child[1];
        node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[7]->child[2];
        node->child[2]->neighbor[7]					=	node->neighbor[7]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[7]->child[0];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[7]->child[1];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[7]->child[2];
        node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[7]->child[3];

//	Update neighbor count.

        node->child[0]->nNeighbor					=	node->child[0]->nNeighbor+2;
        node->child[2]->nNeighbor					=	node->child[2]->nNeighbor+2;
    }
}


void LET::getUpwardEquivalent_acc_cache(index_t nSurface, vector<Point> &location, vector<Point> &upwardCheckSurface,
                                  VectorXd &charge, VectorXd &density) {
    /*
     *  1. calculate right hand side's potential for this node(leaf node).
     */

    MatrixXd rhs_K; kernelEvaluate(upwardCheckSurface, location, rhs_K);
    cache.push_back(rhs_K);
    density += rhs_K * charge;
}

void LET::getUpwardEquivalent_acc_fast(index_t nSurface, vector<Point> &location, vector<Point> &upwardCheckSurface,
                                  VectorXd &charge, VectorXd &density) {
    /*
     *  1. calculate right hand side's potential for this node(leaf node).
     */

    density += cache[cacheIndex] * charge;
    (cacheIndex)++;
}



void LET::getUpwardEquivalent_acc(index_t nSurface, vector<Point> &location, vector<Point> &upwardCheckSurface,
                                  VectorXd &charge, VectorXd &density) {
    /*
     *  1. calculate right hand side's potential for this node(leaf node).
     */

    MatrixXd rhs_K; kernelEvaluate(upwardCheckSurface, location, rhs_K);

    density += rhs_K * charge;
}


void LET::getUpwardEquivalent_inv(index_t nSurface,vector<Point> &upwardEquivalentSurface,
                                  vector<Point> &upwardCheckSurface, VectorXd &rhs, VectorXd &density) {
    /*
     *  2. invert equivalent-check matrix
     */
    density = VectorXd::Zero(nSurface);

    MatrixXd equivalentCheck; kernelEvaluate(upwardCheckSurface, upwardEquivalentSurface, equivalentCheck);

    pseudoInverse(equivalentCheck, rhs, density);

}


void LET::getUpwardEquivalent_inv_cache(index_t nSurface,vector<Point> &upwardEquivalentSurface,
                                  vector<Point> &upwardCheckSurface, VectorXd &rhs, VectorXd &density) {
    /*
     *  2. invert equivalent-check matrix
     */
    density = VectorXd::Zero(nSurface);

    MatrixXd equivalentCheck; kernelEvaluate(upwardCheckSurface, upwardEquivalentSurface, equivalentCheck);

    pseudoInverse(equivalentCheck, rhs, density);

}


void LET::getUpwardEquivalent_inv_fast(index_t nSurface,vector<Point> &upwardEquivalentSurface,
                                  vector<Point> &upwardCheckSurface, VectorXd &rhs, VectorXd &density) {
    /*
     *  2. invert equivalent-check matrix
     */
	MatrixXd equivalentCheck;
    pseudoInverse(equivalentCheck, rhs, density);

}


/*
 * Since there is a single inversion applied to all interaction list. Here output should be just local potential.
 */
void LET::getDownwardEquivalent_acc(index_t nSurface, vector<Point> &location, vector<Point> &downCheckSurface, VectorXd &charge, VectorXd &density) {

    /*
     *  for M2L part.
     */
    MatrixXd rhs_K; kernelEvaluate(downCheckSurface, location, rhs_K);
    density += rhs_K * charge;

}

void LET::getDownwardEquivalent_acc_cache(index_t nSurface, vector<Point> &location, vector<Point> &downCheckSurface, VectorXd &charge, VectorXd &density) {

    /*
     *  for M2L part.
     */
    MatrixXd rhs_K; kernelEvaluate(downCheckSurface, location, rhs_K);

    cache.push_back(rhs_K);
    density += rhs_K * charge;

}


void LET::getDownwardEquivalent_acc_fast(index_t nSurface, vector<Point> &location, vector<Point> &downCheckSurface, VectorXd &charge, VectorXd &density) {

    /*
     *  for M2L part.
     */
    density += cache[cacheIndex] * charge;
    (cacheIndex)++;

}


void LET::getDownwardEquivalent_inv(index_t nSurface, vector<Point> &downwardEquivalentSurface, vector<Point> &downCheckSurface,
                                    VectorXd &rhs, VectorXd &density) {

    density = VectorXd::Zero(nSurface);

    MatrixXd equivalentCheck; kernelEvaluate(downCheckSurface, downwardEquivalentSurface, equivalentCheck);

    pseudoInverse(equivalentCheck, rhs, density);

}

void LET::getDownwardEquivalent_inv_cache(index_t nSurface, vector<Point> &downwardEquivalentSurface, vector<Point> &downCheckSurface,
                                    VectorXd &rhs, VectorXd &density) {

    density = VectorXd::Zero(nSurface);

    MatrixXd equivalentCheck; kernelEvaluate(downCheckSurface, downwardEquivalentSurface, equivalentCheck);

    //cache.push_back(equivalentCheck);

    pseudoInverse(equivalentCheck, rhs, density);
}

void LET::getDownwardEquivalent_inv_fast(index_t nSurface, vector<Point> &downwardEquivalentSurface, vector<Point> &downCheckSurface,
                                    VectorXd &rhs, VectorXd &density) {

	MatrixXd equivalentCheck;
    pseudoInverse(equivalentCheck, rhs, density);

}


void LET::downPass(LET_Node *&node, VectorXd &potential) {
    if (!node->isEmpty) {
        if (node->isLeaf) {
            /*
             *  downpass to leaf, accumulate all contributions from U list.
             */
            MatrixXd neighborMatrix;
            for (index_t k = 0; k < 8; ++k) {
                if (node->neighbor[k]!= nullptr) {
                    if (!node->neighbor[k]->isEmpty) {
                        // each time reset matrix.
                    	if (option == Option ::Cache) {
                    		kernelEvaluate(node->location, node->neighbor[k]->location, neighborMatrix);
                    		cache.push_back(neighborMatrix);
                    	}
                    	else {
                    		neighborMatrix = cache[cacheIndex];
                    		cacheIndex++;
                    	}

                        getCharge(node->neighbor[k]);
                        node->potential+= neighborMatrix * node->neighbor[k]->charge;
                    }
                }
            }
            /*
             * from downward equivalent surface to locations.
             */

            MatrixXd K;
            if (option == Option::Cache) {
            	kernelEvaluate(node->location, node->shiftedDownEquivalentSurface, K);
            	cache.push_back(K);
            }
            else {
            	K = cache[cacheIndex];
            	cacheIndex++;
            }
            node->potential += K * node->downwardEquivalent;

            /*
             * self contribution in U list
             */

            if (option == Option::Cache) {
				kernelEvaluate(node->location, node->location, K);
				cache.push_back(K);
            }
            else {
            	K = cache[cacheIndex];
            	cacheIndex++;
            }


            node->potential += K * node->charge;

            transferPotential(node, potential);
        }
        else {

            bool computePotential = false;
            for (index_t k = 0; k < 8; ++k) {
                if (node->neighbor[k] != nullptr) {
                    if (!node->neighbor[k]->isEmpty) {
                        if (node->neighbor[k]->isLeaf) {
                            MatrixXd neighborMatrix;
                            if (option == Option::Cache) {
                            	kernelEvaluate(node->location, node->neighbor[k]->location, neighborMatrix);
                            	cache.push_back(neighborMatrix);
                            }
                            else {
                            	neighborMatrix = cache[cacheIndex];
                            	cacheIndex++;
                            }
                            getCharge(node->neighbor[k]);
                            node->potential += neighborMatrix * node->neighbor[k]->charge;
                            computePotential = true;
                        }
                    }
                }
            }
            /*
             *  V list, M2L
             */

            interactionList(node);

            /*
             *  from parent to child, L2L
             */
            parentToChild(node);

            if (computePotential) {
                transferPotential(node, potential);
            }
            for (index_t k = 0; k < 4; ++k) {
                downPass(node->child[k], potential);
            }
        }
    }
}

void LET::parentToChild(LET_Node *&node) {
    for (index_t k = 0; k < 4; ++k) {
        if (!node->child[k]->isEmpty) {
            VectorXd contribution; contribution = VectorXd::Zero(nSurface);
            VectorXd result; result = VectorXd::Zero(nSurface);
            if (option == Option::Cache) {
				getDownwardEquivalent_acc_cache(nSurface, node->shiftedDownEquivalentSurface, node->child[k]->shiftedDownCheckSurface, node->downwardEquivalent, contribution);
				getDownwardEquivalent_inv_cache(nSurface, node->child[k]->shiftedDownEquivalentSurface, node->child[k]->shiftedDownCheckSurface, contribution, result);
            }
            else {
            	getDownwardEquivalent_acc_fast(nSurface, node->shiftedDownEquivalentSurface, node->child[k]->shiftedDownCheckSurface, node->downwardEquivalent, contribution);
            	getDownwardEquivalent_inv_fast(nSurface, node->child[k]->shiftedDownEquivalentSurface, node->child[k]->shiftedDownCheckSurface, contribution, result);
            }
            node->child[k]->downwardEquivalent += result;
        }

    }
}



void LET::interactionList(LET_Node *&node) {

    for (index_t k = 0 ; k < 4; ++k) {
        if (!node->child[k]->isEmpty) {
            /*
             * interaction list. aka V list.
             */
            VectorXd contribution; contribution = VectorXd::Zero(nSurface);
            for (index_t interactionIndex = 0; interactionIndex < node->child[k]->nInteraction; ++interactionIndex) {

            	if (option == Option::Cache) {
            		getDownwardEquivalent_acc_cache(nSurface, node->child[k]->interaction[interactionIndex]->shiftedUpEquivalentSurface,
                                          node->child[k]->shiftedDownCheckSurface, node->child[k]->interaction[interactionIndex]->upwardEquivalent, contribution);
            	}
            	else {
            		getDownwardEquivalent_acc_fast(nSurface, node->child[k]->interaction[interactionIndex]->shiftedUpEquivalentSurface,
            		                                          node->child[k]->shiftedDownCheckSurface, node->child[k]->interaction[interactionIndex]->upwardEquivalent, contribution);
            	}
            }
            if (option == Option::Cache) {
            	getDownwardEquivalent_inv_cache(nSurface, node->child[k]->shiftedDownEquivalentSurface,
                                      node->child[k]->shiftedDownCheckSurface, contribution, node->child[k]->downwardEquivalent);
            }
            else {
            	getDownwardEquivalent_inv_fast(nSurface, node->child[k]->shiftedDownEquivalentSurface,
                                      node->child[k]->shiftedDownCheckSurface, contribution, node->child[k]->downwardEquivalent);
            }
        }

    }
}


void LET::transferPotential(LET_Node *&node, VectorXd &potential) {
    for (size_t k = 0; k < node->N; ++k) {
        potential.row((size_t)node->index(k)) += node->potential.row(k);
    }
}


void LET::getLocalSurface(const Identity id, Point &center, Point& radius, vector<Point>& surface) {
    surface.resize(nSurface);
    if (id == Identity::UE) {
        for (index_t k = 0; k < nSurface; ++k) {
            surface[k].x = center.x + radius.x * standardUpEquivalentSurface[k].x;
            surface[k].y = center.y + radius.y * standardUpEquivalentSurface[k].y;
        }
    }
    else if (id == Identity::UC) {
        for (index_t k = 0; k < nSurface; ++k) {
            surface[k].x = center.x +  radius.x * standardUpCheckSurface[k].x;
            surface[k].y = center.y +  radius.y * standardUpCheckSurface[k].y;
        }
    }
    else if (id == Identity::DE) {
        for (index_t k = 0; k < nSurface; ++k) {
            surface[k].x = center.x + radius.x * standardDownEquivalentSurface[k].x;
            surface[k].y = center.y + radius.y * standardDownEquivalentSurface[k].y;
        }
    }
    else if (id == Identity::DC) {
        for (index_t k = 0; k < nSurface; ++k) {
            surface[k].x = center.x + radius.x * standardDownCheckSurface[k].x;
            surface[k].y = center.y + radius.y * standardDownCheckSurface[k].y;
        }
    }
    else {
        cout << "Wrong Identity Surface found.\n" << endl;
        exit(0);
    }
}


void LET::getCharge(LET_Node *&node) {
    if (node->chargeComputed == true) {
        return;
    }
    else {
        node->chargeComputed = true;
        node->charge = VectorXd::Zero(node->N);
        for (size_t k = 0; k < node->N; ++k) {
            node->charge.row(k) = chargeTree.row((size_t)node->index(k));
        }
    }

}

void LET::getStandardSurface(const index_t nSurface, scalar_t radius, vector<Point> &surface) {
//    scalar_t theta = 2 * M_PI/nSurface;
//    for (index_t k = 0; k < nSurface; ++k) {
//        surface[k] = {radius * cos(k *theta), radius * sin(k * theta)};
//    }
    index_t each_side = nSurface/4;
    for (index_t k = 0; k < each_side; ++k) {
        surface[k] = {radius, -radius + k * (2 * radius)/each_side};
        surface[k + each_side] = {radius - k * (2 * radius)/each_side, radius};
        surface[k + 2 * each_side] = {-radius, radius - k * (2 * radius)/each_side};
        surface[k + 3 * each_side] = {-radius + k * (2 * radius)/each_side,-radius};
    }

//    index_t each_side = nSurface/4;
//    for (index_t k = 0; k < each_side; ++k) {
//        surface[k] = {radius, -radius * cos(k * 2 * M_PI/each_side)};
//        surface[k + each_side] = {radius * cos(k * 2 * M_PI/each_side), radius};
//        surface[k + 2 * each_side] = {-radius, radius * cos(k * 2 * M_PI/each_side)};
//        surface[k + 3 * each_side] = {-radius * cos(k * 2 * M_PI/each_side),-radius};
//    }

}


void LET::getCenterRadius(const vector<Point>& location, Point& center, Point& radius){
    double maxX;
    double maxY;
    double minX;
    double minY;
    maxAndMinCoordinates(location, maxX, maxY, minX, minY);
    center.x = 0.5*(maxX + minX);
    center.y = 0.5*(maxY + minY);
    radius.x = 0.5*(maxX - minX);
    radius.y = 0.5*(maxY - minY);
}

void LET::maxAndMinCoordinates(const vector<Point> &vec, scalar_t &maxX, scalar_t &maxY, scalar_t &minX,
                                scalar_t &minY) {
    maxX = vec[0].x; maxY = vec[0].y;
    minX = maxX; minY = maxY;
    for (size_t i = 0; i < vec.size(); i++) {
        if (vec[i].x > maxX) { maxX = vec[i].x; }
        if (vec[i].y > maxY) { maxY = vec[i].y; }
        if (vec[i].x < minX) { minX = vec[i].x; }
        if (vec[i].y < minY) { minY = vec[i].y; }
    }
}


/*
 * K (x , y)
 */
void LET::kernelEvaluate(vector<Point> &targetLocation, vector<Point> &sourceLocation, MatrixXd &K) {
    size_t targetSize = targetLocation.size();
    size_t sourceSize = sourceLocation.size();

    K = MatrixXd::Zero(targetSize, sourceSize);

    for (size_t targetIndex = 0; targetIndex < targetSize; ++targetIndex) {
        for (size_t sourceIndex = 0; sourceIndex < sourceSize; ++sourceIndex) {
            K(targetIndex, sourceIndex) = kernelFunction(targetLocation[targetIndex], sourceLocation[sourceIndex]);
        }
    }
}

void LET::pseudoInverse(MatrixXd &K, VectorXd &rhs, VectorXd& lhs) {

	if (option != Option::Cache) {
		lhs = cache[cacheIndex] * rhs;
		cacheIndex++;
		return;
	}

    double tolerance=1.e-8; // choose your tolerance wisely!

    JacobiSVD<MatrixXd> svd(K, ComputeThinU | ComputeThinV);
    MatrixXd u = svd.matrixU();
    MatrixXd v = svd.matrixV();
    VectorXd s = svd.singularValues();

    for (size_t i=0; i< s.rows(); ++i) {
        s(i) = s(i) > tolerance ? 1.0/s(i) : 0.;
    }


    cache.push_back(v*s.asDiagonal()*u.transpose());


    lhs = v*s.asDiagonal()*u.transpose() * rhs;

}

LET::~LET() {
    if(root!=nullptr){
        delete root;
        root = nullptr;
    }
}
