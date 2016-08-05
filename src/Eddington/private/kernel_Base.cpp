/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	kernel_Base.cpp
*/

#include"kernel_Base.hpp"
#include <omp.h>

#define scal  0.2

void kernel_Base::calculate_Potential_cache_svd(H2_2D_Node*& node, MatrixXd& potential,H2_2D_Tree& tree) {
    if(!node->isEmpty){
		if(node->isLeaf){
			MatrixXd tempK;

			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
					kernel_2D(node->N , node->location, node->neighbor[k]->N, node->neighbor[k]->location, tempK);
					// cached, copied.

					int row = int(tempK.rows() * scal);
					JacobiSVD<MatrixXd> svd(tempK, ComputeThinU | ComputeThinV);
					MatrixXd u = svd.matrixU().leftCols(row);
					MatrixXd v = svd.matrixV().leftCols(row);
					VectorXd s = svd.singularValues().head(row);
					cache.push_back(u);
					cache.push_back(v);
					sigma.push_back(s);
					tree.get_Charge(node->neighbor[k]);
					node->potential+=tempK*node->neighbor[k]->charge;

					}
				}
			}
            //			Potential from Chebyshev nodes
			node->potential+=node->R*node->nodePotential;
            //			Self potential
			kernel_2D(node->N , node->location, node->N , node->location, tempK);
			//cached, copied
			const int row = int(tempK.rows() *scal );
			JacobiSVD<MatrixXd> svd(tempK, ComputeThinU | ComputeThinV);
			MatrixXd u = svd.matrixU().leftCols(row);
			MatrixXd v = svd.matrixV().leftCols(row);
			VectorXd s = svd.singularValues().head(row);
			cache.push_back(u);
			cache.push_back(v);
			sigma.push_back(s);


			node->potential+=tempK *node->charge;
			tranfer_Potential_To_Potential_Tree(node, potential);
		}
		else{
			bool computePotential	=	false;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
					if(node->neighbor[k]->isLeaf){
						MatrixXd tempK;
						kernel_2D(node->N, node->location, node->neighbor[k]->N, node->neighbor[k]->location, tempK);
						//cached, copied

						const int row = int(tempK.rows() *scal );
						JacobiSVD<MatrixXd> svd(tempK, ComputeThinU | ComputeThinV);
						MatrixXd u = svd.matrixU().leftCols(row);
						MatrixXd v = svd.matrixV().leftCols(row);
						VectorXd s = svd.singularValues().head(row);
						cache.push_back(u);
						cache.push_back(v);
						sigma.push_back(s);

						tree.get_Charge(node->neighbor[k]);
						node->potential+=tempK *node->neighbor[k]->charge;
						computePotential	=	true;
					}
				}
				}
			}
			calculate_NodePotential_From_Wellseparated_Clusters_cache_svd(node,tree.rank,tree.nChebNodes);
			transfer_NodePotential_To_Child(node,tree.R);
			if(computePotential){
				tranfer_Potential_To_Potential_Tree(node, potential);
			}
			for(unsigned short k=0;k<4;++k){
				calculate_Potential_cache_svd(node->child[k], potential,tree);
			}
		}
	}
}

void kernel_Base::calculate_Potential_cache_svd(H2_2D_Tree& tree, double* potential) {
    MatrixXd potentialMatrix;
    potentialMatrix = MatrixXd::Zero(tree.N,tree.m);
    set_Tree_Potential_Zero(tree.root);
//    std::cout << "Calculating potential caching..." << std::endl;
    calculate_Potential_cache_svd(tree.root,potentialMatrix,tree);
    Map<MatrixXd>(potential, potentialMatrix.rows(), potentialMatrix.cols()) = potentialMatrix;
//    std::cout << "Calculated potential caching at storage of "<< cache.size() << "." << std::endl;
}

void kernel_Base::calculate_NodePotential_From_Wellseparated_Clusters_cache_svd(
		H2_2D_Node*& node, unsigned short rank,unsigned short nChebNodes) {
	MatrixXd tempK = MatrixXd::Zero(rank, rank);
	for(unsigned short k=0; k<4; ++k){
		if(!node->child[k]->isEmpty){
			for(unsigned short i=0; i<node->child[k]->nInteraction; ++i){
				if (node->child[k]->interaction[i] != NULL && !node->child[k]->interaction[i]->isEmpty) {
					kernel_Cheb_2D(nChebNodes,node->child[k]->scaledCnode,nChebNodes,node->child[k]->interaction[i]->scaledCnode,tempK);

					const int row = int(tempK.rows() * scal);
					JacobiSVD<MatrixXd> svd(tempK, ComputeThinU | ComputeThinV);
					MatrixXd u = svd.matrixU().leftCols(row);
					MatrixXd v = svd.matrixV().leftCols(row);
					VectorXd s = svd.singularValues().head(row);
					cache.push_back(u);
					cache.push_back(v);
					sigma.push_back(s);

					node->child[k]->nodePotential += tempK*
							node->child[k]->interaction[i]->nodeCharge;

				}
			}
		}
	}
}


void kernel_Base::calculate_Potential_fast_svd(H2_2D_Node*& node, MatrixXd& potential,H2_2D_Tree& tree, int*& index) {
    if(!node->isEmpty){
		if(node->isLeaf){
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
                    //	Potential from neighbors
					tree.get_Charge(node->neighbor[k]);
					node->potential+=
							cache[(*index) * 2] *
							sigma[(*index)].asDiagonal() *
							cache[(*index) * 2 + 1].transpose()*node->neighbor[k]->charge;
					(*index) = (*index) + 1;
				}
				}
			}
            //			Potential from Chebyshev nodes
			node->potential+=node->R*node->nodePotential;

			node->potential+=(cache[(*index) * 2] *
					sigma[(*index)].asDiagonal() *
					cache[(*index) * 2 + 1].transpose())*node->charge;

			(*index) = (*index) + 1;


			tranfer_Potential_To_Potential_Tree(node, potential);
		}
		else{
			bool computePotential	=	false;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
					if(node->neighbor[k]->isLeaf){
						tree.get_Charge(node->neighbor[k]);
						node->potential+=(cache[(*index) * 2] *
								sigma[(*index)].asDiagonal() *
								cache[(*index) * 2 + 1].transpose() )* node->neighbor[k]->charge;

						(*index) = (*index) + 1;
						computePotential	=	true;
					}
				}
				}
			}
			// inc index by 1
			calculate_NodePotential_From_Wellseparated_Clusters_fast_svd(node,tree.rank,tree.nChebNodes, index);
			transfer_NodePotential_To_Child(node,tree.R);
			if(computePotential){
				tranfer_Potential_To_Potential_Tree(node, potential);
			}
			for(unsigned short k=0;k<4;++k){
				calculate_Potential_fast_svd(node->child[k], potential,tree, index);
			}
		}
	}
}
void kernel_Base::calculate_Potential_fast_svd(H2_2D_Tree& tree, double* potential) {
	int index = 0;
	int* index_ptr = &index;
    MatrixXd potentialMatrix;
    potentialMatrix = MatrixXd::Zero(tree.N,tree.m);
    set_Tree_Potential_Zero(tree.root);
    calculate_Potential_fast_svd(tree.root, potentialMatrix,tree, index_ptr);
    Map<MatrixXd>(potential, potentialMatrix.rows(), potentialMatrix.cols()) = potentialMatrix;

}
void kernel_Base::calculate_NodePotential_From_Wellseparated_Clusters_fast_svd(
		H2_2D_Node*& node, unsigned short rank,unsigned short nChebNodes, int*& index) {

	for(unsigned short k=0; k<4; ++k){
		if(!node->child[k]->isEmpty){
			for(unsigned short i=0; i<node->child[k]->nInteraction; ++i){
				if (node->child[k]->interaction[i] != NULL && !node->child[k]->interaction[i]->isEmpty) {
					node->child[k]->nodePotential	+= ( cache[(*index) * 2] *
						sigma[(*index)].asDiagonal() *
						cache[(*index) * 2 + 1].transpose() )*node->child[k]->interaction[i]->nodeCharge;

					(*index) = (*index) + 1;

				}
			}
		}
	}
}




void kernel_Base::calculate_Potential_cache(H2_2D_Node*& node, MatrixXd& potential,H2_2D_Tree& tree) {
    if(!node->isEmpty){
		if(node->isLeaf){
			MatrixXd tempK;

			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
					kernel_2D(node->N , node->location, node->neighbor[k]->N, node->neighbor[k]->location, tempK);
					// cached, copied.
					cache.push_back(tempK);

					tree.get_Charge(node->neighbor[k]);
					node->potential+=tempK*node->neighbor[k]->charge;
					}
				}
			}
            //			Potential from Chebyshev nodes
			node->potential+=node->R*node->nodePotential;
            //			Self potential
			kernel_2D(node->N , node->location, node->N , node->location, tempK);
			//cached, copied
			cache.push_back(tempK);

			node->potential+=tempK*node->charge;

			tranfer_Potential_To_Potential_Tree(node, potential);
		}
		else{
			bool computePotential	=	false;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
					if(node->neighbor[k]->isLeaf){
						MatrixXd tempK;
						kernel_2D(node->N, node->location, node->neighbor[k]->N, node->neighbor[k]->location, tempK);
						//cached, copied
						cache.push_back(tempK);
						tree.get_Charge(node->neighbor[k]);
						node->potential+=tempK*node->neighbor[k]->charge;
						computePotential	=	true;
					}
				}
				}
			}
			calculate_NodePotential_From_Wellseparated_Clusters_cache(node,tree.rank,tree.nChebNodes);
			transfer_NodePotential_To_Child(node,tree.R);
			if(computePotential){
				tranfer_Potential_To_Potential_Tree(node, potential);
			}
			for(unsigned short k=0;k<4;++k){
				calculate_Potential_cache(node->child[k], potential,tree);
			}
		}
	}
}

void kernel_Base::calculate_Potential_cache(H2_2D_Tree& tree, double* potential){
    MatrixXd potentialMatrix;
    potentialMatrix = MatrixXd::Zero(tree.N,tree.m);
    set_Tree_Potential_Zero(tree.root);
//    std::cout << "Calculating potential caching..." << std::endl;
    calculate_Potential_cache(tree.root,potentialMatrix,tree);
    Map<MatrixXd>(potential, potentialMatrix.rows(), potentialMatrix.cols()) = potentialMatrix;
//    std::cout << "Calculated potential caching at storage of "<< cache.size() << "." << std::endl;
}

void kernel_Base::calculate_Potential_fast(H2_2D_Node*& node, MatrixXd& potential,H2_2D_Tree& tree, int*& index) {
    if(!node->isEmpty){
		if(node->isLeaf){
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
                    //	Potential from neighbors
					tree.get_Charge(node->neighbor[k]);
					node->potential+=cache[*index]*node->neighbor[k]->charge;
					(*index) = (*index) + 1;
				}
				}
			}
            //			Potential from Chebyshev nodes
			node->potential+=node->R*node->nodePotential;

			node->potential+=cache[*index]*node->charge;
			(*index) = (*index) + 1;


			tranfer_Potential_To_Potential_Tree(node, potential);
		}
		else{
			bool computePotential	=	false;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
					if(node->neighbor[k]->isLeaf){
						tree.get_Charge(node->neighbor[k]);
						node->potential+=cache[*index]*node->neighbor[k]->charge;
						(*index) = (*index) + 1;
						computePotential	=	true;
					}
				}
				}
			}
			// inc index by 1
			calculate_NodePotential_From_Wellseparated_Clusters_fast(node,tree.rank,tree.nChebNodes, index);
			transfer_NodePotential_To_Child(node,tree.R);
			if(computePotential){
				tranfer_Potential_To_Potential_Tree(node, potential);
			}
			for(unsigned short k=0;k<4;++k){
				calculate_Potential_fast(node->child[k], potential,tree, index);
			}
		}
	}
}

void kernel_Base::calculate_Potential_fast(H2_2D_Tree& tree, double* potential){
	int index = 0;
	int* index_ptr = &index;
    MatrixXd potentialMatrix;
    potentialMatrix = MatrixXd::Zero(tree.N,tree.m);
    set_Tree_Potential_Zero(tree.root);
    calculate_Potential_fast(tree.root, potentialMatrix,tree, index_ptr);
    Map<MatrixXd>(potential, potentialMatrix.rows(), potentialMatrix.cols()) = potentialMatrix;
}

void kernel_Base::calculate_Potential(H2_2D_Node*& node, MatrixXd& potential,H2_2D_Tree& tree){
    if(!node->isEmpty){
		if(node->isLeaf){
			MatrixXd tempK;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
					kernel_2D(node->N , node->location, node->neighbor[k]->N, node->neighbor[k]->location, tempK);
                    //	Potential from neighbors
					tree.get_Charge(node->neighbor[k]);
					node->potential+=tempK*node->neighbor[k]->charge;
				}
				}
			}
            //			Potential from Chebyshev nodes
			node->potential+=node->R*node->nodePotential;
            //			Self potential
			kernel_2D(node->N , node->location, node->N , node->location, tempK);
			node->potential+=tempK*node->charge;
            
			tranfer_Potential_To_Potential_Tree(node, potential);
		}
		else{
			bool computePotential	=	false;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isEmpty){
						if(node->neighbor[k]->isLeaf){
							MatrixXd tempK;
							kernel_2D(node->N, node->location, node->neighbor[k]->N, node->neighbor[k]->location, tempK);
							tree.get_Charge(node->neighbor[k]);
							node->potential+=tempK*node->neighbor[k]->charge;
							computePotential	=	true;
						}
					}
				}
			}
			calculate_NodePotential_From_Wellseparated_Clusters(node,tree.rank,tree.nChebNodes);
			transfer_NodePotential_To_Child(node,tree.R);
			if(computePotential){
				tranfer_Potential_To_Potential_Tree(node, potential);
			}
			for(unsigned short k=0;k<4;++k){
				calculate_Potential(node->child[k], potential,tree);
			}
		}
	}
}

void kernel_Base::set_Tree_Potential_Zero(H2_2D_Node* node){
    if (node) {
        node->potential     =   MatrixXd::Zero(node->potential.rows(),node->potential.cols());
        node->nodePotential =   MatrixXd::Zero(node->nodePotential.rows(),node->nodePotential.cols());
        for (unsigned short k=0; k<4; ++k) {
            set_Tree_Potential_Zero(node->child[k]);
        }
    }
}

//	Calculates potential;
void kernel_Base::calculate_Potential(H2_2D_Tree& tree, double* potential){
    MatrixXd potentialMatrix;
    potentialMatrix = MatrixXd::Zero(tree.N,tree.m);
    set_Tree_Potential_Zero(tree.root);
    calculate_Potential(tree.root,potentialMatrix,tree);
    Map<MatrixXd>(potential, potentialMatrix.rows(), potentialMatrix.cols()) = potentialMatrix;
}

//	Obtains Chebyshev node potential from well separated clusters;
void kernel_Base::calculate_NodePotential_From_Wellseparated_Clusters(H2_2D_Node*& node, unsigned short rank,unsigned short nChebNodes){
	MatrixXd K = MatrixXd::Zero(rank, rank);
	for(unsigned short k=0; k<4; ++k){
		if(!node->child[k]->isEmpty){
			for(unsigned short i=0; i<node->child[k]->nInteraction; ++i){
				if (node->child[k]->interaction[i] != NULL && !node->child[k]->interaction[i]->isEmpty) {
					kernel_Cheb_2D(nChebNodes,node->child[k]->scaledCnode,nChebNodes,node->child[k]->interaction[i]->scaledCnode,K);

					node->child[k]->nodePotential	=	node->child[k]->nodePotential+K*node->child[k]->interaction[i]->nodeCharge;

				}
			}
		}
	}
}

void kernel_Base::calculate_NodePotential_From_Wellseparated_Clusters_cache(H2_2D_Node*& node, unsigned short rank,unsigned short nChebNodes){
	MatrixXd K = MatrixXd::Zero(rank, rank);
	for(unsigned short k=0; k<4; ++k){
		if(!node->child[k]->isEmpty){
			for(unsigned short i=0; i<node->child[k]->nInteraction; ++i){
				if (node->child[k]->interaction[i] != NULL && !node->child[k]->interaction[i]->isEmpty) {
					kernel_Cheb_2D(nChebNodes,node->child[k]->scaledCnode,nChebNodes,node->child[k]->interaction[i]->scaledCnode,K);

					//cached, copied.
					cache.push_back(K);


					node->child[k]->nodePotential	+= K*node->child[k]->interaction[i]->nodeCharge;

				}
			}
		}
	}
}

void kernel_Base::calculate_NodePotential_From_Wellseparated_Clusters_fast(
		H2_2D_Node*& node, unsigned short rank,unsigned short nChebNodes, int*& index){
//	MatrixXd K = MatrixXd::Zero(rank, rank);
	for(unsigned short k=0; k<4; ++k){
		if(!node->child[k]->isEmpty){
		for(unsigned short i=0; i<node->child[k]->nInteraction; ++i){
            if (node->child[k]->interaction[i] != NULL && !node->child[k]->interaction[i]->isEmpty) {
                node->child[k]->nodePotential	+=
                		cache[*index]*node->child[k]->interaction[i]->nodeCharge;
                (*index) = (*index) + 1;

            }
		}
	}
	}
}

//	Tranfers potential from node to final potential matrix when needed;
void kernel_Base::tranfer_Potential_To_Potential_Tree(H2_2D_Node*& node, MatrixXd& potential){
	for(unsigned long k=0; k<node->N; ++k){
			potential.row((size_t)node->index(k))+=node->potential.row(k);
	}
}

//	Evaluate kernel at Chebyshev nodes;
void kernel_Base::kernel_Cheb_2D(const unsigned short& M, const vector<Point>& xVec, const unsigned short& N, const vector<Point>& yVec, MatrixXd& K){
	vector<Point> xNew;
	vector<Point> yNew;
	K		=	MatrixXd::Zero(M*M,N*N);
	for(unsigned short j=0;j<M;++j){
		for(unsigned short i=0;i<M;++i){
            Point newPoint;
            newPoint.x     =   xVec[i].x;
            newPoint.y     =   xVec[j].y;
			xNew.push_back(newPoint);
		}
	}
	for(unsigned short j=0;j<N;++j){
		for(unsigned short i=0;i<N;++i){
            Point newPoint;
            newPoint.x     =   yVec[i].x;
            newPoint.y     =   yVec[j].y;
            yNew.push_back(newPoint);
		}
	}
	kernel_2D(M*M, xNew, N*N, yNew, K);
}

//	Tranfers potential from Chebyshev node of parent to Chebyshev node of children;
void kernel_Base::transfer_NodePotential_To_Child(H2_2D_Node*& node, MatrixXd R[]){
	for(unsigned short k=0;k<4;++k){
		if(!node->child[k]->isEmpty){
			node->child[k]->nodePotential	=	node->child[k]->nodePotential+R[k]*node->nodePotential;
		}
	}
}


void kernel_Base::kernel_2D(const unsigned long M, const vector<Point>& x, const unsigned long N, const vector<Point>& y, MatrixXd& kernel) {
	kernel	=	MatrixXd::Zero(M,N);
	for(unsigned long i=0;i<M;++i){
		for(unsigned long j=0;j<N;++j){
            kernel(i,j) =   kernel_Func(x[i],y[j]);
        }
    }
}

