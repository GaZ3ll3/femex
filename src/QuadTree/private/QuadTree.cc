/*
 * QuadTree.cpp
 *
 *  Created on: Sep 18, 2015
 *      Author: lurker
 */

#include "QuadTree.h"

QuadTree::QuadTree(const vector<scalar_t>& _position_, scalar_t _size_):
	position(_position_),
	size(_size_),
	parent(nullptr),
	intensity(0.),
	status(Status::UNSET){}

QuadTree::~QuadTree() {
	for (auto child : children) {
		delete child;
	}
}

/*
 * add photon to quadtree
 */
void QuadTree::addPhoton(Photon* _photon_) noexcept {
	photons.push_back(_photon_);
}

/*
 * clear photons in a quadtree
 */
void QuadTree::clearPhotons() noexcept {
	for (auto photon : photons) {
		delete photon;
	}
	photons.clear();
}

/*
 * recursively populate quadtree structure.
 */
void QuadTree::populate() noexcept {
	if      (getPhotonSize() == 0) { setStatus(Status::EMPTY); return; }
	else if (getPhotonSize() == 1) { setStatus(Status::LEAF);  return; }
	else {
		/*
		 * reach a branch or at root
		 */
		if (getStatus() != Status::ROOT) {
			setStatus(Status::BRANCH);
		}

		for (auto i = 0; i < 4; i++) {
			vector<double> child_position {
				position[0] + size/2.0 * (i & 1),
				position[1] + size/2.0 * ((i >> 1) & 1)
			};

			auto child = new QuadTree(child_position, size/2.0);
			child->setParent(this);
			children.push_back(child);
		}

		for(auto photon : photons) {
			/*
			 * distributes photons
			 */
			size_t index = 0;
			for (auto j = 0; j < 2; j++) {
				int temp = (photon->position[j] - position[j])/(size/2.0);
				/*
				 * the only case that temp == 2 is photon locates at boundary.
				 * 1. avoid putting photons(nodes) on boundary
				 * 2. if necessary, boundary photons will be considered inside.
				 */
				temp = (temp == 2) ? temp - 1 : temp;
				index |= (temp << j);
			}
			children[index]->addPhoton(photon);
		}

		for(auto child : children) {
			/*
			 * calculate intensity and center
			 *
			 * slow top bottom swiping
			 */
			child->updateAttribute();
			child->populate();
		}
	}
}

/*
 * update the averaged intensity over a quadtree.
 */
void QuadTree::updateAttribute() noexcept {
	if (getPhotonSize() != 0) {
		this->intensity = 0.;
		for(auto photon : photons) {
			this->intensity += photon->intensity;
		}
		this->intensity /= getPhotonSize();
	}
}

QuadTree* QuadTree::getParent() noexcept {
	return parent;
}

size_t QuadTree::getPhotonSize() noexcept {
	return photons.size();
}

Status QuadTree::getStatus() noexcept {
	return status;
}

void QuadTree::setParent(QuadTree* _parent_) noexcept {
	parent = _parent_;
}

void QuadTree::setStatus(Status _status_) noexcept {
	status = _status_;
}






