/*
 * Cell.cpp
 *
 *  Created on: Sep 13, 2015
 *      Author: lurker
 */

#include "Cell.h"

Cell::Cell(const std::vector<double>& pos, double sz) :
	position(pos),
	size(sz),
	parent(nullptr),
	intensity(0.),
	status(CellStatus::UNSET){center.resize(2);};

Cell::~Cell() {
	for(auto child : children) {
		delete child;
	}
}

void Cell::deleteParticles() noexcept {
	for(auto p : particles) {
		delete p;
	}
}
/*
 * split grid into 4 equal parts until empty or only one particle inside.
 *
 */
void Cell::split() noexcept {
	if (getNumParticles() == 0) {
		setStatus(CellStatus::EMPTY);
		return;
	}
	else if (getNumParticles() == 1) {
		setStatus(CellStatus::LEAF);
		return;
	}
	else {

		if (getStatus() != CellStatus::ROOT) {
			setStatus(CellStatus::BRANCH);
		}

		/*
		 * split until single photon recursively.
		 */
		for (auto i = 0; i < 4; i++) {
			std::vector<double> pos {
				this->position[0] + this->size/2.0 * (i & 1),
				this->position[1] + this->size/2.0 * ((i >> 1) & 1)
			};

			auto child = new Cell(pos, this->size/2.0);
			child->setParent(this);
			children.push_back(child);
		}

		for(auto photon : particles) {
			/*
			 * split photons
			 */
			size_t index = 0;
			for (auto j = 0; j < 2; j++) {
				int temp = (photon->position[j] - this->position[j])/(this->size/2.0);
				/*
				 * the only case that temp == 2 is locating on right most edge. This case
				 * is considered as inside. A zero measure difference will result in 4 sections.
				 *
				 */
				temp = (temp == 2) ? temp - 1 : temp;
				index |= (temp << j);
			}
			this->children[index]->addParticle(photon);
		}

		for(auto child : children) {
			/*
			 * calculate intensity and center
			 *
			 * slow top bottom swiping
			 */
			child->updateCenterIntensity();
			child->split();
		}
	}
}

void Cell::updateCenterIntensity() noexcept {
	if (this->getNumParticles() != 0) {
		this->intensity = 0.;
		for(auto photon : this->particles) {
			this->intensity += photon->intensity;
			this->center[0] += photon->intensity * photon->position[0];
			this->center[1] += photon->intensity * photon->position[1];
		}
		this->intensity /= this->getNumParticles();
		this->center[0] /= this->getNumParticles();
		this->center[1] /= this->getNumParticles();
	}
}

size_t Cell::getNumParticles() const noexcept {
	return particles.size();
}

void Cell::addParticle(Photon* p) noexcept{
	particles.push_back(p);
}

CellStatus Cell::getStatus() const noexcept {
	return status;
}

void Cell::setStatus(CellStatus status) noexcept{
	this->status = status;
}

void Cell::clearParticles() noexcept {
	particles.clear();
}

void Cell::setSize(double size) noexcept {
	this->size = size;
}

double Cell::getIntensity() noexcept {
	return this->intensity;
}

double Cell::getSize() const noexcept{
	return size;
}

Cell* Cell::getParent() const noexcept {
	return parent;
}

void Cell::setParent(Cell* parent) noexcept {
	this->parent = parent;
}

std::vector<Cell*>& Cell::getChildren() noexcept {
	return this->children;
}

std::vector<double>& Cell::getCenter() noexcept {
	return this->center;
}

/*
 * test function
 */
inline double eval(std::vector<double>& x, std::vector<double>& y) {
	return exp(-0.1 * sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1])));
}

inline double visit(Photon* photon, Cell* cell) {
//	std::cout << "this cell has " << cell->getNumParticles() << " particles, intensity "
//			<< cell->getIntensity() << " at "<<cell->getCenter() << std::endl;

	return eval(photon->position, cell->getCenter()) * cell->getSize() * cell->getSize() /
			sqrt(
					(photon->position[0] - cell->getCenter()[0]) *
					(photon->position[0] - cell->getCenter()[0]) +
					(photon->position[1] - cell->getCenter()[1]) *
					(photon->position[1] - cell->getCenter()[1])
				);
}

inline void trasverse(Photon* photon, Cell* cell, double* ptr, size_t n) {
	double distance = sqrt(
			(photon->position[0] - cell->getCenter()[0]) *
			(photon->position[0] - cell->getCenter()[0]) +
			(photon->position[1] - cell->getCenter()[1]) *
			(photon->position[1] - cell->getCenter()[1])
		);
	if (cell != nullptr && (cell->getStatus()== CellStatus::LEAF ||
			cell->getSize()/distance < 0.2) ) {
		auto result = visit(photon, cell);
		for (auto _photon : cell->particles) {
			if (photon->id != _photon->id) {
				ptr[photon->id * n + _photon->id] = result/cell->getNumParticles();
			}
		}
	}
	else{
		for(auto child : cell->getChildren()) {
			if (child->getStatus() != CellStatus::EMPTY) {
				trasverse(photon, child, ptr, n);
			}
		}
	}
}



template class mexplus::Session<Cell>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 1);

	auto posptr = mxGetPr(prhs[0]);
	std::vector<double> pos {posptr[0], posptr[1]};
	auto szptr = mxGetPr(prhs[1]);

	output.set(0, Session<Cell>::create(new Cell(
			pos,
			*szptr)));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Cell* root = Session<Cell>::get(input.get(0));

	root->deleteParticles();

	Session<Cell>::destroy(input.get(0));
}

MEX_DEFINE(import) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])  {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 0);

	auto particles = mxGetPr(prhs[1]);
	auto numofparticle = mxGetN(prhs[1]);

	auto root = Session<Cell>::get(input.get(0));

	for (int i = 0; i < numofparticle; i++) {
		/*
		 * photon is (x, y, value) form
		 */
		auto ph = new Photon(particles[3 * i], particles[3 *  i + 1], particles[3 * i + 2]);
		ph->id = i;
		root->addParticle(ph);
	}

	root->updateCenterIntensity();
}

MEX_DEFINE(split) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);

	auto root = Session<Cell>::get(input.get(0));
	root->setStatus(CellStatus::ROOT);
	root->split();
}

MEX_DEFINE(buildmatrix) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);

	auto root = Session<Cell>::get(input.get(0));

	auto num = root->getNumParticles();
	if (num == 0) {
		std::cout << "empty cell, stopped" << std::endl;
		return;
	}
	plhs[0] = mxCreateNumericMatrix(num, num,mxDOUBLE_CLASS, mxREAL);
	auto Radptr  = mxGetPr(plhs[0]);

	for (auto photon : root->particles) {
		trasverse(photon, root, Radptr, num);
	}
}

}

MEX_DISPATCH
