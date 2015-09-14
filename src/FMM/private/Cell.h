/*
 * Cell.h
 *
 *  Created on: Sep 13, 2015
 *      Author: lurker
 */

#ifndef CELL_H_
#define CELL_H_

#include <cstdlib>
#include <vector>
#include <queue>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <mexplus.h>
#include <pprint.h>

#include "utils.h"

using namespace std;
using namespace mexplus;

class Photon {
public:
	Photon(double x, double y) {
		position[0] = x;
		position[1] = y;
	};
	virtual ~Photon() {

	};
	double position[2];
};

enum class CellStatus {
	ROOT,
	BRANCH,
	LEAF,
	EMPTY,
	UNSET
};

class Cell {
public:
	Cell(const std::vector<double>& pos, double sz);

	virtual ~Cell();

	void deleteParticles();

	void split() noexcept;


	size_t getNumParticles() const noexcept;

	void addParticle(Photon* p) noexcept;

	CellStatus getStatus() const noexcept;

	void setStatus(CellStatus status) noexcept;

	void clearParticles() noexcept;

	void setSize(double size) noexcept;

	double getSize() const noexcept;

	Cell* getParent() const noexcept;

	void setParent(Cell* parent) noexcept;


private:
	std::vector<double> position;
	double size;
	std::vector<Photon*> particles;
	std::vector<Cell*> children;
	CellStatus status;

	double intensity;
	Cell* parent;
};

#endif /* CELL_H_ */
