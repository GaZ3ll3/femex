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


class Integral {

};

class Photon {
public:
	Photon(double _x, double _y, double _intensity) :
		position {_x, _y},
		intensity(_intensity),
		id(0){};
	virtual ~Photon() {

	};
	std::vector<double> position;
	double intensity;
	int id;
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

	void deleteParticles() noexcept;

	void split() noexcept;

	size_t getNumParticles() const noexcept;

	void addParticle(Photon* p) noexcept;

	CellStatus getStatus() const noexcept;

	void setStatus(CellStatus status) noexcept;

	void clearParticles() noexcept;

	void updateCenterIntensity() noexcept;

	double getIntensity() noexcept;

	std::vector<double>& getCenter() noexcept;

	void setSize(double size) noexcept;

	double getSize() const noexcept;

	Cell* getParent() const noexcept;
	std::vector<Cell*>& getChildren() noexcept;

	void setParent(Cell* parent) noexcept;

	friend void buildmatrix(Cell* cell) noexcept;

	std::vector<Photon*> particles;
	std::vector<double> position;


private:

	double size;
	std::vector<Cell*> children;
	CellStatus status;

	double intensity;
	std::vector<double> center;
	Cell* parent;
};

#endif /* CELL_H_ */
