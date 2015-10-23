/*
 * quadtree.h
 *
 *  Created on: Sep 18, 2015
 *      Author: lurker
 */
#ifndef QUADTREE_H
#define QUADTREE_H


#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string.h>
#include <memory>
#include <cassert>

//#include "tbb/tbb.h"

using namespace std;
//using namespace tbb;

typedef double scalar_t;
typedef double attribute_t;
typedef int32_t ord_t;


/*
 * Status of quadtree
 */
enum class Status {
    ROOT,
    BRANCH,
    LEAF,
    EMPTY,
    UNSET
};
/*
 *  point structure.
 *  @members:
 *
 *  @id : the ordering of each point
 *  @x  : x coordinate
 *  @y  : y coordinate
 *  @attribute: attribute at the point
 */
typedef struct point{
    ord_t         id;
    scalar_t      x;
    scalar_t      y;
    attribute_t   attribute;
    attribute_t   value;

    point(scalar_t _x, scalar_t _y, attribute_t _attribute) :
            x(_x), y(_y), attribute(_attribute), id(0), value(0){}
    point(scalar_t _x, scalar_t _y) :
            x(_x), y(_y), attribute(0.0f), id(0), value(0) {}
    // non-copyable
    point(const point&) = delete;
    point& operator=(const point&) = delete;
    ~point() {}
} point;

/*
 * quadtree structure
 * @members:
 *
 * @x : x coordinate of bottom left corner
 * @y : y coordinate of bottom right corner
 * @children : pointers to 4 sub-branches
 * @parent : parent quadtree pointer (non-owning)
 * @points : points located at this quadtree branch
 * @length : side length of this branch
 * @attribute : averaged attribute value from all its points.
 * @status : status of current quadtree
 */
typedef struct quadtree {
    scalar_t                       x;
    scalar_t                       y;
    vector<shared_ptr<quadtree>>   children;
    quadtree*                      parent;
    vector<shared_ptr<point>>      points;
    scalar_t                       length;
    attribute_t                    attribute;
    attribute_t                    value;
    Status                         status;


    /*
     * constructor
     * @params
     *
     * @length: length of side of square
     * @x_start: x coordinate of bottom left corner
     * @y_start: y coordinate of bottom left corner
     */
    explicit quadtree(scalar_t _length, scalar_t _x_start, scalar_t _y_start) noexcept :
            length(_length),
            x(_x_start),
            y(_y_start),
            attribute(0.0f),
            value(0.0f),
            parent(nullptr),
            status(Status::UNSET){}
    // non-copyable
    quadtree(const quadtree&) = delete;
    quadtree& operator=(const quadtree&) = delete;

    /*
     * delete invoke destructor of vector, which will reduce the reference
     * of shared_ptr, and eventually handles all the memory.
     */
    ~quadtree() noexcept {}

    /*
     * addPoint: add point to points
     * @params
     *
     * @point_ptr: pointer to a point
     */
    void addPoint(shared_ptr<point> ptr) noexcept ;

    /*
     * populate: populate a quadtree
     * @no params
     */
    void populate() noexcept ;

    /*
     * set status of quadtree, only used for root.
     * @params
     *
     * @status: status to be assigned.
     */
    void setStatus(Status _status) noexcept ;

    /*
     * update attribute on each quadtree. It takes the
     * averaged attribute of all points inside current
     * quadtree.
     * @no params
     */
    void updateAttribute() noexcept ;

} quadtree;



#endif //QUADTREE_H
