//
// Created by lurker on 9/19/15.
//

#ifndef QUADTREE_TREECODE_H
#define QUADTREE_TREECODE_H
#define EPS 1e-2

#include "quadtree.h"


/*
 * static arrays for integration
 */
static std::vector<scalar_t> X {
        -0.90617984593866396370032134655048,
        -0.53846931010568310771446931539685,
        0,
        0.53846931010568310771446931539685,
        0.90617984593866396370032134655048

};
static std::vector<scalar_t> W {
        0.23692688505618908489935847683228,
        0.47862867049936647090291330641776,
        0.56888888888888888888888888888889,
        0.47862867049936647090291330641776,
        0.23692688505618908489935847683228
};

typedef int level_t;
/*
 * treecode structure
 * @members
 *
 * @x : x coordinate of left bottom corner
 * @y : y coordinate of left bottom corner
 * @length : length of side
 * @size : number of points per slice, it is 2^n
 * @root : root pointer of quadtree
 *
 */
typedef struct treecode {
    scalar_t  x;
    scalar_t  y;
    scalar_t  length;
    level_t   size;
    quadtree*  root;

    /*
     * treecode constructor
     * @params
     *
     */
    explicit treecode(scalar_t _x_start, scalar_t _y_start,
              scalar_t _length, level_t _level)  noexcept :
            x(_x_start), y(_y_start), length(_length), size(1 << _level){
        ord_t order = 0;
        root = new quadtree(length, x, y);
        root->setStatus(Status::ROOT);
        for (auto i = 1; i <= size; i++) {
            for (auto j = 1; j <= size; j++) {
                // assign points with
                auto ptr = make_shared<point>(
                        x + (2.0 * i - 1.) * length/(2.0 * size),
                        y + (2.0 * j - 1.) * length/(2.0 * size)
                );
                ptr->id = order++;
                root->addPoint(ptr);
            }
        }
    }
    /*
     * treecode destructor
     */
    ~treecode() noexcept {delete root;}

    /*
     * return attribute at grid containing (x0, y0)
     */
    attribute_t getAttribute(scalar_t x0, scalar_t y0) noexcept ;

    /*
     * set attribute at the grid containing (x0, y0)
     */
    void setAttribute(attribute_t attr, scalar_t x0, scalar_t y0) noexcept ;

} treecode;


/*
 * external functions
 */

/*
 *  return distance between two locations
 *  @params
 *
 *  @x0 : x coordinate
 *  @y0 : y coordinate
 *  @x1 : x coordinate
 *  @y1 : y coordinate
 */
inline scalar_t distance(scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {
    return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
}

inline scalar_t integral_helper(treecode *tree, scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {
    return distance(x0, y0, x1, y1) * tree->getAttribute((x0 + x1)/2 , (y0 + y1)/2);
}

/*
 * return line integral between two locations
 *
 * @params
 *
 * @tree : treecode structure
 * @x0 : x coordinate
 * @y0 : y coordinate
 * @x1 : x coordinate
 * @y1 : y coordinate
 */
inline scalar_t integral(treecode *tree, scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {

    if (distance(x0, y0, x1, y1) < EPS) return integral_helper(tree, x0, y0, x1, y1);
    else {
        auto xm = (x0 + x1)/2;
        auto ym = (y0 + y1)/2;
        return integral(tree, x0, y0, xm, ym) +
               integral(tree, xm, ym, x1, y1);
    }
}

/*
 * auxiliary function for evaluation
 */
inline scalar_t eval_helper(treecode *tree, scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {
    return exp(-integral(tree, x0, y0, x1, y1)) / distance(x0, y0, x1, y1);
}
/*
 * return evaluation of kernel over a grid
 * @params
 *
 * @tree : treecode pointer
 * @x0 : x coordinate
 * @y0 : y coordinate
 * @branch_ptr : target quadtree
 */
inline scalar_t eval(treecode *tree, scalar_t x0, scalar_t y0, quadtree* branch_ptr) noexcept {
    auto sum = scalar_t(0.0);
    auto center_x = branch_ptr->x +  branch_ptr->length/2;
    auto center_y = branch_ptr->y +  branch_ptr->length/2;
    for (auto i = 0; i < X.size(); i++) {
        auto x_ = center_x +  X[i] * branch_ptr->length/2;
        for (auto j = 0; j < X.size(); j++) {
            sum += eval_helper(tree, x0, y0,
                               x_,
                               center_y + X[j] * branch_ptr->length/2) * W[i] * W[j] / 4;
        }
    }
    return sum * branch_ptr->length * branch_ptr->length;
}
/*
 * traversal quadtree
 *
 * @params
 *
 * @tree : treecode structure pointer
 * @theta : parameter to control accuracy
 * @point_ptr : point pointer
 * @branch_ptr: target quadtree
 * @n : matrix size
 * @matrix_ptr : matrix pointer
 */
inline void traversal(treecode *tree, scalar_t& theta,
                      shared_ptr<point> point_ptr, quadtree* branch_ptr,
                      ord_t& n, scalar_t* matrix_ptr) noexcept {

    auto d = distance(point_ptr->x, point_ptr->y, branch_ptr->x, branch_ptr->y);

    if (branch_ptr->status == Status::LEAF || branch_ptr->length/ d < theta) {
        // if reaches a leaf or in far field
        auto ret = eval(tree, point_ptr->x, point_ptr->y, branch_ptr);
        for (auto point : branch_ptr->points) {
            if (point->id != point_ptr->id) {
                matrix_ptr[n * point->id + point_ptr->id] = ret / branch_ptr->points.size();
            }
            else {
                matrix_ptr[n * point->id + point_ptr->id] =
                        2 * scalar_t(M_PI) *
                        (1 - exp(-branch_ptr->attribute * branch_ptr->length))/branch_ptr->attribute;
            }
        }
    }
    else {
        // else recursively traverse.
        for (auto child : branch_ptr->children) {
            if (child->status != Status::EMPTY) {
                traversal(tree, theta, point_ptr, child.get(), n, matrix_ptr);
            }
        }
    }
}


#endif //QUADTREE_TREECODE_H
