//
// Created by lurker on 9/19/15.
//

#ifndef QUADTREE_TREECODE_H
#define QUADTREE_TREECODE_H
//#define EPS 1e-2

#include <unordered_map>
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
//
//static std::vector<scalar_t> X {0};
//static std::vector<scalar_t> W {2.0};

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
    vector<unordered_map<quadtree*, scalar_t>> interactions;

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
        /*
         * initialize interaction table
         */
        interactions.resize(root->points.size());
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

inline size_t getRow(treecode *tree, scalar_t y) noexcept {
	return size_t(floor(tree->size * y));
}

inline size_t getCol(treecode *tree, scalar_t x) noexcept {
	return size_t(floor(tree->size * x));
}

inline scalar_t integral(treecode *tree, scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {
	auto side = tree->size;
	auto col0 = getRow(tree, x0);
	auto col1 = getRow(tree, x1);
	auto row0 = getCol(tree, y0);
	auto row1 = getCol(tree, y1);

	// 9 cases
	if ((row0 == row1) && (col0 == col1)) {
		return integral_helper(tree, x0, y0, x1, y1);
	}
	else if ((row0 == row1 + 1) && (col0 == col1)) {
		auto row0 = getCol(tree, y0);
		scalar_t ybar = scalar_t(row0)/side;
		scalar_t xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
		return integral_helper(tree, x0, y0, xbar, ybar) + integral_helper(tree, xbar, ybar, x1, y1);
	}
	else if ((row0 == row1 - 1) && (col0 == col1)) {
		auto row1 = getCol(tree, y1);
		scalar_t ybar = scalar_t(row1)/side;
		scalar_t xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
		return integral_helper(tree, x0, y0, xbar, ybar) + integral_helper(tree, xbar, ybar, x1, y1);
	}
	else if ((col0 == col1 + 1) && (row0 == row1)) {
		auto col0 = getCol(tree, x0);
		scalar_t xbar = scalar_t(col0)/side;
		scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
		return integral_helper(tree, x0, y0, xbar, ybar) + integral_helper(tree, xbar, ybar, x1, y1);
	}
	else if ((col0 == col1 - 1) && (row0 == row1)) {
		auto col1 = getCol(tree, x1);
		scalar_t xbar = scalar_t(col1)/side;
		scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
		return integral_helper(tree, x0, y0, xbar, ybar) + integral_helper(tree, xbar, ybar, x1, y1);
	}
	else if ((col0 == col1 + 1) && (row0 == row1 + 1)) {
		scalar_t xbar = scalar_t(col0)/side;
		scalar_t ybar2 = scalar_t(row0)/side;
		scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
		scalar_t xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

		if (xbar < xbar2) {
			return integral_helper(tree, x1, y1, xbar, ybar) + \
					integral_helper(tree, xbar, ybar, xbar2, ybar2) + integral_helper(tree, xbar2, ybar2, x0, y0);
		}
		else {
			return integral_helper(tree, x1, y1, xbar2, ybar2) + \
					integral_helper(tree, xbar, ybar, xbar2, ybar2) + integral_helper(tree, xbar, ybar, x0, y0);

		}
	}
	else if ((col0 == col1 + 1) && (row0 == row1 - 1)) {
		scalar_t xbar = scalar_t(col0)/side;
		scalar_t ybar2 = scalar_t(row1)/side;
		scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
		scalar_t xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

		if (xbar < xbar2) {
			return integral_helper(tree, x1, y1, xbar, ybar) + \
					integral_helper(tree, xbar, ybar, xbar2, ybar2) + integral_helper(tree, xbar2, ybar2, x0, y0);
		}
		else {
			return integral_helper(tree, x1, y1, xbar2, ybar2) + \
					integral_helper(tree, xbar, ybar, xbar2, ybar2) + integral_helper(tree, xbar, ybar, x0, y0);

		}

	}
	else if ((col0 == col1 - 1) && (row0 == row1 + 1)) {
		scalar_t xbar = scalar_t(col1)/side;
		scalar_t ybar2 = scalar_t(row0)/side;
		scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
		scalar_t xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

		if (xbar > xbar2) {
			return integral_helper(tree, x1, y1, xbar, ybar) + \
					integral_helper(tree, xbar, ybar, xbar2, ybar2) + integral_helper(tree, xbar2, ybar2, x0, y0);
		}
		else {
			return integral_helper(tree, x1, y1, xbar2, ybar2) + \
					integral_helper(tree, xbar, ybar, xbar2, ybar2) + integral_helper(tree, xbar, ybar, x0, y0);
		}
	}
	else if ((col0 == col1 - 1) && (row0 == row1 - 1)) {
		scalar_t xbar = scalar_t(col1)/side;
		scalar_t ybar2 = scalar_t(row1)/side;
		scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
		scalar_t xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

		if (xbar > xbar2) {
			return integral_helper(tree, x1, y1, xbar, ybar) + \
					integral_helper(tree, xbar, ybar, xbar2, ybar2) + integral_helper(tree,  xbar2, ybar2, x0, y0);
		}
		else {
			return integral_helper(tree, x1, y1, xbar2, ybar2) + \
					integral_helper(tree, xbar, ybar, xbar2, ybar2) + integral_helper(tree, xbar, ybar, x0, y0);

		}
	}
	else {
		scalar_t xm = (x0 + x1)/2;
		scalar_t ym = (y0 + y1)/2;
		return integral(tree, x0, y0, xm, ym) + integral(tree, xm, ym, x1, y1);
	}
}


inline scalar_t integral_approximate(treecode *tree, scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {

    if (distance(x0, y0, x1, y1) < 1.0/tree->size) return integral_helper(tree, x0, y0, x1, y1);
    else {
        auto xm = (x0 + x1)/2;
        auto ym = (y0 + y1)/2;
        return integral_approximate(tree, x0, y0, xm, ym) +
               integral_approximate(tree, xm, ym, x1, y1);
    }
}

/*
 * auxiliary function for evaluation
 */
inline scalar_t eval_helper(treecode *tree, scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {
    return exp(-integral(tree, x0, y0, x1, y1)) / distance(x0, y0, x1, y1);
}

inline scalar_t integral_corner(scalar_t sigma_t, scalar_t length) {
    scalar_t ret = 0.;
    scalar_t mid = (sqrt(2) + 1.0) * length/4.0;
    scalar_t rad = (sqrt(2) - 1.0) * length/4.0;
    for (int i = 0; i < X.size(); ++i) {
        scalar_t val = (mid + rad * X[i]);
        ret += exp(-sigma_t * val) * (M_PI/2.0 - 2 * acos(length/2.0/val)) * W[i];
    }
    ret *= 4.0 * rad;
    return ret;
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
//    auto sum = scalar_t(0.0);
//    auto center_x = branch_ptr->x +  branch_ptr->length/2;
//    auto center_y = branch_ptr->y +  branch_ptr->length/2;
//    for (auto i = 0; i < X.size(); i++) {
//        auto x_ = center_x +  X[i] * branch_ptr->length/2;
//        for (auto j = 0; j < X.size(); j++) {
//            sum += eval_helper(tree, x0, y0,
//                               x_,
//                               center_y + X[j] * branch_ptr->length/2) * W[i] * W[j] / 4;
//        }
//    }
//    return sum * branch_ptr->length * branch_ptr->length;
	    auto center_x = branch_ptr->x +  branch_ptr->length/2;
	    auto center_y = branch_ptr->y +  branch_ptr->length/2;

	    return eval_helper(tree, x0, y0,
						   center_x,
						   center_y) * branch_ptr->length * branch_ptr->length;
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
        for (auto& point : branch_ptr->points) {
            if (point->id != point_ptr->id) {
                matrix_ptr[n * point->id + point_ptr->id] = ret / branch_ptr->points.size();
            }
            else {
                matrix_ptr[n * point->id + point_ptr->id] =
                        /*
                         * todo: corner parts missing. needs erf function to be implemented.
                         */
                        2 * scalar_t(M_PI) *
                        (1 - exp(-branch_ptr->attribute * branch_ptr->length)/2.0)/branch_ptr->attribute +
                        integral_corner(branch_ptr->attribute, branch_ptr->length);
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

/*
 * import rhs into tree, take point's value slot.
 *
 * @treecode treecode structure pointer
 * @n number of points at finest level
 * @rhs array of values to assign to treecode
 *
 *
 * top-bottom to leaf. assign values to quadtree
 *
 */
inline void traversal_down(quadtree *root,scalar_t* rhs) {
    for(auto& child : root->children) {
        if (child->status == Status::LEAF) {
            child->value = rhs[child->points[0]->id];
        }
        else {
            traversal_down(child.get(), rhs);
        }
    }
}


inline bool traversal_check(quadtree *root){
	auto ret = 0.;
	auto ret_ = true;
	for(auto& child : root->children) {
		ret += child->value;
	}
	if (ret != 4 * root->value) return false;
	else return true;
}


/*
 * bottom-top post-order traversal of tree, averaging
 */
inline scalar_t traversal_up(quadtree *root) {
    scalar_t ret = 0.;
    if (root->status != Status::LEAF) {
        for (auto& child : root->children) {
            ret += traversal_up(child.get());
        }
        ret /= 4.0;
    }
    else {
    	ret = root->value;
    }
    root->value = ret;
    return ret;
}


/*
 * traversal quadtree, without matrix, calculate iteration.
 *
 * @params
 *
 * @tree : treecode structure pointer
 * @theta : parameter to control accuracy
 * @point_ptr : point pointer
 * @branch_ptr: target quadtree
 * @n : matrix size
 * @lhs : vector result
 * @rhs : input source term
 */
inline void traversal (treecode *tree, scalar_t& theta,
                       shared_ptr<point> point_ptr, quadtree* branch_ptr,
                       ord_t& n, scalar_t* lhs, scalar_t* rhs) noexcept {
    /*
     * distance from current node to target branch
     */
    auto d = distance(point_ptr->x, point_ptr->y, branch_ptr->x, branch_ptr->y);
    /*
     * pre-order traversal quadtree.
     */
    if (branch_ptr->status == Status::LEAF || branch_ptr->length / d < theta) {
        // if target branch is a leaf or at far field.
        // todo: finish iteration
        auto ret = eval(tree, point_ptr->x, point_ptr->y, branch_ptr);
        tree->interactions[point_ptr->id][branch_ptr] = ret;

        if (branch_ptr->points[0]->id != point_ptr->id) {
        	lhs[point_ptr->id] += ret * branch_ptr->value;

        }
        else {
            lhs[point_ptr->id] +=
                    (2 * scalar_t(M_PI) *
                     (1 - exp(-branch_ptr->attribute * branch_ptr->length/2.0))
                     /branch_ptr->attribute +
                     integral_corner(branch_ptr->attribute, branch_ptr->length))
                     * rhs[point_ptr->id];
        }
    }
    else {
        for (auto child : branch_ptr->children) {
            if (child->status != Status::EMPTY) {
                traversal(tree, theta, point_ptr, child.get(), n, lhs, rhs);
            }
        }
    }
}

inline void fast_traversal (treecode *tree, scalar_t& theta,
                       shared_ptr<point> point_ptr, quadtree* branch_ptr,
                       ord_t& n, scalar_t* lhs, scalar_t* rhs) noexcept {
    /*
     * distance from current node to target branch
     */
    auto d = distance(point_ptr->x, point_ptr->y, branch_ptr->x, branch_ptr->y);
    /*
     * pre-order traversal quadtree.
     */
    if (branch_ptr->status == Status::LEAF || branch_ptr->length / d < theta) {
        // if target branch is a leaf or at far field.
        // todo: finish iteration
        auto ret = tree->interactions[point_ptr->id][branch_ptr];

        if (branch_ptr->points[0]->id != point_ptr->id) {

        	lhs[point_ptr->id] += ret * branch_ptr->value;
        }
        else {
            lhs[point_ptr->id] +=
                    (2 * scalar_t(M_PI) *
                     (1 - exp(-branch_ptr->attribute * branch_ptr->length/2.0))
                     /branch_ptr->attribute +
                     integral_corner(branch_ptr->attribute, branch_ptr->length))
                     * rhs[point_ptr->id];
        }
    }
    else {
        for (auto child : branch_ptr->children) {
            if (child->status != Status::EMPTY) {
                fast_traversal(tree, theta, point_ptr, child.get(), n, lhs, rhs);
            }
        }
    }
}


#endif //QUADTREE_TREECODE_H
