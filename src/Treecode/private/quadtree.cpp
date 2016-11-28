//
// Created by lurker on 9/18/15.
//


#include "quadtree.h"

void quadtree::addPoint(shared_ptr<point> ptr) noexcept {
    points.push_back(ptr);
}

void quadtree::populate() noexcept  {
    if (points.size() == 0) {
        // no points located in current quadtree branch
        status = Status::EMPTY; return;
    }
    else if (points.size() == 1) {
        // just a single point located in current branch
        status = Status::LEAF; return;

    }
    else {
        if (status != Status::ROOT) {
            status = Status::BRANCH;
        }
        else {
            updateAttribute();
        }
        // split quadtree when more than one points located
        for (auto i = 0 ; i < 4; ++i){
            auto child = make_shared<quadtree>(length / 2.0, x + length / 2.0 * (i & 1),
                                               y + length / 2.0 * ((i >> 1) & 1));
            child->parent = this;
            children.push_back(child);
        };

        // points are not located at edges
        for (auto point : points){
            auto index = int((point->x - x) / (length / 2.0)) |
                    ((int((point->y - y) / (length / 2.0))<< 1));
            children[index]->addPoint(point);
        };

        // recursively populate children quadtree
        for (auto child : children) {
            child->updateAttribute();
            child->populate();
        }

    }
}

void quadtree::setStatus(Status _status) noexcept  {status = _status; }

void quadtree::updateAttribute() noexcept {
    attribute = 0.;
    for (auto point : points) {
        attribute += point->attribute;
    }
    attribute /= points.size();
}