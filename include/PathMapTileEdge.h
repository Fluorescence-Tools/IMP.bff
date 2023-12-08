/**
 *  \file IMP/bff/PathMapTileEdge.h
 *  \brief Tile edges used in path search by PathMap
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PATHMAPTILEEDGE_H
#define IMPBFF_PATHMAPTILEEDGE_H

#include <IMP/bff/bff_config.h>

#include <vector>

#include <IMP/bff/PathMap.h>
#include <IMP/bff/PathMapTile.h>

IMPBFF_BEGIN_NAMESPACE

class PathMap;
class PathMapTile;


class PathMapTileEdge{

friend class PathMapTile;
friend class PathMap;


protected:

    int tile_idx; /// The tile the edge is pointing to
    float length; /// the path length / cost of going to the tile

public:

    /// Length of an edge (usually cartesian distance between tiles)
    float get_length() const{
        return length;
    }

    /*!
     *
     * @param edge_target
     * @param edge_cost
     */
    PathMapTileEdge(
            int edge_target = -1,
            float edge_cost = std::numeric_limits<float>::max()
    ) :
            tile_idx(edge_target), length(edge_cost){}
};


IMPBFF_END_NAMESPACE

#endif //IMPBFF_PATHMAPTILEEDGE_H