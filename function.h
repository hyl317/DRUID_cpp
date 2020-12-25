#ifndef FUNCTION_H
#define FUNCTION_H

#include <algorithm>
#include <string>
#include <math.h>
#include <boost/graph/adjacency_list.hpp>
#include "tools.h"

int getRelfromK(double ibd1, double ibd2, double bkg, double tot_genome, int maxDeg);

struct sample{
    std::string id;
    // could add something else, like age, sex if available
};

enum close_relationship {PC, FS, GP, AV};
struct relationship{
    close_relationship rel;
    Vertex older = NULL;
    // could add something else later
};

// use std::list to store the list of vertices for fast add and removal of vertices
// use std::vector to store the out-going edge for each vertex for fast traversal and small memory overhead
// bidirectionalS makes the graph to have both in_edges() and out_edges() functions, not sure this is useful or not yet
using Pedigree = boost::adjacency_list<
    boost::listS, boost::vecS, boost::undirectedS,
    sample, relationship>;
using Vertex = boost::graph_traits<Pedigree>::vertex_descriptor;
using Edge = boost::graph_traits<Pedigree>::edge_descriptor;

void build_graph(Pedigree &pedigree, 
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<std::string, std::string>, int> &results,
    double tot_genome, double bkg_sharing, int maxDeg);

bool is_avunc(const std::string &fs1, const std::string &fs2, const std::string &avunc, 
        const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
        const std::map<std::string, std::map<int, double>*> &snpmap);


#endif