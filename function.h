#ifndef FUNCTION_H
#define FUNCTION_H

#include <algorithm>
#include <string>
#include <math.h>
#include <boost/graph/adjacency_list.hpp>
#include "tools.h"


enum close_relationship {PC, FS, GP, AV};

struct sample{
    std::string id;
    // could add something else, like age, sex if available
};

struct relationship{
    close_relationship rel;
    std::size_t older;
    bool polarized = false; 
    // could add something else later
};

// use std::list to store the list of vertices for fast add and removal of vertices
// use std::vector to store the out-going edge for each vertex for fast traversal and small memory overhead
// bidirectionalS makes the graph to have both in_edges() and out_edges() functions, not sure this is useful or not yet
using Pedigree = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS,
    sample, relationship>;
using Vertex = boost::graph_traits<Pedigree>::vertex_descriptor;
using Edge = boost::graph_traits<Pedigree>::edge_descriptor;

inline std::pair<Vertex, Vertex> make_pair_v
(Vertex u, Vertex v){
  return u < v? std::make_pair(u, v) : std::make_pair(v, u);
}

int getRelfromK(double ibd1, double ibd2, double bkg, double tot_genome, int maxDeg);

void build_graph(Pedigree &pedigree, 
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results,
    std::map<Vertex, Vertex> &twins,
    double tot_genome, double bkg_sharing, int maxDeg);

bool is_avunc(const std::string &fs1, const std::string &fs2, const std::string &avunc, 
        const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
        const std::map<std::string, std::map<int, double>*> &snpmap);

bool checkAvunc(const std::vector<std::string> &full_sibs, const std::string &avunc, 
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap);

void run_druid(Pedigree &pedigree, 
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results,
    FileOrGZ<FILE *> &logFile,
    double tot_genome, double bkg_sharing, int maxDeg);

struct ConnInfo{
    // this struct keeps track of the information we need for combining segments for every non-trivial (size > 1) connected components
    // all relationships are from the point view of the full-sib generation
    std::vector<Vertex> gp; // grandparents if any
    std::vector<Vertex> av; // AV if any
    std::vector<Vertex> p; // parents if any
    std::vector<Vertex> fs; // full siblings
};

bool isGP(Vertex u, const Pedigree &pedgiree);
bool isAV(Vertex u, const Pedigree &pedgiree);
bool isP(Vertex u, const Pedigree &pedgiree);
bool isFS(Vertex u, const Pedigree &pedgiree);


void write_output(const std::map<std::pair<Vertex, Vertex>, int> &results, 
    const std::string &prefix, const Pedigree &pedigree, const std::map<Vertex, Vertex> &twins);

#endif