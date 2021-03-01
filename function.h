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

struct EqPair{
    bool operator()(const std::pair<Vertex, Vertex> &p1, const std::pair<Vertex, Vertex> &p2) const
    {return p1.first == p2.first && p1.second == p2.second;}
};

struct HashPair{
    size_t operator()(const std::pair<Vertex, Vertex> &key) const
    {return std::hash<Vertex>()(key.first) ^ std::hash<Vertex>()(key.second);}
};

using PairIBD = std::map<std::pair<Vertex, Vertex>, Pair*>;
//using PairIBD = std::unordered_map<std::pair<Vertex, Vertex>, Pair*, HashPair, EqPair>;

int getRelfromK(double K, int maxDeg);

void build_graph(Pedigree &pedigree, PairIBD &allsegs, const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results, std::map<Vertex, Vertex> &twins, 
    const std::map<std::string, int> &id2index, double tot_genome, double bkg_sharing, int maxDeg);

bool is_avunc(const std::string &fs1, const std::string &fs2, const std::string &avunc, 
        const PairIBD &allsegs, const std::map<std::string, int> &id2index,
        const std::map<std::string, std::map<int, double>*> &snpmap);

bool checkAvunc(const std::vector<Vertex> &full_sibs, const Vertex avunc, 
    const PairIBD &allsegs, const std::map<std::string, int> &id2index, 
    const std::map<std::string, std::map<int, double>*> &snpmap);

void run_druid(Pedigree &pedigree, const PairIBD &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results,
    const std::map<std::string, int> &id2index,
    FileOrGZ<FILE *> &logFile,
    double tot_genome, double bkg_sharing, int maxDeg);

struct ConnInfo{
    // this struct keeps track of the information we need for combining segments for every non-trivial (size > 1) connected components
    // all relationships are from the point view of the full-sib generation
    std::vector<Vertex> gp1; // grandparents if any
    std::vector<Vertex> gp2;
    std::vector<Vertex> av1; // AV if any
    std::vector<Vertex> av2;
    std::vector<Vertex> p; // parents if any
    std::vector<Vertex> fs; // full siblings
};

void postorder(const std::vector<Vertex> &components, const Pedigree &pedigree, std::vector<Vertex> &ordered);
void preorder(const std::vector<Vertex> &components, const Pedigree &pedigree, std::vector<Vertex> &ordered);
void grabCloseRelatives_o(const Vertex u, ConnInfo &con, const Pedigree &pedigree);
void grabCloseRelatives(const Vertex &u, ConnInfo &con, const Pedigree &pedigree);
bool isFS2Everyone(const Vertex &u, const std::vector<Vertex> &fs, const Pedigree &pedigree);
void combineIBD(const ConnInfo &con1, const ConnInfo &con2, std::unordered_set<Vertex> &checked1,
    std::unordered_set<Vertex> &checked2, std::map<std::pair<Vertex, Vertex>, int> &results,
    double bkg_sharing, int max_deg);

bool isSingleton(const ConnInfo &con);
void printConnInfo(const ConnInfo &con, const Pedigree &pedigree); // for debugging

double UnionIbdOverTwoSets(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2,
    const std::map<std::string, int> &id2index, const PairIBD &allsegs); // return the total length of combined IBD

int oneVSpedigree(Vertex u, const ConnInfo &con, std::unordered_set<Vertex> &visited,
    const PairIBD &allsegs, std::map<std::pair<Vertex, Vertex>, int> &results,
    const Pedigree &pedigree, const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int max_deg);

std::pair<int, int> pedigreeVSpedigree(const ConnInfo &con1, const ConnInfo &con2,
    std::unordered_set<Vertex> &visited1, std::unordered_set<Vertex> &visited2,
    const PairIBD &allsegs, const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results,
    const Pedigree &pedigree, const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg);

std::pair<bool, Vertex> polarizeUnpolarPC(Vertex v1, int d1, double k1, Vertex v2, int d2, double k2);

void PCpairVSone(Vertex d, const std::pair<Vertex, Vertex> &pc, 
    const PairIBD &allsegs, std::map<std::pair<Vertex, Vertex>, int> &results, int maxDeg);

int PCpairVSpedigree(const std::pair<Vertex, Vertex> &pc, const ConnInfo &con,
    std::unordered_set<Vertex> &visited, const PairIBD &allsegs,
    std::map<std::pair<Vertex, Vertex>, int> &results, const Pedigree &pedigree,
    const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg);

void updateSibsetByTheirParent(int index_p, 
    const std::vector<Vertex> &fs, const std::vector<Vertex> &parents,
    const ConnInfo &con2,
    std::unordered_set<Vertex> &visited1, std::unordered_set<Vertex> &visited2,
    const PairIBD &allsegs, std::map<std::pair<Vertex, Vertex>, int> &results,
    const Pedigree &pedigree, const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg);
// the abvoe function infers fs's relationship to con2 by using fs's parent.

void updateSibsetByTheirGrandParent(int index_gp, int index_av,
    const ConnInfo &con1, const ConnInfo &con2,
    std::unordered_set<Vertex> &visited1, std::unordered_set<Vertex> &visited2,
    const PairIBD &allsegs, std::map<std::pair<Vertex, Vertex>, int> &results,
    const Pedigree &pedigree, const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg);
// the abvoe function infers samples in con1's relationship to con2 by using con1's grandparent.



inline int resetRelationship(int ref, int offset, int maxDeg)
{
    int deg = ref == -1 ? -1 : ref + offset;
    if (deg > maxDeg){deg = -1;}
    return deg;
}

void setRelationshipBetweenTwoSets(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2, 
    std::map<std::pair<Vertex, Vertex>, int> &results, int deg);

void setRelationshipBetweenOneSampleAndSet(Vertex u, const std::vector<Vertex> &set,
    std::map<std::pair<Vertex, Vertex>, int> &results, int deg);

inline void setDeg(Vertex u, Vertex v, int deg, std::map<std::pair<Vertex, Vertex>, int> &results)
{
    if (deg == -1){return;}
    else{results[make_pair_v(u, v)] = deg;}
}

inline int getDeg(Vertex u, Vertex v, const std::map<std::pair<Vertex, Vertex>, int> &results)
{
    auto it = results.find(make_pair_v(u, v));
    if (it == results.end()){return -1;}
    else{return it->second;}
}

double averageKinship(Vertex u, const std::vector<Vertex> &set, 
    const std::map<std::pair<Vertex, Vertex>, Pair*> &allsegs);

double averageKinshipBetweenTwoSets(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2,
    const std::map<std::pair<Vertex, Vertex>, Pair*> &allsegs);

bool includeAunts(const std::vector<Vertex> &aunts, const std::vector<Vertex> &sibs,
    Vertex d, const PairIBD &allsegs);

double minKinshipBetweenTwoSibset(const std::vector<Vertex> &sib1, 
    const std::vector<Vertex> &sib2, const PairIBD &allsegs);

double maxKinshipBetweenTwoSibset(const std::vector<Vertex> &sib1,
    const std::vector<Vertex> &sib2, const PairIBD &allsegs);

int whichAV2Include(const std::vector<Vertex> &av11, const std::vector<Vertex> &av12, 
    const std::vector<Vertex> &sib1, const std::vector<Vertex> &sib2, const PairIBD &allsegs);

bool includeAunts(const std::vector<Vertex> &aunts, const std::vector<Vertex> &sibs, double min_ks1s2, const PairIBD &allsegs);

int whichParent2Include(const std::vector<Vertex> &parents, 
    const std::vector<Vertex> &sib1, const std::vector<Vertex> &sib2, const PairIBD &allsegs);

bool includeParent(Vertex p, double max_ks1s2, const std::vector<Vertex> &sibs, const PairIBD &allsegs);

double IBD0011(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2,
    const std::map<std::string, std::map<int, double>*> &snpmap, 
    const PairIBD &allsegs, const std::map<std::string, int> &id2index);

void grabChildren(Vertex p, std::vector<Vertex> &children, const Pedigree &pedigree);

void propagate(Vertex u, Vertex v, std::unordered_set<Vertex> &visited1, 
    std::unordered_set<Vertex> &visited2, int maxDeg, const Pedigree &pedigree, 
    std::map<std::pair<Vertex, Vertex>, int> &results);

void propagateAlongPedigree(const ConnInfo &con1, const ConnInfo &con2, int av1, int av2,
    std::unordered_set<Vertex> &visited1, std::unordered_set<Vertex> &visited2,
    const Pedigree &pedigree, std::map<std::pair<Vertex, Vertex>, int> &results, int maxDeg);

inline double getTg(int numAV, int numSib){
    if(numAV == 0){return 1.0 - pow(0.5, numSib);}
    else{
        return 1.0 - pow(0.5, numAV) + pow(0.5, numAV+1)*(1.0 - pow(0.5, numSib));
    }
}

inline double calc_bkg_sharing(int numAV, int numSibs, double bkg_sharing)
{
    return 2.0*bkg_sharing*(1.0 - pow(0.5, numAV) + pow(0.5, numAV+1)*(1.0 - pow(0.5, numSibs))) + bkg_sharing*(1.0 - pow(0.5, numSibs)); 
}

inline double calc_bkg_sharing(int numAV1, int numSib1, int numAV2, int numSib2, double bkg_sharing)
{
    return 4.0*bkg_sharing*(1.0 - pow(0.5, numAV1) + pow(0.5, numAV1+1)*(1.0 - pow(0.5, numSib1)))
                    *(1.0 - pow(0.5, numAV2) + pow(0.5, numAV2+1)*(1.0 - pow(0.5, numSib2))) + 
                    2.0*bkg_sharing*(1.0 - pow(0.5, numAV1) + pow(0.5, numAV1+1)*(1.0 - pow(0.5, numSib1)))*(1.0 - pow(0.5, numSib2)) + 
                    2.0*bkg_sharing*(1.0 - pow(0.5, numAV2) + pow(0.5, numAV2+1)*(1.0 - pow(0.5, numSib2)))*(1.0 - pow(0.5, numSib1)) + 
                    bkg_sharing*(1.0 - pow(0.5, numSib1))*(1.0 - pow(0.5, numSib2));
}

std::pair<bool, std::size_t> findUnpolarizedPC(Vertex u, const Pedigree &pedigree); 


bool isGP(Vertex u, const Pedigree &pedgiree);
bool isAV(Vertex u, const Pedigree &pedgiree);
bool isP(Vertex u, const Pedigree &pedgiree);
bool isFS(Vertex u, const Pedigree &pedgiree);


void write_output(const std::map<std::pair<Vertex, Vertex>, int> &results, 
    const std::string &prefix, const Pedigree &pedigree, const std::map<Vertex, Vertex> &twins);

#endif
