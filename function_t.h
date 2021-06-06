#ifndef FUNCTION_T_H
#define FUNCTION_T_H

#include "function.h"

void run_druid_t(Pedigree &pedigree, const PairIBD &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results, const chromMap &id2index, 
    int numThread, FileOrGZ<FILE *> &logFile, double tot_genome, double bkg_sharing, int maxDeg);

void processPairs(int thread_index, int numThread,  const Pedigree &pedigree,
    const std::map<int, std::shared_ptr<std::vector<Vertex>>> &comp_map, const PairIBD &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap, const chromMap &id2index,
    std::map<std::pair<Vertex, Vertex>, int> &results, double tot_genome, double bkg_sharing, int maxDeg);

#endif