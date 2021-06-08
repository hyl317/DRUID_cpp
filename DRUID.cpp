#include <string>
#include <chrono>
#include "tools.h"
#include "function.h"
#include "function_t.h"

int main(int argc, char **argv){

    if (argc <= 1){
        print_help();
    }
    
    std::string segFile, ibd12, bimFile, NeFile;
    std::string prefix;
    std::string exSamples;
    int maxDeg = 11;
    int threads = 1;
    double minIBD = 2.0;
    double bkg_sharing = 0.0;
    parse_command_line(argc, argv, segFile, ibd12, bimFile, NeFile, exSamples, prefix, maxDeg, threads, minIBD, bkg_sharing);

    std::string logFileName = prefix + ".log";
    FILE *logFile;
    logFile = fopen(logFileName.c_str(), "w");
    if(logFile == nullptr){
        fprintf(stderr, "cannot open %s for writing output\n", logFileName.c_str());
        return 1;
    }

    fprintf(logFile, "Running DRUID...\n");

    fprintf(logFile, "Reading bimFile: %s\n", bimFile.c_str());
    std::map<std::string, std::map<int, double>*> snpmap;
    chromMap id2index;
    Eigen::VectorXd chrLens = readBimFile(bimFile, snpmap, id2index, logFile);

    if (NeFile.length() > 0){
        fprintf(logFile, "Correcting for recent demography...\n");
        fprintf(logFile, "\tReading Ne File: %s\n", NeFile.c_str());
        if (minIBD != 0.0){fprintf(logFile, "\tminimum IBD length detected: %lf\n", minIBD);}
        else{
            fprintf(logFile, "\tMinimum IBD length not provided. Use default value 2.0cM.\n");
            minIBD = 2.0;
        }
        bkg_sharing = calc_bkg_sharing(NeFile, chrLens, minIBD);
        fprintf(logFile, "\tExpected background sharing: %lf\n", bkg_sharing);
    }else{
        fprintf(logFile, "Using background sharing value: %lf\n", bkg_sharing);
    }

    fflush(logFile);
    auto t1 = std::chrono::high_resolution_clock::now();
    Pedigree pedigree;
    PairIBD allsegs;
    std::map<char *, Vertex, cmp_str> id2Vertex;
    std::map<std::pair<Vertex, Vertex>, int> results;
    std::map<Vertex, Vertex> twins;
    // read input
    readInput(ibd12, segFile, exSamples, pedigree, allsegs, id2Vertex, id2index, 
        twins, results, snpmap, bkg_sharing, chrLens.sum(), maxDeg, logFile);
    int numSample = id2Vertex.size();
    auto t2 = std::chrono::high_resolution_clock::now();
    auto d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    fprintf(logFile, "Finished reading inputs and building graph from %d samples for analysis, takes %lfs\n", numSample, d/1e6);
    fprintf(logFile, "Maximum Relatedness Reported: degree %d\n", maxDeg);
    fprintf(logFile, "Identifying clusters of close relatives...\n");
    fflush(logFile);

    t1 = std::chrono::high_resolution_clock::now();
    if (threads == 1){
        run_druid(pedigree, allsegs, snpmap, results, id2index, logFile, chrLens.sum(), bkg_sharing, maxDeg);
    }else{
        // TODO
        run_druid_t(pedigree, allsegs, snpmap, results, id2index, threads, logFile, chrLens.sum(), bkg_sharing, maxDeg);
    }
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    fprintf(logFile, "Anaylzing pairwise connected components done, takes %lfs\n", d/1e6);
    
    // testing IBD0011
    // std::vector<Vertex> set1;
    // set1.push_back(id2Vertex["ped4_D5_1_g3-b1-i1"]);
    // set1.push_back(id2Vertex["ped4_D5_1_g3-b2-i1"]);
    // set1.push_back(id2Vertex["ped4_D5_1_g3-b3-i1"]);

    // std::vector<Vertex> set2;
    // set2.push_back(id2Vertex["ped4_D5_1_g3-b4-i1"]);
    // set2.push_back(id2Vertex["ped4_D5_1_g3-b5-i1"]);
    // set2.push_back(id2Vertex["ped4_D5_1_g3-b6-i1"]);
    // double tmp = IBD0011(set1, set2, snpmap, allsegs_v);
    // fprintf(stdout, "ibd0011: %lf\n", tmp);

    // //another test
    // set1.clear();
    // set2.clear();
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b1-i1"]);
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b2-i1"]);
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b3-i1"]);
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b4-i1"]);
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b5-i1"]);

    // set2.push_back(id2Vertex["ped2_D3_1_g3-b6-i1"]);
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b7-i1"]);
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b8-i1"]);
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b9-i1"]);
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b10-i1"]);
    // tmp = IBD0011(set1, set2, snpmap, allsegs_v);
    // fprintf(stdout, "ibd0011: %lf\n", tmp);

    // // another test, this set should have nearly no IBD0011 region because their ungenotyped parents are first-cousin
    // set1.clear();
    // set2.clear();
    // set1.push_back(id2Vertex["ped2_D5_1_g4-b1-i1"]);
    // set1.push_back(id2Vertex["ped2_D5_1_g4-b2-i1"]);
    // set1.push_back(id2Vertex["ped2_D5_1_g4-b3-i1"]);
    // set1.push_back(id2Vertex["ped2_D5_1_g4-b4-i1"]);
    // set1.push_back(id2Vertex["ped2_D5_1_g4-b5-i1"]);

    // set2.push_back(id2Vertex["ped2_D5_1_g4-b6-i1"]);
    // set2.push_back(id2Vertex["ped2_D5_1_g4-b7-i1"]);
    // set2.push_back(id2Vertex["ped2_D5_1_g4-b8-i1"]);
    // set2.push_back(id2Vertex["ped2_D5_1_g4-b9-i1"]);
    // set2.push_back(id2Vertex["ped2_D5_1_g4-b10-i1"]);
    // tmp = IBD0011(set1, set2, snpmap, allsegs_v);
    // fprintf(stdout, "ibd0011: %lf\n", tmp);

    // // another test, this set should alaso have nearly no IBD0011 region
    // set1.clear();
    // set2.clear();
    // set1.push_back(id2Vertex["ped2_D4_1_g3-b2-i1"]);
    // set1.push_back(id2Vertex["ped2_D4_1_g3-b3-i1"]);
    // set1.push_back(id2Vertex["ped2_D4_1_g3-b4-i1"]);
    // set1.push_back(id2Vertex["ped2_D4_1_g3-b5-i1"]);
    // set1.push_back(id2Vertex["ped2_D4_1_g3-b6-i1"]);

    // set2.push_back(id2Vertex["ped2_D4_1_g4-b1-i1"]);
    // set2.push_back(id2Vertex["ped2_D4_1_g4-b2-i1"]);
    // set2.push_back(id2Vertex["ped2_D4_1_g4-b3-i1"]);
    // set2.push_back(id2Vertex["ped2_D4_1_g4-b4-i1"]);
    // set2.push_back(id2Vertex["ped2_D4_1_g4-b5-i1"]);
    // tmp = IBD0011(set1, set2, snpmap, allsegs_v);
    // fprintf(stdout, "ibd0011: %lf\n", tmp);
    // end of test



    // writing output, finishing up
    fprintf(logFile, "Writitng to output file: %s\n", std::string(prefix + ".DRUID").c_str());
    fflush(logFile);
    t1 = std::chrono::high_resolution_clock::now();
    write_output(results, prefix, pedigree, twins);
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    fprintf(logFile, "Writing to output done, takes %lfs\n", d/1e6);

    // clean up
    for(auto it = snpmap.begin(); it != snpmap.end(); it++){
        delete it->second;
    }

    // for(auto it = allsegs.begin(); it != allsegs.end(); it++){
    //     Pair *p = it->second;
    //     for (auto it2 = p->ibd1_map->begin(); it2 != p->ibd1_map->end(); it2++){
    //         delete it2->second;
    //     }
    //     for(auto it3 = p->ibd2_map->begin(); it3 != p->ibd2_map->end(); it3++){
    //         delete it3->second;
    //     }
	//     delete p->ibd1_map;
	//     delete p->ibd2_map;
    //     delete p;
    // }
    fclose(logFile);


    return 0;
}
