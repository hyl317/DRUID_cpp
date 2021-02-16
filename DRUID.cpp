#include <string>
#include <chrono>
#include "tools.h"
#include "function.h"
#include "function_t.h"

int main(int argc, char **argv){

    if (argc <= 1){
        print_help();
    }

    std::string ibdFile, bimFile, NeFile;
    std::string prefix;
    std::string exSamples;
    int maxDeg = 11;
    int threads = 1;
    double minIBD = 2.0;
    parse_command_line(argc, argv, ibdFile, bimFile, NeFile, exSamples, prefix, maxDeg, threads, minIBD);

    std::string logFileName = prefix + ".log";
    FileOrGZ<FILE *> logFile;
    bool ret = logFile.open(logFileName.c_str(), "w");
    if(!ret){
        fprintf(stderr, "cannot open %s for writing output\n", logFileName.c_str());
        return 1;
    }

    logFile.printf("Running DRUID...\n");

    logFile.printf("Reading bimFile: %s\n", bimFile.c_str());
    auto snpmap = std::map<std::string, std::map<int, double>*>();
    Eigen::VectorXd chrLens = readBimFile(bimFile, snpmap, logFile);

    logFile.printf("Reading IBD segment file: %s\n", ibdFile.c_str());
    auto t1 = std::chrono::high_resolution_clock::now();
    PairIBD allsegs;
    std::map<std::string, Vertex> id2Vertex;
    boost::object_pool<Pair> p_pair;
    boost::object_pool<ibdMapType> p_ibdmap;
    if (exSamples.length() == 0){
        readIBDFile(ibdFile, allsegs, id2Vertex, p_pair, p_ibdmap);
    }else{
        // exclude some samples
        readIBDFile_ex(ibdFile, allsegs, id2Vertex, exSamples, logFile);
    }
    int numSample = id2Vertex.size();
    auto t2 = std::chrono::high_resolution_clock::now();
    auto d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    logFile.printf("\tFinished reading segments from %d samples for analysis, takes %lfs\n", numSample, d/1e6);

    double bkg_sharing = 0.0;
    if (NeFile.length() > 0){
        logFile.printf("Correcting for recent demography...\n");
        logFile.printf("\tReading Ne File: %s\n", NeFile.c_str());
        if (minIBD != 0.0){logFile.printf("\tminimum IBD length detected: %lf\n", minIBD);}
        else{
            logFile.printf("\tMinimum IBD length not provided. Use default value 2.0cM.\n");
            minIBD = 2.0;
        }
        bkg_sharing = calc_bkg_sharing(NeFile, chrLens, minIBD);
        logFile.printf("\tExpected background sharing: %lf\n", bkg_sharing);
    }

    logFile.printf("Maximum Relatedness Reported: degree %d\n", maxDeg);
    logFile.printf("Identifying clusters of close relatives...\n");
    t1 = std::chrono::high_resolution_clock::now();
    // make a graph and add vertices properties to it
    Pedigree pedigree = Pedigree(numSample);
    for(auto it = id2Vertex.begin(); it != id2Vertex.end(); it++){
        pedigree[it->second].id = it->first;
    }

    // add edges between close relatives
    std::map<std::pair<Vertex, Vertex>, int> results;
    std::map<Vertex, Vertex> twins;
    build_graph(pedigree, allsegs, snpmap, results, twins, chrLens.sum(), bkg_sharing, maxDeg);
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    logFile.printf("Building graph done, takes %lfs\n", d/1e6);
    t1 = std::chrono::high_resolution_clock::now();
    if (threads == 1){
        run_druid(pedigree, allsegs, snpmap, results, logFile, chrLens.sum(), bkg_sharing, maxDeg);
    }else{
        run_druid_t(pedigree, allsegs, snpmap, results, threads, logFile, chrLens.sum(), bkg_sharing, maxDeg);
    }
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    logFile.printf("Anaylzing pairwise connected components done, takes %lfs\n", d/1e6);
    
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
    logFile.printf("Writitng to output file: %s\n", std::string(prefix + ".DRUID").c_str());
    t1 = std::chrono::high_resolution_clock::now();
    write_output(results, prefix, pedigree, twins);
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    logFile.printf("Writing to output done, takes %lfs\n", d/1e6);

    // clean up
    t1 = std::chrono::high_resolution_clock::now();
    for(auto it = snpmap.begin(); it != snpmap.end(); it++){
        delete it->second;
    }

    for(auto it = allsegs.begin(); it != allsegs.end(); it++){
        Pair *p = it->second;
        for (auto it2 = p->ibd1_map->begin(); it2 != p->ibd1_map->end(); it2++){
            delete it2->second;
        }
        for(auto it3 = p->ibd2_map->begin(); it3 != p->ibd2_map->end(); it3++){
            delete it3->second;
        }
	    delete p->ibd1_map;
	    delete p->ibd2_map;
        delete p;

    }
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    logFile.printf("Clean up heap space takes %lfs\n", d/1e6);
    logFile.close();


    return 0;
}
