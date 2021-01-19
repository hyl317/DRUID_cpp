#include <string>
#include "tools.h"
#include "function.h"

int main(int argc, char **argv){

    // test interval union
    // ibdSegments set1;
    // set1.push_back(std::make_pair(0.079263, 17.0939));
    // set1.push_back(std::make_pair(36.9644, 53.8065));
    // set1.push_back(std::make_pair(54.3264, 80.3289));
    // set1.push_back(std::make_pair(83.5778, 197.2790));

    // ibdSegments set2;
    // set2.push_back(std::make_pair(36.9644, 53.8065));
    // set2.push_back(std::make_pair(54.3264, 83.6337));
    // set2.push_back(std::make_pair(194.4780, 203.9290));

    // ibdSegments set3;
    // interval_union(set1, set2, set3);
    // std::for_each(set3.begin(), set3.end(), [&](const std::pair<double, double> interval){fprintf(stdout, "[%lf, %lf]\n", interval.first, interval.second);});
    // return 0;

    //end of test

    if (argc <= 1){
        print_help();
    }

    std::string ibdFile, bimFile, NeFile;
    std::string prefix;
    std::string exSamples;
    int maxDeg = 11;
    double minIBD = 2.0;
    parse_command_line(argc, argv, ibdFile, bimFile, NeFile, exSamples, prefix, maxDeg, minIBD);

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
    auto allsegs = std::map<std::pair<std::string, std::string>, Pair*>();
    std::set<std::string> inds;
    if (exSamples.length() == 0){
        readIBDFile(ibdFile, allsegs, inds);
    }else{
        // exclude some samples
        readIBDFile_ex(ibdFile, allsegs, inds, exSamples, logFile);
    }
    int numSample = inds.size();
    logFile.printf("\tFinished reading segments from %d samples for analysis\n", numSample);

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
    // make a graph and add vertices properties to it
    Pedigree pedigree = Pedigree(numSample);
    auto pair = boost::vertices(pedigree);
    Pedigree::vertex_descriptor v_head = *(pair.first);
    Pedigree::vertex_descriptor v_tail = *(pair.second);
    auto it2 = inds.begin();
    for(auto it = v_head; it != v_tail; it++){
        assert(it2 != inds.end());
        pedigree[it].id = *it2;
        it2++;
    }
    assert(it2 == inds.end());

    auto vertex_property_map = boost::get(&sample::id, pedigree);
    // make a map from sampleID(string) to Vertex
    std::map<std::string, Vertex> id2Vertex;
    boost::graph_traits<Pedigree>::vertex_iterator vi, vi_end;
    for(boost::tie(vi, vi_end) = boost::vertices(pedigree); vi != vi_end; vi++){
        id2Vertex.insert(std::make_pair(vertex_property_map[*vi], *vi));
        //fprintf(stdout, "%d: %s\n", *vi, vertex_property_map[*vi].c_str());
    }

    std::map<std::pair<Vertex, Vertex>, Pair*> allsegs_v;
    for(auto it = allsegs.begin(); it != allsegs.end(); it++){
        std::pair<std::string, std::string> p = it->first;
        assert(id2Vertex.find(p.first) != id2Vertex.end());
        assert(id2Vertex.find(p.second) != id2Vertex.end());
        Vertex u = id2Vertex[p.first];
        Vertex v = id2Vertex[p.second];
        allsegs_v[make_pair_v(u,v)] = it->second;
    }

    // add edges between close relatives
    std::map<std::pair<Vertex, Vertex>, int> results;
    std::map<Vertex, Vertex> twins;
    build_graph(pedigree, allsegs_v, snpmap, results, twins, chrLens.sum(), bkg_sharing, maxDeg);
    run_druid(pedigree, allsegs_v, snpmap, results, logFile, chrLens.sum(), bkg_sharing, maxDeg);

    // for testing AV detection
    // is_avunc("801120", "801121", "801113", allsegs, snpmap);
    // is_avunc("801120", "801121", "801114", allsegs, snpmap);
    // is_avunc("801120", "801121", "801115", allsegs, snpmap);
    // is_avunc("801120", "801121", "801118", allsegs, snpmap);

    // is_avunc("801120", "801122", "801113", allsegs, snpmap);
    // is_avunc("801120", "801122", "801114", allsegs, snpmap);
    // is_avunc("801120", "801122", "801115", allsegs, snpmap);
    // is_avunc("801120", "801122", "801118", allsegs, snpmap);

    // is_avunc("801121", "801122", "801113", allsegs, snpmap);
    // is_avunc("801121", "801122", "801114", allsegs, snpmap);
    // is_avunc("801121", "801122", "801115", allsegs, snpmap);
    // is_avunc("801121", "801122", "801118", allsegs, snpmap);
    // end of test

    // test UnionIBDover2sets
    // std::vector<Vertex> set1;
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b1-i1"]);
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b2-i1"]);
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b3-i1"]);
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b4-i1"]);
    // set1.push_back(id2Vertex["ped2_D3_1_g3-b5-i1"]);

    // std::vector<Vertex> set2;
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b6-i1"]);
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b7-i1"]);
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b8-i1"]);
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b9-i1"]);
    // set2.push_back(id2Vertex["ped2_D3_1_g3-b10-i1"]);

    // double tmp = UnionIbdOverTwoSets(set1, set2, allsegs_v);
    // fprintf(stdout, "combined ibd1 length: %lf\n", tmp);
    // double t1 = getTg(0, 5);
    // double t2 = getTg(0, 5);
    // double k1 = (tmp/chrLens.sum())/(t1*t2);
    // fprintf(stdout, "k1 is %lf\n", k1);
    // double k2 = (IBD0011(set1, set2, snpmap, allsegs_v)/chrLens.sum())/(t1*t2);
    // fprintf(stdout, "k2 is %lf\n", k2);
    // int deg = getRelfromK(k1/4.0, maxDeg);
    // fprintf(stdout, "Unioned IBD length for ped2_D3: %lf, estimated deg is %d\n", tmp, deg);


    // for (int i = 0; i < 100; i++){
    //    std::random_shuffle(set1.begin(), set1.end());
    //    std::random_shuffle(set2.begin(), set2.end());
    //    double tmp = UnionIbdOverTwoSets(set1, set2, allsegs_v);
    //    fprintf(stdout, "Unioned IBD length: %lf\n", tmp);
    // }
    // end of test
    
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
    //logFile.printf("\tNumber of edges: %d\n", boost::num_edges(pedigree));
    write_output(results, prefix, pedigree, twins);
    logFile.close();

    // clean up
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

    return 0;
}
