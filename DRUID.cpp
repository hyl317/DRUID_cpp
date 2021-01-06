#include <string>
#include "tools.h"
#include "function.h"

int main(int argc, char **argv){
    if (argc <= 1){
        print_help();
    }

    std::string ibdFile, bimFile, NeFile;
    std::string prefix;
    int maxDeg = 11;
    double minIBD = 0.0;
    parse_command_line(argc, argv, ibdFile, bimFile, NeFile, prefix, maxDeg, minIBD);

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
    readIBDFile(ibdFile, allsegs, inds);
    int numSample = inds.size();
    logFile.printf("\tFinished reading segments for %d samples\n", numSample);

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
    build_graph(pedigree, allsegs, snpmap, results, twins, chrLens.sum(), bkg_sharing, maxDeg);
    run_druid(pedigree, allsegs, snpmap, results, logFile, chrLens.sum(), bkg_sharing, maxDeg);

    // for testing purpose
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

    // writing output, finishing up
    logFile.printf("Writitng to output file: %s\n", std::string(prefix + ".DRUID").c_str());
    logFile.printf("\tNumber of edges: %d\n", boost::num_edges(pedigree));
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
