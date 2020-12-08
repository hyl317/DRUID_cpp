#include <string>
#include "tools.h"

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
    int numSample = readIBDFile(ibdFile, allsegs, logFile);

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
    // Identify close relatives to build connected components







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
        delete p;

    }

    return 0;
}