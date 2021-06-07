#include <thread>
#include "function_t.h"
#include <boost/graph/connected_components.hpp>


void run_druid_t(Pedigree &pedigree, const PairIBD &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results, const chromMap &id2index,
    int numThread, FILE *logFile, double tot_genome, double bkg_sharing, int maxDeg)
{
    fprintf(logFile, "Start the primary DRUID algorithm with %d threads...\n", numThread);
    std::vector<int> components(boost::num_vertices(pedigree));
    int num_components = boost::connected_components(pedigree, &components[0]);
    fprintf(logFile, "\tnumber of connected components: %d\n", num_components);

    std::map<int, std::shared_ptr<std::vector<Vertex>>> comp_map;
    boost::graph_traits<Pedigree>::vertex_iterator vi, vi_end;
    for(boost::tie(vi, vi_end) = boost::vertices(pedigree); vi != vi_end; vi++){
        int comp_index = components[*vi];
        if (comp_map.find(comp_index) == comp_map.end()){
            comp_map.insert(std::make_pair(comp_index, std::shared_ptr<std::vector<Vertex>>(new std::vector<Vertex>())));
        }
        comp_map[comp_index]->push_back(*vi);
    }

    for(int i = 0; i < num_components; i++){
        auto ordered = std::shared_ptr<std::vector<Vertex>>(new std::vector<Vertex>());
        postorder(*(comp_map.find(i)->second), pedigree, *ordered);
        comp_map[i] = ordered;
    }

    // using more than 1 threads to analyze pairs of connected components
    // each thread maintain their own updated "results" map
    std::map<int, std::map<std::pair<Vertex, Vertex>, int>*> updated_results;
    for(int i = 0; i < numThread; i++){
        updated_results.insert(std::make_pair(i, new std::map<std::pair<Vertex, Vertex>, int>()));
    }

    std::vector<std::thread> vecThreads;
    for(int i = 0; i < numThread; i++){
        std::thread t(processPairs, i, numThread, std::ref(pedigree), std::ref(comp_map), std::ref(allsegs), std::ref(snpmap), std::ref(id2index), std::ref(*updated_results[i]), tot_genome, bkg_sharing, maxDeg);
        vecThreads.push_back(std::move(t));
    }

    for(std::thread &t : vecThreads){
        if(t.joinable()){t.join();}
    }

    // update results
    for(auto it = updated_results.begin(); it != updated_results.end(); it++){
        for(auto it2 = it->second->begin(); it2 != it->second->end(); it2++){
            results[it2->first] = it2->second;
        }
        delete it->second;
    } 
}

void processPairs(int thread_index, int numThread,  const Pedigree &pedigree,
    const std::map<int, std::shared_ptr<std::vector<Vertex>>> &comp_map, const PairIBD &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap, const chromMap &id2index,
    std::map<std::pair<Vertex, Vertex>, int> &results, double tot_genome, double bkg_sharing, int maxDeg)
{   
    int num_components = comp_map.size();
    int offset = thread_index == 0 ? numThread : thread_index;
    for(int i = 0; i < num_components; i++){
        for(int j = i + offset; j < num_components; j += numThread){
            std::unordered_set<Vertex> visited1;
            for(Vertex u : *(comp_map.find(i)->second)){
                if (visited1.find(u) != visited1.end()){continue;}
                ConnInfo con1;
                grabCloseRelatives(u, con1, pedigree);
                bool isSingleton1 = isSingleton(con1);
                bool hasPC1 = false;
                Vertex w1;
                if (isSingleton1){std::tie(hasPC1, w1) = findUnpolarizedPC(u, pedigree);}
		        std::unordered_set<Vertex> visited2;
                //fprintf(stdout, "---------------------------s--------------------------------\n");
                //fprintf(stdout, "con1:\n");
                //printConnInfo(con1, pedigree);
                for(Vertex v : *(comp_map.find(j)->second)){
                    if(visited2.find(v) != visited2.end()){continue;}
                    ConnInfo con2;
                    grabCloseRelatives(v, con2, pedigree);
                    //fprintf(stdout, "con2:\n");
                    //printConnInfo(con2, pedigree);
                    // analyzing the two ConnInfo component
                    bool isSingleton2 = isSingleton(con2);
                    bool hasPC2 = false;
                    Vertex w2;
                    if (isSingleton2){std::tie(hasPC2, w2) = findUnpolarizedPC(v, pedigree);}
                    if (isSingleton1 && isSingleton2){
                        visited1.insert(u);
                        visited2.insert(v);
                        if (hasPC1){
                            auto pc = std::make_pair(u, w1);
                            PCpairVSone(v, pc, allsegs, results, maxDeg);
                            visited1.insert(w1);
                        }
                        if (hasPC2){
                            auto pc = std::make_pair(v, w2);
                            PCpairVSone(u, pc, allsegs, results, maxDeg);
                            visited2.insert(w2);
                        }
                    }
                    else if(isSingleton1 && !isSingleton2){
                        visited1.insert(u);
                        if (!hasPC1){oneVSpedigree(u, con2, visited2, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);}
                        else{
                            auto pc = std::make_pair(u, w1);
                            visited1.insert(w1);
                            PCpairVSpedigree(pc, con2, visited2, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
                        }
                    }else if(!isSingleton1 && isSingleton2){
                        visited2.insert(v);
                        if (!hasPC2){oneVSpedigree(v, con1, visited1, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);}
                        else{
                            auto pc = std::make_pair(v, w2);
                            PCpairVSpedigree(pc, con1, visited1, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
                        }
                    }else{
                        //printConnInfo(con1, pedigree);
                        //printConnInfo(con2, pedigree);
                        pedigreeVSpedigree(con1, con2, visited1, visited2, allsegs, snpmap, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
                    }

                }
                //fprintf(stdout, "--------------------------e---------------------------------\n");
            }
        }
    }
}
