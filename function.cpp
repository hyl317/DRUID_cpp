#include <math.h>
#include <memory>
#include "function.h"
#include "constants.h"

int getRelfromK(double K, int maxDeg){
    //double K = std::max((ibd1/4.0 + ibd2/2.0 - bkg/4.0)/tot_genome, 0.0); 
    if (K == 0){return -1;}
    int deg = (int)(-log2(K) + 0.5) - 1; // 四舍五入
    if (deg <= maxDeg){return deg;}
    else{return -1;}
}

void build_graph(Pedigree &pedigree, 
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
    std::map<std::pair<std::string, std::string>, int> &results,
    double tot_genome, double bkg_sharing, int maxDeg)
{
    auto vertex_property_map = boost::get(&sample::id, pedigree);
    boost::graph_traits<Pedigree>::vertex_iterator vi1, vi_end1;
    boost::graph_traits<Pedigree>::vertex_iterator vi2, vi_end2;
    using Vertex = boost::graph_traits<Pedigree>::vertex_descriptor;
    std::map<Vertex, std::unique_ptr<std::vector<Vertex>>> first_degs;
    std::map<Vertex, std::unique_ptr<std::vector<Vertex>>> second_degs;
    for(boost::tie(vi1, vi_end1) = boost::vertices(pedigree); vi1 != vi_end1; vi1++){
        std::string id1 = vertex_property_map[*vi1];
        for(boost::tie(vi2, vi_end2) = boost::vertices(pedigree); vi2 != vi_end2; vi2++){
            std::string id2 = vertex_property_map[*vi2];
            if (id1 >= id2){continue;} // avoid analyzing pairs twice
            auto it = allsegs.find(std::make_pair(id1, id2)); 
            // some pairs may not have ibd segments, so need to check here
            if(it == allsegs.end()){continue;};
            Pair p = *(it->second);
            double ibd1 = p.ibd1_tot;
            double ibd2 = p.ibd2_tot;
            double K = std::max((ibd1/4.0 + ibd2/2.0 - bkg_sharing/4.0)/tot_genome, 0.0);
            int deg = getRelfromK(K, maxDeg);
            results.insert(std::make_pair(std::make_pair(id1, id2), deg));
            //std::cout << id1 << "\t" << id2 << "\t" << deg << std::endl;
            // store first and second degree pairs' sample names for later use
            if (deg == 1 || deg == 2){
                auto &map = deg == 1? first_degs : second_degs;
                if(map.find(*vi1) == map.end()){
                    map[*vi1] = std::unique_ptr<std::vector<Vertex>>(new std::vector<Vertex>());
                }
                map[*vi1]->push_back(*vi2);
                //if(map.find(id2) == map.end()){
                //    map[id2] = std::unique_ptr<std::vector<std::string>>(new std::vector<std::string>());
                //}
                //map[id2]->push_back(id1);
            }
        }
    }

    // identify full-sib and parent-offspring relationship and add edges to the graph
    // test if our map is correct
    // for(auto it = first_degs.begin(); it != first_degs.end(); it++){
    //     Vertex u = it->first;
    //     for(Vertex v : *(it->second)){
    //         std::string id1 = vertex_property_map[u];
    //         std::string id2 = vertex_property_map[v];
    //         std::cout << id1 << "\t" << id2 << std::endl;
    //     }
    // }

}