#include <math.h>
#include "function.h"

int getRelfromK(double K, int maxDeg){
    //double K = std::max((ibd1/4.0 + ibd2/2.0 - bkg/4.0)/tot_genome, 0.0); 
    if (K == 0){return -1;}
    int deg = (int)(-log2(K) + 0.5) - 1; // 四舍五入
    if (deg <= maxDeg){return deg;}
    else{return -1;}
}

void build_graph(Pedigree &pedigree, 
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs, 
    double tot_genome, double bkg_sharing, int maxDeg)
{
    auto vertex_property_map = boost::get(&sample::id, pedigree);
    boost::graph_traits<Pedigree>::vertex_iterator vi1, vi_end1;
    boost::graph_traits<Pedigree>::vertex_iterator vi2, vi_end2;
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
            std::cout << id1 << "\t" << id2 << "\t" << deg << std::endl;
        }
    }
}