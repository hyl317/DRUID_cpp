#include <math.h>
#include <memory>
#include <numeric>
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
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<std::string, std::string>, int> &results,
    double tot_genome, double bkg_sharing, int maxDeg)
{
    auto vertex_property_map = boost::get(&sample::id, pedigree);
    boost::graph_traits<Pedigree>::vertex_iterator vi1, vi_end1;
    boost::graph_traits<Pedigree>::vertex_iterator vi2, vi_end2;
    using Vertex = boost::graph_traits<Pedigree>::vertex_descriptor;
    std::vector<std::pair<Vertex, Vertex>> twins;
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
            }else if (deg == 0){
                twins.push_back(std::make_pair(*vi1, *vi2));
            }
        }
    }

    // remove one of the twins
    std::cout << "number of vertex: " << boost::num_vertices(pedigree) << std::endl;
    for(auto twin : twins){boost::remove_vertex(twin.second, pedigree);}
    std::cout << "number of vertex: " << boost::num_vertices(pedigree) << std::endl;

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

bool is_avunc(const std::string &fs1, const std::string &fs2, const std::string &avunc, 
        const std::map<std::pair<std::string, std::string>, Pair*> &allsegs, 
        const std::map<std::string, std::map<int, double>*> &snpmap)
{
    auto fs_p = allsegs.find(make_pair(fs1, fs2));
    assert(fs_p != allsegs.end());
    Pair &full_sib_pair = *(fs_p->second);
    auto fa_p1 = allsegs.find(make_pair(fs1, avunc));
    assert(fa_p1 != allsegs.end());
    Pair &sib1_avunc_pair = *(fa_p1->second);
    auto fa_p2 = allsegs.find(make_pair(fs2, avunc));
    assert(fa_p2 != allsegs.end());
    Pair &sib2_avunc_pair = *(fa_p2->second);

    double ibd011_tot = 0.0;
    double ibd0_tot = 0.0;
    // iterate over chromosomes for which sib1 share ibds with the putative avunc
    for(auto it = sib1_avunc_pair.ibd1_map->begin(); it != sib1_avunc_pair.ibd1_map->end(); it++){
        std::string chr_name = it->first;
        auto it2 = sib2_avunc_pair.ibd1_map->find(chr_name);
        // if sib2 don't share any ibd1 region with the avunc on this chromosome, then there won't be any ibd110 region on this chromosome
        if (it2 == sib2_avunc_pair.ibd1_map->end()){continue;}
        else{
            ibdSegments ibd11;
            interval_intersection(*(it->second), *(it2->second), ibd11);
            // find ibd0 region in the full-sib pair
            // to find ibd0 region, first take the union of ibd1 and ibd2 region
            // then take the complement of their union
            ibdSegments ibd1or2;
            if (full_sib_pair.ibd1_map->find(chr_name) != full_sib_pair.ibd1_map->end()
                && full_sib_pair.ibd2_map->find(chr_name) == full_sib_pair.ibd2_map->end()){
                ibd1or2 = *(full_sib_pair.ibd1_map->find(chr_name)->second);
            }else if (full_sib_pair.ibd1_map->find(chr_name) == full_sib_pair.ibd1_map->end()
                && full_sib_pair.ibd2_map->find(chr_name) != full_sib_pair.ibd2_map->end()){
                ibd1or2 = *(full_sib_pair.ibd2_map->find(chr_name)->second);
            }else if (full_sib_pair.ibd1_map->find(chr_name) != full_sib_pair.ibd1_map->end()
                && full_sib_pair.ibd2_map->find(chr_name) != full_sib_pair.ibd2_map->end()){
                ibdSegments &ibd1 = *(full_sib_pair.ibd1_map->find(chr_name)->second);
                ibdSegments &ibd2 = *(full_sib_pair.ibd2_map->find(chr_name)->second);
                interval_union(ibd1, ibd2, ibd1or2);
            }
            
            ibdSegments ibd0;
            interval_complement(ibd1or2, ibd0, 
                snpmap.find(chr_name)->second->begin()->second, 
                (--snpmap.find(chr_name)->second->end())->second);

            // take the intersection of ibd0 and ibd11 to get ibd011
            ibdSegments ibd011;
            interval_intersection(ibd11, ibd0, ibd011);
            std::vector<double> segLengths;
            std::for_each(ibd011.begin(), ibd011.end(), 
                [&](const std::pair<double, double> &p)
                {segLengths.push_back(p.second - p.first);});
            ibd011_tot += std::accumulate(segLengths.begin(), segLengths.end(), decltype(segLengths)::value_type(0));
        }
    }
    fprintf(stdout, "ibd011 region: %lf\n", ibd011_tot);
    return ibd011_tot >= AVUNC_011;
}