#include <math.h>
#include <memory>
#include <numeric>
#include <queue>
#include <chrono>
#include "function.h"
#include "constants.h"
#include <boost/graph/connected_components.hpp>


int getRelfromK(double K, int maxDeg){
    if (K == 0){return -1;}
    int deg = (int)(-log2(K) + 0.5) - 1; // 四舍五入
    if (deg <= maxDeg){return deg;}
    else{return -1;}
}

void build_graph(Pedigree &pedigree, PairIBD &allsegs, 
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results, std::map<Vertex, Vertex> &twins,
    const std::map<std::string, int> &id2index, double tot_genome, double bkg_sharing, int maxDeg)
{   
    auto vertex_property_map = boost::get(&sample::id, pedigree);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::set<std::pair<Vertex, Vertex>> pcs;
    std::map<Vertex, std::shared_ptr<std::unordered_set<Vertex>>> fs_degs;
    std::map<Vertex, std::shared_ptr<std::unordered_set<Vertex>>> second_degs;
    for(auto it = allsegs.begin(); it != allsegs.end(); it++){
        Vertex u, v;
        boost::tie(u, v) = it->first;
        Pair &p = *(it->second);
        double ibd1 = p.ibd1_tot;
        double ibd2 = p.ibd2_tot;
        double K = std::max((ibd1/4.0 + ibd2/2.0 - bkg_sharing/4.0)/tot_genome, 0.0);
        p.kin = K;
        int deg = getRelfromK(K, maxDeg);
        setDeg(u, v, deg, results);
        // store first and second degree pairs' sample names for later use
        if (deg == 1 || deg == 2){
            bool isFS = true;
            if (deg == 1){
                if (ibd2/tot_genome >= FULL_SIB_MIN_IBD2){
                    boost::add_edge(u, v, pedigree);
                    pedigree[boost::edge(u, v, pedigree).first].rel = FS;
                }else{
                    isFS = false;
                    pcs.insert(std::make_pair(u, v));
                }
            }else if(deg == 2){
                // check for possibility of DC, if so, no need to consider this pair for AV
                // therefore no need to add them to second_deg
                if (ibd2/tot_genome >= DC_MIN_IBD2){continue;}
            }

            auto &map = deg == 1? fs_degs : second_degs;
            if (isFS || deg == 2){
                if(map.find(u) == map.end()){
                    map[u] = std::unique_ptr<std::unordered_set<Vertex>>(new std::unordered_set<Vertex>());
                }
                map[u]->insert(v);
                if(map.find(v) == map.end()){
                    map[v] = std::unique_ptr<std::unordered_set<Vertex>>(new std::unordered_set<Vertex>());
                }
                map[v]->insert(u);
            }
        }else if (deg == 0){twins.insert(std::make_pair(u, v));}
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    fprintf(stdout, "settingn up allsegs_v takes %lf\n", d/1e6);

    t1 = std::chrono::high_resolution_clock::now();
    // add missing full-sib edges (and break incorrect ones)
    std::vector<int> components(boost::num_vertices(pedigree));
    int num_components = boost::connected_components(pedigree, &components[0]);
    std::map<int, std::shared_ptr<std::vector<Vertex>>> connected_comp_map;
    boost::graph_traits<Pedigree>::vertex_iterator vi, vi_end;
    for(boost::tie(vi, vi_end) = boost::vertices(pedigree); vi != vi_end; vi++){
        int comp_index = components[*vi];
        if (connected_comp_map.find(comp_index) == connected_comp_map.end()){
            connected_comp_map.insert(std::make_pair(comp_index, std::shared_ptr<std::vector<Vertex>>(new std::vector<Vertex>())));
        }
        connected_comp_map[comp_index]->push_back(*vi);
    }

    for(auto it = connected_comp_map.begin(); it != connected_comp_map.end(); it++){
        int size = it->second->size();
        if(size == 1){continue;}
        std::vector<Vertex> &sibs = *(it->second);
        std::map<Vertex, int> conn_map;
        for(Vertex u : sibs){
            int counter = 0;
            for(Vertex v : sibs){
                if (u == v){continue;}
                if (boost::edge(u, v, pedigree).second){counter++;}
            }
            conn_map.insert(std::make_pair(u, counter));
        }
            
        for(Vertex u : sibs){
            if (conn_map[u] >= (size-1)/2.0){
                for(Vertex v : sibs){
                    if (u == v){continue;}
                    else if (!boost::edge(u, v, pedigree).second){
                        boost::add_edge(u, v, pedigree);
                        pedigree[boost::edge(u, v, pedigree).first].rel = FS;
                        fs_degs[u]->insert(v);
                        fs_degs[v]->insert(u);
                        // since we have promoted (u,v) as FS, remove this pair from second_deg/pc pair if they were such inferred previously
                        if (second_degs.find(u) != second_degs.end()){second_degs[u]->erase(v);}
                        if (second_degs.find(v) != second_degs.end()){second_degs[v]->erase(u);}
                        if (pcs.find(std::make_pair(u, v)) != pcs.end()){pcs.erase(std::make_pair(u,v));}
                        if (pcs.find(std::make_pair(v, u)) != pcs.end()){pcs.erase(std::make_pair(v,u));}
                        results[make_pair_v(u, v)] = 1;
                        //std::cout << vertex_property_map[u] << " and " << vertex_property_map[v] << " is now a FS pair" << std::endl;
                    }
                }
            }else{
                for(Vertex v : sibs){
                    if (u == v){continue;}
                    else if(boost::edge(u, v, pedigree).second){
                        boost::remove_edge(u, v, pedigree);
                        fs_degs[u]->erase(v);
                        fs_degs[v]->erase(u);
                        // I assume they should be considered 2nd if not full-sib? Or could they be PO? And should we consider the possibility of AV?
                        results[make_pair_v(u, v)] = 2;
                        //std::cout << vertex_property_map[u] << " and " << vertex_property_map[v] << " is no longer a FS pair" << std::endl;
                    }
                }
            }
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    fprintf(stdout, "setting up FS edges takes %lf\n", d/1e6);

    t1 = std::chrono::high_resolution_clock::now();
    // add PC edges
    for(auto pc: pcs){
        boost::add_edge(pc.first, pc.second, pedigree);
        pedigree[boost::edge(pc.first, pc.second, pedigree).first].rel = PC;
        //std::cout << "add a PC edge between " << vertex_property_map[pc.first] << " and " << vertex_property_map[pc.second] << std::endl;
    }

    // polarize Parent-child relationship when we can
    for(auto pc : pcs){
        Vertex v1 = pc.first;
        Vertex v2 = pc.second;
        //std::cout << "polarizing v1 and v2: " << vertex_property_map[v1] << "\t" << vertex_property_map[v2] << std::endl;
        // if neither v1 nor v2 has any full-sibs, we can't do anything
        if (fs_degs.find(v1) == fs_degs.end() && fs_degs.find(v2) == fs_degs.end()){continue;}
        bool v2IsParent = false;
        bool v1IsParent = false;
        if (fs_degs.find(v1) != fs_degs.end()){
            // iterate over v1's full-sib pair to check if they form a parent-child pair with v2
            // if so, then v2 must be the parent
            int counter_formPCwithv2 = 0;
            for(auto fs : *fs_degs[v1]){
                bool isConnected = boost::edge(fs, v2, pedigree).second;
                //std::cout << "check v1's full-sib: " << vertex_property_map[fs] << " is connected? "<< isConnected << std::endl;
                if (isConnected){
                    Edge e = boost::edge(fs, v2, pedigree).first;
                    //std::cout << pedigree[e].rel << std::endl;
                    if (pedigree[e].rel == PC){counter_formPCwithv2++;}
                }
            }
            v2IsParent = counter_formPCwithv2 > fs_degs[v1]->size()*0.5;
            // now assign v2 as the parent of all full-sibs of v1
            if (v2IsParent){
                pedigree[boost::edge(v1, v2, pedigree).first].older = v2;
                pedigree[boost::edge(v1, v2, pedigree).first].polarized = true;
                for (auto fs : *fs_degs[v1]){
                    // don't want to add an edge twice becuase I used std::Vector to represent edge list
                    if (!boost::edge(fs, v2, pedigree).second){
                        boost::add_edge(fs, v2, pedigree);
                        //std::cout << "add a missed PC edge between " << vertex_property_map[fs] << " and " << vertex_property_map[v2] << std::endl;
                    }
                    pedigree[boost::edge(fs, v2, pedigree).first].rel = PC;
                    pedigree[boost::edge(fs, v2, pedigree).first].older = v2;
                    pedigree[boost::edge(fs, v2, pedigree).first].polarized = true;
                }
            }else if(!v2IsParent && counter_formPCwithv2 > 0){
                // get rid of false positive parent edge
                boost::remove_edge(v1, v2, pedigree);
                for (auto fs : *fs_degs[v1]){
                    if (boost::edge(fs, v2, pedigree).second){boost::remove_edge(fs, v2, pedigree);}
                }
            }
        }

        if (fs_degs.find(v2) != fs_degs.end()){
            // iterate over v2's full-sib pair to check if they form a parent-child pair with v1
            // if so, then v1 must be the parent
            int counter_formPCwithv1 = 0;
            for(auto fs : *fs_degs[v2]){
                bool isConnected = boost::edge(fs, v1, pedigree).second;
                //std::cout << "check v2's full-sib: " << vertex_property_map[fs] << " is connected? "<< isConnected << std::endl;
                if (isConnected){
                    Edge e = boost::edge(fs, v1, pedigree).first;
                    //std::cout << pedigree[e].rel << std::endl;
                    if (pedigree[e].rel == PC){counter_formPCwithv1++;}
                }
            }
            v1IsParent = counter_formPCwithv1 > 0.5*fs_degs[v2]->size();
            // now assign v1 as the parent of all full-sibs of v2
            if (v1IsParent){
                pedigree[boost::edge(v1, v2, pedigree).first].older = v1;
                pedigree[boost::edge(v1, v2, pedigree).first].polarized = true;
                for (auto fs : *fs_degs[v2]){
                    // don't want to add an edge twice becuase I used std::Vector to represent edge list
                    if (!boost::edge(fs, v1, pedigree).second){
                        boost::add_edge(fs, v1, pedigree);
                        //std::cout << "add a missed PC edge between " << vertex_property_map[fs] << " and " << vertex_property_map[v1] << std::endl;
                    }
                    pedigree[boost::edge(fs, v1, pedigree).first].rel = PC;
                    pedigree[boost::edge(fs, v1, pedigree).first].older = v1;
                    pedigree[boost::edge(fs, v1, pedigree).first].polarized = true;
                }
            }else if (!v1IsParent && counter_formPCwithv1 > 0){
                boost::remove_edge(v1, v2, pedigree);
                for(auto fs : *fs_degs[v2]){
                    if (boost::edge(fs, v1, pedigree).second){boost::remove_edge(fs, v1, pedigree);}
                }
            }
        }
        // if(v1IsParent && v2IsParent){
        //     std::cout << "something wrong with v1 and v2: " << vertex_property_map[v1] << "\t" << vertex_property_map[v2] << std::endl;
        //     return;
        // }
        assert(!(v1IsParent && v2IsParent));
    }
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    fprintf(stdout, "setting up Parent-offspring edges takes %lf\n", d/1e6);

    // construct second-degree edges
    t1 = std::chrono::high_resolution_clock::now();
    std::unordered_set<Vertex> checked_sibs;
    for(auto it = second_degs.begin(); it != second_degs.end(); it++){
        // here we assume that it->first is the younger generation and test if it forms AV relationship with its putative sescond relaatives
        // actually the younger generation could also be it->second
        // but this is not a problem here because
        // every sample in it->second also has a entry in the second_degs map where they are it->first
        Vertex focal_ind = it->first;
        if (checked_sibs.find(focal_ind) != checked_sibs.end()){continue;}
        if (fs_degs.find(focal_ind) == fs_degs.end()){
            // no full-siblings to this focal individual, no way to determine if the second is avuncular or not
            continue;
        }
        
        std::vector<Vertex> full_sib_vertex_set;
        full_sib_vertex_set.push_back(focal_ind);
        for(Vertex v : *fs_degs[focal_ind]){full_sib_vertex_set.push_back(v);}

        // take the union of second_deg of all the first-sibs
        // I did this because sometimes a person's AV might be classified as 3rd by kinship coefficient
        // this way we can potentially recover such AV pair if it is 2nd to that focal ind's full-sibs
        std::unordered_set<Vertex> second_deg_relatives;
        for(auto fs : full_sib_vertex_set){
            if (second_degs.find(fs) != second_degs.end()){
                for(auto sec : *second_degs[fs]){
                    second_deg_relatives.insert(sec);
                }
            }
        }

        for(Vertex avunc_candidate : second_deg_relatives){
            //std::string avunc_candidate_id = vertex_property_map[avunc_candidate];
            //std::cout << "checking " << avunc_candidate_id << std::endl;
            if (checkAvunc(full_sib_vertex_set, avunc_candidate, allsegs, id2index, snpmap)){
                //std::cout << avunc_candidate_id << " is a AV!" << std::endl;
                // add second-deg edges to the graph
                // what should we do if there is already an edge in between? Let's for now overwrite the previously edge
                for(Vertex v : full_sib_vertex_set){
                    if (!boost::edge(v, avunc_candidate, pedigree).second) {boost::add_edge(v, avunc_candidate, pedigree);}
                    //std::cout << vertex_property_map[v] << " and " << vertex_property_map[avunc_candidate] << " is now a AV pair" << std::endl;
                    pedigree[boost::edge(v, avunc_candidate, pedigree).first].rel = AV;
                    pedigree[boost::edge(v, avunc_candidate, pedigree).first].older = avunc_candidate;
                    pedigree[boost::edge(v, avunc_candidate, pedigree).first].polarized = true;
                    checked_sibs.insert(v);
                }
            }
        }

    }
    t2 = std::chrono::high_resolution_clock::now();
    d = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    fprintf(stdout, "setting up AV edges takes %lf\n", d/1e6);

    // remove all edges coming out of one of the twins
    // don't want to remove its vertex cuz otherwise it would mess up with vertex index of all others
    for(auto twin : twins){
        std::vector<Vertex> remove;
        auto out_edge_iter_pair = boost::out_edges(twin.first, pedigree);
        for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
            Vertex u = boost::source(*it, pedigree);
            Vertex v = boost::target(*it, pedigree);
            Vertex toRemove = u == twin.first? v : u;
            remove.push_back(toRemove);
        }
        for(Vertex v : remove){
            boost::remove_edge(twin.first, v, pedigree);
        }
    }

}

bool is_avunc(const Vertex fs1, const Vertex fs2, const Vertex avunc, const PairIBD &allsegs, 
    const std::map<std::string, int> &id2index, const std::map<std::string, std::map<int, double>*> &snpmap)
{
    auto fs_p = allsegs.find(make_pair_v(fs1, fs2));
    assert(fs_p != allsegs.end());
    Pair &full_sib_pair = *(fs_p->second);
    auto fa_p1 = allsegs.find(make_pair_v(fs1, avunc));
    assert(fa_p1 != allsegs.end());
    Pair &sib1_avunc_pair = *(fa_p1->second);
    auto fa_p2 = allsegs.find(make_pair_v(fs2, avunc));
    assert(fa_p2 != allsegs.end());
    Pair &sib2_avunc_pair = *(fa_p2->second);

    double ibd011_tot = 0.0;
    // iterate over chromosomes for which sib1 share ibds with the putative avunc
    for(auto it = id2index.begin(); it != id2index.end(); it++){
        std::string chr_name = it->first;
        int index = it->second;
        auto it1 = sib1_avunc_pair.ibd1_map[index];
        auto it2 = sib2_avunc_pair.ibd1_map[index];
        // if sib2 don't share any ibd1 region with the avunc on this chromosome, then there won't be any ibd110 region on this chromosome
        if (it1 == nullptr || it2 == nullptr){continue;}
        else{
            ibdSegments ibd11;
            interval_intersection(*it1, *it2, ibd11);
            // find ibd0 region in the full-sib pair
            // to find ibd0 region, first take the union of ibd1 and ibd2 region
            // then take the complement of their union
            ibdSegments ibd1or2;
            if (full_sib_pair.ibd1_map[index] != nullptr
                && (full_sib_pair.ibd2_map == nullptr || full_sib_pair.ibd2_map[index] == nullptr)){
                ibd1or2 = *(full_sib_pair.ibd1_map[index]);
            }else if (full_sib_pair.ibd1_map[index] == nullptr
                && full_sib_pair.ibd2_map != nullptr && full_sib_pair.ibd2_map[index] != nullptr){
                ibd1or2 = *(full_sib_pair.ibd2_map[index]);
            }else if (full_sib_pair.ibd1_map[index] != nullptr
                && full_sib_pair.ibd2_map != nullptr && full_sib_pair.ibd2_map[index] != nullptr){
                ibdSegments &ibd1 = *(full_sib_pair.ibd1_map[index]);
                ibdSegments &ibd2 = *(full_sib_pair.ibd2_map[index]);
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
    return ibd011_tot >= AVUNC_011;
}

bool checkAvunc(const std::vector<Vertex> &full_sibs, const Vertex avunc, 
    const PairIBD &allsegs, const std::map<std::string, int> &id2index, 
    const std::map<std::string, std::map<int, double>*> &snpmap)
{
    // check if avunc is the avunc of the given set of full_sibs
    int num_sibs = full_sibs.size();
    for(int i = 0; i < num_sibs; i++){
        for(int j = i+1; j < num_sibs; j++){
            if(is_avunc(full_sibs[i], full_sibs[j], avunc, allsegs, id2index, snpmap)){return true;}
        }
    }
    return false;
}

void run_druid(Pedigree &pedigree, const PairIBD &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results, const std::map<std::string, int> &id2index,
    FileOrGZ<FILE *> &logFile, double tot_genome, double bkg_sharing, int maxDeg)
{   
    logFile.printf("Start the primary DRUID algorithm...\n");
    std::vector<int> components(boost::num_vertices(pedigree));
    int num_components = boost::connected_components(pedigree, &components[0]);
    logFile.printf("\tnumber of connected components: %d\n", num_components);

    std::map<int, std::shared_ptr<std::vector<Vertex>>> comp_map;
    boost::graph_traits<Pedigree>::vertex_iterator vi, vi_end;
    for(boost::tie(vi, vi_end) = boost::vertices(pedigree); vi != vi_end; vi++){
        int comp_index = components[*vi];
        if (comp_map.find(comp_index) == comp_map.end()){
            comp_map.insert(std::make_pair(comp_index, std::shared_ptr<std::vector<Vertex>>(new std::vector<Vertex>())));
        }
        comp_map[comp_index]->push_back(*vi);
    }

    auto vertex_property_map = boost::get(&sample::id, pedigree);
    for(int i = 0; i < num_components; i++){
        auto ordered = std::shared_ptr<std::vector<Vertex>>(new std::vector<Vertex>());
        preorder(*(comp_map.find(i)->second), pedigree, *ordered);
        comp_map[i] = ordered;
        // test ordering
        //for(auto it = ordered->begin(); it != ordered->end(); it++){
        //    fprintf(stdout, "%s\n", vertex_property_map[*it].c_str());
        //}
        //fprintf(stdout, "\n\n");
    }
    //return;

    // test
    // for(int i = 0; i < num_components; i++){
    //    ConnInfo con;
    //    Vertex u = (*comp_map[i])[0];
    //    fprintf(stdout, "focal ind: %s\n", vertex_property_map[u].c_str());
    //    grabCloseRelatives_o(u, con, pedigree);
    //    if(!isSingleton(con)){printConnInfo(con, pedigree);}
    // }
    // return;
    //end of test

    for(int i = 0; i < num_components; i++){
        for(int j = i+1; j < num_components; j++){
            std::unordered_set<Vertex> visited1;
            for(Vertex u : *comp_map[i]){
                if (visited1.find(u) != visited1.end()){continue;}
                ConnInfo con1;
                grabCloseRelatives_o(u, con1, pedigree);
                bool isSingleton1 = isSingleton(con1);
                bool hasPC1 = false;
                Vertex w1;
                if (isSingleton1){std::tie(hasPC1, w1) = findUnpolarizedPC(u, pedigree);}
		        std::unordered_set<Vertex> visited2;
                //fprintf(stdout, "---------------------------s--------------------------------\n");
                //fprintf(stdout, "con1:\n");
                //printConnInfo(con1, pedigree);
                for(Vertex v : *comp_map[j]){
                    if(visited2.find(v) != visited2.end()){continue;}
                    ConnInfo con2;
                    grabCloseRelatives_o(v, con2, pedigree);
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
                        int av2 = -1;
                        visited1.insert(u);
                        if (!hasPC1){
                            av2 = oneVSpedigree(u, con2, visited2, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
                        }else{
                            auto pc = std::make_pair(u, w1);
                            visited1.insert(w1);
                            av2 = PCpairVSpedigree(pc, con2, visited2, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
                        }
                        propagateAlongPedigree(con1, con2, -1, av2, visited1, visited2, pedigree, results, maxDeg);
                    }else if(!isSingleton1 && isSingleton2){
                        int av1 = -1;
                        visited2.insert(v);
                        if (!hasPC2){
                            av1 = oneVSpedigree(v, con1, visited1, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
                        }else{
                            auto pc = std::make_pair(v, w2);
                            av1 = PCpairVSpedigree(pc, con1, visited1, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
                        }
                        propagateAlongPedigree(con1, con2, av1, -1, visited1, visited2, pedigree, results, maxDeg);
                    }else{
                        //printConnInfo(con1, pedigree);
                        //printConnInfo(con2, pedigree);
                        std::pair<int, int> aunts = pedigreeVSpedigree(con1, con2, visited1, visited2, allsegs, snpmap, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
                        propagateAlongPedigree(con1, con2, aunts.first, aunts.second, visited1, visited2, pedigree, results, maxDeg);
                    }

                }
                //fprintf(stdout, "--------------------------e---------------------------------\n");
            }
        }
    }

    return;
}


void postorder(const std::vector<Vertex> &components, const Pedigree &pedigree, std::vector<Vertex> &ordered)
{   

    // find the youngest generation
    std::queue<Vertex> youngest;
    for(Vertex u : components){
        bool isYounger = true;
        auto out_edge_iter_pair = boost::out_edges(u, pedigree);
        for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
            Edge e = *it;
            if(pedigree[e].polarized && pedigree[e].older == u){isYounger = false;}
        }
        if(isYounger){youngest.push(u);}
    }

    std::unordered_set<Vertex> checked;
    while(!youngest.empty()){
        Vertex u = youngest.front();
        youngest.pop();
        if (checked.find(u) != checked.end()){continue;}
        ordered.push_back(u);
        auto out_edge_iter_pair = boost::out_edges(u, pedigree);
        for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
            Edge e = *it;
            if(pedigree[e].polarized && pedigree[e].older != u){
                Vertex s = boost::source(e, pedigree);
                Vertex t = boost::target(e, pedigree);
                Vertex v = u == s? t : s;
                youngest.push(v);
            }
        }
        checked.insert(u);
    }

    assert(components.size() == ordered.size());
}

void preorder(const std::vector<Vertex> &components, const Pedigree &pedigree, std::vector<Vertex> &ordered)
{
    // find the oldest generation
    std::queue<Vertex> oldest;
    for(Vertex u : components){
        bool isOlder = true;
        auto out_edge_iter_pair = boost::out_edges(u, pedigree);
        for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
            Edge e = *it;
            if(pedigree[e].polarized && pedigree[e].older != u){isOlder = false;}
        }
        if(isOlder){oldest.push(u);}
    }

    std::unordered_set<Vertex> checked;
    while(!oldest.empty()){
        Vertex u = oldest.front();
        oldest.pop();
        if (checked.find(u) != checked.end()){continue;}
        ordered.push_back(u);
        auto out_edge_iter_pair = boost::out_edges(u, pedigree);
        for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
            Edge e = *it;
            if(pedigree[e].polarized && pedigree[e].older == u){
                Vertex s = boost::source(e, pedigree);
                Vertex t = boost::target(e, pedigree);
                Vertex v = u == s? t : s;
                oldest.push(v);
            }
        }
        checked.insert(u);
    }

    assert(components.size() == ordered.size());
}

bool isFS2Everyone(const Vertex &u, const std::vector<Vertex> &fs, const Pedigree &pedigree)
{
    for(Vertex v : fs){
        if(!boost::edge(u, v, pedigree).second){return false;}
        else{
            // edge exists, check if it's a FS edge
            Edge e = boost::edge(u, v, pedigree).first;
            if (pedigree[e].rel != FS){return false;}
        }
    }
    return true;
}

void grabCloseRelatives_o(const Vertex u, ConnInfo &con, const Pedigree &pedigree)
{
    auto vertex_property_map = boost::get(&sample::id, pedigree);
    std::vector<Vertex> children;
    std::vector<Vertex> fs;
    std::vector<Vertex> av;
    auto out_edge_iter_pair = boost::out_edges(u, pedigree);
    for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
        Edge e = *it;
        Vertex s = boost::source(e, pedigree);
        Vertex t = boost::target(e, pedigree);
        Vertex other = s == u ? t : s;
        if (pedigree[e].rel == PC && pedigree[e].polarized && pedigree[e].older == u){    
            children.push_back(other);
        }else if (pedigree[e].rel == FS){fs.push_back(other);}
        else if (pedigree[e].rel == AV && pedigree[e].polarized && pedigree[e].older == u){av.push_back(other);}
    }
    if (children.empty()){
        // check if u's full-sib has children 
        // (this is the same as checking u's AV, 
        // but I believe fulls ib and PC are less error-prone than AV detection, 
        // so I prefer to check fs's children first)
        for(Vertex sib : fs){
            //fprintf(stdout, "checking full-sib %s\n", vertex_property_map[sib].c_str());
            out_edge_iter_pair = boost::out_edges(sib, pedigree);
            for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
                Edge e = *it;
                Vertex s = boost::source(e, pedigree);
                Vertex t = boost::target(e, pedigree);
                Vertex other = s == sib ? t : s;
                //fprintf(stdout, "the other end is %s, edge type is %d\n", vertex_property_map[other].c_str(), pedigree[e].rel);
                if (pedigree[e].rel == PC && pedigree[e].polarized && pedigree[e].older == sib){
                    grabCloseRelatives(other, con, pedigree);
                    //assert(con.av1.contains(u) || con.av2.contains(u));
                    return;
                }
            }
        }

        if (!av.empty()){grabCloseRelatives(av[0], con, pedigree);}
        else {grabCloseRelatives(u, con, pedigree);}
    }else{
        for(Vertex c : children){
            out_edge_iter_pair = boost::out_edges(c, pedigree);
            for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
                Edge e = *it;
                if (pedigree[e].rel == PC && pedigree[e].polarized && pedigree[e].older == c){
                    Vertex s = boost::source(e, pedigree);
                    Vertex t = boost::target(e, pedigree);
                    Vertex other = s == c ? t : s;
                    grabCloseRelatives(other, con, pedigree);
                    //assert(con.gp1.contains(u) || con.gp2.contains(u));
                    return;
                }
            }
        }
        grabCloseRelatives(children[0], con, pedigree);
        //assert(con.p.contains(u));
    }
}

void grabCloseRelatives(const Vertex &u, ConnInfo &con, const Pedigree &pedigree)
{   
    con.fs.push_back(u);
    // first we will find u's full-sibslings, aunts/uncles, parents
    auto out_edge_iter_pair = boost::out_edges(u, pedigree);
    for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
        Edge e = *it;
        Vertex s = boost::source(e, pedigree);
        Vertex t = boost::target(e, pedigree);
        Vertex other = s == u ? t : s;
        if (pedigree[e].rel == PC && pedigree[e].polarized && pedigree[e].older == other){
            con.p.push_back(other);
        }else if (pedigree[e].rel == FS){
            con.fs.push_back(other);
        }else if (pedigree[e].rel == AV && pedigree[e].polarized && pedigree[e].older == other){
            if (con.av1.empty()){con.av1.push_back(other);}
            else{
                // if the other vertex form full-sib to everyone in av1, add it to av1
                if(isFS2Everyone(other, con.av1, pedigree)){con.av1.push_back(other);}
                else{
                    if (con.av2.empty()){con.av2.push_back(other);}
                    else if (isFS2Everyone(other, con.av2, pedigree)){con.av2.push_back(other);}
                }
            }
        }
    }

    // now check if we have grandparents
    std::unordered_set<Vertex> added;
    for(Vertex av : con.av1){
        auto out_edge_iter_pair = boost::out_edges(av, pedigree);
        for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
            Edge e = *it;
            Vertex s = boost::source(e, pedigree);
            Vertex t = boost::target(e, pedigree);
            Vertex other = s == av ? t : s;
            if (pedigree[e].rel == PC && pedigree[e].polarized && pedigree[e].older == other){
                if (added.find(other) == added.end()){
                    con.gp1.push_back(other);
                    added.insert(other);
                }
            }
        }
    }

    for(Vertex av : con.av2){
        auto out_edge_iter_pair = boost::out_edges(av, pedigree);
        for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
            Edge e = *it;
            Vertex s = boost::source(e, pedigree);
            Vertex t = boost::target(e, pedigree);
            Vertex other = s == av ? t : s;
            if (pedigree[e].rel == PC && pedigree[e].polarized && pedigree[e].older == other){
                if (added.find(other) == added.end()){
                    con.gp2.push_back(other);
                    added.insert(other);
                }
            }
        }
    }
}

std::pair<bool, std::size_t> findUnpolarizedPC(Vertex u, const Pedigree &pedigree)
{
    auto out_edge_iter_pair = boost::out_edges(u, pedigree);
    for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
        Edge e = *it;
        if (pedigree[e].rel == PC && !pedigree[e].polarized){
            Vertex s = boost::source(e, pedigree);
            Vertex t = boost::target(e, pedigree);
            Vertex other = s == u ? t : s;
            return std::make_pair(true, other);
        }
    }
    return std::make_pair(false, 0);
}

bool isSingleton(const ConnInfo &con)
{
    return (con.gp1.empty() && con.gp2.empty() 
    && con.av1.empty() && con.av2.empty() && con.p.empty() &&
    con.fs.size() == 1);
}

double UnionIbdOverTwoSets(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2, 
    const std::map<std::string, int> &id2index, const PairIBD &allsegs)
{
    ibdMapType currUnion;
    for(auto s1 : set1){
        for(auto s2 : set2){
            auto p = allsegs.find(make_pair_v(s1, s2));
            if (p == allsegs.end()){continue;}
            else{
                const Pair &pair = *(p->second);
                for(auto it = id2index.begin(); it != id2index.end(); it++){
                    int index = it->second;
                    if (pair.ibd1_map[index] == nullptr){continue;}
                    std::string chrName = it->first;
                    if (currUnion.find(chrName) == currUnion.end()){
                        currUnion.insert(std::make_pair(chrName, new ibdSegments()));
                    }
                    ibdSegments *dest = new ibdSegments();
                    interval_union(*(pair.ibd1_map[index]), *(currUnion.find(chrName)->second), *dest);
                    ibdSegments *prev_ptr = currUnion.find(chrName)->second;
                    currUnion[chrName] = dest;
                    delete prev_ptr;
                }

            }
        }
    }

    double tot_length = 0.0;
    for(auto it = currUnion.begin(); it != currUnion.end(); it++){
        ibdSegments *tmp = it->second;
        std::vector<double> segLengths;
        std::for_each(tmp->begin(), tmp->end(), 
                [&](const std::pair<double, double> &p)
                {segLengths.push_back(p.second - p.first);});
        tot_length += std::accumulate(segLengths.begin(), segLengths.end(), decltype(segLengths)::value_type(0));
    }

    // clean up
    for(auto it = currUnion.begin(); it != currUnion.end(); it++){
        delete it->second;
    }

    return tot_length;
}

double averageKinship(Vertex u, const std::vector<Vertex> &set, const PairIBD &allsegs)
{   
    assert(!set.empty());
    double average = 0.0;
    for(Vertex s : set){
        auto it = allsegs.find(make_pair_v(u, s));
        if (it != allsegs.end()){average += it->second->kin;}
    }
    return average/((double) set.size());
}

double averageKinshipBetweenTwoSets(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2, const PairIBD &allsegs)
{
    assert(!set1.empty());
    assert(!set2.empty());
    double average = 0.0;
    for(Vertex s1 : set1){
        for(Vertex s2 : set2){
            auto it = allsegs.find(make_pair_v(s1, s2));
            if (it != allsegs.end()){average += it->second->kin;}
        }
    }
    return average/((double) set1.size()*set2.size());
}

void inferFStoSingleDistantRelative(Vertex d, const std::vector<Vertex> &fs,
    std::unordered_set<Vertex> &visited, const PairIBD &allsegs, 
    std::map<std::pair<Vertex, Vertex>, int> &results, 
    const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg)
{
    std::vector<Vertex> set1;
    set1.push_back(d);
    int numSibs = fs.size();
    double Tg = getTg(0, numSibs);
    double k1 = (UnionIbdOverTwoSets(set1, fs, id2index, allsegs) - calc_bkg_sharing(0, numSibs, bkg_sharing))/(Tg*tot_genome);
    double K = std::max(0.0, k1/4.0);
    int deg = resetRelationship(getRelfromK(K, maxDeg), 1, maxDeg);
    for(Vertex v : fs){
        //results[make_pair_v(d, v)] = deg;
        setDeg(d, v, deg, results);
        visited.insert(v);
    }
}

int oneVSpedigree(Vertex u, const ConnInfo &con, std::unordered_set<Vertex> &visited,
    const PairIBD &allsegs, std::map<std::pair<Vertex, Vertex>, int> &results,
    const Pedigree &pedigree, const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg)
{   
    if(con.av1.empty() && con.av2.empty() && con.p.empty()){
        inferFStoSingleDistantRelative(u, con.fs, visited, allsegs, results, id2index, bkg_sharing, tot_genome, maxDeg);
        return -1;
    }else if(con.av1.empty() && con.av2.empty() && !con.p.empty()){
        // use full-sibs's parents if applicable
        double maxK = 0.0;
        for(Vertex v : con.fs){
            auto it = allsegs.find(make_pair_v(u, v));
            if(it == allsegs.end()){continue;}
            else if(it->second->kin > maxK){maxK = it->second->kin;}
        }
        if (con.p.size() == 1){
            auto it = allsegs.find(make_pair_v(con.p[0], u));
            if (it != allsegs.end() && it->second->kin > maxK){
                int deg = resetRelationship(getRelfromK(it->second->kin, maxDeg), 1, maxDeg);
                for(Vertex v : con.fs){
                    //results[make_pair_v(u,v)] = deg;
                    setDeg(u, v, deg, results);
                    visited.insert(v);
                }
                visited.insert(con.p[0]);
                //fprintf(stdout, "insert %d into visited, inside oneVSpedigree\n", con.p[0]);
            }else{
                inferFStoSingleDistantRelative(u, con.fs, visited, allsegs, results, id2index, bkg_sharing, tot_genome, maxDeg);
            }
        }else if(con.p.size() ==2){
            std::vector<Vertex> sib1 {u};
            int index = whichParent2Include(con.p, con.fs, sib1, allsegs);
            if (index != -1){
                int deg = resetRelationship(getRelfromK(allsegs.find(make_pair_v(con.p[index], u))->second->kin, maxDeg), 1, maxDeg);
                for(Vertex v : con.fs){
                    //results[make_pair_v(u, v)] = deg;
                    setDeg(u, v, deg, results);
                    visited.insert(v);
                }
                visited.insert(con.p[index]);
                //fprintf(stdout, "insert %d to visited, inside oneVSpedigree\n");
            }else{
                // no strong evidence to choose among the two parents, this implies that these two connected components are mostly unrelated
                //fprintf(stdout, "no parents satisfy the criterion, abandon this branch\n");
                std::for_each(con.fs.begin(), con.fs.end(), [&](Vertex fs){visited.insert(fs);});
                visited.insert(con.p[0]);
                visited.insert(con.p[1]);
            }
        }
        return -1;
    }else{
        bool includeAunts1 = includeAunts(con.av1, con.fs, u, allsegs); // if con.av1 is empty this still works
        bool includeAunts2 = includeAunts(con.av2, con.fs, u, allsegs);
        int index_av = -1;
        if (includeAunts1 && includeAunts2){
            // use the set of aunts with higher average kinship to d
            double average1 = averageKinship(u, con.av1, allsegs);
            double average2 = averageKinship(u, con.av2, allsegs);
            index_av = average1 >= average2 ? 1 : 2;
        }else if(!includeAunts1 && includeAunts2){index_av = 2;}
        else if(includeAunts1 && !includeAunts2){index_av = 1;}

        if (index_av != -1){
            // check if we can use grandparents
            const std::vector<Vertex> &av2use = index_av == 1 ? con.av1 : con.av2; 
            const std::vector<Vertex> &gp2check = index_av == 1 ? con.gp1 : con.gp2;
            std::vector<Vertex> sib1 {u};
            int index_gp = whichParent2Include(gp2check, av2use, sib1, allsegs);

            if (index_gp != -1){
                // can use grandparents for inference
                //fprintf(stdout, "use grandparents for inference\n");
                int deg_gp = getRelfromK(allsegs.find(make_pair_v(gp2check[index_gp], u))->second->kin, maxDeg);
                int deg_av = resetRelationship(deg_gp, 1, maxDeg);
                int deg_fs = resetRelationship(deg_gp, 2, maxDeg);
                for(Vertex a : av2use){
                    //results[make_pair_v(a, u)] = deg_av;
                    setDeg(a, u, deg_av, results);
                    visited.insert(a);
                }
                for(Vertex fs : con.fs){
                    //results[make_pair_v(fs, u)] = deg_fs;
                    setDeg(fs, u, deg_fs, results);
                    visited.insert(fs);
                }
                visited.insert(gp2check[index_gp]);
            }else{
                //fprintf(stdout, "no grandparents satisfy the criterion; use AV set only\n");
                std::vector<Vertex> set1;
                set1.push_back(u);
                std::vector<Vertex> set2;
                set2.insert(set2.end(), con.fs.begin(), con.fs.end());
                set2.insert(set2.end(), av2use.begin(), av2use.end());
                int numAV = av2use.size();
                int numSibs = con.fs.size();
                double Tg = getTg(numAV, numSibs);
                double k1 = (UnionIbdOverTwoSets(set1, set2, id2index, allsegs) - calc_bkg_sharing(numAV, numSibs, tot_genome))/(Tg*tot_genome);
                double K = std::max(0.0, k1/4.0);
                int deg_gp = getRelfromK(K, maxDeg);
                int deg_av = resetRelationship(deg_gp, 1, maxDeg);
                for(Vertex a : av2use){
                    //results[make_pair_v(a, u)] = deg_av;
                    setDeg(a, u, deg_av, results);
                    visited.insert(a);
                }
                int deg_fs = resetRelationship(deg_gp, 2, maxDeg);
                for(Vertex fs : con.fs){
                    //results[make_pair_v(fs, u)] = deg_fs;
                    setDeg(fs, u, deg_fs, results);
                    visited.insert(fs);
                }
            }

            // if we have parents of this sibling set, set the parent that belong to this parent to be of the same degree as AV set
            for(Vertex p : con.p){
                if(isFS2Everyone(p, av2use, pedigree)){
                    //results[make_pair_v(p, u)] = results[make_pair_v(av2use[0], u)];
                    setDeg(p, u, getDeg(av2use[0], u, results), results);
                    visited.insert(p);
                }
            }

        }else{
            // if no aunts/uncle sets satisfy the criterion, we can only use full-sibs
            //fprintf(stdout, "no AV satisfy the criterion, use fs only\n");
            inferFStoSingleDistantRelative(u, con.fs, visited, allsegs, results, id2index, bkg_sharing, tot_genome, maxDeg);
        }
        return index_av;
    }
}


bool includeAunts(const std::vector<Vertex> &aunts, const std::vector<Vertex> &sibs, Vertex d, const PairIBD &allsegs)
{
    // include the set of aunts if max k{a,d} > min k{s,d}
    // this is used for oneVSpedigree
    // Vertex d is the distant relative of interest
    // this return FALSE if the aunts vector is empty
    double min_ksd = 1.0;
    for(Vertex v : sibs){
        auto it = allsegs.find(make_pair_v(v, d));
        if (it == allsegs.end()){
            min_ksd = 0.0;
            break;
        }else{
            double k = it->second->kin;
            if (k < min_ksd){min_ksd = k;}
        }
    }

    double max_kad = 0.0;
    for(Vertex v : aunts){
        auto it = allsegs.find(make_pair_v(v, d));
        if (it == allsegs.end()){continue;}
        else{
            double k = it->second->kin;
            if(k > max_kad){max_kad = k;}
        }
    }

    return max_kad > min_ksd;
}

double minKinshipBetweenTwoSibset(const std::vector<Vertex> &sib1, const std::vector<Vertex> &sib2, const PairIBD &allsegs)
{
    // find min k_{s1, s2}
    double min_ks1s2 = 1.0;
    for(Vertex u : sib1){
        for(Vertex v : sib2){
            auto it = allsegs.find(make_pair_v(u, v));
            if (it == allsegs.end()){
                return 0.0;
            }else if(it->second->kin < min_ks1s2){min_ks1s2 = it->second->kin;}
        }
    }
    return min_ks1s2;
}

double maxKinshipBetweenTwoSibset(const std::vector<Vertex> &sib1,
    const std::vector<Vertex> &sib2, const PairIBD &allsegs)
{
    // find max k_{s1, s2}
    double max_ks1s2 = 0.0;
    for(Vertex u : sib1){
        for(Vertex v : sib2){
            auto it = allsegs.find(make_pair_v(u, v));
            if (it != allsegs.end() && it->second->kin > max_ks1s2){max_ks1s2 = it->second->kin;}
        }
    }
    return max_ks1s2;
}

bool includeAunts(const std::vector<Vertex> &aunts, const std::vector<Vertex> &sibs,
    double min_ks1s2, const PairIBD &allsegs)
{
    double max_kas2 = 0.0;
    for(Vertex a : aunts){
        for(Vertex s : sibs){
            auto it = allsegs.find(make_pair_v(a, s));
            if (it == allsegs.end()){continue;}
            else if (it->second->kin > max_kas2){max_kas2 = it->second->kin;}
        }
    }
    return max_kas2 > min_ks1s2;
}

int whichAV2Include(const std::vector<Vertex> &av11, const std::vector<Vertex> &av12, 
    const std::vector<Vertex> &sib1, const std::vector<Vertex> &sib2, const PairIBD &allsegs)
{
    if (av11.empty() && av12.empty()){return -1;}
    else{
        double min_ks1s2 = minKinshipBetweenTwoSibset(sib1, sib2, allsegs);
        if(av11.empty() && !av12.empty()){
            return includeAunts(av12, sib2, min_ks1s2, allsegs) ? 2 : -1;
        }else if (!av11.empty() && av12.empty()){
            return includeAunts(av11, sib2, min_ks1s2, allsegs) ? 1 : -1;
        }else{
            bool isGreater1 = includeAunts(av11, sib2, min_ks1s2, allsegs);
            bool isGreater2 = includeAunts(av12, sib2, min_ks1s2, allsegs);
            if (isGreater1 && !isGreater2){return 1;}
            else if (!isGreater1 && isGreater2){return 2;}
            else if (!isGreater1 && !isGreater2){return -1;}
            else{
                double average1 = averageKinshipBetweenTwoSets(av11, sib2, allsegs);
                double average2 = averageKinshipBetweenTwoSets(av12, sib2, allsegs);
                return average1 >= average2 ? 1 : 2;
            }
        }
    }
}

bool includeParent(Vertex p, double max_ks1s2, const std::vector<Vertex> &sibs, const PairIBD &allsegs)
{  
    // return true if max k_{p, s2} > max_ks1s2
    for(Vertex s2 : sibs){
        auto it = allsegs.find(make_pair_v(p, s2));
        if (it != allsegs.end() && it->second->kin > max_ks1s2){return true;}
    }
    return false;
}


int whichParent2Include(const std::vector<Vertex> &parents, 
    const std::vector<Vertex> &sib1, const std::vector<Vertex> &sib2, const PairIBD &allsegs)
{
    // parents is the parents of sib1
    if (parents.size() == 0){return -1;}
    else{
        double max_ks1s2 = maxKinshipBetweenTwoSibset(sib1, sib2, allsegs);
        // how to determine if a parent should be included?
        // use max k_{p, s2} > max k_{s1,s2}
        if (parents.size() == 1){
            return includeParent(parents[0], max_ks1s2, sib2, allsegs) ? 0 : -1;
        }else{
            bool include1 = includeParent(parents[0], max_ks1s2, sib2, allsegs);
            bool include2 = includeParent(parents[1], max_ks1s2, sib2, allsegs);
            if (!include1 && !include2){return -1;}
            else if (!include1 && include2){return 1;}
            else if (include1 && !include2){return 0;}
            else{
                // choose the parent with higher average kinship coefficient with sib2
                double average1 = averageKinship(parents[0], sib2, allsegs);
                double average2 = averageKinship(parents[1], sib2, allsegs);
                return average1 > average2 ? 0 : 1;
            }
        }
    }
}

void updateSibsetByTheirParent(int index, 
    const std::vector<Vertex> &fs, const std::vector<Vertex> &parents, 
    const ConnInfo &con2,
    std::unordered_set<Vertex> &visited1, std::unordered_set<Vertex> &visited2,
    const PairIBD &allsegs, std::map<std::pair<Vertex, Vertex>, int> &results,
    const Pedigree &pedigree, const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg)
{
    //fprintf(stdout, "updateSibsetByTheirParent\n");
    Vertex p2use = parents[index];
    visited1.insert(p2use); // samples in fs have already been added to visited1 before this function is called from pedigreeVSpedigree
    //fprintf(stdout, "insert %d\n", p2use);
    oneVSpedigree(p2use, con2, visited2, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
    // update fs's relationship to con2
    if (!con2.gp1.empty()){
        for(Vertex gp : con2.gp1){
            //int base_deg = results[make_pair_v(p2use, gp)];
            int base_deg = getDeg(p2use, gp, results);
            int reset_deg = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(gp, fs, results, reset_deg);
        }
    }

    if (!con2.gp2.empty()){
        for(Vertex gp : con2.gp2){
            //int base_deg = results[make_pair_v(p2use, gp)];
            int base_deg = getDeg(p2use, gp, results);
            int reset_deg = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(gp, fs, results, reset_deg);
        }
    }

    if (!con2.av1.empty()){
        for(Vertex av : con2.av1){
            //int base_deg = results[make_pair_v(p2use, av)];
            int base_deg = getDeg(p2use, av, results);
            int reset_deg = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(av, fs, results, reset_deg);
        }
    }

    if (!con2.av2.empty()){
        for(Vertex av : con2.av2){
            //int base_deg = results[make_pair_v(p2use, av)];
            int base_deg = getDeg(p2use, av, results);
            int reset_deg = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(av, fs, results, reset_deg);
        }
    }

    if (!con2.p.empty()){
        for(Vertex p : con2.p){
            //int base_deg = results[make_pair_v(p2use, p)];
            int base_deg = getDeg(p2use, p, results);
            int reset_deg = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(p, fs, results, reset_deg);
        }
    }

    for(Vertex sib2 : con2.fs){
        //int base_deg = results[make_pair_v(p2use, sib2)];
        int base_deg = getDeg(p2use, sib2, results);
        int reset_deg = resetRelationship(base_deg, 1, maxDeg);
        setRelationshipBetweenOneSampleAndSet(sib2, fs, results, reset_deg);
    }
}


std::pair<int, int> pedigreeVSpedigree(const ConnInfo &con1, const ConnInfo &con2,
    std::unordered_set<Vertex> &visited1, std::unordered_set<Vertex> &visited2,
    const PairIBD &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results,
    const Pedigree &pedigree, const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg)
{   
    // How to properly handle parents and grandparents in this setting?
    int index_av1 = whichAV2Include(con1.av1, con1.av2, con1.fs, con2.fs, allsegs);
    int index_av2 = whichAV2Include(con2.av1, con2.av2, con2.fs, con1.fs, allsegs);
    std::vector<Vertex> set1;
    std::vector<Vertex> set2;
    set1.insert(set1.end(), con1.fs.begin(), con1.fs.end());
    set2.insert(set2.end(), con2.fs.begin(), con2.fs.end());
    std::for_each(con1.fs.begin(), con1.fs.end(), [&](const Vertex fs){visited1.insert(fs);});
    std::for_each(con2.fs.begin(), con2.fs.end(), [&](const Vertex fs){visited2.insert(fs);});
    int numSib1 = con1.fs.size();
    int numSib2 = con2.fs.size();
    int numAV1 = 0;
    int numAV2 = 0;
    // use aunts/uncle if we can
    // if no aunts/uncle can be used, check if we have parents

    Vertex p1, p2;
    bool p1_chosen = false;
    bool p2_chosen = false;
    if(index_av1 != -1){
        const std::vector<Vertex> &av2use = index_av1 == 1 ? con1.av1 : con1.av2;
        set1.insert(set1.end(), av2use.begin(), av2use.end());
        numAV1 = av2use.size();

        // check which parent, if any, is FS to the selected AV set
        for(Vertex p : con1.p){
            if (isFS2Everyone(p, av2use, pedigree)){
                p1 = p;
                p1_chosen = true;
            }
        }

        std::for_each(av2use.begin(), av2use.end(), [&](const Vertex a){visited1.insert(a);});
        // check if we can use grandparents
        const std::vector<Vertex> &gp2use = index_av1 == 1 ? con1.gp1 : con1.gp2;
        int index = whichParent2Include(gp2use, av2use, con2.fs, allsegs);
        if (index != -1){
            //fprintf(stdout, "grandparent to use: %d\n", index);
            updateSibsetByTheirGrandParent(index, index_av1, con1, con2, visited1, visited2, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
            return std::make_pair(index_av1, index_av2);
        }
    }else{
        int index = whichParent2Include(con1.p, con1.fs, con2.fs, allsegs);
        if (index != -1){
            //fprintf(stdout, "use parents from con1\n");
            updateSibsetByTheirParent(index, con1.fs, con1.p, con2, visited1, visited2, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
            return std::make_pair(index_av1, index_av2);
        }else if (index == -1 && con1.p.size() == 2){
            // these two connected components are most likely unrelated
            // we don't add con2.p to visited2 here because
            // although con1 is not related to con2.fs, it might be related to
            // one of con2.p who is married into the pedigree instead of being a descendant
            visited1.insert(con1.p[0]);
            visited1.insert(con1.p[1]);
            //fprintf(stdout, "neither parents of con1 fits. Abandon this branch\n");
        }
    }

    if(index_av2 != -1){
        const std::vector<Vertex> &av2use = index_av2 == 1 ? con2.av1 : con2.av2;
        set2.insert(set2.end(), av2use.begin(), av2use.end());
        numAV2 = av2use.size();

        for(Vertex p : con2.p){
            if(isFS2Everyone(p, av2use, pedigree)){
                p2 = p;
                p2_chosen = true;
            }
        }

        std::for_each(av2use.begin(), av2use.end(), [&](const Vertex a){visited2.insert(a);});
        // check if we can use grandparents
        const std::vector<Vertex> &gp2use = index_av2 == 1 ? con2.gp1 : con2.gp2;
        int index = whichParent2Include(gp2use, av2use, con1.fs, allsegs);
        if (index != -1){
            //fprintf(stdout, "grandparent to use: %d\n", index);
            updateSibsetByTheirGrandParent(index, index_av2, con2, con1, visited2, visited1, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
            return std::make_pair(index_av1, index_av2);
        }
    }else{
        int index = whichParent2Include(con2.p, con2.fs, con1.fs, allsegs);
        if (index != -1){
            //fprintf(stdout, "use parents from con2\n");
            updateSibsetByTheirParent(index, con2.fs, con2.p, con1, visited2, visited1, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
            return std::make_pair(index_av1, index_av2);
        }else if (index == -1 && con2.p.size() == 2){
            visited2.insert(con2.p[0]);
            visited2.insert(con2.p[1]);
            //fprintf(stdout, "neither parents of con2 fits. Abandon this branch\n");
        }
    }

    double Tg1 = getTg(numAV1, numSib1);
    double Tg2 = getTg(numAV2, numSib2);
    double bkg = calc_bkg_sharing(numAV1, numSib1, numAV2, numSib2, bkg_sharing);
    double unioned = UnionIbdOverTwoSets(set1, set2, id2index, allsegs);
    double k1 = (unioned - bkg)/(tot_genome*Tg1*Tg2);
    int deg = getRelfromK(std::max(0.0, k1/4.0), maxDeg);
    if (deg >= 0 && deg <= 3){
        // need to account for IBD2 in the grandparent/parent generation
        double t1, t2, k2;
        std::vector<Vertex> left;
        std::vector<Vertex> right;
        if (index_av1 != -1){
            left = index_av1 == 1 ? con1.av1 : con1.av2;
            //fprintf(stdout, "select av set %d for IBD0011\n", index_av1);
        }else{left = con1.fs;}
        if (index_av2 != -1){
            right = index_av2 == 1 ? con2.av1 : con2.av2;
            //fprintf(stdout, "select av set %d for IBD0011\n", index_av2);
        }else{right = con2.fs;}
        t1 = 1.0 - pow(0.5, left.size());
        t2 = 1.0 - pow(0.5, right.size());
        k2 = IBD0011(left, right, snpmap, allsegs, id2index)/(tot_genome*t1*t2);
        k1 = unioned/(tot_genome*Tg1*Tg2) - bkg/tot_genome;
        //printConnInfo(con1, pedigree);
        //printConnInfo(con2, pedigree);
        //fprintf(stdout, "k1: %lf, k2: %lf, k=%lf\n", k1, k2, k1/4.0 + k2/2.0);
        k1 -= k2;
        deg = getRelfromK(std::max(0.0, k1/4.0 + k2/2.0), maxDeg);
    }

    if (p1_chosen){visited1.insert(p1);}
    if (p2_chosen){visited2.insert(p2);}
    // update result map
    if (index_av1 != -1 && index_av2 != -1){
        // the inferred deg is between gp and gp for both pedigrees
        int deg_av2av = resetRelationship(deg, 2, maxDeg);
        const std::vector<Vertex> &av2use1 = index_av1 == 1 ? con1.av1 : con1.av2;
        const std::vector<Vertex> &av2use2 = index_av2 == 1 ? con2.av1 : con2.av2;
        setRelationshipBetweenTwoSets(av2use1, av2use2, results, deg_av2av);
        int deg_av2fs = resetRelationship(deg, 3, maxDeg);
        setRelationshipBetweenTwoSets(av2use1, con2.fs, results, deg_av2fs);
        setRelationshipBetweenTwoSets(av2use2, con1.fs, results, deg_av2fs);
        int deg_fs2fs = resetRelationship(deg, 4, maxDeg);
        setRelationshipBetweenTwoSets(con1.fs, con2.fs, results, deg_fs2fs);
        if (p1_chosen){
            setRelationshipBetweenOneSampleAndSet(p1, con2.fs, results, deg_av2fs);
            setRelationshipBetweenOneSampleAndSet(p1, av2use2, results, deg_av2av);
        }
        if (p2_chosen){
            setRelationshipBetweenOneSampleAndSet(p2, con1.fs, results, deg_av2fs);
            setRelationshipBetweenOneSampleAndSet(p2, av2use1, results, deg_av2av);
        }
        if (p1_chosen && p2_chosen){
            //results[make_pair_v(p1, p2)] = deg_av2av;
            setDeg(p1, p2, deg_av2av, results);
        }

    }else if (index_av1 != -1 && index_av2 == -1){
        // the inferred deg is between pedigree 1's gp to pedigree 2's parents
        int deg_av2fs = resetRelationship(deg, 2, maxDeg);
        const std::vector<Vertex> &av2use1 = index_av1 == 1 ? con1.av1 : con1.av2;
        setRelationshipBetweenTwoSets(av2use1, con2.fs, results, deg_av2fs);
        int deg_fs2fs = resetRelationship(deg, 3, maxDeg);
        setRelationshipBetweenTwoSets(con1.fs, con2.fs, results, deg_fs2fs);
        // no need to check p2 here because, con2 doesn't have AV set and
        // in that case, we checked if con2 has suitable parent to use, and if there is,
        // we would have called resetSibsetByTheirParents()
        if (p1_chosen){setRelationshipBetweenOneSampleAndSet(p1, con2.fs, results, deg_av2fs);}
    }else if (index_av1 == -1 && index_av2 != -1){
         // the inferred deg is between pedigree 1's gp to pedigree 2's parents
        int deg_av2fs = resetRelationship(deg, 2, maxDeg);
        const std::vector<Vertex> &av2use2 = index_av2 == 1 ? con2.av1 : con2.av2;
        setRelationshipBetweenTwoSets(av2use2, con1.fs, results, deg_av2fs);
        int deg_fs2fs = resetRelationship(deg, 3, maxDeg);
        setRelationshipBetweenTwoSets(con1.fs, con2.fs, results, deg_fs2fs);
        if (p2_chosen){setRelationshipBetweenOneSampleAndSet(p2, con1.fs, results, deg_av2fs);}
    }else{
        // inferred deg is between pedigree 1's parents to pedigree 2's parents
        int deg_fs2fs = resetRelationship(deg, 2, maxDeg);
        setRelationshipBetweenTwoSets(con1.fs, con2.fs, results, deg_fs2fs);
    }
    return std::make_pair(index_av1, index_av2);
}

void updateSibsetByTheirGrandParent(int index_gp, int index_av,
    const ConnInfo &con1, const ConnInfo &con2,
    std::unordered_set<Vertex> &visited1, std::unordered_set<Vertex> &visited2,
    const PairIBD &allsegs, std::map<std::pair<Vertex, Vertex>, int> &results,
    const Pedigree &pedigree, const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg)
{
    const std::vector<Vertex> &avset2use = index_av == 1 ? con1.av1 : con1.av2;
    const std::vector<Vertex> &gpset2use = index_av == 1 ? con1.gp1 : con1.gp2;
    Vertex gp2use = gpset2use[index_gp];
    visited1.insert(gp2use);
    oneVSpedigree(gp2use, con2, visited2, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);

    // determine if con1 has any parent that is in this pedigree
    Vertex p;
    bool useP = false;
    for(Vertex _p : con1.p){
        if (isFS2Everyone(_p, avset2use, pedigree)){
            p = _p;
            useP = true;
        }
    }
    if (useP){visited1.insert(p);}
    if (!con2.gp1.empty()){
        for(Vertex gp : con2.gp1){
            //int base_deg = results[make_pair_v(gp2use, gp)];
            int base_deg = getDeg(gp2use, gp, results);
            int deg2av = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(gp, avset2use, results, deg2av);
            if (useP){
                //results[make_pair_v(p, gp)] = deg2av;
                setDeg(p, gp, deg2av, results);
            }
            setRelationshipBetweenOneSampleAndSet(gp, con1.fs, results, resetRelationship(base_deg, 2, maxDeg));
        }
    }

    if (!con2.gp2.empty()){
        for(Vertex gp : con2.gp2){
            //int base_deg = results[make_pair_v(gp2use, gp)];
            int base_deg = getDeg(gp2use, gp, results);
            int deg2av = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(gp, avset2use, results, deg2av);
            if (useP){
                //results[make_pair_v(p, gp)] = deg2av;
                setDeg(p, gp, deg2av, results);
            }
            setRelationshipBetweenOneSampleAndSet(gp, con1.fs, results, resetRelationship(base_deg, 2, maxDeg));
        }
    }

    if (!con2.av1.empty()){
        for(Vertex av : con2.av1){
            //int base_deg = results[make_pair_v(gp2use, av)];
            int base_deg = getDeg(gp2use, av, results);
            int deg2av = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(av, avset2use, results, deg2av);
            if (useP){
                //results[make_pair_v(p, av)] = deg2av;
                setDeg(p, av, deg2av, results);
            }
            setRelationshipBetweenOneSampleAndSet(av, con1.fs, results, resetRelationship(base_deg, 2, maxDeg));
        }
    }

    if (!con2.av2.empty()){
        for(Vertex av : con2.av2){
            //int base_deg = results[make_pair_v(gp2use, av)];
            int base_deg = getDeg(gp2use, av, results);
            int deg2av = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(av, avset2use, results, deg2av);
            if (useP){
                //results[make_pair_v(p, av)] = deg2av;
                setDeg(p, av, deg2av, results);
            }
            setRelationshipBetweenOneSampleAndSet(av, con1.fs, results, resetRelationship(base_deg, 2, maxDeg));
        }
    }

    if (!con2.p.empty()){
        for(Vertex _p : con2.p){
            //int base_deg = results[make_pair_v(gp2use, _p)];
            int base_deg = getDeg(gp2use, _p, results);
            int deg2av = resetRelationship(base_deg, 1, maxDeg);
            setRelationshipBetweenOneSampleAndSet(_p, avset2use, results, deg2av);
            if (useP){
                //results[make_pair_v(p, _p)] = deg2av;
                setDeg(p, _p, deg2av, results);
            }
            setRelationshipBetweenOneSampleAndSet(_p, con1.fs, results, resetRelationship(base_deg, 2, maxDeg));
        }
    }

    for(Vertex sib2 : con2.fs){
        //int base_deg = results[make_pair_v(gp2use, sib2)];
        int base_deg = getDeg(gp2use, sib2, results);
        int deg2av = resetRelationship(base_deg, 1, maxDeg);
        setRelationshipBetweenOneSampleAndSet(sib2, avset2use, results, deg2av);
        if (useP){
            //results[make_pair_v(p, sib2)] = deg2av;
            setDeg(p, sib2, deg2av, results);
        }
        setRelationshipBetweenOneSampleAndSet(sib2, con1.fs, results, resetRelationship(base_deg, 2, maxDeg));
    }

}

void IBD0011_uniDirection(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2,
    ibdMapType &dest, const std::map<std::string, std::map<int, double>*> &snpmap, 
    const PairIBD &allsegs, const std::map<std::string, int> &id2index)
{
    int numSib1 = set1.size();
    int numSib2 = set2.size();
    for(int i = 0; i < numSib1; i++){
        Vertex u1 = set1[i];
        for(int j = i+1; j < numSib1; j++){
            ibdMapType ibd0Map;
            Vertex u2 = set1[j];
            std::pair<Vertex, Vertex> p = make_pair_v(u1, u2);
            ibdSegments **ibd1 = allsegs.find(p)->second->ibd1_map; // since p is a pair of full-sibs, they for sure will have segments, so no need to check allsegs.find(p) != allsegs.end()
            ibdSegments **ibd2 = allsegs.find(p)->second->ibd2_map;
            // first, find IBD0 region between u1 and u2
            // do this chromosome by chromosome
            for(auto it = id2index.begin(); it != id2index.end(); it++){
                int index = it->second;
                std::string chrName = it->first;
                auto ibd0 = new ibdSegments();
                bool isEmpty_ibd1 = ibd1[index] == nullptr;
                bool isEmpty_ibd2 = ibd2[index] == nullptr;
                double chr_s = snpmap.find(chrName)->second->begin()->second; 
                double chr_e = (--snpmap.find(chrName)->second->end())->second;
                if (!isEmpty_ibd1 && !isEmpty_ibd2){
                    ibdSegments ibdUnion;
                    interval_union(*(ibd1[index]), *(ibd2[index]), ibdUnion);
                    interval_complement(ibdUnion, *ibd0, chr_s, chr_e);
                }else if (!isEmpty_ibd1 && isEmpty_ibd2){
                    interval_complement(*(ibd1[index]), *ibd0, chr_s, chr_e);
                }else if (isEmpty_ibd1 && !isEmpty_ibd2){
                    interval_complement(*(ibd2[index]), *ibd0, chr_s, chr_e);
                }else{
                    ibd0->push_back(std::make_pair(chr_s, chr_e));
                }
                ibd0Map.insert(std::make_pair(chrName, ibd0));
            }

            for(int m = 0; m < numSib2; m++){
                Vertex v1 = set2[m];
                for(int n = 0; n < numSib2 && n != m; n++){
                    Vertex v2 = set2[n];
                    // now identify ibd0011 region, also chrom by chrom
                    auto it1 = allsegs.find(make_pair_v(u1, v1));
                    auto it2 = allsegs.find(make_pair_v(u2, v2));
                    if (it1 == allsegs.end() || it2 == allsegs.end()){continue;}
                    ibdSegments **pair1 = it1->second->ibd1_map;
                    ibdSegments **pair2 = it2->second->ibd1_map;
                    for(auto it = id2index.begin(); it != id2index.end(); it++){
                        int index = it->second;
                        std::string chrName = it->first;
                        if (pair1[index] == nullptr || pair2[index] == nullptr || ibd0Map.find(chrName) == ibd0Map.end()){continue;}
                        else{
                            ibdSegments ibd11;
                            interval_intersection(*(pair1[index]), *(pair2[index]), ibd11);
                            if(dest.find(chrName) == dest.end()){
                                dest.insert(std::make_pair(chrName, new ibdSegments));
                            }
                            ibdSegments ibd0011;
                            interval_intersection(ibd11, *(ibd0Map.find(chrName)->second), ibd0011);
                            // add the new ibd0011 region to the big union
                            ibdSegments *updated_0011Union = new ibdSegments();
                            interval_union(ibd0011, *(dest.find(chrName)->second), *updated_0011Union);
                            ibdSegments *prev_ptr = dest.find(chrName)->second;
                            dest[chrName] = updated_0011Union;
                            delete prev_ptr; 
                        }
                    }
                }
            }
            for(auto it = ibd0Map.begin(); it != ibd0Map.end(); it++){delete it->second;}
        }
    }
}

double IBD0011(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2,
    const std::map<std::string, std::map<int, double>*> &snpmap, 
    const PairIBD &allsegs, const std::map<std::string, int> &id2index)
{
    ibdMapType union1; // stores IBD0011 region for each of the chromosome
    IBD0011_uniDirection(set1, set2, union1, snpmap, allsegs, id2index);
    ibdMapType union2;
    IBD0011_uniDirection(set2, set1, union2, snpmap, allsegs, id2index);

    double ibd0011_tot = 0.0;
    for(auto it = snpmap.begin(); it != snpmap.end(); it++){
        std::string chr = it->first;
        auto iter1 = union1.find(chr);
        auto iter2 = union2.find(chr);
        std::vector<double> segLengths;
        if (iter1 != union1.end() && iter2 != union2.end()){
            ibdSegments dest;
            interval_union(*(iter1->second), *(iter2->second), dest);
            std::for_each(dest.begin(), dest.end(), 
                [&](const std::pair<double, double> &p)
                {segLengths.push_back(p.second - p.first);});
        }else if (iter1 != union1.end() || iter2 != union2.end()){
            ibdSegments &nonEmpty = iter1 != union1.end() ? *(iter1->second) : *(iter2->second);
            std::for_each(nonEmpty.begin(), nonEmpty.end(), 
                [&](const std::pair<double, double> &p)
                {segLengths.push_back(p.second - p.first);});
        }else{continue;}
        ibd0011_tot += std::accumulate(segLengths.begin(), segLengths.end(), decltype(segLengths)::value_type(0));
    }

    // clean up
    for(auto it = union1.begin(); it != union1.end(); it++){delete it->second;}
    for(auto it = union2.begin(); it != union2.end(); it++){delete it->second;}
    return ibd0011_tot;
}



// for debugging
void printConnInfo(const ConnInfo &con, const Pedigree &pedigree)
{
    auto vertex_property_map = boost::get(&sample::id, pedigree);
    if(!con.gp1.empty()){
        fprintf(stdout, "grandparents1:\n");
        for(Vertex u : con.gp1){
            fprintf(stdout, "%s\n", vertex_property_map[u].c_str());
        }
    }

    if(!con.gp2.empty()){
        fprintf(stdout, "grandparents2:\n");
        for(Vertex u : con.gp2){
            fprintf(stdout, "%s\n", vertex_property_map[u].c_str());
        }
    }

    if(!con.av1.empty()){
        fprintf(stdout, "av1:\n");
        for(Vertex u : con.av1){
            fprintf(stdout, "%s\n", vertex_property_map[u].c_str());
        }
    }

    if(!con.av2.empty()){
        fprintf(stdout, "av2:\n");
        for(Vertex u : con.av2){
            fprintf(stdout, "%s\n", vertex_property_map[u].c_str());
        }
    }

    if(!con.p.empty()){
        fprintf(stdout, "parents:\n");
        for(Vertex u : con.p){
            fprintf(stdout, "%s\n", vertex_property_map[u].c_str());
        }
    }

    if(!con.fs.empty()){
        fprintf(stdout, "full-sibs:\n");
        for(Vertex u : con.fs){
            fprintf(stdout, "%s\n", vertex_property_map[u].c_str());
        }
    }

    fprintf(stdout, "\n");

}

bool isFS(Vertex u, const Pedigree &pedigree)
{
    // either u has no polarized edges, or in polarized edges it is the younger generation
    // actually this is problematic if the pedigree spans more than three generations
    // TODO: HANDLE PEDIGREES SPANNING MORE THAN THREE GENERATIONS
    auto out_edge_iter_pair = boost::out_edges(u, pedigree);
    for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
        if (pedigree[*it].polarized){
            if(pedigree[*it].older == u){return false;}
        }
        else if (pedigree[*it].rel != FS){return false;}
    }
    return true;
}

bool isAV(Vertex u, const Pedigree &pedigree){
    // u must have full-sibs as its younger generation
    bool foundFS = false;
    auto out_edge_iter_pair = boost::out_edges(u, pedigree);
    for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
        if (pedigree[*it].rel == AV && pedigree[*it].older == u){
            Vertex s = boost::source(*it, pedigree);
            Vertex t = boost::target(*it, pedigree);
            Vertex v = s == u? t:s;
            if (isFS(v, pedigree)){foundFS = true;}
        }
    }
    return foundFS;
}

bool isP(Vertex u, const Pedigree &pedigree){
    bool foundFS = false;
    auto out_edge_iter_pair = boost::out_edges(u, pedigree);
    for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
        if (pedigree[*it].rel == PC && pedigree[*it].polarized && pedigree[*it].older == u){
            Vertex s = boost::source(*it, pedigree);
            Vertex t = boost::target(*it, pedigree);
            Vertex v = s == u? t:s;
            if (!isFS(v, pedigree)){return false;}
            else{foundFS = true;}
        }
    }
    return foundFS;
}

bool isGP(Vertex u, const Pedigree &pedigree){
    bool foundAVorP = false;
    auto out_edge_iter_pair = boost::out_edges(u, pedigree);
    for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
        if (pedigree[*it].rel == PC && pedigree[*it].polarized && pedigree[*it].older == u){
            Vertex s = boost::source(*it, pedigree);
            Vertex t = boost::target(*it, pedigree);
            Vertex v = s == u? t:s;
            if (isP(v, pedigree) || isAV(v, pedigree)){foundAVorP = true;}
            else{return false;}
        }
    }
    return foundAVorP;
}

void write_output(const std::map<std::pair<Vertex, Vertex>, int> &results, 
    const std::string &prefix, const Pedigree &pedigree, const std::map<Vertex, Vertex> &twins)
{
  std::string outFileName = prefix + ".DRUID";
    FileOrGZ<FILE *> outFile;
    bool ret = outFile.open(outFileName.c_str(), "w");
    if(!ret){
        fprintf(stderr, "cannot open %s for writing output\n", outFileName.c_str());
        exit(1);
    }

    auto vertex_property_map = boost::get(&sample::id, pedigree);
    std::map<int, std::string> rel2string = {
        {0, "PC"},
        {1, "FS"},
        {2, "GP"},
        {3, "AV"}
    };

    boost::graph_traits<Pedigree>::vertex_iterator vi1, vi_end1;
    boost::graph_traits<Pedigree>::vertex_iterator vi2, vi_end2;
    for(boost::tie(vi1, vi_end1) = boost::vertices(pedigree); vi1 != vi_end1; vi1++){
        std::string id1 = vertex_property_map[*vi1];
        Vertex u = twins.find(*vi1) == twins.end()? *vi1 : twins.find(*vi1)->second;
        for(boost::tie(vi2, vi_end2) = boost::vertices(pedigree); vi2 != vi_end2; vi2++){
            std::string id2 = vertex_property_map[*vi2];
            Vertex v = twins.find(*vi2) == twins.end()? *vi2 : twins.find(*vi2)->second;
            if (*vi1 >= *vi2){continue;} // avoid reporting pairs twice
            else if (u == v){
                // u == v means we have twins
                outFile.printf("%s\t%s\t%s\n", id1.c_str(), id2.c_str(), "MZ");
            }
            else if (boost::edge(u, v, pedigree).second){
                if (!pedigree[boost::edge(u, v, pedigree).first].polarized){
                    outFile.printf("%s\t%s\t%s\n", id1.c_str(), id2.c_str(), rel2string[pedigree[boost::edge(u, v, pedigree).first].rel].c_str());
                }else{
                    std::string older = vertex_property_map[pedigree[boost::edge(u, v, pedigree).first].older];
                    outFile.printf("%s\t%s\t%s\t%s\n", id1.c_str(), id2.c_str(), rel2string[pedigree[boost::edge(u, v, pedigree).first].rel].c_str(), older.c_str());
                }
            }else{
                auto p = results.find(make_pair_v(u, v));
                // pairs that don't have segments shared is not stored in the results map, so need to check this
                int deg = p == results.end()? -1 : p->second;
                outFile.printf("%s\t%s\t%d\n", id1.c_str(), id2.c_str(), deg);
            }
        }
    }
    outFile.close();
}

void setRelationshipBetweenTwoSets(const std::vector<Vertex> &set1, const std::vector<Vertex> &set2, 
    std::map<std::pair<Vertex, Vertex>, int> &results, int deg)
{
    // set the degree between all pairs of set1 and set2 to be deg
    // assume deg <= maxDeg (i.e, deg is a valid integer)
    if (set1.empty() || set2.empty()){return;}
    for(Vertex u : set1){
        for(Vertex v : set2){
            //results[make_pair_v(u, v)] = deg;
            setDeg(u, v, deg, results);
        }
    }
}

void setRelationshipBetweenOneSampleAndSet(Vertex u, const std::vector<Vertex> &set,
    std::map<std::pair<Vertex, Vertex>, int> &results, int deg)
{   
    if (set.empty()){return;}
    for(Vertex v : set){
        //results[make_pair_v(u, v)] = deg;
        setDeg(u, v, deg, results);
    }
}

void PCpairVSone(Vertex d, const std::pair<Vertex, Vertex> &pc, const PairIBD &allsegs,
    std::map<std::pair<Vertex, Vertex>, int> &results, int maxDeg)
{
    // deal with an unpolarized PC pair with a single putative distant relative
    auto it1 = results.find(make_pair_v(pc.first, d));
    int d1 = it1 == results.end() ? -1 : it1->second;
    it1 = results.find(make_pair_v(pc.second, d));
    int d2 = it1 == results.end() ? -1 : it1->second;
    auto it2 = allsegs.find(make_pair_v(pc.first, d));
    double k1 = it2 == allsegs.end() ? 0.0 : it2->second->kin;
    it2 = allsegs.find(make_pair_v(pc.second, d));
    double k2 = it2 == allsegs.end() ? 0.0 : it2->second->kin;
    std::pair<bool, Vertex> tuple = polarizeUnpolarPC(pc.first, d1, k1, pc.second, d2, k2);

    if (tuple.first){
        Vertex p = tuple.second;
        Vertex c = p == pc.first ? pc.second : pc.first;
        //results[make_pair_v(c, d)] = resetRelationship(results[make_pair_v(p, d)], 1, maxDeg);
        setDeg(c, d, resetRelationship(getDeg(p, d, results), 1, maxDeg), results);
    }

}



std::pair<bool, Vertex> polarizeUnpolarPC(Vertex v1, int d1, double k1, Vertex v2, int d2, double k2)
{
    // determine in a parent-offspring pair (v1, v2), which one is the older one
    // d1, k1 is the degree estimate and kinship coefficient of v1 to a putative distant relative (could be a real sample, or a reconstructed ungenotpyed parent/grandparent)
    // and similarly d2, k2
    // return (true, v1/v2) if v1, v2 is inferred to be the parent
    // return (false, 0) if neither satisfies the criterion

    if (d1 != -1 && k1 > 0.0 && (0.75 - 0.025*d1)*k1 > k2 && k1 < MAX_PARENT_CHILD_KINSHIP_RATIO*k2){
        // test if v1 is the parent
        return std::make_pair(true, v1);
    }else if (d2 != -1 && k2 > 0.0 && (0.75 - 0.025*d2)*k2 > k1 && k2 < MAX_PARENT_CHILD_KINSHIP_RATIO*k1){
        return std::make_pair(true, v2);
    }else{return std::make_pair(false, 0);}
}

int PCpairVSpedigree(const std::pair<Vertex, Vertex> &pc, const ConnInfo &con,
    std::unordered_set<Vertex> &visited, const PairIBD &allsegs,
    std::map<std::pair<Vertex, Vertex>, int> &results, const Pedigree &pedigree,
    const std::map<std::string, int> &id2index, double bkg_sharing, double tot_genome, int maxDeg)
{
    // test
    //fprintf(stdout, "PCpairVSpedigree\n");
    //auto vertex_property_map = boost::get(&sample::id, pedigree);
    //fprintf(stdout, "%s\n", vertex_property_map[pc.first].c_str());
    //fprintf(stdout, "%s\n", vertex_property_map[pc.second].c_str());
    // end of test

    bool aunt11 = includeAunts(con.av1, con.fs, pc.first, allsegs);
    bool aunt12 = includeAunts(con.av1, con.fs, pc.second, allsegs);
    bool aunt21 = includeAunts(con.av2, con.fs, pc.first, allsegs);
    bool aunt22 = includeAunts(con.av2, con.fs, pc.second, allsegs);
    bool aunt1 = aunt11 || aunt12;
    bool aunt2 = aunt21 || aunt22;
    int index_av = -1;
    if (aunt1 && !aunt2){index_av = 1;}
    else if (!aunt1 && aunt2){index_av = 2;}
    else if (aunt1 && aunt2){
        double average11 = averageKinship(pc.first, con.av1, allsegs);
        double average12 = averageKinship(pc.second, con.av1, allsegs);
        double average21 = averageKinship(pc.first, con.av2, allsegs);
        double average22 = averageKinship(pc.second, con.av2, allsegs);
        double average1 = std::max(average11, average12);
        double average2 = std::max(average21, average22);
        index_av = average1 > average2 ? 1 : 2;
    }

    std::vector<Vertex> s1 {pc.first};
    std::vector<Vertex> s2 {pc.second};
    if (index_av == -1){
        // check if we can use full-sibs' parent
        int p2use = whichParent2Include(con.p, con.fs, s1, allsegs);
        if (p2use == -1){
            p2use = whichParent2Include(con.p, con.fs, s2, allsegs);
        }
        int d1, d2;
        double k1, k2;
        if (p2use != -1){
            Vertex p = con.p[p2use];
            visited.insert(p);
            d1 = getDeg(pc.first, p, results);
            d2 = getDeg(pc.second, p, results);
            auto it2 = allsegs.find(make_pair_v(pc.first, p));
            k1 = it2 == allsegs.end()? 0.0 : it2->second->kin;
            it2 = allsegs.find(make_pair_v(pc.second, p));
            k2 = it2 == allsegs.end()? 0.0 : it2->second->kin;
        }else{
            // use reconstructed parents to calculate d1, d2, k1, k2
            double tg = getTg(0, con.fs.size());
            double tot_ibd1 = std::max(0.0, UnionIbdOverTwoSets(s1, con.fs, id2index, allsegs) - calc_bkg_sharing(0, con.fs.size(), bkg_sharing));
            k1 = (tot_ibd1/(tot_genome*4.0))/tg;
            d1 = getRelfromK(k1, maxDeg);
            double tot_ibd2 = std::max(0.0, UnionIbdOverTwoSets(s2, con.fs, id2index, allsegs) -calc_bkg_sharing(0, con.fs.size(), bkg_sharing));
            k2 = (tot_ibd2/(tot_genome*4.0))/tg;
            d2 = getRelfromK(k2, maxDeg);
        }

        std::pair<bool, Vertex> tuple = polarizeUnpolarPC(pc.first, d1, k1, pc.second, d2, k2);
        if (!tuple.first){
            oneVSpedigree(pc.first, con, visited, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
            oneVSpedigree(pc.second, con, visited, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
        }else{
            Vertex p = tuple.second;
            Vertex c = p == pc.first ? pc.second : pc.first;
            int deg_pp = pc.first == p ? d1 : d2;
            for(auto fs : con.fs){
                //results[make_pair_v(fs, p)] = resetRelationship(deg_pp, 1, maxDeg);
                //results[make_pair_v(fs, c)] = resetRelationship(deg_pp, 2, maxDeg);
                setDeg(fs, p, resetRelationship(deg_pp, 1, maxDeg), results);
                setDeg(fs, c, resetRelationship(deg_pp, 2, maxDeg), results);
                visited.insert(fs);
            }
        }
    }else{
        // check if we can use full-sib's grandparent
        const std::vector<Vertex> &av2use = index_av == 1 ? con.av1 : con.av2;
        const std::vector<Vertex> &gp2check = index_av == 1 ? con.gp1 : con.gp2;
        int index_gp = whichParent2Include(gp2check, av2use, s1, allsegs);
        if (index_gp == -1){
            index_gp = whichParent2Include(gp2check, av2use, s2, allsegs);
        }
        int d1, d2;
        double k1, k2;
        if (index_gp != -1){
            Vertex gp = gp2check[index_gp];
            visited.insert(gp);
            d1 = getDeg(pc.first, gp, results);
            d2 = getDeg(pc.second, gp, results);
            auto it2 = allsegs.find(make_pair_v(pc.first, gp));
            k1 = it2 == allsegs.end()? 0.0 : it2->second->kin;
            it2 = allsegs.find(make_pair_v(pc.second, gp));
            k2 = it2 == allsegs.end()? 0.0 : it2->second->kin;
        }else{
            // use reconstructed grandparents to calculate d1, d2, k1, k2
            std::vector<Vertex> rhs;
            std::for_each(con.fs.begin(), con.fs.end(), [&](Vertex v){rhs.push_back(v);});
            std::for_each(av2use.begin(), av2use.end(), [&](Vertex v){rhs.push_back(v);});
            double tg = getTg(av2use.size(), con.fs.size());
            double tot_ibd1 = std::max(0.0, UnionIbdOverTwoSets(s1, rhs, id2index, allsegs) - calc_bkg_sharing(av2use.size(), con.fs.size(), bkg_sharing));
            k1 = (tot_ibd1/(tot_genome*4.0))/tg;
            d1 = getRelfromK(k1, maxDeg);
            double tot_ibd2 = std::max(0.0, UnionIbdOverTwoSets(s2, rhs, id2index, allsegs) -calc_bkg_sharing(av2use.size(), con.fs.size(), bkg_sharing));
            k2 = (tot_ibd2/(tot_genome*4.0))/tg;
            d2 = getRelfromK(k2, maxDeg);
        }

        std::pair<bool, Vertex> tuple = polarizeUnpolarPC(pc.first, d1, k1, pc.second, d2, k2);
        if (!tuple.first){
            oneVSpedigree(pc.first, con, visited, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
            oneVSpedigree(pc.second, con, visited, allsegs, results, pedigree, id2index, bkg_sharing, tot_genome, maxDeg);
        }else{
            Vertex p = tuple.second;
            Vertex c = p == pc.first ? pc.second : pc.first;
            int deg_pvsgp = pc.first == p ? d1 : d2;
            for(auto av : av2use){
                //results[make_pair_v(av, p)] = resetRelationship(deg_pvsgp, 1, maxDeg);
                //results[make_pair_v(av, c)] = resetRelationship(deg_pvsgp, 2, maxDeg);
                setDeg(av, p, resetRelationship(deg_pvsgp, 1, maxDeg), results);
                setDeg(av, c, resetRelationship(deg_pvsgp, 2, maxDeg), results);
                visited.insert(av);
            }
            
            for(auto fs : con.fs){
                //results[make_pair_v(fs, p)] = resetRelationship(deg_pvsgp, 2, maxDeg);
                //results[make_pair_v(fs, c)] = resetRelationship(deg_pvsgp, 3, maxDeg);
                setDeg(fs, p, resetRelationship(deg_pvsgp, 2, maxDeg), results);
                setDeg(fs, c, resetRelationship(deg_pvsgp, 3, maxDeg), results);
                visited.insert(fs);
            }
        }
    }
    return index_av;
    
}

void grabChildren(Vertex p, std::vector<Vertex> &children, const Pedigree &pedigree)
{
    // store all children of p in the children vector
    auto out_edge_iter_pair = boost::out_edges(p, pedigree);
    for(auto it = out_edge_iter_pair.first; it != out_edge_iter_pair.second; it++){
        Edge e = *it;
        Vertex s = boost::source(e, pedigree);
        Vertex t = boost::target(e, pedigree);
        Vertex other = s == p ? t : s;
        if (pedigree[e].rel == PC && pedigree[e].polarized && pedigree[e].older == p){    
            children.push_back(other);
        }
    }

}

void propagate(Vertex u, Vertex v, std::unordered_set<Vertex> &visited1, 
    std::unordered_set<Vertex> &visited2, int maxDeg, const Pedigree &pedigree, 
    std::map<std::pair<Vertex, Vertex>, int> &results)
{
    int baseDeg = getDeg(u, v, results);
    // set relationship between u and v's descendents given that u,v are baseDeg degrees related
    //if (baseDeg >= maxDeg){return;} Don't return here because I want to add all u,v's descendants to visited1, visited2
    std::vector<Vertex> children_u;
    std::vector<Vertex> children_v;
    grabChildren(u, children_u, pedigree);
    grabChildren(v, children_v, pedigree);
    setRelationshipBetweenOneSampleAndSet(u, children_v, results, resetRelationship(baseDeg, 1, maxDeg));
    setRelationshipBetweenOneSampleAndSet(v, children_u, results, resetRelationship(baseDeg, 1, maxDeg));
    setRelationshipBetweenTwoSets(children_u, children_v, results, resetRelationship(baseDeg, 2, maxDeg));
    for(Vertex child_u : children_u){
        for(Vertex child_v : children_v){
            propagate(child_u, child_v, visited1, visited2, maxDeg, pedigree, results);
            visited1.insert(child_u);
            visited2.insert(child_v);
        }
    }  

}

void propagateAlongPedigree(const ConnInfo &con1, const ConnInfo &con2, int av1, int av2,
    std::unordered_set<Vertex> &visited1, std::unordered_set<Vertex> &visited2,
    const Pedigree &pedigree, std::map<std::pair<Vertex, Vertex>, int> &results, int maxDeg)
{   
    // no need to check con1.p and con2.p because their descendants are just con1.fs, con2.fs
    // and no need to check the gp generation since their descendants are just the AVs 
    // check AV
    std::vector<Vertex> set1;
    std::vector<Vertex> set2;
    copy(con1.fs.begin(), con1.fs.end(), std::back_inserter(set1));
    copy(con2.fs.begin(), con2.fs.end(), std::back_inserter(set2));
    if (av1 != -1){
        const std::vector<Vertex> av2use = av1 == 1 ? con1.av1 : con1.av2;
        copy(av2use.begin(), av2use.end(), std::back_inserter(set1));
    }
    if (av2 != -1){
        const std::vector<Vertex> av2use = av2 == 1 ? con2.av1 : con2.av2;
        copy(av2use.begin(), av2use.end(), std::back_inserter(set2));
    }

    for(Vertex u : set1){
        for(Vertex v : set2){
            propagate(u, v, visited1, visited2, maxDeg, pedigree, results);
        }
    }

}
