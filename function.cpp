#include <math.h>
#include <memory>
#include <numeric>
#include <queue>
#include "function.h"
#include "constants.h"
#include <boost/graph/connected_components.hpp>


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
    std::map<std::pair<Vertex, Vertex>, int> &results,
    std::map<Vertex, Vertex> &twins,
    double tot_genome, double bkg_sharing, int maxDeg)
{
    auto vertex_property_map = boost::get(&sample::id, pedigree);
    boost::graph_traits<Pedigree>::vertex_iterator vi1, vi_end1;
    boost::graph_traits<Pedigree>::vertex_iterator vi2, vi_end2;
    std::set<std::pair<Vertex, Vertex>> pcs;
    std::map<Vertex, std::shared_ptr<std::unordered_set<Vertex>>> fs_degs;
    std::map<Vertex, std::shared_ptr<std::unordered_set<Vertex>>> second_degs;
    for(boost::tie(vi1, vi_end1) = boost::vertices(pedigree); vi1 != vi_end1; vi1++){
        std::string id1 = vertex_property_map[*vi1];
        //std::cout << id1 << "\t" << *vi1 << std::endl;
        for(boost::tie(vi2, vi_end2) = boost::vertices(pedigree); vi2 != vi_end2; vi2++){
            std::string id2 = vertex_property_map[*vi2];
            if (*vi1 >= *vi2){continue;} // avoid analyzing pairs twice
            auto it = allsegs.find(make_pair_str(id1, id2)); 
            // some pairs may not have ibd segments, so need to check here
            if(it == allsegs.end()){continue;};
            Pair &p = *(it->second);
            double ibd1 = p.ibd1_tot;
            double ibd2 = p.ibd2_tot;
            double K = std::max((ibd1/4.0 + ibd2/2.0 - bkg_sharing/4.0)/tot_genome, 0.0);
            p.kin = K;
            int deg = getRelfromK(K, maxDeg);
            //std::cout << id1 << "\t" << id2 << "\t" << K << "\t" << deg << std::endl;
            results.insert(std::make_pair(std::make_pair(*vi1, *vi2), deg));
            // store first and second degree pairs' sample names for later use
            if (deg == 1 || deg == 2){
                bool isFS = true;
                if (deg == 1){
                    if (ibd2/tot_genome >= FULL_SIB_MIN_IBD2){
                        boost::add_edge(*vi1, *vi2, pedigree);
                        pedigree[boost::edge(*vi1, *vi2, pedigree).first].rel = FS;
                        //std::cout << id1 << " and " << id2  << " is inferred to be FS" << std::endl;
                    }else{
                        isFS = false;
                        pcs.insert(std::make_pair(*vi1, *vi2));
                    }
                }else if(deg == 2){
                    // check for possibility of DC, if so, no need to consider this pair for AV
                    // therefore no need to add them to second_deg
                    if (ibd2/tot_genome >= DC_MIN_IBD2){continue;}
                }

                auto &map = deg == 1? fs_degs : second_degs;
                if (isFS || deg == 2){
                    if(map.find(*vi1) == map.end()){
                        map[*vi1] = std::unique_ptr<std::unordered_set<Vertex>>(new std::unordered_set<Vertex>());
                    }
                    map[*vi1]->insert(*vi2);
                    if(map.find(*vi2) == map.end()){
                        map[*vi2] = std::unique_ptr<std::unordered_set<Vertex>>(new std::unordered_set<Vertex>());
                    }
                    map[*vi2]->insert(*vi1);
                }
                
            }else if (deg == 0){
                twins.insert(std::make_pair(*vi1, *vi2));
            }
        }
    }

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
                        // let's update the results map when we actually write the output, or shall we?
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
        // if neither v1 nor v2 has any full-sibs, we can't do anything
        if (fs_degs.find(v1) == fs_degs.end() && fs_degs.find(v2) == fs_degs.end()){continue;}
        bool v2IsParent = false;
        bool v1IsParent = false;
        if (fs_degs.find(v1) != fs_degs.end()){
            // iterate over v1's full-sib pair to check if they form a parent-child pair with v2
            // if so, then v2 must be the parent
            for(auto fs : *fs_degs[v1]){
                bool isConnected = boost::edge(fs, v2, pedigree).second;
                //std::cout << "check v1's full-sib: " << vertex_property_map[fs] << " is connected? "<< isConnected << std::endl;
                if (isConnected){
                    Edge e = boost::edge(fs, v2, pedigree).first;
                    //std::cout << pedigree[e].rel << std::endl;
                    if (pedigree[e].rel == PC){v2IsParent=true;}
                }
            }

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
            }
        }

        if (fs_degs.find(v2) != fs_degs.end()){
            // iterate over v2's full-sib pair to check if they form a parent-child pair with v1
            // if so, then v1 must be the parent
            for(auto fs : *fs_degs[v2]){
                bool isConnected = boost::edge(fs, v1, pedigree).second;
                //std::cout << "check v2's full-sib: " << vertex_property_map[fs] << " is connected? "<< isConnected << std::endl;
                if (isConnected){
                    Edge e = boost::edge(fs, v1, pedigree).first;
                    //std::cout << pedigree[e].rel << std::endl;
                    if (pedigree[e].rel == PC){v1IsParent=true;}
                }
            }

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
            }
        }
        assert(!(v1IsParent && v2IsParent));
    }

    // construct second-degree edges
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
        
        const std::unordered_set<Vertex> &fs_to_focal_ind = *(fs_degs[focal_ind]);
        // consider all pairs of full-siblings of the focal individual
        std::unordered_set<Vertex> full_sib_vertex_set = fs_to_focal_ind; // copy
        full_sib_vertex_set.insert(focal_ind);
        std::vector<std::string> full_sib_id_set(full_sib_vertex_set.size());
        std::transform(full_sib_vertex_set.begin(), full_sib_vertex_set.end(),
            full_sib_id_set.begin(), [&](Vertex v){return vertex_property_map[v];});

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
            std::string avunc_candidate_id = vertex_property_map[avunc_candidate];
            //std::cout << "checking " << avunc_candidate_id << std::endl;
            if (checkAvunc(full_sib_id_set, avunc_candidate_id, allsegs, snpmap)){
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

bool is_avunc(const std::string &fs1, const std::string &fs2, const std::string &avunc, 
        const std::map<std::pair<std::string, std::string>, Pair*> &allsegs, 
        const std::map<std::string, std::map<int, double>*> &snpmap)
{
    auto fs_p = allsegs.find(make_pair_str(fs1, fs2));
    assert(fs_p != allsegs.end());
    Pair &full_sib_pair = *(fs_p->second);
    auto fa_p1 = allsegs.find(make_pair_str(fs1, avunc));
    assert(fa_p1 != allsegs.end());
    Pair &sib1_avunc_pair = *(fa_p1->second);
    auto fa_p2 = allsegs.find(make_pair_str(fs2, avunc));
    assert(fa_p2 != allsegs.end());
    Pair &sib2_avunc_pair = *(fa_p2->second);

    double ibd011_tot = 0.0;
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
    return ibd011_tot >= AVUNC_011;
}

bool checkAvunc(const std::vector<std::string> &full_sibs, const std::string &avunc, 
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap)
{
    // check if avunc is the avunc of the given set of full_sibs
    int num_sibs = full_sibs.size();
    for(int i = 0; i < num_sibs; i++){
        for(int j = i+1; j < num_sibs; j++){
            if(is_avunc(full_sibs[i], full_sibs[j], avunc, allsegs, snpmap)){return true;}
        }
    }
    return false;
}

void run_druid(Pedigree &pedigree, 
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
    const std::map<std::string, std::map<int, double>*> &snpmap,
    std::map<std::pair<Vertex, Vertex>, int> &results,
    FileOrGZ<FILE *> &logFile,
    double tot_genome, double bkg_sharing, int maxDeg)
{   
    logFile.printf("Start the primary DRUID algorithm...\n");
    std::vector<int> components(boost::num_vertices(pedigree));
    int num_components = boost::connected_components(pedigree, &components[0]);
    logFile.printf("\tnumber of connected components: %d\n", num_components);

    //std::map<int, int> comp_size_map;
    std::map<int, std::shared_ptr<std::vector<Vertex>>> comp_map;
    boost::graph_traits<Pedigree>::vertex_iterator vi, vi_end;
    for(boost::tie(vi, vi_end) = boost::vertices(pedigree); vi != vi_end; vi++){
        int comp_index = components[*vi];
        if (comp_map.find(comp_index) == comp_map.end()){
            //comp_size_map.insert(std::make_pair(comp_index, 0));
            comp_map.insert(std::make_pair(comp_index, std::shared_ptr<std::vector<Vertex>>(new std::vector<Vertex>())));
        }
        //comp_size_map[comp_index]++;
        comp_map[comp_index]->push_back(*vi);
    }

    auto vertex_property_map = boost::get(&sample::id, pedigree);
    for(int i = 0; i < num_components; i++){
        auto ordered = std::shared_ptr<std::vector<Vertex>>(new std::vector<Vertex>());
        postorder(*(comp_map.find(i)->second), pedigree, *ordered);
        comp_map[i] = ordered;
    }

    for(int i = 0; i < num_components; i++){
        for(int j = i+1; j < num_components; j++){
            std::unordered_set<Vertex> visited1;
            std::unordered_set<Vertex> visited2;
            for(Vertex u : *comp_map[i]){
                if (visited1.find(u) != visited1.end()){continue;}
                ConnInfo con1;
                grabCloseRelatives(u, con1, pedigree);
                bool isSingleton1 = isSingleton(con1);
                for(Vertex v : *comp_map[j]){
                    if(visited2.find(v) != visited2.end()){continue;}
                    ConnInfo con2;
                    grabCloseRelatives(v, con2, pedigree);
                    // analyzing the two ConInfo component
                    bool isSingleton2 = isSingleton(con2);
                    if (isSingleton1 && isSingleton2){
                        visited1.insert(u);
                        visited2.insert(v);
                        continue;}
                    else if(isSingleton1 && !isSingleton2){
                        visited1.insert(u);
                        oneVSpedigree(u, con2, visited2, allsegs, results, pedigree, bkg_sharing, tot_genome, maxDeg);
                    }else if(!isSingleton1 && isSingleton2){
                        visited2.insert(v);
                        oneVSpedigree(v, con1, visited1, allsegs, results, pedigree, bkg_sharing, tot_genome, maxDeg);
                    }else{

                    }

                }
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

bool isSingleton(const ConnInfo &con)
{
    return (con.gp1.empty() && con.gp2.empty() 
    && con.av1.empty() && con.av2.empty() && con.p.empty() &&
    con.fs.size() == 1);
}

double UnionIbdOverTwoSets(const std::vector<std::string> &set1, const std::vector<std::string> &set2,
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs)
{
    ibdMapType currUnion;
    for(auto s1 : set1){
        for(auto s2 : set2){
            //fprintf(stdout, "%s\t%s\n", s1.c_str(), s2.c_str());
            auto p = allsegs.find(make_pair_str(s1, s2));
            if (p == allsegs.end()){continue;}
            else{
                const Pair &pair = *(p->second);
                for(auto it = pair.ibd1_map->begin(); it != pair.ibd1_map->end(); it++){
                    std::string chrName = it->first;
                    //fprintf(stdout, "chrom Name: %s\n", chrName.c_str());
                    if (currUnion.find(chrName) == currUnion.end()){
                        currUnion.insert(std::make_pair(chrName, new ibdSegments()));
                    }
                    ibdSegments *dest = new ibdSegments();
                    interval_union(*(it->second), *(currUnion.find(chrName)->second), *dest);
                    ibdSegments *prev_ptr = currUnion.find(chrName)->second;
                    currUnion[chrName] = dest;
                    delete prev_ptr;
                    //currUnion.insert(std::make_pair(chrName, dest));
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


void oneVSpedigree(Vertex u, const ConnInfo &con, std::unordered_set<Vertex> &visited,
    const std::map<std::pair<std::string, std::string>, Pair*> &allsegs,
    std::map<std::pair<Vertex, Vertex>, int> &results, 
    const Pedigree &pedigree, double bkg_sharing, double tot_genome, int maxDeg)
{   
    auto vertex_property_map = boost::get(&sample::id, pedigree);
    if(con.av1.empty() && con.av2.empty() && con.p.empty()){
        std::vector<std::string> set1;
        std::vector<std::string> set2;
        set1.push_back(vertex_property_map[u]);
        for(Vertex v : con.fs){set2.push_back(vertex_property_map[v]);}
        double unionLength = UnionIbdOverTwoSets(set1, set2, allsegs);
        int numSibs = con.fs.size();
        unionLength = std::max(0.0, unionLength - 2*bkg_sharing*(1.0 - pow(0.5, numSibs)));
        double K = (unionLength/tot_genome)/(1 - pow(0.5, numSibs));
        int tmp = getRelfromK(K, maxDeg);
	    int deg = tmp == -1 ?  -1 : tmp+1;
	    if(deg > maxDeg){deg = -1;}
        //fprintf(stdout, "new degree is %d\n", deg);
        for(Vertex v : con.fs){
            results[make_pair_v(u,v)] = deg;
            visited.insert(v);
        }
    }
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

    fprintf(stdout, "\n\n");

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
                outFile.printf("%s\t%s\t%s\n", id1.c_str(), id2.c_str(), rel2string[pedigree[boost::edge(u, v, pedigree).first].rel].c_str());
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
