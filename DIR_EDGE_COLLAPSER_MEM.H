#define DIR_EDGE_COLLAPSER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <bitset>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <memory>
#include <queue>

using Vertex = int;
using Filtration_value = double;
using Edge = std::tuple<Vertex, Vertex, Filtration_value>;
constexpr Filtration_value infinity = std::numeric_limits<double>::infinity();

// Memory-efficient vertex set using bitset
constexpr size_t MAX_VERTICES = 10000;
using Vertex_set = std::bitset<MAX_VERTICES>;
using Neighborhood = std::unordered_map<Vertex, Filtration_value>;
using Vertex_set_tuple = std::tuple<Vertex_set, Vertex_set, Vertex_set>;

class SetOps {
public:
    static Vertex_set setUnion(const Vertex_set& set1, const Vertex_set& set2) {
        return set1 | set2;
    }

    static Vertex_set setIntersection(const Vertex_set& set1, const Vertex_set& set2) {
        return set1 & set2;
    }

    static bool isSubset(const Vertex_set& set1, const Vertex_set& set2) {
        return (set1 & ~set2).none();
    }
};

// Custom hash for pair
struct PairHash {
    std::size_t operator()(const std::pair<Vertex, Vertex>& p) const {
        auto h1 = std::hash<Vertex>{}(p.first);
        auto h2 = std::hash<Vertex>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

// Cache for neighborhood computations - adapted for bitsets
class NeighborhoodCache {
private:
    struct CacheEntry {
        Filtration_value max_filt;
        Vertex_set neighbors;
    };
    
    mutable std::unordered_map<Vertex, std::vector<CacheEntry>> outgoing_cache;
    mutable std::unordered_map<Vertex, std::vector<CacheEntry>> incoming_cache;
    
public:
    void clear() {
        outgoing_cache.clear();
        incoming_cache.clear();
    }
    
    void invalidate(Vertex v) {
        outgoing_cache.erase(v);
        incoming_cache.erase(v);
    }
    
    const Vertex_set* get_cached_outgoing(Vertex v, Filtration_value max_filt) const {
        auto it = outgoing_cache.find(v);
        if (it != outgoing_cache.end()) {
            for (const auto& entry : it->second) {
                if (entry.max_filt == max_filt) {
                    return &entry.neighbors;
                }
            }
        }
        return nullptr;
    }
    
    void cache_outgoing(Vertex v, Filtration_value max_filt, Vertex_set neighbors) const {
        outgoing_cache[v].push_back({max_filt, std::move(neighbors)});
    }
    
    const Vertex_set* get_cached_incoming(Vertex v, Filtration_value max_filt) const {
        auto it = incoming_cache.find(v);
        if (it != incoming_cache.end()) {
            for (const auto& entry : it->second) {
                if (entry.max_filt == max_filt) {
                    return &entry.neighbors;
                }
            }
        }
        return nullptr;
    }
    
    void cache_incoming(Vertex v, Filtration_value max_filt, Vertex_set neighbors) const {
        incoming_cache[v].push_back({max_filt, std::move(neighbors)});
    }
};

class FilteredGraph {
private:
    std::unordered_map<Vertex, std::vector<std::pair<Vertex, Filtration_value>>> outgoing_edges;
    std::unordered_map<Vertex, std::vector<std::pair<Vertex, Filtration_value>>> incoming_edges;
    
    // Cache for neighborhood computations
    mutable NeighborhoodCache cache;
    
public:
    // Quick edge lookup (made public for access)
    std::unordered_map<Vertex, std::unordered_map<Vertex, Filtration_value>> edge_map;
    
public:
    void add_edge(Vertex u, Vertex v, Filtration_value t) {
        // Update cache invalidation
        cache.invalidate(u);
        cache.invalidate(v);
        
        // Update edge storage
        auto& out_edges = outgoing_edges[u];
        auto& in_edges = incoming_edges[v];
        
        // Check if edge exists and update
        bool found = false;
        for (auto& [vertex, filt] : out_edges) {
            if (vertex == v) {
                filt = t;
                found = true;
                break;
            }
        }
        if (!found) {
            out_edges.emplace_back(v, t);
        }
        
        found = false;
        for (auto& [vertex, filt] : in_edges) {
            if (vertex == u) {
                filt = t;
                found = true;
                break;
            }
        }
        if (!found) {
            in_edges.emplace_back(u, t);
        }
        
        // Update quick lookup
        edge_map[u][v] = t;
    }
    
    void remove_edge(Vertex u, Vertex v) {
        cache.invalidate(u);
        cache.invalidate(v);
        
        if (auto it = outgoing_edges.find(u); it != outgoing_edges.end()) {
            auto& edges = it->second;
            edges.erase(std::remove_if(edges.begin(), edges.end(),
                [v](const auto& p) { return p.first == v; }), edges.end());
            if (edges.empty()) outgoing_edges.erase(it);
        }
        
        if (auto it = incoming_edges.find(v); it != incoming_edges.end()) {
            auto& edges = it->second;
            edges.erase(std::remove_if(edges.begin(), edges.end(),
                [u](const auto& p) { return p.first == u; }), edges.end());
            if (edges.empty()) incoming_edges.erase(it);
        }
        
        if (auto it = edge_map.find(u); it != edge_map.end()) {
            it->second.erase(v);
            if (it->second.empty()) edge_map.erase(it);
        }
    }
    
    Vertex_set get_outgoing_neighbors_set(Vertex u, Filtration_value max_filtration) const {
        // Check cache first
        if (const auto* cached = cache.get_cached_outgoing(u, max_filtration)) {
            return *cached;
        }
        
        Vertex_set neighbors;
        auto it = outgoing_edges.find(u);
        if (it != outgoing_edges.end()) {
            for (const auto& [v, t] : it->second) {
                if (t <= max_filtration && v < MAX_VERTICES) {
                    neighbors.set(v);
                }
            }
        }
        
        // Cache the result
        cache.cache_outgoing(u, max_filtration, neighbors);
        return neighbors;
    }
    
    Vertex_set get_incoming_neighbors_set(Vertex v, Filtration_value max_filtration) const {
        // Check cache first
        if (const auto* cached = cache.get_cached_incoming(v, max_filtration)) {
            return *cached;
        }
        
        Vertex_set neighbors;
        auto it = incoming_edges.find(v);
        if (it != incoming_edges.end()) {
            for (const auto& [u, t] : it->second) {
                if (t <= max_filtration && u < MAX_VERTICES) {
                    neighbors.set(u);
                }
            }
        }
        
        // Cache the result
        cache.cache_incoming(v, max_filtration, neighbors);
        return neighbors;
    }
    
    void get_neighborhoods(Vertex u, Vertex v, Filtration_value max_filtration,
                          Vertex_set_tuple& current_neighborhood,
                          std::vector<std::tuple<Filtration_value, Vertex, int>>& future_neighborhood,
                          std::vector<std::tuple<Filtration_value, std::pair<std::pair<Vertex, int>, std::pair<Vertex, int>>>>& future_edges) const {
        
        // Get all neighbors at once
        const auto out_u = outgoing_edges.find(u);
        const auto out_v = outgoing_edges.find(v);
        const auto in_u = incoming_edges.find(u);
        const auto in_v = incoming_edges.find(v);
        
        std::unordered_map<Vertex, std::pair<Filtration_value, int>> common_map;
        
        // Process common incoming
        if (in_u != incoming_edges.end() && in_v != incoming_edges.end()) {
            for (const auto& [w, t1] : in_u->second) {
                if (w >= MAX_VERTICES) continue;
                for (const auto& [w2, t2] : in_v->second) {
                    if (w == w2) {
                        Filtration_value f = std::max(t1, t2);
                        if (f <= max_filtration) {
                            std::get<0>(current_neighborhood).set(w);
                        } else {
                            future_neighborhood.emplace_back(f, w, 0);
                        }
                        common_map[w] = {f, 0};
                        break;
                    }
                }
            }
        }
        
        // Process common outgoing
        if (out_u != outgoing_edges.end() && out_v != outgoing_edges.end()) {
            for (const auto& [w, t1] : out_u->second) {
                if (w >= MAX_VERTICES) continue;
                for (const auto& [w2, t2] : out_v->second) {
                    if (w == w2) {
                        Filtration_value f = std::max(t1, t2);
                        if (f <= max_filtration) {
                            std::get<1>(current_neighborhood).set(w);
                        } else {
                            future_neighborhood.emplace_back(f, w, 1);
                        }
                        common_map[w] = {f, 1};
                        break;
                    }
                }
            }
        }
        
        // Process outgoing-incoming
        if (out_u != outgoing_edges.end() && in_v != incoming_edges.end()) {
            for (const auto& [w, t1] : out_u->second) {
                if (w >= MAX_VERTICES) continue;
                for (const auto& [w2, t2] : in_v->second) {
                    if (w == w2) {
                        Filtration_value f = std::max(t1, t2);
                        if (f <= max_filtration) {
                            std::get<2>(current_neighborhood).set(w);
                        } else {
                            future_neighborhood.emplace_back(f, w, 2);
                        }
                        common_map[w] = {f, 2};
                        break;
                    }
                }
            }
        }
        
        for (const auto& [w1, info1] : common_map) {
            for (const auto& [w2, info2] : common_map) {
                if (w1 >= w2) continue;
                
                // Check if edge exists
                if (auto it1 = edge_map.find(w1); it1 != edge_map.end()) {
                    if (auto it2 = it1->second.find(w2); it2 != it1->second.end()) {
                        if (it2->second > max_filtration) {
                            future_edges.emplace_back(it2->second,
                                std::make_pair(std::make_pair(w1, info1.second), 
                                              std::make_pair(w2, info2.second)));
                        }
                    }
                }
                if (auto it1 = edge_map.find(w2); it1 != edge_map.end()) {
                    if (auto it2 = it1->second.find(w1); it2 != it1->second.end()) {
                        if (it2->second > max_filtration) {
                            future_edges.emplace_back(it2->second,
                                std::make_pair(std::make_pair(w2, info2.second), 
                                              std::make_pair(w1, info1.second)));
                        }
                    }
                }
            }
        }
        
        std::sort(future_edges.begin(), future_edges.end());
    }
    
    int is_connected(const Vertex& u, const Vertex& v, Filtration_value max_filtration) const {
        if (auto it1 = edge_map.find(u); it1 != edge_map.end()) {
            if (auto it2 = it1->second.find(v); it2 != it1->second.end()) {
                if (it2->second <= max_filtration) return 1;
            }
        }
        if (auto it1 = edge_map.find(v); it1 != edge_map.end()) {
            if (auto it2 = it1->second.find(u); it2 != it1->second.end()) {
                if (it2->second <= max_filtration) return 0;
            }
        }
        return -1;
    }
    
    void print_edges() const {
        std::cout << "Edges in the graph (with filtration values):" << std::endl;
        for (const auto& [u, edges] : outgoing_edges) {
            for (const auto& [v, t] : edges) {
                std::cout << "(" << u << " -> " << v << ", filtration: " << t << ")" << std::endl;
            }
        }
    }
    
    // Get neighbors for compatibility
    Neighborhood get_outgoing_neighbors(Vertex u, Filtration_value max_filtration) const {
        Neighborhood neighbors;
        if (auto it = outgoing_edges.find(u); it != outgoing_edges.end()) {
            for (const auto& [v, t] : it->second) {
                if (t <= max_filtration) {
                    neighbors[v] = t;
                }
            }
        }
        return neighbors;
    }
    
    Neighborhood get_incoming_neighbors(Vertex v, Filtration_value max_filtration) const {
        Neighborhood neighbors;
        if (auto it = incoming_edges.find(v); it != incoming_edges.end()) {
            for (const auto& [u, t] : it->second) {
                if (t <= max_filtration) {
                    neighbors[u] = t;
                }
            }
        }
        return neighbors;
    }
    
    std::pair<Vertex_set, std::vector<std::pair<Vertex, Filtration_value>>> 
    get_common_neighbors_split(const Neighborhood& nbr1, const Neighborhood& nbr2, 
                              Filtration_value max_filtration) const {
        Vertex_set curr_common_nbr;
        std::vector<std::pair<Vertex, Filtration_value>> future_common_nbr;
        
        for (const auto& [u, t1] : nbr1) {
            if (auto it = nbr2.find(u); it != nbr2.end()) {
                Filtration_value f = std::max(t1, it->second);
                if (f <= max_filtration && u < MAX_VERTICES) {
                    curr_common_nbr.set(u);
                } else if (u < MAX_VERTICES) {
                    future_common_nbr.emplace_back(u, f);
                }
            }
        }
        return {curr_common_nbr, future_common_nbr};
    }
};

class CoreDirFlagFiltration {
private:
    std::vector<Edge> edges;
    FilteredGraph graph;
    
    std::unordered_map<std::pair<Vertex, Vertex>, std::vector<std::pair<Vertex, int>>, 
                      PairHash> dominator_cache;
    
public:
    CoreDirFlagFiltration(const std::vector<Edge>& input_edges) : edges(input_edges) {
        // Add all edges to the graph
        for (const auto& e : edges) {
            Vertex u = std::get<0>(e);
            Vertex v = std::get<1>(e);
            Filtration_value t = std::get<2>(e);
            graph.add_edge(u, v, t);
        }
    }
    
    void process_edges(std::vector<Edge>& Ec) {
        // Sort edges in descending order of their filtration values
        std::sort(edges.begin(), edges.end(), [](const Edge& e1, const Edge& e2) {
            return std::get<2>(e1) > std::get<2>(e2);
        });
        
        // Process each edge
        for (const auto& e : edges) {
            Vertex u = std::get<0>(e);
            Vertex v = std::get<1>(e);
            Filtration_value t_init = std::get<2>(e);
            
            Vertex_set_tuple e_nbd_ioc_tuple;
            std::vector<std::tuple<Filtration_value, Vertex, int>> e_nbd_later;
            std::vector<std::tuple<Filtration_value, std::pair<std::pair<Vertex, int>, std::pair<Vertex, int>>>> edges_later;
            
            graph.get_neighborhoods(u, v, t_init, e_nbd_ioc_tuple, e_nbd_later, edges_later);
            
            std::vector<std::pair<Vertex, int>> dominators;
            is_dominated(u, v, e_nbd_ioc_tuple, dominators, t_init);
            
            if (!dominators.empty()) {
                // Process with domination logic
                process_dominated_edge(u, v, t_init, e_nbd_ioc_tuple, e_nbd_later, 
                                     edges_later, dominators, Ec);
            } else {
                Ec.emplace_back(u, v, t_init);
            }
        }
        
        std::cout << "Final " << Ec.size() << " filtered directed edges (Ec):" << std::endl;
    }
    
private:
    void process_dominated_edge(Vertex u, Vertex v, Filtration_value t_init,
                               Vertex_set_tuple& e_nbd_ioc_tuple,
                               std::vector<std::tuple<Filtration_value, Vertex, int>>& e_nbd_later,
                               std::vector<std::tuple<Filtration_value, std::pair<std::pair<Vertex, int>, std::pair<Vertex, int>>>>& edges_later,
                               std::vector<std::pair<Vertex, int>>& dominators,
                               std::vector<Edge>& Ec) {
        
        auto cmpareT = [](const auto& a, const auto& b) { return std::get<0>(a) > std::get<0>(b); };
        std::make_heap(e_nbd_later.begin(), e_nbd_later.end(), cmpareT);
        std::make_heap(edges_later.begin(), edges_later.end(), cmpareT);
        
        auto e_nbd_later_begin = e_nbd_later.begin();
        auto e_nbd_later_end = e_nbd_later.end();
        auto edges_later_begin = edges_later.begin();
        auto edges_later_end = edges_later.end();
        
        while (true) {
            if ((e_nbd_later_begin == e_nbd_later_end) && (edges_later_begin == edges_later_end)) {
                graph.remove_edge(u, v);
                break;
            }
            
            bool process_neighbor = false;
            bool process_edge = false;
            Filtration_value t_prime;
            
            // Determine what to process next
            if (e_nbd_later_begin != e_nbd_later_end && edges_later_begin == edges_later_end) {
                process_neighbor = true;
                t_prime = std::get<0>(*e_nbd_later_begin);
            } else if (e_nbd_later_begin == e_nbd_later_end && edges_later_begin != edges_later_end) {
                process_edge = true;
                t_prime = std::get<0>(*edges_later_begin);
            } else {
                Filtration_value t_nbr = std::get<0>(*e_nbd_later_begin);
                Filtration_value t_edge = std::get<0>(*edges_later_begin);
                if (t_nbr <= t_edge) {
                    process_neighbor = true;
                    t_prime = t_nbr;
                } else {
                    process_edge = true;
                    t_prime = t_edge;
                }
            }
            
            if (process_neighbor) {
                Vertex w = std::get<1>(*e_nbd_later_begin);
                int type = std::get<2>(*e_nbd_later_begin);
                
                if (still_dominated(e_nbd_ioc_tuple, dominators, w, type, t_prime)) {
                    std::pop_heap(e_nbd_later_begin, e_nbd_later_end--, cmpareT);
                    // Update neighborhood
                    switch (type) {
                        case 0: std::get<0>(e_nbd_ioc_tuple).set(w); break;
                        case 1: std::get<1>(e_nbd_ioc_tuple).set(w); break;
                        case 2: std::get<2>(e_nbd_ioc_tuple).set(w); break;
                    }
                } else {
                    Ec.emplace_back(u, v, t_prime);
                    graph.add_edge(u, v, t_prime);
                    break;
                }
            } else if (process_edge) {
                Vertex x = std::get<0>(std::get<1>(*edges_later_begin)).first;
                int type_x = std::get<0>(std::get<1>(*edges_later_begin)).second;
                Vertex y = std::get<1>(std::get<1>(*edges_later_begin)).first;
                int type_y = std::get<1>(std::get<1>(*edges_later_begin)).second;
                
                if (still_dominated_edge(dominators, x, type_x, y, type_y, t_prime)) {
                    std::pop_heap(edges_later_begin, edges_later_end--, cmpareT);
                } else {
                    Ec.emplace_back(u, v, t_prime);
                    graph.add_edge(u, v, t_prime);
                    break;
                }
            }
        }
    }
    
    bool is_3_simplex(const Vertex& dom, const int& domType, const Vertex& w, 
                     const int& type_w, Filtration_value& t_prime, int& con) const {
        con = graph.is_connected(dom, w, t_prime);
        if (con == -1) return false;
        
        static const bool simplex_table[3][3][2] = {
            // domType = 0 (in)
            {{true, true}, {false, true}, {false, true}},
            // domType = 1 (out)
            {{true, false}, {true, true}, {true, false}},
            // domType = 2 (center)
            {{true, false}, {false, true}, {true, true}}
        };
        
        if (domType >= 0 && domType < 3 && type_w >= 0 && type_w < 3 && con >= 0 && con < 2) {
            return simplex_table[domType][type_w][con];
        }
        return false;
    }
    
    bool is_3_simplex_edge(const Vertex& x, const int& type_x, 
                          const Vertex& y, const int& type_y) const {
        if (type_x == type_y) return true;
        if (type_x == 0 && (type_y == 1 || type_y == 2)) return true;
        if (type_x == 2 && type_y == 1) return true;
        return false;
    }
    
    bool still_dominated_edge(std::vector<std::pair<Vertex, int>>& dominators,
                             const Vertex& x, const int& type_x,
                             const Vertex& y, const int& type_y,
                             Filtration_value& t_prime_edge) const {
        if (dominators.empty()) return false;
        
        // Early return if edge involves a dominator
        for (const auto& [dom, domType] : dominators) {
            if (dom == x || dom == y) return true;
        }
        
        if (!is_3_simplex_edge(x, type_x, y, type_y)) return false;
        
        // Check connections and update dominators
        auto it = dominators.begin();
        while (it != dominators.end()) {
            const auto& [dom, domType] = *it;
            int dxcon = graph.is_connected(dom, x, t_prime_edge);
            int dycon = graph.is_connected(dom, y, t_prime_edge);
            
            if (dxcon == 1 && dycon == 0) {
                it = dominators.erase(it);
            } else {
                ++it;
            }
        }
        
        return !dominators.empty();
    }
    
    bool still_dominated(const Vertex_set_tuple& e_nbd_ioc_tuple,
                        std::vector<std::pair<Vertex, int>>& dominators,
                        const Vertex& w, const int& type_w,
                        Filtration_value& t_prime) const {
        if (dominators.empty()) return false;
        
        auto it = dominators.begin();
        while (it != dominators.end()) {
            const auto& [dom, domType] = *it;
            int dwcon;
            
            if (!is_3_simplex(dom, domType, w, type_w, t_prime, dwcon)) {
                it = dominators.erase(it);
                continue;
            }
            
            bool still_valid = true;
            
            // Check against appropriate neighborhoods based on type
            const Vertex_set* neighborhoods[] = {
                &std::get<0>(e_nbd_ioc_tuple),  // incoming
                &std::get<1>(e_nbd_ioc_tuple),  // outgoing
                &std::get<2>(e_nbd_ioc_tuple)   // center
            };
            
            for (int nbr_type = 0; nbr_type < 3; ++nbr_type) {
                // Iterate over bitset
                for (size_t nbr = 0; nbr < MAX_VERTICES; ++nbr) {
                    if (!(*neighborhoods[nbr_type])[nbr]) continue;
                    if (nbr == dom) continue;
                    
                    int nbr_w_con;
                    if (is_3_simplex(nbr, nbr_type, w, type_w, t_prime, nbr_w_con)) {
                        int nbr_dom_con = graph.is_connected(nbr, dom, t_prime);
                        if (nbr_dom_con == dwcon && dwcon != nbr_w_con) {
                            still_valid = false;
                            break;
                        }
                    }
                }
                if (!still_valid) break;
            }
            
            if (!still_valid) {
                it = dominators.erase(it);
            } else {
                ++it;
            }
        }
        
        return !dominators.empty();
    }
    
    bool split_edge_condition(const Vertex_set& N_x, const Vertex_set& N_y) const {
        if (N_x.none() || N_y.none()) return true;
        
        for (size_t x = 0; x < MAX_VERTICES; ++x) {
            if (!N_x[x]) continue;
            if (auto it = graph.edge_map.find(x); it != graph.edge_map.end()) {
                for (size_t y = 0; y < MAX_VERTICES; ++y) {
                    if (!N_y[y]) continue;
                    if (it->second.count(y)) return false;
                }
            }
        }
        return true;
    }
    
    void is_dominated(const Vertex& u, const Vertex& v, const Vertex_set_tuple& e_nbd_tuple,
                     std::vector<std::pair<Vertex, int>>& dominators,
                     Filtration_value max_filtration) const {
        // Check cache first
        auto key = std::make_pair(u, v);
        if (auto it = dominator_cache.find(key); it != dominator_cache.end()) {
            dominators = it->second;
            return;
        }
        
        is_in_dominated(u, v, e_nbd_tuple, dominators, max_filtration);
        is_center_dominated(u, v, e_nbd_tuple, dominators, max_filtration);
        is_out_dominated(u, v, e_nbd_tuple, dominators, max_filtration);
        
        // Cache result
        const_cast<decltype(dominator_cache)&>(dominator_cache)[key] = dominators;
    }
    
    void is_in_dominated(const Vertex& u, const Vertex& v, const Vertex_set_tuple& e_nbd_tuple,
                        std::vector<std::pair<Vertex, int>>& dominators,
                        Filtration_value max_filtration) const {
        const auto& commonIncoming = std::get<0>(e_nbd_tuple);
        const auto& commonOutgoing = std::get<1>(e_nbd_tuple);
        const auto& outgoingIncoming = std::get<2>(e_nbd_tuple);
        
        // Precompute neighborhoods
        const auto incomingU = graph.get_incoming_neighbors_set(u, max_filtration);
        const auto incomingV = graph.get_incoming_neighbors_set(v, max_filtration);
        const auto incomingUV = SetOps::setIntersection(incomingU, incomingV);
        const auto unionOutgoing = SetOps::setUnion(commonOutgoing, outgoingIncoming);
        
        for (size_t vPrime = 0; vPrime < MAX_VERTICES; ++vPrime) {
            if (!commonIncoming[vPrime]) continue;
            
            const auto incomingVPrime = graph.get_incoming_neighbors_set(vPrime, max_filtration);
            const auto outgoingVPrime = graph.get_outgoing_neighbors_set(vPrime, max_filtration);
            
            // Check conditions
            if (!SetOps::isSubset(unionOutgoing, outgoingVPrime)) continue;
            
            auto unionIncoming = incomingVPrime;
            unionIncoming |= outgoingVPrime;
            unionIncoming.set(vPrime);
            if (!SetOps::isSubset(commonIncoming, unionIncoming)) continue;
            
            const auto commonIncomingVpUV = SetOps::setIntersection(incomingVPrime, incomingUV);
            const auto outgoingIncomingVpUV = SetOps::setIntersection(outgoingVPrime, incomingUV);
            
            if (!split_edge_condition(outgoingIncomingVpUV, commonIncomingVpUV)) continue;
            
            dominators.emplace_back(vPrime, 0);
        }
    }
    
    void is_center_dominated(const Vertex& u, const Vertex& v, const Vertex_set_tuple& e_nbd_tuple,
                            std::vector<std::pair<Vertex, int>>& dominators,
                            Filtration_value max_filtration) const {
        const auto& commonIncoming = std::get<0>(e_nbd_tuple);
        const auto& commonOutgoing = std::get<1>(e_nbd_tuple);
        const auto& outgoingIncoming = std::get<2>(e_nbd_tuple);
        
        const auto incomingV = graph.get_incoming_neighbors_set(v, max_filtration);
        const auto outgoingU = graph.get_outgoing_neighbors_set(u, max_filtration);
        const auto outgoingIncomingUV = SetOps::setIntersection(outgoingU, incomingV);
        
        for (size_t vPrime = 0; vPrime < MAX_VERTICES; ++vPrime) {
            if (!outgoingIncoming[vPrime]) continue;
            
            const auto incomingVPrime = graph.get_incoming_neighbors_set(vPrime, max_filtration);
            const auto outgoingVPrime = graph.get_outgoing_neighbors_set(vPrime, max_filtration);
            
            if (!SetOps::isSubset(commonOutgoing, outgoingVPrime)) continue;
            if (!SetOps::isSubset(commonIncoming, incomingVPrime)) continue;
            
            auto unionVPrime = incomingVPrime;
            unionVPrime |= outgoingVPrime;
            unionVPrime.set(vPrime);
            if (!SetOps::isSubset(outgoingIncoming, unionVPrime)) continue;
            
            const auto UVp_V = SetOps::setIntersection(outgoingIncomingUV, outgoingVPrime);
            const auto U_VpV = SetOps::setIntersection(outgoingIncomingUV, incomingVPrime);
            
            if (!split_edge_condition(UVp_V, U_VpV)) continue;
            
            dominators.emplace_back(vPrime, 2);
        }
    }
    
    void is_out_dominated(const Vertex& u, const Vertex& v, const Vertex_set_tuple& e_nbd_tuple,
                         std::vector<std::pair<Vertex, int>>& dominators,
                         Filtration_value max_filtration) const {
        const auto& commonIncoming = std::get<0>(e_nbd_tuple);
        const auto& commonOutgoing = std::get<1>(e_nbd_tuple);
        const auto& outgoingIncoming = std::get<2>(e_nbd_tuple);
        
        const auto outgoingU = graph.get_outgoing_neighbors_set(u, max_filtration);
        const auto outgoingV = graph.get_outgoing_neighbors_set(v, max_filtration);
        const auto outgoingUV = SetOps::setIntersection(outgoingU, outgoingV);
        
        for (size_t vPrime = 0; vPrime < MAX_VERTICES; ++vPrime) {
            if (!commonOutgoing[vPrime]) continue;
            
            const auto incomingVPrime = graph.get_incoming_neighbors_set(vPrime, max_filtration);
            const auto outgoingVPrime = graph.get_outgoing_neighbors_set(vPrime, max_filtration);
            
            const auto unionInNeighbors = SetOps::setUnion(commonIncoming, outgoingIncoming);
            if (!SetOps::isSubset(unionInNeighbors, incomingVPrime)) continue;
            
            auto unionVPrime = incomingVPrime;
            unionVPrime |= outgoingVPrime;
            unionVPrime.set(vPrime);
            if (!SetOps::isSubset(commonOutgoing, unionVPrime)) continue;
            
            const auto outgoingIncomingUVVp = SetOps::setIntersection(outgoingUV, incomingVPrime);
            const auto commonOutgoingUVVp = SetOps::setIntersection(outgoingVPrime, outgoingUV);
            
            if (!split_edge_condition(commonOutgoingUVVp, outgoingIncomingUVVp)) continue;
            
            dominators.emplace_back(vPrime, 1);
        }
    }
};
