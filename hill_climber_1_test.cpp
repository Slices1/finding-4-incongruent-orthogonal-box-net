#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;

// ==========================================
// 1. DATA STRUCTURES & MATH
// ==========================================
struct Point { int x, y; };

struct Vec3 {
    int x, y, z;
    bool operator<(const Vec3& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }
    bool operator==(const Vec3& o) const { return x == o.x && y == o.y && z == o.z; }
};

Vec3 v_add(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 v_sub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 v_mul(int s, Vec3 v) { return {s * v.x, s * v.y, s * v.z}; }
Vec3 v_neg(Vec3 v) { return {-v.x, -v.y, -v.z}; }

struct Node { int id; Vec3 C, N, u, v; };
struct Edge { int to, out_dir, in_dir; };
struct EdgeData { int p_node, p_dir, c_node, c_in; };

const Point DIRS[4] = {{0, 1}, {1, 0}, {0, -1}, {-1, 0}};

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

// ==========================================
// 2. TOPOLOGY GENERATOR
// ==========================================
struct BoxGraph {
    int W, H, D, total_squares;
    vector<vector<Edge>> adj;
    vector<pair<int, int>> all_edges;

    BoxGraph(int w, int h, int d) : W(w), H(h), D(d) {
        total_squares = 2 * (W * H + H * D + D * W);
        adj.resize(total_squares);
        build();
    }

    void build() {
        vector<Node> nodes;
        int id = 0;
        
        auto add_face = [&](int wc, int hc, auto center_func, Vec3 N, Vec3 u, Vec3 v) {
            for (int i = 0; i < wc; ++i) {
                for (int j = 0; j < hc; ++j) {
                    nodes.push_back({id++, center_func(i, j), N, u, v});
                }
            }
        };

        add_face(W, H, [](int x, int y) -> Vec3 { return {x*2+1, y*2+1, 0}; }, {0,0,-2}, {2,0,0}, {0,-2,0});
        add_face(W, H, [this](int x, int y) -> Vec3 { return {x*2+1, y*2+1, D*2}; }, {0,0,2}, {2,0,0}, {0,2,0});
        add_face(W, D, [](int x, int z) -> Vec3 { return {x*2+1, 0, z*2+1}; }, {0,-2,0}, {2,0,0}, {0,0,2});
        add_face(W, D, [this](int x, int z) -> Vec3 { return {x*2+1, H*2, z*2+1}; }, {0,2,0}, {-2,0,0}, {0,0,2});
        add_face(H, D, [](int y, int z) -> Vec3 { return {0, y*2+1, z*2+1}; }, {-2,0,0}, {0,-2,0}, {0,0,2});
        add_face(H, D, [this](int y, int z) -> Vec3 { return {W*2, y*2+1, z*2+1}; }, {2,0,0}, {0,2,0}, {0,0,2});

        map<pair<Vec3, Vec3>, int> nodes_map;
        for (auto& n : nodes) nodes_map[{n.C, n.N}] = n.id;

        for (auto& n : nodes) {
            Vec3 dirs[4] = {n.v, n.u, v_neg(n.v), v_neg(n.u)};
            for (int d = 0; d < 4; ++d) {
                Vec3 vec = dirs[d];
                Vec3 c_test = v_add(n.C, vec);
                Vec3 in_vec;
                int neighbor_id;

                if (nodes_map.count({c_test, n.N})) {
                    neighbor_id = nodes_map[{c_test, n.N}];
                    in_vec = v_neg(vec);
                } else {
                    Vec3 half_vec = {vec.x/2, vec.y/2, vec.z/2};
                    Vec3 half_N = {n.N.x/2, n.N.y/2, n.N.z/2};
                    Vec3 corner_test = v_sub(v_add(n.C, half_vec), half_N);
                    
                    if (!nodes_map.count({corner_test, vec})) {
                        cout << "\nCRITICAL ERROR: Mathematical hole in topology at Box " << W << "x" << H << "x" << D << endl;
                        exit(1);
                    }
                    neighbor_id = nodes_map[{corner_test, vec}];
                    in_vec = n.N;
                }

                auto& nb = nodes[neighbor_id];
                int back_dir = -1;
                if (in_vec == nb.v) back_dir = 0;
                else if (in_vec == nb.u) back_dir = 1;
                else if (in_vec == v_neg(nb.v)) back_dir = 2;
                else if (in_vec == v_neg(nb.u)) back_dir = 3;

                adj[n.id].push_back({neighbor_id, d, back_dir});
                if (n.id < neighbor_id) all_edges.push_back({n.id, neighbor_id});
            }
        }
    }
};

// ==========================================
// 3. POLYOMINO MATH & OVERLAPS
// ==========================================
void transform_net(const vector<Point>& pts, int form, vector<Point>& res) {
    res.clear();
    for (auto& p : pts) {
        if (form == 0) res.push_back({p.x, p.y});
        else if (form == 1) res.push_back({-p.y, p.x});
        else if (form == 2) res.push_back({-p.x, -p.y});
        else if (form == 3) res.push_back({p.y, -p.x});
        else if (form == 4) res.push_back({-p.x, p.y});
        else if (form == 5) res.push_back({p.x, -p.y});
        else if (form == 6) res.push_back({p.y, p.x});
        else if (form == 7) res.push_back({-p.y, -p.x});
    }
}

int overlap_counts[600][600] = {0}; 

int get_max_overlap(const vector<Point>& A, const vector<Point>& B, vector<Point>& transB, vector<Point>& touched) {
    int max_overlap = 0;
    for (int form = 0; form < 8; ++form) {
        transform_net(B, form, transB);
        touched.clear();
        
        for (auto& pa : A) {
            for (auto& pb : transB) {
                int dx = pa.x - pb.x + 300;
                int dy = pa.y - pb.y + 300;
                
                if (dx < 0 || dx >= 600 || dy < 0 || dy >= 600) continue; 
                
                if (overlap_counts[dx][dy] == 0) touched.push_back({dx, dy});
                overlap_counts[dx][dy]++;
                if (overlap_counts[dx][dy] > max_overlap) {
                    max_overlap = overlap_counts[dx][dy];
                }
            }
        }
        for(auto& p : touched) overlap_counts[p.x][p.y] = 0;
    }
    return max_overlap;
}

// ==========================================
// 4. TREE & NET MANAGEMENT
// ==========================================
struct NetState {
    BoxGraph* box;
    vector<pair<int, int>> tree_edges;
    vector<Point> coords;

    vector<pair<int, int>> available;
    vector<vector<int>> t_adj;
    vector<bool> vis;
    vector<int> parent;
    vector<int> q;
    vector<pair<int, int>> cycle;
    vector<vector<Edge>> edge_adj;
    vector<bool> visited_coords;
    vector<pair<int, int>> q_coords;
    vector<Point> pos_map;

    NetState(BoxGraph* b) : box(b) {
        t_adj.resize(box->total_squares);
        vis.resize(box->total_squares);
        parent.resize(box->total_squares);
        edge_adj.resize(box->total_squares);
        visited_coords.resize(box->total_squares);
        pos_map.resize(box->total_squares);
    }

    bool build_coords(int start_node) {
        coords.clear();
        for(auto& v : edge_adj) v.clear();
        for (auto& e : tree_edges) {
            for(auto& ge : box->adj[e.first]) if(ge.to == e.second) edge_adj[e.first].push_back(ge);
            for(auto& ge : box->adj[e.second]) if(ge.to == e.first) edge_adj[e.second].push_back(ge);
        }

        fill(visited_coords.begin(), visited_coords.end(), false);
        q_coords.clear();
        
        q_coords.push_back({start_node, 0});
        visited_coords[start_node] = true;
        pos_map[start_node] = {0, 0};
        coords.push_back({0, 0});

        int head = 0;
        while(head < q_coords.size()) {
            int curr = q_coords[head].first;
            int c_offset = q_coords[head].second;
            Point c_pos = pos_map[curr];
            head++;

            for (auto& edge : edge_adj[curr]) {
                if (!visited_coords[edge.to]) {
                    int abs_dir = (edge.out_dir + c_offset) % 4;
                    Point n_pos = {c_pos.x + DIRS[abs_dir].x, c_pos.y + DIRS[abs_dir].y};
                    
                    for (auto& c : coords) if (c.x == n_pos.x && c.y == n_pos.y) return false;
                    
                    coords.push_back(n_pos);
                    pos_map[edge.to] = n_pos;
                    visited_coords[edge.to] = true;
                    
                    int in_abs_dir = (abs_dir + 2) % 4;
                    int n_offset = (in_abs_dir - edge.in_dir + 4) % 4;
                    q_coords.push_back({edge.to, n_offset});
                }
            }
        }
        return coords.size() == box->total_squares;
    }

    bool build_initial_net(vector<bool>& visited_dfs, vector<Point>& occ, 
                           vector<pair<Point, int>>& states, vector<EdgeData>& frontier) {
        if (occ.size() == box->total_squares) return true;

        for (int i = 0; i < frontier.size(); ++i) {
            auto ed = frontier[i];
            if (visited_dfs[ed.c_node]) continue;

            auto p_state = states[ed.p_node];
            int abs_dir = (ed.p_dir + p_state.second) % 4;
            Point c_pos = {p_state.first.x + DIRS[abs_dir].x, p_state.first.y + DIRS[abs_dir].y};

            bool collision = false;
            for(auto& p : occ) if(p.x == c_pos.x && p.y == c_pos.y) { collision = true; break; }
            if (collision) continue;

            int in_abs = (abs_dir + 2) % 4;
            int c_offset = (in_abs - ed.c_in + 4) % 4;

            visited_dfs[ed.c_node] = true;
            occ.push_back(c_pos);
            states[ed.c_node] = {c_pos, c_offset};
            tree_edges.push_back({ed.p_node, ed.c_node});

            vector<EdgeData> next_frontier;
            for (int j = 0; j < frontier.size(); ++j) {
                if (j != i && !visited_dfs[frontier[j].c_node]) next_frontier.push_back(frontier[j]);
            }
            for (auto& nbr : box->adj[ed.c_node]) {
                if (nbr.out_dir != ed.c_in && !visited_dfs[nbr.to]) {
                    next_frontier.push_back({ed.c_node, nbr.out_dir, nbr.to, nbr.in_dir});
                }
            }

            if (build_initial_net(visited_dfs, occ, states, next_frontier)) return true;

            tree_edges.pop_back();
            occ.pop_back();
            visited_dfs[ed.c_node] = false;
        }
        return false;
    }

    void generate_seed() {
        tree_edges.clear();
        vector<bool> visited_dfs(box->total_squares, false);
        vector<Point> occ;
        vector<pair<Point, int>> states(box->total_squares);
        vector<EdgeData> frontier;

        visited_dfs[0] = true;
        occ.push_back({0,0});
        states[0] = {{0,0}, 0};

        for (auto& nbr : box->adj[0]) {
            frontier.push_back({0, nbr.out_dir, nbr.to, nbr.in_dir});
        }

        if (!build_initial_net(visited_dfs, occ, states, frontier)) {
            cout << "FATAL ERROR: Topology connected, but failed to unfold." << endl;
            exit(1);
        }
        build_coords(0); 

        // Apply random mutations so the plateau breaker generates a scrambled shape!
        for(int k=0; k<100; ++k) { mutate(); }
    }

    bool mutate() {
        available.clear();
        for (auto& e : box->all_edges) {
            bool found = false;
            for (auto& te : tree_edges) {
                if ((te.first == e.first && te.second == e.second) ||
                    (te.first == e.second && te.second == e.first)) {
                    found = true; break;
                }
            }
            if (!found) available.push_back(e);
        }

        auto added_edge = available[rng() % available.size()];
        
        for (auto& v : t_adj) v.clear();
        for (auto& e : tree_edges) {
            t_adj[e.first].push_back(e.second);
            t_adj[e.second].push_back(e.first);
        }

        fill(vis.begin(), vis.end(), false);
        q.clear();
        q.push_back(added_edge.first);
        vis[added_edge.first] = true;

        int head = 0;
        while(head < q.size()) {
            int u = q[head++];
            if(u == added_edge.second) break;
            for(int v : t_adj[u]) {
                if(!vis[v]) { vis[v] = true; parent[v] = u; q.push_back(v); }
            }
        }

        cycle.clear();
        int curr = added_edge.second;
        while(curr != added_edge.first) {
            cycle.push_back({parent[curr], curr});
            curr = parent[curr];
        }

        auto removed_edge = cycle[rng() % cycle.size()];
        
        for (auto& e : tree_edges) {
            if ((e.first == removed_edge.first && e.second == removed_edge.second) ||
                (e.first == removed_edge.second && e.second == removed_edge.first)) {
                e = added_edge;
                break;
            }
        }
        
        if (build_coords(0)) return true; 
        
        // Internal Rollback
        for (auto& e : tree_edges) {
            if ((e.first == added_edge.first && e.second == added_edge.second) ||
                (e.first == added_edge.second && e.second == added_edge.first)) {
                e = removed_edge;
                break;
            }
        }
        return false;
    }
};

// ==========================================
// 5. HILL CLIMBING ENGINE
// ==========================================

// UPDATED for 2 boxes: only 1 pairwise comparison instead of 6 pairs.
int evaluate_fitness(vector<NetState>& nets, vector<Point>& transB, vector<Point>& touched) {
    int total_score = 0;
    for(int i = 0; i < 2; ++i) {
        for(int j = i+1; j < 2; ++j) {
            total_score += get_max_overlap(nets[i].coords, nets[j].coords, transB, touched);
        }
    }
    return total_score;
}

void print_success(vector<NetState>& nets) {
    cout << "\n==================================================" << endl;
    cout << " SUCCESS! FOUND A UNIVERSALLY SHARED 22-SQUARE NET!" << endl;
    cout << "==================================================" << endl;

    for (int i = 0; i < 2; ++i) {
        cout << "\n--------------------------------------------------" << endl;
        cout << "BOX " << i+1 << " (" << nets[i].box->W << "x" << nets[i].box->H << "x" << nets[i].box->D << ")" << endl;
        cout << "--------------------------------------------------" << endl;

        int min_x = 1000, min_y = 1000, max_x = -1000, max_y = -1000;
        for (auto p : nets[i].coords) {
            min_x = min(min_x, p.x); min_y = min(min_y, p.y);
            max_x = max(max_x, p.x); max_y = max(max_y, p.y);
        }

        cout << "2D ASCII Pattern (Numbers are 3D Node IDs. Adjacent squares imply a fold):\n" << endl;
        for (int y = max_y; y >= min_y; --y) { // Start from max_y down so it prints correctly visually
            for (int x = min_x; x <= max_x; ++x) {
                int node_id = -1;
                for(int n = 0; n < nets[i].box->total_squares; ++n) {
                    if (nets[i].pos_map[n].x == x && nets[i].pos_map[n].y == y) {
                        node_id = n; 
                        break;
                    }
                }

                if (node_id != -1) {
                    if (node_id < 10) cout << "[" << node_id << " ]";
                    else cout << "[" << node_id << "]";
                } else {
                    cout << "    "; 
                }
            }
            cout << endl;
        }
        
        cout << "\nHinges (3D Edges kept intact for folding):" << endl;
        for (auto e : nets[i].tree_edges) {
            cout << e.first << "-" << e.second << "  ";
        }
        cout << endl;
    }
}

int main() {
    cout << "Initializing 3D Box Graphs (Testing 2-Box System)..." << endl;
    BoxGraph b1(1, 1, 5);
    BoxGraph b2(1, 2, 3);

    vector<NetState> nets = {{&b1}, {&b2}};
    
    cout << "Generating deterministic valid starting nets (22 squares)..." << endl;
    for (int i = 0; i < 2; ++i) nets[i].generate_seed();

    vector<Point> transB; transB.reserve(100);
    vector<Point> touched; touched.reserve(5000);

    int current_score = evaluate_fitness(nets, transB, touched);
    
    // Max possible score for 1 pair of 22-square boxes is exactly 22.
    int max_possible_score = 22; 
    
    cout << "Starting Fitness: " << current_score << " / " << max_possible_score << endl;
    cout << "Beginning Stochastic Hill Climb...\n" << endl;

    long long iterations = 0;
    long long successful_mutations = 0;
    long long iterations_since_improvement = 0; 
    
    vector<pair<int, int>> backup_tree;
    vector<Point> backup_coords;

    while (current_score < max_possible_score) {
        iterations++;
        int target = rng() % 2; // Only pick between 0 and 1
        
        backup_tree = nets[target].tree_edges;
        backup_coords = nets[target].coords;
        
        int attempt = 0;
        while (!nets[target].mutate()) {
            attempt++;
            if (attempt > 200000) {
                cout << "[Warning] Box " << target << " struggling to find valid mutation..." << endl;
                attempt = 0; 
            }
        }

        int new_score = evaluate_fitness(nets, transB, touched);

        if (new_score >= current_score) {
            if (new_score > current_score) {
                iterations_since_improvement = 0; 
            } else {
                iterations_since_improvement++;
            }

            current_score = new_score;
            successful_mutations++;
            
            if (new_score > current_score || successful_mutations % 1000 == 0) {
                cout << "[Iter: " << iterations << "] New Fitness: " << current_score 
                     << " / " << max_possible_score << endl;
            }
        } else {
            nets[target].tree_edges = backup_tree;
            nets[target].coords = backup_coords;
            iterations_since_improvement++;
        }

        if (iterations_since_improvement > 300000) {
            cout << "\n*** PLATEAU DETECTED *** Nuking Box " << target << " to break the deadlock!\n" << endl;
            nets[target].generate_seed(); 
            current_score = evaluate_fitness(nets, transB, touched);
            iterations_since_improvement = 0;
        }
    }

    print_success(nets);
    return 0;
}