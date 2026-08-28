import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ==========================================
# 1. CONFIGURATION
# ==========================================
WIDTH  = 1
HEIGHT = 5
DEPTH  = 5

ENABLE_VISUALS = True  # Set to False to run in terminal-only mode (blazing fast)

TOTAL_SQUARES = 2 * (WIDTH * HEIGHT + HEIGHT * DEPTH + DEPTH * WIDTH)

# ==========================================
# 2. DYNAMIC 3D TOPOLOGY GENERATOR
# ==========================================
def generate_cuboid(W, H, D):
    nodes = []
    node_id = 0
    
    def add_face_nodes(w_count, h_count, center_func, normal, u, v):
        nonlocal node_id
        for i in range(w_count):
            for j in range(h_count):
                nodes.append({
                    'id': node_id, 'C': center_func(i, j),
                    'N': normal, 'u': u, 'v': v
                })
                node_id += 1

    add_face_nodes(W, H, lambda x, y: (x+0.5, y+0.5, 0), (0,0,-1), (1,0,0), (0,-1,0))
    add_face_nodes(W, H, lambda x, y: (x+0.5, y+0.5, D), (0,0,1),  (1,0,0), (0,1,0))
    add_face_nodes(W, D, lambda x, z: (x+0.5, 0, z+0.5), (0,-1,0), (1,0,0), (0,0,1))
    add_face_nodes(W, D, lambda x, z: (x+0.5, H, z+0.5), (0,1,0), (-1,0,0), (0,0,1))
    add_face_nodes(H, D, lambda y, z: (0, y+0.5, z+0.5), (-1,0,0), (0,-1,0), (0,0,1))
    add_face_nodes(H, D, lambda y, z: (W, y+0.5, z+0.5), (1,0,0),  (0,1,0), (0,0,1))

    nodes_map = {(n['C'][0], n['C'][1], n['C'][2], n['N'][0], n['N'][1], n['N'][2]): n for n in nodes}
    graph = {}
    faces_3d = {}
    
    def v_add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
    def v_sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def v_mul(s, v): return (s*v[0], s*v[1], s*v[2])
    def v_neg(v):    return (-v[0], -v[1], -v[2])

    for n in nodes:
        nid = n['id']
        C, N, u, v = n['C'], n['N'], n['u'], n['v']
        
        c0 = v_sub(v_sub(C, v_mul(0.5, u)), v_mul(0.5, v))
        c1 = v_sub(v_add(C, v_mul(0.5, u)), v_mul(0.5, v))
        c2 = v_add(v_add(C, v_mul(0.5, u)), v_mul(0.5, v))
        c3 = v_add(v_sub(C, v_mul(0.5, u)), v_mul(0.5, v))
        faces_3d[nid] = [c0, c1, c2, c3]
        
        graph[nid] = {}
        dirs = [v, u, v_neg(v), v_neg(u)]
        
        for local_dir_idx, vec in enumerate(dirs):
            c_test = v_add(C, vec)
            pos_norm = (c_test[0], c_test[1], c_test[2], N[0], N[1], N[2])
            
            if pos_norm in nodes_map:
                neighbor = nodes_map[pos_norm]
                in_vec = v_neg(vec)
            else:
                c_test = v_sub(v_add(C, v_mul(0.5, vec)), v_mul(0.5, N))
                pos_norm = (c_test[0], c_test[1], c_test[2], vec[0], vec[1], vec[2])
                neighbor = nodes_map[pos_norm]
                in_vec = N
                
            nu, nv = neighbor['u'], neighbor['v']
            if in_vec == nv: back_dir = 0
            elif in_vec == nu: back_dir = 1
            elif in_vec == v_neg(nv): back_dir = 2
            elif in_vec == v_neg(nu): back_dir = 3
            
            graph[nid][local_dir_idx] = (neighbor['id'], back_dir)
            
    return graph, faces_3d

graph, faces_3d = generate_cuboid(WIDTH, HEIGHT, DEPTH)
DIRS = {0: (0, 1), 1: (1, 0), 2: (0, -1), 3: (-1, 0)}

# ==========================================
# 3. DEDUPLICATION
# ==========================================
def get_canonical_net(coords):
    forms = []
    transformations = [
        lambda x, y: (x, y), lambda x, y: (-y, x), lambda x, y: (-x, -y), lambda x, y: (y, -x),
        lambda x, y: (-x, y), lambda x, y: (x, -y), lambda x, y: (y, x), lambda x, y: (-y, -x)
    ]
    for transform in transformations:
        transformed = [transform(x, y) for x, y in coords]
        min_x = min(x for x, y in transformed)
        min_y = min(y for x, y in transformed)
        shifted = tuple(sorted((x - min_x, y - min_y) for x, y in transformed))
        forms.append(shifted)
    return min(forms)

# ==========================================
# 4. OPTIMIZED REDELMEIER DFS 
# ==========================================
def dfs_unfold(visited, net_coords, node_states, available_edges, tree_edges):
    if len(visited) == TOTAL_SQUARES:
        yield ("SUCCESS", None, visited, net_coords, tree_edges)
        return

    for i in range(len(available_edges)):
        edge = available_edges[i]
        p_node, p_dir, c_node, c_in = edge
        
        if c_node in visited:
            continue 
            
        p_pos, p_offset = node_states[p_node]
        
        abs_dir = (p_dir + p_offset) % 4
        dx, dy = DIRS[abs_dir]
        c_pos = (p_pos[0] + dx, p_pos[1] + dy)
        incoming_abs_dir = (abs_dir + 2) % 4
        c_offset = (incoming_abs_dir - c_in) % 4
        
        if c_pos in net_coords:
            yield ("COLLISION", c_pos, visited, net_coords, tree_edges)
            continue 
            
        visited.add(c_node)
        net_coords[c_pos] = c_node
        node_states[c_node] = (c_pos, c_offset)
        tree_edges.add((p_node, c_node))
        tree_edges.add((c_node, p_node))
        
        yield ("ADDED", c_pos, visited, net_coords, tree_edges)
        
        new_edges = []
        for d in range(4):
            if d != c_in:
                n, n_in = graph[c_node][d]
                if n not in visited:
                    new_edges.append((c_node, d, n, n_in))
                    
        next_available = available_edges[i+1:] + new_edges
        
        yield from dfs_unfold(visited, net_coords, node_states, next_available, tree_edges)
        
        visited.remove(c_node)
        del net_coords[c_pos]
        del node_states[c_node]
        tree_edges.remove((p_node, c_node))
        tree_edges.remove((c_node, p_node))
        yield ("BACKTRACK", c_pos, visited, net_coords, tree_edges)

# ==========================================
# 5. VISUALIZATION ENGINE 
# ==========================================
class SimState:
    fast_forward = False

def run_simulation():
    output_dir = "net_generator_optimised_proof_of_concept_images"
    
    if ENABLE_VISUALS:
        os.makedirs(output_dir, exist_ok=True)
        plt.ion()
        fig = plt.figure(figsize=(12, 6))
        ax3d = fig.add_subplot(121, projection='3d')
        ax2d = fig.add_subplot(122)
        fig.canvas.manager.set_window_title("Press SPACE to toggle Fast-Forward")
        
        sim_state = SimState()
        def on_key(event):
            if event.key == ' ':
                sim_state.fast_forward = not sim_state.fast_forward
                print(f"\n>>> FAST FORWARD: {'ON' if sim_state.fast_forward else 'OFF'} <<<\n")
        fig.canvas.mpl_connect('key_press_event', on_key)
    else:
        # If visuals are disabled, naturally fast-forward through everything
        sim_state = SimState()
        sim_state.fast_forward = True
    
    visited = {0}
    net_coords = {(0, 0): 0}
    node_states = {0: ((0, 0), 0)}
    initial_edges = [(0, d, graph[0][d][0], graph[0][d][1]) for d in range(4)]
    tree_edges = set()
    
    dfs_gen = dfs_unfold(visited, net_coords, node_states, initial_edges, tree_edges)
    
    unique_nets_found = set()
    stats = {"attempts": 0, "collisions": 0, "duplicates": 0, "unique": 0}
    
    max_dim = max(WIDTH, HEIGHT, DEPTH)
    pad = TOTAL_SQUARES // 2 + 1
    
    print(f"Starting Optimized Search on {WIDTH}x{HEIGHT}x{DEPTH} Box...")
    if ENABLE_VISUALS:
        print("Keep the plot window focused and PRESS SPACEBAR to fast-forward.\n")
    else:
        print("Visuals DISABLED. Running in terminal-only mode...\n")
    
    for status, pos, current_visited, current_net_coords, current_tree_edges in dfs_gen:
        stats["attempts"] += 1
        
        # Keep UI responsive during intense fast-forwarding
        if ENABLE_VISUALS and sim_state.fast_forward and stats["attempts"] % 25 == 0:
            plt.pause(0.001)

        # Logic Logging
        if status == "COLLISION":
            stats["collisions"] += 1
            if not sim_state.fast_forward: print(f"[{stats['attempts']:03d}] INVALID: Overlap detected")
        elif status == "SUCCESS":
            coords = list(current_net_coords.keys())
            canonical_id = get_canonical_net(coords)
            
            if canonical_id in unique_nets_found:
                stats["duplicates"] += 1
                if not sim_state.fast_forward: print(f"[{stats['attempts']:03d}] DUPLICATE: Tree maps to existing shape.")
            else:
                stats["unique"] += 1
                unique_nets_found.add(canonical_id)
                print(f"[{stats['attempts']:03d}] UNIQUE: Found net #{stats['unique']}!")

        # Skip rendering if visuals are disabled entirely
        if not ENABLE_VISUALS:
            continue

        # Skip drawing if we are fast-forwarding and it's not a brand new unique net
        if sim_state.fast_forward and (status != "SUCCESS" or canonical_id in unique_nets_found and len(unique_nets_found) > 0):
            continue

        # --- RENDER FRAME ---
        ax3d.clear()
        ax2d.clear()
        
        # 3D Render
        ax3d.set_title(f"3D Box ({WIDTH}x{HEIGHT}x{DEPTH})\nUnique Nets: {stats['unique']}")
        ax3d.set_xlim([0, max_dim]); ax3d.set_ylim([0, max_dim]); ax3d.set_zlim([0, max_dim])
        ax3d.view_init(elev=20, azim=45)
        
        edge_corners = {0: (3, 2), 1: (1, 2), 2: (0, 1), 3: (0, 3)}
        
        for face_id, corners in faces_3d.items():
            if face_id in current_visited:
                poly = Poly3DCollection([corners], alpha=0.8, facecolor='blue', edgecolor='none')
                ax3d.add_collection3d(poly)
                for local_dir in range(4):
                    neighbor, _ = graph[face_id][local_dir]
                    is_hinge = (face_id, neighbor) in current_tree_edges
                    ecolor = 'black' if is_hinge else 'red'
                    lw = 3 if is_hinge else 2
                    
                    idx1, idx2 = edge_corners[local_dir]
                    p1, p2 = corners[idx1], corners[idx2]
                    ax3d.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color=ecolor, linewidth=lw)
            else:
                poly = Poly3DCollection([corners], alpha=0.2, facecolor='whitesmoke', edgecolor='lightgray')
                ax3d.add_collection3d(poly)
                
        # 2D Render
        ax2d.set_title(f"Status: {status} (SPACE: {'FFWD' if sim_state.fast_forward else 'NORM'})\nUnique Found: {stats['unique']}")
        ax2d.set_xlim([-pad, pad]); ax2d.set_ylim([-pad, pad])
        ax2d.set_aspect('equal')
        ax2d.grid(True, linestyle=':', alpha=0.6)
        
        for (nx, ny), node_id in current_net_coords.items():
            face_color = 'lightgreen' if status == "SUCCESS" else 'lightgray'
            rect = Rectangle((nx - 0.5, ny - 0.5), 1, 1, facecolor=face_color, edgecolor='black')
            ax2d.add_patch(rect)
            ax2d.text(nx, ny, str(node_id), ha='center', va='center')
        
        if status == "COLLISION":
            rect = Rectangle((pos[0] - 0.5, pos[1] - 0.5), 1, 1, facecolor='red', edgecolor='black', alpha=0.5)
            ax2d.add_patch(rect)
            
        # UI Timing & Saves
        if status == "SUCCESS" and canonical_id in unique_nets_found:
            filename = os.path.join(output_dir, f"net_{stats['unique']:03d}.png")
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            
        if sim_state.fast_forward:
            plt.pause(0.001)
        else:
            plt.pause(0.45)

    print("\n" + "="*50)
    print(f"SEARCH COMPLETE for {WIDTH}x{HEIGHT}x{DEPTH}")
    print(f"Total Operations: {stats['attempts']}")
    print(f"Total Collisions: {stats['collisions']}")
    print(f"Total Duplicates: {stats['duplicates']}")
    print(f"Unique Nets:      {stats['unique']}")
    print("="*50)

    if ENABLE_VISUALS:
        plt.ioff()
        plt.show()

if __name__ == "__main__":
    run_simulation()