import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# 1. DEFINE THE 3D GRAPH FOR A 1x1x1 CUBE
graph = {
    0: {0: (1, 2), 1: (2, 2), 2: (3, 2), 3: (4, 2)},  # Bottom
    1: {0: (5, 2), 1: (2, 3), 2: (0, 0), 3: (4, 1)},  # Front
    2: {0: (5, 1), 1: (3, 3), 2: (0, 1), 3: (1, 1)},  # Right
    3: {0: (5, 0), 1: (4, 3), 2: (0, 2), 3: (2, 1)},  # Back
    4: {0: (5, 3), 1: (1, 3), 2: (0, 3), 3: (3, 1)},  # Left
    5: {0: (3, 0), 1: (2, 0), 2: (1, 0), 3: (4, 0)},  # Top
}

faces_3d = {
    0: [(0,0,0), (1,0,0), (1,1,0), (0,1,0)],  # Bottom
    1: [(0,0,0), (1,0,0), (1,0,1), (0,0,1)],  # Front
    2: [(1,0,0), (1,1,0), (1,1,1), (1,0,1)],  # Right
    3: [(1,1,0), (0,1,0), (0,1,1), (1,1,1)],  # Back
    4: [(0,1,0), (0,0,0), (0,0,1), (0,1,1)],  # Left
    5: [(0,0,1), (1,0,1), (1,1,1), (0,1,1)],  # Top
}

DIRS = {0: (0, 1), 1: (1, 0), 2: (0, -1), 3: (-1, 0)}


# 2. DEDUPLICATION (Canonicalization)
def get_canonical_net(coords):
    """
    Takes a list of (x,y) coordinates and returns a unique, standardized tuple 
    that represents the shape regardless of its position, rotation, or reflection.
    """
    forms = []
    # The 8 symmetries of a 2D grid (rotations and reflections)
    transformations = [
        lambda x, y: (x, y),          # Original
        lambda x, y: (-y, x),         # Rotate 90
        lambda x, y: (-x, -y),        # Rotate 180
        lambda x, y: (y, -x),         # Rotate 270
        lambda x, y: (-x, y),         # Reflect across Y-axis
        lambda x, y: (x, -y),         # Reflect across X-axis
        lambda x, y: (y, x),          # Reflect across Y=X
        lambda x, y: (-y, -x)         # Reflect across Y=-X
    ]
    
    for transform in transformations:
        # Apply the transformation
        transformed = [transform(x, y) for x, y in coords]
        
        # Shift to origin (bottom-left at 0,0)
        min_x = min(x for x, y in transformed)
        min_y = min(y for x, y in transformed)
        shifted = tuple(sorted((x - min_x, y - min_y) for x, y in transformed))
        
        forms.append(shifted)
        
    # Return the lexicographically smallest representation as the canonical ID
    return min(forms)


# 3. THE BRANCHING DFS GENERATOR
def dfs_unfold(visited, net_coords, node_states):
    if len(visited) == 6:
        yield ("SUCCESS", None, visited, net_coords)
        return

    # Find ALL unvisited neighbors of ALL currently visited squares
    possible_moves = []
    for v_node in visited:
        v_pos, v_offset = node_states[v_node]
        for local_dir in range(4):
            neighbor, back_dir = graph[v_node][local_dir]
            if neighbor not in visited:
                possible_moves.append((v_node, v_pos, v_offset, local_dir, neighbor, back_dir))
                
    for move in possible_moves:
        v_node, v_pos, v_offset, local_dir, neighbor, back_dir = move
        
        abs_dir = (local_dir + v_offset) % 4
        dx, dy = DIRS[abs_dir]
        next_pos = (v_pos[0] + dx, v_pos[1] + dy)
        
        incoming_abs_dir = (abs_dir + 2) % 4
        next_offset = (incoming_abs_dir - back_dir) % 4
        
        if next_pos in net_coords:
            continue # Prune this path (Collision)
            
        visited.add(neighbor)
        net_coords[next_pos] = neighbor
        node_states[neighbor] = (next_pos, next_offset)
        
        # We are skipping yielding "ADDED" here to let the script run lightning fast
        yield from dfs_unfold(visited, net_coords, node_states)
        
        visited.remove(neighbor)
        del net_coords[next_pos]
        del node_states[neighbor]


# 4. VISUALIZATION AND SAVING ENGINE
def run_simulation():
    # Setup Output Directory
    output_dir = "net_generator_proof_of_concept_images"
    os.makedirs(output_dir, exist_ok=True)
    
    plt.ion()
    fig = plt.figure(figsize=(12, 6))
    ax3d = fig.add_subplot(121, projection='3d')
    ax2d = fig.add_subplot(122)
    
    # Initialize the DFS State
    visited = {0}
    net_coords = {(0, 0): 0}
    node_states = {0: ((0, 0), 0)}
    dfs_gen = dfs_unfold(visited, net_coords, node_states)
    
    unique_nets_found = set()
    
    print(f"Starting search... Saving images to ./{output_dir}/")
    
    for status, pos, current_visited, current_net_coords in dfs_gen:
        if status == "SUCCESS":
            # Extract coordinates and generate the unique fingerprint
            coords = list(current_net_coords.keys())
            canonical_id = get_canonical_net(coords)
            
            # If we haven't seen this exact shape before:
            if canonical_id not in unique_nets_found:
                unique_nets_found.add(canonical_id)
                net_number = len(unique_nets_found)
                
                # --- UPDATE THE GUI ---
                ax3d.clear()
                ax2d.clear()
                
                # Draw 3D Cube
                ax3d.set_title(f"3D Cube\nUnique Nets Found: {net_number}/11")
                ax3d.set_xlim([0, 1]); ax3d.set_ylim([0, 1]); ax3d.set_zlim([0, 1])
                ax3d.view_init(elev=20, azim=45)
                for face_id, corners in faces_3d.items():
                    poly = Poly3DCollection([corners], alpha=0.8, facecolor='blue', edgecolor='black')
                    ax3d.add_collection3d(poly)
                    
                # Draw 2D Net
                ax2d.set_title(f"Unique Net #{net_number}")
                ax2d.set_xlim([-4, 4]); ax2d.set_ylim([-4, 4])
                ax2d.set_aspect('equal')
                ax2d.grid(True, linestyle=':', alpha=0.6)
                
                for (nx, ny), node_id in current_net_coords.items():
                    rect = Rectangle((nx - 0.5, ny - 0.5), 1, 1, facecolor='lightgreen', edgecolor='black')
                    ax2d.add_patch(rect)
                    ax2d.text(nx, ny, str(node_id), ha='center', va='center')
                
                # --- SAVE THE IMAGE ---
                filename = os.path.join(output_dir, f"cube_net_{net_number:02d}.png")
                plt.savefig(filename, dpi=150, bbox_inches='tight')
                print(f"Found and saved: {filename}")
                
                plt.pause(0.1) # Brief pause so you can see it flash on screen
                
                # Stop when we find all 11 known cube nets
                if net_number == 11:
                    print("\nSuccess! All 11 unique nets found.")
                    break

    plt.ioff()
    plt.show()

if __name__ == "__main__":
    run_simulation()