import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# 1. DEFINE THE 3D GRAPH FOR A 1x1x1 CUBE
graph = {
    0: {0: (1, 2), 1: (2, 2), 2: (3, 2), 3: (4, 2)},
    1: {0: (5, 2), 1: (2, 3), 2: (0, 0), 3: (4, 1)},
    2: {0: (5, 1), 1: (3, 3), 2: (0, 1), 3: (1, 1)},
    3: {0: (5, 0), 1: (4, 3), 2: (0, 2), 3: (2, 1)},
    4: {0: (5, 3), 1: (1, 3), 2: (0, 3), 3: (3, 1)},
    5: {0: (3, 0), 1: (2, 0), 2: (1, 0), 3: (4, 0)},
}

faces_3d = {
    0: [(0,0,0), (1,0,0), (1,1,0), (0,1,0)], 1: [(0,0,0), (1,0,0), (1,0,1), (0,0,1)],
    2: [(1,0,0), (1,1,0), (1,1,1), (1,0,1)], 3: [(1,1,0), (0,1,0), (0,1,1), (1,1,1)],
    4: [(0,1,0), (0,0,0), (0,0,1), (0,1,1)], 5: [(0,0,1), (1,0,1), (1,1,1), (0,1,1)],
}

DIRS = {0: (0, 1), 1: (1, 0), 2: (0, -1), 3: (-1, 0)}


# 2. DEDUPLICATION (Canonicalization)
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


# 3. THE BRANCHING DFS GENERATOR
def dfs_unfold(visited, net_coords, node_states):
    if len(visited) == 6:
        yield ("SUCCESS", None, visited, net_coords)
        return

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
            yield ("COLLISION", next_pos, visited, net_coords)
            continue 
            
        visited.add(neighbor)
        net_coords[next_pos] = neighbor
        node_states[neighbor] = (next_pos, next_offset)
        
        # Yield the addition so we can see it on screen
        yield ("ADDED", next_pos, visited, net_coords)
        
        yield from dfs_unfold(visited, net_coords, node_states)
        
        visited.remove(neighbor)
        del net_coords[next_pos]
        del node_states[neighbor]
        
        # Yield the backtrack so we can see it shrink
        yield ("BACKTRACK", next_pos, visited, net_coords)


# 4. VISUALIZATION AND SAVING ENGINE
def run_simulation():
    output_dir = "net_generator_proof_of_concept_images"
    os.makedirs(output_dir, exist_ok=True)
    
    plt.ion()
    fig = plt.figure(figsize=(12, 6))
    ax3d = fig.add_subplot(121, projection='3d')
    ax2d = fig.add_subplot(122)
    
    visited = {0}
    net_coords = {(0, 0): 0}
    node_states = {0: ((0, 0), 0)}
    dfs_gen = dfs_unfold(visited, net_coords, node_states)
    
    unique_nets_found = set()
    stats = {"collisions": 0, "duplicates": 0, "unique": 0, "attempts": 0}
    
    print(f"Starting visual search... Saving images to ./{output_dir}/\n")
    
    for status, pos, current_visited, current_net_coords in dfs_gen:
        stats["attempts"] += 1
        
        # Logging logic
        if status == "COLLISION":
            stats["collisions"] += 1
            print(f"[{stats['attempts']:04d}] INVALID: Overlap detected")
        elif status == "SUCCESS":
            coords = list(current_net_coords.keys())
            canonical_id = get_canonical_net(coords)
            
            if canonical_id in unique_nets_found:
                stats["duplicates"] += 1
                print(f"[{stats['attempts']:04d}] DUPLICATE: Valid 6-square net, but already found.")
            else:
                stats["unique"] += 1
                unique_nets_found.add(canonical_id)
                print(f">>> [{stats['attempts']:04d}] UNIQUE: Found net #{stats['unique']}! <<<")

        # --- ALWAYS UPDATE GUI ON EVERY STEP ---
        ax3d.clear()
        ax2d.clear()
        
        # Draw 3D Cube (highlight visited faces)
        ax3d.set_title(f"3D Cube\nUnique Nets: {stats['unique']}/11")
        ax3d.set_xlim([0, 1]); ax3d.set_ylim([0, 1]); ax3d.set_zlim([0, 1])
        ax3d.view_init(elev=20, azim=45)
        for face_id, corners in faces_3d.items():
            color = 'blue' if face_id in current_visited else 'whitesmoke'
            alpha = 0.8 if face_id in current_visited else 0.2
            poly = Poly3DCollection([corners], alpha=alpha, facecolor=color, edgecolor='black')
            ax3d.add_collection3d(poly)
            
        # Draw 2D Net
        ax2d.set_title(f"Status: {status} | Attempt: {stats['attempts']}")
        ax2d.set_xlim([-4, 4]); ax2d.set_ylim([-4, 4])
        ax2d.set_aspect('equal')
        ax2d.grid(True, linestyle=':', alpha=0.6)
        
        # Draw the placed squares
        for (nx, ny), node_id in current_net_coords.items():
            face_color = 'lightgreen' if status == "SUCCESS" else 'lightgray'
            rect = Rectangle((nx - 0.5, ny - 0.5), 1, 1, facecolor=face_color, edgecolor='black')
            ax2d.add_patch(rect)
            ax2d.text(nx, ny, str(node_id), ha='center', va='center')
        
        # Highlight collision overlap in red
        if status == "COLLISION":
            rect = Rectangle((pos[0] - 0.5, pos[1] - 0.5), 1, 1, facecolor='red', edgecolor='black', alpha=0.5)
            ax2d.add_patch(rect)
        
        # Dynamic pause logic and saving
        if status == "SUCCESS":
            if canonical_id in unique_nets_found:
                # Save the image only when a brand new unique net is displayed
                filename = os.path.join(output_dir, f"cube_net_{stats['unique']:02d}.png")
                plt.savefig(filename, dpi=150, bbox_inches='tight')
                
            plt.pause(1.5) # Longer pause to admire completed nets (even duplicates)
            
            if stats["unique"] == 11:
                print("\n" + "="*50)
                print(f"SUCCESS! All 11 unique nets found.")
                print(f"Total Collisions Hit: {stats['collisions']}")
                print(f"Total Duplicates Generated: {stats['duplicates']}")
                print(f"Total Unique Nets: {stats['unique']}")
                print("="*50)
                break
        else:
            plt.pause(0.4) # Slow enough to watch the DFS build step-by-step

    plt.ioff()
    plt.show()

if __name__ == "__main__":
    run_simulation()