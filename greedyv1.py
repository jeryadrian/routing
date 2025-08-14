# Optimize Temporary Storage Sites (TSS) Locations
# This script finds optimal locations for temporary storage sites for material origins
# using network analysis. It minimizes total transport effort (tonne-kilometers) 
# considering road network constraints. Materials only move from origins to TSS.

import geopandas as gpd  # For handling geospatial data
import igraph           # For network/graph operations
import pandas as pd     # For data manipulation
from shapely.geometry import Point, LineString
import numpy as np
import os
import time

# Configuration and File Paths
print("--- 1. Configuration ---")

# Input GeoPackage Files:
# - project.gpkg: Material origin points with tonnage information
# - road.gpkg: Road network for routing
# - grid.gpkg: Grid of potential TSS locations
# Adjust file paths as needed
ORIGIN_GPKG = 'input/project.gpkg'
ROAD_NETWORK_GPKG = 'input/road.gpkg'
GRID_INPUT_GPKG = 'input/cch_grid.gpkg'

# Output Directory and Files
# Results include optimal TSS locations, routes, debug information,
# and optimization logs in both GeoPackage and CSV formats
# Adjust output paths as needed
OUTPUT_DIR = 'output/greedyv1_output'
OPTIMAL_TSS_GPKG = os.path.join(OUTPUT_DIR, 'optimal_tss_locations.gpkg')
CANDIDATE_NODES_GPKG = os.path.join(OUTPUT_DIR, 'candidate_tss_nodes_debug.gpkg')
LOADED_GRID_GPKG = os.path.join(OUTPUT_DIR, 'loaded_grid_debug.gpkg')
GRAPH_NODES_GPKG_OPT = os.path.join(OUTPUT_DIR, 'graph_nodes_opt_debug.gpkg')
GRAPH_EDGES_GPKG_OPT = os.path.join(OUTPUT_DIR, 'graph_edges_opt_debug.gpkg')
LOG_CSV_PATH = os.path.join(OUTPUT_DIR, 'optimization_log.csv')
FINAL_ROUTES_CSV = os.path.join(OUTPUT_DIR, 'final_origin_routes.csv')
FINAL_ROUTES_GPKG = os.path.join(OUTPUT_DIR, 'final_origin_routes.gpkg')

# Analysis Parameters
# Adjust these parameters based on your analysis needs
TARGET_CRS = 'EPSG:28992'          # Dutch coordinate system (Rijksdriehoekstelsel)
ORIGIN_ID_COL = 'name'             # Column name for origin point identifiers
ORIGIN_TONNAGE_COL = 'tonnage'     # Column name for material quantities
MAX_SNAP_DISTANCE = 500            # Maximum distance (meters) to snap points to network
NUM_TSS_TO_LOCATE = 17             # Number of temporary storage sites to optimize
DISTANCE_DIVISOR = 1000.0          # Convert distances from meters to kilometers

# Constants for facility type identification in outputs
TSS_TYPE_VALUE = 'TSS'  # Temporary Storage Site

# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Created output directory: {OUTPUT_DIR}")

# Global variables for network graph construction
node_id_counter = 0            # Unique identifier for each node
coords_to_node_id = {}         # Maps coordinates to node IDs for deduplication

print("\n--- 2. Defining Helper Functions ---")

def get_node_id_global(coords, graph_obj, precision=6):
    """
    Get or create a unique node ID for given coordinates in the graph.
    Ensures each unique location has exactly one node in the graph.
    
    Args:
        coords: (x,y) coordinate tuple
        graph_obj: iGraph object to add vertices to
        precision: Decimal places for coordinate rounding (default: 6)
    
    Returns:
        node_id: Integer ID of the node
    """
    global node_id_counter, coords_to_node_id
    rounded_coords = (round(coords[0], precision), round(coords[1], precision))
    node_id = coords_to_node_id.get(rounded_coords)
    if node_id is None:
        # Create new node if coordinates not seen before
        node_id = node_id_counter
        coords_to_node_id[rounded_coords] = node_id
        graph_obj.add_vertex(
            original_nodeID=node_id,
            x=coords[0],
            y=coords[1],
            geometry=Point(coords)
        )
        node_id_counter += 1
    return node_id

def snap_points_to_nodes(points_gdf, nodes_gdf, nodes_sindex, max_snap_dist=250, point_id_col=None):
    """
    Snap points (origins, TSS candidates, etc.) to their nearest network nodes.
    This ensures all facilities are connected to the road network.
    
    Args:
        points_gdf: GeoDataFrame containing points to snap
        nodes_gdf: GeoDataFrame of network nodes to snap to
        nodes_sindex: Spatial index of nodes for efficient nearest neighbor search
        max_snap_dist: Maximum allowed snapping distance in meters
        point_id_col: Column name containing point identifiers for reporting
    
    Returns:
        GeoDataFrame with added columns:
        - nearest_node: ID of closest network node
        - snap_dist: Distance to nearest node
    """
    print(f"  Snapping {len(points_gdf)} points...")
    if points_gdf.empty:
        points_gdf['nearest_node'] = pd.NA; points_gdf['snap_dist'] = np.nan; return points_gdf
    
    # Ensure consistent coordinate systems
    target_crs_gdf = nodes_gdf.crs
    if points_gdf.crs != target_crs_gdf:
         print(f"    Warn: Reprojecting points {point_id_col or ''} to {target_crs_gdf}")
         points_gdf = points_gdf.to_crs(target_crs_gdf)
    
    try:
        # Find nearest nodes within maximum distance
        nearest_output = nodes_sindex.nearest(points_gdf.geometry, return_distance=False, max_distance=max_snap_dist, return_all=False)
        # Handle different output formats from spatial index nearest neighbor search
        # The output format can vary based on the version of the spatial index library
        if isinstance(nearest_output, tuple) and len(nearest_output) == 2:
            # Standard format: (point_indices, node_indices)
            point_indices, node_indices = nearest_output
        else:
            # Handle various alternative output formats
            if len(points_gdf) == len(nearest_output[1]):
                # Format: ([indices], [matched_indices])
                point_indices = points_gdf.index.values[nearest_output[0]]
                node_indices = nearest_output[1]
            else:
                if len(nearest_output) == len(points_gdf) and not isinstance(nearest_output[0], np.ndarray):
                    # Format: [matched_indices] for all points
                    point_indices = points_gdf.index.values
                    node_indices = nearest_output
                elif len(nearest_output[0]) == 0:
                    # No matches found
                    point_indices = np.array([])
                    node_indices = np.array([])
                else:
                    # Default case: assume first format
                    point_indices = nearest_output[0]
                    node_indices = nearest_output[1]

        # Process the nearest neighbor results
        if point_indices.size == 0:
            # No points could be snapped within the maximum distance
            print(f"    Warn: No points found within {max_snap_dist}m during snapping for {point_id_col or 'points'}.")
            snapped_gdf = points_gdf.copy()
            snapped_gdf['nearest_node'] = pd.NA
            snapped_gdf['snap_dist'] = np.nan
            snapped_gdf['node_idx'] = pd.NA
        else:
            # Create DataFrame with point-to-node mappings
            nearest_df = pd.DataFrame({'point_idx': point_indices, 'node_idx': node_indices})
            
            # Handle index type mismatches between points and results
            if points_gdf.index.dtype != nearest_df['point_idx'].dtype:
                if not all(idx in points_gdf.index for idx in nearest_df['point_idx'].unique()):
                    nearest_df['point_idx'] = points_gdf.index[nearest_df['point_idx']]
            
            # Join the nearest node information with original points
            nearest_df = nearest_df.set_index('point_idx')
            snapped_gdf = points_gdf.join(nearest_df)

        # Get node IDs and calculate snapping distances
        # Convert node indices to actual node IDs from the network
        snapped_gdf['nearest_node'] = snapped_gdf['node_idx'].apply(
            lambda idx: nodes_gdf.iloc[int(idx)]['nodeID'] if pd.notna(idx) else pd.NA
        )
        
        # Get geometries of the snapped nodes
        snapped_nodes_geom = snapped_gdf['node_idx'].apply(
            lambda idx: nodes_gdf.iloc[int(idx)].geometry if pd.notna(idx) else None
        )
        snapped_nodes_geoseries = gpd.GeoSeries(snapped_nodes_geom, crs=points_gdf.crs)
        
        # Calculate actual distances between original points and their snapped locations
        snapped_gdf['snap_dist'] = points_gdf.geometry.distance(snapped_nodes_geoseries, align=True)
        snapped_gdf.loc[snapped_gdf['nearest_node'].isna(), 'snap_dist'] = np.nan
        
        # Report snapping results
        num_snapped = snapped_gdf['nearest_node'].notna().sum()
        print(f"    → Snapped {num_snapped} points within {max_snap_dist}m.")
        if num_snapped < len(points_gdf):
            unsnapped_indices = points_gdf.index.difference(snapped_gdf.dropna(subset=['nearest_node']).index)
            print(f"    Warn: Failed to snap {len(unsnapped_indices)} points (Indices: {unsnapped_indices.tolist()[:10]}...).")
        
        return snapped_gdf.drop(columns=['node_idx'], errors='ignore')
    except Exception as e:
        print(f"ERROR during snapping: {e}")
        points_gdf['nearest_node'] = pd.NA
        points_gdf['snap_dist'] = np.nan
        return points_gdf

def geom_from_nodes(node_list, graph_nodes_dict):
    """
    Create a LineString geometry from a list of node IDs, representing a path in the network.
    
    Args:
        node_list: List of node IDs defining the path
        graph_nodes_dict: Dictionary mapping node IDs to their coordinates and properties
    
    Returns:
        LineString geometry connecting the nodes, or None if invalid input/error
    """
    if node_list is None or not isinstance(node_list, list) or len(node_list) < 2:
        return None
    try:
        # Extract coordinates for each node in the path
        coords = [(graph_nodes_dict[n]['x'], graph_nodes_dict[n]['y']) 
                 for n in node_list if n in graph_nodes_dict]
        # Create LineString if we have at least 2 valid points
        return LineString(coords) if len(coords) > 1 else None
    except Exception as e:
        print(f"    Warn: Error creating geometry from node list {node_list}: {e}")
        return None


# --- 3. Data Loading and Initial Processing ---
print("\n--- 3. Data Loading and Initial Processing ---")
start_time_script = time.time()
try:
    # Load all required geospatial datasets (removed CCH loading)
    origins_gdf = gpd.read_file(ORIGIN_GPKG)      # Material origin points
    roads_gdf = gpd.read_file(ROAD_NETWORK_GPKG)  # Road network
    input_grid_gdf = gpd.read_file(GRID_INPUT_GPKG)  # Potential TSS locations grid
    
    print(f"  → Loaded {len(origins_gdf)} origins, "
          f"{len(roads_gdf)} roads, {len(input_grid_gdf)} grid points.")
    
    # Validate required columns exist
    if ORIGIN_ID_COL not in origins_gdf.columns:
        raise ValueError(f"Missing required column {ORIGIN_ID_COL} in origins file")
    # Validate all required columns exist in input files
    if ORIGIN_TONNAGE_COL not in origins_gdf.columns:
        raise ValueError(f"Missing tonnage column {ORIGIN_TONNAGE_COL} in origins file")
    
    # Clean and validate tonnage data
    origins_gdf[ORIGIN_TONNAGE_COL] = pd.to_numeric(origins_gdf[ORIGIN_TONNAGE_COL], errors='coerce')
    origins_gdf = origins_gdf.dropna(subset=[ORIGIN_TONNAGE_COL]).copy()
    print(f"  → Using {len(origins_gdf)} origins after tonnage validation.")
    
    # Ensure all layers use the same coordinate reference system
    print(f"Reprojecting layers to target CRS ({TARGET_CRS})...")
    gdfs_to_reproject = {
        'origins': origins_gdf,
        'roads': roads_gdf,
        'input_grid': input_grid_gdf
    }
    
    # Reproject each layer if needed
    for name, gdf_obj in gdfs_to_reproject.items():
        if gdf_obj.crs is None:
            print(f"Warning: CRS for {name} is None. Assuming {TARGET_CRS}.")
            gdf_obj.crs = TARGET_CRS
        if str(gdf_obj.crs).upper() != TARGET_CRS.upper():
            gdfs_to_reproject[name] = gdf_obj.to_crs(TARGET_CRS)
    
    # Update references to reprojected dataframes
    origins_gdf, roads_gdf, input_grid_gdf = list(gdfs_to_reproject.values())
except Exception as e:
    print(f"ERROR: Data loading/validation: {e}")
    exit()


# --- 4. Process Road Network and Build Graph ---
print("\n--- 4. Process Road Network and Build Graph ---")
try:
    # Clean and prepare the road network
    # 1. Explode multipart geometries into single parts
    roads_exploded = roads_gdf.explode(index_parts=False).reset_index(drop=True)
    # 2. Keep only valid LineString geometries
    roads_lines = roads_exploded[roads_exploded.geometry.geom_type == 'LineString'].copy()
    # 3. Remove null, empty, or too short geometries
    roads_cleaned = roads_lines[
        roads_lines.geometry.notnull() & 
        ~roads_lines.geometry.is_empty & 
        (roads_lines.geometry.length > 0.001)
    ].copy()
    roads_cleaned['length'] = roads_cleaned.geometry.length

    # Densify road segments to ensure adequate network detail
    MAX_SEGMENT_LENGTH = 25  # Maximum distance between nodes in meters
    print(f"Densifying road geometries (max segment length: {MAX_SEGMENT_LENGTH}m)...")
    roads_cleaned.geometry = roads_cleaned.geometry.segmentize(max_segment_length=MAX_SEGMENT_LENGTH)

    # Initialize network graph and tracking variables
    graph = igraph.Graph(directed=False)  # Undirected graph for road network
    edges_list = []  # Track edges for later geometry creation
    print("Building network graph...")
    start_time_graph = time.time()
    node_id_counter = 0
    coords_to_node_id = {}

    # Convert road segments to graph edges
    for _, row in roads_cleaned.iterrows():
        if row.geometry is None or row.geometry.is_empty:
            continue
        
        # Process each segment of the road
        coords = list(row.geometry.coords)
        for i in range(len(coords) - 1):
            s_coords, e_coords = coords[i], coords[i+1]
            if s_coords == e_coords:  # Skip zero-length segments
                continue
                
            # Get or create nodes for segment endpoints
            s_id = get_node_id_global(s_coords, graph)
            e_id = get_node_id_global(e_coords, graph)
            if s_id == e_id:  # Skip self-loops
                continue
                
            # Create edge geometry and calculate length
            s_geom = LineString([s_coords, e_coords])
            s_len = s_geom.length
            # Add valid edges to the graph (non-zero length, no duplicates)
            if s_len > 0:
                # Check if edge already exists to avoid duplicates
                if graph.get_eid(s_id, e_id, directed=False, error=False) == -1:
                    # Store edge data for later export
                    edges_list.append({
                        'geometry': s_geom,
                        'length': s_len,
                        'node_start': s_id,
                        'node_end': e_id
                    })
                    # Add edge to graph with length attribute for routing
                    graph.add_edge(s_id, e_id, length=s_len)

    # Report graph statistics
    print(f"  → Graph: {graph.vcount()} nodes, {graph.ecount()} edges.")
    print(f"  → Graph building time: {time.time() - start_time_graph:.2f} seconds.")
    if graph.vcount() == 0:
        raise ValueError("Empty graph - no valid road segments found.")

    # Create GeoDataFrames for visualization and debugging
    # Extract node data from graph
    nodes_data = [{
        'nodeID': v['original_nodeID'],
        'x': v['x'],
        'y': v['y'],
        'geometry': v['geometry']
    } for v in graph.vs]
    
    # Convert to GeoDataFrames
    nodes_gdf = gpd.GeoDataFrame(nodes_data, geometry='geometry', crs=TARGET_CRS)
    edges_gdf = gpd.GeoDataFrame(edges_list, geometry='geometry', crs=TARGET_CRS)
    
    # Save debug files for visualization
    print(f"Saving graph nodes debug file to {GRAPH_NODES_GPKG_OPT}...")
    nodes_gdf.to_file(GRAPH_NODES_GPKG_OPT, driver='GPKG')
    print(f"Saving graph edges debug file to {GRAPH_EDGES_GPKG_OPT}...")
    edges_gdf.to_file(GRAPH_EDGES_GPKG_OPT, driver='GPKG')
    
    # Create spatial index and node lookup for efficient operations
    nodes_sindex = nodes_gdf.sindex  # For spatial queries
    graph_nodes_dict = {v['original_nodeID']: {
        'x': v['x'],
        'y': v['y'],
        'geometry': v['geometry']
    } for v in graph.vs}  # For quick node lookups
except Exception as e:
    print(f"ERROR: Graph building: {e}")
    exit()


# --- 5. Generate Candidate TSS Locations from Input Grid ---
print("\n--- 5. Generate Candidate TSS Locations from Input Grid ---")
try:
    # Save input grid for visualization/debugging
    print(f"Saving loaded grid debug file to {LOADED_GRID_GPKG}...")
    input_grid_gdf.to_file(LOADED_GRID_GPKG, driver='GPKG')
    
    # Check if grid points have IDs for reference
    grid_id_col = 'id' if 'id' in input_grid_gdf.columns else None
    # Snap grid points to network nodes
    snapped_grid_gdf = snap_points_to_nodes(
        input_grid_gdf,
        nodes_gdf,
        nodes_sindex,
        MAX_SNAP_DISTANCE,
        grid_id_col
    )
    
    # Extract unique network nodes that correspond to grid points
    # These will be our candidate locations for TSS
    candidate_nodes_set = set(
        snapped_grid_gdf.dropna(subset=['nearest_node'])['nearest_node'].astype(int)
    )
    print(f"  → Identified {len(candidate_nodes_set)} unique candidate nodes from input grid.")
    
    if not candidate_nodes_set:
        raise ValueError("No candidate nodes from grid - check snapping distance or grid coverage")
    
    # Create GeoDataFrame of candidate nodes for visualization
    candidate_nodes_gdf = nodes_gdf[nodes_gdf['nodeID'].isin(candidate_nodes_set)].copy()
    print(f"Saving candidate nodes debug file to {CANDIDATE_NODES_GPKG}...")
    candidate_nodes_gdf.to_file(CANDIDATE_NODES_GPKG, driver='GPKG')
except Exception as e:
    print(f"ERROR: Candidate generation: {e}")
    exit()


# --- 6. Snap Origins to Network ---
print("\n--- 6. Snap Origins to Network ---")
try:
    # Snap origin points to nearest network nodes
    origins_snapped_gdf = snap_points_to_nodes(
        origins_gdf,
        nodes_gdf,
        nodes_sindex,
        MAX_SNAP_DISTANCE,
        ORIGIN_ID_COL
    )
    
    # Filter out origins that couldn't be snapped to network
    origins_final = origins_snapped_gdf.dropna(subset=['nearest_node']).copy()
    origins_final['nearest_node'] = origins_final['nearest_node'].astype(int)
    
    if origins_final.empty:
        raise ValueError("No origins could be snapped to network - check snapping distance")
    
    origin_map = dict(zip(origins_final[ORIGIN_ID_COL], origins_final['nearest_node']))
    origin_tonnage_map = dict(zip(origins_final[ORIGIN_ID_COL], origins_final[ORIGIN_TONNAGE_COL]))
except Exception as e: 
    print(f"ERROR: Origin snapping: {e}")
    exit()


# --- 7. Precompute Shortest Path Distances (Sequential) ---
print("\n--- 7. Precompute Shortest Path Distances (Sequential) ---")
precompute_start_time = time.time()
dist_origin_to_candidate = {}
print("  Calculating distances FROM Origins to TSS candidates...")
origin_node_list = list(origin_map.values())

for i, origin_node in enumerate(origin_node_list):
    if not (0 <= origin_node < graph.vcount()):
        dist_origin_to_candidate[origin_node] = {}
        continue
    try:
        distances_matrix = graph.shortest_paths_dijkstra(source=origin_node, target=None, weights='length', mode=igraph.OUT)
        all_distances_from_origin = distances_matrix[0]
        dist_origin_to_candidate[origin_node] = {
            cn: all_distances_from_origin[cn]
            for cn in candidate_nodes_set
            if cn < len(all_distances_from_origin) and all_distances_from_origin[cn] != float('inf')
        }
    except Exception as e:
        print(f"    Error precomputing distances for origin_node {origin_node}: {e}")
        dist_origin_to_candidate[origin_node] = {}
    if (i + 1) % 50 == 0 or (i+1) == len(origin_node_list): 
        print(f"    Processed origin distances for {i+1}/{len(origin_node_list)} origins...")

print(f"  → Distance precomputation time: {time.time() - precompute_start_time:.2f} seconds.")


# --- 8. Optimization Helper Function: Calculates Cost (Origin to TSS Only) ---
print("\n--- 8. Defining Optimization Helper Function (Origin to TSS Only) ---")
def calculate_total_tonne_km_origin_to_tss(
                             active_tss_nodes, origin_map, origin_tonnage_map,
                             dist_origin_to_candidate, distance_divisor):
    """
    Calculate total transport cost for origins going directly to TSS locations.
    Materials are delivered to their nearest TSS and remain there.
    
    Returns:
        tuple: (total_tonne_km, origins_unserved, total_tonnage_served, 
                total_distance_km, origins_served_count)
    """
    total_tonne_km = 0.0
    origins_unserved = 0
    total_tonnage_served = 0.0
    total_distance_km = 0.0
    origins_served_count = 0

    for origin_name, origin_node in origin_map.items():
        tonnage = origin_tonnage_map.get(origin_name, 0)

        # Find nearest TSS for this origin
        min_dist_to_tss = float('inf')
        if origin_node in dist_origin_to_candidate and active_tss_nodes:
            for tss_node in active_tss_nodes:
                current_dist = dist_origin_to_candidate[origin_node].get(tss_node, float('inf'))
                if current_dist < min_dist_to_tss:
                    min_dist_to_tss = current_dist

        if min_dist_to_tss != float('inf') and tonnage > 0:
            # Origin can reach a TSS
            distance_km = min_dist_to_tss / distance_divisor
            tonne_km_contribution = tonnage * distance_km
            
            total_tonne_km += tonne_km_contribution
            total_tonnage_served += tonnage
            total_distance_km += distance_km
            origins_served_count += 1
        else:
            # Origin cannot reach any TSS or has no tonnage
            if tonnage > 0:
                origins_unserved += 1
            
    return (total_tonne_km, origins_unserved, total_tonnage_served, 
            total_distance_km, origins_served_count)

print("  → `calculate_total_tonne_km_origin_to_tss` function defined.")


# --- 9. Optimization Algorithm: Greedy Heuristic (PRIORITIZE TONNAGE) ---
print("\n--- 9. Running Optimization Algorithm (Greedy - Prioritize Tonnage to TSS) ---")
optimization_start_time = time.time()

def find_optimal_tss_locations_greedy(
                                  candidate_nodes, num_tss_to_locate,
                                  origin_map, origin_tonnage_map,
                                  dist_origin_to_candidate, distance_divisor):
    chosen_tss_nodes = set()
    available_candidates = set(candidate_nodes)
    results_log = []

    # Calculate baseline (no TSS available)
    baseline_tk, baseline_unserved, baseline_tonnage_served, baseline_dist_km, baseline_origins_served = calculate_total_tonne_km_origin_to_tss(
        chosen_tss_nodes, origin_map, origin_tonnage_map,
        dist_origin_to_candidate, distance_divisor)
    
    print(f"Baseline (0 TSS):")
    print(f"  Total Tonne-km: {baseline_tk:.2f}")
    print(f"  Total Distance (km): {baseline_dist_km:.2f}")
    print(f"  Origins served: {baseline_origins_served}")
    print(f"  Tonnage served: {baseline_tonnage_served:.2f}")
    print(f"  Origins unserved: {baseline_unserved}")

    results_log.append({
        'num_tss': 0, 
        'tonne_km': baseline_tk, 
        'added_node': None,
        'origins_unserved': baseline_unserved, 
        'tonnage_served': baseline_tonnage_served,
        'total_distance_km': baseline_dist_km,
        'origins_served': baseline_origins_served
    })

    for i in range(num_tss_to_locate):
        iter_time = time.time()
        print(f"\nSelecting TSS #{i+1} of {num_tss_to_locate}...")
        
        best_candidate = None
        best_tonnage_served = -1.0
        best_tk = float('inf')
        best_unserved = float('inf')
        best_distance_km = 0.0
        best_origins_served = -1

        if not available_candidates: 
            print("  No more candidates.")
            break
            
        current_set = chosen_tss_nodes.copy()
        candidates_to_evaluate = list(available_candidates)

        for eval_idx, candidate in enumerate(candidates_to_evaluate):
            temp_set = current_set | {candidate}
            new_tk, new_unserved, new_tonnage_served, new_distance_km, new_origins_served = calculate_total_tonne_km_origin_to_tss(
                temp_set, origin_map, origin_tonnage_map,
                dist_origin_to_candidate, distance_divisor)

            update_best = False

            # Primary criterion: maximize tonnage served
            if new_tonnage_served > best_tonnage_served + 1e-6: 
                update_best = True
            elif abs(new_tonnage_served - best_tonnage_served) < 1e-6:
                # Secondary criterion: minimize tonne-km (transport cost)
                if new_tk < best_tk - 1e-6: 
                    update_best = True

            if update_best:
                best_tonnage_served = new_tonnage_served
                best_tk = new_tk
                best_unserved = new_unserved
                best_distance_km = new_distance_km
                best_origins_served = new_origins_served
                best_candidate = candidate

            if (eval_idx + 1) % 200 == 0 or (eval_idx + 1) == len(candidates_to_evaluate):
                print(f"    Evaluated {eval_idx + 1}/{len(candidates_to_evaluate)} candidates...")

        if best_candidate is not None:
            chosen_tss_nodes.add(best_candidate)
            available_candidates.remove(best_candidate)
            
            print(f"  → Added Node {best_candidate} as TSS #{len(chosen_tss_nodes)}.")
            print(f"  Total Tonne-km with {len(chosen_tss_nodes)} TSS: {best_tk:.2f}")
            print(f"  Total Distance (km): {best_distance_km:.2f}")
            print(f"  Origins served: {best_origins_served}")
            print(f"  Tonnage served: {best_tonnage_served:.2f}")
            print(f"  Origins unserved: {best_unserved}")
            
            results_log.append({
                'num_tss': len(chosen_tss_nodes),
                'tonne_km': best_tk,
                'added_node': best_candidate,
                'origins_unserved': best_unserved,
                'tonnage_served': best_tonnage_served,
                'total_distance_km': best_distance_km,
                'origins_served': best_origins_served
            })
        else:
            print("  No suitable candidate found in this iteration.")
            break  # Exit the optimization loop if no improvement is possible

    print(f"\nGreedy Optimization Finished. Selected {len(chosen_tss_nodes)} TSS.")
    print(f"  → Optimization algorithm time: {time.time() - optimization_start_time:.2f} seconds.")
    
    return chosen_tss_nodes, results_log

# Initialize result tracking variables
selected_tss_nodes = set()
results_log = []

# Run the optimization algorithm
selected_tss_nodes, results_log = find_optimal_tss_locations_greedy(
    candidate_nodes_set,
    NUM_TSS_TO_LOCATE,
    origin_map,
    origin_tonnage_map,
    dist_origin_to_candidate,
    DISTANCE_DIVISOR
)


# --- 10. Save Optimization Results ---
print("\n--- 10. Save Optimization Results ---")
try:
    if selected_tss_nodes:
        # Convert selected TSS nodes to GeoDataFrame
        tss_data = []
        for i, node in enumerate(selected_tss_nodes):
            node_info = graph_nodes_dict[node]
            tss_data.append({
                'id': i,
                'node_id': node,
                'type': TSS_TYPE_VALUE,
                'geometry': node_info['geometry']
            })
        
        tss_points = gpd.GeoDataFrame(tss_data, geometry='geometry', crs=TARGET_CRS)
        
        # Save optimization results
        print(f"Saving optimal TSS locations to {OPTIMAL_TSS_GPKG}...")
        tss_points.to_file(OPTIMAL_TSS_GPKG, driver='GPKG')
    else:
        print("No TSS locations selected.")
    
    # Save optimization log
    if results_log:
        log_df = pd.DataFrame(results_log)
        print(f"Saving optimization log to {LOG_CSV_PATH}...")
        log_df.to_csv(LOG_CSV_PATH, index=False)
    
    print("Successfully saved optimization results.")
except Exception as e:
    print(f"ERROR: Failed to save optimization results: {e}")
    exit()


# --- 11. Perform Final Routing & Save Detailed Routes ---
print("\n--- 11. Perform Final Routing & Save Detailed Routes ---")
final_routing_start_time = time.time()
final_route_results = []

print(f"Generating detailed routes for {len(origin_map)} origins to {len(selected_tss_nodes)} selected TSS locations...")

processed_final_routes = 0
for origin_name, origin_node in origin_map.items():
    tonnage = origin_tonnage_map.get(origin_name, 0)
    route_data = {
        'origin': origin_name, 
        'tonnage': tonnage, 
        'status': 'No_Path_Found',
        'destination_type': None, 
        'destination_loc': None, 
        'distance_km': np.nan,
        'tonne_km': np.nan, 
        'route_nodes': None
    }

    if not (0 <= origin_node < graph.vcount()):
        route_data['status'] = 'Origin_Node_Not_In_Graph'
        final_route_results.append(route_data)
        continue

    try:
        # Calculate distances from origin to all nodes
        all_distances_from_origin = graph.distances(
            source=origin_node, target=None, weights='length', mode=igraph.OUT)[0]

        # Find nearest TSS
        min_dist_to_tss = float('inf')
        nearest_tss_node = None
        path_to_tss = None
        
        if selected_tss_nodes:
            for tss_node in selected_tss_nodes:
                if tss_node < len(all_distances_from_origin) and all_distances_from_origin[tss_node] != float('inf'):
                    dist = all_distances_from_origin[tss_node]
                    if dist < min_dist_to_tss:
                        min_dist_to_tss = dist
                        nearest_tss_node = tss_node
                        try:
                            path_to_tss = graph.get_shortest_paths(
                                origin_node, to=tss_node, weights='length', mode=igraph.OUT, output='vpath')[0]
                        except Exception: 
                            path_to_tss = None

        if nearest_tss_node is not None and path_to_tss and min_dist_to_tss != float('inf'):
            distance_km = min_dist_to_tss / DISTANCE_DIVISOR
            tonne_km = tonnage * distance_km
            
            route_data.update({
                'status': 'Routed_to_TSS',
                'destination_type': TSS_TYPE_VALUE,
                'destination_loc': f"TSS_{nearest_tss_node}",
                'distance_km': distance_km,
                'tonne_km': tonne_km,
                'route_nodes': path_to_tss
            })
        else:
            route_data['status'] = 'No_TSS_Reachable'

    except Exception as e:
        print(f"Error routing origin {origin_name} (node {origin_node}): {e}")
        route_data['status'] = 'Routing_Error'

    final_route_results.append(route_data)
    processed_final_routes += 1
    if processed_final_routes % 50 == 0 or processed_final_routes == len(origin_map):
        print(f"  Generated final routes for {processed_final_routes}/{len(origin_map)} origins...")

print(f"  → Final route generation time: {time.time() - final_routing_start_time:.2f} seconds.")


# --- 12. Save Final Route Results ---
print("\n--- 12. Saving Final Route Results ---")
try:
    if final_route_results:
        final_routes_df = pd.DataFrame(final_route_results)
        
        # Create GeoDataFrame with route geometries
        final_routes_gdf = gpd.GeoDataFrame(
            final_routes_df.drop(columns=['route_nodes']),
            geometry=final_routes_df['route_nodes'].apply(lambda x: geom_from_nodes(x, graph_nodes_dict)),
            crs=TARGET_CRS
        )
        
        # Keep only routes with valid geometries
        final_routes_gdf = final_routes_gdf[
            final_routes_gdf.geometry.notnull() & ~final_routes_gdf.geometry.is_empty
        ].copy()
        
        print(f"Saving final {len(final_routes_gdf)} routes to {FINAL_ROUTES_GPKG}...")
        final_routes_gdf.to_file(FINAL_ROUTES_GPKG, driver='GPKG')
        
        print(f"Saving final {len(final_routes_df)} routes details to {FINAL_ROUTES_CSV}...")
        final_routes_df.to_csv(FINAL_ROUTES_CSV, index=False)
    else:
        print("No final route results to save.")
except Exception as e:
    print(f"ERROR: Saving final route results: {e}")


# --- 13. Final Summary ---
print("\n--- 13. Final Summary ---")
total_time_end = time.time()
print(f"Total script execution time: {total_time_end - start_time_script:.2f} seconds.")

if results_log:
    final_result = results_log[-1]
    
    print(f"\nOptimization Results:")
    print(f"  Selected {len(selected_tss_nodes)} optimal TSS locations")
    print(f"  Final tonne-km: {final_result['tonne_km']:.2f}")
    print(f"  Total distance (km): {final_result['total_distance_km']:.2f}")
    print(f"  Origins served: {final_result['origins_served']}")
    print(f"  Tonnage served: {final_result['tonnage_served']:.2f}")
    print(f"  Origins unserved: {final_result['origins_unserved']}")
    
    # Calculate service coverage
    total_origins = len(origin_map)
    service_rate = (final_result['origins_served'] / total_origins * 100) if total_origins > 0 else 0
    print(f"  Service coverage: {service_rate:.1f}%")
    
    total_tonnage = sum(origin_tonnage_map.values())
    tonnage_coverage = (final_result['tonnage_served'] / total_tonnage * 100) if total_tonnage > 0 else 0
    print(f"  Tonnage coverage: {tonnage_coverage:.1f}%")

print(f"\nResults saved to {OUTPUT_DIR}:")
print(f"  - Optimal TSS locations: {OPTIMAL_TSS_GPKG}")
print(f"  - Optimization log: {LOG_CSV_PATH}")
print(f"  - Detailed routes: {FINAL_ROUTES_GPKG}")
print(f"  - Route data: {FINAL_ROUTES_CSV}")

print("\nOptimization completed! Materials flow from origins directly to TSS locations.")