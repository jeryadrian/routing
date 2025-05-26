# Optimize Temporary Storage Sites (TSS) Locations
# This script finds optimal locations for temporary storage sites between material origins
# and circular construction hubs (CCH) using network analysis. It minimizes total 
# transport effort (tonne-kilometers) considering road network constraints.

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
# - cch.gpkg: Circular Construction Hub locations (final destinations)
# - road.gpkg: Road network for routing
# - grid.gpkg: Grid of potential TSS locations
# - Adjust file paths as needed
ORIGIN_GPKG = 'input/project.gpkg'
CCH_GPKG = 'input/5cch.gpkg' 
ROAD_NETWORK_GPKG = 'input/road.gpkg'
GRID_INPUT_GPKG = 'input/grid.gpkg'

# Output Directory and Files
# Results include optimal TSS locations, routes, debug information,
# and optimization logs in both GeoPackage and CSV formats
# Adjust output paths as needed
OUTPUT_DIR = 'output/greedy_output'
OPTIMAL_TSS_GPKG = os.path.join(OUTPUT_DIR, 'optimal_tss_locations.gpkg')
CANDIDATE_NODES_GPKG = os.path.join(OUTPUT_DIR, 'candidate_tss_nodes_debug.gpkg')
LOADED_GRID_GPKG = os.path.join(OUTPUT_DIR, 'loaded_grid_debug.gpkg')
GRAPH_NODES_GPKG_OPT = os.path.join(OUTPUT_DIR, 'graph_nodes_opt_debug.gpkg')
GRAPH_EDGES_GPKG_OPT = os.path.join(OUTPUT_DIR, 'graph_edges_opt_debug.gpkg')
DEBUG_SNAPPED_CCH_GPKG = os.path.join(OUTPUT_DIR, 'snapped_cch_debug.gpkg')
LOG_CSV_PATH = os.path.join(OUTPUT_DIR, 'optimization_log.csv')
FINAL_ROUTES_CSV = os.path.join(OUTPUT_DIR, 'final_origin_routes.csv')
FINAL_ROUTES_GPKG = os.path.join(OUTPUT_DIR, 'final_origin_routes.gpkg')

# Analysis Parameters
# Adjust these parameters based on your analysis needs
TARGET_CRS = 'EPSG:28992'          # Dutch coordinate system (Rijksdriehoekstelsel)
ORIGIN_ID_COL = 'name'             # Column name for origin point identifiers
ORIGIN_TONNAGE_COL = 'tonnage'     # Column name for material quantities
CCH_ID_COL = 'loc'                 # Column name for CCH identifiers
MAX_SNAP_DISTANCE = 500            # Maximum distance (meters) to snap points to network
NUM_TSS_TO_LOCATE = 10             # Number of temporary storage sites to optimize
DISTANCE_DIVISOR = 1000.0          # Convert distances from meters to kilometers

# Constants for facility type identification in outputs
TSS_TYPE_VALUE = 'TSS'  # Temporary Storage Site
CCH_TYPE_VALUE = 'CCH'  # Circular Construction Hub

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
    Snap points (origins, CCHs, etc.) to their nearest network nodes.
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
    # Load all required geospatial datasets
    origins_gdf = gpd.read_file(ORIGIN_GPKG)      # Material origin points
    cch_gdf = gpd.read_file(CCH_GPKG)            # Circular Construction Hubs
    roads_gdf = gpd.read_file(ROAD_NETWORK_GPKG)  # Road network
    input_grid_gdf = gpd.read_file(GRID_INPUT_GPKG)  # Potential TSS locations grid
    
    print(f"  → Loaded {len(origins_gdf)} origins, {len(cch_gdf)} CCH, "
          f"{len(roads_gdf)} roads, {len(input_grid_gdf)} grid points.")
    
    # Validate required columns exist
    if ORIGIN_ID_COL not in origins_gdf.columns:
        raise ValueError(f"Missing required column {ORIGIN_ID_COL} in origins file")
    # Validate all required columns exist in input files
    if ORIGIN_TONNAGE_COL not in origins_gdf.columns:
        raise ValueError(f"Missing tonnage column {ORIGIN_TONNAGE_COL} in origins file")
    if CCH_ID_COL not in cch_gdf.columns:
        raise ValueError(f"Missing ID column {CCH_ID_COL} in CCH file")
    
    # Clean and validate tonnage data
    origins_gdf[ORIGIN_TONNAGE_COL] = pd.to_numeric(origins_gdf[ORIGIN_TONNAGE_COL], errors='coerce')
    origins_gdf = origins_gdf.dropna(subset=[ORIGIN_TONNAGE_COL]).copy()
    print(f"  → Using {len(origins_gdf)} origins after tonnage validation.")
    
    # Ensure all layers use the same coordinate reference system
    print(f"Reprojecting layers to target CRS ({TARGET_CRS})...")
    gdfs_to_reproject = {
        'origins': origins_gdf,
        'cch': cch_gdf,
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
    origins_gdf, cch_gdf, roads_gdf, input_grid_gdf = list(gdfs_to_reproject.values())
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


# --- 6. Snap Origins and CCHs to Network ---
print("\n--- 6. Snap Origins and CCHs to Network ---")
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
    cch_snapped_gdf = snap_points_to_nodes(cch_gdf, nodes_gdf, nodes_sindex, MAX_SNAP_DISTANCE, CCH_ID_COL)
    cch_final = cch_snapped_gdf.dropna(subset=['nearest_node']).copy()
    cch_final['nearest_node'] = cch_final['nearest_node'].astype(int)
    if cch_final.empty: raise ValueError("No CCHs snapped.")
    print(f"Saving snapped CCH debug file to {DEBUG_SNAPPED_CCH_GPKG}...")
    if not cch_final.empty: cch_final.to_file(DEBUG_SNAPPED_CCH_GPKG, driver='GPKG')
    origin_map = dict(zip(origins_final[ORIGIN_ID_COL], origins_final['nearest_node']))
    origin_tonnage_map = dict(zip(origins_final[ORIGIN_ID_COL], origins_final[ORIGIN_TONNAGE_COL]))
    cch_nodes_set = set(cch_final['nearest_node'])
    node_to_cch_name_map = dict(zip(cch_final['nearest_node'], cch_final[CCH_ID_COL]))
except Exception as e: print(f"ERROR: Origin/CCH snapping: {e}"); exit()


# --- 7. Precompute Shortest Path Distances (Sequential) ---
print("\n--- 7. Precompute Shortest Path Distances (Sequential) ---")
precompute_start_time = time.time()
dist_origin_to_candidate = {}; dist_origin_to_cch = {}; dist_candidate_to_nearest_cch = {}
print("  Calculating distances FROM Origins...")
origin_node_list = list(origin_map.values())

for i, origin_node in enumerate(origin_node_list):
    if not (0 <= origin_node < graph.vcount()):
        dist_origin_to_candidate[origin_node] = {}; dist_origin_to_cch[origin_node] = {}
        continue
    try:
        distances_matrix = graph.shortest_paths_dijkstra(source=origin_node, target=None, weights='length', mode=igraph.OUT)
        all_distances_from_origin = distances_matrix[0]
        dist_origin_to_candidate[origin_node] = {
            cn: all_distances_from_origin[cn]
            for cn in candidate_nodes_set
            if cn < len(all_distances_from_origin) and all_distances_from_origin[cn] != float('inf')
        }
        dist_origin_to_cch[origin_node] = {
            ccn: all_distances_from_origin[ccn]
            for ccn in cch_nodes_set
            if ccn < len(all_distances_from_origin) and all_distances_from_origin[ccn] != float('inf')
        }
    except Exception as e:
        print(f"    Error precomputing distances for origin_node {origin_node}: {e}")
        dist_origin_to_candidate[origin_node] = {}; dist_origin_to_cch[origin_node] = {}
    if (i + 1) % 50 == 0 or (i+1) == len(origin_node_list): print(f"    Processed origin distances for {i+1}/{len(origin_node_list)} origins...")

print("  Calculating distances FROM Candidates TO nearest CCH...")
candidate_node_list = list(candidate_nodes_set)
for i, cand_node in enumerate(candidate_node_list):
    if not (0 <= cand_node < graph.vcount()):
        dist_candidate_to_nearest_cch[cand_node] = float('inf')
        continue
    min_d = float('inf')
    try:
        distances_matrix_cand = graph.shortest_paths_dijkstra(source=cand_node, target=None, weights='length', mode=igraph.OUT)
        all_distances_from_cand = distances_matrix_cand[0]
        for cch_node in cch_nodes_set:
            if cch_node < len(all_distances_from_cand) and all_distances_from_cand[cch_node] != float('inf'):
                min_d = min(min_d, all_distances_from_cand[cch_node])
        dist_candidate_to_nearest_cch[cand_node] = min_d
    except Exception as e:
        print(f"    Error precomputing distances for candidate_node {cand_node}: {e}")
        dist_candidate_to_nearest_cch[cand_node] = float('inf')
    if (i + 1) % 100 == 0 or (i+1) == len(candidate_node_list): print(f"    Processed candidate distances for {i+1}/{len(candidate_node_list)} candidates...")

if dist_candidate_to_nearest_cch:
    can_reach_cch_count = sum(1 for d in dist_candidate_to_nearest_cch.values() if d is not None and d != float('inf'))
    print(f"\n  Check: {can_reach_cch_count} out of {len(dist_candidate_to_nearest_cch)} candidate nodes can reach a CCH.")
    if can_reach_cch_count == 0 and candidate_nodes_set: print("  CRITICAL WARNING: NO candidate nodes can reach any CCH.")
print(f"  → Distance precomputation time: {time.time() - precompute_start_time:.2f} seconds.")


# --- 8. Optimization Helper Function: Calculates Cost ---
print("\n--- 8. Defining Optimization Helper Function (Strict First Leg - V5 Scaled) ---")
def calculate_total_tonne_km_strict_first_leg_v5_scaled(
                             active_tss_nodes, origin_map, origin_tonnage_map,
                             dist_origin_to_candidate, dist_origin_to_cch,
                             dist_candidate_to_nearest_cch, distance_divisor):
    total_tonne_km = 0.0; origins_unserved = 0
    origins_using_tss_count = 0; tonnage_using_tss_sum = 0.0

    tk_origin_to_tss_leg = 0.0
    tk_tss_to_cch_leg = 0.0
    tk_origin_to_cch_direct_leg = 0.0

    sum_dist_origin_to_tss_km = 0.0
    sum_dist_tss_to_cch_km = 0.0
    sum_dist_origin_to_cch_direct_km = 0.0

    for origin_name, origin_node in origin_map.items():
        tonnage = origin_tonnage_map.get(origin_name, 0)

        dist_o_direct_cch = float('inf')
        if origin_node in dist_origin_to_cch and dist_origin_to_cch[origin_node]:
            o_cch_distances = list(dist_origin_to_cch[origin_node].values())
            if o_cch_distances: dist_o_direct_cch = min(o_cch_distances)

        dist_o_leg1_tss = float('inf'); chosen_tss_for_leg1 = None
        if origin_node in dist_origin_to_candidate and active_tss_nodes:
            for active_tss_node_iter in active_tss_nodes:
                current_dist_o_to_active_tss = dist_origin_to_candidate[origin_node].get(active_tss_node_iter, float('inf'))
                if current_dist_o_to_active_tss < dist_o_leg1_tss:
                    dist_o_leg1_tss = current_dist_o_to_active_tss
                    chosen_tss_for_leg1 = active_tss_node_iter

        final_journey_dist_m = float('inf')
        current_origin_tonne_km_contribution = 0.0

        if dist_o_direct_cch <= dist_o_leg1_tss:
            if dist_o_direct_cch != float('inf'):
                final_journey_dist_m = dist_o_direct_cch
                journey_km = final_journey_dist_m / distance_divisor
                sum_dist_origin_to_cch_direct_km += journey_km
                if tonnage > 0:
                    current_origin_tonne_km_contribution = tonnage * journey_km
                    tk_origin_to_cch_direct_leg += current_origin_tonne_km_contribution
        else:
            if chosen_tss_for_leg1 is not None and dist_o_leg1_tss != float('inf'):
                dist_leg2_tss_cch = dist_candidate_to_nearest_cch.get(chosen_tss_for_leg1, float('inf'))
                if dist_leg2_tss_cch != float('inf'):
                    final_journey_dist_m = dist_o_leg1_tss + dist_leg2_tss_cch
                    leg1_km = dist_o_leg1_tss / distance_divisor
                    leg2_km = dist_leg2_tss_cch / distance_divisor
                    sum_dist_origin_to_tss_km += leg1_km
                    sum_dist_tss_to_cch_km += leg2_km
                    if tonnage > 0:
                        current_leg1_tk = tonnage * leg1_km
                        current_leg2_tk = tonnage * leg2_km
                        tk_origin_to_tss_leg += current_leg1_tk
                        tk_tss_to_cch_leg += current_leg2_tk
                        current_origin_tonne_km_contribution = current_leg1_tk + current_leg2_tk
                        origins_using_tss_count += 1
                        tonnage_using_tss_sum += tonnage

        if tonnage > 0:
            if final_journey_dist_m != float('inf'):
                total_tonne_km += current_origin_tonne_km_contribution
            else:
                origins_unserved += 1
            
    return (total_tonne_km, origins_unserved, origins_using_tss_count, tonnage_using_tss_sum,
            tk_origin_to_tss_leg, tk_tss_to_cch_leg, tk_origin_to_cch_direct_leg,
            sum_dist_origin_to_tss_km, sum_dist_tss_to_cch_km, sum_dist_origin_to_cch_direct_km)
print("  → `calculate_total_tonne_km_strict_first_leg_v5_scaled` function defined (with tonne-km and total distance breakdowns).")


# --- 9. Optimization Algorithm: Greedy Heuristic (PRIORITIZE TONNAGE) ---
print("\n--- 9. Running Optimization Algorithm (Greedy - Prioritize Tonnage Using TSS) ---")
optimization_start_time = time.time()
def find_optimal_locations_greedy_prioritize_tonnage(
                                  candidate_nodes, num_tss_to_locate,
                                  origin_map, origin_tonnage_map,
                                  dist_origin_to_candidate, dist_origin_to_cch,
                                  dist_candidate_to_nearest_cch, distance_divisor):
    chosen_tss_nodes = set(); available_candidates = set(candidate_nodes); results_log = []

    baseline_tk, baseline_unserved, baseline_origins_using_tss, baseline_tonnage_using_tss, \
    baseline_tk_o_tss, baseline_tk_tss_cch, baseline_tk_o_cch_direct, \
    baseline_dist_o_tss_km, baseline_dist_tss_cch_km, baseline_dist_o_cch_direct_km = calculate_total_tonne_km_strict_first_leg_v5_scaled(
        chosen_tss_nodes, origin_map, origin_tonnage_map,
        dist_origin_to_candidate, dist_origin_to_cch, dist_candidate_to_nearest_cch, distance_divisor)
    
    baseline_total_dist_km = baseline_dist_o_cch_direct_km + baseline_dist_o_tss_km + baseline_dist_tss_cch_km
    
    print(f"Baseline (0 TSS):")
    print(f"  Total Tonne-km: {baseline_tk:.2f}")
    print(f"  Total Route Distance (km): {baseline_total_dist_km:.2f}")
    print(f"  Breakdown:")
    print(f"    - Sum Dist (Direct O-CCH): {baseline_dist_o_cch_direct_km:.2f} km")
    print(f"    - Sum Dist (Via TSS, O-TSS leg): {baseline_dist_o_tss_km:.2f} km")
    print(f"    - Sum Dist (Via TSS, TSS-CCH leg): {baseline_dist_tss_cch_km:.2f} km")
    print(f"  ({baseline_unserved} unserved, {baseline_origins_using_tss} origins, {baseline_tonnage_using_tss:.2f}t using TSS)")

    results_log.append({'num_tss': 0, 'tonne_km': baseline_tk, 'added_node': None,
                        'origins_unserved': baseline_unserved, 'origins_using_tss': baseline_origins_using_tss,
                        'tonnage_using_tss': baseline_tonnage_using_tss,
                        'tk_o_tss': baseline_tk_o_tss, 'tk_tss_cch': baseline_tk_tss_cch,
                        'tk_o_cch_direct': baseline_tk_o_cch_direct,
                        'dist_o_tss_km': baseline_dist_o_tss_km,
                        'dist_tss_cch_km': baseline_dist_tss_cch_km,
                        'dist_o_cch_direct_km': baseline_dist_o_cch_direct_km })

    for i in range(num_tss_to_locate):
        iter_time = time.time(); print(f"\nSelecting TSS #{i+1} of {num_tss_to_locate}...")
        best_cand_iter = None; best_iter_tonnage_using_tss = -1.0; best_iter_cand_dist_to_cch = float('inf')
        best_iter_tk_for_best_cand = float('inf'); best_iter_unserved_for_best_cand = float('inf'); best_iter_origins_count_for_best = -1
        
        best_iter_tk_o_tss = 0.0
        best_iter_tk_tss_cch = 0.0
        best_iter_tk_o_cch_direct = 0.0
        best_iter_dist_o_tss_km = 0.0
        best_iter_dist_tss_cch_km = 0.0
        best_iter_dist_o_cch_direct_km = 0.0

        if not available_candidates: print("  No more candidates."); break
        current_eval_set_base = chosen_tss_nodes.copy()
        candidates_to_evaluate = list(available_candidates)

        for eval_idx, cand_being_tested in enumerate(candidates_to_evaluate):
            temp_set_for_eval = current_eval_set_base | {cand_being_tested}
            new_tk, new_unserved, new_origins_using_tss, new_tonnage_using_tss, \
            current_tk_o_tss, current_tk_tss_cch, current_tk_o_cch_direct, \
            current_dist_o_tss_km, current_dist_tss_cch_km, current_dist_o_cch_direct_km = calculate_total_tonne_km_strict_first_leg_v5_scaled(
                temp_set_for_eval, origin_map, origin_tonnage_map,
                dist_origin_to_candidate, dist_origin_to_cch, dist_candidate_to_nearest_cch, distance_divisor)

            cand_dist_to_cch = dist_candidate_to_nearest_cch.get(cand_being_tested, float('inf'))
            update_best = False

            if new_tonnage_using_tss > best_iter_tonnage_using_tss + 1e-6: update_best = True
            elif abs(new_tonnage_using_tss - best_iter_tonnage_using_tss) < 1e-6:
                if cand_dist_to_cch < best_iter_cand_dist_to_cch - 1e-6: update_best = True
                elif abs(cand_dist_to_cch - best_iter_cand_dist_to_cch) < 1e-6:
                    if new_tk < best_iter_tk_for_best_cand - 1e-6: update_best = True

            if update_best:
                best_iter_tonnage_using_tss = new_tonnage_using_tss
                best_iter_cand_dist_to_cch = cand_dist_to_cch
                best_iter_tk_for_best_cand = new_tk
                best_iter_unserved_for_best_cand = new_unserved
                best_iter_origins_count_for_best = new_origins_using_tss
                best_cand_iter = cand_being_tested
                best_iter_tk_o_tss = current_tk_o_tss
                best_iter_tk_tss_cch = current_tk_tss_cch
                best_iter_tk_o_cch_direct = current_tk_o_cch_direct
                best_iter_dist_o_tss_km = current_dist_o_tss_km
                best_iter_dist_tss_cch_km = current_dist_tss_cch_km
                best_iter_dist_o_cch_direct_km = current_dist_o_cch_direct_km

            if (eval_idx + 1) % 200 == 0 or (eval_idx + 1) == len(candidates_to_evaluate):
                print(f"    Evaluated {eval_idx + 1}/{len(candidates_to_evaluate)} candidates...")

        if best_cand_iter is not None:
            chosen_tss_nodes.add(best_cand_iter)
            available_candidates.remove(best_cand_iter)
            selected_tss_nodes.add(best_cand_iter)
            best_iter_total_dist_km = best_iter_dist_o_cch_direct_km + best_iter_dist_o_tss_km + best_iter_dist_tss_cch_km
            print(f"  → Added Node {best_cand_iter} as TSS #{len(chosen_tss_nodes)} (Prioritized Tonnage).")
            print(f"  Total Tonne-km with {len(chosen_tss_nodes)} TSS: {best_iter_tk_for_best_cand:.2f}")
            print(f"  Total Route Distance (km): {best_iter_total_dist_km:.2f}")
            print(f"  Breakdown:")
            print(f"    - Sum Dist (Direct O-CCH): {best_iter_dist_o_cch_direct_km:.2f} km")
            print(f"    - Sum Dist (Via TSS, O-TSS leg): {best_iter_dist_o_tss_km:.2f} km")
            print(f"    - Sum Dist (Via TSS, TSS-CCH leg): {best_iter_dist_tss_cch_km:.2f} km")
            print(f"  ({best_iter_unserved_for_best_cand} unserved, {best_iter_origins_count_for_best} origins using TSS, {best_iter_tonnage_using_tss:.2f}t using TSS)")
            
            # Record metrics for this iteration's best candidate
            # This helps track the optimization progress and final solution quality
            tkm_results.append(best_iter_tk_for_best_cand)
            tonnage_results.append(best_iter_tonnage_using_tss)
            candidate_points_used.append(best_cand_iter)
            
            # Print progress feedback for monitoring long-running optimizations
            print(f"Iteration {len(candidate_points_used)}")
            print(f"Best TKm: {best_iter_tk_for_best_cand:,.2f}")
            print(f"Tonnage through TSS: {best_iter_tonnage_using_tss:,.2f}")
            
            # Break if we've reached the desired number of TSS locations
            if len(chosen_tss_nodes) >= num_tss_to_locate:
                break
        else:
            print("  No suitable candidate found in this iteration.")
            break  # Exit the optimization loop if no improvement is possible

    print(f"\nGreedy Prioritize Tonnage Finished. Selected {len(chosen_tss_nodes)} TSS.")
    
    final_tk_overall = 0.0; final_unserved_overall = 0; final_n_tss_overall = 0
    final_origins_using_tss_overall = 0; final_tonnage_using_tss_overall = 0.0
    final_tk_o_tss_overall = 0.0; final_tk_tss_cch_overall = 0.0; final_tk_o_cch_direct_overall = 0.0
    final_dist_o_tss_km_overall = 0.0; final_dist_tss_cch_km_overall = 0.0; final_dist_o_cch_direct_km_overall = 0.0

    if results_log:
        final_log_entry = results_log[-1]
        final_tk_overall = final_log_entry['tonne_km']
        final_unserved_overall = final_log_entry['origins_unserved']
        final_n_tss_overall = final_log_entry['num_tss']
        final_origins_using_tss_overall = final_log_entry.get('origins_using_tss', 0)
        final_tonnage_using_tss_overall = final_log_entry.get('tonnage_using_tss', 0.0)
        final_tk_o_tss_overall = final_log_entry.get('tk_o_tss', 0.0)
        final_tk_tss_cch_overall = final_log_entry.get('tk_tss_cch', 0.0)
        final_tk_o_cch_direct_overall = final_log_entry.get('tk_o_cch_direct', 0.0)
        final_dist_o_tss_km_overall = final_log_entry.get('dist_o_tss_km', 0.0)
        final_dist_tss_cch_km_overall = final_log_entry.get('dist_tss_cch_km', 0.0)
        final_dist_o_cch_direct_km_overall = final_log_entry.get('dist_o_cch_direct_km', 0.0)

    print(f"  → Optimization algorithm time: {time.time() - optimization_start_time:.2f} seconds.")
    return chosen_tss_nodes, final_tk_overall, pd.DataFrame(results_log)

# Initialize variables that will be assigned by the optimization function
selected_tss_nodes = set()
final_min_tonne_km = float('inf') 
results_log_df = pd.DataFrame()

# Initialize result tracking variables
tkm_results = []  # Track transport-kilometer results for each iteration
tonnage_results = []  # Track tonnage routed through TSS for each iteration
candidate_points_used = []  # Store selected TSS locations

# Set the target number of TSS locations to find
num_tss = 3  # Can be adjusted based on requirements

# Main optimization loop to find best TSS locations
while True:
    best_iter_tk_for_best_cand = float('inf')
    best_iter_tonnage_using_tss = 0
    
    # Run the optimization algorithm
    find_optimal_locations_greedy_prioritize_tonnage(
        candidate_nodes_set,
        NUM_TSS_TO_LOCATE,
        origin_map,
        origin_tonnage_map,
        dist_origin_to_candidate,
        dist_origin_to_cch,
        dist_candidate_to_nearest_cch,
        DISTANCE_DIVISOR
    )
    
    if len(tkm_results) >= NUM_TSS_TO_LOCATE:
        print("Reached target number of TSS locations.")
        break
    elif not tkm_results or tkm_results[-1] >= final_min_tonne_km:
        print("No further improvement possible.")
        break
    
    final_min_tonne_km = tkm_results[-1]


# --- 10. Save Optimization Results ---
print("\n--- 10. Save Optimization Results ---")
try:
    # Convert selected TSS nodes to GeoDataFrame
    tss_points = gpd.GeoDataFrame(
        geometry=[graph_nodes_dict[node]['geometry'] for node in selected_tss_nodes],
        crs=TARGET_CRS
    )
    
    # Add metadata columns
    tss_points['id'] = range(len(tss_points))
    tss_points['type'] = TSS_TYPE_VALUE
    tss_points['tkm'] = tkm_results
    tss_points['tonnage'] = tonnage_results
    
    # Save optimization results
    print(f"Saving optimal TSS locations to {OPTIMAL_TSS_GPKG}...")
    tss_points.to_file(OPTIMAL_TSS_GPKG, driver='GPKG')
    
    # Save optimization log
    log_df = pd.DataFrame({
        'iteration': range(len(tkm_results)),
        'tonne_km': tkm_results,
        'tonnage_using_tss': tonnage_results
    })
    print(f"Saving optimization log to {LOG_CSV_PATH}...")
    log_df.to_csv(LOG_CSV_PATH, index=False)
    
    print("Successfully saved all optimization results.")
except Exception as e:
    print(f"ERROR: Failed to save optimization results: {e}")
    exit()


# --- 11. Perform Final Routing & Save Detailed Routes (Updated) ---
print("\n--- 11. Perform Final Routing & Save Detailed Routes ---")
final_routing_start_time = time.time()
final_route_results = []

sum_tonne_km_origin_tss = 0.0
sum_tonne_km_tss_cch = 0.0
sum_tonne_km_origin_cch_direct = 0.0
# You could add distance sums here too if needed for cross-checking Sec 13 with Sec 11
# sum_dist_km_origin_tss_sec11 = 0.0
# sum_dist_km_tss_cch_sec11 = 0.0
# sum_dist_km_origin_cch_direct_sec11 = 0.0


final_active_tss_nodes = selected_tss_nodes
print(f"Generating detailed routes for {len(origin_map)} origins using the {len(final_active_tss_nodes)} selected TSS locations...")

processed_final_routes = 0
for origin_name, origin_node in origin_map.items():
    tonnage = origin_tonnage_map.get(origin_name, 0)
    route_data = {'origin': origin_name, 'tonnage': tonnage, 'status': 'No_Path_Found',
                  'leg1_type': None, 'leg1_loc': None, 'leg1_dist': np.nan,
                  'leg2_type': None, 'leg2_loc': None, 'leg2_dist': np.nan,
                  'total_dist': np.nan, 'total_tonne_km': np.nan, 'route_nodes': None}

    if not (0 <= origin_node < graph.vcount()):
        route_data['status']='Origin_Node_Not_In_Graph'; final_route_results.append(route_data); continue

    try:
        all_distances_from_origin = graph.distances(
            source=origin_node, target=None, weights='length', mode=igraph.OUT)[0]

        min_dist_o_direct_cch = float('inf'); nearest_direct_cch_node = None; path_o_direct_cch = None
        for cch_node in cch_nodes_set:
            if cch_node < len(all_distances_from_origin) and all_distances_from_origin[cch_node] != float('inf'):
                dist = all_distances_from_origin[cch_node]
                if dist < min_dist_o_direct_cch:
                    min_dist_o_direct_cch = dist
                    nearest_direct_cch_node = cch_node
                    try:
                        path_o_direct_cch = graph.get_shortest_paths(
                            origin_node, to=cch_node, weights='length', mode=igraph.OUT, output='vpath')[0] 
                    except Exception: path_o_direct_cch = None 

        min_dist_o_leg1_tss = float('inf'); chosen_leg1_tss_node = None; path_o_leg1_tss = None
        if final_active_tss_nodes:
            for tss_node in final_active_tss_nodes:
                if tss_node < len(all_distances_from_origin) and all_distances_from_origin[tss_node] != float('inf'):
                    dist = all_distances_from_origin[tss_node]
                    if dist < min_dist_o_leg1_tss:
                        min_dist_o_leg1_tss = dist
                        chosen_leg1_tss_node = tss_node
                        try:
                            path_o_leg1_tss = graph.get_shortest_paths(
                                origin_node, to=tss_node, weights='length', mode=igraph.OUT, output='vpath')[0]
                        except Exception: path_o_leg1_tss = None

        if min_dist_o_direct_cch <= min_dist_o_leg1_tss:
            if nearest_direct_cch_node is not None and path_o_direct_cch: 
                current_total_dist_km = min_dist_o_direct_cch / DISTANCE_DIVISOR
                current_total_tonne_km = (tonnage * min_dist_o_direct_cch) / DISTANCE_DIVISOR
                route_data.update({
                    'status': 'Direct_CCH', 'leg1_type': CCH_TYPE_VALUE,
                    'leg1_loc': node_to_cch_name_map.get(nearest_direct_cch_node, f"CCH_{nearest_direct_cch_node}"),
                    'leg1_dist': current_total_dist_km,
                    'total_dist': current_total_dist_km,
                    'total_tonne_km': current_total_tonne_km,
                    'route_nodes': path_o_direct_cch })
                if pd.notna(current_total_tonne_km) and tonnage > 0:
                    sum_tonne_km_origin_cch_direct += current_total_tonne_km
                # if pd.notna(current_total_dist_km): sum_dist_km_origin_cch_direct_sec11 += current_total_dist_km

        else:
            if chosen_leg1_tss_node is not None and path_o_leg1_tss:
                dist_leg2 = float('inf'); path_leg2 = None; nearest_cch_from_tss = None
                try:
                    all_distances_from_tss = graph.distances(
                        source=chosen_leg1_tss_node, target=None, weights='length', mode=igraph.OUT)[0]
                    for cch_node in cch_nodes_set:
                        if cch_node < len(all_distances_from_tss) and all_distances_from_tss[cch_node] != float('inf'):
                            dist = all_distances_from_tss[cch_node]
                            if dist < dist_leg2:
                                dist_leg2 = dist
                                nearest_cch_from_tss = cch_node
                                try:
                                    path_leg2 = graph.get_shortest_paths(
                                        chosen_leg1_tss_node, to=cch_node, weights='length', mode=igraph.OUT, output='vpath')[0]
                                except Exception: path_leg2 = None
                except Exception as e_tss_route: print(f"    Error routing from TSS {chosen_leg1_tss_node}: {e_tss_route}")

                if nearest_cch_from_tss is not None and path_leg2:
                    full_path_nodes = path_o_leg1_tss + path_leg2[1:]
                    total_dist_m = min_dist_o_leg1_tss + dist_leg2
                    leg1_dist_km = min_dist_o_leg1_tss / DISTANCE_DIVISOR
                    leg2_dist_km = dist_leg2 / DISTANCE_DIVISOR
                    current_total_dist_km = total_dist_m / DISTANCE_DIVISOR
                    current_total_tonne_km = (tonnage * total_dist_m) / DISTANCE_DIVISOR
                    
                    route_data.update({
                        'status': 'Via_TSS', 'leg1_type': TSS_TYPE_VALUE,
                        'leg1_loc': f"TSS_{chosen_leg1_tss_node}",
                        'leg1_dist': leg1_dist_km,
                        'leg2_type': CCH_TYPE_VALUE,
                        'leg2_loc': node_to_cch_name_map.get(nearest_cch_from_tss, f"CCH_{nearest_cch_from_tss}"),
                        'leg2_dist': leg2_dist_km,
                        'total_dist': current_total_dist_km,
                        'total_tonne_km': current_total_tonne_km,
                        'route_nodes': full_path_nodes })
                    
                    if tonnage > 0:
                        if pd.notna(leg1_dist_km): sum_tonne_km_origin_tss += tonnage * leg1_dist_km
                        if pd.notna(leg2_dist_km): sum_tonne_km_tss_cch += tonnage * leg2_dist_km
                    # if pd.notna(leg1_dist_km): sum_dist_km_origin_tss_sec11 += leg1_dist_km
                    # if pd.notna(leg2_dist_km): sum_dist_km_tss_cch_sec11 += leg2_dist_km
                else:
                    route_data['status'] = 'Via_TSS_No_CCH_Path'
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
        final_routes_gdf = gpd.GeoDataFrame(
            final_routes_df.drop(columns=['route_nodes']),
            geometry=final_routes_df['route_nodes'].apply(lambda x: geom_from_nodes(x, graph_nodes_dict)),
            crs=TARGET_CRS
        )
        final_routes_gdf = final_routes_gdf[final_routes_gdf.geometry.notnull() & ~final_routes_gdf.geometry.is_empty].copy()
        print(f"Saving final {len(final_routes_gdf)} routes to {FINAL_ROUTES_GPKG}...")
        final_routes_gdf.to_file(FINAL_ROUTES_GPKG, driver='GPKG')
        print(f"Saving final {len(final_routes_df)} routes details to {FINAL_ROUTES_CSV}...")
        final_routes_df.to_csv(FINAL_ROUTES_CSV, index=False)
    else:
        print("No final route results to save.")
except Exception as e:
    print(f"ERROR: Saving final route results: {e}")


# --- 13. Final Message ---
print("\n--- 13. Finished ---")
total_time_end = time.time()
print(f"Total script execution time: {total_time_end - start_time_script:.2f} seconds.")

if not results_log_df.empty:
    final_log_entry = results_log_df.iloc[-1]
    num_actually_selected_tss = len(selected_tss_nodes) 
    
    final_tk_total_from_log = final_log_entry['tonne_km']
    final_origins_using_tss_count = final_log_entry['origins_using_tss']
    final_tonnage_using_tss_val = final_log_entry['tonnage_using_tss']
    final_unserved_count = final_log_entry['origins_unserved']
    
    final_dist_o_tss_km_log = final_log_entry.get('dist_o_tss_km', 0.0)
    final_dist_tss_cch_km_log = final_log_entry.get('dist_tss_cch_km', 0.0)
    final_dist_o_cch_direct_km_log = final_log_entry.get('dist_o_cch_direct_km', 0.0)
    final_total_dist_km_log = final_dist_o_cch_direct_km_log + final_dist_o_tss_km_log + final_dist_tss_cch_km_log

    print(f"Found {num_actually_selected_tss} optimal TSS locations.")
    print(f"Final State (based on log for {final_log_entry['num_tss']} TSS):")
    print(f"  Total Tonne-km: {final_tk_total_from_log:.2f}")
    print(f"  Total Route Distance (km): {final_total_dist_km_log:.2f}")
    print(f"  Breakdown:")
    print(f"    - Sum Dist (Direct O-CCH): {final_dist_o_cch_direct_km_log:.2f} km")
    print(f"    - Sum Dist (Via TSS, O-TSS leg): {final_dist_o_tss_km_log:.2f} km")
    print(f"    - Sum Dist (Via TSS, TSS-CCH leg): {final_dist_tss_cch_km_log:.2f} km")
    print(f"  ({final_unserved_count} unserved, {final_origins_using_tss_count} origins using TSS, {final_tonnage_using_tss_val:.2f}t using TSS).")

elif final_min_tonne_km != float('inf'): # Check if it was ever updated from its initial 'inf'
    print(f"No TSS selected or optimization log DataFrame was empty. Overall Tonne-km (baseline or last valid): {final_min_tonne_km:.2f}.")
else:
    print("Optimization results or log not available, or optimization did not run successfully to produce a tonne-km value.")

print(f"\nCheck {OUTPUT_DIR} for optimization log, selected TSS locations, and detailed final routes.")

print("\nOptimization completed!")
print(f"Found {len(candidate_points_used)} TSS locations")
print(f"Final TKm: {tkm_results[-1]:,.2f}")
print(f"Final tonnage through TSS: {tonnage_results[-1]:,.2f}")
