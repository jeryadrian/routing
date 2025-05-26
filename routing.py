"""
This script performs optimized routing calculations using the python-igraph library.
It finds the shortest paths from origins to Circular Construction Hubs (CCH),
potentially via Temporary Storage Sites (TSS), based on a given road network.
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import os
import time
import igraph
from shapely.geometry import LineString, MultiLineString
import warnings

# Suppress UserWarnings, often from GeoPandas/Shapely operations
warnings.filterwarnings('ignore', category=UserWarning)

print("--- 1. Configuration & Constants ---")

# Define input GeoPackage file paths for origins, TSS, CCH, and the road network.
ORIGIN_GPKG = 'input/project.gpkg'
TSS_GPKG = 'input/8tss.gpkg'
CCH_GPKG = 'input/3cch.gpkg'
ROAD_NETWORK_GPKG = 'input/road.gpkg'

# Define output directory and file paths for routing results.
OUTPUT_DIR = 'output/routing_output'
RESULTS_CSV = os.path.join(OUTPUT_DIR, 'routing_results.csv')
RESULTS_GPKG = os.path.join(OUTPUT_DIR, 'routing_results.gpkg')

# Define constants for Coordinate Reference System (CRS) and facility types.
TARGET_CRS = 'EPSG:28992'
TSS_TYPE_VALUE = 'TSS'
CCH_TYPE_VALUE = 'CCH'

# Create the output directory if it does not already exist.
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Created output directory: {OUTPUT_DIR}")

print("\n--- 2. Load Input Data ---")
try:
    # Load geographical data from GeoPackage files into GeoDataFrames.
    origins_gdf = gpd.read_file(ORIGIN_GPKG)
    tss_gdf = gpd.read_file(TSS_GPKG)
    cch_gdf = gpd.read_file(CCH_GPKG)
    roads_gdf = gpd.read_file(ROAD_NETWORK_GPKG)
    print(f"  → Loaded {len(origins_gdf)} origins, {len(tss_gdf)} TSS, {len(cch_gdf)} CCH, {len(roads_gdf)} roads.")
except Exception as e:
    print(f"ERROR: Failed to load input files: {e}")
    exit()

# Define required column names for input GeoDataFrames.
required_cols = {'origins': 'name', 'tss': 'loc', 'cch': 'loc'}

# Verify that all required columns exist in the loaded GeoDataFrames.
for df, name in [(origins_gdf, 'origins'), (tss_gdf, 'tss'), (cch_gdf, 'cch')]:
    if required_cols[name] not in df.columns:
        print(f"ERROR: {name} GeoPackage missing required column '{required_cols[name]}'")
        exit()

print("\n--- 3. Build Graph ---")
try:
    # Initialize an undirected igraph Graph.
    g = igraph.Graph(directed=False)
    
    # Dictionaries to map coordinates to node IDs and vice-versa.
    coords_to_node = {}
    node_to_coords = {} 
    node_counter = 0
    
    def get_line_coords(geom):
        """
        Extracts coordinates from LineString or MultiLineString geometry objects.
        Handles both single and multiple line geometries.
        """
        if isinstance(geom, LineString):
            return [coord for coord in geom.coords]
        elif isinstance(geom, MultiLineString):
            return [coord for line in geom.geoms for coord in line.coords]
        else:
            raise ValueError(f"Unsupported geometry type: {type(geom)}")
    
    print("Building network nodes...")
    # First pass: Iterate through road geometries to collect all unique coordinates
    # and assign them a unique node ID in the graph.
    for idx, row in roads_gdf.iterrows():
        try:
            coords = get_line_coords(row.geometry)
            for coord in coords:
                if coord not in coords_to_node:
                    coords_to_node[coord] = node_counter
                    node_to_coords[node_counter] = coord 
                    node_counter += 1
                    g.add_vertex() # Add a new vertex to the graph for each unique coordinate
        except Exception as e:
            print(f"Warning: Skipping road segment due to geometry error: {e}")
            continue
    
    # Lists to store edge tuples (node_id, node_id) and their corresponding weights (lengths).
    edge_tuples = []
    edge_weights = []
    
    print("Building network edges...")
    # Second pass: Iterate through road geometries to create edges between connected nodes.
    # Edge weights are calculated as the Euclidean distance between consecutive coordinates.
    for idx, row in roads_gdf.iterrows():
        try:
            # Handle both LineString and MultiLineString geometries.
            if isinstance(row.geometry, MultiLineString):
                for line in row.geometry.geoms:
                    coords = list(line.coords)
                    for i in range(len(coords)-1):
                        n1 = coords_to_node[coords[i]]
                        n2 = coords_to_node[coords[i+1]]
                        dist = ((coords[i+1][0] - coords[i][0])**2 + (coords[i+1][1] - coords[i][1])**2)**0.5
                        edge_tuples.append((n1, n2))
                        edge_weights.append(dist)
            else:  # LineString
                coords = list(row.geometry.coords)
                for i in range(len(coords)-1):
                    n1 = coords_to_node[coords[i]]
                    n2 = coords_to_node[coords[i+1]]
                    dist = ((coords[i+1][0] - coords[i][0])**2 + (coords[i+1][1] - coords[i][1])**2)**0.5
                    edge_tuples.append((n1, n2))
                    edge_weights.append(dist)
        except Exception as e:
            print(f"Warning: Skipping road segment {idx} due to error: {e}")
            continue
    
    # Add all collected edges to the graph and assign their calculated lengths as weights.
    g.add_edges(edge_tuples)
    g.es['length'] = edge_weights
    
    print(f"  → Graph: {g.vcount()} nodes, {g.ecount()} edges")
    
    print("\n--- 4. Setup Routing ---")
    
    def snap_to_network(points_gdf):
        """
        Snaps geographical points (origins, TSS, CCH) to the nearest node in the road network graph.
        Returns a list of the nearest node IDs for each point.
        """
        snapped = []
        for _, point in points_gdf.iterrows():
            min_dist = float('inf')
            nearest_node = None
            pt_coords = (point.geometry.x, point.geometry.y)
            
            # Iterate through all network nodes to find the closest one.
            for coord, node_id in coords_to_node.items():
                dist = ((coord[0] - pt_coords[0])**2 + (coord[1] - pt_coords[1])**2)**0.5
                if dist < min_dist:
                    min_dist = dist
                    nearest_node = node_id
                    
            snapped.append(nearest_node)
        return snapped

    # Snap all origin, TSS, and CCH points to their nearest network nodes.
    origins_gdf['nearest_node'] = snap_to_network(origins_gdf)
    tss_gdf['nearest_node'] = snap_to_network(tss_gdf)
    cch_gdf['nearest_node'] = snap_to_network(cch_gdf)

    # Filter out any points that could not be snapped to the network (e.g., if nearest_node is None).
    origins_gdf_snapped = origins_gdf.dropna(subset=['nearest_node'])
    tss_gdf_snapped = tss_gdf.dropna(subset=['nearest_node'])
    cch_gdf_snapped = cch_gdf.dropna(subset=['nearest_node'])

    print(f"  → Snapped {origins_gdf_snapped['nearest_node'].count()} origins, "
          f"{tss_gdf_snapped['nearest_node'].count()} TSS, "
          f"{cch_gdf_snapped['nearest_node'].count()} CCH to the network.")

    # Create dictionaries to map facility names/locations to their snapped node IDs.
    origin_map = dict(zip(origins_gdf_snapped[required_cols['origins']], origins_gdf_snapped['nearest_node']))
    tss_map = dict(zip(tss_gdf_snapped[required_cols['tss']], tss_gdf_snapped['nearest_node']))
    cch_map = dict(zip(cch_gdf_snapped[required_cols['cch']], cch_gdf_snapped['nearest_node']))
    
    # Create reverse lookup dictionaries to get facility names from node IDs.
    node_to_tss = {v: k for k, v in tss_map.items()}
    node_to_cch = {v: k for k, v in cch_map.items()}
    
    # Create sets of node IDs for TSS and CCH locations for efficient lookup.
    tss_nodes = list(tss_map.values())
    cch_nodes = list(cch_map.values())
    all_dests = list(set(tss_nodes) | set(cch_nodes)) # Combined set of all possible destination nodes.
except Exception as e:
    print(f"ERROR: Failed to setup routing: {e}")
    exit()

print("\n--- 5. Routing ---")
print("Computing shortest paths...")
results = []

# Ensure there are valid destinations (TSS or CCH) to route to.
if not all_dests:
    print("ERROR: No valid destinations (TSS or CCH) found!")
    exit()

print(f"Processing {len(origin_map)} origins to {len(tss_nodes)} TSS and {len(cch_nodes)} CCH locations")

routing_start_time = time.time()
total_origins = len(origin_map)

# Iterate through each origin to find the optimal route to a CCH.
for i, (orig_name, orig_node) in enumerate(origin_map.items(), 1):
    # Print progress updates periodically.
    if i % max(1, total_origins // 20) == 0: 
        elapsed = time.time() - routing_start_time
        estimated_total = (elapsed / i) * total_origins
        remaining = estimated_total - elapsed
        print(f"Processing origin {i}/{total_origins} ({(i/total_origins)*100:.1f}%)... "
              f"Elapsed: {elapsed:.1f}s, Remaining: {remaining:.1f}s")
    
    try:
        # Calculate shortest distances from the current origin to all other nodes in the graph.
        dists = g.distances(source=orig_node, target=None, weights='length')[0]
        
        # Initialize a dictionary to store the best direct CCH route details.
        best_cch_route = {
            'total_dist': float('inf'),
            'leg1_type': None, 'leg1_loc': None, 'leg1_dist': None,
            'leg2_type': None, 'leg2_loc': None, 'leg2_dist': None,
            'route_nodes': None, 'status': 'NoReach_TSS_CCH'
        }

        # --- Option 1: Direct path to CCH ---
        # Find all reachable CCHs from the origin and their distances.
        cch_dists = [(n, dists[n]) for n in cch_nodes if not np.isinf(dists[n])]
        if cch_dists:
            # Select the closest CCH and its distance.
            min_cch_node, min_cch_dist = min(cch_dists, key=lambda x: x[1])
            # Get the actual shortest path (sequence of nodes) to this closest CCH.
            path_direct_cch = g.get_shortest_paths(orig_node, to=min_cch_node, weights='length')[0]
            
            # Update best_cch_route with direct CCH details.
            best_cch_route.update({
                'total_dist': min_cch_dist,
                'leg1_type': CCH_TYPE_VALUE,
                'leg1_loc': node_to_cch[min_cch_node],
                'leg1_dist': min_cch_dist,
                'route_nodes': path_direct_cch,
                'status': 'Direct_CCH'
            })

        # --- Option 2: Path via TSS to CCH ---
        # Initialize variables to track the best route via a TSS.
        best_tss_cch_dist = float('inf')
        best_tss_cch_path1 = None
        best_tss_cch_path2 = None
        best_tss_node_found = None
        best_cch_node_via_tss = None

        # Find all reachable TSS from the origin and their distances.
        reachable_tss_from_orig = [(n, dists[n]) for n in tss_nodes if not np.isinf(dists[n])]

        if reachable_tss_from_orig:
            # Sort reachable TSS by distance from origin to prioritize closer TSS for efficiency.
            reachable_tss_from_orig.sort(key=lambda x: x[1])

            # Iterate through each reachable TSS to find the best two-leg route (origin -> TSS -> CCH).
            for tss_node, tss_dist_from_orig in reachable_tss_from_orig:
                # Get the path from origin to the current TSS.
                path_orig_to_tss = g.get_shortest_paths(orig_node, to=tss_node, weights='length')[0]
                if not path_orig_to_tss: 
                    continue

                # Calculate distances from the current TSS to all other nodes.
                dists_from_tss = g.distances(source=tss_node, target=None, weights='length')[0]
                # Find all reachable CCHs from the current TSS and their distances.
                reachable_cch_from_tss = [(n, dists_from_tss[n]) for n in cch_nodes if not np.isinf(dists_from_tss[n])]

                if reachable_cch_from_tss:
                    # Select the closest CCH from the current TSS.
                    min_cch_node_from_tss, min_cch_dist_from_tss = min(reachable_cch_from_tss, key=lambda x: x[1])
                    # Get the path from the current TSS to its closest CCH.
                    path_tss_to_cch = g.get_shortest_paths(tss_node, to=min_cch_node_from_tss, weights='length')[0]
                    
                    if path_tss_to_cch: 
                        # Calculate the total distance for this two-leg route.
                        current_total_dist = tss_dist_from_orig + min_cch_dist_from_tss
                        # If this route is better than the previously found best TSS-CCH route, update.
                        if current_total_dist < best_tss_cch_dist:
                            best_tss_cch_dist = current_total_dist
                            best_tss_cch_path1 = path_orig_to_tss
                            best_tss_cch_path2 = path_tss_to_cch
                            best_tss_node_found = tss_node
                            best_cch_node_via_tss = min_cch_node_from_tss

        # --- Determine closest CCH and TSS from origin ---
        # Find the closest direct CCH from the origin.
        min_cch_dist_from_orig = float('inf')
        closest_cch_node_from_orig = None
        if cch_nodes:
            cch_dists_from_orig = [(n, dists[n]) for n in cch_nodes if not np.isinf(dists[n])]
            if cch_dists_from_orig:
                closest_cch_node_from_orig, min_cch_dist_from_orig = min(cch_dists_from_orig, key=lambda x: x[1])

        # Find the closest TSS from the origin.
        min_tss_dist_from_orig = float('inf')
        closest_tss_node_from_orig = None
        if tss_nodes:
            tss_dists_from_orig = [(n, dists[n]) for n in tss_nodes if not np.isinf(dists[n])]
            if tss_dists_from_orig:
                closest_tss_node_from_orig, min_tss_dist_from_orig = min(tss_dists_from_orig, key=lambda x: x[1])

        chosen_route = None

        # --- Decision Logic: Choose based on closest facility type ---
        # Compare the shortest direct CCH distance with the shortest distance to any TSS.
        if min_cch_dist_from_orig <= min_tss_dist_from_orig:
            # If a direct CCH is closer or equally close, choose the direct CCH route.
            if closest_cch_node_from_orig is not None:
                path_direct_cch = g.get_shortest_paths(orig_node, to=closest_cch_node_from_orig, weights='length')[0]
                chosen_route = {
                    'origin': orig_name,
                    'leg1_type': CCH_TYPE_VALUE,
                    'leg1_loc': node_to_cch[closest_cch_node_from_orig],
                    'leg1_dist': min_cch_dist_from_orig,
                    'leg2_type': None,
                    'leg2_loc': None,
                    'leg2_dist': 0.0,
                    'total_dist': min_cch_dist_from_orig,
                    'route_nodes': path_direct_cch,
                    'status': 'Direct_CCH'
                }
        else: 
            # If the closest TSS is strictly closer than any direct CCH, consider routing via TSS.
            if closest_tss_node_from_orig is not None:
                # Calculate path from this closest TSS to its closest CCH.
                dists_from_closest_tss = g.distances(source=closest_tss_node_from_orig, target=None, weights='length')[0]
                cch_dists_from_closest_tss = [(n, dists_from_closest_tss[n]) for n in cch_nodes if not np.isinf(dists_from_closest_tss[n])]

                if cch_dists_from_closest_tss:
                    min_cch_node_from_tss, min_cch_dist_from_tss = min(cch_dists_from_closest_tss, key=lambda x: x[1])
                    path_orig_to_tss = g.get_shortest_paths(orig_node, to=closest_tss_node_from_orig, weights='length')[0]
                    path_tss_to_cch = g.get_shortest_paths(closest_tss_node_from_orig, to=min_cch_node_from_tss, weights='length')[0]

                    if path_orig_to_tss and path_tss_to_cch:
                        # Concatenate the two path segments (origin->TSS and TSS->CCH).
                        full_path = []
                        full_path.extend(path_orig_to_tss)
                        # Avoid duplicating the intermediate TSS node if it's the last node of the first path and first of the second.
                        if path_orig_to_tss[-1] == path_tss_to_cch[0]:
                            full_path.extend(path_tss_to_cch[1:])
                        else:
                            full_path.extend(path_tss_to_cch)

                        chosen_route = {
                            'origin': orig_name,
                            'leg1_type': TSS_TYPE_VALUE,
                            'leg1_loc': node_to_tss[closest_tss_node_from_orig],
                            'leg1_dist': min_tss_dist_from_orig,
                            'leg2_type': CCH_TYPE_VALUE,
                            'leg2_loc': node_to_cch[min_cch_node_from_tss],
                            'leg2_dist': min_cch_dist_from_tss,
                            'total_dist': min_tss_dist_from_orig + min_cch_dist_from_tss,
                            'route_nodes': full_path,
                            'status': 'Via_TSS'
                        }
        
        # If no valid route (direct CCH or via TSS) was found for the origin.
        if chosen_route is None:
            chosen_route = {
                'origin': orig_name,
                'status': 'NoReach_TSS_CCH', # Status indicating no reachable CCH or TSS.
                'leg1_type': None, 'leg1_loc': None, 'leg1_dist': None,
                'leg2_type': None, 'leg2_loc': None, 'leg2_dist': None,
                'total_dist': None, 'route_nodes': None
            }
        
        results.append(chosen_route)
            
    except Exception as e:
        # Handle any errors during routing for a specific origin.
        print(f"Error routing for {orig_name}: {str(e)}")
        results.append({
            'origin': orig_name,
            'status': 'Error', # Status indicating an error occurred.
            'leg1_type': None,
            'leg1_loc': None,
            'leg1_dist': None,
            'leg2_type': None,
            'leg2_loc': None,
            'leg2_dist': None,
            'total_dist': None,
            'route_nodes': None
        })

print("--- 7 & 8: Save results ---")
# Check if any results were generated before attempting to save.
if not results:
    print("WARNING: No routing results to save!")
    exit()

# Convert the list of routing results into a Pandas DataFrame.
df = pd.DataFrame(results)
print(f"Processing {len(df)} routing results...")

# Merge tonnage information from the origins GeoDataFrame if the 'tonnage' column exists.
if 'tonnage' in origins_gdf.columns:
    ton = origins_gdf[[required_cols['origins'],'tonnage']]
    df = df.merge(ton, how='left', left_on='origin', right_on=required_cols['origins'])
    # Drop the redundant origin column if it's different from the merged 'origin' column.
    if required_cols['origins'] in df.columns and required_cols['origins']!='origin': 
        df.drop(columns=[required_cols['origins']], inplace=True)
else:
    df['tonnage']=None # If no tonnage column, add it with None values.

def geom_from_nodes(nlist):
    """
    Generates a LineString geometry from a list of node IDs.
    Uses the node_to_coords map to retrieve coordinates for each node.
    """
    if not nlist or len(nlist) < 2:
        return None
    pts = [node_to_coords[n] for n in nlist]
    return LineString(pts)

# Apply the geom_from_nodes function to create geometry objects for each route.
df['geometry']=df['route_nodes'].apply(geom_from_nodes)
# Convert the DataFrame to a GeoDataFrame, setting the geometry column and CRS.
res_gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=TARGET_CRS)
# Save the results to a CSV file (without index) and a GeoPackage file.
df.to_csv(RESULTS_CSV, index=False, float_format='%.3f')
res_gdf.to_file(RESULTS_GPKG, driver='GPKG')
print("Done.")
