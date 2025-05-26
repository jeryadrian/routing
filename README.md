# Transport Network Optimization for Temporary Storage Sites

This project implements a [greedy algorithm](https://en.wikipedia.org/wiki/Greedy_algorithm) to find optimal locations for temporary storage sites (TSS) between material origin points and circular construction hubs (CCH). The primary goal is to maximize the tonnage of materials routed through TSS, while also minimizing total transport effort measured in tonne-kilometers and considering road network constraints.

## Project Structure

```markdown
.
├── greedy.py             # Script for optimizing TSS locations
├── routing.py            # Script for direct routing calculations
├── tkm_calculator.ipynb  # Analysis notebook for t-km calculations
├── requirements.txt      # Python package dependencies
├── input/                # Input data files
│   ├── project.gpkg      # Material origin points with tonnage
│   ├── cch.gpkg          # Circular Construction Hub locations
│   ├── road.gpkg         # Road network for routing
│   ├── grid.gpkg         # Grid of potential TSS locations (used by greedy.py)
│   ├── tss.gpkg          # Pre-defined TSS locations (used by routing.py)
│   ├── cch.gpkg          # Pre-defined CCH locations (used by routing.py)
└── output/               # Generated results and debug files
```

## Key Features

- Network-based optimization using Python's [igraph](https://igraph.org/) library
- Support for both direct routes (Origin → CCH) and indirect routes (Origin → TSS → CCH)
- Automatic snapping of points to road network nodes
- Detailed distance and tonne-kilometer calculations
- Comprehensive debug outputs and visualizations

## Routing Mechanism

The core of the routing functionality relies on the `igraph` library, a powerful tool for creating and manipulating graphs. The road network data (`road.gpkg`) is transformed into a graph where road intersections and significant points become 'nodes' (vertices) and road segments become 'edges'. Each edge is assigned a 'length' attribute, representing the real-world distance of that road segment.

Once the graph is constructed, shortest path distances between any two points on the network (e.g., an origin and a CCH, or an origin and a TSS) are calculated using [Dijkstra's algorithm](https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm), which is efficiently implemented in `igraph`. These calculated distances are then used by the greedy optimization algorithm to iteratively select the best TSS locations, prioritizing those that maximize tonnage routed through them while minimizing overall transport effort.

## Script Usage and Purpose

This project includes two main Python scripts, `greedy.py` and `routing.py`, each serving a distinct purpose based on the input data they utilize:

### `greedy.py` (TSS Location Optimization)

- **Purpose:** This script is designed to *optimize the placement* of a specified number of Temporary Storage Sites (TSS). It identifies the best locations for these sites from a grid of potential candidates to minimize the total transport effort (tonne-kilometers) between material origins and Circular Construction Hubs (CCH).
- **Input Requirements:**
  - `project.gpkg`: Material origin points with tonnage.
  - `cch.gpkg`: Circular Construction Hub locations.
  - `road.gpkg`: The road network.
  - `grid.gpkg`: A grid of *potential* TSS locations from which the optimal sites will be chosen.
- **Output:** Generates optimal TSS locations, detailed final routes based on these optimized sites, and an optimization log.

### `routing.py` (Direct Routing Calculation)

- **Purpose:** This script is for *calculating shortest paths* from origins to Circular Construction Hubs (CCH), with the option to route via *pre-defined* Temporary Storage Sites (TSS). It performs direct routing calculations given all facility locations are already known.
- **Input Requirements:**
  - `input/project.gpkg`: Material origin points.
  - `input/tss.gpkg`: Specific, *pre-defined* Temporary Storage Site (TSS) locations.
  - `input/cch.gpkg`: Specific, *pre-defined* Circular Construction Hub (CCH) locations.
  - `input/road.gpkg`: The road network.
- **Output:** Generates routing results (paths and distances) for each origin to its chosen destination (direct CCH or via a pre-defined TSS).

## Detailed Calculation Differences

While both `greedy.py` and `routing.py` perform routing calculations, their core logic and objectives differ significantly, especially in how they handle Temporary Storage Sites (TSS) and determine optimal routes.

### 1. Primary Objective

- **`greedy.py` (Optimization):** Its main goal is to *find and select* the optimal TSS locations from a set of candidates. It's an optimization script that determines *where* TSS should be placed to minimize overall transport effort and maximize tonnage routed through them.
- **`routing.py` (Calculation):** Its main goal is to *calculate routes* and associated tonne-kilometers using *pre-defined* TSS locations. It assumes TSS locations are already given and focuses on finding the best path through them or directly to CCHs.

### 2. TSS Handling

- **`greedy.py`:** Implements a **greedy heuristic** to iteratively choose a specified number of TSS from a grid of potential sites. The selection is based on which candidate site provides the most benefit (e.g., highest tonnage routed through TSS, lowest overall tonne-km).
- **`routing.py`:** Does **not perform any TSS optimization or selection**. It simply uses the TSS locations provided in its input file (`input/tss.gpkg`).

### 3. Route Decision Logic for an Origin

- **`greedy.py`:** For each origin, it evaluates both a direct route to the closest CCH and a two-leg route via any *active* TSS (origin → active TSS → closest CCH from that TSS). It then selects the route that results in the **absolute minimum total distance** for that origin.
- **`routing.py`:** For each origin, it first identifies the *closest direct CCH* and the *closest pre-defined TSS*. The decision is then made based on which of these two facilities is *closer in the first leg*:
  - If the closest direct CCH is closer than or equally close to the closest TSS, it routes directly to that CCH.
  - If the closest TSS is strictly closer than any direct CCH, it routes via that closest TSS to its own closest CCH.
    This "first-leg proximity" rule can lead to different route choices compared to `greedy.py`'s "overall minimum total distance" rule.

### 4. Precomputation of Distances

- **`greedy.py`:** Precomputes and stores all relevant shortest path distances (origin to candidate TSS, origin to CCH, candidate TSS to CCH) in dictionaries before the main optimization loop. This speeds up the iterative evaluation of candidate TSS.
- **`routing.py`:** Calculates shortest path distances on-the-fly for each origin and for each TSS-to-CCH segment within its main routing loop. It does not rely on extensive precomputed global distance matrices.

## Requirements

- Python 3.8+
- Required packages listed in requirements.txt:
  - geopandas>=0.13.0
  - pandas>=2.0.0
  - numpy>=1.24.0
  - shapely>=2.0.0
  - igraph>=0.10.0

## Usage

1. Place input files in the `input/` directory:
   - Material origin points (GeoPackage with tonnage data)
   - CCH locations (GeoPackage)
   - Road network (GeoPackage)
   - Potential TSS locations grid (GeoPackage)

2. Run the optimization:

   ```bash
   python greedy.py
   ```

3. Analyze results using the Jupyter notebook:

   ```bash
   jupyter notebook tkm_calculator.ipynb
   ```

## Input Data Format

### project.gpkg

- Point geometries of material origins
- Required columns:
  - `name`: Unique identifier for each origin
  - `tonnage`: Material quantity in tonnes

### cch.gpkg

- Point geometries of circular construction hubs
- Required columns:
  - `loc`: Unique identifier for each CCH

### road.gpkg

- LineString geometries of the road network
- Used for routing and distance calculations

### grid.gpkg

- Point geometries of potential TSS locations
- Forms the search space for optimization

## Output Files

This project generates various output files depending on the script executed:

### From `greedy.py` (Optimization)

- `optimal_tss_locations.gpkg`: Selected optimal TSS locations.
- `final_origin_routes.gpkg`: Detailed routing results for origins based on optimized TSS.
- `final_origin_routes.csv`: Route statistics and metrics in CSV format.
- Various debug files (e.g., `candidate_tss_nodes_debug.gpkg`, `graph_nodes_opt_debug.gpkg`) for visualization and validation of intermediate steps.

### From `routing.py` (Direct Routing)

- `routing_results.csv`: Routing results (paths and distances) for each origin.
- `routing_results.gpkg`: GeoPackage containing the geometries of the calculated routes.

## Analysis Tools

The `tkm_calculator.ipynb` Jupyter notebook is provided for comprehensive analysis of routing results. It can process the output CSV files from both `greedy.py` (`final_origin_routes.csv`) and `routing.py` (`routing_results.csv`). The notebook provides:

- Detailed transport-kilometer calculations.
- Route statistics and summaries by CCH and TSS locations.
- Data validation and verification of calculated distances and tonne-kilometers.
- Visual analysis of results (if plotting capabilities are utilized within the notebook).

## Algorithm Overview

The greedy algorithm in `greedy.py`:

1. Loads and validates input data
2. Builds a network graph from road data
3. Snaps all points (origins, CCHs, grid) to network
4. Precomputes shortest path distances
5. Iteratively selects optimal TSS locations
6. Generates detailed routing solutions
7. Saves results and debug information
