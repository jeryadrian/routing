# Transport Network Optimization for Temporary Storage Sites

This project implements a network-based optimization algorithm to find optimal locations for temporary storage sites (TSS) between material origin points and circular construction hubs (CCH). The primary goal is to maximize the tonnage of materials routed through TSS, while also minimizing total transport effort measured in tonne-kilometers and considering road network constraints.

## Project Structure

```markdown
.
├── greedy.py             # Main optimization algorithm
├── tkm_calculator.ipynb  # Analysis notebook for transport-kilometer calculations
├── requirements.txt      # Python package dependencies
├── input/                # Input data files
│   ├── project.gpkg      # Material origin points with tonnage
│   ├── cch.gpkg          # Circular Construction Hub locations
│   ├── road.gpkg         # Road network for routing
│   └── grid.gpkg         # Grid of potential TSS locations
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

Once the graph is constructed, shortest path distances between any two points on the network (e.g., an origin and a CCH, or an origin and a TSS) are calculated using Dijkstra's algorithm, which is efficiently implemented in `igraph`. This allows the optimization algorithm to accurately determine the transport effort for various routing scenarios.

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

The optimization produces several output files:

- `optimal_tss_locations.gpkg`: Selected TSS locations
- `final_origin_routes.gpkg`: Detailed routing results
- `final_origin_routes.csv`: Route statistics and metrics
- Various debug files for visualization and validation

## Analysis Tools

The `tkm_calculator.ipynb` notebook provides:

- Detailed transport-kilometer calculations
- Route statistics and summaries
- Data validation and verification
- Visual analysis of results

## Algorithm Overview

The greedy algorithm in `greedy.py`:

1. Loads and validates input data
2. Builds a network graph from road data
3. Snaps all points (origins, CCHs, grid) to network
4. Precomputes shortest path distances
5. Iteratively selects optimal TSS locations
6. Generates detailed routing solutions
7. Saves results and debug information
