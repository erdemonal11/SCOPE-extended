# SCOPE-Extended

SCOPE-Extended is a literature review tool for detecting chemical patterns across scientific publications, enhanced with publication date analytics. It visualizes chemicals resulting from a query search in an interactive density plot (mass over logP), while providing temporal analysis of the literature.

## Features

- Search Europe PMC database with custom queries
- Extract ChEBI chemical identifiers from publications
- Visualize chemical compounds in interactive heatmap plots
- Analyze distribution of chemical properties (mass, logP)
- **NEW**: Track and visualize publication dates of papers containing specific chemicals
- Compare chemical distributions across different search queries
- Interactive controls for blur, saturation, and visualization modes

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/SCOPE-Extended.git
cd SCOPE-Extended
```

2. Download the required ChEBI property files:
```bash
python download_files.py
```

## Usage

### 1. Create a query file
Create a text file with search queries in the following format:
```
<output-tag-1>, <query-1>
<output-tag-2>, <query-2>
```

For example:
```
CE-binding, "capillary electrophoresis" AND ("protein-ligand binding" OR "binding constant" OR "affinity measurement")
```
The search queries should follow the syntax used on the Europe PMC site.

### 2. Search for publications
```bash
python search_query.py -i <path-to-query-file>
```
This will retrieve publications matching your queries, extract ChEBI annotations, and save the results in the results folder.

**Note**: For large queries with many hits, this process may take a considerable amount of time. Consider starting with more specific queries when testing.

### 3. Create data tables
```bash
python make_table.py -i results -t folder
```
Or for a specific result file:
```bash
python make_table.py -i <path-to-result-file> -t file
```
This will generate tables containing chemical properties and publication year information.

### 4. Visualize the results
```bash
python visualize_multiplot.py -i tables -o <plot-name>
```

Optional parameters:
- `-xmin`: Minimum logP value (default: -5)
- `-xmax`: Maximum logP value (default: 10)
- `-ymax`: Maximum mass in Da (default: 1600)
- `-class`: ChEBI class ID to highlight

## Interactive Visualization

The visualization includes several interactive features:

- Term selection: Switch between different search terms
- Blur control: Adjust Gaussian blur for density visualization
- Saturation: Modify the intensity scaling
- Color schemes: Toggle between Viridis and grayscale
- TFIDF: Toggle between count and TFIDF normalization
- Statistics button: View detailed statistics including publication year analysis
- Metadata button: View query metadata

When hovering over a hexagon in the plot, a tooltip displays the top ChEBI compounds in that region.

## Publication Date Analysis

The new publication date analysis feature shows:

- Year range (earliest to latest publication)
- Mean and median publication years
- Distribution of publications per year in tabular format

This temporal information helps identify trends in chemical research literature.

## Directories

- `files/`: ChEBI property files (downloaded by download_files.py)
- `results/`: ChEBI IDs extracted from publications
- `metadata/`: Query metadata
- `tables/`: Processed data tables ready for visualization
- `plots/`: Generated visualization HTML files
- `publication_dates/`: Publication year information

## Requirements

- Python 3.6+
- Required Python packages:
  - pandas
  - numpy
  - bokeh
  - requests
  - tqdm

## Troubleshooting

If you encounter connection issues when querying the Europe PMC API:

- Check your internet connection
- Verify that the API is accessible from your location
- Try using a more specific query to reduce the number of results
- Increase the timeout settings in the script if necessary
- The script includes automatic retries for failed connections 