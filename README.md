# SCOPE-Extended

A comprehensive text mining tool for analyzing chemical patterns in scientific literature with temporal analysis

SCOPE-Extended is an enhanced version of the original SCOPE repository (Search and Chemical Ontology Plotting Environment), specifically optimized for protein-ligand binding studies. This tool was developed for and utilized in the research paper "A Text Mining Approach on Protein–Ligand Binding Studies Based on Analytical Techniques".

## Overview

SCOPE-Extended combines natural language processing with chemical informatics to extract and visualize chemical patterns across scientific publications. It leverages Europe PMC's database to detect ChEBI chemical identifiers and presents them in interactive density plots, plotting mass against logP values. The tool has been enhanced with publication date analytics to provide temporal insights into chemical research trends.

## Key Enhancements Over Original SCOPE

- **Temporal Analysis**: Track publication dates of papers containing specific chemicals
- **Improved Performance**: Optimized query processing for faster results
- **Enhanced Visualizations**: More interactive control options and better rendering
- **Specialized Filters**: Tailored for protein-ligand binding studies
- **Stability Improvements**: Fixed bugs from the original repository
- **User Experience**: Streamlined interface and improved documentation
- **Drug-Likeness Analysis**: New features for analyzing drug-like properties
- **Extended Chemical Properties**: Additional molecular property calculations
- **Smart Protein Family Prediction**: AI-based protein family assignment
- **Enhanced Chemical Type Classification**: Improved chemical type determination

## Features

- Literature Mining: Search Europe PMC database with custom queries
- Chemical Entity Recognition: Extract ChEBI chemical identifiers from publications
- Property Visualization: Interactive heatmap plots of chemical compounds
- Chemical Distribution Analysis: Analyze mass and logP distributions
- Publication Timeline: Track and visualize publication dates for chemical research
- Comparative Analysis: Compare chemical distributions across different search queries
- Interactive Controls: Adjust blur, saturation, and visualization modes
- Publication Metrics: Statistics on research trends over time
- Drug-Likeness Rules: Analyze compounds against Lipinski's Rule of Five and Veber's Rule
- Extended Properties: Calculate H-bond donors/acceptors, TPSA, and rotatable bonds
- Protein Family Prediction: Smart assignment of protein families based on techniques
- Chemical Type Classification: Enhanced classification of compounds into ligands, buffers, and solvents

## Installation

1. Clone the repository:
```bash
git clone https://github.com/erdemonal11/SCOPE-extended.git
cd SCOPE-Extended
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Download the required ChEBI property files:
```bash
python download_files.py
```

## Workflow

### 1. Create a Query File
Create a text file with search queries in the following format:
```
<output-tag-1>, <query-1>
<output-tag-2>, <query-2>
```

Example query for protein-ligand binding studies:
```
CE-binding, "capillary electrophoresis" AND ("protein-ligand binding" OR "binding constant" OR "affinity measurement")
MS-binding, "mass spectrometry" AND ("protein-ligand binding" OR "binding constant" OR "affinity measurement")
SPR-binding, "surface plasmon resonance" AND ("protein-ligand binding" OR "binding constant" OR "affinity measurement")
```

The search queries should follow the Europe PMC syntax.

### 2. Search for Publications
```bash
python search_query.py -i <path-to-query-file>
```

This will:
- Retrieve publications matching your queries
- Extract ChEBI annotations from the publications
- Save the results in the results folder

**Note**: For large queries, this process may take time. Consider starting with specific queries when testing.

### 3. Generate Data Tables
```bash
python make_table.py -i results -t folder
```

Or for a specific result file:
```bash
python make_table.py -i <path-to-result-file> -t file
```

This creates tables containing chemical properties and publication year information.

### 4. Visualize the Results
```bash
python visualize_multiplot.py -i tables -o <plot-name>
```

Optional parameters:
- `-xmin`: Minimum logP value (default: -5)
- `-xmax`: Maximum logP value (default: 10)
- `-ymax`: Maximum mass in Da (default: 1600)
- `-class`: ChEBI class ID to highlight
- `-years`: Filter by year range (e.g., "2010-2022")

### 5. Generate Additional Figures
```bash
python generate_figures.py
```

This will create additional visualizations including:
- Protein family distribution
- Chemical type analysis
- Molecular property correlations
- Drug-likeness analysis

## Interactive Visualization Features

The visualization interface includes:

- Term Selection: Switch between different search terms
- Blur Control: Adjust Gaussian blur for density visualization
- Saturation: Modify the intensity scaling
- Color Schemes: Toggle between Viridis and grayscale
- Normalization: Toggle between count and TFIDF normalization
- Statistics Panel: View detailed statistics including publication year analysis
- Metadata Panel: View query metadata
- Chemical Information: Hover tooltips showing top ChEBI compounds in each region
- Export Options: Save visualizations as PNG or interactive HTML files
- Drug-Likeness Highlighting: Visualize compounds that follow drug-likeness rules
- Extended Property Analysis: View additional molecular properties
- Protein Family Distribution: Analyze distribution across protein families
- Chemical Type Analysis: Compare different chemical types

## Publication Date Analysis

The temporal analysis feature reveals:

- Year Range: Earliest to latest publication containing specific chemicals
- Central Tendency: Mean and median publication years
- Temporal Distribution: Visualize publications per year
- Trend Analysis: Identify emerging or declining research topics

## Directory Structure

- `files/`: ChEBI property files
- `results/`: ChEBI IDs extracted from publications
- `metadata/`: Query metadata
- `tables/`: Processed data tables
- `plots/`: Generated visualization HTML files
- `publication_dates/`: Publication year information
- `examples/`: Example query files and results
- `docs/`: Additional documentation

## Requirements

- Python 3.6+
- Required packages:
  - pandas
  - numpy
  - bokeh
  - requests
  - tqdm
  - scikit-learn
  - matplotlib
  - rdkit (for molecular property calculations)
  - seaborn (for enhanced visualizations)

## Application in Research

SCOPE-Extended was successfully applied in the paper "A Text Mining Approach on Protein–Ligand Binding Studies Based on Analytical Techniques" to:

- Identify common chemicals used in protein-ligand binding studies
- Track the evolution of analytical techniques over time
- Compare chemical distributions across different analytical methods
- Discover trends in ligand properties based on analytical technique preference
- Analyze drug-likeness patterns in protein-ligand binding studies
- Study protein family preferences for different chemical types
- Investigate temporal trends in chemical property distributions

## Troubleshooting

If you encounter issues:

- Connection Problems: Check your internet connection and Europe PMC API accessibility
- Memory Errors: Use more specific queries to reduce result size
- Visualization Issues: Update bokeh to the latest version
- Missing Data: Verify that the ChEBI property files were downloaded correctly
- Timeout Errors: Increase the timeout settings in the script
- RDKit Errors: Ensure RDKit is properly installed for molecular property calculations
- Property Calculation Issues: Check SMILES string validity for problematic compounds

## Acknowledgments

- Original SCOPE repository by [ReinV](https://github.com/ReinV/SCOPE)
- Europe PMC for providing the API
- ChEBI for the chemical ontology database
- RDKit for molecular property calculations

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 