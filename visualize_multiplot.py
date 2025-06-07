#!/usr/bin/python
import time
import numpy as np
import pandas as pd
import argparse
from math import exp
from math import sqrt
from datetime import datetime
import os
from collections import Counter

# BOKEH
from bokeh import events
from bokeh.io import output_file, show
from bokeh.models import CustomJS, HoverTool, ColumnDataSource, Slider, CheckboxGroup, RadioGroup, Button, MultiSelect, Div
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.transform import log_cmap
from bokeh.util.hex import axial_to_cartesian
from bokeh.util.hex import cartesian_to_axial
from bokeh.layouts import column, row
from bokeh.palettes import Viridis256, Greys256

def import_table(file):
    '''
    This function imports the pkl files from the tables folder, and returns a pandas dataframe.
    '''
    table = pd.read_pickle(file)
    return table

def create_array(table):
    '''
    This function receives a dataframe with logP and Mass values for every ChEBI identifier.
    It returns two numpy arrays: for mass and logP.
    Fixed to ensure both arrays have the same length.
    '''
    # Filter for rows that have valid (non-NaN) values for both Mass and logP
    valid_data = table.dropna(subset=['Mass', 'logP'])
    
    if len(valid_data) == 0:
        print(f"Warning: No compounds with both Mass and logP data found!")
        return np.array([]), np.array([])
    
    # Create lists from the filtered data
    x = [float(logP) for logP in valid_data.logP]
    y = [float(mass) for mass in valid_data.Mass]
    
    print(f"Creating arrays for {len(x)} compounds with complete Mass/logP data")
    
    return np.asarray(x), np.asarray(y)

def hexbin(df, x, y, size, aspect_scale, orientation):
    '''
    This function receives x and y coordinate arrays and converts these into q and r hexagon coordinates by calling Bokeh's "cartesian_to_axial" function.
    The q and r coordinates are added to the dataframe, and this dataframe is returned.
    '''
    q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=aspect_scale)

    df.loc[:,'q'] = q
    df.loc[:,'r'] = r

    return df

def add_tooltip_columns(df, table):
    '''
    Enhanced tooltip with drug-likeness information
    '''
    table = table.drop(['Class', 'logP', 'Mass'], axis=1, errors='ignore')
    table = table.reset_index()

    # Define tooltip size
    TOOLTIP_COUNT = 3
    columns = table.columns
    tooltip_columns = {column+str(i):[] for i in range(1, TOOLTIP_COUNT+1) for column in columns}

    # Extract ChEBI identifiers from dataframe
    chebi_ids = [ids if isinstance(ids, list) else [] for ids in df.ChEBI]
    list_for_df = []

    # Use chebi identifiers to look up information in the "table" dataframe
    for ids in chebi_ids:
        rows = table[table.ChEBI.isin(ids)]

        # Sort and select most frequent ChEBI identifiers
        rows = rows.sort_values(by='Count', ascending=False)
        values_nested = rows[0:TOOLTIP_COUNT].values.tolist()

        # Unnest information from the table,
        values_unnested = [item for sublist in values_nested for item in sublist]
        while len(values_unnested) < (TOOLTIP_COUNT*len(columns)):
            values_unnested.append("-")

        list_for_df.append(values_unnested)

    df_tooltip = pd.DataFrame(list_for_df, columns = tooltip_columns)

    df = df.join(df_tooltip, how='left')

    return df

def get_blur(x,y,sigma_x,sigma_y):
    '''
    This function receives x, y values and sigma x, sigma y values and returns the calculated blur value.
    See https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    '''
    return exp(-0.5*(x*x/sigma_x/sigma_x + y*y/sigma_y/sigma_y))

def get_rows(q, r, counts, tfidf, kernel, blur_max, step_size):
    '''
    For every hexagon (=row in the dataframe), counts will be distributed to surrounding hexagons in a gaussian manner (blur).
    This function receives information for one hexagon so that it can distribute its counts to other hexagons.
    To distribute the counts to other hexagons/rows, we use the "kernel".
    We create new rows for other hexagons, which means we might create many rows with the same hexagonal coordinates.
    Newly created rows will be returned
    '''
    # initiate original row
    rows = [[q, r, counts, tfidf] + [counts for i in np.arange(0, blur_max+step_size, step_size)] + [tfidf for i in np.arange(0, blur_max+step_size, step_size)]]

    # use kernel to distribute counts to other rows
    for coords_new, blur_factors in kernel.items():
        q_new = q + coords_new[0]
        r_new = r + coords_new[1]

        # create new row, use kernels coordinates and calculated blur values
        new_row = [q_new, r_new, 0, 0] + list(map(lambda n: n * counts, blur_factors)) + list(map(lambda n: n*tfidf, blur_factors))
        rows.append(new_row)
    return rows

def construct_kernel(blur_max, step_size):
    '''
    This function receives the maximum blur and step size to construct the kernel.
    The kernel is a dictionary that uses the coordinates as keys, and the blur values as values.
    The blur values depend on sd_x, so the kernel will return a list of blur values from 0 to blur_max.
    Blur values are calculated by "get_blur()"".
    The kernel will be used by "get_rows()" to distribute counts to surrounding hexagons.
    '''
    coordinates_to_distance = {(-5,2):(7.5, sqrt(3)/2),(-5,3):(7.5, sqrt(3)/2),
    (-4,1):(6, sqrt(3)),(-4,2):(6, 0),(-4,3):(6,sqrt(3)),
    (-3,0):(4.5, 3*sqrt(3)/2),(-3,1):(4.5, sqrt(3)/2),(-3,2):(4.5, sqrt(3)/2),(-3,3):(4.5, 3*sqrt(3)/2),
    (-2,-1):(3, 2*sqrt(3)),(-2,0):(3, sqrt(3)),(-2,1):(3, 0),(-2,2):(3, sqrt(3)),(-2,3):(3, 2*sqrt(3)),
    (-1,-1):(1.5, 3*sqrt(3)/2),(-1,0):(1.5, sqrt(3)/2),(-1,1):(1.5, sqrt(3)/2),(-1,2):(1.5, 3*sqrt(3)/2),
    (0,-2):(0, 2*sqrt(3)),(0,-1):(0, sqrt(3)),(0,1):(0, sqrt(3)),(0,2):(0, 2*sqrt(3)),
    (1,-2):(1.5, 3*sqrt(3)/2),(1,-1):(1.5, sqrt(3)/2),(1,0):(1.5, sqrt(3)/2),(1,1):(1.5, 3*sqrt(3)/2),(2,-3):(3, 2*sqrt(3)),
    (2,-2):(3, sqrt(3)),(2,-1):(3, 0),(2,0):(3, sqrt(3)),(2,1):(3, 2*sqrt(3)),
    (3,-3):(4.5, 3*sqrt(3)/2),(3,-2):(4.5, sqrt(3)/2),(3,-1):(4.5, sqrt(3)/2),(3,0):(4.5, 3*sqrt(3)/2),
    (4,-3):(6, sqrt(3)),(4,-2):(6, 0),(4,-1):(6, sqrt(3)),
    (5,-3):(7.5, sqrt(3)/2),(5,-2):(7.5, sqrt(3)/2)}

    kernel = {}
    for key, distance in coordinates_to_distance.items():
        kernel[key] = [0]
        for sd_x in np.arange(step_size, blur_max+step_size, step_size):
            kernel[key].append(get_blur(distance[0], distance[1], sd_x, sd_x/2))

    return kernel

def add_gaussian_blur(df, blur_max, step_size):
    '''
    Function:
    This function adds gaussian blur to the plot by going through all hexagons and applying the kernel to the neighbouring hexagons.
    To speed up the process, the pandas dataframe is put in a python dictionary, so to quickly find the neighbouring hexagon coordinates using the coordinates as keys.
    After bluring, the dictionary is put in a pandas dataframe again, and this dataframe is returned.

    Columns:
    sd_x values for calculating blur values are used to name a column with counts that result from blurring with that specific sd_x.
    This makes selecting correct sd_x column easy with the slider code, because the slider returns values that represent sd_x values.
    ColumnDataSource does not accept intergers as column names, so sd_x column names are changed to string.
    sd_x columns contain lists, in a list the first value is the normal counts, the second value is the tfidf count.
    '''

    kernel = construct_kernel(blur_max, step_size)

    columns = ['q', 'r', 'Count', 'TFIDF'] + [str(sd_x) for sd_x in np.arange(0, (blur_max+step_size), step_size)] + ['%s_tfidf' % str(sd_x) for sd_x in np.arange(0, (blur_max+step_size), step_size)]

    df_blur = pd.concat([pd.DataFrame(get_rows(q, r, counts, tfidf, kernel, blur_max, step_size),
        columns=columns) for q, r, counts, tfidf in zip(df.q, df.r, df.Count, df.TFIDF)], ignore_index=True)

    df_blur = df_blur.groupby(['q', 'r'], as_index=False).agg(sum)

    df_joined = df_blur.merge(df.loc[:,['q', 'r', 'ChEBI']], on=['q', 'r'], how='outer')

    return df_joined

def check_for_ids(id_list, class_id):
    '''
    This functions checks if the Class identifier is in the "id_list".
    This ChEBI ID-specific list contains other ChEBI identifiers of chemicals that are of higher level hierarchical class.
    '''
    check = any(str(class_id) == str(id) for id in id_list)
    return check

def analyze_publication_years(table):
    '''
    Analyzes publication years from the table and returns statistics.
    '''
    years_data = []
    
    # Extract all years from the PublicationYears column if it exists
    if 'PublicationYears' in table.columns:
        for years_list in table['PublicationYears']:
            if isinstance(years_list, list) and years_list:
                years_data.extend(years_list)
    
    if not years_data:
        return {
            'has_data': False,
            'min_year': 'N/A',
            'max_year': 'N/A',
            'mean_year': 'N/A',
            'median_year': 'N/A',
            'year_counts': {}
        }
    
    # Calculate statistics
    min_year = min(years_data)
    max_year = max(years_data)
    mean_year = np.mean(years_data)
    median_year = np.median(years_data)
    
    # Count publications per year
    year_counts = Counter(years_data)
    year_counts = {str(year): count for year, count in sorted(year_counts.items())}
    
    return {
        'has_data': True,
        'min_year': min_year,
        'max_year': max_year,
        'mean_year': round(mean_year, 2),
        'median_year': median_year,
        'year_counts': year_counts
    }

def analyze_drug_likeness(table):
    '''
    Analyze drug-likeness statistics from the table
    '''
    drug_stats = {}
    
    if 'Lipinski_RO5' in table.columns:
        ro5_compliant = table['Lipinski_RO5'].sum()
        ro5_total = table['Lipinski_RO5'].count()
        drug_stats['lipinski_compliant'] = ro5_compliant
        drug_stats['lipinski_percentage'] = round((ro5_compliant / ro5_total) * 100, 1) if ro5_total > 0 else 0
    
    if 'Veber_Rule' in table.columns:
        veber_compliant = table['Veber_Rule'].sum()
        veber_total = table['Veber_Rule'].count()
        drug_stats['veber_compliant'] = veber_compliant
        drug_stats['veber_percentage'] = round((veber_compliant / veber_total) * 100, 1) if veber_total > 0 else 0
    
    if 'Drug_Like' in table.columns:
        drug_like = table['Drug_Like'].sum()
        drug_total = table['Drug_Like'].count()
        drug_stats['drug_like'] = drug_like
        drug_stats['drug_like_percentage'] = round((drug_like / drug_total) * 100, 1) if drug_total > 0 else 0
    
    return drug_stats

def create_class_source(table, size, ratio, orientation, class_id):
    '''
    This function finds all chemicals belonging to class "class_id" as defined by the ChEBI ontology.
    It returns a dataframe filtered for these chemicals.
    '''
    if class_id == None:
        # if no class id is given, then class source should be empty
        df = pd.DataFrame().reindex_like(table)

    else:
        table_class_only = table[[check_for_ids(id_list, class_id) for id_list in table["Class"]]]

        # create array
        x, y = create_array(table_class_only)
        q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=ratio)
        df = table_class_only.reset_index()
        df.loc[:,"q"] = q
        df.loc[:,"r"] = r

        df = df.groupby(['q', 'r']).agg({'Count': 'sum', 'TFIDF': 'sum', 'ChEBI': list}).reset_index()
        df = df.drop(columns='ChEBI')
    # Create plot source
    source = ColumnDataSource(df)

    return source

def create_drug_like_source(table, size, ratio, orientation, rule_type='Drug_Like'):
    '''
    Create a source for highlighting drug-like compounds
    '''
    if rule_type not in table.columns:
        # Return empty source if rule column doesn't exist
        df = pd.DataFrame(columns=['q', 'r', 'Count', 'TFIDF'])
        return ColumnDataSource(df)
    
    # Filter for compounds that meet the drug-likeness rule
    drug_like_table = table[table[rule_type] == True]
    
    if len(drug_like_table) == 0:
        df = pd.DataFrame(columns=['q', 'r', 'Count', 'TFIDF'])
        return ColumnDataSource(df)
    
    # Create hexagonal coordinates
    x, y = create_array(drug_like_table)
    if len(x) == 0 or len(y) == 0:
        df = pd.DataFrame(columns=['q', 'r', 'Count', 'TFIDF'])
        return ColumnDataSource(df)
        
    q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=ratio)
    df = drug_like_table.reset_index()
    df.loc[:,"q"] = q
    df.loc[:,"r"] = r

    # Group by hexagon coordinates
    df = df.groupby(['q', 'r']).agg({'Count': 'sum', 'TFIDF': 'sum'}).reset_index()
    
    return ColumnDataSource(df)

def create_data_source(table, term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE):
    '''
    This function receives a query-specific "table" dataframe, and constructs a source for the hexagonal plot.
    Fixed to handle missing Mass/logP data properly.
    '''
    import pandas as pd
    import numpy as np
    from bokeh.models import ColumnDataSource
    from bokeh.util.hex import cartesian_to_axial
    
    # Filter table for compounds with both Mass and logP data
    table_filtered = table.dropna(subset=['Mass', 'logP']).copy()
    
    if len(table_filtered) == 0:
        print(f"Warning: No compounds with complete Mass/logP data for {term}")
        # Return empty source
        empty_df = pd.DataFrame(columns=['q', 'r', 'Count', 'TFIDF'])
        return ColumnDataSource(empty_df), f'No data available for {term}'
    
    print(f"Using {len(table_filtered)} compounds with complete Mass/logP data out of {len(table)} total")
    
    # create array with mass and logP values from filtered table
    x = [float(logP) for logP in table_filtered.logP]
    y = [float(mass) for mass in table_filtered.Mass]
    
    if len(x) == 0 or len(y) == 0:
        print(f"Warning: Empty arrays generated for {term}")
        empty_df = pd.DataFrame(columns=['q', 'r', 'Count', 'TFIDF'])
        return ColumnDataSource(empty_df), f'No plottable data for {term}'

    # create hexagonal coordinates from mass and logP values
    q, r = cartesian_to_axial(np.asarray(x), np.asarray(y), size, orientation=orientation, aspect_scale=ratio)

    # create dataframe with hexagonal coordinates as rows (row = hexagon)
    # IMPORTANT: Use the filtered table and reset index to match coordinate arrays
    df = table_filtered.reset_index(drop=False)  # Keep the original index as a column
    df = df.iloc[:len(q)]  # Ensure we only take as many rows as we have coordinates
    df.loc[:,"q"] = q
    df.loc[:,"r"] = r

    # sum rows with identical coordinates together
    # Group by hexagonal coordinates and aggregate
    agg_dict = {'Count': 'sum', 'ChEBI': list}
    if 'TFIDF' in df.columns:
        agg_dict['TFIDF'] = 'sum'
    else:
        df['TFIDF'] = df['Count']  # Use Count as TFIDF if not available
        agg_dict['TFIDF'] = 'sum'
    
    df_grouped = df.groupby(['q', 'r']).agg(agg_dict).reset_index()
    
    # Add gaussian blur - simplified version
    try:
        df_blurred = add_gaussian_blur(df_grouped, BLUR_MAX, BLUR_STEP_SIZE)
        df_final = add_tooltip_columns(df_blurred, table_filtered)
    except Exception as e:
        print(f"Warning: Advanced features failed ({e}), using basic version")
        # Fallback to basic version without blur and advanced tooltips
        df_final = df_grouped.copy()
        
        # Add basic blur columns
        for sd_x in np.arange(0, BLUR_MAX + 0.25, 0.25):
            df_final[str(sd_x)] = df_final['Count']
            if f'{sd_x}_tfidf' not in df_final.columns:
                df_final[f'{sd_x}_tfidf'] = df_final['TFIDF']
    
    # Clean up and ensure required columns exist
    if 'ChEBI' in df_final.columns:
        df_final = df_final.drop(columns='ChEBI', errors='ignore')
    
    df_final.loc[:,"Count_total"] = df_final.loc[:,"Count"]
    
    # Ensure all required numeric columns are present and valid
    numeric_columns = ['Count', 'TFIDF', 'Count_total']
    for col in numeric_columns:
        if col in df_final.columns:
            df_final[col] = pd.to_numeric(df_final[col], errors='coerce').fillna(0)

    # plot title and source
    title = f'Hexbin plot for {len(x)} annotated chemicals with query {term}'
    source = ColumnDataSource(df_final)
    return source, title

def create_stats_description(table):
    '''
    Enhanced function with drug-likeness and extended property statistics
    '''
    total_count = table.Count.sum()
    
    # Basic molecular properties
    stats_listed = [
        'Enhanced Chemical Property Statistics',
        f'Total amount of chemicals: {total_count}',
    ]
    
    # Core properties
    properties = {
        'logP': 'LogP',
        'Mass': 'Molecular Weight (Da)',
        'HBD': 'Hydrogen Bond Donors',
        'HBA': 'Hydrogen Bond Acceptors',
        'PSA': 'Polar Surface Area (Ų)',
        'RotBonds': 'Rotatable Bonds'
    }
    
    for prop_col, prop_name in properties.items():
        if prop_col in table.columns:
            prop_data = pd.to_numeric(table[prop_col], errors='coerce')
            prop_data = prop_data.dropna()
            
            if len(prop_data) > 0:
                # Weight by count for distribution stats
                weighted_data = np.repeat(prop_data, table.loc[prop_data.index, 'Count'])
                
                stats_listed.extend([
                    f'{prop_name} mean: {weighted_data.mean():.3f}',
                    f'{prop_name} std dev: {weighted_data.std():.3f}',
                    f'{prop_name} median: {np.median(weighted_data):.3f}'
                ])
    
    # Drug-likeness statistics
    drug_stats = analyze_drug_likeness(table)
    if drug_stats:
        stats_listed.append('<br><b>Drug-likeness Analysis:</b>')
        
        if 'lipinski_compliant' in drug_stats:
            stats_listed.append(f"Lipinski RO5 compliant: {drug_stats['lipinski_compliant']} ({drug_stats['lipinski_percentage']}%)")
        
        if 'veber_compliant' in drug_stats:
            stats_listed.append(f"Veber rule compliant: {drug_stats['veber_compliant']} ({drug_stats['veber_percentage']}%)")
        
        if 'drug_like' in drug_stats:
            stats_listed.append(f"Overall drug-like: {drug_stats['drug_like']} ({drug_stats['drug_like_percentage']}%)")
    
    # Publication year analysis
    pub_years_stats = analyze_publication_years(table)
    
    if pub_years_stats['has_data']:
        stats_listed.append('<br><b>Publication Timeline:</b>')
        stats_listed.append(f"Year range: {pub_years_stats['min_year']} - {pub_years_stats['max_year']}")
        stats_listed.append(f"Mean publication year: {pub_years_stats['mean_year']}")
        stats_listed.append(f"Median publication year: {pub_years_stats['median_year']}")
        
        # Add year distribution table
        stats_listed.append('<br><b>Publications per year:</b>')
        year_dist_html = '<table border="1" cellpadding="3" style="border-collapse: collapse;">'
        year_dist_html += '<tr><th>Year</th><th>Publications</th></tr>'
        
        sorted_years = sorted(pub_years_stats['year_counts'].keys())
        for year in sorted_years:
            count = pub_years_stats['year_counts'][year]
            year_dist_html += f'<tr><td>{year}</td><td>{count}</td></tr>'
        
        year_dist_html += '</table>'
        stats_listed.append(year_dist_html)
    else:
        stats_listed.append('<br><b>Publication Timeline:</b> Data not available')
    
    stats_description = return_html(stats_listed)
    return stats_description

def return_html(metadata):
    '''
    This function receives the metadata. The metadata is put in the html string and this string returned.
    '''
    html_content = f"""
    <HTML>
    <HEAD>
    <TITLE>{metadata[0]}</TITLE>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.4; }}
        h1 {{ color: #444; }}
        .stats-container {{ margin-top: 20px; }}
        .drug-stats {{ margin-top: 15px; background-color: #f0f8ff; padding: 10px; border-radius: 5px; }}
        .year-stats {{ margin-top: 15px; }}
        table {{ margin-top: 10px; }}
        th {{ background-color: #e6e6e6; }}
    </style>
    </HEAD>
    <BODY BGCOLOR="FFFFFF">
    <h1>{metadata[0]}</h1>
    <div class="stats-container">
    """
    
    # Add all the metadata items
    for i in range(1, len(metadata)):
        html_content += f"<div>{metadata[i]}</div>\n"
    
    html_content += """
    </div>
    </BODY>
    </HTML>
    """
    return html_content

def get_tables(files):
    '''
    Enhanced function to import tables with extended property support
    '''
    tables = dict()
    for file in files:
        table = import_table(file)
        term = file.split('/')[1].split('_table')[0]
        metadata_file = f'metadata/{term}.txt'
        
        try:
            with open(metadata_file, 'r') as f:
                metadata = f.readlines()
        except:
            print(f"Warning: Metadata file for {term} not found or couldn't be read.")
            metadata = [f"Metadata for {term}", "No metadata available"]
        
        # Add information about available extended properties
        extended_props = ['HBD', 'HBA', 'PSA', 'RotBonds', 'Role', 'Lipinski_RO5', 'Veber_Rule', 'Drug_Like']
        available_props = [prop for prop in extended_props if prop in table.columns]
        
        if available_props:
            metadata.append(f"Extended properties available: {', '.join(available_props)}\n")
        
        # Check for publication years data
        if 'PublicationYears' in table.columns:
            pub_years_stats = analyze_publication_years(table)
            if pub_years_stats['has_data']:
                year_info = f"Publication years range: {pub_years_stats['min_year']} - {pub_years_stats['max_year']}\n"
                if not any("publication years range" in line.lower() for line in metadata):
                    metadata.append(year_info)
        
        tables[term] = {'table': table, 'metadata': metadata}

    return tables

def get_files(folder):
    '''
    This function receives the input folder and returns a list of file paths in this folder.
    '''
    files = [f'{folder}/{file}' for file in os.listdir(folder) if '.pkl' in file]
    return files

def return_JS_code(widget):
    '''
    Enhanced JavaScript code with drug-likeness highlighting support
    '''
    if widget == 'multi_select_class':
        code = """
            var source_data = source.data;
            var source_class_data = source_class.data;
            var term_to_source = term_to_source;
            var term = cb_obj.value[0];
            var f = slider1.value
            var sd_x = slider2.value;
            var p = p;
            var mapper = mapper;
            var class_hex = class_hex;
            var checkbox_class = checkbox_class
            var term_to_class = term_to_class;
            var term_to_stats = term_to_stats;
            var stats_button = stats_button;

            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf
            if (checkbox.active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // select new_data
            var new_data = term_to_source[term]['source'].data
            var new_class = term_to_class[term]['source'].data

            // replace current data
            for (var key in new_data) {
                source_data[key] = [];
                for (i=0;i<new_data[key].length;i++) {
                    source_data[key].push(new_data[key][i]);
                }
            }

            // class
            class_hex.visible = false

            for (var key in new_class){
                source_class_data[key] = [];
                for (i=0;i<new_class[key].length;i++) {
                    source_class_data[key].push(new_class[key][i])
                }
            }

            class_hex.source = source_class_data
            checkbox_class.active = []

            // new title
            var title = term_to_source[term]['title']
            p.title.text = title

            // apply blur and saturation
            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }

            // maximum value for fill_color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // Apply changes
            source_class.change.emit();
            source.change.emit();
        """
        
    elif widget == 'drug_likeness':
        code = """
            var source_drug_data = source_drug.data;
            var term_to_drug_sources = term_to_drug_sources;
            var term = multi_select.value[0];
            var active = cb_obj.active;
            var drug_hex = drug_hex;

            if (active.length > 0) {
                // Show drug-like highlighting
                var drug_data = term_to_drug_sources[term]['source'].data;
                
                for (var key in drug_data) {
                    source_drug_data[key] = [];
                    for (i=0; i<drug_data[key].length; i++) {
                        source_drug_data[key].push(drug_data[key][i]);
                    }
                }
                
                drug_hex.visible = true;
            } else {
                // Hide drug-like highlighting
                drug_hex.visible = false;
            }
            
            source_drug.change.emit();
        """
    
    elif widget == 'multi_select':
        code = """
            var source_data = source.data;
            var term_to_source = term_to_source;
            var term = cb_obj.value[0];
            var f = slider1.value
            var sd_x = slider2.value;
            var p = p;
            var mapper = mapper;

            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf
            if (checkbox.active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // select new_data
            var new_data = term_to_source[term]['source'].data

            // replace current data
            for (var key in new_data) {
                source_data[key] = [];
                for (i=0;i<new_data[key].length;i++) {
                    source_data[key].push(new_data[key][i]);
                }
            }

            // new title
            var title = term_to_source[term]['title']
            p.title.text = title

            // apply blur and saturation
            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }

            // maximum value for fill_color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // apply changes
            source.change.emit();
        """
        
    elif widget == 'tooltips':
        code = """
            <style>
            table, td, th {
              border-collapse: collapse;
              border: 1px solid #dddddd;
              padding: 2px;
              table-layout: fixed;
              height: 20px;
            }

            tr:nth-child(even) {
              background-color: #dddddd;
            }
            </style>

            <table>
              <col width="100">
              <col width="80">
              <col width="80">
              <col width="305">
              <tr>
                <th>Total Counts</th>
                <th>@Count_total</th>
                <th>@TFIDF</th>
                <th>($x, $y)</th>
              </tr>
              <tr style="color: #fff; background: black;">
                <th>ChEBI ID</th>
                <th>Count</th>
                <th>TFIDF</th>
                <th>Name</th>
              </tr>
              <tr>
                <th>@ChEBI1</th>
                <th>@Count1</th>
                <th>@TFIDF1</th>
                <th>@Names1</th>
              </tr>
              <tr>
                <th>@ChEBI2</th>
                <th>@Count2</th>
                <th>@TFIDF2</th>
                <th>@Names2</th>
              </tr>
              <tr>
                <th>@ChEBI3</th>
                <th>@Count3</th>
                <th>@TFIDF3</th>
                <th>@Names3</th>
              </tr>
            </table>
            """

    elif widget == 'slider2':
        code = """
            var source_data = source.data;
            var f = slider1.value
            var mapper = mapper;
            var checkbox = checkbox;
            var sd_x = cb_obj.value;

            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf
            if (checkbox.active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // apply blur and saturation
            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }

            // maximum value for fill_color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // apply changes
            source.change.emit();
        """

    elif widget == 'slider1':
        code = """
            var mapper = mapper
            var source_data = source.data;
            var f = cb_obj.value;
            var checkbox = checkbox;
            var sd_x = slider2.value;
            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf
            if (checkbox.active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // apply scaling
            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }

            // maximum value for fill color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // apply changes
            source.change.emit();
            """

    elif widget == 'rbg':
        code = """
            var active = cb_obj.active;
            var p = p;
            var Viridis256 = Viridis256;
            var Greys256 = Greys256;

            if (active == 1){
            mapper.transform.palette = Greys256
            p.background_fill_color = '#000000'
            }

            if (active == 0){
            mapper.transform.palette = Viridis256
            p.background_fill_color = '#440154'
            }
            """

    elif widget == 'checkbox':
        code = """
            var source_data = source.data;
            var active = cb_obj.active
            var f = slider1.value
            var mapper = mapper
            var sd_x = slider2.value;

            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf
            if (active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // apply scaling
            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }

            // maximum value for fill color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // apply changes
            source.change.emit();
            """

    elif widget == 'hover':
        code = """
            var tooltips = document.getElementsByClassName("bk-tooltip");
            for (var i = 0, len = tooltips.length; i < len; i ++) {
                tooltips[i].style.top = ""; // unset what bokeh.js sets
                tooltips[i].style.left = "";
                tooltips[i].style.bottom = "150px";
                tooltips[i].style.left = "575px";
                tooltips[i].style.width = "500px";
            }
            """

    elif widget == 'button':
        code = """
            var term_to_metadata = term_to_metadata
            var term = multi_select.value[0]
            var metadata = term_to_metadata[term]
            var wnd = window.open("about:blank", "", "_blank");
            wnd.document.write(metadata)
            """
            
    elif widget == 'stats':
        code = """
            var term_to_stats = term_to_stats
            var term = multi_select.value[0]
            var stats_description = term_to_stats[term]
            var wnd = window.open("about:blank", "", "_blank");
            wnd.document.write(stats_description)
            """

    elif widget == 'class':
        code = """
            var term_to_class = term_to_class;
            var multi_select = multi_select;
            var class_hex = class_hex;
            var active = cb_obj.active;

            if (active.length == 1) {
                class_hex.visible = true
            } else {
                class_hex.visible = false
            }
            """
            
    return code

def plot(tables, output_filename, xmin, xmax, ymin, ymax, class_id):
    '''
    Enhanced plot function with drug-likeness highlighting
    '''
    if not os.path.isdir('plots'):
        os.mkdir('plots')
    file_name = f'plots/{output_filename}.html'
    output_file(file_name)
    
    # Blur and saturation values
    BLUR_MAX = 4
    BLUR_STEP_SIZE = 0.25
    SATURATION_MAX = 5
    SATURATION_STEP_SIZE = 0.25

    # Hexagon plot properties
    SIZE_HEXAGONS = 10
    orientation = 'flattop'
    ratio = ((ymax-ymin) / (xmax-xmin))
    size = SIZE_HEXAGONS / ratio
    hexagon_height = sqrt(3) * size
    hexagon_height = hexagon_height*ratio

    # Make figure
    p = figure(x_range=[xmin, xmax], y_range=[ymin-(hexagon_height/2), ymax],
               tools="wheel_zoom,reset,save", background_fill_color='#440154')

    p.grid.visible = False
    p.xaxis.axis_label = "log(P)"
    p.yaxis.axis_label = "mass in Da"
    p.xaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_font_style = 'normal'

    # Source for widgets
    term_to_source = dict()
    term_to_class = dict()
    term_to_drug_sources = dict()
    term_to_metadata = dict()
    term_to_stats = dict()

    options = []
    
    # Check if any tables have drug-likeness data
    has_drug_data = any('Drug_Like' in tables[term]['table'].columns for term in tables.keys())
    
    # Loop for plot sources
    for term in tables.keys():
        print(f'Sourcing {term}')

        # Add term to the options for the multiplot
        options.append((term, term))

        # Get table
        table = tables[term]['table']

        # Plot sources
        source, title = create_data_source(table, term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE)
        source_class = create_class_source(table, size, ratio, orientation, class_id)
        
        # Drug-likeness sources
        drug_sources = {}
        if 'Drug_Like' in table.columns:
            drug_sources['Drug_Like'] = create_drug_like_source(table, size, ratio, orientation, 'Drug_Like')
        if 'Lipinski_RO5' in table.columns:
            drug_sources['Lipinski_RO5'] = create_drug_like_source(table, size, ratio, orientation, 'Lipinski_RO5')
        if 'Veber_Rule' in table.columns:
            drug_sources['Veber_Rule'] = create_drug_like_source(table, size, ratio, orientation, 'Veber_Rule')
        
        stats_description = create_stats_description(table)
        
        term_to_source[term] = {'source': source, 'title': title}
        term_to_class[term] = {'source': source_class, 'show_class': True}
        term_to_drug_sources[term] = {'source': drug_sources.get('Drug_Like', ColumnDataSource())}
        term_to_stats[term] = stats_description

        # Metadata
        metadata = return_html(tables[term]['metadata'])
        term_to_metadata[term] = metadata

    # Make default source for plot
    default_term = list(tables.keys())[0]
    table = tables[default_term]['table']
    source, title = create_data_source(table, default_term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE)
    p.title.text = title
    metadata = tables[default_term]['metadata']
    metadata = return_html(metadata)

    # Color mapper
    mapper = linear_cmap('Count', 'Viridis256', 0, max(source.data['Count']))
    hex = p.hex_tile(q="q", r="r", size=size, line_color=None, source=source, aspect_scale=ratio, orientation=orientation,
           fill_color=mapper)

    # Class highlighting (existing functionality)
    if class_id:
        source_class = create_class_source(table, size, ratio, orientation, class_id)
        class_hex = p.hex_tile(q='q', r="r", size=size, line_color=None, source=source_class, aspect_scale=ratio,orientation=orientation,
            fill_color='#ff007f')
        class_hex.visible = False

    # Drug-likeness highlighting
    if has_drug_data and 'Drug_Like' in table.columns:
        source_drug = create_drug_like_source(table, size, ratio, orientation, 'Drug_Like')
        drug_hex = p.hex_tile(q='q', r="r", size=size, line_color='#00ff00', line_width=2, 
                             source=source_drug, aspect_scale=ratio, orientation=orientation,
                             fill_color='#00ff00', fill_alpha=0.3)
        drug_hex.visible = False

    # HOVER
    TOOLTIPS = return_JS_code('tooltips')
    code_callback_hover = return_JS_code('hover')
    callback_hover = CustomJS(code=code_callback_hover)
    hover = HoverTool(renderers=[hex], tooltips=TOOLTIPS, callback=callback_hover, show_arrow=False)
    p.add_tools(hover)

    # Widgets
    slider1 = Slider(start=1, end=SATURATION_MAX, value=1, step=SATURATION_STEP_SIZE, title="Saturation", width=100)
    slider2 = Slider(start=0, end=BLUR_MAX, value=0, step=BLUR_STEP_SIZE, title="Blur", width=100)
    multi_select = MultiSelect(title=output_filename, value=[default_term], options=options, width=100, height=200)
    checkbox = CheckboxGroup(labels=["TFIDF"], active=[])
    radio_button_group = RadioGroup(labels=["Viridis256", "Greys256"], active=0)
    button = Button(label="Metadata", button_type="default", width=100)
    stats = Button(label="Statistics", button_type="default", width=100)
    
    # Drug-likeness widgets
    drug_widgets = []
    if has_drug_data:
        checkbox_druglike = CheckboxGroup(labels=["Highlight Drug-like"], active=[])
        drug_widgets.append(checkbox_druglike)
    
    if class_id:
        checkbox_class = CheckboxGroup(labels=[f"Show {class_id}"], active=[])

    # JavaScript code
    code_callback_slider1 = return_JS_code('slider1')
    code_callback_slider2 = return_JS_code('slider2')
    code_callback_ms = return_JS_code('multi_select')
    code_callback_checkbox = return_JS_code('checkbox')
    code_callback_rbg = return_JS_code('rbg')
    code_callback_button = return_JS_code('button')
    code_callback_stats = return_JS_code('stats')
    
    if class_id:
        code_callback_ms = return_JS_code('multi_select_class')
        code_callback_class = return_JS_code('class')

    # Drug-likeness callbacks
    if has_drug_data:
        code_callback_druglike = return_JS_code('drug_likeness')

    # Callbacks
    callback_slider1 = CustomJS(args={'source': source, 'mapper': mapper, 'slider2': slider2, 'checkbox': checkbox}, code=code_callback_slider1)
    callback_slider2 = CustomJS(args={'source': source, 'mapper': mapper, 'slider1': slider1, 'checkbox': checkbox}, code=code_callback_slider2)
    callback_ms = CustomJS(args={'source': source, 'term_to_source': term_to_source, 'slider1': slider1, 'slider2': slider2, 'checkbox': checkbox, 'p': p, 'mapper': mapper}, code=code_callback_ms)
    callback_checkbox = CustomJS(args={'source': source, 'slider2': slider2, 'mapper': mapper, 'slider1': slider1}, code=code_callback_checkbox)
    callback_radio_button_group = CustomJS(args={'p': p, 'multi_select': multi_select, 'mapper': mapper, 'term_to_class': term_to_class, 'Viridis256': Viridis256, 'Greys256': Greys256}, code=code_callback_rbg)
    callback_button = CustomJS(args={'term_to_metadata': term_to_metadata, 'multi_select': multi_select}, code=code_callback_button)
    callback_stats = CustomJS(args={'term_to_stats': term_to_stats, 'multi_select': multi_select}, code=code_callback_stats)
    
    if class_id:
        callback_ms = CustomJS(args={'source': source, 'term_to_source': term_to_source, 'slider1': slider1, 'slider2': slider2,
        'checkbox': checkbox, 'p': p, 'mapper': mapper, 'source_class': source_class, 'term_to_class': term_to_class, 'class_hex': class_hex,
        'checkbox_class': checkbox_class}, code=code_callback_ms)
        callback_class = CustomJS(args={'multi_select': multi_select, 'term_to_class': term_to_class, 'class_hex': class_hex}, code=code_callback_class)

    if has_drug_data:
        callback_druglike = CustomJS(args={'source_drug': source_drug, 'term_to_drug_sources': term_to_drug_sources, 'multi_select': multi_select, 'drug_hex': drug_hex}, code=code_callback_druglike)

    # On change
    slider1.js_on_change('value', callback_slider1)
    slider2.js_on_change('value', callback_slider2)
    multi_select.js_on_change("value", callback_ms)
    checkbox.js_on_change('active', callback_checkbox)
    radio_button_group.js_on_change('active', callback_radio_button_group)
    button.js_on_event(events.ButtonClick, callback_button)
    stats.js_on_event(events.ButtonClick, callback_stats)
    
    if class_id:
        checkbox_class.js_on_change('active', callback_class)
    
    if has_drug_data:
        checkbox_druglike.js_on_change('active', callback_druglike)

    # Layout
    widgets_column = [slider1, slider2, checkbox, radio_button_group]
    
    if has_drug_data:
        widgets_column.extend(drug_widgets)
    
    if class_id:
        widgets_column.append(checkbox_class)
    
    widgets_column.extend([button, stats])
    
    layout = row(multi_select, p, column(*widgets_column))
    show(layout)

def parser():
    parser = argparse.ArgumentParser(description='Enhanced visualization with drug-likeness analysis')
    parser.add_argument('-i', required=True, metavar='input_folder', dest='input_folder', help='Input folder with table files')
    parser.add_argument('-o', required=True, metavar='output_filename', dest='output_filename', help='Output HTML filename')
    parser.add_argument('-xmin', required=False, metavar='xmin', dest='xmin', help='X axis minimum (logP), default is -5')
    parser.add_argument('-xmax', required=False, metavar='xmax', dest='xmax', help='X axis maximum (logP), default is 10')
    parser.add_argument('-ymax', required=False, metavar='ymax', dest='ymax', help='Y axis maximum (mass in Da), default is 1600')
    parser.add_argument('-class', required=False, metavar='class_id', dest='class_id', help='Class to highlight in the plot')
    arguments = parser.parse_args()
    return arguments

def main():
    start_time = datetime.now()

    args = parser()
    folder = args.input_folder
    output_filename = args.output_filename
    xmin = args.xmin
    xmax = args.xmax
    ymax = args.ymax

    # Default settings
    if not xmin:
        xmin = -5
    else:
        xmin = float(args.xmin)
    if not xmax:
        xmax = 10
    else:
        xmax = float(args.xmax)
    if not ymax:
        ymax = 1600
    else:
        ymax = float(args.ymax)
    ymin = 0  # cannot be negative

    # Create required directories if they don't exist
    required_dirs = ['plots', 'tables', 'metadata', 'publication_dates']
    for directory in required_dirs:
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")

    files = get_files(folder)
    if not files:
        print(f"Error: No .pkl files found in {folder}")
        print("Make sure to run make_table.py first to generate the required files.")
        return
    
    print(f"Found {len(files)} table files to process")
    
    tables = get_tables(files)
    
    # Report on available features
    feature_summary = {}
    for term, data in tables.items():
        table = data['table']
        features = []
        
        if 'Drug_Like' in table.columns:
            features.append('Drug-likeness')
        if 'Lipinski_RO5' in table.columns:
            features.append('Lipinski RO5')
        if 'Veber_Rule' in table.columns:
            features.append('Veber Rule')
        if 'PublicationYears' in table.columns:
            features.append('Publication Timeline')
        
        extended_props = [col for col in ['HBD', 'HBA', 'PSA', 'RotBonds', 'Role'] if col in table.columns]
        if extended_props:
            features.append(f"Extended Properties ({', '.join(extended_props)})")
        
        feature_summary[term] = features
    
    print("\nEnhanced Features Available:")
    for term, features in feature_summary.items():
        if features:
            print(f"  {term}: {', '.join(features)}")
        else:
            print(f"  {term}: Basic properties only")

    plot(tables, output_filename, xmin, xmax, ymin, ymax, args.class_id)

    print(f"\n✓ Enhanced visualization completed in {datetime.now() - start_time}")
    print(f"Interactive plot saved as: plots/{output_filename}.html")
    print("\nNew features include:")
    print("- Drug-likeness rule highlighting")
    print("- Extended chemical property analysis")
    print("- Enhanced publication timeline statistics")
    print("- Improved filtering transparency")

if __name__ == '__main__':
    main()