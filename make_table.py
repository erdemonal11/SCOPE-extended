#!/usr/bin/python

import argparse
import math
import sys
import pandas as pd
import numpy as np
import os
from pathlib import PurePath

def import_properties():
    '''
    This function reads a file and adds the information + id in a dictionary. This is added to another dictionary as value, and the term that describes the
    information (e.g. "Names") as key.
    '''
    FOLDER = 'files'
    files = os.listdir(FOLDER)
    data = dict()
    for file in files:
        path = os.path.join(FOLDER, file)
        key = file.split('2')[1].split('_')[0]
        if '.pkl' in file:
            df = pd.read_pickle(path)
            df['ChEBI'] = df['ChEBI'].astype(int)
            df = df.set_index('ChEBI')
        else:
            df = pd.read_csv(path, sep='\t', header=None, names=['ChEBI', 'Info'], index_col='ChEBI', dtype={"ChEBI": "int", "Info": "str"})
        id_to_info = df.to_dict()
        data[key] = id_to_info
    return data

def get_info(data, key, id):
    '''
    This function recieves:
        - a ChEBI identifier
        - a dictionary key indicating the type of information
        - data dictionary containing information of every ChEBI identifier
    The function returns specific information retrieved from the data dictionary (e.g. Name or Mass)
    '''
    try:
        info = data[key]['Info'][id]
        # Ensure we're returning a scalar value
        if isinstance(info, (list, dict, pd.Series)):
            # Convert to string if it's a complex structure
            info = str(info)
    except:
        info = float('NaN')
    return info

def import_publication_dates(term):
    '''
    This function imports publication dates for a specific search term.
    Returns a dictionary mapping publication IDs to years.
    '''
    pub_dates_dict = {}
    pub_dates_file = os.path.join('publication_dates', f'{term}_pub_dates.tsv')
    
    if os.path.exists(pub_dates_file):
        try:
            pub_dates_df = pd.read_csv(pub_dates_file, sep='\t')
            pub_dates_dict = dict(zip(pub_dates_df['PublicationID'], pub_dates_df['PublicationYear']))
        except Exception as e:
            print(f"Error reading publication dates file: {str(e)}")
    else:
        print(f"Publication dates file not found for term: {term}")
    
    return pub_dates_dict

def make_table(data, df_results, publication_dates=None):
    '''
    This function recieves the data of all the ChEBI files in the files folder and the ids of the query search.
    It returns a dictionary of the query ids and their properties from the ChEBI files, if those properties are there.
    e.g. if there is no logP value, the id is not added to the dictionary that is returned.
    '''
    table = df_results.copy()  # Use copy to avoid modifying the original
    
    # create column lists
    for key in data.keys():
        # Process one column at a time to handle potential errors
        try:
            column_list = []
            for id in df_results.index.to_list():
                column_list.append(get_info(data, key, id))
                
            # Convert to series first to handle mixed types better
            table[key] = pd.Series(column_list, index=df_results.index)
        except Exception as e:
            print(f"Error processing column '{key}': {str(e)}")
            # Skip columns that cause errors
            continue

    # Handle special cases for Mass and logP
    if 'Mass' in table.columns:
        table['Mass'] = pd.to_numeric(table['Mass'].replace('-', np.nan), errors='coerce')
    if 'logP' in table.columns:
        table['logP'] = pd.to_numeric(table['logP'].replace('-', np.nan), errors='coerce')
    
    # Add publication years if available
    if publication_dates:
        # Process the publications for each ChEBI ID
        pub_years = {}
        for chebi_id, row in df_results.iterrows():
            # Find all publication IDs for this ChEBI ID
            pub_ids = []
            if isinstance(df_results, pd.DataFrame) and 'Publication' in df_results.columns:
                pub_ids = [row['Publication']]
            elif 'Publication' in df_results.reset_index().columns:
                pub_id_series = df_results.reset_index().loc[df_results.reset_index()['ChEBI'] == chebi_id, 'Publication']
                if not pub_id_series.empty:
                    pub_ids = pub_id_series.tolist()
            
            # Get years for these publication IDs
            years = [int(publication_dates.get(pid)) for pid in pub_ids if pid in publication_dates and publication_dates.get(pid)]
            if years:
                pub_years[chebi_id] = years
        
        # Add to the table
        if pub_years:
            table['PublicationYears'] = pd.Series({idx: pub_years.get(idx, []) for idx in table.index})
    
    table = table.dropna()
    table = table.sort_values(by='Count', ascending=False)
    return table

def write_to_file(table, term):
    '''
    This function writes the table in a .pkl file for easy importation into the visualization script,
    and a .tsv for own inspection. The file is named with the (shortest) query search term.
    '''
    # Make sure the tables directory exists
    os.makedirs('tables', exist_ok=True)
    
    path = 'tables/'+term
    table.to_csv(path+'_table.tsv', sep='\t')
    table.to_pickle(path+'_table.pkl')

def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=True, metavar='input', dest='input', help='[i] to select input folder or input file from the results folder ')
    parser.add_argument('-t', required=True, metavar='type', dest='type', help='[t] to select type of input: file or folder')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parser()
    input_type = args.type
    input = args.input
    if input_type == 'file':
        results = [input]
    elif input_type == 'folder':
        files = os.listdir(input)
        results = [input+'/'+file for file in files]
    else:
        sys.exit('Error: please give \'file\' or \'folder\' as input type')

    #gather properties
    data = import_properties()

    for result in results:
        try:
            term = PurePath(result).parts[1].split('_ChEBI_IDs.tsv')[0]
            print('making table for %s' % term)

            # Import publication dates if available
            publication_dates = import_publication_dates(term)

            # import results
            df = pd.read_csv(result, sep = '\t', names=['ChEBI', 'Publication'], dtype={"ChEBI": "int", "Publication": "str"})
            df_results = df.groupby(by=['ChEBI']).count().rename(columns={"Publication": "Count"})

            # make table
            table = make_table(data, df_results, publication_dates)

            # Check if 'idf' column exists before performing normalization
            if 'idf' in table.columns:
                # ensure numeric types for calculation
                table["Count"] = pd.to_numeric(table["Count"], errors='coerce')
                table["idf"] = pd.to_numeric(table["idf"], errors='coerce')
                
                # perform normalization
                table.loc[:,"TFIDF"] = table["Count"].astype(float) * table["idf"].astype(float)
                table.loc[:,"TFIDF"] = table.loc[:,"TFIDF"].round(decimals=0).astype(int)
            else:
                print(f"Warning: 'idf' column not found for {term}, skipping TFIDF calculation")
            
            print(table)

            # write table to file
            write_to_file(table, term)
        except Exception as e:
            print(f"Error processing {result}: {str(e)}")
            continue

if __name__ == '__main__':
    main()