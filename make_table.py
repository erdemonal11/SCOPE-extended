#!/usr/bin/python

import argparse
import math
import sys
import pandas as pd
import numpy as np
import os
from pathlib import PurePath

# Default stop list of common non-ligands (ChEBI IDs)
DEFAULT_STOP_LIST = [
    '15377',  # water
    '26710',  # sodium chloride
    '9754',   # Tris
    '42334',  # HEPES
    '28262',  # DMSO
    '16236',  # ethanol
    '17883',  # methanol
    '15347',  # acetonitrile
    '27385',  # potassium chloride
    '9344',   # phosphate
    '16810',  # hydrogen chloride
    '26836',  # magnesium chloride
]

def load_stop_list(stop_list_path=None):
    """Load stop list from file or use default"""
    if stop_list_path and os.path.exists(stop_list_path):
        with open(stop_list_path, 'r') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return DEFAULT_STOP_LIST

def apply_stop_list_filter(df_results, stop_list):
    """Remove ChEBI IDs in stop list from results"""
    stop_ids = []
    for id_str in stop_list:
        try:
            stop_ids.append(int(id_str))
        except ValueError:
            print(f"Warning: Invalid ChEBI ID in stop list: {id_str}")
    
    initial_count = len(df_results)
    df_results = df_results[~df_results.index.isin(stop_ids)]
    removed_count = initial_count - len(df_results)
    
    if removed_count > 0:
        print(f"Stop list filter removed {removed_count} common non-ligands")
    
    return df_results

def filter_by_roles(df_results, data, exclude_roles=None, include_roles=None):
    """Filter results based on ChEBI roles/classes"""
    if 'Role' not in data:
        if exclude_roles or include_roles:
            print("Warning: Role information not available for filtering")
        return df_results
    
    role_dict = data['Role']['Info']
    initial_count = len(df_results)
    
    if exclude_roles:
        exclude_list = [role.strip().lower() for role in exclude_roles.split(',')]
        exclude_ids = []
        
        for chebi_id, role_info in role_dict.items():
            if isinstance(role_info, str):
                role_lower = role_info.lower()
                if any(ex_role in role_lower for ex_role in exclude_list):
                    exclude_ids.append(chebi_id)
        
        df_results = df_results[~df_results.index.isin(exclude_ids)]
        removed = initial_count - len(df_results)
        print(f"Role exclusion filter removed {removed} chemicals with roles: {exclude_roles}")
    
    if include_roles:
        include_list = [role.strip().lower() for role in include_roles.split(',')]
        include_ids = []
        
        for chebi_id, role_info in role_dict.items():
            if isinstance(role_info, str):
                role_lower = role_info.lower()
                if any(in_role in role_lower for in_role in include_list):
                    include_ids.append(chebi_id)
        
        df_results = df_results[df_results.index.isin(include_ids)]
        kept = len(df_results)
        print(f"Role inclusion filter kept {kept} chemicals with roles: {include_roles}")
    
    return df_results

def calculate_drug_likeness_rules(table):
    """Calculate compliance with drug-likeness rules"""
    
    def check_lipinski_ro5(row):
        """Lipinski's Rule of Five"""
        try:
            mw = pd.to_numeric(row.get('Mass', np.nan), errors='coerce')
            logp = pd.to_numeric(row.get('logP', np.nan), errors='coerce')
            hbd = pd.to_numeric(row.get('HBD', np.nan), errors='coerce')
            hba = pd.to_numeric(row.get('HBA', np.nan), errors='coerce')
            
            if any(pd.isna([mw, logp])):  # Only require Mass and logP for basic rule
                return False
            
            conditions = [
                mw <= 500,     # Molecular weight ≤ 500 Da
                logp <= 5,     # LogP ≤ 5
            ]
            
            # Add HBD/HBA conditions if available
            if not pd.isna(hbd):
                conditions.append(hbd <= 5)      # HBD ≤ 5
            if not pd.isna(hba):
                conditions.append(hba <= 10)     # HBA ≤ 10
            
            return sum(conditions) >= 2  # At least 2 conditions met (relaxed for missing data)
        except:
            return False
    
    def check_veber_rule(row):
        """Veber's rule for oral bioavailability"""
        try:
            psa = pd.to_numeric(row.get('PSA', np.nan), errors='coerce')
            rotbonds = pd.to_numeric(row.get('RotBonds', np.nan), errors='coerce')
            
            if pd.isna(psa) or pd.isna(rotbonds):
                return False
            
            return psa <= 140 and rotbonds <= 10
        except:
            return False
    
    # Apply rules if we have the required properties
    if all(prop in table.columns for prop in ['Mass', 'logP']):
        table['Lipinski_RO5'] = table.apply(check_lipinski_ro5, axis=1)
        compliant_count = table['Lipinski_RO5'].sum()
        print(f"Calculated Lipinski RO5 compliance for {compliant_count} compounds")
    
    if all(prop in table.columns for prop in ['PSA', 'RotBonds']):
        table['Veber_Rule'] = table.apply(check_veber_rule, axis=1)
        compliant_count = table['Veber_Rule'].sum()
        print(f"Calculated Veber rule compliance for {compliant_count} compounds")
    
    if 'Lipinski_RO5' in table.columns and 'Veber_Rule' in table.columns:
        table['Drug_Like'] = table['Lipinski_RO5'] & table['Veber_Rule']
        drug_like_count = table['Drug_Like'].sum()
        print(f"Identified {drug_like_count} drug-like compounds")
    elif 'Lipinski_RO5' in table.columns:
        # Use only Lipinski if Veber data is not available
        table['Drug_Like'] = table['Lipinski_RO5']
        drug_like_count = table['Drug_Like'].sum()
        print(f"Identified {drug_like_count} drug-like compounds (Lipinski RO5 only)")
    
    return table

def import_properties():
    '''
    Enhanced function to read ChEBI property files including roles and additional properties
    '''
    FOLDER = 'files'
    if not os.path.exists(FOLDER):
        print(f"Error: {FOLDER} directory not found. Run download_files.py first.")
        return {}
    
    files = os.listdir(FOLDER)
    data = dict()
    
    # Property mapping for new files
    property_mapping = {
        'role': 'Role',
        'application': 'Application',
        'hbd': 'HBD',
        'hba': 'HBA',
        'psa': 'PSA',
        'rotbonds': 'RotBonds',
        'rotatable': 'RotBonds',
        'smiles': 'SMILES'
    }
    
    for file in files:
        path = os.path.join(FOLDER, file)
        print(f"Processing file: {file}")
        
        # Determine property type
        key = None
        file_lower = file.lower()
        
        for prop_identifier, prop_key in property_mapping.items():
            if prop_identifier in file_lower:
                key = prop_key
                break
        
        if key is None:
            # Use existing logic for standard files
            if '_' in file and '2' in file:
                try:
                    key = file.split('2')[1].split('_')[0]
                except IndexError:
                    key = file.split('.')[0]
            else:
                key = file.split('.')[0]
        
        try:
            if '.pkl' in file:
                df = pd.read_pickle(path)
                if 'ChEBI' in df.columns:
                    df['ChEBI'] = df['ChEBI'].astype(int)
                    df = df.set_index('ChEBI')
            else:
                # Handle different file formats
                try:
                    df = pd.read_csv(path, sep='\t', header=None, 
                                   names=['ChEBI', 'Info'], index_col='ChEBI',
                                   dtype={"ChEBI": "int", "Info": "str"})
                except ValueError:
                    # Try with header
                    df = pd.read_csv(path, sep='\t', index_col=0)
                    if len(df.columns) == 1:
                        df.columns = ['Info']
            
            id_to_info = df.to_dict()
            data[key] = id_to_info
            print(f"Loaded {len(df)} entries for property: {key}")
            
        except Exception as e:
            print(f"Error processing file {file}: {str(e)}")
            continue
    
    return data

def get_info(data, key, id):
    '''
    Enhanced function to retrieve information with better error handling
    '''
    try:
        info = data[key]['Info'][id]
        # Ensure we're returning a scalar value
        if isinstance(info, (list, dict, pd.Series)):
            info = str(info)
        # Handle numeric conversion for specific properties
        if key in ['Mass', 'logP', 'HBD', 'HBA', 'PSA', 'RotBonds']:
            info = pd.to_numeric(info, errors='coerce')
        return info
    except (KeyError, TypeError):
        return np.nan

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

def make_table(data, df_results, publication_dates=None, exclude_roles=None, include_roles=None, use_stop_list=True, stop_list_path=None):
    '''
    Enhanced function to create table with filtering options
    '''
    table = df_results.copy()
    
    # Apply stop list filter first if enabled
    if use_stop_list:
        stop_list = load_stop_list(stop_list_path)
        table = apply_stop_list_filter(table, stop_list)
    
    # Apply role-based filtering
    if exclude_roles or include_roles:
        table = filter_by_roles(table, data, exclude_roles, include_roles)
    
    # Create property columns
    for key in data.keys():
        try:
            column_list = []
            for id in table.index.to_list():
                column_list.append(get_info(data, key, id))
            
            table[key] = pd.Series(column_list, index=table.index)
            
            # Clean up specific numeric columns
            if key in ['Mass', 'logP', 'HBD', 'HBA', 'PSA', 'RotBonds']:
                table[key] = pd.to_numeric(table[key].replace('-', np.nan), errors='coerce')
                
        except Exception as e:
            print(f"Error processing column '{key}': {str(e)}")
            continue
    
    # Add publication years if available
    if publication_dates:
        pub_years = {}
        for chebi_id, row in df_results.iterrows():
            if chebi_id in table.index:  # Only for filtered results
                pub_ids = []
                if isinstance(df_results, pd.DataFrame) and 'Publication' in df_results.columns:
                    pub_ids = [row['Publication']]
                elif 'Publication' in df_results.reset_index().columns:
                    pub_id_series = df_results.reset_index().loc[df_results.reset_index()['ChEBI'] == chebi_id, 'Publication']
                    if not pub_id_series.empty:
                        pub_ids = pub_id_series.tolist()
                
                years = [int(publication_dates.get(pid)) for pid in pub_ids if pid in publication_dates and publication_dates.get(pid)]
                if years:
                    pub_years[chebi_id] = years
        
        if pub_years:
            table['PublicationYears'] = pd.Series({idx: pub_years.get(idx, []) for idx in table.index})
    
    # Calculate drug-likeness rules
    table = calculate_drug_likeness_rules(table)
    
    # Clean up table
    table = table.dropna(subset=['Count'])
    table = table.sort_values(by='Count', ascending=False)
    
    return table

def write_to_file(table, term):
    '''
    Enhanced function to write detailed output files with better error handling
    '''
    os.makedirs('tables', exist_ok=True)
    
    path = 'tables/'+term
    
    try:
        # Create detailed TSV output
        enhanced_table = table.copy()
        enhanced_table.reset_index(inplace=True)
        
        # Define column order for better readability
        base_columns = ['ChEBI', 'Count']
        
        # Add name column if available
        if 'Names' in enhanced_table.columns:
            base_columns.insert(1, 'Names')
        
        # Add core properties
        property_columns = ['Mass', 'logP']
        if 'TFIDF' in enhanced_table.columns:
            property_columns.append('TFIDF')
        
        # Add extended properties if available
        extended_props = ['HBD', 'HBA', 'PSA', 'RotBonds', 'Role', 'Application']
        extended_props = [col for col in extended_props if col in enhanced_table.columns]
        
        # Add drug-likeness columns
        rule_columns = ['Lipinski_RO5', 'Veber_Rule', 'Drug_Like']
        rule_columns = [col for col in rule_columns if col in enhanced_table.columns]
        
        # Add remaining columns
        all_ordered_cols = base_columns + property_columns + extended_props + rule_columns
        remaining_cols = [col for col in enhanced_table.columns if col not in all_ordered_cols]
        final_column_order = all_ordered_cols + remaining_cols
        
        # Write detailed TSV
        enhanced_table[final_column_order].to_csv(path+'_table_detailed.tsv', sep='\t', index=False)
        
        # Write standard pickle file
        table.to_pickle(path+'_table.pkl')
        
        # Write summary statistics
        summary_stats = {
            'total_compounds': len(table),
            'mean_mass': table['Mass'].mean() if 'Mass' in table.columns and not table['Mass'].isna().all() else None,
            'mean_logp': table['logP'].mean() if 'logP' in table.columns and not table['logP'].isna().all() else None,
            'drug_like_count': table['Drug_Like'].sum() if 'Drug_Like' in table.columns else None,
            'lipinski_compliant': table['Lipinski_RO5'].sum() if 'Lipinski_RO5' in table.columns else None
        }
        
        with open(path+'_summary.txt', 'w') as f:
            f.write(f"Summary for {term}\n")
            f.write("="*50 + "\n")
            for key, value in summary_stats.items():
                if value is not None:
                    if isinstance(value, float):
                        f.write(f"{key}: {value:.3f}\n")
                    else:
                        f.write(f"{key}: {value}\n")
                        
        print(f"Files written successfully: {path}_table_detailed.tsv, {path}_table.pkl, {path}_summary.txt")
        
    except Exception as e:
        print(f"Error writing files for {term}: {str(e)}")
        # Try to write at least the pickle file
        try:
            table.to_pickle(path+'_table.pkl')
            print(f"At least pickle file written: {path}_table.pkl")
        except:
            print(f"Failed to write any files for {term}")

def parser():
    parser = argparse.ArgumentParser(description='Enhanced table creation with filtering options')
    parser.add_argument('-i', required=True, metavar='input', dest='input', 
                       help='Input folder or file from the results folder')
    parser.add_argument('-t', required=True, metavar='type', dest='type', 
                       help='Type of input: file or folder')
    
    # Filtering options
    parser.add_argument('--exclude-roles', type=str, 
                       help='Comma-separated roles to exclude (e.g., "solvent,buffer")')
    parser.add_argument('--include-roles', type=str,
                       help='Comma-separated roles to include (e.g., "drug,metabolite")')
    parser.add_argument('--no-stop-list', action='store_true',
                       help='Disable default stop list filtering')
    parser.add_argument('--stop-list', type=str,
                       help='Path to custom stop list file')
    
    return parser.parse_args()

def main():
    args = parser()
    input_type = args.type
    input_path = args.input
    
    if input_type == 'file':
        results = [input_path]
    elif input_type == 'folder':
        if not os.path.exists(input_path):
            sys.exit(f'Error: Input folder {input_path} does not exist')
        files = os.listdir(input_path)
        results = [input_path+'/'+file for file in files if file.endswith('.tsv')]
    else:
        sys.exit('Error: please give \'file\' or \'folder\' as input type')

    if not results:
        sys.exit(f'Error: No TSV files found in {input_path}')

    # Gather properties
    print("Loading ChEBI properties...")
    data = import_properties()
    
    if not data:
        sys.exit("Error: No property data loaded. Check files directory.")
    
    print(f"Loaded properties: {list(data.keys())}")

    for result in results:
        try:
            # Extract term name from file path
            filename = os.path.basename(result)
            if '_ChEBI_IDs.tsv' in filename:
                term = filename.replace('_ChEBI_IDs.tsv', '')
            else:
                term = filename.replace('.tsv', '')
                
            print(f'\nMaking table for {term}')

            # Import publication dates if available
            publication_dates = import_publication_dates(term)

            # Import results
            df = pd.read_csv(result, sep='\t', names=['ChEBI', 'Publication'], 
                           dtype={"ChEBI": "int", "Publication": "str"})
            df_results = df.groupby(by=['ChEBI']).count().rename(columns={"Publication": "Count"})

            print(f"Initial compounds found: {len(df_results)}")

            # Make table with filtering
            table = make_table(data, df_results, publication_dates,
                             exclude_roles=args.exclude_roles,
                             include_roles=args.include_roles,
                             use_stop_list=not args.no_stop_list,
                             stop_list_path=args.stop_list)

            print(f"Final compounds after filtering: {len(table)}")

            # Check if 'idf' column exists before performing normalization
            if 'idf' in table.columns:
                print("Calculating TF-IDF values...")
                # Ensure numeric types for calculation
                table["Count"] = pd.to_numeric(table["Count"], errors='coerce')
                table["idf"] = pd.to_numeric(table["idf"], errors='coerce')
                
                # Perform normalization with proper handling of non-finite values
                tfidf_values = table["Count"].astype(float) * table["idf"].astype(float)
                
                # Handle non-finite values
                tfidf_values = tfidf_values.fillna(0)  # Replace NaN with 0
                tfidf_values = tfidf_values.replace([np.inf, -np.inf], 0)  # Replace inf with 0
                
                # Round and convert to int, but keep as float if there are decimal values
                table.loc[:,"TFIDF"] = tfidf_values.round(decimals=3)
                print("✓ TF-IDF values calculated successfully")
            else:
                print(f"Warning: 'idf' column not found for {term}, skipping TFIDF calculation")
            
            print("\nTop 10 compounds:")
            display_cols = ['Count']
            if 'Names' in table.columns:
                display_cols.insert(0, 'Names')
            if 'Mass' in table.columns:
                display_cols.append('Mass')
            if 'logP' in table.columns:
                display_cols.append('logP')
            if 'TFIDF' in table.columns:
                display_cols.append('TFIDF')
            if 'Drug_Like' in table.columns:
                display_cols.append('Drug_Like')
            
            print(table[display_cols].head(10))

            # Write table to file
            write_to_file(table, term)
            
        except Exception as e:
            print(f"Error processing {result}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\n✅ Table generation completed!")
    print(f"📁 Check the 'tables/' directory for output files")

if __name__ == '__main__':
    main()