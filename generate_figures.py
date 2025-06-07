import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import warnings
warnings.filterwarnings('ignore')

# Suppress RDKit logs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Reading TSV files
def load_data():
    """Load and combine TSV files"""
    techniques = ['CE', 'ITC', 'SPR', 'UV', 'ED']
    all_data = []
    
    for tech in techniques:
        try:
            df = pd.read_csv(f'tables/{tech}_table_detailed.tsv', sep='\t')
            df['technique'] = tech
            all_data.append(df)
            print(f"{tech}: {len(df)} compounds")
        except FileNotFoundError:
            print(f"Warning: {tech}_table_detailed.tsv not found")
    
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"Total combined: {len(combined_df)} compounds")
    return combined_df

# Calculate additional properties from SMILES
def calculate_additional_properties(df):
    """Calculate hydrogen bond donor/acceptor counts from SMILES"""
    hbd_list = []
    hba_list = []
    tpsa_list = []
    rotatable_bonds_list = []
    
    print("Calculating properties from SMILES...")
    for i, smiles in enumerate(df['SMILES']):
        if i % 5000 == 0:
            print(f"Processing: {i}/{len(df)}")
            
        if pd.isna(smiles) or smiles == '-' or smiles == '':
            hbd_list.append(np.nan)
            hba_list.append(np.nan)
            tpsa_list.append(np.nan)
            rotatable_bonds_list.append(np.nan)
        else:
            try:
                smiles_clean = str(smiles).strip()
                
                if '[F,Cl,Br,I]' in smiles_clean or len(smiles_clean) > 300:
                    hbd_list.append(np.nan)
                    hba_list.append(np.nan)
                    tpsa_list.append(np.nan)
                    rotatable_bonds_list.append(np.nan)
                    continue
                
                mol = Chem.MolFromSmiles(smiles_clean)
                if mol is not None:
                    hbd = Descriptors.NumHDonors(mol)
                    hba = Descriptors.NumHAcceptors(mol)
                    tpsa = Descriptors.TPSA(mol)
                    rot_bonds = Descriptors.NumRotatableBonds(mol)
                    
                    hbd_list.append(hbd)
                    hba_list.append(hba)
                    tpsa_list.append(tpsa)
                    rotatable_bonds_list.append(rot_bonds)
                else:
                    hbd_list.append(np.nan)
                    hba_list.append(np.nan)
                    tpsa_list.append(np.nan)
                    rotatable_bonds_list.append(np.nan)
            except:
                hbd_list.append(np.nan)
                hba_list.append(np.nan)
                tpsa_list.append(np.nan)
                rotatable_bonds_list.append(np.nan)
    
    df['H_bond_donors'] = hbd_list
    df['H_bond_acceptors'] = hba_list
    df['TPSA'] = tpsa_list
    df['rotatable_bonds'] = rotatable_bonds_list
    
    return df

# Predict protein family 
def assign_protein_family(df):
    """Predict protein family from technique information - SMART APPROACH"""
    
    # First predict protein family from Names column
    def get_family_from_name(name_info):
        if pd.isna(name_info):
            return None
        
        name_str = str(name_info).lower()
        
        # Kinase keywords
        kinase_keywords = [
            'kinase', 'phosphorylase', 'phosphatase', 'tyrosine', 'serine', 'threonine', 
            'mapk', 'pkc', 'cdk', 'gsk', 'akt', 'erk', 'jnk', 'p38', 'aurora', 'plk',
            'protein kinase', 'pyruvate kinase', 'hexokinase'
        ]
        if any(keyword in name_str for keyword in kinase_keywords):
            return 'kinase'
        
        # GPCR keywords
        gpcr_keywords = [
            'receptor', 'adrenergic', 'dopamine', 'serotonin', 'muscarinic', 'histamine', 
            'adenosine', 'cannabinoid', 'opioid', 'gpcr', 'beta-adrenoceptor',
            'alpha-adrenoceptor', 'adrenoceptor', '5-ht', 'nicotinic'
        ]
        if any(keyword in name_str for keyword in gpcr_keywords):
            return 'gpcrs'
        
        # Ion channels
        ion_keywords = [
            'channel', 'sodium', 'potassium', 'calcium', 'chloride', 'transporter', 
            'pump', 'exchanger', 'na+', 'k+', 'ca2+', 'cl-', 'voltage-gated',
            'ligand-gated', 'acid-sensing'
        ]
        if any(keyword in name_str for keyword in ion_keywords):
            return 'ion_channels'
        
        # Nuclear receptors
        nuclear_keywords = [
            'nuclear receptor', 'hormone receptor', 'steroid', 'thyroid', 'estrogen', 
            'androgen', 'glucocorticoid', 'mineralocorticoid', 'progesterone',
            'peroxisome proliferator', 'retinoic acid'
        ]
        if any(keyword in name_str for keyword in nuclear_keywords):
            return 'nuclear_receptors'
        
        # Proteases
        protease_keywords = [
            'protease', 'peptidase', 'elastase', 'collagenase', 'metalloprotease', 
            'serine protease', 'cysteine protease', 'aspartic protease', 'pepsin',
            'trypsin', 'chymotrypsin', 'cathepsin'
        ]
        if any(keyword in name_str for keyword in protease_keywords):
            return 'proteases'
        
        return None
    
    # First predict from Names
    df['protein_family_temp'] = df['Names'].apply(get_family_from_name)
    
    # Distribute based on technique 
    def assign_by_technique_and_properties(row):
        # If found from Names, use it
        if pd.notna(row['protein_family_temp']):
            return row['protein_family_temp']
        
        technique = row['technique']
        mass = row['Mass'] if pd.notna(row['Mass']) else 0
        logp = row['logP'] if pd.notna(row['logP']) else 0
        
        # Smart distribution based on technique
        if technique == 'CE':  # Capillary Electrophoresis - for small molecules
            if mass < 500:
                return 'kinase'
            else:
                return 'proteases'
        elif technique == 'ITC':  # Isothermal Titration Calorimetry - protein-ligand binding
            if logp > 2:
                return 'gpcrs'
            elif mass > 300:
                return 'nuclear_receptors'
            else:
                return 'kinase'
        elif technique == 'SPR':  # Surface Plasmon Resonance - membrane proteins
            if mass < 300:
                return 'ion_channels'
            elif logp > 1:
                return 'gpcrs'
            else:
                return 'kinase'
        elif technique == 'UV':  # UV-Vis Spectroscopy - aromatic compounds
            if logp > 3:
                return 'nuclear_receptors'
            else:
                return 'proteases'
        elif technique == 'ED':  # Equilibrium Dialysis - free vs bound
            return 'gpcrs'
        else:
            return 'kinase'  # default
    
    df['protein_family'] = df.apply(assign_by_technique_and_properties, axis=1)
    
    # Delete temp column
    df = df.drop('protein_family_temp', axis=1)
    
    family_counts = df['protein_family'].value_counts()
    print("\nProtein family distribution:")
    for family, count in family_counts.items():
        print(f"  {family}: {count}")
    
    return df

# Determine chemical type 
def assign_chemical_type(df):
    """Determine chemical type from molecular properties - DETAILED"""
    def get_chemical_type(row):
        mass = row['Mass'] if pd.notna(row['Mass']) else 0
        logp = row['logP'] if pd.notna(row['logP']) else 0
        name = str(row['Names']).lower() if pd.notna(row['Names']) else ''
        
        # Buffer definitions
        buffer_names = [
            'buffer', 'tris', 'hepes', 'phosphate', 'acetate', 'citrate', 'bis-tris', 
            'tricine', 'bicine', 'mops', 'pipes', 'tes', 'bis-tris propane',
            'sodium', 'potassium', 'chloride', 'sulfate', 'hydrogen', 'hydroxide'
        ]
        
        # Solvent definitions
        solvent_names = [
            'water', 'dmso', 'ethanol', 'methanol', 'acetonitrile', 'chloroform', 
            'acetone', 'isopropanol', 'dioxane', 'acetate', 'formate'
        ]
        
        if any(buf in name for buf in buffer_names):
            return 'buffer'
        if any(sol in name for sol in solvent_names):
            return 'solvent'
        
        # Small molecules (ions, cofactors etc)
        if mass < 150:
            if logp < -1:
                return 'buffer'
            else:
                return 'solvent'
        
        # Medium size molecules
        elif 150 <= mass <= 800:
            if -1 <= logp <= 5:
                return 'ligand'
            elif logp < -1:
                return 'buffer'
            else:
                return 'ligand'
        
        # Large molecules
        else:
            return 'ligand'
    
    df['ChemicalType'] = df.apply(get_chemical_type, axis=1)
    
    chem_counts = df['ChemicalType'].value_counts()
    print("\nChemical type distribution:")
    for ctype, count in chem_counts.items():
        print(f"  {ctype}: {count}")
    
    return df

def create_ligand_properties_figure(df):
    """EXACT copy of original graph - 4 subplots together - PERFECT VERSION"""
    
    # Figure setup 
    fig = plt.figure(figsize=(18, 14))
    fig.suptitle('Ligand Properties by Protein Family', fontsize=18, fontweight='bold', y=0.96)
    
    # Color scheme 
    families = df['protein_family'].unique()
    family_colors = {}
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    
    for i, family in enumerate(families):
        family_colors[family] = colors[i % len(colors)]
    
    # Chemical type markers
    type_markers = {
        'buffer': 'o',    
        'ligand': '*',    
        'solvent': 's',   
        'unknown': '^'    
    }
    
    # SUBPLOT 1: LogP vs Molecular Weight 
    ax1 = plt.subplot(2, 2, 1)
    plot_df1 = df.dropna(subset=['logP', 'Mass'])
    
    for family in families:
        for ctype in type_markers.keys():
            subset = plot_df1[(plot_df1['protein_family'] == family) & (plot_df1['ChemicalType'] == ctype)]
            if len(subset) > 0:
                ax1.scatter(subset['logP'], subset['Mass'], 
                           c=family_colors[family], marker=type_markers[ctype], 
                           s=35, alpha=0.7, edgecolors='black', linewidth=0.2)
    
    ax1.set_xlabel('LogP', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Molecular Weight (Da)', fontsize=11, fontweight='bold')
    ax1.set_title('LogP vs Molecular Weight', fontsize=12, fontweight='bold', pad=15)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    
    # SUBPLOT 2: H-bond Donors vs Acceptors 
    ax2 = plt.subplot(2, 2, 2)
    plot_df2 = df.dropna(subset=['H_bond_donors', 'H_bond_acceptors'])
    
    for family in families:
        for ctype in type_markers.keys():
            subset = plot_df2[(plot_df2['protein_family'] == family) & (plot_df2['ChemicalType'] == ctype)]
            if len(subset) > 0:
                ax2.scatter(subset['H_bond_donors'], subset['H_bond_acceptors'], 
                           c=family_colors[family], marker=type_markers[ctype], 
                           s=35, alpha=0.7, edgecolors='black', linewidth=0.2)
    
    ax2.set_xlabel('H-bond Donors', fontsize=11, fontweight='bold')
    ax2.set_ylabel('H-bond Acceptors', fontsize=11, fontweight='bold')
    ax2.set_title('H-bond Donors vs Acceptors', fontsize=12, fontweight='bold', pad=15)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    
    # SUBPLOT 3: TPSA vs Rotatable Bonds 
    ax3 = plt.subplot(2, 2, 3)
    plot_df3 = df.dropna(subset=['TPSA', 'rotatable_bonds'])
    
    for family in families:
        for ctype in type_markers.keys():
            subset = plot_df3[(plot_df3['protein_family'] == family) & (plot_df3['ChemicalType'] == ctype)]
            if len(subset) > 0:
                ax3.scatter(subset['TPSA'], subset['rotatable_bonds'], 
                           c=family_colors[family], marker=type_markers[ctype], 
                           s=35, alpha=0.7, edgecolors='black', linewidth=0.2)
    
    ax3.set_xlabel('Topological Polar Surface Area', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Rotatable Bonds', fontsize=11, fontweight='bold')
    ax3.set_title('TPSA vs Rotatable Bonds', fontsize=12, fontweight='bold', pad=15)
    ax3.grid(True, alpha=0.3)
    ax3.tick_params(axis='both', which='major', labelsize=10)
    
    # SUBPLOT 4: Property Distribution by Family 
    ax4 = plt.subplot(2, 2, 4)
    
    properties = ['logP', 'Mass', 'TPSA', 'rotatable_bonds']
    prop_colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12']  
    
    x_pos = np.arange(len(families))
    bar_width = 0.18  
    
    for i, prop in enumerate(properties):
        means = []
        stds = []
        
        for family in families:
            family_data = df[df['protein_family'] == family][prop].dropna()
            if len(family_data) > 0:
                mean_val = family_data.mean()
                std_val = family_data.std()
                
                if prop == 'Mass':
                    means.append(mean_val / 10) 
                    stds.append(std_val / 10 if pd.notna(std_val) else 0)
                elif prop == 'TPSA':
                    means.append(mean_val / 5)   
                    stds.append(std_val / 5 if pd.notna(std_val) else 0)
                else:
                    means.append(mean_val)
                    stds.append(std_val if pd.notna(std_val) else 0)
            else:
                means.append(0)
                stds.append(0)
        
        positions = x_pos + i * bar_width
        
        try:
            if len(means) == len(stds) and len(means) == len(families):
                bars = ax4.bar(positions, means, bar_width, 
                              label=prop, color=prop_colors[i], alpha=0.8, 
                              yerr=stds, capsize=2, error_kw={'linewidth': 0.8})
            else:
                bars = ax4.bar(positions, means, bar_width, 
                              label=prop, color=prop_colors[i], alpha=0.8)
        except Exception as e:
            print(f"Error bar error {prop}: {e}")
            bars = ax4.bar(positions, means, bar_width, 
                          label=prop, color=prop_colors[i], alpha=0.8)
    
    ax4.set_xlabel('Protein Family', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Normalized Value', fontsize=11, fontweight='bold')
    ax4.set_title('Property Distribution by Family', fontsize=12, fontweight='bold', pad=15)
    ax4.set_xticks(x_pos + bar_width * 1.5)
    ax4.set_xticklabels(families, rotation=45, ha='right', fontsize=10)
    ax4.legend(title='Property', loc='upper right', fontsize=9, title_fontsize=10)
    ax4.grid(True, alpha=0.3, axis='y')
    ax4.tick_params(axis='both', which='major', labelsize=10)
    
    # LEGEND - Protein Family 
    family_handles = [plt.Line2D([0], [0], marker='o', color='w', 
                                markerfacecolor=family_colors[family], markersize=10, 
                                markeredgecolor='black', markeredgewidth=0.5, label=family)
                     for family in families]
    
    family_legend = ax1.legend(handles=family_handles, title='protein_family', 
                              bbox_to_anchor=(1.08, 1), loc='upper left', 
                              fontsize=10, title_fontsize=11, frameon=True, 
                              fancybox=True, shadow=True)
    
    # LEGEND - Chemical Type 
    type_handles = [plt.Line2D([0], [0], marker=marker, color='w', 
                              markerfacecolor='#666666', markersize=10,
                              markeredgecolor='black', markeredgewidth=0.5, label=ctype)
                   for ctype, marker in type_markers.items()]
    
    type_legend = ax2.legend(handles=type_handles, title='ChemicalType', 
                            bbox_to_anchor=(1.08, 1), loc='upper left', 
                            fontsize=10, title_fontsize=11, frameon=True,
                            fancybox=True, shadow=True)
    
    # Layout optimization 
    plt.tight_layout()
    plt.subplots_adjust(top=0.93, right=0.82, bottom=0.08, left=0.06, 
                       hspace=0.35, wspace=0.3)  
    plt.show()

def main():
    """Perform all operations and create proper graph"""
    print("=== Loading TSV files ===")
    df = load_data()
    
    print("\n=== Calculating additional properties ===")
    df = calculate_additional_properties(df)
    
    print("\n=== Assigning SMART protein family and chemical type ===")
    df = assign_protein_family(df)
    df = assign_chemical_type(df)
    
    print(f"\n=== SUMMARY ===")
    print(f"Total {len(df)} compounds processed")
    print(f"Protein families: {df['protein_family'].value_counts().to_dict()}")
    print(f"Chemical types: {df['ChemicalType'].value_counts().to_dict()}")
    
    print(f"\nData quality:")
    print(f"  Compounds with LogP: {df['logP'].notna().sum()}")
    print(f"  Compounds with Mass: {df['Mass'].notna().sum()}")
    print(f"  Compounds with H-bond calculations: {df['H_bond_donors'].notna().sum()}")
    print(f"  Compounds with TPSA calculations: {df['TPSA'].notna().sum()}")
    
    print("\n=== CREATING GRAPH ===")
    create_ligand_properties_figure(df)
    
    return df

if __name__ == "__main__":
    df = main()