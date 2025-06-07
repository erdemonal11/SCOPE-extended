#!/usr/bin/python

import os
import requests
import json 
import shutil
import tqdm
import collections

# Additional property files to look for
EXTENDED_PROPERTY_FILES = [
    'ChEBI2Role_relation.tsv',
    'ChEBI2Application_relation.tsv', 
    'ChEBI2BiologicalRole_relation.tsv',
    'ChEBI2HBD_relation.tsv',      # Hydrogen Bond Donors
    'ChEBI2HBA_relation.tsv',      # Hydrogen Bond Acceptors  
    'ChEBI2PSA_relation.tsv',      # Polar Surface Area
    'ChEBI2RotBonds_relation.tsv', # Rotatable Bonds
    'ChEBI2SMILES_relation.tsv',   # SMILES strings
    'ChEBI2Pharmacology_relation.tsv',
    'ChEBI2Drug_relation.tsv'
]

def download_files(file_to_link, folder):
    '''
    This function downloads the missing files from the OSF project (https://osf.io/pvwu2/).
    '''
    print(f"Starting download of {len(file_to_link)} files...")
    
    for file, link in file_to_link.items():
        try:
            r = requests.get(link, stream=True)
            r.raise_for_status()  # Raise exception for bad status codes
            
            file_size = int(r.headers.get('Content-Length', 0))
            desc = f'Downloading {file}'
            path = os.path.join(folder, file)

            with tqdm.tqdm.wrapattr(r.raw, "read", total=file_size, desc=desc) as r_raw:
                with open(path, "wb") as f:
                    shutil.copyfileobj(r_raw, f)
            print(f"Successfully downloaded: {file}")
            
        except Exception as e:
            print(f"Error downloading {file}: {str(e)}")
    
    return

def get_response(url):
    '''
    This function returns url response with better error handling.
    '''
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        response_data = json.loads(response.text)
        return response_data
    except requests.exceptions.RequestException as e:
        print(f"Network error: {str(e)}")
        raise
    except json.JSONDecodeError as e:
        print(f"JSON decode error: {str(e)}")
        raise

def check_available_files(url):
    '''
    This function lists all available files in the OSF repository.
    Useful for identifying what property files are actually available.
    '''
    print("Checking available files in repository...")
    try:
        response = get_response(url)
        available_files = []
        
        for obj in response["data"]:
            file_name = obj["attributes"]["name"]
            file_size = obj["attributes"].get("size", "Unknown")
            available_files.append((file_name, file_size))
            
        print(f"Found {len(available_files)} files in repository:")
        for name, size in sorted(available_files):
            print(f"  {name} ({size} bytes)")
            
        return available_files
        
    except Exception as e:
        print(f"Error checking available files: {str(e)}")
        return []

def get_osf_to_rel(url):
    '''
    Enhanced function to map OSF files to releases with better error handling
    '''
    response = get_response(url)
    osf_to_rel = collections.defaultdict(dict)

    for obj in response["data"]:
        file_name = obj["attributes"]["name"]
        download_link = obj["links"]["download"]
        file_size = obj["attributes"].get("size", 0)

        if "_" in file_name:
            rel = file_name.split("_")[1].split(".")[0]
            file = file_name.split("_")[0]
        else:
            rel = ""
            file = file_name.split(".")[0]
            
        osf_to_rel[file]["rel"] = rel
        osf_to_rel[file]["link"] = download_link
        osf_to_rel[file]["name"] = file_name
        osf_to_rel[file]["size"] = file_size

    return osf_to_rel

def get_repo_to_rel(folder):
    '''
    Enhanced function to check existing files with better error handling
    '''
    if not os.path.exists(folder):
        return collections.defaultdict(dict)
        
    files = os.listdir(folder)
    repo_to_rel = collections.defaultdict(dict)

    for file_name in files:
        try:
            if "_" in file_name:
                rel = file_name.split("_")[1].split(".")[0]
                file = file_name.split("_")[0]
            else:
                rel = ""
                file = file_name.split(".")[0]
                
            path = os.path.join(folder, file_name)
            file_size = os.path.getsize(path)

            repo_to_rel[file]["rel"] = rel
            repo_to_rel[file]["path"] = path
            repo_to_rel[file]["size"] = file_size
        except Exception as e:
            print(f"Error processing existing file {file_name}: {str(e)}")
            
    return repo_to_rel

def get_files_to_download(osf_to_rel, repo_to_rel):
    '''
    Enhanced function to determine which files need downloading
    '''
    files_to_download = {}
    priority_files = []
    extended_files = []

    for file in osf_to_rel.keys():
        rel = osf_to_rel[file]['rel']
        link = osf_to_rel[file]['link']
        name = osf_to_rel[file]['name']
        size = osf_to_rel[file]['size']

        should_download = False
        reason = ""

        if file in repo_to_rel.keys():
            # File exists, check if it needs updating
            if rel != repo_to_rel[file]["rel"]:
                should_download = True
                reason = f"newer release ({rel} vs {repo_to_rel[file]['rel']})"
                # Remove old version
                old_path = repo_to_rel[file]["path"]
                try:
                    os.remove(old_path)
                    print(f"Removed outdated file: {old_path}")
                except Exception as e:
                    print(f"Warning: Could not remove old file {old_path}: {str(e)}")
            elif abs(size - repo_to_rel[file]["size"]) > 1000:  # Size difference > 1KB
                should_download = True
                reason = f"size mismatch ({size} vs {repo_to_rel[file]['size']} bytes)"
        else:
            # File doesn't exist
            should_download = True
            reason = "new file"

        if should_download:
            files_to_download[name] = link
            
            # Categorize files for reporting
            if any(ext_file in name for ext_file in EXTENDED_PROPERTY_FILES):
                extended_files.append((name, reason))
            else:
                priority_files.append((name, reason))

    # Report what will be downloaded
    if priority_files:
        print(f"\nCore files to download ({len(priority_files)}):")
        for name, reason in priority_files:
            print(f"  {name} - {reason}")
    
    if extended_files:
        print(f"\nExtended property files to download ({len(extended_files)}):")
        for name, reason in extended_files:
            print(f"  {name} - {reason}")
    
    if not files_to_download:
        print("\nAll files are up to date!")

    return files_to_download

def validate_downloaded_files(folder):
    '''
    Validate that critical files were downloaded successfully
    '''
    if not os.path.exists(folder):
        print(f"Warning: Download folder {folder} does not exist")
        return False
    
    files = os.listdir(folder)
    
    # Check for core required files
    core_files = ['ChEBI2Mass', 'ChEBI2logP', 'ChEBI2Names']
    missing_core = []
    
    for core_file in core_files:
        if not any(core_file in f for f in files):
            missing_core.append(core_file)
    
    if missing_core:
        print(f"Warning: Missing core files: {missing_core}")
        return False
    
    # Report on extended files
    extended_found = []
    for ext_file in EXTENDED_PROPERTY_FILES:
        if any(ext_file.replace('_relation.tsv', '') in f for f in files):
            extended_found.append(ext_file)
    
    print(f"\nDownload validation:")
    print(f"  Core files: {'✓' if not missing_core else '✗'}")
    print(f"  Extended properties found: {len(extended_found)}")
    
    if extended_found:
        print("  Available extended properties:")
        for ext_file in extended_found:
            print(f"    - {ext_file}")
    
    return len(missing_core) == 0

def main():
    # Define constants
    folder = 'files'
    url = 'https://api.osf.io/v2/nodes/pvwu2/files/osfstorage/611252ba847d1304ca38b4d4/'

    print("SCOPE-Extended Enhanced File Downloader")
    print("="*50)
    
    # Create directory if it doesn't exist
    if not os.path.isdir(folder):
        os.mkdir(folder)
        print(f"Created directory: {folder}")

    # Check what files are available (optional, for debugging)
    # available_files = check_available_files(url)

    # Get file mappings
    print("\nAnalyzing repository...")
    osf_to_rel = get_osf_to_rel(url)
    repo_to_rel = get_repo_to_rel(folder)

    print(f"Found {len(osf_to_rel)} files in repository")
    print(f"Found {len(repo_to_rel)} existing local files")

    # Determine what needs downloading
    files_to_download = get_files_to_download(osf_to_rel, repo_to_rel)

    if files_to_download:
        print(f"\nDownloading {len(files_to_download)} files...")
        download_files(files_to_download, folder)
    
    # Validate download
    print("\nValidating downloads...")
    validation_success = validate_downloaded_files(folder)
    
    if validation_success:
        print("\n✓ Download completed successfully!")
        print(f"Files are ready in the '{folder}' directory")
        print("\nNext steps:")
        print("1. Run search_query.py to find publications")
        print("2. Run make_table.py with enhanced filtering options")
        print("3. Run visualize_multiplot.py for interactive analysis")
    else:
        print("\n⚠ Download completed with warnings")
        print("Some files may be missing. The basic functionality should still work.")

if __name__ == '__main__':
    main()