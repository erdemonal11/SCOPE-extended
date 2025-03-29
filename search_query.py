#!/usr/bin/python

import argparse
import urllib.request
import urllib.parse
import json
import concurrent.futures
import csv
import time
import datetime
import math
import os

def construct_url(query, pageSize, cursorMark):
    """
    This function constructs a url that searches through the europepmc website and returns the results in json format.
    """

    query_string = urllib.parse.quote(query)

    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + query_string + "&pageSize=" + str(pageSize) + "&resultType=lite&format=JSON&cursorMark=" + str(cursorMark)

    return url


def get_data(url, timeout, sleep_time):
    """
    This function recieves a url and returns the json data as type dictionary.
    If there is no response within the given timeout, an error message is printed and the connection is retried after the given sleep time.
    """
    connection = False
    while connection == False:
        try:
            response = urllib.request.urlopen(url, timeout=timeout)
            data = json.loads(response.read())
            connection = True
        except:
            print('connection failed')
            time.sleep(sleep_time) # in seconds
    return data


def find_publications_with_tmt(data):
    """
    This function looks for publications with text mined terms.
    It returns a list of tuples with the publication ID as the first value, the publication source as the second value,
    and the publication date as the third value.
    """

    page = []
    for publication in data['resultList']['result']:
        if publication['hasTextMinedTerms'] == 'Y':
            # Extract publication date if available
            pub_date = publication.get('pubYear', None) or publication.get('firstPublicationDate', None)
            publication = (publication['id'], publication['source'], pub_date)
            page.append(publication)

    return page


def search_publications(query, pageSize):
    """
    This function searches the europe pmc site with the query and retrieves all hits (publications).
    Variable 'CursorMark' is used to go through the search result pages until all publcations are retrieved.
    Publications are selected for having text mined terms with function 'find_publications_with_tmt'.
    Their annotations are downloaded with function 'get_annotations' and saved in a dictionary with publication id as key and chemical ids in a list as value.
    After the search is done, this dictionary is returned.
    """
    chebi_dict = dict()
    pub_date_dict = dict()  # Dictionary to store publication dates
    TIMEOUT=60
    SLEEP_TIME=5

    print("searching publications...")
    cursorMark = "*"
    pages = []
    url = construct_url(query, pageSize, cursorMark)
    query_data = get_data(url, TIMEOUT, SLEEP_TIME)
    total_hits = query_data['hitCount']
    print("total hits: %d" % total_hits)

    try:
        nextCursorMark = query_data['nextCursorMark']
    except:
        nextCursorMark = None
    publications = find_publications_with_tmt(query_data)
    chebi_dict, pub_date_dict = get_annotations(publications, chebi_dict, pub_date_dict)

    counter = 1
    while nextCursorMark != None:

        print('%d/%d pages retrieved' % (counter, math.ceil(total_hits/1000)))

        cursorMark = nextCursorMark
        url = construct_url(query, pageSize, cursorMark)
        query_data = get_data(url, TIMEOUT, SLEEP_TIME)
        try:
            nextCursorMark = query_data['nextCursorMark']
        except:
            nextCursorMark = None
        publications = find_publications_with_tmt(query_data)
        chebi_dict, pub_date_dict = get_annotations(publications, chebi_dict, pub_date_dict)

        counter += 1

    return chebi_dict, pub_date_dict


def get_annotations(publications, dict, pub_date_dict):
    """
    This function searches through the publications with text mined terms for annotations of type 'Chemicals'.
    From the ChEBI urls, the ChEBI ID's are extracted and returned as values with the publication ID's as keys in a dictionary.
    Also stores publication dates in a separate dictionary.
    Now includes improved error handling and progress reporting.
    """
    if not publications:
        print("No publications to process")
        return dict, pub_date_dict

    # Construct urls to the annotation data for every publication
    print(f"Processing {len(publications)} publications...")
    urls = []
    for publication in publications:
        id = publication[0]
        source = publication[1]
        pub_date = publication[2]
        
        # Store publication date
        if pub_date:
            pub_date_dict[id] = pub_date
            
        url = "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds=" + str(source) + ":" + str(id) + "&format=JSON"
        urls.append((url, id, source))

    json_data = []
    CONNECTIONS = 5  # Reduced from 10 to put less strain on the API
    TIMEOUT = 15     # Increased timeout
    SLEEP_TIME = 5   # Reduced sleep time for retries

    # Process in smaller batches to avoid overwhelming the API
    batch_size = 50
    total_batches = (len(urls) + batch_size - 1) // batch_size
    
    for batch_idx in range(total_batches):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, len(urls))
        batch_urls = urls[start_idx:end_idx]
        
        print(f"Processing batch {batch_idx+1}/{total_batches} ({start_idx+1}-{end_idx} of {len(urls)})")
        
        # ThreadPoolExecutor allows multiple connections, thereby speeding up the process of downloading annotations
        with concurrent.futures.ThreadPoolExecutor(max_workers=CONNECTIONS) as executor:
            futures = []
            for url_info in batch_urls:
                url, id, source = url_info
                futures.append(executor.submit(get_data, url, TIMEOUT, SLEEP_TIME))
            
            completed = 0
            for future in concurrent.futures.as_completed(futures):
                try:
                    data = future.result()
                    if data:  # Only append if we got valid data
                        json_data.append(data)
                    completed += 1
                    if completed % 10 == 0:
                        print(f"Completed {completed}/{len(batch_urls)} requests in current batch")
                except Exception as e:
                    print(f"Error processing annotation request: {str(e)}")
        
        # Add a small delay between batches
        if batch_idx < total_batches - 1:
            print(f"Batch complete. Pausing for 2 seconds before next batch...")
            time.sleep(2)

    print(f"Retrieved annotation data for {len(json_data)} publications")
    
    # Select the downloaded annotations for type 'Chemicals' and extract ChEBI ID
    chemicals_found = 0
    for data in json_data:
        if not data or len(data) == 0:
            continue
            
        try:
            annotations = data[0]['annotations']
            pub_id = data[0]['extId']
            
            for annotation in annotations:
                if annotation['type'] == 'Chemicals':
                    chebi_url = annotation['tags'][0]['uri']
                    if '_' in chebi_url:  # Make sure we can split the URL properly
                        chebi_id = chebi_url.split('_')[1]
                        chemicals_found += 1
                        
                        # Check if publication ID is already in dictionary, if not, then a list as value needs to be created for potentially multiple chebi id's per publication
                        try:
                            dict[pub_id]
                            dict[pub_id].append(chebi_id)
                        except:
                            dict[pub_id] = [chebi_id]
        except Exception as e:
            print(f"Error processing annotation data: {str(e)}")
            
    print(f"Found {chemicals_found} chemical annotations across {len(dict)} publications")
    return dict, pub_date_dict

def write_results(dict, pub_date_dict, term, query):
    """
    This function writes the ChEBI urls and publication ID's in a seperate csv file and the metadata to a text file.
    Also writes publication dates to a separate file.
    """
    # Create directories if they don't exist
    if not os.path.exists('results'):
        os.makedirs('results')
    if not os.path.exists('metadata'):
        os.makedirs('metadata')
    if not os.path.exists('publication_dates'):
        os.makedirs('publication_dates')

    # Write ChEBI IDs and publication IDs
    file = 'results/'+str(term)+'_ChEBI_IDs.tsv'
    count = 0
    uniques = set()

    with open(file, 'w', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        for pub_id in dict.keys():
            for chebi_id in dict[pub_id]:
                count += 1
                uniques.add(chebi_id)
                writer.writerow([chebi_id, pub_id])

    # Write publication dates
    date_file = 'publication_dates/'+str(term)+'_pub_dates.tsv'
    with open(date_file, 'w', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        writer.writerow(['PublicationID', 'PublicationYear'])
        for pub_id, date in pub_date_dict.items():
            writer.writerow([pub_id, date])

    print('%s query results are written to file' % term)

    # Get metadata
    current_day = datetime.date.today()
    number_of_papers = len(dict.keys())
    number_of_chemicals = count
    number_of_unique_chemicals = len(uniques)
    
    # Publication date statistics
    years = [int(y) for y in pub_date_dict.values() if y and y.isdigit()]
    year_stats = ""
    if years:
        earliest_year = min(years)
        latest_year = max(years)
        year_stats = f"\npublication years range: {earliest_year} - {latest_year}"

    file = 'metadata/'+str(term)+'.txt'
    f = open(file, 'w')
    f.write('metadata for %s\n' % term
    + 'query: %s\n' % query
    + 'search date: %s\n' % current_day
    + 'number of papers: %d\n' % number_of_papers
    + 'number of chemicals: %d\n' % number_of_chemicals
    + 'number of unique chemicals: %d (note: not all chemicals can be plotted due to missing logP values)' % number_of_unique_chemicals
    + year_stats)
    f.close()

def read_input(file):
    '''
    This function reads the input file and returns the query terms in a dictionary.
    Every line in the input file should be a new query.
    In the lines, synonyms are sperated by ", ".
    The dictionary values are lists of the terms, including the first tirm that is used as a key.
    '''
    f = open(file, 'r')
    input = f.readlines()
    queries = dict()

    for line in input:
        query_list = line.split(',', 1)
        term = query_list[0]
        query = query_list[1].strip()
        queries[term] = query

    return queries

def parser():
    parser = argparse.ArgumentParser(description='A script that parses JSON data from EPC publications and extracts annotations of type Chemicals')
    parser.add_argument('-i', required=True, metavar='input_file', dest='input_file', help='[i] to select input file from the queries folder')
    parser.add_argument('-p', required=False, metavar='pageSize', dest='pageSize', help='[p] to select pageSize, not required, defeault=1000')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parser()
    input_file = args.input_file
    #pageSize = int(args.pageSize)
    pageSize = 1000

    queries = read_input(input_file)
    for term in queries.keys():
        query = queries[term]
        print('searching with: %s' % query)
        chebi_dict, pub_date_dict = search_publications(query, pageSize)
        write_results(chebi_dict, pub_date_dict, term, query)
        print('%d publications with text mined terms and annotations of type \'chemical\' found for %s' % (len(chebi_dict.keys()), term) )

if __name__ == '__main__':
    main()