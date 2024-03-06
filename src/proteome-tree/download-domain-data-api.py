import json
import requests
from multiprocessing import Pool

def get_data(api_url):
    response = requests.get(api_url)
     
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Error: {response.status_code}")
 
def url_signalp_phobius(accession):
    return f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/{accession}?extra_features=&amp;format=json"

def url_pfam(accession):
    return f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{accession}?format=json"

def save_data(protein):
    protein_accession = protein['metadata']['accession']
    # Get data
    data1 = get_data(url_signalp_phobius(protein_accession))
    data2 = get_data(url_pfam(protein_accession))
    # Save data
    with open(f"data/pfam/api/{protein_accession}_data1.json", 'w') as outfile:
        json.dump(data1, outfile)
    with open(f"data/pfam/api/{protein_accession}_data2.json", 'w') as outfile:
        json.dump(data2, outfile)

if __name__ == '__main__':
    json_object = json.load(open('data/proteome-tree/export.json'))
    print(len(json_object))
    with(Pool(7)) as p:
        p.map(save_data, json_object)
