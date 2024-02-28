import json
from common import get_taxon
import pandas as pd

json_list = json.load(open('data/pfam/proteome-matching-PF00264.json'))

def Merge(dict1, dict2):
    res = {**dict1, **dict2}
    return res

metadata_json = json.load(open('data/proteome-tree/proteomes_AND_proteome_type_1_2024_02_28.json'))
metadata_dict = {}
for genome in metadata_json['results']:
    if 'buscoReport' in genome['proteomeCompletenessReport']:
        dict1 = genome['proteomeCompletenessReport']['buscoReport']
    else:
        dict1 = {'complete': '', 'completeSingle': '', 'completeDuplicated': '', 'fragmented': '', 'missing': '', 'total': '', 'lineageDb': '', 'score': ''}
    dict2 = genome['proteomeCompletenessReport']['cpdReport']
    merged = Merge(dict1, dict2)
    metadata_dict[genome['id']] = merged

data = {}
for entry in json_list:
    taxid = entry['metadata']['taxonomy']
    acc = entry['metadata']['accession']
    name = entry['metadata']['name']
    count_tyrosinases = entry['extra_fields']['counters']['proteins']
    kingdom = get_taxon(taxid, 'kingdom')
    phylum = get_taxon(taxid, 'phylum')
    class_tax = get_taxon(taxid, 'class')
    order = get_taxon(taxid, 'order')
    family = get_taxon(taxid, 'family')
    genus = get_taxon(taxid, 'genus')
    species = get_taxon(taxid, 'species')
    
    dict1 = {'taxid': taxid, 'name': name, 'count_tyrosinases': count_tyrosinases,
                 'kingdom': kingdom, 'phylum': phylum, 'class': class_tax, 'order': order,
                'family':family, 'genus':genus, 'species': species}
    
    data[acc] = Merge(dict1, metadata_dict[acc])

    
df = pd.DataFrame.from_dict(data).T.sort_values(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

df[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', \
    'count_tyrosinases', 'taxid', 'name', 'complete', 'completeSingle', \
    'completeDuplicated', 'fragmented', 'missing', 'total', 'lineageDb', 'score', \
    'proteomeCount', 'stdCdss', 'averageCdss', 'confidence', 'status'\
    ]].to_csv('data/proteome-tree/proteome-data.tsv', sep='\t')