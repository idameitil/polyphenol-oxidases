from ete3 import NCBITaxa
import pandas as pd

ncbi = NCBITaxa()
def get_taxon(taxid, rank):
    try:
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        taxon = ncbi.get_taxid_translator([ranks2lineage[rank]])[ranks2lineage[rank]]
    except:
        taxon = 'NA'
    return taxon

def read_interproscan_output(filename, acc2domain={}, unique_domains=set()):
    with open(filename, 'r') as infile:
        for line in infile:
            splitted_line = line.split('\t')
            protein_accession, md5, sequence_length, analysis, signature_accession, \
            signature_description, start_location, stop_location, score, status,\
            date, interpro_annotations_accession, interpro_annotations_description \
            = splitted_line
            domain_full_name = f"domain_{signature_accession}_{signature_description}"
            unique_domains.add(domain_full_name)
            if protein_accession not in acc2domain:
                acc2domain[protein_accession] = {domain_full_name: f"{start_location}-{stop_location}"}
            else:
                if domain_full_name not in acc2domain[protein_accession]:
                    acc2domain[protein_accession][domain_full_name] = f"{start_location}-{stop_location}"
                else:
                    acc2domain[protein_accession][domain_full_name] += f",{start_location}-{stop_location}"
    return acc2domain, unique_domains

def make_interproscan_df(acc2domain, unique_domains):
    data = {}
    for acc in acc2domain:
        data[acc] = []
        for domain_full_name in unique_domains:
            if domain_full_name not in acc2domain[acc]:
                data[acc].append('')
            else:
                data[acc].append(acc2domain[acc][domain_full_name])
    interproscan_df = pd.DataFrame.from_dict(data, orient='index', columns=list(unique_domains))
    interproscan_df.index.name = 'protein_accession'
    return interproscan_df