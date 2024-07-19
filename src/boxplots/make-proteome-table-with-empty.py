import sys
import json
import pandas as pd

json_list = json.load(open("data/pfam/proteome-matching-PF00264.json"))


def get_taxon(rank, genome):
    if "taxonLineage" in genome:
        for entry in genome["taxonLineage"]:
            if entry["rank"] == rank:
                return entry["scientificName"]
    return ""


wanted_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
metadata_json = json.load(
    open("data/proteome-tree/proteomes_AND_proteome_type_1_2024_02_28.json")
)
metadata_dict = {}
for genome in metadata_json["results"]:
    if "buscoReport" in genome["proteomeCompletenessReport"]:
        if genome["proteomeCompletenessReport"]["buscoReport"]["score"] < 75:
            continue
        dict1 = genome["proteomeCompletenessReport"]["buscoReport"]
    else:
        continue
    for rank in wanted_ranks:
        taxon = get_taxon(rank, genome)
        dict1[rank] = taxon
    metadata_dict[genome["id"]] = dict1
metadata_df = pd.DataFrame.from_dict(metadata_dict, orient="index")
metadata_df.index.name = "accession"
metadata_df.to_csv("test2.tsv", sep="\t")

data = {}
for entry in json_list:
    acc = entry["metadata"]["accession"]
    name = entry["metadata"]["name"]
    count_tyrosinases = entry["extra_fields"]["counters"]["proteins"]
    data[acc] = {"name": name, "count_tyrosinases": count_tyrosinases}

df_tyrosinase_count = pd.DataFrame.from_dict(data, orient="index")
df_tyrosinase_count.index.name = "accession"

df_merged = pd.merge(
    left=metadata_df, right=df_tyrosinase_count, how="left", on="accession"
)
df_sorted = df_merged.sort_values(
    ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
)

df_sorted.to_csv("data/boxplots/boxplot-data.tsv", sep="\t")
