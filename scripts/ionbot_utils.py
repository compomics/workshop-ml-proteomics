import re
import pandas as pd
from collections import Counter


def get_ionbot_condition(condition,folder,dataset, features=False):
    #read consition results
    ionbot = pd.read_csv("%s/%s%s/ionbot.first.csv"%(folder,dataset,condition))
    ionbot["rank"] = ["first"]*len(ionbot)
    if features:
        tmp = pd.read_csv("%s/%s%s/ionbot.features.csv"%(folder,dataset,condition))
        ionbot = ionbot.merge(tmp,on="ionbot_match_id",how="left")

    ionbot["PSM"] = ionbot["scan"].astype(str) + ionbot["spectrum_file"] +  "|PSM|" +  ionbot["matched_peptide"]
    ionbot["num_matches"] = ionbot.groupby(["spectrum_file","scan"])["scan"].transform("count")
    return ionbot

def get_unique_proteins(protein_series):
    uprotein = []
    for protein in protein_series:
        uprotein.extend(protein.split("||"))
    
    return set(uprotein)

def return_search_metrics(df):

    len_psms = len(df.index)
    len_peptides = len(set(df["matched_peptide"]+"|"+df["modifications"].fillna("")))
    len_proteins = len(get_unique_proteins(df["proteins"]))

    return [len_psms, len_peptides, len_proteins]

def get_all_mods(df):
    mod_list = []

    for ionbot_modifications in df["modifications"].fillna(""):
        mod_names = re.sub('\[.*?\]', '', ionbot_modifications)
        for mod_name in mod_names.split("|"):
            if mod_name == "":
                continue
            try:
                int(mod_name)
            except ValueError:
                mod_list.append(mod_name.lower())

    mod_list = [am for am in mod_list if am != ""]
    all_mods_counter_obj = Counter(mod_list)
    return all_mods_counter_obj