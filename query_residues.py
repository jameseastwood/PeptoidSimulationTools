import requests, json, sys
url = "https://databank.peptoids.org/api"
residues = requests.get(url+"/residues").json()
keycodes = [list(residue.keys())[0] for residue in residues]
values = [list(residue.values())[0] for residue in residues]
residue_dict = {k: v for k, v in zip(keycodes, values)}

name = residue_dict[sys.argv[1]]['short_name']
if name:
    print(name)
