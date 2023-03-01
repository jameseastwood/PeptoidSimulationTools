import argparse
from ccdc.search import SMARTSSubstructure
from ccdc.search import SubstructureSearch
from ccdc.io import MoleculeReader
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--name', help='base name (residue name) to save files under')
parser.add_argument('-s', '--smiles', help='SMILES (or SMARTS) string representing the residue to be analyzed. ex. *N(CCOC)CC(=O)*. Atom order is important!')
parser.add_argument('-x', '--numchi', type=int, help='number of chi dihedrals to analyze')
parser.add_argument('-d', '--db_path', help='path to file defining database for CSD (e.g. gcd file of refcodes)')
args = parser.parse_args(['-n','NSpe','-s','*N([C@@H](C)c1ccccc1)CC(=O)*','-x','2','-d','peptoiddb.gcd',])
args = parser.parse_args()

peptoid_db = MoleculeReader(args.db_path)

def search_and_measure(residue_smiles, numchi, db=peptoid_db):
    #first we define a peptoid backbone substructure
    residue = SMARTSSubstructure(residue_smiles)
    peptoid_bb = SMARTSSubstructure('CC(=O)N(C*)CC(=O)NC')
    substructure_search = SubstructureSearch()
    bb_id = substructure_search.add_substructure(peptoid_bb)
    sub_id = substructure_search.add_substructure(residue)
    substructure_search.add_torsion_angle_measurement('chi1', sub_id, 0, sub_id, 1, sub_id, 2, sub_id, 3)
    substructure_search.add_torsion_angle_measurement('chi2', sub_id, 1, sub_id, 2, sub_id, 3, sub_id, 4)
    if (numchi > 2):
        substructure_search.add_torsion_angle_measurement('chi3', sub_id, 2, sub_id, 3, sub_id, 4, sub_id, 5)    
    substructure_search.add_torsion_angle_measurement('omega', bb_id, 0, bb_id, 1, bb_id, 3, bb_id, 6)
    substructure_search.add_torsion_angle_measurement('phi', bb_id, 1, bb_id, 3, bb_id, 6, bb_id, 7)
    substructure_search.add_torsion_angle_measurement('psi', bb_id, 3, bb_id, 6, bb_id, 7, bb_id, 9)
    hits = substructure_search.search(database=db)
    #print(len(hits))
    df = pd.DataFrame({h.identifier : h.measurements for h in hits})
    return(df.T)

residue_df = search_and_measure(args.smiles, args.numchi)

residue_df.hist(sharex=True, sharey=True, layout=(2,3), legend=False, bins=12)
plt.suptitle(args.name, fontweight='bold')
plt.savefig(args.name+'_hist.png', dpi=300)
plt.close()

pd.plotting.scatter_matrix(residue_df, figsize=(6,6), diagonal='kde', marker='x')
plt.suptitle(args.name, fontweight='bold')
plt.savefig(args.name+'_pairs.png', dpi=300)
plt.close()
