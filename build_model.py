import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--smiles', help='SMILES or SMARTS string representing the residue to be modelled') 
parser.add_argument('-n', '--nterm', help='capping group SMILES for the N-terminus, with resname if applicable', nargs=2)
parser.add_argument('-c', '--cterm', help='capping group SMILES for the C-terminus, with resname if applicable', nargs=2)
parser.add_argument('-f', '--name', help='file name to save under')
args = parser.parse_args(['-s', '*N([C@@H](C)COC)CC(=O)*', '-n', 'CC(=O)*', 'ACE', '-c', '*N(C)C', 'NDM', '-f', 'NsMP'])
args = parser.parse_args()

#Input: SMILES string of residue, SMILES strings of capping groups
#output: mol2 file of capped residue in arbitrary geometry

def substitute_polymer(n, c):
    from rdkit.Chem import AllChem, Draw
    atom_map = {}
    #declare that we have not made the connection yet
    disconnected = True
    #First we find the last connection point on the nterm residue
    i = 1
    while n.GetAtomWithIdx(n.GetNumAtoms()-i).GetSymbol() != '*':
        i += 1
    c_connector = n.GetNumAtoms()-i
    print(c_connector)
    #The connection can only have one neighbor, or it's not a valid connection.
    c_connect = n.GetAtomWithIdx(c_connector).GetNeighbors()[0]
    print(c_connect.GetSymbol(), c_connect.GetIdx())
    #make a writeable copy of the N side, add every atom in the c side to it, map the id in c to the new id. The first connection point is removed, and all bonds are added to the last connection point of the Nterm fragment.
    pw = Chem.RWMol(n) # can't use `with` block because RWMol doesn';t have __enter__ defined
    for atom in c.GetAtoms():
        if disconnected:
            if atom.GetSymbol() == '*':
                atom_map[atom.GetIdx()] = c_connect.GetIdx()
                disconnected = False
        else:
            nid = pw.AddAtom(atom)
            atom_map[atom.GetIdx()] = nid
        #add a bond to the polymer for each bond in the cterm res, calling on the map of atom ids.
    for bond in c.GetBonds():
        print(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
        print(atom_map[bond.GetBeginAtomIdx()], atom_map[bond.GetEndAtomIdx()])
        pw.AddBond(atom_map[bond.GetBeginAtomIdx()], atom_map[bond.GetEndAtomIdx()], bond.GetBondType())
    pw.RemoveAtom(c_connector)
    #print(Chem.MolToSmiles(pw))
    #AllChem.Compute2DCoords(pw)
    #Draw.ShowMol(pw)
    Chem.SanitizeMol(pw)
    return(pw)

def label_sc_atoms(root_atom, label_list, counter=1):
    for a in root_atom.GetNeighbors():
        if a.GetMonomerInfo().GetName():
            pass
        else:
            if a.GetAtomicNum() == 6:
                try:
                    a.GetMonomerInfo().SetName(' C'+label_list[0]+' ')
                except:
                    a.GetMonomerInfo().SetName('{0:^4s}'.format(a.GetSymbol()+str(a.GetIdx())))
                new_list = label_list[1:]
            else:
                a.GetMonomerInfo().SetName(' ' + a.GetSymbol()+root_atom.GetMonomerInfo().GetName()[2:])
                new_list = label_list
            label_sc_atoms(a, new_list)

ni = Chem.AtomPDBResidueInfo()
ni.SetResidueName(args.nterm[1])
ni.SetResidueNumber(1)
ri = Chem.AtomPDBResidueInfo()
#Do some reformatting to make sure this is 3 uppercase characters
ri.SetResidueName('{0:<3}'.format(args.name[args.name.find('N')+1:args.name.find('N')+4].upper()))
ri.SetResidueNumber(2)
ci = Chem.AtomPDBResidueInfo()
ci.SetResidueName(args.cterm[1])
ci.SetResidueNumber(3)

nmol = Chem.MolFromSmiles(args.nterm[0])
[atom.SetMonomerInfo(ni) for atom in nmol.GetAtoms()]
[atom.GetPDBResidueInfo().SetName('{0: <4}'.format(atom.GetSymbol()+str(atom.GetIdx()))) for atom in nmol.GetAtoms()]
rmol = Chem.MolFromSmiles(args.smiles)
[atom.SetMonomerInfo(ri) for atom in rmol.GetAtoms()]
cmol = Chem.MolFromSmiles(args.cterm[0])
[atom.SetMonomerInfo(ci) for atom in cmol.GetAtoms()]
[atom.GetPDBResidueInfo().SetName('{0: <4}'.format(atom.GetSymbol()+str(atom.GetIdx()))) for atom in cmol.GetAtoms()]

m = substitute_polymer(nmol, substitute_polymer(rmol, cmol))
AllChem.Compute2DCoords(m)
img = Chem.Draw.MolToFile(m, args.name+'.png')

peptoid_bb = Chem.MolFromSmarts('[C](=[O])[N:1]([#6:5])[C:2][C:3](=[O:4])[N](C)C')
ind_map_bb = {}
for atom in peptoid_bb.GetAtoms() :
    map_num = atom.GetAtomMapNum()
    if map_num:
        ind_map_bb[map_num-1] = atom.GetIdx()
map_list_bb = [ind_map_bb[x] for x in sorted(ind_map_bb)]

bb_labels = [' N  ', ' CA ', ' C  ', ' O  ', ' CAN']
sc_labels = ['B', 'G', 'E', 'Z']
for match in m.GetSubstructMatches( peptoid_bb ) :
    mas = [match[x] for x in map_list_bb]

for a in range(0, len(mas)):
    m.GetAtomWithIdx(mas[a]).GetPDBResidueInfo().SetName(bb_labels[a])
label_sc_atoms(m.GetAtomWithIdx(mas[-1]), sc_labels)

m1 = Chem.AddHs(m, addResidueInfo=True)
#need to fix H names so prepgen doesn't get confused. They must be unique
i = 1
j = 1
for atom in m1.GetAtoms():
    if (atom.GetPDBResidueInfo().GetResidueName() == args.cterm[1]) & (atom.GetSymbol() == 'H'):
        atom.GetPDBResidueInfo().SetName('{0:^4s}'.format('HC'+str(i)))
        i += 1
    elif (atom.GetPDBResidueInfo().GetResidueName() == args.nterm[1]) & (atom.GetSymbol() == 'H'):
        atom.GetPDBResidueInfo().SetName('{0:^4}'.format('HN'+str(j)))
        j += 1
#AllChem.EmbedMolecule(m1, randomSeed=0xf00d)
Chem.rdDistGeom.EmbedMolecule(m1, enforceChirality=True)
AllChem.MMFFOptimizeMolecule(m1)
Chem.MolToPDBFile(m1, flavor=8,filename=args.name+'.pdb')
#Chem.MolToMolFile(m1, filename=args.name+'.mol2', includeStereo=True)
