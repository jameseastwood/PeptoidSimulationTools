for file in smiles/*.smi
do name=`basename $file '.smi'`
   shortname=`python query_residues.py $name | awk '{print toupper($0)}'`
   smiles=`cat $file`
   python build_model.py -s $smiles -n 'CC(=O)*' ACE -c '*N(C)C' NDM -f $shortname
done
