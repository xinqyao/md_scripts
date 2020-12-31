module load anaconda2

python mmgbsa_entropy.py --GBSA_dir=sample --GBSA_RL_FN=PBSA_RL.dat --GBSA_R_FN=PBSA_R.dat --GBSA_L_FN=PBSA_L.dat --ligand_dir=sample/ligand --ligand_prmtop_FN=top.pdb --ligand_dcd_FN=ligand.dcd --skip_analysis

