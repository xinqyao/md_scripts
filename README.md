# Recipes for Molecular Dynamics Simulation #

This repository contains script and force field files for running and analyzing MD simulations 
in the Hamelberg Lab.

## Before MD ##
Normally, a MD simulation starts with a crystallographic structure from PDB. 
Prior to any processing, ask following questions:

* What is the biological unit/assembly you want to simulate?
   - This is to determine chains and cofactors to be included in the simulation.
   - Sometimes, a PDB file does not give a complete structure of a functional complex 
     due to symmetry. For example, `2DN1` contains only two chains but the system it 
     represents (i.e., hemoglobin) is a tetramer. Use `bio3d` to reconstruct the 
     biological unit, or download it directly from PDB.

* Are there missing residues?
   - If so, determine whether a structural modeling is needed.
   - If missing residues are only on terminii of a chain, you can ignore them.
     Optionally, you may cap the N- and C-terminus using acetyl and methylamide groups, respectively.
   - If missing residues break a chain, you may need to model them using Modeller or Swiss-Model.
     For short continuous missing residues (<5 aa), Modeller without additional template works well.
     For long missing fragments, structural templates from homologs may be needed.

* Any disulfide bonds?
   - If so, keep in mind that additional modifications are needed for both PDB and tleap input file.

* Any essential post-translational modifications (phosphorylation, glycosylation, etc.)?
   - If so, use proper force fields for the modified residues.
   - For phosphorylated residues, use the files in this repository.
   - For glycosylated residues, take a look of http://glycam.org.

* Any engineering mutations introduced only for expression or crystallization? 
   - If so, may consider changing them back to wildtype.

* Any molecules that require a force field out of Amber?
   - Check the website, http://research.bmh.manchester.ac.uk/bryce/amber, which collects force field files
     for a few commonly encountered cofactors, such as ATP, GTP, NADH, etc.
   - Google, check Amber tutorials, or ask colleagues to see if there is a force field available.
   - Use Antechamer to generate a force field (will need to load the GAFF force field in tleap).
   - Optionally, may run Gaussian to generate partial charges.

## Setup ##
First, split the PDB file into protein, ligand, and water (we keep crystal water in case these
water molecules have a functional role). Remove all hydrogen atoms if there is any 
(you can easily do it with `bio3d`).

Check pKa of residues and determine the protonation state of histidines. Recommend use PROPKA embeded in PDB2PQR.
You can use either standalone or web server of PDB2PQR. For standalone, following command should be used:

```
#!bash
pdb2pqr --ff=amber --ffout=amber --chain --nodebump --noopt --ph-calc-method=propka --with-ph=7 infile outfile
```

Check the out.propka file. In most cases, Asp and Glu will be deprotonated (negatively charged), Arg and Lys will 
be protonated (positively charged), and His will be deprotonated (HID or HIE). If other protonation states are 
assigned, pay attention to them. Then, determine HID/HIE for His residues by visually inspecting 
micro-environment of the residue using e.g. PyMol. If the hydrogen is located close to other hydrogen atoms 
associated with carbon or nitrogen, it may need to be moved to the alternative nitrogen of the histidine 
(from HID to HIE or reverse). On the other hand, if there is potential to form hydrogen bonds (e.g. close to 
an oxygen), the proton should be placed accordingly. Once determined the protonation state, change
residue names from HIS to HID/HIE (or other ionizable residues if necessary) in the original PDB file of protein 
(or the file generated by PDB2PQR).

If your protein contains disulfide bonds, change the name of corresponding cysteins from CYS to CYX.

If you plan to run both wildtype and mutants, do mutagenesis now. I recommend use the **Mutagenesis Wizard**
of PyMol because by it you can select different side-chain conformations.

Make a working directory for a specific system and copy all related force field files, structure files (protein,
ligand, and water), and the three '.in' files under it. Modify 'new_prep.in' to match your case (i.e. file names
and options in tleap.in) and then run `./new_prep.in`.

Check the 'prep/log' file to make sure there is no errors and significant warnings. Some warnings like 
'WARNING: The unperturbed charge of the unit: 1.000000 is not zero.' can be ignored.

Also recommend check 'prep/XXX_box.prmtop' and 'prep/XXX_box.inpcrd' files using VMD. Make sure everything 
looks fine, especially there is no weird geometry in cofactors.

## Energy minimization, heating, equilibration ##
Modify 'new_equil.in' and 'new_production.in' files to match your case. In most cases, don't change parameters 
of the input files of `pmemd`.  Then, run

```
#!bash
./new_equil.in
./new_production.in
```

Go to 'equil/' and run

```
#!bash
nohup ./jobe.in >& err.log &
```

## Production ##
After equilibration is done, go to 'production/' and run
```
#!bash
nohup ./jobp.in >& err.log &
```

## Analyzing trajectories ##
 
