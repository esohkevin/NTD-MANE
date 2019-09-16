# Directives on saving TLR4 chain A and B
- `pymol 4A8G.pdb`
- `indicate c. A:B
`
- Click File; Save molecule; Ok
- Remove HOH molecules 
	`grep -v "HOH" file.pdb > newFile.pdb`
	
	Not always adviceable to grep cos one may be interested 
	in crystal water around active sites
	
- Generate topology files with pdb2gmx: You'll be prompted to choose a force field
	This is an important step. Make sure to read literature to understand what
	force field is suitable for the study
	
	`gmx pdb2gmx -f 1AKI.pdb -o 1AKI_processed.gro -water spce`
	
- Define the simulation box (cubic, rhombic dodecahedron - preferably as it allows 
	enough space for solvent molecules)
	
	`gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic`
	
- Fill the box with solvent (water)

	`gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top`
	
	spc216.gro is part gromacs standard installation. One can choose any
- Next we have to add ions with gmx genion. But first we have to create the input
	file for genion using gmxpp, a gromacs preprocesor, and obtain a molecular
	dynamic parameter file (.mdp). This file is normally used to run energy 
	minimizations in MD. One can search for templates for specific projects
	
	`gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr`
	`gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral`
	
- Now we perform energy minimization on the solvated electroneutral system in order to
	'relax' the structure. We'll use grompp again to make a .tpr file and pass it to 
	mdrun
	
	`gmx grompp -f mdout.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr`
	`gmx mdrun -v -deffnm em`
	
- Now we can perform some first pass analysis by looking at the effect of EM on our system:

	`gmx energy -f em.edr -o potential.xvg`
	
	The Etot (Total energy) should be negative indicating that the system is stable and
	ready fo simulation. If positive. it means the system is unstable and one should
	consider going through the process once more
- Next we equilibrate the slvent and ions around the protein to prevent the system from
	collapsing during unrestrained dynamics. Bring the solvent to temp at which simula
	would be done. Equilibrate first at isothermal-isochoric (NVT cons # particles,
	volume and temp) - make take some time to run!
	
	`gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr`
	`gmx mdrun -deffnm nvt`
