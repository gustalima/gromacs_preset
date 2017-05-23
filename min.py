# -*- coding: utf-8 -*-
import os
import sys


protein_list = [s for s in os.listdir('.') if s.endswith('.pdb')]

for protein in protein_list:
    protein = protein.split('.pdb')[0]
    os.mkdir(protein+'minim')

    print "++ SYSTEM SETUP ++"

    #------------------------------------------
    #Water model
    #------------------------------------------
    water_model = 'spce'
    os.system('gmx pdb2gmx -f %s.pdb -o %s_processed.gro -water %s' %(protein, protein, water_model))
    #------------------------------------------
    #------------------------------------------

    #------------------------------------------
    #Box type
    #------------------------------------------
    box_type = 'cubic'

    #------------------------------------------
    #Solute distance
    #------------------------------------------
    dist_sol = 2.0

    os.system('gmx editconf -f %s_processed.gro -o %s_newbox.gro -c -d %s -bt %s'%(protein,protein,dist_sol,box_type))



    #------------------------------------------
    #Solvatation
    #------------------------------------------
    cs = 'spc216'
    os.system('gmx solvate -cp %s_newbox.gro -cs %s.gro -o %s_solv.gro -p %s_topol.top' %(protein,cs,protein,protein))


    #------------------------------------------
    #Definig Ions
    #------------------------------------------

    try:
        os.remove('ions.mdp')
    except OSError:
        pass

    with open('ions.mdp', 'w') as ion:
        ion.write('''
    ; ions.mdp - used as input into grompp to generate ions.tpr
    ; Parameters describing what to do, when to stop and what to save
    integrator	= steep		; Algorithm (steep = steepest descent minimization)
    emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
    emstep      = 0.01      ; Energy step size
    nsteps		= 50000	  	; Maximum number of (minimization) steps to perform

    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
    nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
    cutoff-scheme   = Verlet
    ns_type		    = grid		; Method to determine neighbor list (simple, grid)
    coulombtype	    = PME		; Treatment of long range electrostatic interactions
    rcoulomb	    = 1.0		; Short-range electrostatic cut-off
    rvdw		    = 1.0		; Short-range Van der Waals cut-off
    pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)
    ''')


    os.system('gmx grompp -f ions.mdp -c %s_solv.gro -p %s_topol.top -o %s_ions.tpr'%(protein,protein,protein))

    genion_options = '-pname NA -nname CL -neutral -conc 0.15'

    os.system('gmx genion -s %s_ions.tpr -o %s_solv_ions.gro -p %s_topol.top -pname NA -nname CL -neutral -conc 0.15' %(protein,protein,protein))


    #------------------------------------------
    #Energy minimization
    #------------------------------------------

    rootdir = os.curdir

    try:
        os.remove('minim.mdp')
    except OSError:
        pass

    with open('minim.mdp', 'w') as minim:
        minim.write('''
    ; minim.mdp - used as input into grompp to generate em.tpr
    integrator	= steep		; Algorithm (steep = steepest descent minimization)
    emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
    emstep      = 0.01      ; Energy step size
    nsteps		= 50000	  	; Maximum number of (minimization) steps to perform

    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
    nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
    cutoff-scheme   = Verlet
    ns_type		    = grid		; Method to determine neighbor list (simple, grid)
    coulombtype	    = PME		; Treatment of long range electrostatic interactions
    rcoulomb	    = 1.0		; Short-range electrostatic cut-off
    rvdw		    = 1.0		; Short-range Van der Waals cut-off
    pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)
    ''')

    print "+++ ENERGY MINIZATION +++"

    os.system('gmx grompp -f minim.mdp -c %s_solv_ions.gro -p %s_topol.top -o %s_em.tpr' %(protein,protein,protein))
    os.system('gmx mdrun -v -deffnm %s_em'%protein)



    #------------------------------------------
    #NVT Equilibration
    #------------------------------------------

    try:
        os.remove('nvt.mdp')
    except OSError:
        pass

    with open('nvt.mdp', 'w') as nvt:
        nvt.write('''
    title		= NVT equilibration
    define		= -DPOSRES	; position restrain the protein
    ; Run parameters
    integrator	= md		; leap-frog integrator
    nsteps		= 50000		; 2 * 50000 = 100 ps
    dt		    = 0.002		; 2 fs
    ; Output control
    nstxout		= 500		; save coordinates every 1.0 ps
    nstvout		= 500		; save velocities every 1.0 ps
    nstenergy	= 500		; save energies every 1.0 ps
    nstlog		= 500		; update log file every 1.0 ps
    ; Bond parameters
    continuation	        = no		; first dynamics run
    constraint_algorithm    = lincs	    ; holonomic constraints
    constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
    lincs_iter	            = 1		    ; accuracy of LINCS
    lincs_order	            = 4		    ; also related to accuracy
    ; Neighborsearching
    cutoff-scheme   = Verlet
    ns_type		    = grid		; search neighboring grid cells
    nstlist		    = 10		; 20 fs, largely irrelevant with Verlet
    rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
    rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
    ; Electrostatics
    coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics
    pme_order	    = 4		; cubic interpolation
    fourierspacing	= 0.16	; grid spacing for FFT
    ; Temperature coupling is on
    tcoupl		= V-rescale	            ; modified Berendsen thermostat
    tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
    tau_t		= 0.1	  0.1           ; time constant, in ps
    ref_t		= 300 	  300           ; reference temperature, one for each group, in K
    ; Pressure coupling is off
    pcoupl		= no 		; no pressure coupling in NVT
    ; Periodic boundary conditions
    pbc		= xyz		    ; 3-D PBC
    ; Dispersion correction
    DispCorr	= EnerPres	; account for cut-off vdW scheme
    ; Velocity generation
    gen_vel		= yes		; assign velocities from Maxwell distribution
    gen_temp	= 300		; temperature for Maxwell distribution
    gen_seed	= -1		; generate a random seed
    ''')

    os.system('gmx grompp -f nvt.mdp -c %s_em.gro -p %s_topol.top -o %s_nvt.tpr'%(protein,protein,protein))
    os.system('gmx mdrun -deffnm %s_nvt -v'%protein)


    #------------------------------------------
    #System and pressure equilibration
    #------------------------------------------
    try:
        os.remove('npt.mdp')
    except OSError:
        pass

    with open('npt.mdp', 'w') as npt:
        npt.write('''
    title		= NPT equilibration
    define		= -DPOSRES	; position restrain the protein
    ; Run parameters
    integrator	= md		; leap-frog integrator
    nsteps		= 50000		; 2 * 50000 = 100 ps
    dt		    = 0.002		; 2 fs
    ; Output control
    nstxout		= 500		; save coordinates every 1.0 ps
    nstvout		= 500		; save velocities every 1.0 ps
    nstenergy	= 500		; save energies every 1.0 ps
    nstlog		= 500		; update log file every 1.0 ps
    ; Bond parameters
    continuation	        = yes		; Restarting after NVT
    constraint_algorithm    = lincs	    ; holonomic constraints
    constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
    lincs_iter	            = 1		    ; accuracy of LINCS
    lincs_order	            = 4		    ; also related to accuracy
    ; Neighborsearching
    cutoff-scheme   = Verlet
    ns_type		    = grid		; search neighboring grid cells
    nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
    rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
    rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
    ; Electrostatics
    coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
    pme_order	    = 4		    ; cubic interpolation
    fourierspacing	= 0.16		; grid spacing for FFT
    ; Temperature coupling is on
    tcoupl		= V-rescale	            ; modified Berendsen thermostat
    tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
    tau_t		= 0.1	  0.1	        ; time constant, in ps
    ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
    ; Pressure coupling is on
    pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
    pcoupltype	        = isotropic	            ; uniform scaling of box vectors
    tau_p		        = 2.0		            ; time constant, in ps
    ref_p		        = 1.0		            ; reference pressure, in bar
    compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
    refcoord_scaling    = com
    ; Periodic boundary conditions
    pbc		= xyz		; 3-D PBC
    ; Dispersion correction
    DispCorr	= EnerPres	; account for cut-off vdW scheme
    ; Velocity generation
    gen_vel		= no		; Velocity generation is off
    ''')


    os.system('gmx grompp -f npt.mdp -c %s_nvt.gro -t %s_nvt.cpt -p %s_topol.top -o %s_npt.tpr'%(protein,protein,protein,protein))
    os.system('gmx mdrun -deffnm %s_npt'%protein)

    #---------------------
    #MOVE RESULTS
    #---------------------
    out_list = filter(os.path.isfile, [s for s in os.listdir('.') if not s.endswith('.pdb')])
    for f_ile in out_list:
        os.system('mv %s %s/'%(f_ile,protein))
