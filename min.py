# -*- coding: utf-8 -*-
import os
import sys


protein_list = [s for s in os.listdir('.') if s.endswith('.pdb')]

for protein in protein_list:
    protein=protein.split('.pdb')[0]
    print 'creating directories'
    os.mkdir('result_min_%s'%protein)
    os.system('gmx pdb2gmx -f %s.pdb -o %s_processed.gro -water spce -ff oplsaa' % (protein,protein))
    os.system('mv posre.itp result_min_%s' %protein)
    os.system('mv topol.top result_min_%s' %protein)
    os.system('mv *%s_proc* result_min_%s'%(protein,protein))


    os.system('gmx editconf -f result_min_%s/%s_processed.gro -o result_min_%s/%s_newbox.gro -c -d 2.0 -bt cubic'%(protein,protein,protein,protein))


    #------------------------------------------
    #Solvatation
    #------------------------------------------
    os.system('gmx solvate -cp result_min_%s/%s_newbox.gro -cs spc216.gro -o result_min_%s/%s_solv.gro -p result_min_%s/topol.top' %(protein,protein,protein,protein,protein))


    #------------------------------------------
    #Definig Ions
    #------------------------------------------


    os.system('gmx grompp -f ions.mdp -c result_min_%s/%s_solv.gro -p result_min_%s/topol.top -o result_min_%s/%s_ions.tpr'%(protein,protein,protein,protein,protein))



    os.system('echo 13 | gmx genion -s result_min_%s/%s_ions.tpr -o result_min_%s/%s_solv_ions.gro -p result_min_%s/topol.top -pname NA -nname CL -neutral -conc 0.15' %(protein,protein,protein,protein,protein))


    #------------------------------------------
    #Energy minimization
    #------------------------------------------


    print "+++ ENERGY MINIZATION +++"

    os.system('gmx grompp -f minim.mdp -c result_min_%s/%s_solv_ions.gro -p result_min_%s/topol.top -o result_min_%s/%s_em.tpr' %(protein,protein,protein,protein,protein))
    os.system('gmx mdrun -v -deffnm result_min_%s/%s_em'%(protein,protein))

    os.system("echo '10 0' | gmx energy -f result_min_%s/%s_em.edr -o result_min_%s/%s_potential.xvg" %(protein,protein,protein,protein))

    #------------------------------------------
    #NVT Equilibration
    #------------------------------------------

    os.system('gmx grompp -f nvt.mdp -c result_min_%s/%s_em.gro -p result_min_%s/topol.top -o result_min_%s/%s_nvt.tpr'%(protein,protein,protein,protein,protein))
    os.system('gmx mdrun -deffnm result_min_%s/%s_nvt -v'%(protein,protein))

    os.system("echo '15 0' | gmx energy -f result_min_%s/%s_nvt.edr -o result_min_%s/%s_temperature.xvg" %(protein,protein,protein,protein))

    #------------------------------------------
    #System and pressure equilibration
    #------------------------------------------



    os.system('gmx grompp -f npt.mdp -c result_min_%s/%s_nvt.gro -t result_min_%s/%s_nvt.cpt -p result_min_%s/topol.top -o result_min_%s/%s_npt.tpr'%(protein,protein,protein,protein,protein,protein,protein))
    os.system('gmx mdrun -deffnm result_min_%s/%s_npt'%(protein,protein))

    os.system("echo '16 0' | gmx energy -f result_min_%s/%s_nvt.edr -o result_min_%s/%s_pressure.xvg" %(protein,protein,protein,protein))
    os.system("echo '22 0' | gmx energy -f result_min_%s/%s_nvt.edr -o result_min_%s/%s_density.xvg" %(protein,protein,protein,protein))
