import os

min_result = [s for s in os.listdir('.') if os.path.isdir(s)]
labl = [x.split('result_min_')[1] for x in min_result]
for prot,folder in zip(labl,min_result):

        os.system('gmx grompp -f md.mdp -c %s/%s_npt.gro -t %s/%s_npt.cpt -p %s/topol.top -o %s/%s_md.tpr' %(folder,prot,folder,prot,folder,folder,prot))
        os.system('gmx mdrun -v -s %s/%s_md.tpr -o %s/%s_md -x %s/%s_md -deffnm %s/%s'%(folder,prot,folder,prot,folder,prot,folder,prot,))

		os.system("echo '1 0' | gmx trjconv -s %s/%s_md.tpr -f %s/%s_md.xtc -center -ur compact -pbc mol -o %s/%s_md_center.xtc" %(folder,prot,folder,prot,folder,prot))
		os.system("echo '1 0' | gmx trjconv -s %s/%s_md.tpr -f %s/%s_md_center.xtc -fit rot+trans -o %s/%s_md_center_rot_trans.xtc" %(folder,prot,folder,prot,folder,prot))

		