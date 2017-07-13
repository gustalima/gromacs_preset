import matplotlib.pyplot as plt


with open('rmsd.xvg', 'r') as rmsd_inp:
        raw_data = rmsd_inp.readlines() 

data = []

for i in raw_data: 
        if not any(x in i for x in list("@#&") ):
                data.append(i)
        if "@    title" in i:
                plot_title = i.split("@    title")[1].replace('"', '')
        if "@    xaxis  label" in i:
                xlbl = i.split("@    xaxis  label")[1].replace('"', '')
        if "@    yaxis  label" in i:
                ylbl = i.split("@    yaxis  label")[1].replace('"', '')
        if "#   gmx rms -s" in i:
                legend = i.split('#   gmx rms -s')[1].split('_md.tpr')[0]

data = [[x.split()[0],x.split()[1]] for x in data]
        
plt.plot(*zip(*data), label='RMSD %s' %legend)
plt.ylabel(ylbl)
plt.xlabel(xlbl)
plt.title(plot_title)
plt.legend()
plt.grid(True)

plt.savefig('rmsd_%s.png' %legend)

plt.show()
