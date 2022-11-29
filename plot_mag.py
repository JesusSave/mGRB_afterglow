## Comparision of U3s parameters: Python jump model, Z&K05 Fig 1 results, Z&K05 Modification with adiabatic index  
# version 30,10,2019
# Hallelujah!

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
cm = plt.cm.get_cmap('jet')
import matplotlib.ticker
from scipy.optimize import curve_fit

PATH = './'
DFILE_lc = 'data/lightcurves'
DFILE_dyn = 'data/dynamics'
DFILE_mag = 'data/magnetization'
DFIG = 'figs'

plot_model = False
plot_model = True

plot_jump = True
plot_jump = False

simnum = 4
# (PATH, folder, filename), sigma, (ls, lw, color)
SIMs = [
  [(PATH, DFILE_mag, 'mag'),0.1, ('-',1,1./simnum)]
]
sigma_cmp = 0


fig1, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex=True,figsize=(9,6))

# Two subplots, the axes array is 1-d
mf = matplotlib.ticker.ScalarFormatter(useMathText=True)
mf.set_powerlimits((-2,2))
fig1.gca().yaxis.set_major_formatter(mf)
fig1.tight_layout(rect=[0.02, 0.025, 1, 0.99])

ax1.set_title('(a)',loc='left')
ax1.set_ylabel(r'log($\rm u_{3s}$)')
ax1.set_xlabel(r'log($\rm \gamma_{34}$)')
ax1.axis([1e0,1e3,1e-1,1e2])
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([1e0,1e3])

ax2.set_title('(b)',loc='left')
ax2.set_ylabel(r'log($\rm \gamma_{4s}/\gamma_{34}$)')
ax2.set_xlabel(r'log($\rm \gamma_{34}$)')
ax2.axis([1e0,1e3,1e0,1e2])
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim([1e0,1e3])

ax3.set_title('(c)',loc='left')
ax3.set_ylabel(r'$\rm f_a=(e_3/n_3m_pc^2)/(\gamma_{34}-1)$')
ax3.set_xlabel(r'log($\rm \gamma_{34}$)')
ax3.axis([1e0,1e3,1e-1,2e0])
ax3.set_xscale('log')
ax3.set_xlim([1e0,1e3])

ax4.set_title('(d)',loc='left')
ax4.set_ylabel(r'$\rm f_b=(n_3/n_4)/(4\gamma_{34}+3)$')
ax4.set_xlabel(r'log($\rm \gamma_{34}$)')
ax4.axis([1e0,1e3,1e-1,2e0])
ax4.set_xscale('log')
ax4.set_xlim([1e0,1e3])

ax5.set_title('(e)',loc='left')
ax5.set_ylabel(r'log($\rm p_{b,3}/p_3$)=log($f_c-1$)')
ax5.set_xlabel(r'log($\rm \gamma_{34}$)')
ax5.axis([1e0,1e3,1e-1,1e5])
ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.set_xlim([1e0,1e3])

ax6.set_title('(f)',loc='left')
ax6.set_ylabel(r'log(F)')
ax6.set_xlabel(r'log($\rm \gamma_{34}$)')
ax6.axis([1e0,1e3,5e-1,1e4])
ax6.set_xscale('log')
ax6.set_yscale('log')
ax6.set_xlim([1e0,1e3])


linestyles = ['-','-','-','-','-','-','-','-','-','-']

if (plot_jump):
  for tt in N:
    print ('reading timestep data t =',tt,N.index(tt) ,'/', len(N))
    ls = linestyles[N.index(tt)]
    cnorm = (N.index(tt)+1.0)/len(N)

    g = open(DIR3+'jump_u'+str(tt)+'.dat','r') 
    data = np.array([[float(data) for data in line.split()] for line in g.readlines()])
    g.close()

    gam34 = data[:,0]
    u3s = data[:,1]
    fa = data[:,4]
    fb = data[:,5]
    fc = data[:,6]
    F = fa*fb*fc

    pb3p3 = [x-1.0 for x in fc]

    gam3s = np.sqrt(u3s*u3s+1.0)
    beta34 = np.sqrt(gam34*gam34-1.0)/gam34
    g4sg34 = gam3s+u3s*beta34
  
    if (N.index(tt)==0):
      ax1.plot(gam34,u3s,ls =ls,lw=3*lw,c=cm(cnorm),label='jump (<tx)')
      ax2.plot(gam34,g4sg34,ls =ls,lw=3*lw,c=cm(cnorm))
      ax3.plot(gam34,fa,ls =ls,lw=3*lw,c=cm(cnorm))
      ax4.plot(gam34,fb,ls =ls,lw=3*lw,c=cm(cnorm))
      ax6.plot(gam34,F,ls =ls,lw=3*lw,c=cm(cnorm))
    else:
      ax1.plot(gam34,u3s,ls =ls,c=cm(cnorm))
      ax2.plot(gam34,g4sg34,ls =ls,c=cm(cnorm))
      ax3.plot(gam34,fa,ls =ls,c=cm(cnorm))
      ax4.plot(gam34,fb,ls =ls,c=cm(cnorm))
      ax5.plot(gam34,pb3p3,ls =ls,c=cm(cnorm))
      ax6.plot(gam34,F,ls =ls,c=cm(cnorm))


if (plot_model):
  for k in range(len(SIMs)):
    print ('sigma =',k+1,'/',len(SIMs))

    (PATH, DFILE, filename),sigma, (ls,lw,cnorm) = SIMs[k]
    print ('  MODEL SIM',k+1,'/',len(SIMs))

    print('    reading fluxfs data')
    f = open(PATH+f'/{DFILE}/{filename}.dat','r')
    data = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()

    t = data[:,0]
    gam3 = data[:,1]
    gam34 = data[:,2]
    u3s = data[:,3]
    fa = data[:,4]
    fb = data[:,5]
    fc = data[:,6]

    gam3s = np.sqrt(u3s*u3s+1.0)
    beta34 = np.sqrt(gam34*gam34-1.0)/gam34
    g4sg34 = gam3s+u3s*beta34
    F = fa*fb*fc

    if (sigma==0):
      ax1.plot(gam34,u3s,ls =ls,lw=3*lw,c=cm(cnorm),label='Torch')
      ax2.plot(gam34,g4sg34,ls =ls,lw=3*lw,c=cm(cnorm))
      ax3.plot(gam34,fa,ls =ls,lw=3*lw,c=cm(cnorm))
      ax4.plot(gam34,fb,ls =ls,lw=3*lw,c=cm(cnorm))
      ax6.plot(gam34,F,ls =ls,lw=3*lw,c=cm(cnorm))
    else:
      ax1.plot(gam34,u3s,ls =ls,lw=lw,c=cm(cnorm))
      ax2.plot(gam34,g4sg34,ls =ls,lw=lw,c=cm(cnorm))
      ax3.plot(gam34,fa,ls =ls,lw=lw,c=cm(cnorm))
      ax4.plot(gam34,fb,ls =ls,lw=lw,c=cm(cnorm))
      ax5.plot(gam34,fc-1.,ls =ls,lw=lw,c=cm(cnorm))
      ax6.plot(gam34,F,ls =ls,lw=lw,c=cm(cnorm))


ax1.tick_params(which='both', direction='in', right=True, top = True)
ax2.tick_params(which='both', direction='in', right=True, top = True)
ax3.tick_params(which='both', direction='in', right=True, top = True)
ax4.tick_params(which='both', direction='in', right=True, top = True)
ax5.tick_params(which='both', direction='in', right=True, top = True)
ax6.tick_params(which='both', direction='in', right=True, top = True)

ax1.legend(loc=0,frameon=False,ncol=3)
plt.tight_layout()

fig1.savefig(PATH+f"/{DFIG}/mag_parameters.png",format='png')
fig1.savefig(PATH+f"/{DFIG}/mag_parameters.pdf",format='pdf')
plt.show()


