#!/usr/bin/python3
# -*- coding: utf-8 -*
# Ref: 1507.01984v1, Gao & Meszaros 2015

import numpy as np
import matplotlib.pyplot as plt
cm = plt.cm.get_cmap('jet')

PATH = '.'

plot_ana = True
plot_ana = False

plot_model = False
plot_model = True

plot_jump = True
plot_jump = False

row = 2
col = 4

plt.clf()
plt.minorticks_on()
plt.gcf().set_size_inches(col*4,row*4)

plt1 = plt.subplot(row,col,1)
plt2 = plt.subplot(row,col,2)
plt3 = plt.subplot(row,col,3)
plt4 = plt.subplot(row,col,4)
plt5 = plt.subplot(row,col,5)
plt6 = plt.subplot(row,col,6)
plt7 = plt.subplot(row,col,7)
plt8 = plt.subplot(row,col,8)


def gam2(t):
  if k == 0: # ISM
    if Delt0 > SedL/(gam4**(8./3)): # Gao He 15, xi < 1, thick shell
      #return np.where(t<ttw,(SedL/Delt0)**(3./8)*(4*t/ttw)**(-1./4),(17.0*E0/(1024.0*pi*n*mp*c**5*t**3))**(1./8))
      return np.where(t<ttw,(SedL/Delt0)**(3./8)*(4*t/ttw)**(-1./4),(SedL/Delt0)**(3./8)*(4*ttw/ttw)**(-1./4)*(t/ttw)**(-3./8))
    else:
      return np.where(t<ttn,gam4,(17.0*E0/(1024.0*pi*n*mp*c**5*t**3))**(1./8))
  else:
    if Delt0 > SedL/(gam4**(8./3)): # Wind thick shell
      return np.where(t<ttw, 1.0/np.sqrt(2)*(SedL/Delt0)**(1./4),(17.0*E0/(1024.0*pi*n*mp*c**5*t**3))**(1./8))
    else:
      return np.where(t<ttn,gam4,(17.0*E0/(1024.0*pi*n*mp*c**5*t**3))**(1./8))

def gam3(t):
  if k == 0: # ISM
    if Delt0 > SedL/(gam4**(8./3)): # Gao He 15, xi < 1, thick shell
      return np.where(t<ttw,(SedL/Delt0)**(3./8)*(4*t/ttw)**(-1./4),(SedL/Delt0)**(3./8)*(4*ttw/ttw)**(-1./4)*(t/ttw)**(-7./16))
    else:
      return np.where(t<ttn,gam4,gam4*(t/ttn)**(-2./5))
  else:
    if Delt0 > SedL/(gam4**(8./3)): # Wind thick shell
      return np.where(t<ttw, 1.0/np.sqrt(2)*(SedL/Delt0)**(1./4),1.0/np.sqrt(2)*(SedL/Delt0)**(1./4)*(t/ttw)**(-3./8))
    else:
      return np.where(t<ttn,gam4,gam4*(t/ttn)**(-1./3))

def n3(t):
  if k ==0: #ISM
    if Delt0 > SedL/(gam4**(8./3)): # thick
      return np.where(t<ttw,8.*gam3(t)**3*n/gam4,8*gam3(ttw)**3*n/gam4*(t/ttw)**(-13./16))
  else:
    if Delt0 > SedL/(gam4**(8./3)):
      return np.where(t<ttn,7*n*gam4**2*(t/ttn)**-3,7*n*gam4**2*(t/ttn)**(-6./7))
    else:
      return np.where(t<ttn,7.0*A*gam4**6/l**2*(t/ttn)**-3, 7.0*A*gam4**6/l**2*(t/ttn)**(-7.0/8))

def e3(t):
  if k ==0: #ISM
    if Delt0 > SedL/(gam4**(8./3)): # thick
      return np.where(t<ttw,4.*gam3(t)**2*n*mp*c**2,4*gam3(ttw)**2*n*mp*c**2*(t/ttw)**(-13./12))
      #return np.where(t<ttw,(gam34(t)-1.)*n3(t)*mp*c**2,4*gam3(ttw)**2*n*mp*c**2*(t/ttw)**(-13./12))
    else:
      return np.where(t<ttn,4*gam4**2*n*mp*c**2,4*gam4**2*n*mp*c**2*(t/ttn)**(-8./7))
  else:
    if Delt0 > SedL/(gam4**(8./3)): # Wind thick
      return np.where(t<ttw, 8.0*A*mp*c*c/SedL**(1./2)/Delt0**(3./2)*(t/ttw)**(-2), 8.0*A*mp*c*c/SedL**(1./2)/Delt0**(3./2)*(t/ttw)**(-3./2))
    else:
      return np.where(t<ttn,4.0*A*mp*c*c*gam4**6.0/SedL**2.0*(t/ttn)**(-2.),4.0*A*mp*c*c*gam4**6.0/SedL**2.0*(t/ttn)**(-32./21))

def Ne3(t):
  if k ==0: #ISM
    if Delt0 > SedL/(gam4**(8./3)): # thick
      return np.where(t<ttw,N0*(t/ttw),N0)
    else:
      return np.where(t<ttn,N0*(t/ttn)**(3./2),N0)
  else:
    if Delt0 > SedL/(gam4**(8./3)):
      return np.where(t<ttw,N0*(t/ttw),N0)
    else:
      return np.where(t<ttn,N0*(t/ttn)**(1./2),N0)

gem = lambda t : epseR*e3(t)/(n3(t)*me*c**2)*(p-2)/(p-1)
BR  = lambda t : np.sqrt(8*np.pi*epsBR*e3(t))
gec = lambda t : 6.0*pi*me*c/(sigT*BR(t)**2.0*gam3(t)*t)
num = lambda t : q*BR(t)/(2*np.pi*me*c)*gem(t)**2*gam3(t)
nuc = lambda t : q*BR(t)/(2.0*np.pi*me*c)*gec(t)**2.0*gam3(t)
Fmax = lambda t : Ne3(t)*me*c**2*sigT/3.0/q*BR(t)*gam3(t)/(4.0*pi*DL**2) # Sari 97

gem_f = lambda t : epseF*(p-2)/(p-1)*mp/me*gam2(t)
B_f = lambda t : np.sqrt(32.0*pi*mp*epsBF*n)*gam2(t)*c
gec_f = lambda t : 6*pi*me*c/(sigT*gam2(t)*B_f(t)**2*t)
num_f = lambda t : gam2(t)*gem_f(t)**2*q*B_f(t)/(2*pi*me*c)
nuc_f = lambda t : gam2(t)*gec_f(t)**2*q*B_f(t)/(2*pi*me*c)
R_f = lambda t : 2.0*gam2(t)**2.0*c*t
#R_f = lambda t : (17.0*E0*t/(4*pi*mp*n*c))**(1./4)
Ne2 = lambda t : 4.0*pi*R_f(t)**3*n/3.0
Pmax_f = lambda t : me*c**2*sigT*gam2(t)*B_f(t)/3/q
Fmax_f = lambda t : Ne2(t)*Pmax_f(t)/(4*pi*DL**2)

bet3 = lambda t : np.sqrt(gam3(t)*gam3(t)-1.0)/gam3(t)
gam34 = lambda t : gam3(t)*gam4*(1.0-bet3(t)*bet4)


#FmaxIC(t):
#  return sigT*number(t)/(4*math.pi*radius(t)**2)*Fmax(t)
#numIC(t):
#  return 2*gem(t)**2*num(t)
#nucIC(t):
#  return 2*gec(t)**2*nuc(t)

def rspe(t,u):
  if num(t)<nuc(t): # Fast cooling 
    #print 'fast cooling', u, num(t), nuc(t)
    #print 'Fmax',Fmax(t), 1.1e5*epsBR**(1./2)*E52*n**(1./2)*D28**(-2)
    return np.where(u<num(t),(u/num(t))**(1./3)*Fmax(t),np.where(u<nuc(t),(u/num(t))**(-(p-1.)/2)*Fmax(t),(u/nuc(t))**(-p/2)*(nuc(t)/num(t))**(-(p-1.)/2)*Fmax(t)))
  else:
    return np.where(u<nuc(t),(u/nuc(t))**(1./3)*Fmax(t),np.where(u<num(t),(u/nuc(t))**(-1./2)*Fmax(t),(u/num(t))**(-p/2)*(num(t)/nuc(t))**(-1./2)*Fmax(t)))

def fspe(t,u):
  if num_f(t)<nuc_f(t): 
    return np.where(u<num_f(t),(u/num_f(t))**(1./3)*Fmax_f(t),np.where(u<nuc_f(t),(u/num_f(t))**(-(p-1.)/2)*Fmax_f(t),(u/nuc_f(t))**(-p/2)*(nuc_f(t)/num_f(t))**(-(p-1.)/2)*Fmax_f(t)))
  else:
    return np.where(u<nuc_f(t),(u/nuc_f(t))**(1./3)*Fmax_f(t),np.where(u<num_f(t),(u/nuc_f(t))**(-1./2)*Fmax_f(t),(u/num_f(t))**(-p/2)*(num_f(t)/nuc_f(t))**(-1./2)*Fmax_f(t)))

#------------IC code need check
#
#def rspeIC(t,u):
#  if numIC(t)<nucIC(t):
#    return np.where(u<numIC(t),(u/numIC(t))**(1./3)*FmaxIC(t),np.where(u<nucIC(t),(u/numIC(t))**(-(p-1.)/2)*FmaxIC(t),(u/nucIC(t))**(-p/2)*(nucIC(t)/numIC(t))**(-(p-1.)/2*FmaxIC(t)))*u
#  else:
#    return np.where(u<nucIC(t),u/nucIC(t)**(1./3)*FmaxIC(t),np.where(u<numIC(t),(u/nucIC(t))**(-1./2)*FmaxIC(t),(u/numIC(t))**(-p/2)*(numIC(t)/nucIC(t))**(-1.2)*FmaxIC(t)))*u



c = 299792458e2
hp = 6.626*1e-27
hb = 1.055*1e-27
kb = 1.381*1e-16
g = 6.673*1e-8
me = 9.1093897e-28
mp = 1.6726231e-24
q = 4.803*1e-10
sigT = 6.652*1e-25
pi = np.pi

# The evolution of the characteristic frequencies
p = 2.5
r014 = 1
E53 = 1e-1
E52 = E53*10
gam42 = 3.0
Delt12 = 1e1
epsBR2 = 1
epseR1 = 1
#DLG = 1 

# model parameters
r0 = r014*10**14
E0 = E53*10**53
gam4 = gam42*10**2
bet4 = np.sqrt(gam4*gam4-1.0)/gam4

Delt0 = Delt12*10**12
epseR = epseR1*0.1
epsBR = epsBR2*0.01
epseF = 1*epseR
epsBF = 1*epsBR
k = 0
n0 =n1 = 1.0
A = 3.e34 # torch
D28 = 1e0
DL = D28*10**28
nu_obs = 2e14
#1Gpe=3.085*10**27cm#

if k == 0:
  n = n0
else:
  n = A/r0**k

g2x1 = 2
g2x2 = 16
g2y1 = 4
g2y2 = g2y1+(g2x2-g2x1)*(-3./4)
g2mx1 = (g2x2+g2x1)/2
g2my1 = (g2y2+ g2y1)/2

N0 = E0/(gam4*mp*c**2)
#SedL = ((3-k)*E0/(4*np.pi*A*mp*c**2))**(1./(3-k)) # Gao 2015
SedL = ((3-k)*E0/(4*np.pi*n1*mp*c**2))**(1./(3-k))

ttw = Delt0/c   # RRS Crossing Time for Thick shell #Gao 2015
#ttw = SedL/(2*gam4**(8/3)) # NRS Crossing Time 
#ttn = SedL/(2*c*gam4**(3./8))  #Crossing Time for Thin shell
#ttn = SedL/(2**(2./3)*c*gam4**(8./3)) # Qiang Chen choose this equation
ttn = SedL/(2**(k-2)*c**(k-3)*gam4**(2*k-8)) # Gao 2015
print('Sedov Length = ',SedL)
Reta = SedL/gam4**(2./3)

xi = np.sqrt(SedL/Delt0)*gam4**(-(4.-k)/(3-k)) # Gao 2015
print('xi = %e' % xi)

#if Delt0 > SedL/(gam4**(8./3)):
if xi < 1.0:
  txd = ttw
  print('THICK SHELL, Tx = %e' % ttw)
else: 
  txd = ttn
  print('THIN SHELL, Tx = %e' % ttn)

Rcross = SedL**(3./4)*Delt0**(1./4)
Rdec = SedL/(gam4**(2./3))
print('ZK07 cross time',Rcross, Rdec, (Rcross-Rdec)/c)
#exit()


if (plot_ana):
  t = np.logspace(-2,8,300)
  plt1.plot(t,gam2(t),color ='black',lw=2,label=r'Ana $\sigma=0$')
  plt1.plot([txd,txd],[1e-1,gam2(txd)], color ='black', linewidth=1, linestyle="--",alpha=0.5)
  #plt1.scatter([txd,],[gam2(txd)], 150, color ='red')
  plt1.text(txd*1.1,1.1e0,r'$t_\times$')

  plt2.plot(t,R_f(t),color ='black',lw=2)#,label=r'$\sigma=0$')
  #plt2.axhline(y=Rcross, color='k', lw=1, linestyle='--')
  plt2.plot([1e0,txd],[Rcross,Rcross], color='k', lw=1, linestyle='--',alpha=0.5)
  plt2.text(2e0,Rcross*1.1,r'$R_\Delta$')

  plt3.plot(t,Ne3(t)*mp,color ='black',lw=2)#,label=r'$\sigma=0$')
  plt4.plot(t,n3(t),color ='black',lw=2)#,label=r'$\sigma=0$')
  

if (plot_model):

    print('reading FS dyn data')
    f = open(PATH+'/data/dynamics/dynfs.dat','r')
    dynfs = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    t2 = dynfs[:,0]
    a = int(0)
    b = int(len(t2))
    t2 = dynfs[a:b,0]
    gam2 = dynfs[a:b,1]
    m2 = dynfs[a:b,2]
    r2 = dynfs[a:b,3]
    Ek2 = dynfs[a:b,4]
    U2 = dynfs[a:b,5]
    mom2 = dynfs[a:b,6]
    print('reading RS dyn data')
    f = open(PATH+'/data/dynamics/dynrs.dat','r')
    dynrs = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    t3 = dynrs[:,0]
    a = int(0)
    b = int(len(t3))
    t3 = dynrs[a:b,0]
    gam3 = dynrs[a:b,1]
    m3 = dynrs[a:b,2]
    r3 = dynrs[a:b,3]
    Ek3 = dynrs[a:b,4]
    U3 = dynrs[a:b,5]
    EB3 = dynrs[a:b,6]

    mom3 = dynrs[a:b,7]
    mom4 = dynrs[a:b,8]
    E4 = dynrs[a:b,9]
    print('reading FS density data')
    f = open(PATH+'/data/densities/denfs.dat','r')
    denfs = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    n1 = denfs[a:b,0]
    n2 = denfs[a:b,1]
    p2 = denfs[a:b,2]
    e2 = denfs[a:b,3]
    hgam2 = denfs[a:b,4]
    print('reading RS density data')
    f = open(PATH+'/data/densities/denrs.dat','r')
    denrs = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    n4 = denrs[a:b,0]
    n3 = denrs[a:b,1]
    p3 = denrs[a:b,2]
    e3 = denrs[a:b,3]
    hgam3 = denrs[a:b,4] 
    E2 = Ek2+U2
    E3 = Ek3+U3+EB3

    #c = len(np.argwhere(t3<1e5)) # mask late unrealistic RS
    #U3[c:]=U3[c]
    #E3[c:]=E3[c]

    print('concatenate data')
    tcon = np.unique(np.concatenate((t3,t2)))
    E2c = np.interp(tcon, t2, E2, left=0, right=0)
    E3c = np.interp(tcon, t3, E3, left=0, right=0)
    E4c = np.interp(tcon, t2, E4, left=0, right=0)
    mom2c = np.interp(tcon, t2, mom2, left=0, right=0)
    mom3c = np.interp(tcon, t3, mom3, left=0, right=0)
    mom4c = np.interp(tcon, t3, mom4, left=0, right=0)
    mom2c[mom2c<0]=0
    mom3c[mom3c<0]=0
    mom4c[mom4c<0]=0
    momtot = mom2c+mom3c+mom4c
    E2c[E2c<0]=0
    E3c[E3c<0]=0
    E4c[E4c<0]=0
    Etot = E2c+E3c+E4c

    print('plotting')
    plt1.plot(t2,gam2,ls='-',c='b',label=r'$\gamma_2$')
    plt1.plot(t3,gam3,ls='--',c='m',label=r'$\gamma_3$')
    plt2.plot(t2,r2,ls='-',c='b',label=r'$R$')
    plt3.plot(t3,m3,ls='--',c='m',label=r'$m_3$')
    plt4.plot(t2,n2,ls='-',c='b',label=r'$n_2$')
    plt4.plot(t3,n3,ls='--',c='m',label=r'$n_3$')
    plt5.plot(t2,p2,ls='-',c='b',label=r'$p_2$')
    plt5.plot(t3,p3,ls='--',c='m',label=r'$p_3$')
    plt6.plot(t2,mom2/momtot[0],ls='-',c='b',lw=1,label=r'mom2')
    plt6.plot(t3,mom3/momtot[0],ls='--',c='m',lw=1,label=r'mom3')
    plt6.plot(t3,mom3/momtot[0],ls=':',c='m',lw=1,label=r'mom4')
    plt6.plot(tcon,momtot/momtot[0],ls='-',lw=2,c='k',label=r'momtot')
    plt7.plot(t2,E2/Etot[0],ls='-',c='b',lw=2,label=r'E2')
    plt7.plot(t2,U2/Etot[0],ls='--',c='b',lw=1,label=r'U2')
    plt7.plot(t2,Ek2/Etot[0],ls='-.',c='b',lw=1,label=r'E2k')
    plt7.plot(t3,E3/Etot[0],ls='-',c='m',lw=2,label=r'E3')
    plt7.plot(t3,U3/Etot[0],ls='--',c='m',lw=1,label=r'U3')
    plt7.plot(t3,Ek3/Etot[0],ls='-.',c='m',lw=1,label=r'Ek3')
    plt7.plot(t3,EB3/Etot[0],ls=':',c='m',lw=1,label=r'EB3')
    plt7.plot(t2,E4/Etot[0],ls='-',c='g',lw=2,label=r'E4')
    plt7.plot(tcon,Etot/Etot[0],ls='-',c='k',lw=2,label=r'Etot')
    plt8.plot(t2,hgam2,ls='-',c='g',lw=1,label=r'hgam2')
    plt8.plot(t2,hgam3,ls='--',c='g',lw=1,label=r'hgam3')


plt1.set_title(r'FS Lorentz factor')
plt1.set_ylabel(r'$\rm \gamma_2$')
plt1.axis([1e0,1e8,1e0,3e2])
plt1.set_xscale('log')
plt1.set_yscale('log')
#plt1.set_xlabel('t (s)')
plt1.legend(loc=0,frameon=False,prop={'size': 11},ncol=1)

plt2.set_title(r'radius, k=0')
plt2.set_ylabel(r'R (cm)')
plt2.axis([1e0,1e8,3e15,1e18])
plt2.set_xscale('log')
plt2.set_yscale('log')
#plt2.set_xlabel('t (s)')
plt2.legend(loc=0,frameon=False,prop={'size': 10},ncol=2)

plt3.set_title(r'RS shocked mass')
plt3.set_ylabel(r'$\rm m_3$ (g)')
plt3.axis([1e0,1e8,1e26,1e29])
plt3.set_xscale('log')
plt3.set_yscale('log')
#plt3.set_xlabel('t (s)')
plt3.legend(loc=0,frameon=False,prop={'size': 11},ncol=2)

plt4.set_title(r'RS number density')
plt4.set_ylabel(r'$\rm n_3$ ($\rm cm^{-3}$)')
plt4.axis([1e0,1e8,1e-1,1e6])
#plt4.set_xlim([1e0,1e8])
plt4.set_xscale('log')
plt4.set_yscale('log')
#plt4.set_xlabel('t (s)')
#plt3.grid()
plt4.legend(loc=0,frameon=True,prop={'size': 11},ncol=1)

plt5.set_title(r'pressure')
plt5.set_xscale('log')
plt5.set_yscale('log')
plt5.axis([1e0,1e8,1e-6,3e2])
plt5.legend(loc=0,frameon=True,prop={'size': 11},ncol=1)
plt5.set_xlabel('t (s)')

plt6.set_title(r'mometume conservation')
plt6.set_xscale('log')
plt6.set_yscale('log')
plt6.axis([1e0,1e8,1e-3,3e0])
plt6.legend(loc=0,frameon=True,prop={'size': 11},ncol=1)
plt6.set_xlabel('t (s)')

plt7.set_title(r'energy conservation')
plt7.set_xscale('log')
plt7.set_yscale('log')
plt7.axis([1e0,1e8,1e-3,3e0])
plt7.legend(loc=0,frameon=True,prop={'size': 11},ncol=2)
plt7.set_xlabel('t (s)')

plt8.set_title(r'adiabatic index')
plt8.set_xscale('log')
#plt8.set_yscale('log')
plt8.axis([1e0,1e8,4./3-0.1,5./3+0.1])
plt8.axhline(y=4./3, color='r', linestyle='--',lw=0.5)
plt8.axhline(y=5./3, color='r', linestyle='--',lw=0.5)
plt8.legend(loc=0,frameon=True,prop={'size': 11},ncol=1)
plt8.set_xlabel('t (s)')

plt1.set_title('(a)',loc='left')
plt2.set_title('(b)',loc='left')
plt3.set_title('(c)',loc='left')
plt4.set_title('(d)',loc='left')
plt5.set_title('(e)',loc='left')
plt6.set_title('(f)',loc='left')
plt7.set_title('(g)',loc='left')
plt8.set_title('(h)',loc='left')

plt1.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt2.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt3.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt4.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt5.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt6.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt7.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt8.tick_params(which='both', direction='in', right=True, top = True, labelright='off')

plt.tight_layout() 
plt.savefig(PATH+f"/figs/dyncom.pdf", format='pdf')
plt.savefig(PATH+f"/figs/dyncom.png",format='png')
plt.show()

