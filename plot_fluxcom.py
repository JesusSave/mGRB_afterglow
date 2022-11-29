#!/usr/bin/python3
# -*- coding: utf-8 -*
# Ref: 1507.01984v1, Gao & Meszaros 2015

import numpy as np
import matplotlib.pyplot as plt
cm = plt.cm.get_cmap('jet')
from scipy.integrate import quad

PATH = './'
PATH2 = '../Chen19v5_ISM_radio'
DFILE_lc = 'data/lightcurves'
DFILE_dyn = 'data/dynamics'
DFIG = 'figs'

plot_ana = True
plot_ana = False

plot_model = False
plot_model = True

plot_jump = True
plot_jump = False

plot_fs = False
plot_fs = True

plot_rs = False
plot_rs = True

plot_total_light_curve = False
plot_total_light_curve = True

sigma = 0.1
epsBrs = 0.01


def Gam2(t):
  if k == 0: # ISM
    if delt0 > SedL/(gam4**(8./3)): # Gao He 15, xi < 1, thick shell
      #return np.where(t<ttw,(SedL/delt0)**(3./8)*(4*t/ttw)**(-1./4),(17.0*E0/(1024.0*pi*n*mp*c**5*t**3))**(1./8))
      return np.where(t<ttw,(SedL/delt0)**(3./8)*(4*t/ttw)**(-1./4),(SedL/delt0)**(3./8)*(4*ttw/ttw)**(-1./4)*(t/ttw)**(-3./8))
    else:
      return np.where(t<ttn,gam4,(17.0*E0/(1024.0*pi*n*mp*c**5*t**3))**(1./8))
  else:
    if delt0 > SedL/(gam4**(8./3)): # Wind thick shell
      return np.where(t<ttw, 1.0/np.sqrt(2)*(SedL/delt0)**(1./4),(17.0*E0/(1024.0*pi*n*mp*c**5*t**3))**(1./8))
    else:
      return np.where(t<ttn,gam4,(17.0*E0/(1024.0*pi*n*mp*c**5*t**3))**(1./8))

def Gam3(t):
  if k == 0: # ISM
    if delt0 > SedL/(gam4**(8./3)): # Gao He 15, xi < 1, thick shell
      return np.where(t<ttw,(SedL/delt0)**(3./8)*(4*t/ttw)**(-1./4),(SedL/delt0)**(3./8)*(4*ttw/ttw)**(-1./4)*(t/ttw)**(-7./16))
    else:
      return np.where(t<ttn,gam4,gam4*(t/ttn)**(-2./5))
  else:
    if delt0 > SedL/(gam4**(8./3)): # Wind thick shell
      return np.where(t<ttw, 1.0/np.sqrt(2)*(SedL/delt0)**(1./4),1.0/np.sqrt(2)*(SedL/delt0)**(1./4)*(t/ttw)**(-3./8))
    else:
      return np.where(t<ttn,gam4,gam4*(t/ttn)**(-1./3))

def n3(t):
  if k ==0: #ISM
    if delt0 > SedL/(gam4**(8./3)): # thick
      return np.where(t<ttw,8*Gam3(t)**3*n/gam4,8*Gam3(ttw)**3*n/gam4*(t/ttw)**(-13./16))
  else:
    if delt0 > SedL/(gam4**(8./3)):
      return np.where(t<ttn,7*n*gam4**2*(t/ttn)**-3,7*n*gam4**2*(t/ttn)**(-6./7))
    else:
      return np.where(t<ttn,7.0*A*gam4**6/l**2*(t/ttn)**-3, 7.0*A*gam4**6/l**2*(t/ttn)**(-7.0/8))

def e3(t):
  if k ==0: #ISM
    if delt0 > SedL/(gam4**(8./3)): # thick
      return np.where(t<ttw,4*Gam3(t)**2*n*mp*c**2,4*Gam3(ttw)**2*n*mp*c**2*(t/ttw)**(-13./12))
    else:
      return np.where(t<ttn,4*gam4**2*n*mp*c**2,4*gam4**2*n*mp*c**2*(t/ttn)**(-8./7))
  else:
    if delt0 > SedL/(gam4**(8./3)): # Wind thick
      return np.where(t<ttw, 8.0*A*mp*c*c/SedL**(1./2)/delt0**(3./2)*(t/ttw)**(-2), 8.0*A*mp*c*c/SedL**(1./2)/delt0**(3./2)*(t/ttw)**(-3./2))
    else:
      return np.where(t<ttn,4.0*A*mp*c*c*gam4**6.0/SedL**2.0*(t/ttn)**(-2.),4.0*A*mp*c*c*gam4**6.0/SedL**2.0*(t/ttn)**(-32./21))

def Ne3(t):
  if k ==0: #ISM
    if delt0 > SedL/(gam4**(8./3)): # thick
      return np.where(t<ttw,N0*(t/ttw),N0)
    else:
      return np.where(t<ttn,N0*(t/ttn)**(3./2),N0)
  else:
    if delt0 > SedL/(gam4**(8./3)):
      return np.where(t<ttw,N0*(t/ttw),N0)
    else:
      return np.where(t<ttn,N0*(t/ttn)**(1./2),N0)

gem = lambda t : epseR*e3(t)/(n3(t)*mp*c**2)*(p-2)/(p-1)*mp/me
BR  = lambda t : np.sqrt(8.0*np.pi*epsBR*e3(t))
gec = lambda t : 6.0*pi*me*c/(sigT*BR(t)**2.0*Gam3(t)*t)
num = lambda t : Gam3(t)*gem(t)**2.*q*BR(t)/(2.0*np.pi*me*c)/(1.+z)
nuc = lambda t : Gam3(t)*gec(t)**2.*q*BR(t)/(2.0*np.pi*me*c)/(1.+z)
Pmax = lambda t: me*c**2*sigT/3.0/q*Gam3(t)*BR(t)
Fmax = lambda t : Ne3(t)*Pmax(t)/(4.0*pi*DL**2)*(1.+z) # Sari 97

gem_f = lambda t : epseF*(p-2)/(p-1)*mp/me*(Gam2(t)-1.)
B_f = lambda t : np.sqrt(32.0*pi*mp*epsBF*n)*Gam2(t)*c
gec_f = lambda t : 6.0*pi*me*c/(sigT*Gam2(t)*B_f(t)**2*t)
nuc_f = lambda t : Gam2(t)*gec_f(t)**2*q*B_f(t)/(2.0*np.pi*me*c)/(1.+z)
num_f = lambda t : Gam2(t)*gem_f(t)**2*q*B_f(t)/(2.0*np.pi*me*c)/(1.+z)
R_f = lambda t : 2.0*Gam2(t)**2.0*c*t
#R_f = lambda t : (17.0*E0*t/(4*pi*mp*n*c))**(1./4)
Ne2 = lambda t : 4.0*pi*R_f(t)**3*n/3.0
Pmax_f = lambda t : me*c**2*sigT*Gam2(t)*B_f(t)/3/q
Fmax_f = lambda t : Ne2(t)*Pmax_f(t)/(4*pi*DL**2)*(1.+z)

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


def func_Dc(x,a,b):
  return 10.0*c/H0/np.sqrt(a*(1.0+x)**3+b)

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
delt12 = 1e1
epsBR2 = 1
epseR1 = 1
#DLG = 1 
# model parameters
r0 = r014*10**14
E0 = E53*10**53
gam4 = gam42*10**2
delt0 = delt12*10**12
epseR = epseR1*0.1
epsBR = epsBR2*0.01
epseF = 1*epseR
epsBF = 1*epsBR
k = 0
n1 = 1.0
A = 3.e34 # torch
z = 1.0
H0 = 68.0 # [km/s/Mpc = 10^-1 m/s/pc]
OmegaM = 0.31 
OmegaL = 0.69
Dc = quad(func_Dc, 0.0, z, args=(OmegaM,OmegaL)) # [pc]
# Omega_k = 0, DM = Dc 
# DL = (1+z)DM ! unit pc
print ('Dc',Dc)
D28 = 3.085678*1.e18*(1.0+z)*Dc[0]/1e28
print ('D28',D28) 
DL = D28*10**28
nu_obs = 2e14 # optic
#nu_obs = 5e9 # radio
#1Gpe=3.085*10**27cm#

if k == 0:
  n = n1
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
ttw = delt0/c   # RRS Crossing Time for Thick shell #Gao 2015
#ttw = SedL/(2*gam4**(8/3)) # NRS Crossing Time 
#ttn = SedL/(2*c*gam4**(3./8))  #Crossing Time for Thin shell
#ttn = SedL/(2**(2./3)*c*gam4**(8./3)) # Qiang Chen choose this equation
ttn = SedL/(2**(k-2)*c**(k-3)*gam4**(2*k-8)) # Gao 2015
print('Sedov Length = ',SedL)
Reta = SedL/gam4**(2./3)
xi = np.sqrt(SedL/delt0)*gam4**(-(4.-k)/(3-k)) # Gao 2015
print('xi = %e' % xi)
#if delt0 > SedL/(gam4**(8./3)):
if xi < 1.0:
  txd = ttw
  print('THICK SHELL, Tx = %e' % ttw)
else: 
  txd = ttn
  print('THIN SHELL, Tx = %e' % ttn)

row = 2
col = 3

plt.clf()
plt.minorticks_on()
plt.gcf().set_size_inches(col*4.5,row*4.5)
plt.tight_layout(rect=[0.04, 0.03, 1, 0.99]) 

plt1 = plt.subplot(row,col,1)
plt2 = plt.subplot(row,col,2)
plt3 = plt.subplot(row,col,3)
plt4 = plt.subplot(row,col,4)
plt5 = plt.subplot(row,col,5)
plt6 = plt.subplot(row,col,6)

if (plot_ana):
  t = np.logspace(-2,8,100)
  fr = np.array([rspe(t1,nu_obs) for t1 in t])
  ff = np.array([fspe(t2,nu_obs) for t2 in t])
  plt1.plot(t,1e23*1e3*ff+1e23*1e3*fr,color='black',lw=2,label=r'ana')

if (plot_model):
    print('reading fluxfs data')
    f = open(PATH+f'/data/lightcurves/fluxfs.dat','r')
    fluxfs = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    t2 = fluxfs[:,0]
    a = int(0)
    b = int(len(t2))
    t2 = fluxfs[a:b,0]
    flux_FS = fluxfs[a:b,1]
    print('reading fluxRS data')
    f = open(PATH+f'/data/lightcurves/fluxrs.dat','r')
    fluxrs = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    t3 = fluxrs[:,0]
    a = int(0)
    b = int(len(t3))
    t3 = fluxrs[a:b,0]
    flux_RS = fluxrs[a:b,1]

    print(' plotting')
    if (plot_fs):
      plt1.plot(t2,1e23*1e3*flux_FS,c='b',lw=1,ls='-.',label=r'FS')
    if (plot_rs):
      plt1.plot(t3,1e23*1e3*flux_RS,c='b',lw=1,ls='--',label=r'RS')
    # concatenate RS and FS flux
    tcon = np.unique(np.concatenate((t3,t2)))
    flux_RS = np.interp(tcon, t3, flux_RS, left=0, right=0)
    flux_FS = np.interp(tcon, t2, flux_FS, left=0, right=0)
    # replace nan value with zero
    flux_RS[np.isnan(flux_RS)] = 0
    flux_FS[np.isnan(flux_FS)] = 0
    # replace negative value
    #flux_RS[flux_RS <= 0] = 1e-50
    #flux_FS[flux_FS <= 0] = 1
    if (plot_total_light_curve):
      plt1.plot(tcon,1e23*1e3*flux_RS+1e23*1e3*flux_FS,c='b',lw=2, label=r'tot')

    print('reading FS frequency data')
    f = open(PATH+f'/data/lightcurves/frequencyfs.dat','r')
    frefs = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    t2 = frefs[:,0]
    a = int(0)
    b = int(len(t2))        
    t2 = frefs[a:b,0]
    num2 = frefs[a:b,7]
    nuc2 = frefs[a:b,8]
    nua2 = frefs[a:b,9]
    Pmax2 = frefs[a:b,11]
    B2 = frefs[a:b,12]
    print('reading RS frequency data')
    f = open(PATH+f'/data/lightcurves/frequencyrs.dat','r')
    frers = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    t3 = frers[:,0]
    a = int(0)
    b = int(len(t3))        
    t3 = frers[a:b,0]
    num3 = frers[a:b,7]
    nuc3 = frers[a:b,8]
    nua3 = frers[a:b,9]
    Pmax3 = frers[a:b,11]
    B3 = frers[a:b,12]

    plt2.plot(t2, num2, c='b',lw=1,ls='--',label=r'$\nu_m^f$')  
    plt2.plot(t2, nuc2, c='b',lw=1,ls='-.',label=r'$\nu_c^f$')  
    plt2.plot(t2, nua2, c='b',lw=1,ls=':',label=r'$\nu_a^f$')   
    plt2.plot(t3, num3, c='g',lw=1,ls='--',label=r'$\nu_m^r$')  
    plt2.plot(t3, nuc3, c='g',lw=1,ls='-.',label=r'$\nu_c^r$')  
    plt2.plot(t3, nua3, c='g',lw=1,ls=':',label=r'$\nu_a^r$')
    plt4.plot(t2, B2, c='b',lw=1,ls='--',label=r'B2') 
    plt4.plot(t3, B3, c='g',lw=1,ls=':',label=r'B3') 
    plt5.plot(t2, Pmax2, c='b',lw=1,ls='--',label=r'Pmax2') 
    plt5.plot(t3, Pmax3, c='g',lw=1,ls=':',label=r'Pmax3') 

    print('reading magnetic data')
    f = open(PATH+f'/data/magnetization/mag.dat','r')
    mag = np.array([[float(data) for data in line.split()] for line in f.readlines()])
    f.close()
    t3 = mag[:,0]
    a = int(0)
    b = int(len(t3))  
    t3 = mag[a:b,0]
    fc = mag[a:b,6]
    hgam3 = mag[a:b,7]
    if(sigma==0):
      epsBrs = epsBrs
    else:
      epsBrs = (hgam3-1.)*(fc-1.)
    plt3.plot(t3, epsBrs, c='g',lw=1,ls=':',label=r'$\epsilon_B^r$') 



if (plot_jump):
    print('  Optic section')

plt1.set_title(r'total optic light curve')
plt1.axis([1e0,1e7,1e-4,1e2])
plt1.set_xscale('log')
plt1.set_yscale('log')
plt1.set_ylabel(r'Flux (mJy)')
#plt.grid()
plt1.legend(loc=0,frameon=False,prop={'size': 11},ncol=1)

plt2.set_title(r'characteristic frequency')
plt2.axis([1e0,1e7,1e4,1e20])
#plt2.set_xlim([1e0,1e7])
plt2.set_xscale('log')
plt2.set_yscale('log')
#plt2.set_ylabel(r'Flux (mJy)')
#plt.grid()
plt2.legend(loc=0,frameon=False,prop={'size': 11},ncol=3)

plt3.set_title(r'magnetic energy fraction')
plt3.axis([1e0,1e7,1e-1,1e1])
plt3.set_xscale('log')
plt3.set_yscale('log')
#plt3.set_ylabel(r'Flux (mJy)')
#plt.grid()
plt3.legend(loc=0,frameon=False,prop={'size': 11},ncol=1)

plt4.set_title(r'magnetic field strength')
plt4.axis([1e0,1e7,1e-2,1e1])
plt4.set_xscale('log')
plt4.set_yscale('log')
#plt4.set_ylabel(r'Flux (mJy)')
#plt.grid()
plt4.legend(loc=0,frameon=False,prop={'size': 11},ncol=3)
plt4.set_xlabel('t (s)')

plt5.set_title(r'peak spectral power')
#plt5.set_xlim([1e0,1e7])
plt5.axis([1e0,1e7,1e-23,1e-20])
plt5.set_xscale('log')
plt5.set_yscale('log')
plt5.set_xlabel('t (s)')

plt6.set_xlim([1e0,1e7])
plt6.set_xscale('log')
plt6.set_yscale('log')
plt6.set_xlabel('t (s)')

plt1.set_title('(a)',loc='left')
plt2.set_title('(b)',loc='left')
plt3.set_title('(c)',loc='left')
plt4.set_title('(d)',loc='left')
plt5.set_title('(e)',loc='left')
plt6.set_title('(f)',loc='left')

plt1.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt2.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt3.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt4.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt5.tick_params(which='both', direction='in', right=True, top = True, labelright='off')
plt6.tick_params(which='both', direction='in', right=True, top = True, labelright='off')

plt.savefig(PATH+f"/{DFIG}/flux_radio.pdf", format='pdf')
plt.savefig(PATH+f"/{DFIG}/flux_radio.png",format='png')
plt.show()

