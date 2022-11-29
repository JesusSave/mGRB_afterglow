#===================================================================
#           MHD shock jump condition dynamics evolution
#
# This program use MHD jump conditions to derive GRB afterglow dynamics
# The central assumption is equity of pressure and Lorentz factor at CD
#
# VERSION: v11 (3 Aug 2020) 
# 
# FEATURE:
# - dynamics format is consistent with TORCH
# - include a naive radiatition section
#
# AUTHORS: 
# - Mr. QC chen@camk.edu.pl
# - Prof. Dr. KN
#
# IMPORTANT: 
# - very slow in speed >_<
# - use python3 f-format
# - equaty pressure and Lorentz factor at CD vialate energy conservation
# - light curves is recommaned to plot by TORCH reading dynamics
#
# REMAINED:
# - cannot deal with radiative case, eps2=eps3=0 is required  @_@
#
# REFERENCE:
# - Zhang & Kobayashi 2005, astro-ph/0404140
# - Chen & Liu 2019, in prep.
#===================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.special import kn
from scipy.integrate import quad
cm = plt.cm.get_cmap('jet')

PATH = '.'

def mshock(n4,gam4,gam3,sigma):
  global MAX_ITERATIONS # def can read global data, must declare in case of modifing

  bet3 = np.sqrt(gam3*gam3-1.0)/gam3
  bet4 = np.sqrt(gam4*gam4-1.0)/gam4
  gam34 = gam3*gam4*(1.0-bet3*bet4)
  bet34 = (bet3-bet4)/(1.0-bet4*bet3)

  #calculation of adiabatic index hgam from relativistic temperature T = p / (n m c^2)
  hgam = hgam1
  hgam0 = hgam1

  i = 0
  #while (True):
  for i in range(100):
    i += 1
    '''
    #-----------------------Z&K05's Final Results---------------
    a1  =  hgam*(2.0-hgam)*(gam34-1.0)+2.0
    b1  =  -(gam34+1.0)*((2.0-hgam)*(hgam*gam34**2.0+1.0)+hgam*(hgam-1.0)*gam34)*sigma-(gam34-1.0)*(hgam*(2.0-hgam)*(gam34**2.0-2.0)+(2.0*gam34+3.0))
    c1  =  (gam34+1.0)*(hgam*(1.0-hgam/4.0)*(gam34**2.0-1.0)+1.0)*sigma**2.0+ (gam34**2.0-1.0)*(2.0*gam34-(2.0-hgam)*(hgam*gam34-1.0))*sigma + (gam34+1.0)*(gam34-1.0)**2*(hgam-1.0)**2.0
    d1  =  -(gam34-1.0)*(gam34+1.0)**2.0*(2.0-hgam)**2.0*sigma**2.0/4.0
    '''
    #-------------Chen's Version, identical to Z&K05-------------
    a1  = hgam*(hgam-2.0)*gam34**2.0-2.0*(hgam-1.0)**2.0*gam34+(hgam-1.0)**2.0+1.0
    b1  = hgam*(2.0-hgam)*(1.0+sigma)*gam34**4.0  +((2.0+sigma)*hgam-2.0)*(hgam-1.0)*gam34**3.0  +((1.0+sigma)*hgam**2.0-(3.0*sigma+2.0)*hgam+2.0*sigma-1.0)*gam34**2.0  +(hgam-1.0)*(4.0-hgam*(sigma+4.0))*gam34+2.0*hgam**2.0+(sigma-4.0)*hgam-2.0*sigma+3.0
    c1 = (0.25*(sigma**2.0-4.0*sigma-4.0)*hgam**2.0-(sigma**2.0-2.0*sigma-2.0)*hgam-2.0*sigma-1.0)*gam34**4.0  +(hgam-1.0)*(hgam*(sigma+2.0)-2.0)*gam34**3.0 +(-0.5*sigma*(sigma-2.0)*hgam**2.0+sigma*(2.0*sigma-3.0)*hgam-sigma*(sigma-4.0))*gam34**2.0 +(hgam-1)*(2.0-sigma*hgam-2.0*hgam)*gam34 +0.25*(sigma**2.0+4.0)*hgam**2.0-(sigma**2.0-sigma+2.0)*hgam+(sigma-1.0)**2.0
    d1  = (gam34**2-1.0)**2.0*(2.0-hgam)**2.0*sigma**2.0/4.0

    coeff = [a1,b1,c1,d1]
    u3psqr = np.roots(coeff)
    u3p = np.sqrt(u3psqr[1])*np.sign(-bet34)
    bet3p = u3p/np.sqrt(u3p**2.0+1.0)

    betRS = (bet3-bet3p)/(1.0-bet3*bet3p)
    GRS = 1.0/np.sqrt(1.0-betRS**2.0)
    uRS = GRS*betRS

    bet4p = (bet3p-bet34)/(1.0-bet3p*bet34)
    gam4p = 1.0/np.sqrt(1.0-bet4p**2.0)
    u4p = gam4p*bet4p
    #betRS_v2 = (bet4-bet4p)/(1.0-bet4*bet4p)

    if (u4p==0 ) or (u3p ==0):
      fa = 1.0
      fb = 1.0
      fc = 1.0
    else:
      fa = 1.0 - (gam34 +1.0)/(2.0*(u3p**2.0*gam34 +u3p*np.sqrt(u3p**2.0+1.0)*np.sqrt(gam34*gam34-1.0)))*sigma
      fb = (gam34+np.sign(u3p)*(np.sqrt(u3p**2.0+1.0)/u3p)*np.sqrt(gam34**2.0-1.0))/((hgam*gam34+1.0)/(hgam-1.0))
      fc = 1.0 + 1.0/(2.0*(hgam-1.0))*u4p/u3p/((gam34-1.0)*fa)*sigma
    n3 = n4*(hgam*gam34+1.0)/(hgam-1.0)*fb
    e3 = n3*mp*c*c*(gam34-1.0)*fa
    p3 = (hgam-1.0)*e3*fc
    T3 = p3/(n3*mp*c*c)
    if (T3 > 2.0e-3):
      hgam = 1.0+T3/(3.0*T3+kn(1,1.0/T3)/kn(2,1.0/T3)-1.0)
    else:
      hgam = 5.0/3.0
    if (abs(hgam - hgam0) < 1.0e-8):
      if (i > MAX_ITERATIONS):
        MAX_ITERATIONS = i
      break
    hgam0 = hgam

  return uRS,u3p,gam34,fa,fb,fc,n3,e3,p3,hgam


def loc_mshock():

  gam3_min = 1e0
  gam3_max = gam4
  Ngam3 = 100
  arr_gam3 = np.logspace(np.log10(gam3_min),np.log10(gam3_max),Ngam3)

  arr_dp = []
  arr_p2 = []
  arr_p3 = []

  for gam3 in arr_gam3:

    uRS,u3p,gam34,fa3,fb3,fc3,n3,e3,p3,hgam3 = mshock(n4,gam4,gam3,sigma)
    uFS,u2p,gam21,fa,fb,fc,n2,e2,p2,hgam2 = mshock(n1,gam1,gam3,0.0)
    dp = abs((p2-p3)/p3)
    arr_dp.append(dp)
    arr_p2.append(p2)
    arr_p3.append(p3)

  minp = np.argmin(arr_dp) 
  gam3_loc = arr_gam3[minp]


  if (np.min(dp) > 1e0):
    arr_gam3_refine = arr_gam3

    for levels in range(30):
      if (minp < 5):
        arr_gam3_refine = np.logspace(np.log10(arr_gam3[0]),np.log10(arr_gam3[9]),Ngam3)
      elif (minp > 90):
        arr_gam3_refine = np.logspace(np.log10(arr_gam3[90]),np.log10(arr_gam3[99]),Ngam3) 
      else:
        arr_gam3_refine = np.logspace(np.log10(arr_gam3[minp-5]),np.log10(arr_gam3[minp+5]),Ngam3) 

      arr_dp = []
      arr_p2 = []
      arr_p3 = []
      for gam3 in arr_gam3_refine:
        uRS,u3p,gam34,fa3,fb3,fc3,n3,e3,p3,hgam3 = mshock(n4,gam4,gam3,sigma)
        uFS,u2p,gam21,fa,fb,fc,n2,e2,p2,hgam2 = mshock(n1,gam1,gam3,0.0)
        dp = abs((p2-p3)/p3)
        arr_dp.append(dp)
        arr_p2.append(p2)
        arr_p3.append(p3)
      minp = np.argmin(arr_dp) 
      gam3_loc = arr_gam3_refine[minp]
      arr_gam3 = arr_gam3_refine

      if (np.min(dp)<1e-10):
        break
    #print ('      p2-p3',arr_dp[minp], minp)
  return gam3_loc


def Flux(num, nuc, Fmax, t , u):
  if num < nuc: # Fast cooling 
    return np.where(u<num,(u/num)**(1./3)*Fmax,np.where(u<nuc,(u/num)**(-(p-1.)/2)*Fmax,(u/nuc)**(-p/2)*(nuc/num)**(-(p-1.)/2)*Fmax))
  else:
    return np.where(u<nuc,(u/nuc)**(1./3)*Fmax,np.where(u<num,(u/nuc)**(-1./2)*Fmax,(u/num)**(-p/2)*(num/nuc)**(-1./2)*Fmax))


def fs_after_cross(v,x):
  r  = np.power(10.0,x)   
  gam2 = v[0]
  m2 = v[1]
  t2 = v[2]
  drdx = r*np.log(10.0)
  gam2dx = -4.0*np.pi*r*r*n1r(r)*mp*(gam2*gam2-1.0)/((2.0*gam2-2.0*eps2*gam2+eps2)*m2)*drdx
  m2dx = 4.0*np.pi*(r**2.0)*n1r(r)*mp*drdx
  u2 = np.sqrt(gam2*gam2-1.0)
  beta2 = u2/gam2
  dtdx = (1.0-beta2)/(beta2*c)*drdx
  return [gam2dx, m2dx, dtdx]

def n1r(r):
  if k ==0:
    n1 = n0
  else:
    n1 = Aast*3e35/r**k
  return n1

def func_Dc(x,a,b):
  return 10.0*c/H0/np.sqrt(a*(1.0+x)**3+b)


def adiabatic_index(p,n):
  T = p/(n*mp*c*c)
  h = kn(1,1.0/T)/kn(2,1.0/T)
  hgam = 1.0 + T/(3.0*T + h -1.0)
  return hgam

def adiabatic_index_hydro(e,n,hgam1):
  hgam = hgam1
  hgam0 = hgam1
  i = 0
  for i in range(100):
    i += 1
    p = (hgam-1.0)*e
    T = p/(n*mp*c*c)
    if (T > 2.0e-3):
      hgam = 1.0+T/(3.0*T+kn(1,1.0/T)/kn(2,1.0/T)-1.0)
    else:
      hgam = 5.0/3.0
    if (abs(hgam - hgam0) < 1.0e-8):
      break
    hgam0 = hgam
  return hgam

# CONSTANTS
c = 299792458e2
me = 9.1093897e-28
mp = 1.6726231e-24
q = 4.803*1e-10
sigT = 6.652*1e-25
# SIMULATION RESOLUTION
nstep = 100
# BEGIN KEY PARAMETERS
gam1 = 1.0
gam4 = 3e2
p = 2.5
n1 = 1.0
n0 = 1.0
p1 = 0.0
Aast = 1e-2
p4 = 0.0
sigmas = [0,1e-2,1e-1,1e0,1e1,1e2,1e3]
sigmas = [0.1]
Delta0 = 1.e13
Eiso = 1.0
k = 0

epsBR2 = 1
epseR1 = 1
epseR = epseR1*0.1
epsBR0 = epsBR = epsBR2*0.01
epseF = 1*epseR
epsBF = 1*epsBR
eps2 = 0.0
eps3 = 0.0

z = 1.0
H0 = 68.0
OmegaM = 0.31
OmegaL = 0.69
Dc = quad(func_Dc, 0.0, z, args=(OmegaM,OmegaL))
# Omega_k = 0, DM = Dc 
# DL = (1+z)DM ! unit pc
print ('Dc',Dc)
D28 = 3.085678*1.e18*(1.0+z)*Dc[0]/1e28
print ('D28',D28 )
DL = D28*10**28

nu_obs = 2e14
#nu_obs = 5e9
#1Gpe=3.085*10**27cm#
#nu_obs = c/1928.0/1e-8
DL = D28*10**28

Del12 = Delta0/10**12.
gam2p5 = gam4/10**2.5
E53 = Eiso*10**52/10**53

if (k==0):
  l = ((3.0-k)*Eiso*1e52/(4.0*np.pi*n0*mp*c*c))**(1./(3.0-k))
else:
  l = ((3.0-k)*Eiso*1e52/(4.0*np.pi*Aast*3e35*mp*c*c))**(1./(3.0-k))


rdec = l*gam4**(-2./(3.-k))
rs = Delta0*gam4*gam4
xi = np.sqrt(rdec/rs)

print ('Sedov length', l)
print ('xi', xi)
print ('Deceleration radius', rdec)

if Delta0 > l*gam4**(-8./3)/c:
  THINSHELL = False
  print ('  Thick Shell')
  print ('  Rdec_thick = Deltaa**1/4*l**3/4 =', Delta0**(1./4)*l**(3./4)) 
else:
  THINSHELL = True
  print ('  Thin Shell')
  print ('  Rdec_thin = l/ega**2./3 = ', l/(gam4**(2./(3.-k)))) 

for sigma in sigmas:
  print ('sigma',sigma, sigmas.index(sigma)+1,'/',len(sigmas))
  cnorm = (sigmas.index(sigma)+1.0)/len(sigmas) 


  MAX_ITERATIONS = 0

  # initialize hgam
  T4 = 0.0
  if (T4 > 2.0e-3):
      hgam4 = 1.0+T4/(3.0*T4+kn(1,1.0/T4)/kn(2,1.0/T4)-1.0)
  else:
      hgam4 = 5.0/3.0
  hgam = hgam1 = hgam4 

  Rdec = l*gam4**(-2./(3.0-k))*(1.0+sigma)**(-(3.0-k))
  print ('  Deceleration radius', Rdec)
  R0 = Delta0*gam4

  x0 = np.log10(R0)
  x1 = x0+6.0
  dx = (x1-x0)/nstep

  rs = Delta0*gam4*gam4
  xi = np.sqrt(Rdec/rs)

  t = t0 = R0/c
  bet4 = np.sqrt(gam4*gam4-1.0)/gam4
  tobs = t0*(1+z)*(1.-bet4)
  m2 = 4.0/3.0*np.pi*(R0**3.0)*n1r(R0)*mp
  m3 = 0.0

  mej = Eiso*1e52/(gam4*c*c*(1.0+sigma))
  n40 = mej/(4.0*np.pi*R0*R0*gam4*max(Delta0,R0/gam4**2)*mp)

  print ('  Sedov length', l)
  print ('  spreading radius', rs)
  print ('  xi', xi)
  print ('  Rdec', Rdec)

  f_dynf = open(PATH+f'/data/dynamics/dynfs.dat','w')
  f_dynr = open(PATH+f'/data/dynamics/dynrs.dat','w')
  f_mag = open(PATH+f'/data/magnetization/mag.dat','w')
  f_denf = open(PATH+f'/data/densities/denfs.dat','w')
  f_denr = open(PATH+f'/data/densities/denrs.dat','w')
  f_fref = open(PATH+f'/data/lightcurves/frequencyfs.dat','w')
  f_frer = open(PATH+f'/data/lightcurves/frequencyrs.dat','w')
  f_fluxf = open(PATH+f'/data/lightcurves/fluxfs.dat','w')
  f_fluxr = open(PATH+f'/data/lightcurves/fluxrs.dat','w')

  arr_x = list(np.arange(x0,x1,dx))
  for x in arr_x:
    #print ('  x=',x, arr_x.index(x)+1,'/',len(arr_x))
    r = np.power(10.0,x)    
    dr = r-np.power(10.0,x-dx)
    n4 = mej/(4.0*np.pi*r*r*gam4*max(Delta0,r/gam4**2)*mp)
    n1 = n1r(r)
    #print ('    n4 =',n4)

    if (m3 < mej):
      gam3 = loc_mshock()
      uRS,u3p,gam34,fa,fb,fc,n3,e3,p3,hgam3 = mshock(n4,gam4,gam3,sigma)
      uFS,u2p,gam21,fa2,fb2,fc2,n2,e2,p2,hgam2 = mshock(n1,gam1,gam3,0.0)
      gam2 = gam3

      u2 = np.sqrt(gam2*gam2-1.0)
      beta2 = np.sqrt(u2*u2/(1.0+u2*u2))
      dt = dr/beta2/c
      t += dt
      dtobs = dt*(1.0-beta2)*(1+z)
      tobs += dtobs
      tobs3 = tobs
      
      bet3p = u3p/np.sqrt(u3p*u3p+1.0)
      bet34 = np.sqrt(gam34*gam34-1.0)/gam34
      bet3 = np.sqrt(gam2*gam2-1.0)/gam2
      #betRS = uRS/np.sqrt(uRS*uRS+1.0)      
      betRS  = (bet3-bet3p)/(1.0-bet3*bet3p)

      m3 += 4.0*np.pi*(r**2.0)*(bet4-betRS)/betRS*gam4*n4*mp*dr
      m2 += 4.0*np.pi*(r**2.0)*n1*mp*dr

      gam3tx = gam3
      n3tx = n3
      e3tx = e3
      tx = tobs
      m3tx = m3
      fatx = fa
      fbtx = fb
      fctx = fc
      ftx = n4/n1
      u3ptx = u3p
      turn = arr_x.index(x)

    else:
      if (arr_x.index(x) -turn ==1):
        v_after = odeint(fs_after_cross, [gam3tx,m2,tx],arr_x[turn+1:])
      gam2 = v_after[:,0][arr_x.index(x)-turn-1]
      m2 = v_after[:,1][arr_x.index(x)-turn-1]
      tobs = v_after[:,2][arr_x.index(x)-turn-1]
      u3p = u3ptx
      fa = fatx
      fb = fbtx 
      fc = fctx

      ratio = 10**(x-arr_x[turn])                # ratio = r/rx
      if(k==0.):
        if(THINSHELL):
          g = 2.
          tobs3 = tx*ratio**(2.0*g+1.0)          # t3:(r/rx)**(2g+1)
          gam3 = gam3tx*ratio**(-g)              # gam3:(r/rx)**(-g)
          e3 = e3tx*ratio**(-8.0*(3.0+g)/7.0)
          n3 = n3tx*ratio**(-6.0*(3.0+g)/7.0)
        else: 
          g = 3.50
          tobs3 = tx*ratio**(2.0*g+1.0)           # t3:(r/rx)**(2g+1)
          gam3 = gam3tx*ratio**(-g)               # gam3:(r/rx)**(-g)
          n3 = n3tx*ratio**((2.0*g+1.0)*(-13./16))
          e3 = e3tx*ratio**((2.0*g+1.0)*(-13./12))
      elif(k==2.0): 
        if(THINSHELL):
          g = 1.0
          tobs3= tx*ratio**(2.0*g+1.0)            # t3:(r/rx)**(2g+1)
          gam3 = gam3tx*ratio**(-g)               # gam3:(r/rx)**(-g)
          e3 = e3tx*ratio**(-8.0*(3.0+g)/7.0)
          n3 = n3tx*ratio**(-6.0*(3.0+g)/7.0)
        else:
          g = 1.50
          tobs3 = tx*ratio**(2.0*g+1.0)           # t3:(r/rx)**(2g+1)
          gam3 = gam3tx*ratio**(-g)               # gam3:(r/rx)**(-g)
          n3 = n3tx*ratio**((2.0*g+1.0)*(-9./8))
          e3 = e3tx*ratio**((2.0*g+1.0)*(-3./2))
      else:
        print('cant deal with k=',k)
        exit()

      #if gam3 < 1:
      #  print 'break since gam3 =', gam3
      #  break

      hgam3 = hgam1
      hgam0 = hgam1
      i = 0
      for i in range(100):
        i += 1
        p3 = (hgam3-1.0)*e3*fc
        T3 = p3/(n3*mp*c*c)
        if (T3 > 2.0e-3):
          hgam3 = 1.0+T3/(3.0*T3+kn(1,1.0/T3)/kn(2,1.0/T3)-1.0)
        else:
          hgam3 = 5.0/3.0
        if (abs(hgam3 - hgam0) < 1.0e-8):
          if (i > MAX_ITERATIONS):
            MAX_ITERATIONS = i
          break
        hgam0 = hgam3

      hgam2 = hgam1
      hgam0 = hgam1
      i = 0
      for i in range(100):
        i += 1
        n2 = (hgam2*gam2+1.0)/(hgam2-1.0)*n1r(10**x)
        e2 = (gam2-1.0)*n2*mp*c*c
        p2 = (hgam2-1.0)*e2
        T2 = p2/(n2*mp*c*c)
        if (T2 > 2.0e-3):
          hgam2 = 1.0+T2/(3.0*T2+kn(1,1.0/T2)/kn(2,1.0/T2)-1.0)
        else:
          hgam2 = 5.0/3.0
        if (abs(hgam2 - hgam0) < 1.0e-12):
          if (i > MAX_ITERATIONS):
            MAX_ITERATIONS = i
          #print '    hgam = %g after %d iterations' % (hgam,i)
          break
        hgam0 = hgam2
      print('i',arr_x.index(x))
      print('x',x)
      print('gam2',gam2)
      print('hgam2',hgam2)

    u2 = np.sqrt(gam2*gam2-1.0)
    u3 = np.sqrt(gam3*gam3-1.0)
    u4 = np.sqrt(gam4*gam4-1.0)
    gam34 = gam4*gam3-u4*u3       
    mom2 = m2*u2
    mom3 = m3*u3
    mom4 = (mej-m3)*u4
    Ek2 = (gam2-1.)*m2*c*c 
    Eint2com = (gam2-1.)*m2*c*c
    Eint2 = (1.-eps2)*gam2*Eint2com
    E2 = Ek2+Eint2
    Ek3 = (gam3-1.)*m3*c*c
    Eint3com = fa*(gam34-1.)*m3*c*c
    Eint3 = (1.-eps3)*gam3*Eint3com
    EB3com = (hgam3-1.)*(fc-1.)*Eint3com
    EB3 = gam3*EB3com
    E3 = Ek3+Eint3+EB3
    E4 = (gam4-1.)*(1.+sigma)*(mej-m3)*c*c

    if sigma == 0:
      epsBR = epsBR0
    else:
      epsBR = (hgam3-1)*(fc-1.0)

    gem = epseR*e3/(n3*mp*c**2)*(p-2)/(p-1)*mp/me
    BR  = np.sqrt(8.0*np.pi*e3*epsBR)
    gec = 6.0*np.pi*me*c/(sigT*BR**2.0*gam3*tobs3)
    num = q*BR/(2.0*np.pi*me*c)*gem**2.0*gam3
    nuc = q*BR/(2.0*np.pi*me*c)*gec**2.0*gam3
    Ne3 = m3/mp
    Pmax = me*c**2*sigT/3.0/q*gam3*BR
    Fmax = Ne3*Pmax/(4.0*np.pi*DL**2.0) # Sari et al 97

    gem_f = epseF*(p-2)/(p-1)*mp/me*(gam2-1.)
    Bfs = np.sqrt(32.0*np.pi*mp*epsBF*n1r(r))*gam2*c
    gec_f = 6.0*np.pi*me*c/(sigT*gam2*Bfs**2.0*tobs)
    num_f = gam2*gem_f**2*q*Bfs/(2.0*np.pi*me*c)
    nuc_f = gam2*gec_f**2.0*q*Bfs/(2.0*np.pi*me*c)
    Ne2 = m2/mp
    Pmax_f = me*c**2.0*sigT*gam2*Bfs/3.0/q
    Fmax_f = Ne2*Pmax_f/(4.0*np.pi*DL**2)

    Flux_r = Flux(num, nuc, Fmax, tobs3 , nu_obs)
    Flux_f = Flux(num_f, nuc_f, Fmax_f, tobs , nu_obs)

    f_dynf.write(f'{tobs:g} {gam2:g} {m2:g} {r:g} {Ek2:g} {Eint2:g} {mom2:g}\n')
    f_dynr.write(f'{tobs3:g} {gam3:g} {m3:g} {r:g} {Ek3:g} {Eint3:g} {EB3:g} {mom3:g} {mom4:g} {E4:g}\n')
    f_mag.write(f'{tobs3:g} {gam3:g} {gam34:g} {u3p:g} {fa:g} {fb:g} {fc:g} {hgam3:g}\n')
    f_denf.write(f'{n1:g} {n2:g} {p2:g} {e2:g} {hgam2:g}\n')
    f_denr.write(f'{n4:g} {n3:g} {p3:g} {e3:g} {hgam3:g}\n')
    # RREQUENCY ITEMS: te,gam,R,Ntot*mp,gm,gc,gmax,num,nuc,nua,numax,Pmax,B
    # gmax, nua and numax are dummy!
    f_fref.write(f'{tobs:g} {gam2:g} {r:g} {m2:g} {gem_f:g} {gec_f:g} {gec_f:g} {num_f:g} {nuc_f:g} {nuc_f:g} {nuc_f:g} {Pmax_f:g} {Bfs:g}\n')
    f_frer.write(f'{tobs3:g} {gam3:g} {r:g} {m3:g} {gem:g} {gec:g} {gec:g} {num:g} {nuc:g} {nuc:g} {nuc:g} {Pmax:g} {BR:g}\n')
    f_fluxf.write(f'{tobs:g} %g\n'%(Flux_f))
    f_fluxr.write(f'{tobs3:g} %g\n'%(Flux_r))

  print ('ADIABATIC INDEX MAX_ITERATIONS =',MAX_ITERATIONS)
  f_dynf.close()
  f_dynr.close()
  f_mag.close()
  f_denf.close()
  f_denr.close()
  f_fref.close()
  f_frer.close()
  f_fluxf.close()
  f_fluxr.close()
