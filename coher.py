# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 22:46:44 2019

@author: charc
"""
import scipy.constants as sc
import numpy as np
import matplotlib.pyplot as plt

def coh(lat,natom,maxb,emax,flag):
    """
    Compute Bragg energies and associated 
    structure factors for coherent elastic 
    scattering.
    """
    k=1
    imax=1
    
    # Dictionary containing all the original material parameters from NJOY.
    materials=[
              {"a":2.4573e-8,"c":6.700e-8,"m":12.011,"x":5.5},
              {"a":2.2856e-8,"c":3.5832e-8,"m":9.01,"x":7.53},
              {"a":2.695e-8,"c":4.39e-8,"m":12.5,"x":1},
              {"a":4.04e-8,"m":26.7495,"x":1.495},
              {"a":4.94e-8,"m":207,"x":1},
              {"a":2.86e-8,"m":55.454,"x":12.9}
              ] 
    
    # Define two arrays for the energy and weight components of the bragg 
    # vector that will be built.
    energies=[0.0]*int(maxb/2)
    weights=[0.0]*int(maxb/2)
    
    # Shift physical constants to the units required for calculation
    eV=sc.eV*1e7
    hbar=sc.hbar*1e7
    neutron_mass=sc.neutron_mass*1e3
    atomic_mass=sc.atomic_mass*1e3
    econ=eV*8*neutron_mass/(hbar**2)
    recon=1/econ
    tsqx=econ/20
    
    # Constants hardcoded by NJOY
    twopis=39.478417604362633
    twothd=0.666666666667e0
    eps=0.05
    toler=1e-5
    wint=0
    
    # If the flag is true, use corrected values and redefine material 
    # parameters and constants.
    if flag:
        materials=[
              {"a":2.46e-8,"c":6.71e-8,"m":12.011,"x":5.551},
              {"a":2.226e-8,"c":3.57e-8,"m":9.0122,"x":7.63},
              {"a":2.71e-8,"c":4.402e-8,"m":12.5,"x":11.862},
              {"a":2.856e-8,"m":26.986,"x":1.495},
              {"a":3.571e-8,"m":207,"x":11.115},
              {"a":2.467e-8,"m":55.845,"x":11.22}
              ]    
        
        twopis=(2*np.pi)**2
        twothd=2/3
    
    # Define material parameters as variables for use in later calculations.
    a=materials[lat-1]["a"]
    amsc=materials[lat-1]["m"]
    scoh=materials[lat-1]["x"]/natom
    
    t2=hbar/(2*atomic_mass*amsc)
    ulim=econ*emax
    
    # For all hex materials
    if lat<4:
        c=materials[lat-1]["c"]
        c1=4/(3*a**2)
        c2=1/(c**2)
        scon=scoh*(4*np.pi)**2/(2*a**2*c*np.sqrt(3)*econ)
        phi=ulim/twopis
        i1m=int(a*np.sqrt(phi))
        i1m=i1m+1
        k=0
        for i1 in range(1,i1m):
          l1=i1-1
          i2m=int((l1+np.sqrt(3*(a*a*phi-l1*l1)))/2)
          i2m=i2m+1
          for i2 in range(i1,i2m):
             l2=i2-1
             x=phi-c1*(l1*l1+l2*l2-l1*l2)
             i3m=0
             if x > 0: i3m=int(c*np.sqrt(x))
             i3m=i3m+1
             for i3 in range(1,i3m):
                l3=i3-1
                w1=2
                if l1==l2: w1=1
                w2=2
                if l1==0 or l2 == 0: w2=1
                if l1==0 and l2 == 0: w2=1/2
                w3=2
                if l3==0: w3=1
                tsq=tausq(l1,l2,l3,c1,c2,twothd,twopis,lat)
                if tsq>0 and tsq<=ulim:
                   tau=np.sqrt(tsq)
                   w=np.exp(-tsq*t2*wint)*w1*w2*w3/tau
                   f=w*formf(lat,l1,l2,l3)
                   if k<=0 or tsq<=tsqx: 
                      k=k+1
                      if 2*k>maxb: print('coh' and
                        'storage exceeded ')
                      energies[k-1]=tsq
                      #print(k)
                      weights[k-1]=f 
                   else:
                      i=0
                      idone=0
                      while i<k and idone==0:
                         i=i+1
                         if tsq>=energies[i-1] and tsq<=(1+eps)*energies[i-1]:
                               weights[i-1]=weights[i-1]+f
                               idone=1
                      if idone==0:
                         k=k+1
                         if 2*k>maxb: print('coh' and
                           'storage exceeded ')
                         energies[k-1]=tsq
                         weights[k-1]=f
                tsq=tausq(l1,-l2,l3,c1,c2,twothd,twopis,lat)
                if tsq>0 and tsq<=ulim:
                   tau=np.sqrt(tsq)
                   w=np.exp(-tsq*t2*wint)*w1*w2*w3/tau
                   f=w*formf(lat,l1,-l2,l3)
                   if k<=0 or tsq<=tsqx: 
                      k=k+1
                      if 2*k> maxb: print('coh' and
                           'storage exceeded ')
                      energies[k-1]=tsq
                      weights[k-1]=f
                   else:
                      i=0
                      idone=0
                      while i<k and idone==0:
                         i=i+1
                         if tsq>=energies[i-1] and tsq<(1+eps)*energies[i-1]:
                            weights[i-1]=weights[i-1]+f
                            idone=1
                      if idone==0:
                         k=k+1
                         if 2*k>maxb: print('coh' and
                           'storage exceeded ')
                         energies[k-1]=tsq
                         weights[k-1]=f
        imax=k-1
       
    # For body centered cubic materials    
    elif lat<6:
       c1=3/(a*a)
       scon=scoh*(4*np.pi)**2/(16*a*a*a*econ)
       phi=ulim/twopis
       i1m=15
       k=0
       n=0
       for i1 in range(-i1m,i1m):
          i2m=i1m
          for i2 in range(-i2m,i2m):
             i3m=i1m
             for i3 in range(-i3m,i3m):
                tsq=tausq(i1,i2,i3,c1,0,twothd,twopis,lat)
                if tsq>0.0 and tsq<=ulim:
                   tau=np.sqrt(tsq)
                   w=np.exp(-tsq*t2*wint)/tau
                   f=w*formf(lat,i1,i2,i3)
                   k=k+1
                   if (2*k)> maxb: 
                       print('coh','storage exceeded',' ')
                   energies[k-1]=tsq
                   weights[k-1]=f
       imax=k
    
    else:
        c1=2/(a*a)
        scon=scoh*(4*np.pi)**2/(8*a*a*a*econ)
        phi=ulim/twopis
        i1m=int(a*np.sqrt(phi))
        i1m=15
        k=0
        for i1 in range(-i1m,i1m):
           i2m=i1m
           for i2 in range(-i2m,i2m):
              i3m=i1m
              for i3 in range(-i3m,i3m):
                 tsq=tausq(i1,i2,i3,c1,0,twothd,twopis,lat)
                 print(tsq)
                 if tsq>0 and tsq<=ulim:
                    tau=np.sqrt(tsq)
                    w=np.exp(-tsq*t2*wint)/tau
                    f=w*formf(lat,i1,i2,i3)
                    k=k+1
                    if (2*k)> maxb:
                        print('coh','storage exceeded',' ')
                    energies[k-1]=tsq
                    weights[k-1]=f
        imax=k-1      
    
    for i in range(0,imax):
        jmin=i+1
        for j in range(jmin,k):
            if energies[j]<energies[i]:
                st=energies[i]
                sf=weights[i]
                energies[i]=energies[j]
                weights[i]=weights[j]
                energies[j]=st
                weights[j]=sf
    k+=1
    energies[k]=ulim
    weights[k]=weights[k-1]
    maxb=2*k    
    bel=-1
    j=-1
    for i in range(0,k):
        be=energies[i]*recon
        bs=weights[i]*scon
        if (be-bel)<toler:
            weights[j]=weights[j]+bs
        else:
            j+=1
            energies[j]=be
            weights[j]=bs
            bel=be
    nbe=j
    maxb=2*nbe
    plotting(energies,weights,lat,flag)
    
def tausq(m1,m2,m3,c1,c2,twothd,twopis,lat):
    """
    Evaluate tausq equation
    """
    if lat < 4:
        return (c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*twopis
    elif lat < 6:
        return c1*(m1*m1+m2*m2+m3*m3+twothd*m1*m2+twothd*m1*m3-twothd*m2*m3)*twopis
    else:
        return c1*(m1*m1+m2*m2+m3*m3+m1*m2+m2*m3+m1*m3)*twopis
        

def formf(lat,l1,l2,l3):
    """
    Compute form factors for specified lattice
    """
    c1=7.54
    c2=4.24
    c3=11.31
    if lat==1:
        i=int(l3/2)  
        if (2*i)!=l3: 
            formf=np.sin(np.pi*(l1-l2)/3)**2
        else:
            formf=(6+10*np.cos(2*np.pi*(l1-l2)/3))/4
    elif lat==2:
        formf=1+np.cos(2*np.pi*(2*l1+4*l2+3*l3)/6)
    elif lat==3:
        formf=(1+np.cos(2*np.pi*(2*l1+4*l2+3*l3)/6))*(c1+c2+c3*np.cos(3*np.pi*l3/4))
    elif lat==4 or lat==5: 
        e1=2*np.pi*l1
        e2=2*np.pi*(l1+l2)
        e3=2*np.pi*(l1+l3)
        formf=(1+np.cos(e1)+np.cos(e2)+np.cos(e3))**2+(np.sin(e1)+np.sin(e2)+np.sin(e3))**2
    elif lat==6:
        e1=2*np.pi*(l1+l2+l3)
        formf=(1+np.cos(e1))**2+(np.sin(e1))**2
    return formf
    

def plotting(energies,weights,lat,flag):
    """
    Plot the bragg edges.
    """
    name=['G' 'Be' 'BeO' 'Al' 'Pb' 'Fe']
    fileName=name(lat-1)
    
    zero=energies.index(0)
    energies=energies[0:zero]
    weights=weights[0:zero]
    fullE=np.logspace(-3,0.5,1000)
    fullXS=[0.0]*len(fullE)
    for i in range(len(energies)):
        for j in range(len(fullE)):
            if fullE[j]>energies[i] and j<=len(fullE):
                  fullXS[j]+=np.exp(-4.0*energies[i])*weights[i]/fullE[j]
    plt.plot(fullE, fullXS,label = fileName)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energies (eV)')
    plt.ylabel('Cross Section (b)')
    plt.title(fileName +' Bragg Peaks')
    plt.legend(loc='best')
    #plt.savefig('temp.png')
    plt.show('temp.png')



    