'''
Created on 16 May 2013

Waveguide dielectric simulator to generate S-parameters for given dielectric parameters

@author: Hardie Pienaar
'''

import numpy as np 
import matplotlib.pyplot as plt
import itertools
import readSParamFromFile as fileIO
import DielectricExtractionAlgorithms as dea
import parameterPlots as pp
from matplotlib.patches import Patch

if __name__ == '__main__':
    pass

#Show intro in command line
print("*******************************")
print("Wave Guide MUT Simulator")
print("*******************************")
print("")

#Define free-space constants
eps0 = 8.854187817e-12
mue0 = 4*np.pi*1e-7
c0 = 1/np.sqrt(eps0*mue0)

#Define input parameters
#S-band
# d = 25e-3
# a = 72.136e-3
# b = 34.036e-3
# startFreq = 2.6e9
# stopFreq = 3.2e9
#X-band
d = 5e-3
a = 22.86e-3
b = 10.16e-3
startFreq = 9e9
stopFreq = 9.1e9

epsr = 10
cond = 0;
muer = 1
noFreqPoints = 201

resolution = 1e-3 
d = 5*d

print("Waveguide specified parameters:")
print("Width                : "+str(a*1e3)+"mm")
print("Height               : "+str(b*1e3)+"mm")
print("Sample Holder Length : "+str(d*1e3)+"mm")
print("Epsillon Relative    : "+str(epsr))
print("Conductivity         : "+str(cond)+"S/m")
print("Mue Relative         : "+str(muer))
print("Start Frequency      : "+str(startFreq*1e-9)+"GHz")
print("Stop Frequency       : "+str(stopFreq*1e-9)+"GHz")
print("")

#Calculate Waveguide Parameters
freq = np.linspace(startFreq,stopFreq,noFreqPoints) #Generate frequency points
epsr_comp = epsr -1j*(cond/(freq*2*np.pi)) #Calculate complex epsilon


fc_air = 1/(2*a*np.sqrt(mue0*eps0)) #Calculate air fundamental mode cutoff frequency TE10
fc_MUT = 1/(2*a*np.sqrt(mue0*muer*eps0*epsr)) #Calculate MUT fundamental mode cutoff frequency TE10

lambda0 = c0/freq #Generate free space wavelengths
B0 = (2*np.pi)/lambda0; #Propagation constants for free space

B_air = B0*np.sqrt(1-np.power(fc_air/freq,2)) #Propagation constant in air part of wave guide for positive traveling waves
B_MUT = B0*np.sqrt(1-np.power(fc_MUT/freq,2)) #Propagation constant in MUT part of wave guide for positive traveling waves

Z_w_air = (2*np.pi*freq*mue0/B_air) #Wave impedance in air region
Z_w_MUT = (2*np.pi*freq*mue0*muer/B_MUT) #Wave impedance in MUT region

lambda_air = (2*np.pi)/B_air #Wave length in air region
lambda_MUT = (2*np.pi)/B_MUT #Wave length in air region

eps_prime = np.real(epsr)
eps_doublePrime = -1*np.imag(epsr) #Conductivity/Angular frequency
lossTangent = eps_doublePrime/eps_prime; #Calculate loss tangent

 
print("Waveguide calculated parameters:")
print("Cutoff frequency for air region        : "+str(fc_air*1e-9)+"GHz")
print("Cutoff frequency for MUT region        : "+str(fc_MUT*1e-9)+"GHz")
print("Propagation constant for free space    : "+str(B0[0])+" 1/m to "+str(B0[-1])+" 1/m")
print("Propagation constant for air region    : "+str(B_air[0])+" 1/m to "+str(B_air[-1])+" 1/m")
print("Propagation constant for MUT region    : "+str(B_MUT[0])+" 1/m to "+str(B_MUT[-1])+" 1/m")
print("Wave impedance for air region          : "+str(Z_w_air[0])+" Ohm to "+str(Z_w_air[-1])+" Ohm")
print("Wave impedance for MUT region          : "+str(Z_w_MUT[0])+" Ohm to "+str(Z_w_MUT[-1])+" Ohm")
print("Wave length for free space             : "+str(lambda0[0])+" m to "+str(lambda0[-1])+" m")
print("Wave length for air region             : "+str(lambda_air[0])+" m to "+str(lambda_air[-1])+" m")
print("Wave length for MUT region             : "+str(lambda_MUT[0])+" m to "+str(lambda_MUT[-1])+" m")
print("")

#Calculate E and H Fields for TE10 mode
EyMax = 1 #Unity input field

#Calculate directional propagation constants
Bx = np.pi/a #Propagation constant in x-direction, transverse to E field and propagation direction
By = 0  #Propagation constant in y-direction, In line with E-field
Bc = Bx + By #Bc is the cutoff wave number, Bc = B when Bz = 0

A10 = (eps0*a*EyMax)/np.pi #Calculate constant term
P = np.power((A10*Bc)/(eps0),2)*(B_air[0]/(2*2*np.pi*freq[0]*mue0))*(a/2)*(b/1) #Calculate Power Flow

plt.ion()
fig = plt.figure()
# plt.colorbar()

for phase in np.linspace(0,2*np.pi,100):
    #TODO: Incorporate testing variables
    #TODO: Add interface position to inputs
    zInt = 0.5*lambda_air[0]
    x = np.linspace(0,a,a/resolution)
    zAir1 = np.linspace(0,zInt,zInt/resolution)
    Ey = np.zeros(len(x)*len(zAir1),dtype=complex).reshape(len(x),len(zAir1))
    Hx = np.zeros(len(x)*len(zAir1),dtype=complex).reshape(len(x),len(zAir1))
    Hz = np.zeros(len(x)*len(zAir1),dtype=complex).reshape(len(x),len(zAir1))
    
    #Fields in air region 2D array with [X position] [Z position(In propagation direction)]
    for i in np.arange(0,len(x)-1): 
        Ey[i] = -(A10*np.pi/(eps0*a))*np.sin(np.pi*x[i]/a)*np.exp(-1j*(B_air[0]*zAir1 - phase))
        Hx[i] = A10*((B_air[0]*np.pi)/(2*np.pi*freq[0]*mue0*eps0*a))*np.sin(np.pi*x[i]/a)*np.exp(-1j*(B_air[0]*zAir1 - phase))
        Hz[i] = -1j*(A10/(2*np.pi*freq[0]*mue0*eps0))*np.power(np.pi/a,2)*np.cos(np.pi*x[i]/a)*np.exp(-1j*(B_air[0]*zAir1 - phase))
    
    #Integrate over the fields at the input port to get a numerical estimate of the power
#     Pint = 0
#     for i in np.arange(1,len(x)): 
#         Pint = Pint + 0.5*np.real(-1*Ey[i-1][0]*np.conjugate(Hx[i-1][0]))*(x[i]-x[i-1])*b   
        
    #Calculate Transmitted and reflected fields at Air-MUT interface
    tau = (2*Z_w_MUT[0])/(Z_w_air[0] + Z_w_MUT[0])
    rho = (Z_w_MUT[0] - Z_w_air[0])/(Z_w_air[0] + Z_w_MUT[0])
    EyMaxTrans = tau*np.abs(Ey[len(x)/2,len(zAir1)-1])
    EyMaxRefl = rho*np.abs(Ey[len(x)/2,len(zAir1)-1])
    EyPhaseAtInt = -np.angle(Ey[len(x)/2,len(zAir1)-1])-np.pi
    
    
    #Calculate reflected fields and fields in MUT
    zMUT = np.linspace(0,d,d/resolution)
    Eyt = np.zeros(len(x)*len(zMUT),dtype=complex).reshape(len(x),len(zMUT))
    Hxt = np.zeros(len(x)*len(zMUT),dtype=complex).reshape(len(x),len(zMUT))
    Hzt = np.zeros(len(x)*len(zMUT),dtype=complex).reshape(len(x),len(zMUT))
    
    
    
    #Fields in MUT region 2D array with [X position] [Z position(In propagation direction)]
    #TODO: These fields have have a infinite number of reflections
    A10t = (eps0*epsr*a*EyMaxTrans)/np.pi #Calculate constant term
    for i in np.arange(0,len(x)-1): 
        Eyt[i] = -(A10t*np.pi/(eps0*epsr*a))*np.sin(np.pi*x[i]/a)*np.exp(-1j*(B_MUT[0]*zMUT + EyPhaseAtInt))    
        Hxt[i,:] = A10t*((B_air[0]*np.pi)/(2*np.pi*freq[0]*mue0*eps0*epsr*a))*np.sin(np.pi*x[i]/a)*np.exp(-1j*(B_MUT[0]*zMUT+ EyPhaseAtInt ))
        Hzt[i,:] = -1j*(A10t/(2*np.pi*freq[0]*mue0*eps0*epsr))*np.power(np.pi/a,2)*np.cos(np.pi*x[i]/a)*np.exp(-1j*(B_MUT[0]*zMUT+ EyPhaseAtInt))
    
    #Calculate Transmitted and reflected fields at MUT-Air interface
    tau2 = (2*Z_w_air[0])/(Z_w_air[0] + Z_w_MUT[0])
    rho2 = (Z_w_air[0] - Z_w_MUT[0])/(Z_w_air[0] + Z_w_MUT[0])
    EyMaxTrans2 = tau*np.abs(Eyt[len(x)/2,len(zMUT)-1])
    EyMaxRefl2 = rho*np.abs(Eyt[len(x)/2,len(zMUT)-1])
    EyPhaseAtInt2 = -np.angle(Eyt[len(x)/2,len(zMUT)-1])-np.pi
    
    #Calculate reflected fields and fields in 2nd air region
    zAir2 = np.linspace(0,zInt,zInt/resolution)
    Eyt2 = np.zeros(len(x)*len(zAir2),dtype=complex).reshape(len(x),len(zAir2))
    Hxt2 = np.zeros(len(x)*len(zAir2),dtype=complex).reshape(len(x),len(zAir2))
    Hzt2 = np.zeros(len(x)*len(zAir2),dtype=complex).reshape(len(x),len(zAir2))
        
    #Fields in second air region 2D array with [X position] [Z position(In propagation direction)]
    #TODO: These fields have have a infinite number of transmissions
    A10t2 = (eps0*a*EyMaxTrans2)/np.pi #Calculate constant term
    for i in np.arange(0,len(x)-1): 
        Eyt2[i] = -(A10t2*np.pi/(eps0*a))*np.sin(np.pi*x[i]/a)*np.exp(-1j*(B_air[0]*zAir2 + EyPhaseAtInt2))
        Hxt2[i] = A10t2*((B_air[0]*np.pi)/(2*np.pi*freq[0]*mue0*eps0*a))*np.sin(np.pi*x[i]/a)*np.exp(-1j*(B_air[0]*zAir2+ EyPhaseAtInt2 ))
        Hzt2[i] = -1j*(A10t2/(2*np.pi*freq[0]*mue0*eps0))*np.power(np.pi/a,2)*np.cos(np.pi*x[i]/a)*np.exp(-1j*(B_air[0]*zAir2+ EyPhaseAtInt2))    
        
    
    
    print("Calculated field parameters: ")
    print("Input Voltage                         : "+str(EyMax)+" V/m")
    print("Field Constant at input               : "+str(A10))
    print("Field Constant at MUT trans           : "+str(A10t))
    print("Input Power                           : "+str(P*1e9)+" nW")    
    print("Input Power (Integrated)              : "+str(P*1e9)+" nW")
    print("Reflection coefficient at air-MUT     : "+str(rho))
    print("Transmission coefficient at air MUT   : "+str(tau))
    print("Reflected E-field at air-MUT          : "+str(EyMaxRefl)+" V/m")
    print("Transmitted E-field at air-MUT        : "+str(EyMaxTrans)+" V/m")
    print("E-field phase at air-MUT              : "+str(EyPhaseAtInt*(180/np.pi))+" deg")
    print("Reflection coefficient at MUT-air     : "+str(rho2))
    print("Transmission coefficient at MUT-air   : "+str(tau2))
    print("Reflected E-field at air-MUT          : "+str(EyMaxRefl2)+" V/m")
    print("Transmitted E-field at air-MUT        : "+str(EyMaxTrans2)+" V/m")
    print("E-field phase at air-MUT              : "+str(EyPhaseAtInt2*(180/np.pi))+" deg")
    
    #Append length vectors
    z = np.zeros((len(zAir1) + len(zMUT))+ len(zAir2))
    for i in np.arange(0,len(z)):
        if i < len(zAir1):
            z[i] = zAir1[i]
        elif i >= len(zAir1) and i < len(zAir1)+len(zMUT):
            z[i] = zMUT[i- len(zAir1)] + zAir1[len(zAir1) - 1]
        else:
            z[i] = zAir2[i- len(zAir1) - len(zMUT)] + zAir1[len(zAir1) - 1] + zMUT[len(zMUT)-1]
                    
    E = np.zeros(len(x)*(len(zAir1) + len(zMUT) + len(zAir2)),dtype=complex).reshape(len(x),len(zAir1) + len(zMUT)+ len(zAir2))
    for i in np.arange(0,len(z)):
        if i < len(Ey[0]):
            E[:,i] = Ey[:,i]
        elif i >=  len(Ey[0]) and i < len(Eyt[0]) + len(Ey[0]):
            E[:,i] = Eyt[:,i- len(Ey[0])]
        else:
            E[:,i] = Eyt2[:,i- len(Ey[0])- len(Eyt[0])]        
    
    # plt.figure()
    # plt.plot(z,np.real(E[len(x)/2]))
    # plt.axvline(zInt,linewidth=2,color='black')
    # plt.axvline(zInt+d,linewidth=2,color='black')
    # plt.figure()
    # plt.plot(x,np.real(Ey[:,0]))
    
    X,Y = np.meshgrid(z, x)
    plt.pcolor(X,Y,np.real(E))
    plt.axvline(zInt,linewidth=2,color='black')
    plt.axvline(zInt+d,linewidth=2,color='black')
    plt.xlim(0,z[len(z)-1])
    plt.ylim(0,x[len(x)-1])
       
    plt.draw()
# # call the animator.  blit=True means only re-draw the parts that have changed.
# from matplotlib import animation
# 
# anim = animation.FuncAnimation(fig, animate)    
# anim.save('basic_animation.mp4', fps=30)
