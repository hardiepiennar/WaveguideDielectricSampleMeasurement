'''
Created on May 12, 2013

Parameter plot methods for different kinds of data

@author: Hardie Pienaar
'''
import numpy as np 

def plotSParam(plt,freq,S11,S21,S12,S22):
    plt.subplot(2,2,1)
    p1, = plt.plot(freq,20*np.log10(np.abs(S11)))
    p2, = plt.plot(freq,20*np.log10(np.abs(S22)))
    plt.title("Scattering Parameters")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude [dB]")
    plt.xlim([freq[0], freq[-1]])
    plt.legend([p1,p2],["S11","S22"])
    plt.grid(True)
    plt.subplot(2,2,2)
    p1, = plt.plot(freq,20*np.log10(np.abs(S21)))
    p2, = plt.plot(freq,20*np.log10(np.abs(S12)))
    plt.title("Scattering Parameters")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude [dB]")
    plt.xlim([freq[0], freq[-1]])
    plt.legend([p1,p2],["S21","S12"])
    plt.grid(True)
    plt.subplot(2,2,3)
    p1, = plt.plot(freq,np.angle(S11,True))
    p2, = plt.plot(freq,np.angle(S22,True))
    plt.title("Scattering Parameters")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase [Degrees]")
    plt.xlim([freq[0], freq[-1]])
    plt.legend([p1,p2],["S11","S22"])
    plt.grid(True)
    plt.subplot(2,2,4)
    p1, = plt.plot(freq,np.angle(S21,True))
    p2, = plt.plot(freq,np.angle(S12,True))
    plt.title("Scattering Parameters")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase [Degrees]")
    plt.xlim([freq[0], freq[-1]])
    plt.legend([p1,p2],["S21","S12"])
    plt.grid(True)
    return;

def plotDielParam(plt,freq,epsr,mur):
    plt.subplot(2,2,1)
    plt.plot(freq,np.real(epsr))
    plt.title("epsilon'")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("")
    plt.xlim([freq[0], freq[-1]])
    plt.grid(True)
    plt.subplot(2,2,3)
    plt.plot(freq,-1*np.imag(epsr))
    plt.title("epsilon''")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("")
    plt.xlim([freq[0], freq[-1]])
    plt.grid(True)
    plt.subplot(2,2,2)
    plt.plot(freq,np.real(mur))
    plt.title("mue'")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("")
    plt.xlim([freq[0], freq[-1]])
    plt.grid(True)
    plt.subplot(2,2,4)
    plt.plot(freq,-1*np.imag(mur))
    plt.title("mue''")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("")
    plt.xlim([freq[0], freq[-1]])
    plt.grid(True)
    return;

#Plots loss tangent and returns loss tangent vector
def plotLossTangent(plt,freq,epsr,conductivity):
    w = 2*np.pi*freq
    lossTangent = (w*-1*np.imag(epsr) + conductivity)/(w*np.real(epsr)) 
    
    plt.plot(freq,lossTangent)
    plt.title("Loss Tangent")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("")
    plt.xlim([freq[0], freq[-1]])
    plt.grid(True)
    return lossTangent;