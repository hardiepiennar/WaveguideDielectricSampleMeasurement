'''
Created on May 12, 2013

Wave-guide parameter extraction code using different parameter 
extraction methods on data imported from S-parameter touch-stone
files.

@author: Hardie Pienaar
'''
import numpy as np 
import matplotlib.pyplot as plt
import readSParamFromFile as fileIO
import DielectricExtractionAlgorithms as dea
import parameterPlots as pp

if __name__ == '__main__':
    pass

#Show intro in command line
print("*******************************")
print("Wave guide parameter extraction")
print("*******************************")

#Define free-space constants
eps0 = 8.854187817e-12
mu0 = 4*np.pi*1e-7
c0 = 1/np.sqrt(eps0*mu0)

#Define input parameters
#filename = '../../Measurements/S-band 8 Mei 2013/PERSPEX_MUT_25mm_NTG'
filename = '../../Measurements/S-band 8 Mei 2013/FOAM_MUT_25mm_TG'
d = 25e-3
a = 72.136e-3
b = 34.036e-3


# filename = '../../Measurements/X-band 7 Mei 2013/MUT_5mm_TG'
# d = 5e-3
# a = 22.86e-3
# b = 10.16e-3


#Import S-Parameters
print("Importing S-Parameters")
freq, S11, S21, S12, S22 = fileIO.readSParamFromFile(filename)
print("")

#Plot S-Parameters
print("Plotting S-Parameters")
plt.figure(1)
pp.plotSParam(plt,freq,S11,S21,S12,S22)
print("")

#Run extraction algorithm to find epsilon and mu
print("Running Nicolson Ross Weir Extraction Method")
epsr, mur = dea.nicolsonRossWeir(freq,S11,S21,S12,S22,a,d)
print("Running Nicolson Ross Weir Extraction Method")
epsr2, mur2 = dea.nicolsonRossWeirWithStepWisePhase(freq,S11,S21,S12,S22,a,d)
print("")

#Calculate validity region
#TODO: Add start of valid region
#TODO: incorperate sample length
print("Calculating Valid Data Region")
endOfValiditySample = 0;
for i in range(1,len(freq)):
    if np.real(epsr2[i]) < 0 or -1*np.imag(epsr2[i]) < 0:
        endOfValiditySample = i;
        break;
    else:
        endOfValiditySample = len(freq) - 1
print("Valid region ends at "+ str(freq[endOfValiditySample]/1000000000) +"GHz")
print("")

#Do a straight line fit to the data
print("Fitting lines to material parameters")
order = 0
endOfValiditySample = np.floor(len(freq)/2)
coef = np.polyfit(freq[0:endOfValiditySample],epsr2[0:endOfValiditySample], order)
epsrFit= np.polyval(coef,freq)
coef = np.polyfit(freq[0:endOfValiditySample],mur2[0:endOfValiditySample], order)
murFit= np.polyval(coef,freq)
print("Fitting line with data up to " + str(freq[endOfValiditySample]/1000000000) +"GHz")
print("")

#Plot dielectric parameters with fitted data
print("Plotting Dielectric Parameters")
plt.figure(2)
pp.plotDielParam(plt,freq,epsr,mur)
pp.plotDielParam(plt,freq,epsr2,mur2)
pp.plotDielParam(plt,freq,epsrFit,murFit)
print("Plotting Loss Tangent")
plt.figure(3)
pp.plotLossTangent(plt, freq, epsr, 0) #TODO: Assume a conductivity of zero
pp.plotLossTangent(plt, freq, epsr2, 0) #TODO: Assume a conductivity of zero
lossTangent = pp.plotLossTangent(plt, freq, epsrFit, 0) #TODO: Assume a conductivity of zero
print("Epsilon Relative: " + str(epsrFit[0]))
print("Mue Relative    : " + str(murFit[0]))
print("Loss Tangent    : " + str(lossTangent[0]))

plt.show()
