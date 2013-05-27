'''
Created on May 12, 2013
Dielectric extraction algorithms
@author: Hardie Pienaar
'''
import numpy as np

#Extraction method for permitivity and permeability given S11 and S21
#This method has some restrictions in that it is ambiguous at half wavelength samples
#and also that an assumption is made in the B-field boundary conditions that no surface currents
#are generated. This is only true for low loss materials. Lastly multiple roots are generated
#when taking the log of e_gamma_d, this algorithm takes the fundamental root which might be wrong 
#for materials with high permittivity or permeabillity 
#TODO: Raise a warning if sample length is to long
def nicolsonRossWeir(freq,S11,S21,S12,S22,guideWidth,sampleLength):
    #Define free-space constants
    eps0 = 8.854187817e-12
    mu0 = 4*np.pi*1e-7
    c0 = 1/np.sqrt(eps0*mu0)
    
    #Calculate wave-guide cutoff frequency
    fc = 1/(2*guideWidth*np.sqrt(mu0*eps0)); 
    print("Fundamental mode cutoff: "+str(fc/1000000000)+"GHz")
    print("TE20 mode cutoff: "+str(2*fc/1000000000)+"GHz")
    #Convert to angular frequency
    wc = 2*np.pi*fc;
    w = 2*np.pi*freq;
    
    #Calculate correct gamma
    X = (np.square(S11) - np.square(S21) + 1)/(2*S11)

    gamma = []
    for XIter in X:
        g = XIter + np.sqrt(np.square(XIter) - 1)  
        #print(XIter)
        #print(g)
        if np.abs(g) > 1:
            gamma.append(XIter - np.sqrt(np.square(XIter) -1))
        else: 
            gamma.append(g)
    Gamma = np.array(gamma)
    
    
    
    #Calculate propagation term and refractive index respectfully    
    egamd = (S11 + S21 - Gamma)/(1 - (S11+S21)*Gamma);
    nsqr = -1.*np.square((c0/(w*sampleLength))*(np.log(1/egamd))) + np.square(wc/w);
    
    #Calculate relative dielectric parameters mu and eps
    mur = ((1+Gamma)/(1-Gamma))*np.sqrt((nsqr - np.square(wc/w))/(1 - np.square(wc/w)));
    epsr = nsqr/mur;
    return [epsr,mur];

#Extraction method for permitivity and permeability given S11 and S21
#This method has some restrictions in that it is ambiguous at half wavelength samples
#and also that an assumption is made in the B-field boundary conditions that no surface currents
#are generated. This is only true for low loss materials. Lastly multiple roots are handled.
#This extraction method makes sure the correct root is chosen in the log term by building the
#log term from a abs value and a unwrapped calculated phase term
#TODO: Raise a warning if sample length is to long
def nicolsonRossWeirWithStepWisePhase(freq,S11,S21,S12,S22,guideWidth,sampleLength):
    #Define free-space constants
    eps0 = 8.854187817e-12
    mu0 = 4*np.pi*1e-7
    c0 = 1/np.sqrt(eps0*mu0)
    
    #Calculate wave-guide cutoff frequency
    fc = 1/(2*guideWidth*np.sqrt(mu0*eps0)); 
    print("Fundamental mode cutoff: "+str(fc/1000000000)+"GHz")
    print("TE20 mode cutoff: "+str(2*fc/1000000000)+"GHz")

    #Convert to angular frequency
    wc = 2*np.pi*fc;
    w = 2*np.pi*freq;
    
    #Calculate correct gamma
    X = (np.square(S11) - np.square(S21) + 1)/(2*S11)
    gamma = []
    for XIter in X:
        g = XIter + np.sqrt(np.square(XIter) - 1)  
        if np.abs(g) > 1:
            gamma.append(XIter - np.sqrt(np.square(XIter) -1))
        else: 
            gamma.append(g)
    Gamma = np.array(gamma)
    
     
    z = np.sqrt((np.square(1 + S11) - np.square(S21))/(np.square(1-S11) - np.square(S21)));
    egamd = (1-np.square(S11)+np.square(S21))/(2*np.array(S21)) + (2*np.array(S11))/((z - 1/z)*S21);
   
    #Calculate unwrapped phase
    phi0 = float(np.angle(egamd[0],0))
    phi = [phi0]
    tmp_phi = 0
    for n in range(1,len(freq)):    
        for a in range(1,n):       
            tmp_phi = tmp_phi + float(np.angle(egamd[a]/egamd[a-1],0))        
        phi.append(phi0+tmp_phi)
        tmp_phi=0       
    phi = np.array(phi)
    
    
    #Calculate log of egamd with correct imaginary root
    lnegamd = np.log(abs(egamd)) + 1j*phi; 
    
    #Calculate propagation term and refractive index respectfully    
    nsqr = -1.*np.square((c0/(w*sampleLength))*(lnegamd)) + np.square(wc/w);
    
    #Calculate relative dielectric parameters mu and eps
    mur = ((1+Gamma)/(1-Gamma))*np.sqrt((nsqr - np.square(wc/w))/(1 - np.square(wc/w)));
    epsr = nsqr/mur;
    return [epsr,mur];

