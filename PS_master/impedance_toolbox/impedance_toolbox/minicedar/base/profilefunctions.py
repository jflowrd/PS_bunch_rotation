'''
**Module containing profile functions to fit data.**
'''

from __future__ import division
import numpy as np


def leastSquareResidualFunction(fitParameters, *fittingArgList):
    '''
    * Function to be used for fitting in the minimize function (least square).*
    '''
    
    profileFitFunction = fittingArgList[0]
    profileInputX = fittingArgList[1]
    fittedProfileInputY = fittingArgList[2]
    
    residue = np.sqrt(np.sum((fittedProfileInputY - profileFitFunction(profileInputX, *fitParameters))**2))

    return residue
    
    
def gaussian(profileInputX, *fitParameters):
    '''
    * Gaussian line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    sigma = abs(fitParameters[2])

    lineDensityFunction = amplitude * np.exp(- (profileInputX - bunchPosition)**2 /(2*sigma**2))
    
    return lineDensityFunction


def generalizedGaussian(profileInputX, *fitParameters):
    '''
    * Generalized gaussian line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    alpha = abs(fitParameters[2])
    exponent = abs(fitParameters[3])

    lineDensityFunction = amplitude * np.exp(- (np.abs(profileInputX - bunchPosition) /alpha)**exponent)
    
    return lineDensityFunction


def waterbag(profileInputX, *fitParameters):
    '''
    * Waterbag distribution line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    bunchLength = abs(fitParameters[2])

    lineDensityFunction = np.zeros(len(profileInputX))
    lineDensityFunction[np.abs(profileInputX - bunchPosition) < bunchLength/2] = \
        amplitude * (1 - ((profileInputX[np.abs(profileInputX - bunchPosition) < bunchLength/2] - bunchPosition) / (bunchLength/2))**2)**0.5
    
    return lineDensityFunction
    

def parabolicLine(profileInputX, *fitParameters):
    '''
    * Parabolic line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    bunchLength = abs(fitParameters[2])

    lineDensityFunction = np.zeros(len(profileInputX))
    lineDensityFunction[np.abs(profileInputX - bunchPosition) < bunchLength/2] = \
        amplitude * (1 - ((profileInputX[np.abs(profileInputX - bunchPosition) < bunchLength/2] - bunchPosition) / (bunchLength/2))**2)
    
    return lineDensityFunction
    

def parabolicAmplitude(profileInputX, *fitParameters):
    '''
    * Parabolic in action line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    bunchLength = abs(fitParameters[2])

    lineDensityFunction = np.zeros(len(profileInputX))
    lineDensityFunction[np.abs(profileInputX - bunchPosition) < bunchLength/2] = \
        amplitude * (1 - ((profileInputX[np.abs(profileInputX - bunchPosition) < bunchLength/2] - bunchPosition) / (bunchLength/2))**2)**1.5
    
    return lineDensityFunction
    

def binomialAmplitude2(profileInputX, *fitParameters):
    '''
    * Binomial exponent 2 in action line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    bunchLength = abs(fitParameters[2])

    lineDensityFunction = np.zeros(len(profileInputX))
    lineDensityFunction[np.abs(profileInputX - bunchPosition) < bunchLength/2] = \
        amplitude * (1 - ((profileInputX[np.abs(profileInputX - bunchPosition) < bunchLength/2] - bunchPosition) / (bunchLength/2))**2)**2.5
    
    return lineDensityFunction
    

def binomialAmplitudeN(profileInputX, *fitParameters):
    '''
    * Binomial exponent n in action line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    bunchLength = abs(fitParameters[2])
    exponent = abs(fitParameters[3])
    
    lineDensityFunction = np.zeros(len(profileInputX))
    lineDensityFunction[np.abs(profileInputX - bunchPosition) < bunchLength/2] = \
        amplitude * (1 - ((profileInputX[np.abs(profileInputX - bunchPosition) < bunchLength/2] - bunchPosition) / (bunchLength/2))**2)**(exponent+0.5)
    
    return lineDensityFunction
    
    
def cosine(profileInputX, *fitParameters):
    '''
    * Cosine line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    bunchLength = abs(fitParameters[2])

    lineDensityFunction = np.zeros(len(profileInputX))
    lineDensityFunction[np.abs(profileInputX - bunchPosition) < bunchLength/2] = \
        amplitude * np.cos(np.pi * (profileInputX[np.abs(profileInputX - bunchPosition) < bunchLength/2] - bunchPosition) / bunchLength)
    
    return lineDensityFunction   
    

def cosineSquared(profileInputX, *fitParameters):
    '''
    * Cosine squared line density *
    '''
    
    amplitude = abs(fitParameters[0])
    bunchPosition = abs(fitParameters[1])
    bunchLength = abs(fitParameters[2])

    lineDensityFunction = np.zeros(len(profileInputX))
    lineDensityFunction[np.abs(profileInputX - bunchPosition) < bunchLength/2] = \
        amplitude * np.cos(np.pi * (profileInputX[np.abs(profileInputX - bunchPosition) < bunchLength/2] - bunchPosition) / bunchLength)**2.
    
    return lineDensityFunction 



