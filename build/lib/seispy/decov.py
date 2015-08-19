# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 14:21:38 2014
[RFI, RMS, it]=makeRFitdecon(UIN,WIN,DT,NT,TSHIFT,F0,ITMAX,MINDERR)

In:
UIN = numerator (radial for PdS)
WIN = denominator (vertical component for PdS)
DT = sample interval (s)
NT = number of samples
TSHIFT = Time until beginning of receiver function (s)
F0 = width of gaussian filter
ITMAX = max # iterations
MINDERR = Min change in error required for stopping iterations

Out:
RFI = receiver function
RMS = Root mean square error for predicting numerator after each iteration

@author: xumj
"""

import numpy as np
from obspy.signal.util import nextpow2
from math import pi

def gaussFilter( dt, nft, f0 ):
	df = 1.0/(nft*dt)
	nft21 = 0.5*nft + 1
	f = df*np.arange(0,nft21)
	w = 2*pi*f
	
	gauss = np.zeros([nft,1])
	gauss1 = np.exp(-0.25*(w/f0)**2)/dt
	gauss1.shape=(len(gauss1),1)
	gauss[0:int(nft21)] =gauss1
	gauss[int(nft21):] = np.flipud(gauss[1:int(nft21)-1])
	gauss = gauss[:,0]
	
	return gauss

def gfilter(x, nfft, gauss, dt):
	Xf = np.fft.fft(x, nfft)
	Xf = Xf*gauss*dt
	xnew =  np.fft.ifft(Xf, nfft).real	
	return xnew

def correl( R, W, nfft ):
	x = np.fft.ifft(np.fft.fft(R,nfft)*np.conj(np.fft.fft(W,nfft)), nfft)
	x = x.real	
	return x

def phaseshift( x, nfft, DT, TSHIFT ):
	Xf = np.fft.fft(x, nfft)	
	shift_i = int(TSHIFT/DT)
	p = 2*pi*np.arange(1,nfft+1)*shift_i/(nfft);
	Xf = Xf*np.vectorize(complex)(np.cos(p),-np.sin(p))	
	x = np.fft.ifft(Xf, nfft)/np.cos(2*pi*shift_i/nfft)
	x = x.real	
	return x
	
	
def decovit(UIN,WIN,DT,NT,TSHIFT,F0,ITMAX,MINDERR):	
	print('Iterative Decon (Ligorria & Ammon):\n')  	
	
	RMS = np.zeros([ITMAX,1]) 
	nfft = nextpow2(NT) 
	P0 = np.zeros(nfft) 

	U0 = np.zeros(nfft)
	W0 = np.zeros(nfft)
	
	U0[0:NT]=UIN
	W0[0:NT]=WIN
	
	gaussF = gaussFilter( DT, nfft, F0 )
	
	U = gfilter( U0, nfft, gaussF , DT)
	W = gfilter( W0, nfft, gaussF , DT)
	
	Wf = np.fft.fft( W0, nfft )
	R = U
	
	powerU = np.sum(U**2)
	
	it = 0
	sumsq_i = 1
	d_error = 100*powerU + MINDERR
	maxlag = 0.5*nfft
	print('\tMax Spike Display is '+str((maxlag)*DT))
	
	while np.abs(d_error) > MINDERR  and it < ITMAX:
		
		
		RW= correl(R, W, nfft)
		RW = RW/np.sum(W**2)
		
		i1=np.argmax( np.abs( RW[0:int(maxlag)-1] ) )
		amp = RW[i1]/DT
		
		P0[i1] = P0[i1] + amp
		P = gfilter(P0, nfft, gaussF, DT)
		P = gfilter(P, nfft, Wf, DT)
		
		R = U - P
		sumsq = np.sum( R**2 )/powerU
		RMS[it] = sumsq
		d_error = 100*(sumsq_i - sumsq)
		
		sumsq_i = sumsq
		
		it = it+1
		
	P = gfilter( P0, nfft, gaussF, DT )
	P = phaseshift(P,nfft,DT,TSHIFT)
	RFI=P[0:NT]
	RMS = RMS[0:it-1]
	
	print('\t# iterations: ',it,'\n')
	print('\tFinal RMS: ',float(RMS[it-2]),'\n')
	
	return RFI, RMS, it
	
	
	
	
	
	
	
