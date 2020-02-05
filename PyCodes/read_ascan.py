#
#  Functions to manipulate A-scan waveform data
#
import numpy as np

def read_dat(fname):	# Read A-scan from csv file
	nhead=2
	fp=open(fname,'r')
	k=0
	time=[]
	amp=[]
	for raw in fp:
		item=raw.split(',')
		if k > 1: 
			time.append(float(item[0])*1.e06)
			amp.append(float(item[1]))
		k=k+1
#	print(fname)
	amp=amp-np.mean(amp)
	return time,amp

def ascan_fft(time,amp): # Apply FFT on A-scan
	Nt=len(time)
	dt=time[1]-time[0]
	fmax=1/dt	
	df=1/dt/Nt
	freq=np.linspace(0,fmax,Nt)
	Amp=np.fft.fft(amp)
#	Amp[Nt/2:]=0.0

	return freq,Amp


def butterworth(time,tb,sig):	# generate Butterworth filter 
	arg=(np.asarray(time)-tb)/sig
	Bw=1/(1+arg**4)

	return Bw

def butterworth(time,tb,sig,mexp):	# generate Butterworth filter 
	arg=(np.asarray(time)-tb)/sig
	Bw=1/(1+arg**mexp)

	return Bw

	
