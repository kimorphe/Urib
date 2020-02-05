#!/home/kazushi/anaconda3/bin/python
#  A-scan waveform class 
"""
	A-scan (time-series) signal and basic signal processing methods  
"""

import numpy as np
import matplotlib.pyplot as plt

class ASCAN:
    # Read A-scan from csv file
    def __init__(self,fname):
        nhead=2
        #print(fname)
        fp=open(fname,'r')
        k=0
        time=[]
        amp=[]
        for raw in fp:
            item=raw.strip().split(',')
            if k >= 0: 
                time.append(float(item[0])*1.e06)
                #print(k,item[1])
                amp.append(float(item[1]))
            k=k+1
        amp=amp-np.mean(amp)

        self.time=np.array(time);
        self.amp=np.array(amp);
        self.dt=time[2]-time[1];
        self.Nt=len(time);
        self.fname=fname;

    def fft(self): # Apply FFT on A-scan
        dt=self.dt;
        fmax=1/self.dt
        self.df=1/self.dt/self.Nt
        self.freq=np.linspace(0,fmax,self.Nt)
        self.Amp=np.fft.fft(self.amp)
#       Amp[Nt/2:]=0.0

    def butterworth(self,tb,sig,apply=True):	# generate Butterworth filter
        arg=(np.asarray(self.time)-tb)/sig
        Bw=1/(1+arg**4)
        self.Bw=Bw;
        if apply==True:
            self.amp*=Bw;
        return Bw
    def gdelay(self):
        PI=np.pi; self.fft();
        phi=np.unwrap(np.angle(self.Amp))	# get phase angle
        self.tg=-np.diff(phi)/(2*self.df*PI) #	Group Delay
        self.freq1=self.freq[1:]	# frequency axis data for group delay
        # tgb=np.mean(tg[nf1:nf2])	# mean group delay
    def show_amp(self,ax):
        ax.plot(self.time,self.amp);
        ax.grid(True);
        #ax.set_title(self.fname)
        #ax.set_xlabel("time (micro sec)");
        ax.set_ylabel("amplitude (volt)");

    def show_Amp(self,ax):
        #ax.plot(self.freq,np.abs(self.Amp));
        Amax=np.max(np.abs(self.Amp))
        ax.semilogy(self.freq,np.abs(self.Amp)/Amax);
        ax.grid(True);
        ax.set_xlim([0,10.0]);
        ax.set_ylim([0.001,1]);
        #ax.set_title(self.fname)
        #ax.set_xlabel("frequency  (MHz)");
        ax.set_ylabel("Fourier amplitude");
"""
	def show_tg(self,ax):
		ax.plot(self.freq1,self.tg);
		ax.grid(True);
		ax.set_xlim([0,0.6]);
		ax.set_title(self.fname)
		ax.set_xlabel("frequency  (MHz)");
		ax.set_ylabel("group delay (micro sec)");

	def export(self,fname):
		fp=open(fname,"w");
		item="frequency [MHz],group delay [micro sec]\n"
		fp.write(item)
		Nt2=int(self.Nt*0.5)
		for k in range(Nt2):
			data_string="%lf,%lf\n" % (self.freq1[k],self.tg[k]);
			fp.write(data_string);
		fp.close()
	def interp(self, time):
		ts=self.time[0];
		te=self.time[-1];	
		amp=[];
		for tt in time:
			tm=tt-ts;
			if tm < 0.0:
				amp.append(0.0);				
			elif tm > te:
				amp.append(0.0);
			else:
				i1=int(np.floor(tm/self.dt));
				#print i1	
				i2=i1+1;
				xi=tm-self.dt*i1;	
				eta=1.-xi;
				amp.append( self.amp[i1]*eta+self.amp[i2]*xi);
		return amp

"""	
	
if __name__ == "__main__":

	fname="scope_10.csv"
	f0=0.2; tb=10.0;
	sig=1./f0;
	fname=input("Enter file name (e.g. scope_10.csv)")
	tb=float(input(" Enter approximate flight-time [micro sec] "))


	fig1=plt.figure(); ax1=fig1.add_subplot(111)
	fig2=plt.figure(); ax2=fig2.add_subplot(111)
	fig3=plt.figure(); ax3=fig3.add_subplot(111)

	wv1=ASCAN(fname);	# create A-scan instance by loading data
	wv1.fft();
	wv1.show_amp(ax1);	
	wv1.show_Amp(ax2);	
	plt.show()
	wv1.butterworth(tb,sig,apply=True);	# apply butterworth window
	wv1.fft();
	wv1.show_amp(ax1);	
	wv1.show_Amp(ax2);	
	wv1.gdelay();	# get group delay

	wv1.show_tg(ax3);

	wv1.export("group_delay.csv");	
	plt.show();
