
import numpy as np
import matplotlib.pyplot as plt
import Ascan

class Bscan:
	def __init__(self,nums):
		self.nums=nums;
		dat=np.array([])
		for num in nums:
			fname="70_4_0211/scope_"+str(num)+".csv"
			#fname="70_0_0211/scope_"+str(num)+".csv"
			#fname="70_2_0212/scope_"+str(num)+".csv"
			print(fname)
			asc=Ascan.ASCAN(fname);
			dat=np.append(dat,np.array(asc.amp))
		self.Nx=len(nums);
		self.Nt=len(asc.time);

		dat=np.reshape(dat,(self.Nx,self.Nt))*2.0	# [m/sec]
		print(np.shape(dat));
		self.dat=dat;
		self.time=asc.time;
		self.Nt=len(asc.time);

	def show(self,ax,x1=0,dx=1):
		t1=self.time[0];
		t2=self.time[-1]; x2=x1+dx*(self.Nx-1)
		cax=ax.imshow(self.dat,extent=[t1,t2,x1,x2],aspect="auto",origin="lower",cmap="gray",vmin=-0.03,vmax=0.03);
		ax.set_xlabel("time [micro sec]",fontsize=12);
		ax.set_ylabel("y ",fontsize=12);

		return cax

class Ascans:
	def __init__(self,nums):
		self.nums=nums;
		self.nwv=len(nums)
		awvs=[];
		for num in nums:
			fname="scope_"+str(num)+".csv"
			print(fname)
			asc=Ascan.ASCAN(fname);
			asc.amp*=2.0;
			awvs.append(asc)

		self.awvs=awvs;
		self.Nt=len(asc.amp);

	def showB(self,ax):
		dat=np.array([])
		f0=0.2;
		sig=1/f0;
		tb=15.0
		for k in range(self.nwv):
			#dat=np.append(dat,self.awvs[k].amp)			
			self.awvs[k].butterworth(tb,sig,apply=True)
			self.awvs[k].gdelay();
			#dat=np.append(dat,self.awvs[k].amp)
			dat=np.append(dat,self.awvs[k].tg)

		#dat=np.reshape(dat,(self.nwv,self.Nt))
		dat=np.reshape(dat,(self.nwv,self.Nt-1))
		wv=self.awvs[0];
		t1=wv.time[0];
		t2=wv.time[-1];
		th1=0;
		th2=360;
		cax=ax.imshow(dat,extent=[t1,t2,th1,th2],aspect="auto",vmin=0.,vmax=20.);
	
		return cax	


if __name__=="__main__":

	nums=range(3,404,5); # file numbers

	fig=plt.figure()
	ax=fig.add_subplot(111)

	Bscn=Bscan(nums);   # read B-scan data
	#ascs=Ascans(nums);
	cax=Bscn.show(ax,x1=140,dx=0.25) # show B-scan image
	#cax=ascs.showB(ax)
	#cbar=fig.colorbar(cax,orientation="vertical"); # set colorbar
	ax.set_xlim([16,66]) # set x limit
	fig.savefig("bscan.png",bbox_inches="tight") # save image 
	plt.show()
