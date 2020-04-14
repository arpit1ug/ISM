import numpy as np
import matplotlib.pyplot as plt
from scipy import *
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import csv

global c,W,Y,F,L,num,start,end,wave,Nflux,data,atom
c=3*10**5
L=1144.93
W=2*np.pi*c*10**13/L
Y=3.52*10**8
F=0.083
num=1000
start=L-2
end=L+2

def voigtI(x,u,a):
	integrand=np.exp(-x**2)/((x-u)**2+a**2)
	return integrand


def voigt(l,b,N):
	w=2*np.pi*c*10**13/l
	u=(w-W)/W*c/b
	a=Y/(4*W)*c/b
	WD=W*b/c
	y=2*10**(N-4)*F/np.sqrt(np.pi)*a/WD*quad(voigtI,-np.inf,np.inf,args=(u,a))[0]
	return np.exp(-y)

def profile(b,N):
	l=start
	bin=(end-start)/num
	p=[]
	wave=[]
	while (l<end):
		p.append(voigt(l,b,N))
		wave.append(l)
		l+=bin
	return (wave,p)


#-----reading IGM data--------------------
f=open("hlsp_igm_hst_cos_1es1553_g130m-g160m_v3_spec.dat")
lines=f.readlines()
f.close()

rows, cols = (len(lines), 4) 
data = np.array([[0.]*cols]*rows)
#data 0 = wavelength
#data 1 = flux
#data 2 = error
#data 3 = continum
for l in range(0,len(lines)):
	data[l][0]=float(lines[l].split()[0])
	data[l][1]=float(lines[l].split()[1])
	data[l][2]=float(lines[l].split()[2])
	data[l][3]=float(lines[l].split()[3])

wave=data[:,0]


#------reading atomic lines data-------------
f=open("atom.dat")
lines=f.readlines()
f.close()
del lines[:6] # deleting first 6 lines
del lines[len(lines)-1]
WL=len(lines)
l=0
x1=[]
x2=[]
while (l<WL):
	if(lines[l].split()[0]=='!'): # deleting all ines starting with '!'
		del lines[l]
		WL=WL-1
	else:
		x=lines[l].split()[1]
		if(len(x)>3):
			y=float(lines[l].split()[1])
			if(y<1135.0 or y>1800.0): # deleting lines not required
				x1.append(lines[l])
				del lines[l]
				WL=WL-1
			else:
				l=l+1
		else:
			y=float(lines[l].split()[2])
			if(y<1135.0 or y>1800.0):
				x2.append(lines[l])
				del lines[l]
				WL=WL-1
			else:
				l=l+1

rows, cols = (len(lines), 3)
atom = np.array([[0.]*cols]*rows)
dataID=[]
#dataID = wavelength ID
#atom 0 = rest wavelength
#atom 1 = f
#atom 2 = gamma
for l in range(0,len(lines)):
	x=lines[l].split()[1]
	if(len(x)>3):
		dataID.append(lines[l].split()[0])
		atom[l][0]=float(lines[l].split()[1])
		atom[l][1]=float(lines[l].split()[2])
		atom[l][2]=float(lines[l].split()[3])
	else:
		dataID.append(lines[l].split()[0]+lines[l].split()[1])
		atom[l][0]=float(lines[l].split()[2])
		atom[l][1]=float(lines[l].split()[3])
		atom[l][2]=float(lines[l].split()[4])
#-------------------------------------------------------

def width(wl): # Q2A
	le=wave2index(wl)
	l=le
	m=le
	f1=0
	f2=0
	lo=0
	hi=0
	while(f1<0.95 or f2<0.95):
		if(f1<0.95):
			lo+=1
			l=l-1
		if(f2<0.95):
			hi+=1
			m=m+1
		f1=Nflux[l]
		f2=Nflux[m]
		#print lo,hi,f1,f2
	sum=0
	#print (wave[le-lo],wave[le+hi-1])
	for i in range(le-lo+1,le+hi):
		sum+=(1.0-Nflux[i])*0.03
	if (wl>1200 and wl<1202):
		#print wl
		return sum/2
	else: 
		return sum

def wave2index(l):
	i=max(min(int((l-1135.43)/0.03),20117),0)
	if(l<1135.43 or l>1799.95):
		return None
	while abs(wave[i]-l)>0.02:
		if(wave[i]<l):
			i=i+1
		else:
			i=i-1
	return i

#------------------correcting Ly-alpha-----------------------
i1=wave2index(1173)
i2=wave2index(1245)
l=i1
while(l<i2):
	data[l][3]=1.6*10**(-14)
	l=l+1

i1=wave2index(1213.93)
i2=wave2index(1217.21)
l=i1
while(l<i2):
	data[l][1]=0
	l=l+1

Nflux=np.divide(data[:,1],data[:,3])
#plt.plot(wave,Nflux)
#plt.show()

#------------------------------------------------------------


#-2 A-----------reading lines.txt and finding EW of all lines
f=open("lines.txt")
lines=f.readlines()
f.close()

name=[]
wave_rest=[]
EW=[]
for l in range(0,len(lines)):
	name.append(lines[l].split()[0])
	wave_rest.append(float(lines[l].split()[1]))
	EW.append(width(wave_rest[l]))

#zip(name,wave_rest,EW)
with open('line_output.csv', 'w') as f:
	writer = csv.writer(f, delimiter='\t')
	writer.writerows(zip(name,wave_rest,EW))

#----------------------------------------------------------------

#----2D----------finding Column density of Ni II and Fe II-------
name1=[]
wave_rest1=[]
for i in range(0,len(name)):
	if(name[i]=='FeII' or name[i]=='NiII'):
		name1.append(name[i])
		wave_rest1.append(wave_rest[i])

with open('line2c.csv', 'w') as f:
	writer = csv.writer(f, delimiter='\t')
	writer.writerows(zip(name1,wave_rest1))

		
#----------------------------------------------------------------











def search_atom(wl):
	x=atom[:,0]
	for i in range(0,len(x)):
		if(int(x[i])==wl):
			#print dataID[i]
			return (atom[i][0],atom[i][1],atom[i][2])
		else:
			i+=1

L,F,Y=search_atom(1142)
W=2*np.pi*c*10**13/L
start=L-2
end=L+2
num=200
b=20




i1=wave2index(L-0.5)
i2=wave2index(L+0.5)
NF1=Nflux[i1:i2]
wav1=np.add(0.12,wave)

N=15
NMax=0

for i in range(0,50):
	N=15+i/10.0
	print N
	wl,p1=profile(b,N)
	vmin=np.min(p1)
	dmin=np.min(NF1)
	print vmin,dmin
	print dmin-vmin
	if((dmin-vmin)<0.05 and (dmin-vmin)<0):
		i=i+1
	else:
		Nmax=N
		break

N=Nmax
b=30
N=23.3
wl,p1=profile(b,N)
plt.xlim(L-10,L+10)
plt.plot(wl,p1,"r",label="Voigt Profile b=30Km/s and logN=23.3")
plt.plot(wav1,Nflux,label="Lymen alpha data")
plt.legend()
plt.show()
