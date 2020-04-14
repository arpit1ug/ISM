from scipy import *
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

global c,W,Y,F,L,num,start,end
c=3*10**5
L=1215.67*10**(-10)
W=2*np.pi*c*10**3/L
Y=6.265*10**8
F=0.42
num=1000
start=1212.0
end=1219.0

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


#------------------for different values of N-----------------
p=[]
for i in range(13,21):
	wave,p1=profile(12,i)
	p.append(p1)
	plt.plot(wave,p1,label="logN=%i"%i)

plt.xlabel("wavelength(in $\AA$)")
plt.ylabel("Normalized flux")
plt.legend()
plt.title(r"Lymen-$\alpha$ with b = 12km/s")
plt.xlim(1215,1216.5)
plt.show()

#------------------for different values of b-----------------
y=[]
p=[]
for j in range(1,6):
	wave,p1=profile(10*j,13.5)
	y=p1
	for i in range(10,len(y)-2):
		if (abs(y[i+1]-y[i])>5*abs(y[i]-y[i-1])):
			y[i+1]=0.5*y[i]+0.5*y[i+2]
	p.append(p1)
	e=10*j
	plt.plot(wave,p1,label="b=%ikm/s"%e)

plt.xlim(1215,1216.5)
plt.xlabel("wavelength(in $\AA$)")
plt.ylabel("Normalized flux")
plt.legend()
plt.title(r"Lymen-$\alpha$ with N = 13.5")
plt.show()

#-------------------COG--------------------------------------
Nstart=11.0
Nend=22.0
n=Nstart
bin=(Nend-Nstart)/100
WL=[]
N=[]
while(n<Nend):
	x,y=profile(12,n)
	WL.append(np.log10(np.sum(np.subtract(1,y))*(end-start)/(num*1215.67)))
	N.append(n+np.log10(1215.67*F)+4)
	#N.append(n)
	n+=bin

plt.plot(N,WL,label=r"COG for Lymen-$\alpha$")
plt.xlabel("Log[N$\lambda$f]")
plt.ylabel("Log[W/$\lambda$]")
plt.title(r"Curve of Growth for Lymen-$\alpha$ line with b=12km/s ")
m=1
c=-20-np.log10(1.132)-6.422
x=np.linspace(17,23,100)
y=m*x+c
plt.plot(x,y,':',label="linear region")
m=0.5
c=-18-np.log10(1.88)+1.27378
x=np.linspace(25,29,100)
y=m*x+c
plt.plot(x,y,'-.',label="damped region")
x=np.linspace(22,26,100)
y=6*np.log10(x/(1215.67*F))+4.2
plt.plot(x,y,'--',label="saturated region")
plt.legend()
plt.show()
