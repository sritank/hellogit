import sys
import numpy as np
import scipy
import pylab as plt
import scipy.sparse.linalg as splinalg
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from pdb import set_trace as keyboard
from matplotlib import rc as matplotrc

PI=np.pi

matplotrc('text.latex',preamble='\usepackage{color}')
matplotrc('text',usetex=True)
matplotrc('font',family='serif')

figure_folder='../report'

x_eval=0.
#Collocated stencil

#   Centered
#   l=r=1
recon_order=3
l=2
r=2
n=l+r+1
weights=np.zeros(n)
basefun=np.zeros(n-1)
deriv_order = 3 
polyn_order = n-1
x0=0
x1=PI
a=np.linspace(np.log(6)/np.log(10),np.log(25000)/np.log(10),15)
Nmatrix=10**a
print "Nmatrix=",Nmatrix
Nmatrix=np.ceil(Nmatrix)
DX1 = 10**(np.linspace(0,-8,17))
DX2 = 10**(np.linspace(0,-4,17))
DX3 = 10**(np.linspace(0,-8/3,17))
DX4 = 10**(np.linspace(0,-2,17))
dx1=np.zeros(len(a))
eps=np.zeros(len(a))
for iN,nn in enumerate(Nmatrix):
    N=nn.astype(int)
    L=np.zeros((N,N))
    R=L
    X = np.linspace(x0,x1,N)
    f = np.sin(X)
    dx=X[1]-X[0]
    dx1[iN]=dx
    print "dx=",dx
    D=np.zeros((N,N))
    print "length N = ",N
    for ii in range(2,N-2):
        L[ii,ii-2] = -2/(dx**3)
        L[ii,ii-1] = 4/(dx**3)
        L[ii,ii+1] = -4/(dx**3)
        L[ii,ii+2] = 2/(dx**3)
    L[N-2,N-4] = -2/dx**3
    L[N-2,N-3] = 6/dx**3
    L[N-2,N-2] = -6/dx**3
    L[N-2,N-1] = 2/dx**3
    L[0,0] = 1
    L[1,0] = -3/(2*dx)
    L[1,1] = 2./(dx)
    L[1,2] = -1/(2*dx)

    L[N-1,N-3] = 1/(2*dx)
    L[N-1,N-2] = -2/(dx)
    L[N-1,N-1] = 3/(2*dx)
   #######################################################
#    print "L=",L
#    keyboard()

    y=np.zeros((N,1))
    y3=-1.0*np.cos(X)       #Third derivative of sin(x)
    for i in range(1,N-2):
        y[1+i,0] = y3[i-1]+2*y3[i] + y3[i+1]
    y[1,0] = np.cos(X[0])
    y[0,0] = np.sin(X[0])
    y[N-1,0] = np.cos(X[N-1])
#    print "y=",y
    ysol=splinalg.spsolve(L,y)
    Y = np.transpose(ysol)
#    print "Y=",Y
    eps[iN]=np.power(np.sum(np.power((Y-f),2))/N,0.5)

figure_folder = "../report/figures/"
figure_name = "Pade.jpg"
figure_name2 = "padesolutions.jpg"
figwidth=18
figheight=9
lineWidth=3
textFontSize=28
gcafontSize=12
#keyboard()
#print "D=",D

fig = plt.figure(0, figsize = (figwidth,figheight))
fig1 = fig.add_subplot(1,1,1)
plt.axes(fig1)
fig1.set_xlabel(r"$\frac{1}{\delta x}$ axis",fontsize=textFontSize)
fig1.set_ylabel(r"truncation error($\epsilon$)",fontsize=textFontSize,rotation=90)  
#print "1/dx=",1/dx1
#fig1.plot(X,ysol,'-b',linewidth = lineWidth)
#fig1.plot(X,np.sin(X),'-r',linewidth = lineWidth)
es,=fig1.loglog(1/dx1,eps,'-b',linewidth=lineWidth,label = "RMS error")
x1,=fig1.loglog(1/DX1,DX1,'--r',linewidth=1,label = "dx^1")
x2,=fig1.loglog(1/DX2,DX2**2,'--r',linewidth=2,label = "dx^2")
x3,=fig1.loglog(1/DX3,DX3**3,'--r',linewidth=3,label = "dx^3")
x4,=fig1.loglog(1/DX4,DX4**4,'--r',linewidth=4,label = "dx^4")
          
handles, labels = fig1.get_legend_handles_labels()   
fig1.legend(loc='upper right')
figure_filepath = figure_folder + figure_name
print "saving figure"
plt.tight_layout()
plt.savefig(figure_filepath)
plt.close()


#plotting the solution
fig = plt.figure(0, figsize = (figwidth,figheight))
fig1 = fig.add_subplot(1,1,1)
plt.axes(fig1)
fig1.plot(X,ysol,'-b',linewidth = lineWidth)
fig1.plot(X,np.sin(X),'-r',linewidth = lineWidth)
fig1.grid('on',which = 'both')
fig1.set_xlabel(r"X", fontsize = textFontSize)
fig1.set_ylabel(r"f", fontsize = textFontSize)

figure_filepath = figure_folder + figure_name2
print "saving figure"
plt.tight_layout()
plt.savefig(figure_filepath)
plt.close()


