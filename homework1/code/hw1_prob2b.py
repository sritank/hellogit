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
#Staggered stencil

#   Centered
#   l=r=2
#Third order polynomiali
recon_order=3
l=2
r=2
n=l+r
weights=np.zeros(n)
basefun=np.zeros(n-1)
deriv_order = 3 
polyn_order = n-1
x0=0
x1=PI
a=np.linspace(np.log(6)/np.log(10),np.log(15000)/np.log(10),11)
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
    X = np.linspace(x0,x1,N)
    f = np.sin(X)
    Xst=(X[1:N]+X[0:N-1])/2.
    x=np.zeros(n)
    dx=X[1]-X[0]
    dx1[iN]=dx
    print "dx=",dx
    D = lil_matrix((N,N))
    print "length N = ",N
    for ii in range(((recon_order-1)/2),(N-(recon_order+1)/2)):
        x[0] = X[ii]-(l-1)*dx
        for j in range(1,n):
            x[j]=x[0]+j*dx
        x_eval = Xst[ii]    
        for i,val in enumerate(x):
           basefun = np.zeros(n,)
           basefun[i] = 1.0
           p_coef=np.polyfit(x, basefun, recon_order)
           weight3 = np.polyval(np.polyder(p_coef, 3), x_eval)
           weights[i] = weight3
           D[ii+1,ii-1+i]=weight3
############# For the 2nd row of D##########################
    x=X[0:recon_order+1]
    for i,val in enumerate(x):
        basefun = np.zeros(n,)
        basefun[i] = 1.0
#       keyboard()
        p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
        weight1 = np.polyval(np.polyder(p_coef, 1), Xst[0])
        weights[i] = weight1
        D[1,i]=weight1   
    D[0,0] = 1    
#######################################################


############# For the last row of D##########################
    x=X[N-recon_order-1:N]
    for i,val in enumerate(x):
        basefun = np.zeros(n,)
        basefun[i] = 1.0
#       keyboard()
        p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
        weight1 = np.polyval(np.polyder(p_coef, 1), Xst[len(Xst)-1])

        weights[i] = weight1
        D[N-1,N-recon_order-1+i] = weight1
#######################################################
#    keyboard()
#    print "D=",D
#    keyboard()
    D=D.tocsr()
    y=np.zeros((N,1))
    y1=-1.0*np.cos(Xst)
    for i in range(1,N-2):
        y[1+i,0]=-1.*np.cos(Xst[i])
    y[1,0] = np.cos(Xst[0])
    y[0,0] = np.sin(X[0])
    y[N-1,0] = np.cos(Xst[N-2])
#    print "y=",y
    ysol=splinalg.spsolve(D,y)
    Y = np.transpose(ysol)
#    print "Y=",Y
    eps[iN]=np.power(np.sum(np.power((Y-f),2))/N,0.5)
#    print "eps=",eps

print weights
figure_folder = "../report/figures/"
figure_name = "staggered_prob2b_3rdorder.jpg"
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
fig1.grid('on',which = 'both')
fig1.set_xlabel(r"$\delta x^{-1}", fontsize = textFontSize)
fig1.set_ylabel(r"$\epsilon$", fontsize = textFontSize)

#print "1/dx=",1/dx1
#fig1.plot(X,ysol,'-b',linewidth = lineWidth)
#fig1.plot(X,np.sin(X),'-r',linewidth = lineWidth)
fig1.loglog(1/dx1,eps,'-b',linewidth=lineWidth,label="RMS error")
x1,=fig1.loglog(1/DX1,DX1,'--r',linewidth=1,label = "dx^1")
x2,=fig1.loglog(1/DX2,DX2**2,'--r',linewidth=2,label = "dx^2")
x3,=fig1.loglog(1/DX3,DX3**3,'--r',linewidth=3,label = "dx^3")
x4,=fig1.loglog(1/DX4,DX4**4,'--r',linewidth=4,label = "dx^4")
fig1.legend(loc='upper right')

figure_filepath = figure_folder + figure_name
print "saving figure"
plt.tight_layout()
plt.savefig(figure_filepath)
plt.close()








##########################5th order polynomial#################################

recon_order=5
l=3
r=3
n=l+r
weights=np.zeros(n)
basefun=np.zeros(n-1)
deriv_order = 3 
polyn_order = n-1
x0=0
x1=PI
a=np.linspace(np.log(10)/np.log(10),np.log(10000)/np.log(10),10)
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
    X = np.linspace(x0,x1,N)
    f = np.sin(X)
    Xst=(X[1:N]+X[0:N-1])/2.
    x=np.zeros(n)
    dx=X[1]-X[0]
    dx1[iN]=dx
    print "dx=",dx
    D = np.zeros((N,N))
    print "length N = ",N
    for ii in range(((recon_order-1)/2),(N-(recon_order+1)/2)):
        x[0] = X[ii]-(l-1)*dx
        for j in range(1,n):
            x[j]=x[0]+j*dx
        x_eval = Xst[ii]    
        for i,val in enumerate(x):
           basefun = np.zeros(n,)
           basefun[i] = 1.0
           p_coef=np.polyfit(x, basefun, recon_order)
           weight3 = np.polyval(np.polyder(p_coef, 3), x_eval)
           weights[i] = weight3
           D[ii+1,ii-2+i]=weight3
############# For the 2nd and 3rd row of D##########################
    x=X[0:recon_order+1]
    for i,val in enumerate(x):
        basefun = np.zeros(n,)
        basefun[i] = 1.0
#       keyboard()
        p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
        weight1 = np.polyval(np.polyder(p_coef, 1), Xst[0])
        weights[i] = weight1
        D[1,i]=weight1   
        weight3 = np.polyval(np.polyder(p_coef, 3), Xst[1])
        D[2,i] = weight3

    D[0,0] = 1    
#######################################################


############# For the second last and last row of D##########################
    x=X[N-recon_order-1:N]
    for i,val in enumerate(x):
        basefun = np.zeros(n,)
        basefun[i] = 1.0
#       keyboard()
        p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
        weight1 = np.polyval(np.polyder(p_coef, 1), Xst[len(Xst)-1])

        weight3 = np.polyval(np.polyder(p_coef, 3), Xst[len(Xst)-2])
        D[N-2,N-recon_order-1+i] = weight3
        D[N-1,N-recon_order-1+i] = weight1
        weights[i] = weight1
        
#######################################################
#    keyboard()
#    print "D=",D
#    keyboard()
#    D=D.tocsr()
    y=np.zeros((N,1))
    y1=-1.0*np.cos(Xst)
    for i in range(1,N-2):
        y[1+i,0]=-1.*np.cos(Xst[i])
    y[1,0] = np.cos(Xst[0])
    y[0,0] = np.sin(X[0])
    y[N-1,0] = np.cos(Xst[N-2])
#    print "y=",y
    ysol=splinalg.spsolve(D,y)
    Y = np.transpose(ysol)
#    print "Y=",Y
    eps[iN]=np.power(np.sum(np.power((Y-f),2))/N,0.5)
#    print "eps=",eps

#print weights
figure_folder = "../report/figures/"
figure_name = "staggered_prob2b_5thorder.jpg"
figwidth=18
figheight=9
lineWidth=3
textFontSize=28
gcafontSize=12
#keyboard()
#print "D=",D
#print D[N-2,N-15:N-1]
fig = plt.figure(0, figsize = (figwidth,figheight))
fig1 = fig.add_subplot(1,1,1)
plt.axes(fig1)
fig1.grid('on',which = 'both')
fig1.set_xlabel(r"$\delta x^{-1}", fontsize = textFontSize)
fig1.set_ylabel(r"$\epsilon$", fontsize = textFontSize)
#print "1/dx=",1/dx1
#fig1.plot(X,ysol,'-b',linewidth = lineWidth)
#fig1.plot(X,np.sin(X),'-r',linewidth = lineWidth)
es,= fig1.loglog(1/dx1,eps,'-b',linewidth=lineWidth,label="RMS error")
x1,=fig1.loglog(1/DX1,DX1,'--r',linewidth=1,label = "dx^1")
x2,=fig1.loglog(1/DX2,DX2**2,'--r',linewidth=2,label = "dx^2")
x3,=fig1.loglog(1/DX3,DX3**3,'--r',linewidth=3,label = "dx^3")
x4,=fig1.loglog(1/DX4,DX4**4,'--r',linewidth=4,label = "dx^4")
fig1.legend(loc='upper right')
figure_filepath = figure_folder + figure_name
print "saving figure"
plt.tight_layout()
plt.savefig(figure_filepath)
plt.close()







