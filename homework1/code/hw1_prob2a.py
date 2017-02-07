import sys
import numpy as np
import scipy
import pylab as plt
from pdb import set_trace as keyboard
from matplotlib import rc as matplotrc

PI=np.pi

matplotrc('text.latex',preamble='\usepackage{color}')
matplotrc('text',usetex=True)
matplotrc('font',family='serif')

figure_folder='../report/figures/'

x_eval=0.
#Staggered stencil

#   Centered
#   l=r=1
#Third order polynomiali
recon_order=3
l=2
r=2
n=l+r
a=np.linspace(0,-15,31)
dx=10**a
DX1 = 10**(np.linspace(0,-16,17))
DX2 = 10**(np.linspace(0,-8,17))
DX3 = 10**(np.linspace(0,-5.3,17))
DX4 = 10**(np.linspace(0,-4,17))
weights=np.zeros(n)
basefun=np.zeros(n-1)
deriv_order = 3
polyn_order = n-1
m=np.size(a)
polyder=np.zeros(len(a))
trun_err=np.zeros(len(a))
U=np.zeros((m,n))
f=np.zeros((m,n))
#staggered stencil
# print "U=",U
# print "f=",f
x0=0
x1=1
N=10
X=np.linspace(x0,x1,N)
Xst=(X[1:N]+X[0:N-1])/2
x=np.zeros(n)
dx=X[1]-X[0]
D3=np.zeros((N-1,N))
D1 = np.zeros((N-1,N))
print "length",len(X)
for ii in range((recon_order-1)/2,N-(recon_order+1)/2):
    x[0] = X[ii]-(l-1)*dx
    for j in range(1,n):
        x[j]=x[0]+j*dx
#    print "x=",x
    x_eval = Xst[ii]    
    for i,val in enumerate(x):
#        print i
#        print val
        basefun = np.zeros(n,)
        basefun[i] = 1.0
#       keyboard()
        p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
        weight3 = np.polyval(np.polyder(p_coef, 3), x_eval)
        weight1 = np.polyval(np.polyder(p_coef, 1), x_eval)
        D1[ii,ii-1+i] = weight1
        weights[i] = weight3
        D3[ii,ii-1+i]=weight3
#        polyder[ii] = polyder[ii] + f[ii,i]*weight
#        trun_err[ii]=((polyder[ii]-U[0,l])**2)**0.5


############# For the 1st row of D##########################
x=X[0:4]
for i,val in enumerate(x):
    basefun = np.zeros(n,)
    basefun[i] = 1.0
#       keyboard()
    p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
    weight3 = np.polyval(np.polyder(p_coef, 3), Xst[0])
    weight1 = np.polyval(np.polyder(p_coef, 1), Xst[0])

    weights[i] = weight3
    D3[0,i]=weight3
    D1[0,i] = weight1
#######################################################









############# For the last row of D##########################
x=X[N-4:N]
for i,val in enumerate(x):
    basefun = np.zeros(n,)
    basefun[i] = 1.0
#       keyboard()
    p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
    weight3 = np.polyval(np.polyder(p_coef, 3), Xst[len(Xst)-1])
    weight1 = np.polyval(np.polyder(p_coef, 1), Xst[len(Xst)-1])

    weights[i] = weight3
    D3[N-2,N-recon_order-1+i] = weight3
    D1[N-2 , N-recon_order-1+i] = weight1
#######################################################






#print weights
#print "approx derivative=",polyder
#print "actual derivative=",U[0,l]
figure_folder = "../report/figures/"
figure_name = "staggered_3rd_order_3rd_derivative.jpg"
figwidth=18
figheight=6
lineWidth=3
textFontSize=28
gcafontSize=10

fig = plt.figure(0, figsize=(figwidth,figheight))
fig1=fig.add_subplot(1,1,1)
fig1.spy(D3, markersize = 4)
figure_file_path= figure_folder + figure_name
print "saving figure"
plt.tight_layout()
plt.savefig(figure_file_path)
plt.close()




figure_name = "staggered_3rd_order_1st_derivative.jpg"
figwidth=18
figheight=6
lineWidth=3
textFontSize=28
gcafontSize=10

fig = plt.figure(0, figsize=(figwidth,figheight))
fig1=fig.add_subplot(1,1,1)
fig1.spy(D1, markersize = 4)
figure_file_path= figure_folder + figure_name
print "saving figure"
plt.tight_layout()
plt.savefig(figure_file_path)
plt.close()
#keyboard()
#print "D3=",D3
#print "D1=",D1











############5th order polynomial############

recon_order=5
l=3
r=3
n=l+r
a=np.linspace(0,-15,31)
dx=10**a
DX1 = 10**(np.linspace(0,-16,17))
DX2 = 10**(np.linspace(0,-8,17))
DX3 = 10**(np.linspace(0,-5.3,17))
DX4 = 10**(np.linspace(0,-4,17))
weights=np.zeros(n)
basefun=np.zeros(n-1)
deriv_order = 1
polyn_order = n-1
m=np.size(a)
polyder=np.zeros(len(a))
trun_err=np.zeros(len(a))
U=np.zeros((m,n))
f=np.zeros((m,n))
#staggered stencil
# print "U=",U
# print "f=",f
x0=0
x1=1
N=10
X=np.linspace(0,1,N)
Xst=(X[1:N]+X[0:N-1])/2
x=np.zeros(n)
dx=X[1]-X[0]
D3=np.zeros((N-1,N))
D1 = np.zeros((N-1,N))
print "length",len(X)
for ii in range((recon_order-1)/2,N-(recon_order+1)/2):
    x[0] = X[ii]-(l-1)*dx
    for j in range(1,n):
        x[j]=x[0]+j*dx
    print "x=",x
    x_eval = Xst[ii]    
    for i,val in enumerate(x):
#        print i
#        print val
        basefun = np.zeros(n,)
        basefun[i] = 1.0
#       keyboard()
        p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
        weight3 = np.polyval(np.polyder(p_coef, 3), x_eval)
        weight1 = np.polyval(np.polyder(p_coef, 1), x_eval)
        D1[ii,ii-2+i] = weight1
        weights[i] = weight3
        D3[ii,ii-2+i]=weight3
#        polyder[ii] = polyder[ii] + f[ii,i]*weight
#        trun_err[ii]=((polyder[ii]-U[0,l])**2)**0.5


############# For the 1st row of D##########################
x=X[0:recon_order+1]
for i,val in enumerate(x):
    basefun = np.zeros(n,)
    basefun[i] = 1.0
#       keyboard()
    p_coef=np.polyfit(x, basefun, recon_order)
#        print "i=",i
#        print "length=",len(weights)
#        print np.polyval(np.polyder(p_coef, deriv_order),x_eval)
    weight3 = np.polyval(np.polyder(p_coef, 3), Xst[0])
    weight1 = np.polyval(np.polyder(p_coef, 1), Xst[0])

    weights[i] = weight3
    D3[0,i]=weight3
    D1[0,i] = weight1

    weight3 = np.polyval(np.polyder(p_coef, 3), Xst[1])
    weight1 = np.polyval(np.polyder(p_coef, 1), Xst[1])
    D3[1,i]=weight3
    D1[1,i] = weight1


   
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
    weight3 = np.polyval(np.polyder(p_coef, 3), Xst[len(Xst)-1])
    weight1 = np.polyval(np.polyder(p_coef, 1), Xst[len(Xst)-1])

    weights[i] = weight3
    D3[N-2,N-recon_order-1+i] = weight3
    D1[N-2 , N-recon_order-1+i] = weight1


    weight3 = np.polyval(np.polyder(p_coef, 3), Xst[len(Xst)-2])
    weight1 = np.polyval(np.polyder(p_coef, 1), Xst[len(Xst)-2])

    weights[i] = weight3
    D3[N-3,N-recon_order-1+i] = weight3
    D1[N-3 , N-recon_order-1+i] = weight1
#
#######################################################






#print weights
#print "approx derivative=",polyder
#print "actual derivative=",U[0,l]
figure_folder = "../report/figures/"
figure_name = "staggered_5th_order_3rd_derivative.jpg"
figwidth=18
figheight=6
lineWidth=3
textFontSize=28
gcafontSize=10

fig = plt.figure(0, figsize=(figwidth,figheight))
fig1=fig.add_subplot(1,1,1)
fig1.spy(D3, markersize = 4)
figure_file_path= figure_folder + figure_name
print "saving figure"
plt.tight_layout()
plt.savefig(figure_file_path)
plt.close()




figure_name = "staggered_5th_order_1st_derivative.jpg"
figwidth=18
figheight=6
lineWidth=3
textFontSize=28
gcafontSize=10

fig = plt.figure(0, figsize=(figwidth,figheight))
fig1=fig.add_subplot(1,1,1)
fig1.spy(D1, markersize = 4)
figure_file_path= figure_folder + figure_name
print "saving figure"
plt.tight_layout()
plt.savefig(figure_file_path)
plt.close()
#keyboard()
#print "D3=",D3
#print "D1=",D1


