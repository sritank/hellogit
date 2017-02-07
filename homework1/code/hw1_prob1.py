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
der = (1/(np.cosh(x_eval))**2)*np.sin(5*x_eval+1.5) + 5*np.tanh(x_eval)*np.cos(5*x_eval+1.5)
###############################Collocated vs Staggered############################
#Collocated stencil

LL=[1,2,3,0,0,0]
LL1=[1,2,3,1,1,1]
RR=[1,2,3,1,2,3]
for l,r,ll in zip(LL,RR,LL1):

    n=l+r+1
    a=np.linspace(0,-15,1501)
    dx=10**a
    DX1 = 10**(np.linspace(0,-16,17))
    DX2 = 10**(np.linspace(0,-8,17))
    DX3 = 10**(np.linspace(0,-5.3,17))
    DX4 = 10**(np.linspace(0,-4,17))
    DX5 = 10**(np.linspace(0,-3.2,17))

    weights=np.zeros(n)
    basefun=np.zeros(n-1)
    deriv_order = 1
    polyn_order = n-1
    m=np.size(a)
    polyder=np.zeros(len(a))
    trun_err=np.zeros(len(a))
    U=np.zeros((m,n))
    f=np.zeros((m,n))
    # print "U=",U
    # print "f=",f
    x=np.zeros(n)
    for ii,dxi in enumerate(dx):
        x[0] = -l*dxi
        for j in range(1,n):
            xx = -(l-j)*dxi
    #        print "x=",x,"xx=",xx
            x[j]=xx
        for i,xi in enumerate(x):
            U[ii,i] = (1/(np.cosh(xi))**2)*np.sin(5*xi+1.5) + 5*np.tanh(xi)*np.cos(5*xi+1.5)
            f[ii,i] = np.tanh(xi)*np.sin(xi+1.5)

        for i,val in enumerate(x):

            basefun = np.zeros(n,)
            basefun[i] = 1.0
            p_coef=np.polyfit(x, basefun, polyn_order)

            weight = np.polyval(np.polyder(p_coef, deriv_order), x_eval)
            weights[i] = weight
            polyder[ii] = polyder[ii] + f[ii,i]*weight
            trun_err[ii]=((polyder[ii]-der )**2)**0.5         #Same as taking the absolute value
                








#Staggered stencil
    n=ll+r
    a=np.linspace(0,-20,2001)
    dx1=10**a
    weights=np.zeros(n)
    basefun=np.zeros(n-1)
    deriv_order = 1
    polyn_order = n-1
    m=np.size(a)
    polyder=np.zeros(len(a))
    trun_err_stag=np.zeros(len(a))
    Us=np.zeros((m,n))
    fs=np.zeros((m,n))
    x=np.zeros(n)
    for ii,dxi in enumerate(dx1):
        x[0] = -(ll-1)*dxi-dxi/2
        for j in range(1,n):
            xx = x[0]+j*dxi
            x[j]=xx
        for i,xi in enumerate(x):
            Us[ii,i] = (1/(np.cosh(xi))**2)*np.sin(5*xi+1.5) + 5*np.tanh(xi)*np.cos(5*xi+1.5)
            fs[ii,i] = np.tanh(xi)*np.sin(xi+1.5)

        for i,val in enumerate(x):
            basefun = np.zeros(n,)
            basefun[i] = 1.0
            p_coef=np.polyfit(x, basefun, polyn_order)

            weight = np.polyval(np.polyder(p_coef, deriv_order), x_eval)
            weights[i] = weight
            polyder[ii] = polyder[ii] + fs[ii,i]*weight
            trun_err_stag[ii]=((polyder[ii]-der )**2)**0.5     #Same as taking the absolute value
            
    figure_folder = "../report/figures/"
    figure_name = "staggered"+"l"+str(ll)+"r"+str(r)+"_collocated"+"l"+str(l)+"r"+str(r)+".jpg"
    figwidth=18
    figheight=12
    lineWidth=4
    textFontSize=30
    gcafontSize=60

    fig = plt.figure(0, figsize=(figwidth,figheight))
    fig1 = fig.add_subplot(1,1,1)        
    plt.axes(fig1)
    fig1.grid('on',which='both')
    fig1.set_xlabel(r"$\frac{1}{\delta x}$ axis",fontsize=textFontSize)
    fig1.set_ylabel(r"truncation error($\epsilon$)",fontsize=textFontSize,rotation=90)        
    ec,=fig1.loglog(1/dx,trun_err,'-g',linewidth=lineWidth,label = "eps collocated")    
    es,=fig1.loglog(1/dx1,trun_err_stag,'-b',linewidth=lineWidth, label = "eps staggered")
    x1,=fig1.loglog(1/DX1,DX1,'--r',linewidth=1,label = "dx^1")
    x2,=fig1.loglog(1/DX2,DX2**2,'--r',linewidth=2,label = "dx^2")
    x3,=fig1.loglog(1/DX3,DX3**3,'--r',linewidth=3,label = "dx^3")
    x4,=fig1.loglog(1/DX4,DX4**4,'--r',linewidth=4,label = "dx^4")
    x5,=fig1.loglog(1/DX5,DX5**5,'--r',linewidth=4,label = "dx^5")
    fig1.set_xlabel(r"$\frac{1}{\delta x}$ axis",fontsize=textFontSize)
    fig1.set_ylabel(r"truncation error($\epsilon$)",fontsize=textFontSize,rotation=90)            
    handles, labels = fig1.get_legend_handles_labels()   
    fig1.legend(loc='upper right')




    figure_file_path= figure_folder + figure_name
    print "saving figure"
    plt.tight_layout()
    plt.savefig(figure_file_path)
    plt.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#########################Centered vs Biased###############################    
    
    
    
#Collocated centered

LL=[1,2,3]
LL1=[0,0,0]
RR=[1,2,3]
for l,r,ll in zip(LL,RR,LL1):

    n=l+r+1
    a=np.linspace(0,-15,1501)
    dx=10**a
    weights=np.zeros(n)
    basefun=np.zeros(n-1)
    deriv_order = 1
    polyn_order = n-1
    m=np.size(a)
    polyder=np.zeros(len(a))
    trun_err_cent=np.zeros(len(a))
    U=np.zeros((m,n))
    f=np.zeros((m,n))
    # print "U=",U
    # print "f=",f
    x=np.zeros(n)
    for ii,dxi in enumerate(dx):
        x[0] = -l*dxi
        for j in range(1,n):
            xx = -(l-j)*dxi
    #        print "x=",x,"xx=",xx
            x[j]=xx
        for i,xi in enumerate(x):
            U[ii,i] = (1/(np.cosh(xi))**2)*np.sin(5*xi+1.5) + 5*np.tanh(xi)*np.cos(5*xi+1.5)
            f[ii,i] = np.tanh(xi)*np.sin(xi+1.5)

        for i,val in enumerate(x):

            basefun = np.zeros(n,)
            basefun[i] = 1.0
            p_coef=np.polyfit(x, basefun, polyn_order)

            weight = np.polyval(np.polyder(p_coef, deriv_order), x_eval)
            weights[i] = weight
            polyder[ii] = polyder[ii] + f[ii,i]*weight
            trun_err_cent[ii]=((polyder[ii]-der )**2)**0.5    
            
####################biased collocated##########################            
    n=ll+r+1
    a=np.linspace(0,-15,1501)
    dx=10**a
    weights=np.zeros(n)
    basefun=np.zeros(n-1)
    deriv_order = 1
    polyn_order = n-1
    m=np.size(a)
    polyder=np.zeros(len(a))
    trun_err_bias=np.zeros(len(a))
    U=np.zeros((m,n))
    f=np.zeros((m,n))
    # print "U=",U
    # print "f=",f
    x=np.zeros(n)
    for ii,dxi in enumerate(dx):
        x[0] = -ll*dxi
        for j in range(1,n):
            xx = -(ll-j)*dxi
    #        print "x=",x,"xx=",xx
            x[j]=xx
        for i,xi in enumerate(x):
            U[ii,i] = (1/(np.cosh(xi))**2)*np.sin(5*xi+1.5) + 5*np.tanh(xi)*np.cos(5*xi+1.5)
            f[ii,i] = np.tanh(xi)*np.sin(xi+1.5)

        for i,val in enumerate(x):

            basefun = np.zeros(n,)
            basefun[i] = 1.0
            p_coef=np.polyfit(x, basefun, polyn_order)

            weight = np.polyval(np.polyder(p_coef, deriv_order), x_eval)
            weights[i] = weight
            polyder[ii] = polyder[ii] + f[ii,i]*weight
            trun_err_bias[ii]=((polyder[ii]-der )**2)**0.5                
    
    
    figure_folder = "../report/figures/"
    figure_name = "Collocated_centered"+"l"+str(l)+"r"+str(r)+"_biased"+"l"+str(ll)+"r"+str(r)+".jpg"
 
    
    
    fig = plt.figure(0, figsize=(figwidth,figheight))
    fig1 = fig.add_subplot(1,1,1)        
    plt.axes(fig1)
    fig1.grid('on',which='both')
    fig1.set_xlabel(r"$\frac{1}{\delta x}$ axis",fontsize=textFontSize)
    fig1.set_ylabel(r"truncation error($\epsilon$)",fontsize=textFontSize,rotation=90)    
    
    ec,=fig1.loglog(1/dx,trun_err_cent,'-g',linewidth=lineWidth,label = "eps centered")    
    es,=fig1.loglog(1/dx,trun_err_bias,'-b',linewidth=lineWidth, label = "eps biased")
    x1,=fig1.loglog(1/DX1,DX1,'--r',linewidth=1,label = "dx^1")
    x2,=fig1.loglog(1/DX2,DX2**2,'--r',linewidth=2,label = "dx^2")
    x3,=fig1.loglog(1/DX3,DX3**3,'--r',linewidth=3,label = "dx^3")
    x4,=fig1.loglog(1/DX4,DX4**4,'--r',linewidth=4,label = "dx^4")
    x5,=fig1.loglog(1/DX5,DX5**5,'--r',linewidth=4,label = "dx^5")
    handles, labels = fig1.get_legend_handles_labels()   
    fig1.legend(loc='upper right')




    figure_file_path= figure_folder + figure_name
    print "saving figure"
#    fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(figure_file_path)
    plt.close()   
    
    
    
    
    
    
    
    
    
########Staggered##########
    
LL=[1,2,3]
LL1=[1,1,1]
RR=[1,2,3]
for l,r,ll in zip(LL,RR,LL1):    
##########centered#################
    n=l+r
    a=np.linspace(0,-20,2001)
    dx1=10**a
    weights=np.zeros(n)
    basefun=np.zeros(n-1)
    deriv_order = 1
    polyn_order = n-1
    m=np.size(a)
    polyder=np.zeros(len(a))
    trun_err_stag_cent=np.zeros(len(a))
    Us=np.zeros((m,n))
    fs=np.zeros((m,n))
    x=np.zeros(n)
    for ii,dxi in enumerate(dx1):
        x[0] = -(l-1)*dxi-dxi/2
        for j in range(1,n):
            xx = x[0]+j*dxi
            x[j]=xx
        for i,xi in enumerate(x):
            Us[ii,i] = (1/(np.cosh(xi))**2)*np.sin(5*xi+1.5) + 5*np.tanh(xi)*np.cos(5*xi+1.5)
            fs[ii,i] = np.tanh(xi)*np.sin(xi+1.5)

        for i,val in enumerate(x):
            basefun = np.zeros(n,)
            basefun[i] = 1.0
            p_coef=np.polyfit(x, basefun, polyn_order)

            weight = np.polyval(np.polyder(p_coef, deriv_order), x_eval)
            weights[i] = weight
            polyder[ii] = polyder[ii] + fs[ii,i]*weight
            trun_err_stag_cent[ii]=((polyder[ii]-der )**2)**0.5
    
    
    



#######biased################
    n=ll+r
    a=np.linspace(0,-20,2001)
    dx1=10**a
    weights=np.zeros(n)
    basefun=np.zeros(n-1)
    deriv_order = 1
    polyn_order = n-1
    m=np.size(a)
    polyder=np.zeros(len(a))
    trun_err_stag_bias=np.zeros(len(a))
    Us=np.zeros((m,n))
    fs=np.zeros((m,n))
    x=np.zeros(n)
    for ii,dxi in enumerate(dx1):
        x[0] = -(ll-1)*dxi-dxi/2
        for j in range(1,n):
            xx = x[0]+j*dxi
            x[j]=xx
        for i,xi in enumerate(x):
            Us[ii,i] = (1/(np.cosh(xi))**2)*np.sin(5*xi+1.5) + 5*np.tanh(xi)*np.cos(5*xi+1.5)
            fs[ii,i] = np.tanh(xi)*np.sin(xi+1.5)

        for i,val in enumerate(x):
            basefun = np.zeros(n,)
            basefun[i] = 1.0
            p_coef=np.polyfit(x, basefun, polyn_order)

            weight = np.polyval(np.polyder(p_coef, deriv_order), x_eval)
            weights[i] = weight
            polyder[ii] = polyder[ii] + fs[ii,i]*weight
            trun_err_stag_bias[ii]=((polyder[ii]-der )**2)**0.5
    

    figure_folder = "../report/figures/"
    figure_name = "Staggered_centered"+"l"+str(l)+"r"+str(r)+"_biased"+"l"+str(ll)+"r"+str(r)+".jpg"
    
    
    fig = plt.figure(0, figsize=(figwidth,figheight))
    fig1 = fig.add_subplot(1,1,1)        
    plt.axes(fig1)
    fig1.grid('on',which='both')
    fig1.set_xlabel(r"$\frac{1}{\delta x}$ axis",fontsize=textFontSize)
    fig1.set_ylabel(r"truncation error($\epsilon$)",fontsize=textFontSize,rotation=90)    

    
    fig1.loglog(1/dx1,trun_err_stag_cent,'-g',linewidth=lineWidth,label = "eps centered")    
    fig1.loglog(1/dx1,trun_err_stag_bias,'-b',linewidth=lineWidth, label = "eps biased")
    fig1.loglog(1/DX1,DX1,'--r',linewidth=1,label = "dx^1")
    fig1.loglog(1/DX2,DX2**2,'--r',linewidth=2,label = "dx^2")
    fig1.loglog(1/DX3,DX3**3,'--r',linewidth=3,label = "dx^3")
    fig1.loglog(1/DX4,DX4**4,'--r',linewidth=4,label = "dx^4")
    fig1.loglog(1/DX5,DX5**5,'--r',linewidth=4,label = "dx^5")
    fig1.legend(loc='upper right')   



    figure_file_path= figure_folder + figure_name
    print "saving figure"
#    fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(figure_file_path)
    plt.close()   
    
        

