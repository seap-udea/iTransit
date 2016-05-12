#!/usr/bin/env python
#############################################################
# EXTRENAL LIBRARIES
#############################################################
from matplotlib import use
use('Agg')
import matplotlib.pylab as plt
import matplotlib.cm as cm
import numpy as np
from sys import argv,exit,stderr
from os import system
from commands import getoutput as System
import json
from scipy.misc import *
from scipy.optimize import minimize,bisect
from scipy.stats import linregress
from scipy import ndimage
import re
import itertools as it
from PIL import Image,ImageDraw


#############################################################
# MACROS
#############################################################
norm=np.linalg.norm
TAB="\t"

#############################################################
# CONSTANTS
#############################################################
DEG=np.pi/180
RAD=180/np.pi

#############################################################
# ROUTINES
#############################################################
def parsePhp(file):
    f=open(file,"r")
    config=dict()
    for line in f:
        line=line.strip()
        if '<' in line or\
           '>' in line:
            continue
        parts=line.split("=");
        var=parts[0][1:]
        value=parts[1]
        value=value.replace("\"","");
        value=value.replace("\'","");
        value=value.replace(";","");
        config[var]=value
    return config

def detectBorder(imagen,threshold=50,ptol=3):

    #Get image properties
    ch,cw=imagen.shape

    xs=np.arange(cw)
    ys=np.arange(ch)
    rs=[]

    #Get points in horizontal sections
    period=1
    for i in xrange(ch):
        if (i%period):continue
        cond=imagen[i,:]>=threshold
        row=xs[cond]
        if len(row)==0:continue
        xmin=row[0]
        xmax=row[-1]
        if (xmin>ptol):rs+=[[xmin,i]]
        if (cw-xmax)>ptol:rs+=[[xmax,i]]
    #Get points in vertical sections
    for j in xrange(cw):
        if (j%period):continue
        cond=(imagen[:,j]>=threshold)
        col=ys[cond]
        if len(col)==0:continue
        ymin=col[0]
        ymax=col[-1]
        if (ymin>ptol):rs+=[[j,ymin]]
        if (ch-ymax)>ptol:rs+=[[j,ymax]]
    rs=np.array(rs)

    #Sort according to abcisa
    rslist=rs.tolist()
    rslist.sort(key=lambda row:row[1])
    rs=np.array(rslist)
    #Remove unique
    rs=np.vstack({tuple(row) for row in rs})

    #Compute the centroid
    xc=rs[:,0].mean()
    yc=rs[:,1].mean()

    #Sort according to angle with respect to centroid
    rslist=rs.tolist()
    rslist.sort(key=lambda row:np.arctan2((row[1]-yc),(row[0]-xc)))
    rs=np.array(rslist)
    
    #Return dict
    border=dict(rs=rs,centroid=[xc,yc])

    return border

def radiusPoints(l,args):
    rs=args["rs"]
    rl=args["rl"]
    rm=args["rm"]

    xc=rm[0]+l*rl[0]
    yc=rm[1]+l*rl[1]

    Rs=np.sqrt((rs[:,0]-xc)**2+(rs[:,1]-yc)**2)
    Rsmean=Rs.mean()
    dRs=Rs.std()
    return Rsmean,dRs

def quadraticEquation(a,b,c):
    disc=b**2-4*a*c
    if disc<0:
        print "Complex roots"
        return 0,0
    else:
        disc=np.sqrt(disc)
        x1=(-b+disc)/(2*a)
        x2=(-b-disc)/(2*a)
        return x1,x2

def centralAngle(r0,q0,r):
    # QUADRATIC EQUATION COEFFICIENTS
    a=r**2/np.sin(q0)**2
    b=-2*r0*r
    c=(r0**2+r**2-a)
    
    # SOLUTION TO QUADRATIC EQUATION
    cq1,cq2=quadraticEquation(a,b,c)
    q1=np.arccos(cq1);q2=np.arccos(cq2)
    
    # SORTED OUTPUT
    qmin=min(q1,q2)
    qmax=max(q1,q2)
    return qmin,qmax

def tRatio(q2,params):

    rs=params["rs"]
    ts=params["ts"]
    verbose=params["verbose"]

    # RADII AND TIMES
    r0,r1,r2=rs
    t0,t1,t2=ts

    # DIRECTION
    direction=np.sign(t2)
    if verbose:print "Direction = ",direction

    if verbose:print "Testing q2 = %.4f"%(q2*RAD)

    # CALCULATE CONDITIONS FOR POSITION 2
    d2=np.sqrt(r0**2+r2**2-2*r0*r2*np.cos(q2))
    if verbose:print "t2 = %.5f, d2 = %.6f, q2 = %.4f"%(t2,d2,q2*RAD)

    # CALCULATE q0 FOR THIS CONFIGURATION
    q0=np.arcsin(np.sin(q2)*r2/d2)
    if verbose:print "q0 = %.4f"%(q0*RAD)

    # CALCULATE CONDITION FOR POSITION 1
    if q0==0:
        #T1 HAPPENS BEFORE THAN T2
        if direction*(t1-t2)<0:
            d1=r0-r1
            q1=0.0
        else:
            d1=r0+r1
            q1=np.pi
    else:
        """
        # VERIFICATION OF THE INVERSE:
        # QUADRATIC EQUATION COEFFICIENTS
        a=r2**2/np.sin(q0)**2
        b=-2*r0*r2
        c=(r0**2+r2**2-a)

        # SOLUTION TO QUADRATIC EQUATION
        cq2_1,cq2_2=quadraticEquation(a,b,c)
        if verbose:print np.arccos(cq2_1)*RAD,np.arccos(cq2_2)*RAD
        exit(0)
        """

        # QUADRATIC EQUATION COEFFICIENTS
        a=r1**2/np.sin(q0)**2
        b=-2*r0*r1
        c=(r0**2+r1**2-a)

        # SOLUTION TO QUADRATIC EQUATION
        cq1_1,cq1_2=quadraticEquation(a,b,c)
        q1_1=np.arccos(cq1_1);q1_2=np.arccos(cq1_2)

        qmin=min(q1_1,q1_2)
        qmax=max(q1_1,q1_2)

        if direction*(t1-t2)<0:
            q1=qmin
        else:
            q1=qmax

        d1=np.sqrt(r0**2+r1**2-2*r0*r1*np.cos(q1))
        if verbose:print "t1 = %.5f, d1 = %.6f, q1 = %.4f"%(t1,d1,q1*RAD)

    f=(d2/d1)-(t2/t1)

    # COMPARE DISTANCES WITH TIMES
    if verbose:print "(d2/d1) = %.4f = (t2/t1) = %.4f ? => f = %.4f "%(d2/d1,t2/t1,f)
    if verbose:raw_input()

    return f,q0,q1,q2,d1,d2
    
def tRatioBisection(q2,params):
    values=tRatio(q2,params)
    return values[0]

def tRatioMinimize(q2,params):
    values=tRatio(q2,params)
    return abs(values[0])

def tdSlope(ql,params):
    
    rs=params["rs"]
    ts=params["ts"]
    verbose=params["verbose"]

    nimages=len(rs)
    qs=np.zeros_like(rs)
    ds=np.zeros_like(qs)

    qs[-1]=ql
    ds[-1]=np.sqrt(rs[0]**2+rs[-1]**2-2*rs[0]*rs[-1]*np.cos(qs[-1]))
    qs[+0]=np.arcsin(np.sin(qs[-1])*rs[-1]/ds[-1])

    direction=np.sign(ts[-1])
    if verbose:print "\nDirection of images: ",direction
    if verbose:print "\nInclination angle q0: ",qs[0]*RAD
    if verbose:print "\nProperties of last point ql,dl: ",qs[-1]*RAD,ds[-1]

    for i in xrange(1,nimages-1):

        if verbose:print "\nPoint %d:"%i
        if verbose:print "\t","ti = %.3f, ri = %.3f, ti-tl = %.3f"%(ts[i],rs[i],ts[i]-ts[-1])

        # QUADRATIC COEFFICIENTS
        a=rs[i]**2/np.sin(qs[0])**2
        b=-2*rs[0]*rs[i]
        c=(rs[0]**2+rs[i]**2-a)

        # SOLUTION TO QUADRATIC EQUATION
        c1,c2=quadraticEquation(a,b,c)
        q1=np.arccos(c1)
        q2=np.arccos(c2)
        if verbose:print "\t","Solution: ",q1*RAD,q2*RAD

        # CHOOSE A SOLUTION ACCORDING TO TIME
        if direction*(ts[+i]-ts[-1])<0:qs[i]=q1
        else:qs[i]=q2
        if verbose:print "\t","Solution chosen: ",qs[i]*RAD

        # CALCULATE DISTANCE TRAVERSED
        ds[i]=np.sqrt(rs[0]**2+rs[i]**2-2*rs[0]*rs[i]*np.cos(qs[i]))
        if verbose:print "\t","Distance: ",ds[i]

    if verbose:print "\nSummary:"
    if verbose:print "\t","Times:",ts
    if verbose:print "\t","Angles:",qs*RAD
    if verbose:print "\t","Distances:",ds

    # FIT VALUES
    it=ts.argsort()
    ts_s=ts[it]
    ds_s=ds[it]

    # LINEAR REGRESSION
    m,b,r,p,s=linregress(ts_s,ds_s)
    logp=np.log10(p)
    tms=np.linspace(0,ts_s[-1],100)
    dms=m*tms+b
    if verbose:print "\nLinear Regression:"
    if verbose:print "\t","r = ",r
    if verbose:print "\t","log(p) = ",logp

    # IMPACT PARAMETER
    B=rs[0]*np.sin(qs[0])
    if verbose:print "\n","Impact parameter: ",B

    return qs,ds,B,m,b,r,logp,s

def tdSlopeMinimize(ql,params):
    qs,ds,B,m,b,r,logp,s=tdSlope(ql,params)
    return logp

def cropCoord(x,w):
    if x<=1:return 0,'l'
    elif x>=w:return w,'r'
    else:return x,'c'
    
def roundFloat(x):
    return int(np.round(x))

def sortAbscisas(rs,xc,yc):
    rslist=rs.tolist()
    rslist.sort(key=lambda row:row[1])
    rs=np.array(rslist)

    #Sort according to angle with respect to centroid                                                   rslist=rs.tolist()
    rslist.sort(key=lambda row:np.arctan2((row[1]-yc),(row[0]-xc)))
    rs=np.array(rslist)

    return rs

def fitCircle(x,params):
    rs=params["rs"]
    xc=x[0]
    yc=x[1]
    R=x[2]

    f=(np.sqrt((rs[:,0]-xc)**2+(rs[:,1]-yc)**2)-R)**2
    return f.sum()

def findCircleProperties(rs):
    npoints=len(rs)
    params=dict(rs=rs)
    xcenter=rs[:,0].mean()
    ycenter=rs[:,1].mean()
    R=max((rs[:,0].max()-rs[:,0].min())/2,
          (rs[:,1].max()-rs[:,1].min())/2)
    solution=minimize(fitCircle,[xcenter,ycenter,R],args=(params,))
    x=solution["x"]
    xcenter=x[0]
    ycenter=x[1]
    R=x[2]
    dR=np.sqrt(fitCircle(x,params))/npoints
    return xcenter,ycenter,R,dR

def rotatePoint(r,theta):
    cost=np.cos(theta)
    sint=np.sin(theta)
    R=np.array([[cost,+sint],[-sint,cost]])
    rp=R.dot(r)
    return rp

def php2json(includefile):
    out=System("""php -r 'include("%s");echo json_encode($report);'"""%includefile)
    dictionary=json.loads(out)
    return dictionary

def imgFigure(Data):
    h,w=Data[:,:,0].shape
    fig=plt.figure(figsize=(w/100.,h/100.))
    ax=fig.add_axes([0,0,1,1])
    ax.imshow(Data,extent=[-1,w,h,1])
    ax.axis('off')
    ax.set_xlim((0,w))
    ax.set_ylim((h,0))
    return fig,ax

def boxFig(r,params):

    verbose=1

    #Subimage
    S=params["S"]
    
    #Size
    h,w=S.shape
    
    #Variables
    x=r[0]
    y=r[1]
    R=r[2]

    if (x**2+y**2)>(w**2+h**2)-R:return 1e100

    if verbose:print "Testing %.0f,%.0f, R = %.2f"%(x,y,R)
    if verbose:print "S:"%(S)

    #Mask
    B=255*np.ones((h,w))
    X,Y=np.meshgrid(np.arange(w),np.arange(h))
    cond=np.sqrt((X-x)**2+(Y-y)**2)<=R
    cond=(np.abs(X-x)<R)*(np.abs(Y-y)<R)
    
    if verbose:print "Cond:\n",cond
    
    B[cond]=0
    if verbose:print "B:\n",B

    #Difference
    D=np.abs(S-B)

    if verbose:print "D:\n",D

    f=D.sum()
    if verbose:
        print "Result:",f
        raw_input()

    return f

def dict2json(dic,jfile):
    f=open(jfile,"w")
    f.write(json.dumps(dic))
    f.close()

def json2dict(jfile):
    f=open(jfile,"r")
    string=f.readline()
    dic=json.loads(string)
    return dic
    
def fileProperties(mfile):
    name=System("basename %s"%mfile)
    mdir=System("dirname %s"%mfile)
    parts=name.split(".")
    ext=parts[-1]
    fname="".join(parts[:-1])
    return mdir,name,fname,ext
    
    
