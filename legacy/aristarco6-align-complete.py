from aristarchus import *
#############################################################
#INPUTS
#############################################################
obsdir=argv[1]
images=argv[2]

#############################################################
#PARAMETERS
#############################################################

#TOLERANCE TO DETERMINE IF A BORDER CORRESPOND TO A CIRCLE
DRTOL=5E-2
NTHRES=3

#############################################################
#ANALYZE IMAGES
#############################################################
print>>stderr,"ARISTARCHUS CAMPAIGN 6 ANALYSIS\n"

limages=images.split(",")
nimages=len(limages)
ipos=np.arange(nimages)
print>>stderr,"Number of images: ",nimages

times=[]
APs=[]
rms=[]
rps=[]
images=[]
i=0
minres=1e100
for image in limages:

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #SPLIT NAME
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parts=image.split(".")
    ext=parts[-1]
    fname="".join(parts[:-1])
    images+=[dict(img=image,name=fname,ext=ext)]
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #READ PHP FILE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    php="%s/%s.php"%(obsdir,fname)
    config=parsePhp(php)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #READ IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Data=imread("%s/%s"%(obsdir,image))
    Mono=Data[:,:,0]
    h,w=Mono.shape
    X,Y=np.meshgrid(np.arange(w),np.arange(h))
    maxval=Mono.max()
    print>>stderr, "\nImage %d: '%s', resolution %d x %d..."%(i,image,h,w)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #GET BORDER OF THE SUN
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xthresmin=0;dRR=1e100

    #DETERMINE THE OPTIMAL THRESHOLD
    for xthres in np.linspace(2.0,10.0,NTHRES):
        border=detectBorder(Mono,threshold=maxval/xthres)
        rs=border["rs"]
        xc,yc=border["centroid"]    
        Rs=np.sqrt((rs[:,0]-xc)**2+(rs[:,1]-yc)**2)
        Rmean=Rs.mean()
        Rstd=Rs.std()
        dRr=Rstd/(1.*Rmean)
        if dRr<dRR:
            xthresmin=xthres
            bordermin=border
            xcenter=xc
            ycenter=yc
            R=Rmean
            dR=Rstd
            dRR=dRr
        #print>>stderr, "Threshold = %d, Rmean,Rstd,dR = "%(maxval/xthres),Rmean,Rstd,dRr

    nborder=rs.shape[0]
    print>>stderr, "\t"*1,"Optimal solution:"
    print>>stderr, "\t"*2,"Threshold = maxval / %.2f"%(xthresmin)
    print>>stderr, "\t"*2,"R = %.2f +/- %.2f (%.5f)"%(R,dR,dRR)
    print>>stderr, "\t"*2,"Center = (%d,%d)"%(xcenter,ycenter)
    
    #"""
    plt.figure(figsize=(8,8))
    plt.plot(rs[:,0],rs[:,1],'ro',ms=5,mec='none')
    #"""

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #IF SOLAR DISK IS NOT COMPLETE FIND CENTER
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qcomplete=True
    if dRR>DRTOL:qcomplete=False
    if qcomplete:
        print>>stderr,"\t"*1,"The Sun is complete"
    if not qcomplete:
        print>>stderr,"\t"*1,"The Sun has been chopped"
        x1=rs[nborder/10*0,0]
        y1=rs[nborder/10*0,1]
        x2=rs[nborder-1,0]
        y2=rs[nborder-1,1]
        
        #Secant midpoint
        xs=(x1+x2)/2.;ys=(y1+y2)/2.

        #Direction
        rl=np.cross([0,0,1],[x2-x1,y2-y1,0])
        rl=rl/norm(rl)
        rl=[rl[0],rl[1]]
        args=dict(rs=rs,rl=rl,rm=[xs,ys])

        dR=1E100
        for l in np.linspace(-w,w,100):
            Rr,dRr=radiusPoints(l,args)
            if dRr<dR:
                dR=dRr
                R=Rr
                lmin=l

        print>>stderr,"\t"*1,"Radial dispersion at solution = ",R,dR
        dRR=dR/(1.*R)
        xcenter=xs+lmin*rl[0]
        ycenter=ys+lmin*rl[1]

        #"""
        plt.plot([x1],[y1],'bs',ms=10)
        plt.plot([x2],[y2],'rs',ms=10)
        plt.plot([xs],[ys],'gs',ms=10)
        plt.plot([x1,x2],[y1,y2],'k-')
        plt.axhline(ycenter)
        plt.axvline(xcenter)
        #"""

        print>>stderr,"\t"*1,"After recalculation:"
        print>>stderr,"\t"*2,"R = %.2f +/- %.2f (%.5f)"%(R,dR,dRR)
        print>>stderr,"\t"*2,"Center = (%d,%d)"%(xcenter,ycenter)

    #"""
    plt.axhline(yc,color='r')
    plt.axvline(xc,color='r')
    extreme=max(w,h)
    plt.xlim((0,extreme))
    plt.ylim((0,extreme))
    plt.savefig("tmp/c.png")
    #"""

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #CROP IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dp=0
    if qcomplete:
        print>>stderr,"\t","Cropping a complete image"
        xmin=int(np.round(xcenter-R))-dp;xmax=int(np.round(xcenter+R))+dp
        ymin=int(np.round(ycenter-R))-dp;ymax=int(np.round(ycenter+R))+dp
        reshape=False
    else:
        print>>stderr,"\t","Cropping a partial image"
        cond=Mono>=maxval/xthresmin
        xs=X[cond]
        ys=Y[cond]
        xmin=xs.min();xmax=xs.max()
        ymin=ys.min();ymax=ys.max()
        reshape=False

    xmin,xmincorner=cropCoord(xmin,w)
    ymin,ymincorner=cropCoord(ymin,h)
    xmax,xmaxcorner=cropCoord(xmax,w)
    ymax,ymaxcorner=cropCoord(ymax,h)

    #Cropped image
    Crop=Data[ymin:ymax,xmin:xmax,:]
    crop="%s/%s-crop-result.%s"%(obsdir,fname,ext)
    images[i]["crop"]=crop
    plt.imsave(crop,Crop)
    ch,cw=Crop[:,:,0].shape
    xcropcenter=roundFloat(xcenter-xmin)
    ycropcenter=roundFloat(ycenter-ymin)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #ADD CANVAS TO AN INCOMPLETE IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if not qcomplete:
        #Canvas image

        #FIND EXTREME OF COMPLETE IMAGE
        cxmin=-int(round((xcropcenter-R)))
        cymin=-int(round((ycropcenter-R)))

        print>>stderr,"Center:",xcropcenter,ycropcenter
        print>>stderr,"Corner before:",cxmin,cymin

        if cxmin<0:cxmin=0
        if cymin<0:cymin=0

        #print>>stderr,"Corner after:",cxmin,cymin

        iR=int(round(R))
        compw=2*iR
        comph=2*iR

        #print>>stderr,"Size before:",compw,comph

        compw=max(compw,cxmin+iR+100)
        comph=max(comph,cymin+iR+100)

        #print>>stderr,"Size after:",compw,comph
        #print>>stderr,"Size:",compw,comph

        #CREATE BLANK IMAGE
        Complete=np.asarray(Image.new('RGB',(compw,comph)))
        Complete.flags.writeable=True

        #print>>stderr,"Shape:",Complete.shape

        #SAVE IMAGE INTO BLANK
        Complete[cymin:cymin+ch,cxmin:cxmin+cw,:]=Crop[:,:,:3]
        imsave(crop,Complete)
        Crop=Complete

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #FIND MERCURY POSITION RESPECT TO CENTER
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ptime=config["time"].split(":")
    time=float(ptime[0])+float(ptime[1])/60.+float(ptime[2])/3600.

    #DISTANCE TO CENTER
    ppos=config["posmercury"].split(",")
    xm=float(ppos[0]);ym=float(ppos[1])
    rm=np.sqrt((xm-xcenter)**2+(ym-ycenter)**2)/R
    
    #APPARENT SIZE OF MERCURY
    Dp=float(ppos[2])
    rp=Dp/2.0/R

    #APPARENT POSITION ANGLE
    """
    AP : Angle between the rightmost point of the Sun and the line towards
    Mercury meaured in the clockwise direction.
    """
    AP=np.arctan2((ym-ycenter),(xm-xcenter))*RAD
    print>>stderr,"\t"*1,"Mercury position : r = %.5f, AP = %.2f deg"%(rm,AP)

    times+=[time]
    rms+=[rm]
    rps+=[rp]
    APs+=[AP]

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #ROTATE IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rotated=ndimage.rotate(Crop,AP,reshape=reshape)

    if Rotated.shape[2]>3:
        cond=Rotated[:,:,3]!=255
        Rotated[cond,3]=255

    if qcomplete:
        rMono=Rotated[:,:,0]
        rh,rw=rMono.shape
        X,Y=np.meshgrid(np.arange(rw),np.arange(rh))
        cond=rMono>=maxval/xthresmin
        xs=X[cond]
        ys=Y[cond]
        xmin=xs.min();xmax=xs.max()
        ymin=ys.min();ymax=ys.max()
        Rotated=Rotated[ymin:ymax,xmin:xmax,:]

    rh,rw=Rotated[:,:,0].shape

    res=rw*rh
    if res<minres:
        hcommon=rh
        wcommon=rw
        minres=res

    rotated="%s/%s-rotated-result.%s"%(obsdir,fname,ext)
    plt.imsave(rotated,Rotated)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #SAVE PHP CONFIGURATION FILE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp=open("%s/%s-align.php"%(obsdir,fname),"w")
    fp.write("""<?php
$center='%d,%d';
$cropceneter='%d,%d';
$R='%.2f';
$dR='%.2f';
$tm='%.8f';
$rm='%.5f';
$AP='%.2f';
?>"""%(xcenter,ycenter,xcropcenter,ycropcenter,R,dR,time,rm,AP))
    fp.close()
    i+=1

times=np.array(times)
rms=np.array(rms)
rps=np.array(rps)
APs=np.array(APs)

#############################################################
#SORT IMAGES ACCORDING TO RADIUS
#############################################################
irm=rms.argsort()[::-1]
rms_s=rms[irm]
rps_s=rps[irm]
times_s=times[irm]
APs_s=APs[irm]
images=[images[i] for i in irm]

print>>stderr, "\nPoint order:",irm

rs=rms[irm]
print>>stderr, "Radii:",rs

ts=np.zeros(nimages)
for i in xrange(1,nimages):ts[i]=times_s[i]-times_s[0]
print>>stderr, "Times:",ts

#############################################################
#SEARCH FOR A SOLUTION
#############################################################
params=dict(ts=ts,rs=rs,verbose=0)
solution=minimize(tdSlopeMinimize,[45*DEG],args=(params,))
ql=solution["x"]

qs,ds,B,m,b,r,logp,s=tdSlope(ql,params)
print>>stderr, "Angles: ",qs*RAD

if logp<-2:
    status="Success"
else:
    status="Failed"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#PLOT SOLUTION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
it=ts.argsort()
ts_s=ts[it]
ds_s=ds[it]
tms=np.linspace(0,ts_s[-1],100)
dms=m*tms+b
fig=plt.figure()
ax=fig.gca()
ax.plot(ts,ds,"rs",ms=20,mec='none')
ax.plot(tms,dms,"r-",label=r"Linear fit, $\dot\theta$ = %.4f $\theta_\odot$/hour"%(m))

ax.grid()
ax.legend(loc='best')
ax.set_xlabel("Time from most external position (hours)")
ax.set_ylabel(r"Distance between points (apparent solar radii, $\theta_\odot$)")
fig.savefig("%s/alignment-result.png"%obsdir)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#TRANSIT DURATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1.0
b=-2*rs[0]*np.cos(qs[0])
c=-(1-rs[0]**2)
d1,d2=quadraticEquation(a,b,c)
dT=np.abs(d1-d2)
print>>stderr, "Traverse distance (R): ",dT
tT=dT/m
print>>stderr, "Transit duration (hours): ",tT

#############################################################
#ROTATE IMAGES
#############################################################
j=0
for i in irm:
    
    # Get image information
    image=images[i]
    print>>stderr, "File: ",image["name"]

    rotated="%s/%s-rotated-result.%s"%(obsdir,image["name"],image["ext"])
    Rotated=imread(rotated)
    final="%s/%s-final-result.%s"%(obsdir,image["name"],image["ext"])

    # Adjusting to a common size
    print>>stderr, "\t","Ajusting image to common resolution %sx%s"%(hcommon,wcommon)
    Rotated=imresize(Rotated,(hcommon,wcommon))
    plt.imsave(rotated,Rotated)

    if j==0:
        Alignment=Rotated
        system("cp %s %s"%(rotated,final))
        j+=1
        continue

    print>>stderr, "\t","Rotating to final position image %d, r = %.2f, q = %.2f..."%(i,rs[i],qs[i]*RAD)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #LOAD ROTATED IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Final=ndimage.rotate(Rotated,qs[i]*RAD,reshape=False)
    cond=Final[:,:,3]!=255
    Final[cond,3]=255
    plt.imsave(final,Final)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #CALCULATE RESULTING IMAGE5
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Alignment=np.minimum(Alignment,Final)
    j+=1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#FINAL ROTATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Alignment=ndimage.rotate(Alignment,qs[0]*RAD,reshape=False)
cond=Alignment[:,:,3]!=255
Alignment[cond,3]=255
alignment="%s/image-alignment-result.%s"%(obsdir,image["ext"])
plt.imsave(alignment,Alignment)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#HTML
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table=""
k=1
for i in irm:
    image=images[i]
    fname=image["name"]
    ext=image["ext"]
    crop="%s/%s-crop-result.%s"%(obsdir,fname,ext)
    table+="""
<!-- CROPPED IMAGE -->
<tr>
<td class="value">%d</td>
<td>
<img src='%s' height='200px'/>
</td>
<!-- BASIC PROPERTIES -->
<td class="value">%.4f</td>
<td class="value">%.2f</td>
<td class="value">%.2f</td><td class="value">%.1f</td>
<td class="value">%.1f</td>
</tr>
"""%(k,crop,times_s[i],ts[i],rms_s[i],APs_s[i],qs[i]*RAD)
    k+=1
    
    
imageresult="%s/image-alignment-result.png"%obsdir

html="""
<center style='font-size:1.2em'>
<a href=\"JavaScript:void(null)\" onclick=\"$('#detail').toggle()\">View/Hide Analysis Results</a>
</center>

<div id='detail' style='display:none'>
<h3>Analysis results</h3>

<h4>Result of the alignment procedure</h4>

<p>
The result was: <b>%s</b>
</p>

This is the resulting image:

<div style="text-align:center">
<a href='%s' target='_blank'>
<img src='%s' width="60%%"/>
</a>
</div>

<p>
After analysing the image we determined the following observable
transit parameters:
</p>

<ul>
<li>Radius of the planet: r/R = %.5f</li>
<li>Impact parameter: b/R = %.5f</li>
<li>Apparent angular velocity: v/R = %.5f h<sup>-1</sup></li>
<li>Transit duration: t<sub>T</sub> = %.4f h</li>
</ul>

<h4>Angles</h4>

<center>

<style>
td.value{
vertical-align:top;
text-align:center;
}
</style>

<table border=1px style="margin-left:10%%;width:60%%;">
<tr>
<th width=0%%>#</th>
<th width=0%%>Cropped image</th>
<th width=30%%>Time</sub></th>
<th width=30%%>t-t<sub>1</sub></th>
<th width=30%%>r/R</th>
<th width=30%%>PA (<sup>o</sup>)</th>
<th width=30%%>&theta; (<sup>o</sup>)</th>
</tr>
%s
</table>
</center>

<p>

Where R is the solar apparent radius and AP is the position angle with
respect to the apparent west side (rightmost side of the image)
measured in the clockwise direction.

</p>

<h4>Fitting</h4>

Result of the fit:

<div style="text-align:center">
<img src='%s/alignment-result.png' width="60%%"/>
</div>

The fitting resulting parameters were:

<ul>
<li>m = %.5f</li>
<li>b = %.5f</li>
<li>r = %.5f</li>
<li>s = %.5f</li>
<li>log p = %.1f</li>
</ul>

<h4>Other download</h4>

<a href='%s/cmd.log' target='_blank'>Command</a> | 
<a href='%s/error.log' target='_blank'>Output</a><br/>
</div>
"""%(status,
     imageresult,imageresult,
     rps.mean(),B,m,tT,
     table,
     obsdir,
     m,b,r,s,logp,
     obsdir,obsdir)

print html
