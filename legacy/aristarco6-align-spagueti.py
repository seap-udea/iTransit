from aristarchus import *
#############################################################
#INPUTS
#############################################################
typealignment=argv[1]
obsdir=argv[2]
images=argv[3]

#############################################################
#PARAMETERS
#############################################################
#TOLERANCE TO DETERMINE IF A BORDER CORRESPOND TO A CIRCLE
DRTOL=5E-2
NTHRES=3

#############################################################
#AT LEAST DOWNLOAD ERROR
#############################################################
html="""
<center>
<a href='%s/cmd.log' target='_blank'>Command</a> | 
<a href='%s/error.log' target='_blank'>Error output</a><br/>
</center>
"""%(obsdir,obsdir)
print html

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
APSs=[]
Rs=[]
dRs=[]
rms=[]
rmercs=[]
rspots=[]
rsps=[]
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
    qcomplete=True
    if dRR>DRTOL:qcomplete=False
    
    nborder=rs.shape[0]
    print>>stderr, "\t"*1,"Optimal threshold:"
    print>>stderr, "\t"*2,"Threshold = maxval / %.2f"%(xthresmin)
    print>>stderr, "\t"*2,"Rimage = %.2f +/- %.2f (%.5f)"%(R,dR,dRR)
    print>>stderr, "\t"*2,"Center Image = (%d,%d)"%(xcenter,ycenter)
    print>>stderr, "\t"*2,"Is solar disk complete? : ",qcomplete
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #CALCULATE CENTER AND RADIUS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xcenter,ycenter,R,dR=findCircleProperties(rs)
    print>>stderr,"\t"*1,"Fit result:"
    print>>stderr,"\t"*2,"R = %.2f +/- %.2f"%(R,dR)
    print>>stderr,"\t"*2,"Center = (%d,%d)"%(xcenter,ycenter)

    """
    #PLOT RESULT OF FITTING
    plt.figure(figsize=(8,8))
    plt.plot(rs[:,0],rs[:,1],'ro',ms=5,mec='none')
    for q in np.linspace(0,2*np.pi,1000):
        plt.plot([xcenter+R*np.cos(q)],[ycenter+R*np.sin(q)],'ko',ms=1,zorder=100)
    plt.xlim(xcenter-R,xcenter+R)
    plt.ylim(ycenter-R,ycenter+R)
    plt.savefig("tmp/c.png")
    #"""

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #FIND MERCURY POSITION RESPECT TO CENTER
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ptime=config["time"].split(":")
    time=float(ptime[0])+float(ptime[1])/60.+float(ptime[2])/3600.

    #DISTANCE TO CENTER
    ppos=config["posmercury"].split(",")
    xm=float(ppos[0]);ym=float(ppos[1])
    rm=np.sqrt((xm-xcenter)**2+(ym-ycenter)**2)/R
    rmerc=np.array([xm,ym])
    
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

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #FIND SUNSPOT POSITION RESPECT TO CENTER
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #DISTANCE TO CENTER
    if 'Use' not in config["possunspot"]:
        spos=config["possunspot"].split(",")
        xsp=float(spos[0]);ysp=float(spos[1])
        rsp=np.sqrt((xsp-xcenter)**2+(ysp-ycenter)**2)/R
        rspot=np.array([xsp,ysp])
        #APPARENT POSITION ANGLE OF THE SUNSPOT
        APS=np.arctan2((ysp-ycenter),(xsp-xcenter))*RAD
        print>>stderr,"\t"*1,"Sunspot position : r = %.5f, AP = %.2f deg"%(rsp,APS)
    else:
        rsp=-1
        APS=0
        rspot=no.array([0,0])

    times+=[time]
    Rs+=[R]
    dRs+=[dR]
    rms+=[rm]
    rsps+=[rsp]
    rps+=[rp]
    APs+=[AP]
    APSs+=[APS]

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #CROP IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    crop="%s/%s-crop-result.%s"%(obsdir,fname,ext)
    iR=roundFloat(R)
    ch=cw=2*iR
    
    # Composite image with white disk
    White=Image.new('RGBA',(cw,ch),"white")
    Crop=Image.new('RGBA',(cw,ch),"black")
    draw=ImageDraw.Draw(Crop)
    draw.ellipse((0,0)+Crop.size,fill=255)
    Crop=Image.composite(Crop,White,Crop)
    print>>stderr,"Size of canvas image:",ch,cw

    # Distance from corner to center of crop image
    cymin=iR-roundFloat(ycenter)
    cxmin=iR-roundFloat(xcenter)

    print>>stderr, "Position of corner (x,y): ",cxmin,cymin
    
    rxmin=0;rxmax=w
    rymin=0;rymax=h

    if cxmin<0:
        rxmin=-cxmin
        rxmax=rxmin+cw
        cxmin=0
        
    if cymin<0:
        rymin=-cymin
        rymax=rymin+ch
        cymin=0

    cxmax=cxmin+(rxmax-rxmin)
    cymax=cymin+(rymax-rymin)
        
    print>>stderr,"Before: Subimage (%d,%d,%d,%d) [%d,%d] -> Image (%d,%d,%d,%d) [%d,%d]"%(rxmin,rxmax,rymin,rymax,
                                                                                           rxmax-rxmin,rymax-rymin,
                                                                                           cxmin,cxmax,cymin,cymax,
                                                                                           cxmax-cxmin,cymax-cymin)
    dx=0;dy=0
    if cxmax>=cw:dx+=cxmax-cw
    if cymax>=ch:dy+=cymax-ch
    if rxmax>=w:dx+=rxmax-w
    if rymax>=h:dy+=rymax-h
    print>>stderr,"Corrections (dx,dy):",dx,dy

    rxmax-=dx;rymax-=dy
    cxmax-=dx;cymax-=dy

    print>>stderr,"After: Subimage (%d,%d,%d,%d) [%d,%d] -> Image (%d,%d,%d,%d) [%d,%d]"%(rxmin,rxmax,rymin,rymax,
                                                                                          rxmax-rxmin,rymax-rymin,
                                                                                          cxmin,cxmax,cymin,cymax,
                                                                                          cxmax-cxmin,cymax-cymin)

    Crop=np.asarray(Crop)
    Crop.flags.writeable=True
    Crop[cymin:cymax,cxmin:cxmax,:3]=Data[rymin:rymax,rxmin:rxmax,:3]
    imsave(crop,Crop)

    res=cw*ch
    if res<minres:
        Rcommon=R
        hcommon=ch
        wcommon=cw
        minres=res

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # MERCURY POSITION IN THE CROPPED IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print "Before:",rmerc
    rmerc=np.array([cxmin,cymin])+(rmerc-np.array([rxmin,rymin]))
    print "After:",rmerc
    rspot=np.array([cxmin,cymin])+(rspot-np.array([rxmin,rymin]))


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #ROTATE IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print>>stderr,"Rotating cropped image by an angle %.2f"%AP
    Rotated=ndimage.rotate(Crop,AP,reshape=False)

    # ROTATE MERCURY AND SPOT POSITION
    rcen=np.array([cw/2.,ch/2.])
    rmerc=rcen+rotatePoint(rmerc-rcen,AP*DEG)
    rspot=rcen+rotatePoint(rspot-rcen,AP*DEG)
    
    if Rotated.shape[2]>3:
        cond=Rotated[:,:,3]!=255
        Rotated[cond,3]=255

    rotated="%s/%s-rotated-result.%s"%(obsdir,fname,ext)
    plt.imsave(rotated,Rotated)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #SAVE PHP CONFIGURATION FILE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print>>stderr,"Mercury after rotation:",rmerc
    print>>stderr,"Sunspot after rotation:",rspot
    rmercs+=[rmerc]
    rspots+=[rspot]

    fp=open("%s/%s-align.php"%(obsdir,fname),"w")
    fp.write("""<?php
$center='%d,%d';
$R='%.2f';
$dR='%.2f';
$tm='%.8f';
$rm='%.5f';
$AP='%.2f';
$rsp='%.5f';
$APS='%.2f';
    ?>"""%(xcenter,ycenter,R,dR,time,rm,AP,rsp,APS))
    fp.close()
    i+=1

times=np.array(times)
Rs=np.array(Rs)
dRs=np.array(dRs)
rms=np.array(rms)
rsps=np.array(rsps)
rps=np.array(rps)
APs=np.array(APs)
APSs=np.array(APSs)
rmercs=np.array(rmercs)
rspots=np.array(rspots)
qs=np.zeros_like(rms)
ds=np.zeros_like(rms)

print>>stderr,"Mercury positions:\n",rmercs

#//////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////
# AUTOMATIC ALIGNMENT
#//////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////
if typealignment=='auto':

    #############################################################
    #SORT IMAGES ACCORDING TO RADIUS
    #############################################################
    irm=rms.argsort()[::-1]
    Rs_s=Rs[irm]
    dRs_s=dRs[irm]
    rms_s=rms[irm]
    rps_s=rps[irm]
    times_s=times[irm]
    APs_s=APs[irm]
    images_s=[images[i] for i in irm]

    rmercs_s=rmercs[irm]
    rspots_s=rspots[irm]


    print>>stderr,"\nPoint order:",irm

    print>>stderr,"Sorted:"
    print>>stderr,"Radii:",rms_s
 
    ts_s=np.zeros(nimages)
    for i in xrange(1,nimages):ts_s[i]=times_s[i]-times_s[0]
    print>>stderr,"Times:",ts_s

    print>>stderr,"Mercury sorted:\n",rmercs_s
    print>>stderr,"Spot sorted:\n",rmercs_s

    #############################################################
    #SEARCH FOR A SOLUTION
    #############################################################
    params=dict(ts=ts_s,rs=rms_s,verbose=0)
    solution=minimize(tdSlopeMinimize,[45*DEG],args=(params,))
    ql=solution["x"]

    qs_s,ds_s,B,m,b,r,logp,s=tdSlope(ql,params)
    print>>stderr,"Angles: ",qs_s*RAD

    #Convert from time ordered angles to original sorting
    for i in xrange(nimages):
        qs[irm[i]]=qs_s[i]
        ds[irm[i]]=ds_s[i]

    if logp<-2:
        status="Success"
    else:
        status="Failed"

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #PLOT SOLUTION
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    it=ts_s.argsort()
    ts_sort=ts_s[it]
    ds_sort=ds_s[it]
    tms=np.linspace(ts_sort.min(),ts_sort.max(),100)
    dms=m*tms+b
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(ts_sort,ds_sort,"rs",ms=20,mec='none')
    ax.plot(tms,dms,"r-",label=r"Linear fit, $\dot\theta$ = %.4f $\theta_\odot$/hour"%(m))
    ax.grid()
    ax.legend(loc='best')
    ax.set_xlabel("Time from reference position (hours)")
    ax.set_ylabel(r"Distance from reference position (apparent solar radii, $\theta_\odot$)")
    fig.savefig("%s/alignment-result.png"%obsdir)

    print>>stderr, "Unsorted:"
    print>>stderr, "Radii:",rms
    print>>stderr, "Angles: ",qs*RAD

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #TRANSIT DURATION
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a=1.0
    b=-2*rms_s[0]*np.cos(qs_s[0])
    c=-(1-rms_s[0]**2)
    d1,d2=quadraticEquation(a,b,c)
    dT=np.abs(d1-d2)
    print>>stderr, "Traverse distance (R): ",dT
    tT=dT/m
    print>>stderr, "Transit duration (hours): ",tT
    qrotatefirst=False
    
#//////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////
# SUNSPOT ALIGNMENT
#//////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////
else:
    status='Success'

    irm=np.arange(nimages)
    
    it=times.argsort()

    images_s=[images[i] for i in it]
    APs_s=APs[it]
    Rs_s=Rs[it]
    dRs_s=dRs[it]
    times_s=times[it]
    ts_s=np.zeros(nimages)
    for i in xrange(1,nimages):ts_s[i]=times_s[i]-times_s[0]
    qs_s=APSs[it]*DEG-APs[it]*DEG
    rms_s=rms[it]
    rmercs_s=rmercs[it]
    rspots_s=rmercs[it]

    # DETERMINE PROPERTIES
    
    m=0
    b=0
    r=0
    s=0
    logp=0
    B=0
    tT=0

    qrotatefirst=True

#############################################################
#ROTATE IMAGES
#############################################################

#CREATE ALIGNMENT BACKGROUND
Alignment=Image.new('RGBA',(wcommon,hcommon),"white")
Background=Image.new('RGBA',(wcommon,hcommon),"black")

j=0
for i in xrange(nimages):
    
    # Get image information
    image=images_s[i]
    print>>stderr, "File: ",image["name"]

    # Position of Mercury and the Sun
    rmerc=rmercs_s[i]
    rspot=rspots_s[i]

    rotated="%s/%s-rotated-result.%s"%(obsdir,image["name"],image["ext"])
    Rotated=imread(rotated)
    rh,rw=Rotated[:,:,0].shape
    final="%s/%s-final-result.%s"%(obsdir,image["name"],image["ext"])

    # Adjusting to a common size
    print>>stderr, "\t","Ajusting image to common resolution %sx%s"%(hcommon,wcommon)
    Rotated=imresize(Rotated,(hcommon,wcommon))
    imsave(rotated,Rotated)

    # Position in scaled image
    sx=(1.*wcommon)/rw;sy=(1.*hcommon)/rh
    rmerc[0]*=sx;rmerc[1]*=sy
    rspot[0]*=sx;rspot[1]*=sy
    rcen=np.array([rw/2.,rh/2.])

    if j==0 and not qrotatefirst:
        system("cp %s %s"%(rotated,final))
        Alignment=Image.fromarray(np.minimum(np.asarray(Alignment),np.asarray(Rotated)))
        rmercs_s[i]=rmerc
        rspots_s[i]=rspot
        print>>stderr,"Mercury after final rotation:",rmerc
        print>>stderr,"Spot after final rotation:",rspot
        j+=1
        continue

    print>>stderr, "\t","Rotating to final position image %d, r = %.2f, q = %.2f..."%(i,rms_s[i],qs_s[i]*RAD)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #ROTATE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Final=ndimage.rotate(Rotated,qs_s[i]*RAD,reshape=False)
    cond=Final[:,:,3]!=255
    Final[cond,3]=255
    plt.imsave(final,Final)

    # Rotate Mercury and Spot
    rmerc=rcen+rotatePoint(rmerc-rcen,qs_s[i])
    rspot=rcen+rotatePoint(rspot-rcen,qs_s[i])
    print>>stderr,"Mercury after final rotation:",rmerc
    print>>stderr,"Spot after final rotation:",rspot

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #CALCULATE RESULTING IMAGE5
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Alignment=Image.fromarray(np.minimum(np.asarray(Alignment),np.asarray(Final)))

    rmercs_s[i]=rmerc
    rspots_s[i]=rspot

    j+=1

alignment="%s/image-alignment-norot-result.%s"%(obsdir,image["ext"])
Alignment.save(alignment)
ah,aw=Alignment.size
rcen=np.array([aw/2.,ah/2.])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#CORD PROPERTIES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
it=times_s.argsort()

times_sort=times_s[it]
rmerc_sort=rmercs_s[it]

mp,bp,rp,pp,sp=linregress(rmerc_sort[:,0]-rcen[0],rmerc_sort[:,1]-rcen[1])
qimg=np.arctan(mp)
print "Cord:",mp,bp,rp,pp,sp,qimg*RAD,qs_s[0]*RAD

#==============================
#CALCULATE IMPACT PARAMETER
#==============================
B=np.abs(bp)/np.sqrt(mp**2+1)/Rcommon
print "Impact parameter:",B

#==============================
#CALCULATE COORD. OF CONTACTS
#==============================
a=-mp
b=1
c=bp
r=Rcommon

disc=r**2*(a**2+b**2)-c**2

xm=(a*c+b*np.sqrt(disc))/(a**2+b**2)+rcen[0]
ym=(b*c-a*np.sqrt(disc))/(a**2+b**2)+rcen[1]
rcont1=np.array([xm,ym])

xM=(a*c-b*np.sqrt(disc))/(a**2+b**2)+rcen[0]
yM=(b*c+a*np.sqrt(disc))/(a**2+b**2)+rcen[1]
rcont2=np.array([xM,yM])

#====================================
#CALCULATE DISTANCE TO CLOSEST POINTS
#====================================
dini=np.linalg.norm(rmerc_sort[0]-rcont1)
dlast=np.linalg.norm(rmerc_sort[-1]-rcont1)
if dini<dlast:
    dcont1=dini/Rcommon
    tcont1=times_sort[0]
else:
    dcont1=dlast/Rcommon
    tcont1=times_sort[-1]

dini=np.linalg.norm(rmerc_sort[0]-rcont2)
dlast=np.linalg.norm(rmerc_sort[-1]-rcont2)
if dini<dlast:
    dcont2=dini/Rcommon
    tcont2=times_sort[0]
else:
    dcont2=dlast/Rcommon
    tcont2=times_sort[-1]

print tcont1,dcont1
print tcont2,dcont2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#FINAL ROTATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alignment="%s/image-alignment-result.%s"%(obsdir,image["ext"])
Alignment=Alignment.rotate(qimg*RAD)

for i in xrange(nimages):
    
    # Position of Mercury and the Sun
    rmerc=rmercs_s[i]
    rspot=rspots_s[i]
    
    rmerc=rcen+rotatePoint(rmerc-rcen,qs_s[0])
    rspot=rcen+rotatePoint(rspot-rcen,qs_s[0])
    rmercs_s[i]=rmerc
    rspots_s[i]=rspot

# IN CASE OF SPOT ALIGNMENT
if typealignment=='spot':
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #CALCULATE DISTANCE BETWEEN POINTS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    it=times_s.argsort()
    j=0
    ds_s=np.zeros_like(times_s)
    for i in it:
        rmerc=rmercs_s[i]
        t=ts_s[i]
        if j==0:
            ds_s[i]=0.0
            rini=rmerc
        if j>0:
            ds_s[i]=np.linalg.norm(rmerc-rini)/Rcommon

        print t,rmerc,ds_s[i],ds_s[i]
        j+=1

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #SOLUTION ANGULAR MOTION
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m,b,r,p,s=linregress(ts_s[it],ds_s[it])
    logp=np.log(p)

    ts_sort=ts_s[it]
    ds_sort=ds_s[it]
    tms=np.linspace(ts_sort.min(),ts_sort.max(),100)
    dms=m*tms+b
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(ts_sort,ds_sort,"rs",ms=20,mec='none')
    ax.plot(tms,dms,"r-",label=r"Linear fit, $\dot\theta$ = %.4f $\theta_\odot$/hour"%(m))
    ax.grid()
    ax.legend(loc='best')
    ax.set_xlabel("Time from reference position (hours)")
    ax.set_ylabel(r"Distance from reference position (apparent solar radii, $\theta_\odot$)")
    fig.savefig("%s/alignment-result.png"%obsdir)

print>>stderr,"Mercury final position:\n",rmercs_s
print>>stderr,"Spot final position:\n",rspots_s
Alignment=Image.composite(Alignment,Background,Alignment)
Alignment.save(alignment)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#FINAL REPORT
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
report="%s/final-report.php"%(obsdir)
fp=open(report,"w")
fp.write("""<?php
$report=array(
""");
for i in xrange(nimages):

    image=images_s[i]
    name=image["name"]

    time=times_s[i]
    tm=ts_s[i]
    rm=rms_s[i]
    q=qs_s[i]
    rmerc=rmercs_s[i]
    rspot=rspots_s[i]

    fp.write("""
              '%s'=>array(
                          'time'=>'%s',
                          'rm'=>%.5f,
                          'tm'=>%.5f,
                          'q'=>%.2f,
                          'rmerc'=>array(%.2f,%.2f),
                          'rspot'=>array(%.2f,%.2f)
                         ),
    """%(name,time,rm,tm,q,rmerc[0],rmerc[1],rspot[0],rspot[1]))
fp.write(");\n?>");
fp.close()

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
"""%(k,crop,times_s[i],ts_s[i],rms_s[i],APs_s[i],qs[i]*RAD)
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
</div>
"""%(status,
     imageresult,imageresult,
     rps.mean(),B,m,tT,
     table,
     obsdir,
     m,b,r,s,logp)

print html
