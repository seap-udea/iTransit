from aristarchus import *
#############################################################
# CONSTANTS
#############################################################

#############################################################
# INPUTS
#############################################################
typealignment=argv[1]
imagefiles=argv[2:]
nimages=len(imagefiles)

#############################################################
# READ PROPERTIES
#############################################################
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

minres=1e100
for imagefile in imagefiles:

    print TAB*0,"Reading image:",imagefile

    #FILE PROPERTIES
    obsdir, image,fname,ext=fileProperties(imagefile)
    images+=[dict(file=imagefile,img=image,name=fname,ext=ext)]

    jfile="%s/scratch/%s-crop.json"%(obsdir,fname)
    conf=json2dict(jfile)

    times+=[conf["time"]]
    Rs+=[conf["R"]]
    dRs+=[conf["dR"]]

    rmercs+=[conf["rmerc"]]
    rms+=[conf["rm"]]
    APs+=[conf["AP"]]
    rps+=[conf["rp"]]

    rspots+=[conf["rspot"]]
    rsps+=[conf["rsp"]]
    APSs+=[conf["APS"]]

    
    crop="%s/scratch/%s-crop.png"%(obsdir,fname)
    cw,ch=Image.open(crop).size
    print TAB*1,"Image resolution (w,h) : ",cw,ch

    res=cw*ch
    if res<minres:
        Rcommon=conf["R"]
        hcommon=ch
        wcommon=cw
        minres=res

print TAB*0,"Minimum resolution (w,h) : ",wcommon,hcommon

times=np.array(times)
Rs=np.array(Rs)
dRs=np.array(dRs)

rmercs=np.array(rmercs)
rms=np.array(rms)
APs=np.array(APs)
rps=np.array(rps)

APSs=np.array(APSs)
rspots=np.array(rspots)
rsps=np.array(rsps)

#############################################################
#NEW PROPERTIES
#############################################################
#Subscript s stand for sorted arrays

# Time from reference
ts_s=np.zeros_like(rms)

# Angle from reference
qs=np.zeros_like(rms)
qs_s=np.zeros_like(rms)
# Distance from reference
ds=np.zeros_like(rms)
ds_s=np.zeros_like(rms)

if typealignment=="auto":

    #############################################################
    #AUTOALIGNMENT
    #############################################################

    #==================================================
    #SORT IMAGES ACPATHING TO RADIUS
    #==================================================
    irm=rms.argsort()[::-1]
    images_s=[images[i] for i in irm]

    times_s=times[irm]
    for i in xrange(1,nimages):ts_s[i]=times_s[i]-times_s[0]

    Rs_s=Rs[irm]
    dRs_s=dRs[irm]

    rmercs_s=rmercs[irm]

    rms_s=rms[irm]
    APs_s=APs[irm]
    rps_s=rps[irm]

    rspots_s=rspots[irm]

    #==================================================
    #SEARCH FOR A SOLUTION
    #==================================================
    print TAB*0,"Searching for an automatic alignment solution..."
    params=dict(ts=ts_s,rs=rms_s,verbose=0)
    solution=minimize(tdSlopeMinimize,[45*DEG],args=(params,))
    print TAB*1,"Solution status:",solution["success"]

    ql=solution["x"]
    qs_s,ds_s,b,m,y0,r,logp,s=tdSlope(ql,params)

    print TAB*0,"Motion properties from automatic alignment:"
    print TAB*1,"Angular velocity (R/h) = ",np.abs(m)
    print TAB*1,"Fit r = ",r
    print TAB*1,"log(p) = ",logp
    print TAB*1,"Impact parameter, b = ",b

    print TAB*0,"Solution:"
    print TAB*1,"Times: ",times_s
    print TAB*1,"Angles: ",qs_s*RAD

    qimg=qs_s[0]
    qs_s[0]=0.0

    qs_s=qs_s+APs_s*DEG
    print TAB*1,"Position angles: ",APs_s
    print TAB*1,"Final angles: ",qs_s*RAD

    #STORE ANGLES AND DISTANCES IN UNSORTED ARRAYS
    for i in xrange(nimages):
        qs[irm[i]]=qs_s[i]
        ds[irm[i]]=ds_s[i]
    
    #==================================================
    #PLOT SOLUTION
    #==================================================
    it=ts_s.argsort()
    ts_t=ts_s[it]
    ds_t=ds_s[it]

    tms=np.linspace(ts_t.min(),ts_t.max(),100)
    dms=m*tms+y0

    fig=plt.figure()
    ax=fig.gca()
    ax.plot(ts_t,ds_t,"rs",ms=20,mec='none')
    ax.plot(tms,dms,"r-",label=r"Linear fit, $\dot\theta$ = %.4f $\theta_\odot$/hour"%(m))
    ax.grid()
    ax.legend(loc='best')
    ax.set_xlabel("Time from reference position (hours)")
    ax.set_ylabel(r"Distance from reference position (apparent solar radii, $\theta_\odot$)")
    fig.savefig("%s/scratch/alignment-result.png"%obsdir)

else:

    #############################################################
    #SPOT ALIGNMENT
    #############################################################
    print TAB*0,"Using the sunspot to align..."

    #==================================================
    #SORT IMAGES ACPATHING TO TIME
    #==================================================
    it=times.argsort()
    images_s=[images[i] for i in it]

    times_s=times[it]
    for i in xrange(1,nimages):ts_s[i]=times_s[i]-times_s[0]

    Rs_s=Rs[it]
    dRs_s=dRs[it]

    rmercs_s=rmercs[it]
    rms_s=rms[it]

    APs_s=APs[it]
    rps_s=rps[it]

    rspots_s=rspots[it]
    
    #==================================================
    #ROTATION ANGLE
    #==================================================
    qs_s=APSs[it]*DEG-APs[it]*DEG
    status="success"

    print TAB*1,"Solution status:",status

    print TAB*0,"Solution:"
    print TAB*1,"Times: ",times_s
    print TAB*1,"Angles: ",qs_s*RAD

#############################################################
#ROTATE IMAGES
#############################################################
Aligned=Image.new('RGBA',(wcommon,hcommon),"white")
Background=Image.new('RGBA',(wcommon,hcommon),"black")

rcen=np.array([wcommon/2.,hcommon/2.])
j=0
for i in xrange(nimages):

    #IMAGE INFORMATION
    image=images_s[i]
    fname=image["name"]
    
    #POSITION OF MERCURY AND SUNSPOT
    rmerc=rmercs_s[i]
    rspot=rspots_s[i]
    
    #OPEN CROP IMAGE
    crop="%s/scratch/%s-crop.png"%(obsdir,fname)
    Crop=Image.open(crop)
    cw,ch=Crop.size
    sx=(1.*wcommon)/cw;sy=(1.*hcommon)/ch
    
    #RESIZE IMAGE
    Resize=Crop.resize((wcommon,hcommon))
    
    #ADJUST POSITION OF MERCURY & SUNSPOT
    rmerc[0]*=sx;rmerc[1]*=sy
    rspot[0]*=sx;rspot[1]*=sy
    
    #ROTATE IMAGE
    rotated="%s/scratch/%s-rotated.%s"%(obsdir,fname,ext)
    Rotated=Resize.rotate(qs_s[i]*RAD)
    Rotated=Image.composite(Rotated,Background,Rotated)
    Rotated.save(rotated)
    
    #ADJUST POSITION OF MERCURY & SUNSPOT
    rmerc=rcen+rotatePoint(rmerc-rcen,qs_s[i])
    rspot=rcen+rotatePoint(rspot-rcen,qs_s[i])

    #STORING MERCURY AND SUNSPOT POSITION
    rmercs_s[i]=rmerc
    rspots_s[i]=rspot

    #BUILD ALIGNED
    Aligned=Image.fromarray(np.minimum(np.asarray(Aligned),np.asarray(Rotated)))

    j+=1

aligned="%s/scratch/aligned.%s"%(obsdir,ext)
Aligned.save(aligned)

#############################################################
#PATH PROPERTIES
#############################################################
it=times_s.argsort()
times_sort=times_s[it]
rmerc_sort=rmercs_s[it]
mp,bp,rp,pp,sp=linregress(rmerc_sort[:,0]-rcen[0],rmerc_sort[:,1]-rcen[1])
qimg=np.arctan(mp)

print TAB*0,"Cord properties:"
print TAB*1,"Inclination angle: ",qimg*RAD
print TAB*1,"Correlation coefficient: ",rp
print TAB*1,"log(p-value): ",np.log(pp)

#############################################################
#ALIGNMENT DEBUG
#############################################################
debug="%s/scratch/aligned-debug.%s"%(obsdir,ext)
fig,ax=imgFigure(np.asarray(Aligned))

xs=np.linspace(0,cw,100)
ys=rcen[1]+mp*(xs-rcen[0])+bp
cond=np.sqrt((xs-rcen[0])**2+(ys-rcen[1])**2)<=Rcommon

ax.plot(rmercs_s[:,0],rmercs_s[:,1],'b+',ms=5,mfc='none')
ax.plot(rmercs_s[:,0],rmercs_s[:,1],'bo',ms=10,mfc='none')
ax.plot(rspots_s[:,0],rspots_s[:,1],'ro',ms=5,mfc='none')
ax.plot(xs[cond],ys[cond],'k--')

fig.savefig(debug)

#############################################################
#FINAL ALIGNMENT
#############################################################
aligned="%s/scratch/final-aligned.%s"%(obsdir,ext)
Aligned=Aligned.rotate(qimg*RAD)
Aligned=Image.composite(Aligned,Background,Aligned)
Aligned.save(aligned)

#ROTATE MERCURY AND SUNSPOT POSITION
rmercs_s=np.array([rcen+rotatePoint(rmercs_s[i]-rcen,qimg) for i in xrange(nimages)])
rspots_s=np.array([rcen+rotatePoint(rspots_s[i]-rcen,qimg) for i in xrange(nimages)])

#FIT PATH
it=times_s.argsort()
times_sort=times_s[it]
rmerc_sort=rmercs_s[it]
mp,bp,rp,pp,sp=linregress(rmerc_sort[:,0]-rcen[0],rmerc_sort[:,1]-rcen[1])
qimg=np.arctan(mp)

#############################################################
#ALIGNMENT DEBUG
#############################################################
debug="%s/scratch/final-aligned-debug.%s"%(obsdir,ext)
fig,ax=imgFigure(np.asarray(Aligned))

xs=np.linspace(0,cw,100)
ys=rcen[1]+mp*(xs-rcen[0])+bp
cond=np.sqrt((xs-rcen[0])**2+(ys-rcen[1])**2)<=Rcommon

ax.plot(rmercs_s[:,0],rmercs_s[:,1],'b+',ms=5,mfc='none')
ax.plot(rmercs_s[:,0],rmercs_s[:,1],'bo',ms=10,mfc='none')
ax.plot(rspots_s[:,0],rspots_s[:,1],'ro',ms=5,mfc='none')
ax.plot(xs[cond],ys[cond],'k--')

fig.savefig(debug)

#############################################################
# SAVE ALIGNMENT RESULT
#############################################################
it=times_s.argsort()
jfile="%s/scratch/aligned.json"%(obsdir)
properties=dict(
    #Images
    imagefiles=[images_s[i]["file"] for i in it],
    typealignment=typealignment,
    size=[wcommon,hcommon],
    rcen=rcen.tolist(),
    times=times_s[it].tolist(),

    #Mercury positions
    rmercs=rmercs_s[it].tolist(),

    #Sunspot positions
    rspots=rspots_s[it].tolist(),
    
    #Coord equation coefficients
    pathcoef=[mp,bp]
)
dict2json(properties,jfile)
