from aristarchus import *
#############################################################
# CONSTANTS
#############################################################
DRTOL=5E-2
NTHRES=3

#############################################################
# INPUTS
#############################################################
imagefile=argv[1]

image=System("basename %s"%imagefile)
obsdir=System("dirname %s"%imagefile)
parts=image.split(".")
ext=parts[-1]
fname="".join(parts[:-1])
system("mkdir -p %s/scratch"%obsdir)

print "Analysing image:",fname

#############################################################
# READ PHP FILE
#############################################################
php="%s/%s.php"%(obsdir,fname)
config=parsePhp(php)

#TIME
hour=config["time"]
ptime=hour.split(":")
time=float(ptime[0])+float(ptime[1])/60.+float(ptime[2])/3600.
print TAB,"Time: ",hour,time

#############################################################
# READ IMAGE
#############################################################
Data=imread("%s/%s"%(obsdir,image))
Mono=Data[:,:,0]
h,w=Mono.shape
X,Y=np.meshgrid(np.arange(w),np.arange(h))
maxval=Mono.max()
print 1*TAB,"Image resolution (w,h): ",w,h

#############################################################
# GET SOLAR BORDER
#############################################################
xthresmin=0;dRR=1e100

#==================================================
#OPTIMAL THRESHOLD FOR MINIMUM RADIUS DISPERSION
#==================================================
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
    
print 1*TAB,"Optimal threshold:"
print 2*TAB,"Optimal threshold = maxval / %.2f"%(xthresmin)
print 2*TAB,"Rimage = %.2f +/- %.2f (%.5f)"%(R,dR,dRR)
print 2*TAB,"Center Image = (%d,%d)"%(xcenter,ycenter)
print 2*TAB,"Is solar disk complete? : ",qcomplete

#==================================================
#CALCULATE CENTER AND RADIUS
#==================================================
xcenter,ycenter,R,dR=findCircleProperties(rs)
rcenter=np.array([xcenter,ycenter])
print TAB*1,"Fit result:"
print TAB*2,"R = %.2f +/- %.2f"%(R,dR)
print TAB*2,"Center = (%d,%d)"%(xcenter,ycenter)

#==================================================
#FIND MERCURY POSITION RESPECT TO CENTER
#==================================================
#POSITION
ppos=config["posmercury"].split(",")
xm=float(ppos[0]);ym=float(ppos[1])
rmerc=np.array([xm,ym])

#DISTANCE
rm=np.sqrt((xm-xcenter)**2+(ym-ycenter)**2)/R

#APPARENT SIZE
Dp=float(ppos[2])
rp=Dp/2.0/R

#APPARENT POSITION ANGLE
AP=np.arctan2((ym-ycenter),(xm-xcenter))*RAD

print TAB,"Mercury position"
print TAB*2,"rmerc = ",rmerc
print TAB*2,"r = %.5f, AP = %.2f deg"%(rm,AP)

#==================================================
#FIND SUNSPOT POSITION RESPECT TO CENTER
#==================================================
if 'Use' not in config["possunspot"]:
    #POSITION
    spos=config["possunspot"].split(",")
    xsp=float(spos[0]);ysp=float(spos[1])
    rspot=np.array([xsp,ysp])

    #DISTANCE
    rsp=np.sqrt((xsp-xcenter)**2+(ysp-ycenter)**2)/R

    #APPARENT POSITION ANGLE OF THE SUNSPOT
    APS=np.arctan2((ysp-ycenter),(xsp-xcenter))*RAD
else:
    rspot=no.array([0,0])
    rsp=0
    APS=0

print TAB,"Sunpot position"
print TAB*2,"rspot = ",rspot
print TAB*2,"r = %.5f, AP = %.2f deg"%(rsp,APS)

#==================================================
#CROP FITTING
#==================================================
debug="%s/scratch/%s-fit-debug.png"%(obsdir,fname)
fig,ax=imgFigure(Data)
ax.plot(rs[::,0],rs[::,1],'g--',ms=1,mec='none')
ax.plot(rcenter[0],rcenter[1],'g+',ms=20)
ax.plot(rmerc[0],rmerc[1],'b+',ms=20)
ax.plot(rspot[0],rspot[1],'r+',ms=20)
fig.savefig(debug)

#############################################################
# CROP
#############################################################
crop="%s/scratch/%s-crop.%s"%(obsdir,fname,ext)

#SIZE OF CROPPED IMAGE
iR=roundFloat(R)
ch=cw=2*iR

#COMPOSITE IMAGE WITH WHITE DISK
White=Image.new('RGBA',(cw,ch),"white")
Crop=Image.new('RGBA',(cw,ch),"black")
draw=ImageDraw.Draw(Crop)
draw.ellipse((0,0)+Crop.size,fill=255)
Crop=Image.composite(Crop,White,Crop)

# Distance from corner to center of crop image
cymin=iR-roundFloat(ycenter)
cxmin=iR-roundFloat(xcenter)

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

dx=0;dy=0
if cxmax>=cw:dx+=cxmax-cw
if cymax>=ch:dy+=cymax-ch
if rxmax>=w:dx+=rxmax-w
if rymax>=h:dy+=rymax-h

rxmax-=dx;rymax-=dy
cxmax-=dx;cymax-=dy

Crop=np.asarray(Crop)
Crop.flags.writeable=True
Crop[cymin:cymax,cxmin:cxmax,:3]=Data[rymin:rymax,rxmin:rxmax,:3]
imsave(crop,Crop)

#==================================================
#RESULTS OF FITTING
#==================================================
#POSITION IN CROPPED IMAGE
rmerc=np.array([cxmin,cymin])+(rmerc-np.array([rxmin,rymin]))
rspot=np.array([cxmin,cymin])+(rspot-np.array([rxmin,rymin]))

#==================================================
#CROP DEBUG
#==================================================
debug="%s/scratch/%s-crop-debug.png"%(obsdir,fname)
fig,ax=imgFigure(Crop)
ax.plot(rmerc[0],rmerc[1],'b+',ms=20)
ax.plot(rspot[0],rspot[1],'r+',ms=20)
fig.savefig(debug)

#############################################################
#SAVE PROPERTIES
#############################################################
jfile="%s/scratch/%s-crop.json"%(obsdir,fname)
properties=dict(
    #Basic information
    hour=hour,
    time=time,
    #Size
    size=[w,h],
    rcenter=[cw/2,ch/2],
    R=R,
    dR=R,
    #Mercury
    rmerc=rmerc.tolist(),
    rm=rm,
    AP=AP,
    rp=rp,
    #Sunspot
    rspot=rspot.tolist(),
    rsp=rsp,
    APS=APS
)
dict2json(properties,jfile)
