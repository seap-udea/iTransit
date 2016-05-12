from aristarchus import *
#############################################################
#INPUTS
#############################################################
img=argv[1]
coords=argv[2]

#############################################################
#READ IMAGE
#############################################################
data=imread(img)
image=data[:,:,1]
h,w=image.shape

#############################################################
#GET SUBIMAGE
#############################################################
coords=eval("np.array(["+coords+"])")
coords[0::2]*=w
coords[1::2]*=h
xmin=min(coords[1],coords[3]);xmax=max(coords[1],coords[3])
ymin=min(coords[0],coords[2]);ymax=max(coords[0],coords[2])
subimage=image[xmin:xmax,ymin:ymax]

sh,sw=subimage.shape
if sw<2 and sh<2:
    print "%d,%d,%d,%d"%(roundFloat(coords[0]),roundFloat(coords[1]),0,0)

#############################################################
#FIND SPOT PIXELS
#############################################################
maxval=subimage.max()
meanval=subimage.mean()
stdval=subimage.std()
js=np.arange(sw)
xs=[];ys=[]

threshold=0.7*maxval
ts=[]
values=[]
for k in xrange(1):
    ts+=[threshold]
    nex=0
    minex=1e100
    for i in xrange(sh):
        line=subimage[i,:]
        condin=line<threshold
        condex=line>=threshold

        limits=js[condin]
        if len(limits)>0:
            ys+=[i]*len(limits)
            xs+=limits.tolist()

        exterior=js[condex]
        if len(exterior):
            minline=subimage[i,exterior].mean()
            values+=[subimage[i,exterior].tolist()]
            if minline<minex:minex=minline
    threshold=minex

xs=np.array(xs)
ys=np.array(ys)
if len(xs)==0 or len(ys)==0:
    print "No object found"
    exit(1)

#############################################################
#CENTROID
#############################################################
xmean=np.mean(xs)
ymean=np.mean(ys)
xext=max(xs)-min(xs)
yext=max(ys)-min(ys)

xmerc=xmean+ymin
ymerc=ymean+xmin

#############################################################
#OUTPUT
#############################################################
print "%d,%d,%.2f"%(roundFloat(xmerc),roundFloat(ymerc),np.mean([xext,yext]))
