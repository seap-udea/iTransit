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
rpos=np.unravel_index(subimage.argmin(),subimage.shape)

subimage=image[rpos[0]-5:rpos[0]+5,rpos[1]-5:rpos[1]+5]
rpos=np.unravel_index(subimage.argmin(),subimage.shape)

params=dict(S=subimage)
x=[rpos[0],rpos[1],1]
solution=minimize(boxFig,x,args=(params,))

print solution

r=solution["x"]
rpos=np.array([r[0],r[1]])+np.array([xmin,ymin])
rpos+=np.array([xmin,ymin])

R=r[2]
print "%.0f,%.0f,%.2f"%(rpos[0],rpos[1],R)

