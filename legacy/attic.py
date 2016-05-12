    #Border of the "white" area
    contour=plt.contour(cond,levels=[0.5])
    border=[]
    R=w/2
    for path in contour.collections[0].get_paths():
        for points in path.vertices:
            print points
            if abs(points[0]-xm)>R or abs(points[1]-ym)>R:
                continue
            border+=[points.tolist()]
    border=np.array(border)
    
    plt.figure()
    plt.plot(border[:,0],border[:,1])
    plt.savefig("tmp/c.png")

    #Border of the "white" area
    contour=plt.contour(cond,levels=[0.5])
    border=[]
    R=w/2
    for path in contour.collections[0].get_paths():
        for points in path.vertices:
            if abs(points[0]-xm)>R or abs(points[1]-ym)>R:
                continue
            border+=[points.tolist()]
    border=np.array(border)
    plt.figure(figsize=(6,6))
    plt.plot(border[:,0],border[:,1],'ko',ms=0.1)
    ext=max(w,h)
    plt.xlim((0,ext))
    plt.ylim((0,ext))
    plt.savefig("tmp/c.png")
    break


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #CROP IMAGE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #"White" area
    cond=Mono>=maxval/xthresmin
    xs=X[cond]
    ys=Y[cond]
    
    #Mean-Min-Max
    xm=xs.mean();ym=ys.mean()
    xmin=xs.min();xmax=xs.max()
    ymin=ys.min();ymax=ys.max()
    print xmin,xmax,ymin,ymax
    exit(0)
    

nex=0
meanex=0
for i in xrange(sh):
    line=subimage[i,:]
    condin=line<(meanval-2*stdval)
    condex=line>=(meanval-2*stdval)

    limits=js[condin]
    if len(limits)>0:
        ys+=[i]*len(limits)
        xs+=limits.tolist()

    exterior=js[condex]
    if len(exterior):
        meanex+=subimage[i,exterior].sum()
        nex+=len(exterior)
meanex/=(1.0*nex)

147,242,0,0 [[139, 138, 146], [138, 136, 153], [131, 76, 107], [102, 72], [118, 104, 121], [130, 137, 139], [139, 138, 146], [138, 136, 153], [131, 76, 107], [102, 72], [118, 104, 121], [130, 137, 139], [139, 138, 146], [138, 136, 153], [131, 76, 107], [102, 72], [118, 104, 121], [130, 137, 139]] [54.454260679724818, 72, 72] 72 153 85.8937970065 107.1 117.333333333 [ 147.055  147.055  147.055] [ 242.1538  242.1538  242.1538]



        rsort=sortAbscisas(rs,xcenter,ycenter)
        rstart=rsort[0]
        rend=rsort[-1]
        
        # GENERATE ALL POINTS
        minstart=1e100
        minend=1e100
        for q in np.linspace(0,2*np.pi,1000):
            x=xcenter+R*np.cos(q)
            y=ycenter+R*np.sin(q)
            dstart=np.sqrt((x-rstart[0])**2+(y-rstart[1])**2)
            dend=np.sqrt((x-rend[0])**2+(y-rend[1])**2)
            if dstart<=minstart:
                xstart=[x,y]
                minstart=dstart
            if dend<=minend:
                xend=[x,y]
                minstart=dstart

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

        print>>stderr,"\t"*1,"After recalculation:"
        print>>stderr,"\t"*2,"R = %.2f +/- %.2f (%.5f)"%(R,dR,dRR)
        print>>stderr,"\t"*2,"Center = (%d,%d)"%(xcenter,ycenter)
        plt.plot([x1],[y1],'bs',ms=10)
        plt.plot([x2],[y2],'rs',ms=10)
        plt.plot([xs],[ys],'gs',ms=10)
        plt.plot([x1,x2],[y1,y2],'k-')
        plt.axhline(ycenter)
        plt.axvline(xcenter)

        print>>stderr,"\t"*1,"After recalculation:"
        print>>stderr,"\t"*2,"R = %.2f +/- %.2f"%(R,dR)
        print>>stderr,"\t"*2,"Center = (%d,%d)"%(xcenter,ycenter)



        #========================================
        #SOLAR DISK IS COMPLETE
        #========================================
        print>>stderr,"\t","Cropping a complete image"

        #Minimum and maximum determined by radius
        xmin=int(np.round(xcenter-R))-dp;xmax=int(np.round(xcenter+R))+dp
        ymin=int(np.round(ycenter-R))-dp;ymax=int(np.round(ycenter+R))+dp
        xmin,xmincorner=cropCoord(xmin,w)
        ymin,ymincorner=cropCoord(ymin,h)
        xmax,xmaxcorner=cropCoord(xmax,w)
        ymax,ymaxcorner=cropCoord(ymax,h)

        #Crop data
        Crop=Data[ymin:ymax,xmin:xmax,:]
        plt.imsave(crop,Crop)
        ch,cw=Crop[:,:,0].shape

        #Center of cropped image
        xcropcenter=roundFloat(xcenter-xmin)
        ycropcenter=roundFloat(ycenter-ymin)

        #Option for rotation
        reshape=False
