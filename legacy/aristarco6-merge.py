from aristarchus import *
#############################################################
#INPUTS
#############################################################
obsdir=argv[1]
images=argv[2]
limages=images.split(",")
nimages=len(limages)

wcommon=1000
hcommon=1000

j=0
alpha=0.4
Alignment=Image.new('RGBA',(wcommon,hcommon),"black")
Background=Image.new('RGBA',(wcommon,hcommon),"black")
for image in limages:
    parts=image.split(".")
    ext=parts[-1]
    fname="".join(parts[:-1])
    Rotated=imread("%s/%s-rotated-result.%s"%(obsdir,fname,ext))
    Rotated=imresize(Rotated,(hcommon,wcommon))
    Alignment=Image.blend(Alignment,Image.fromarray(Rotated),alpha)

alignment="%s/image-alignment-result.%s"%(obsdir,ext)
Alignment=Image.composite(Alignment,Background,Alignment)
Alignment.save(alignment)
