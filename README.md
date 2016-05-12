iTransit
========

This package of python sripts is intended to analyze images of a
planetary transit (Venus or Mercury transit).

We start with a set of N images of the same transit.  Images could
have different resolutions and could correspond to different regions
of the solar disk.  It is mandatory however that the images include at
least a fraction of the solar disk.

Step 1. Enclose position of the planet and a sunspot
----------------------------------------------------

The first step is try to locate the planet and at least one sunspot in
the image.  For that purpose we need to have an initial guess of the
region in the image where the planet or sunspot lie.

The region should be described in the form: 

    ymin,xmin,ymax,xmax

Hereafter 'x' and 'y' correspond to the column and row in the image
(with x=0,y=0 the upper-left corner).

Step 2. Locate the planet and sunspot
-------------------------------------

Run the command:

    $ python itransit-locate.py <image> ymin,xmin,ymax,xmax

The script will return 3 numbers:

    * The abcisa (x)
    * The ordinate (y)
    * The radius of the spot (r)

Once the coordinates of Mercury and the sunspot have been determined,
fill a configuration file in php with the following information:

     <?php
     $time='09:44:38';
     $mercury='0.5033,0.359,0.6,0.4066';
     $sunspot='0.4483,0.4557,0.56,0.5164';
     $posmercury='289,203,2.00';
     $possunspot='260,254,4.50';
     ?>

You need to repeat this procedure for all the images in your set.

Step 3. Crop the images
-----------------------

For each image run:

    $ python itransit-crop.py <image>.<ext>

The script:

    * Find the border of the solar disk.

    * Find the center of the solar disk.

    * Create a crop image with a white solar disk and overimposed the
      fraction of the actual disk in the image.

The script will create a folder 'scratch' including the following files:

    * <image>-fit-debug.png: superposition of solar disk with the
      coordinates of the border, center and the position of the planet
      and the sunspot.

    * <image>-crop.png: cropped image.

    * <image>-crop-debug.png: cropped image with information
      superimposed.


    * <image>.json: information about the image:

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

You will need to repeat this for all images in an observation.

    $ for img in *.png;do python itransit-crop.py $img;done

Step 4. Align several images
----------------------------

Once you have cropped several images you can align them into a single one.

    $ python itransit-alignment.py <typealignment> <targetdirectory> <image1> <image2> <image3> ...

Where:

    * Typealignment:
    
      - auto: the alignment of the images is performed with a
        particular algorithm devised by the Aristarchus Campaigns
        technical team.

      - spot: the alignment of the images is achieved using as a
        reference one sunspot visible during the transit.

Inputs:
	
    * Crop images.
    * Json configurations of the crop images.

The script attempt to align the images using the two methods mentioned
above.

The output of the alignment process is:

    * alignment-plot.png: alignment plot, i.e. distance between points
      vs. time.  In the ideal case of a prefect alignment the plot
      would be a straight line with slope equal to the angular
      velocity of the planet.

    * <image>-rotated.png: rotated images. A version of the image
      rotated in such a way that all images align following the
      algorithm chosen.

    * aligned.png: Put rotated images together. Composition between images
      through a minimum filter.  Each pixel of each image is compared
      and only the minimum value among all images is the final one.

    * aligned-debug.png: Overlay information about the alignment, position of Mercury,
      sunspot and the path of the planet across the solar disk.

    * aligned-final.png: Final aligned image, i.e. an image where the
      path is horizontal.


    * aligned-final-debug.png.  The same as aligned-debug.png but for
      the result of the final alignment.

    * aligned.json. Configuration of the final alignment.

        properties=dict(
            
            #Images
            imagefiles: list with image files aligned.
            typealignment: type of alignment.
            size: size (width, height) of the aligned image.
            rcen: center of the aligned image.
            R: radius of the final image.
            times: times of the images.
            
            #Mercury positions
            rmercs: position of mercury in all images in the reference frame of the aligned image.
            
            #Sunspot positions
            rspots: position of spots in all images in the reference frame of the aligned image.
            
            #Coord equation coefficients
            pathcoef: slope m and intercept b
        )

      List are sorted according to time.

Example
-------

We provide with the package a set of images of the transit of Mercury
of 2006.  The images has been especially prepared to combine images
with different orientations, fractions of the solar disk and
resolutions.  The position of Mercury and a sunspot has been already
determined.

To perform the alignment run:

   $ for img in examples/*.png;do python itransit-crop.py $img;done

   $ python itransit-alignment.py auto examples examples/*.png

