Aristarchus Campaign
********************

Mercury Transit Image Alignment
===============================

Using this set of scripts you will be able to perform a full set of
analysis on the images of a Mercury or Venus Transit.

We start with a set of N images of the same transit.  Images could
have different resolutions and could correspond to different regions
of the solar disk.

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

    $ python aristarco6-locate.py <image> ymin,xmin,ymax,xmax

