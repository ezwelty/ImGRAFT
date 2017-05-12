This is a highly modified version of [ImGRAFT](https://github.com/grinsted/ImGRAFT), the "Image GeoRectification And Feature Tracking" toolbox for MATLAB by Alexandra Messerli and Aslak Grinsted:

> Alexandra Messerli and Aslak Grinsted (2015), Image GeoRectification And Feature Tracking toolbox: ImGRAFT, Geoscientific Instrumentation, Methods and Data Systems, 4, 23-34, doi:[10.5194/gi-4-23-2015](http://dx.doi.org/10.5194/gi-4-23-2015).

New features include:

* Class definitions for image and digital elevation model (DEM) objects
* Faster and more precise undistort routine
* Optimize over multiple images at once
* Optimize using "ground control lines"
* Solve for relative camera rotation without a DEM
* Intersect camera rays with a DEM using a discrete ray traversal algorithm

Be warned, however, that the code currently has some toolbox dependencies, which the original toolbox does not.
