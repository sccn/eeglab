EEGLAB channel location formats are 
- ".loc" for polar coordinates
- ".xyz" for 3-D cartesian coordinates
- ".sph" for 3-D spherical coordinates
- ".txt" for all coordinate formats
(see help readlocs under Matlab for a description of these formats)

Upon reading a file, all coordinate frames are converted into one 
another (polar, 3-D spherical, 3-D cartesian). Polar coordinates 
can for instance be converted to 3-D spherical (assuming radius one), 
then to 3-D cartesian. 

