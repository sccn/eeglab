EEGLAB channel location formats are 
- ".loc" for polar coordinates
- ".xyz" for 3-D carthesian coordinates
- ".sph" for 3-D spherical coordinates
- ".txt" for all coordinate formats
(see help readlocs under Matlab for a description of these formats)

Upon reading a file, all coordinate frames are converted into one 
another (polar, 3-D spherical, 3-D carthesian). Polar coordinates 
can for instance be converted to 3-D spherical (asuming radius one), 
then to 3-D carthesian. 

