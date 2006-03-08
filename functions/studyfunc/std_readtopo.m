% std_readtopo() - Given the ALLEEG structure, a specific EEG dataset index, 
% and a specific component, the function returns the scalp map of that ICA component. 
% The scalp map grid of the dataset ICA components is assumed to be saved in a float 
% file, the EEG dataset include a pointer to this file. If such a float file doesn't exist,
% you can use the std_scalp() function to create it, or use the pre - clustering functions
% that call it: pop_preclust, eeg_preclust & eeg_createdata.  
% Along with the scalp map grid of the selected ICA component the function returns  
% the two axis grid points vectors (x and y). 
%
% Usage:    
%   >> [grid, y, x ] = std_readtopo(ALLEEG, abset, component);  
%   This functions returns the ICA component scalp map grid. 
%   The information is loaded from a float file, which a pointer 
%   to is saved in the EEG dataset. The float file was created
%   by the pre - clustering function std_scalp. 
%
% Inputs:
%   ALLEEG     - an EEGLAB data structure, which holds EEG sets (can also be one EEG set). 
%                      ALLEEG must contain the dataset of interest (the setind).
%   setind         -  [integer] an index of an EEG dataset in the ALLEEG
%                      structure, for which to get the component ERP.
%   component - [integer] a component index in the selected EEG dataset for which 
%                      an ERP will be returned. 
%
% Outputs:
%   grid          - the scalp map grid of the requested ICA component in the
%                      selected dataset. This grid is an interpolated Cartesian grid 
%                      of the component scalp map (the output of the topoplot function). 
%   x              - the x axis points of the interpolated grid, for plotting purposes.  
%   y              - the y axis points of the interpolated grid, for plotting purposes.  
%
%  See also  std_scalp, pop_preclust, eeg_preclust, eeg_createdata           
%
% Authors:  Hilit Serby, SCCN, INC, UCSD, February, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $
% Revision 1.2  2006/03/07 22:14:40  arno
% use fullfile
%

function [grid, yi, xi ] = std_readtopo(ALLEEG, abset, comp)

grid = [];
yi = [];
xi = [];
d = ALLEEG(abset).etc.icascalpparams; %the grid dimension 
if iscell(d)
    d = d{1};
end
try
    scalp = floatread( fullfile( ALLEEG(abset).filepath, ALLEEG(abset).etc.icascalp), [d d+2],[],d*(d+2)*(comp-1));
    grid = scalp(:,1:d);
    yi = scalp(:,d+1);
    xi = scalp(:,d+2).';
catch
    warndlg2(['std_readtopo: file '  ALLEEG(abset).etc.icascalp ' was not found in path ' ALLEEG(abset).filepath], 'Abort - computing scalp map centroid' ); 
    return;
end