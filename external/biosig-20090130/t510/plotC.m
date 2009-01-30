function plotC(input_struct)
% Creates different plots such as ERDS maps, average maps, and so on.
%
% This function creates different plots, depending on the content of the input
% parameter. It is a wrapper function for the specific plot functions in order
% to allow the user to call the same plot function for all kinds of different
% plots (ERDS maps, average maps, spectra, heart rate, ...).
%
% Usage:
%   plotC(input_struct);
%
% Input parameters:
%   input_struct ... Input structure (as obtained by the calculation functions)

% Copyright by Clemens Brunner
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:51 $
% E-Mail: clemens.brunner@tugraz.at

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if (strcmp(input_struct{1}.datatype, 'ERDS map'))
    plotErdsMap(input_struct);
end;

if (strcmp(input_struct{1}.datatype, 'Average map'))
    plotAveMap(input_struct);
end;

if (strcmp(input_struct{1}.datatype, 'Raw trials'))
    plotRawTrials(input_struct);
end;

if (strcmp(input_struct{1}.datatype, 'Spectrum'))
    plotSpectrum(input_struct);
end;

if (strcmp(input_struct{1}.datatype, 'Heart rate'))
    plotHR(input_struct);
end;

if (strcmp(input_struct{1}.datatype, 'Heart rate relative'))
    plotHRr(input_struct);
end;

if (strcmp(input_struct{1}.datatype, 'Heart rate continuous'))
    plotHRwhole(input_struct);
end;

if (strcmp(input_struct{1}.datatype, 'Heart rate spectrum'))
    plotHRSpectrum(input_struct);
end;

if (strcmp(input_struct{1}.datatype, 'PLV'))
    plotPLV(input_struct);
end;