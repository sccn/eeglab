% plugin_pathlast()   - This function is meant to be applied to functions 
%                       that overshadow matlab built-in functions. Given a
%                       function  and a plugin option, this function will check
%                       if the given function belong to the plugin selected.
%                       If so, the function will be moved down to the
%                       overshadowed function.
%
% Usage: >> plugin_pathlast('sopen', 'biosig');
% Inputs:
%   functocheck         - Function to be checked and moved down in the path
%   pluginkeyfunction   - {'biosig', 'fieldtrip'}
% 
% Optional input:
%  <subfolder path>     - If the function 'functocheck' is located in a
%                         subfolder in the plugin folder, then the relative 
%                         path to it should be provided.
%
% Author: Arnaud Delorme
%         Ramon Martinez-Cancino
%
% See also: findduplicatefunctions()
%
% Copyright (C) 2017 Arnaud Delorme, SCCN
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

function plugin_pathlast(functocheck,pluginkeyfunction,varargin)

% List of key functions to easy locate the plugin installed
switch pluginkeyfunction
    case 'biosig'
        keyfunction = 'sopen';
    case 'fieldtrip'
        keyfunction = 'ft_wizard';
end

functocheckpath = fileparts( which(functocheck) );
keyfunctionpath      = fileparts( which(keyfunction) );

% In case the function is located in one of the plugin's subfolder
if ~isempty(varargin)
    keyfunctionpath = fullfile(keyfunctionpath,varargin{1});
end
    
if strcmp(functocheckpath,keyfunctionpath)
    warning('off','MATLAB:dispatcher:nameConflict');
    rmpath(functocheckpath);
    str2doublepat2 = fileparts( which(functocheck) );
    addpath(functocheckpath,'-begin');
    addpath(str2doublepat2,'-begin');
    warning('on','MATLAB:dispatcher:nameConflict');
end