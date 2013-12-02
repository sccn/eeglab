% std_savedat() - save measure for computed data
%
% Usage: std_savedat( filename, structure);
%
% Authors:  Arnaud Delorme, SCCN, INC, UCSD, 2006-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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

function std_savedat( tmpfile, structure)

    delims = find( tmpfile == '.');
    if ~isfield(structure, 'datafile') && ~isfield(structure, 'datafiles') 
        structure.datafile = [ tmpfile(1:delims(end)-1) '.set' ];
    end;
    
    % fix reading problem (bug 764)
    tmpfile2  = which(tmpfile);
    if isempty(tmpfile2), tmpfile2 = tmpfile; end;    
    tmpfile = tmpfile2;
    
    eeglab_options;
    if option_saveversion6
        try
            save('-v6' , tmpfile, '-struct', 'structure');
        catch
            fields = fieldnames(structure);
            for i=1:length(fields)
                eval([ fields{i} '=structure.'  fields{i} ';']);
            end;
            save('-mat', tmpfile, fields{:});
        end;
    else
        save('-v7.3' , tmpfile, '-struct', 'structure');
    end;
