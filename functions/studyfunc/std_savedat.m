% std_savedat() - save measure for computed data
%
% Usage: std_savedat( filename, structure);
%
% Authors:  Arnaud Delorme, SCCN, INC, UCSD, 2006-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function std_savedat( tmpfile, structure)

    delims = find( tmpfile == '.');
    if ~isfield(structure, 'datafile') && ~isfield(structure, 'datafiles') 
        structure.datafile = [ tmpfile(1:delims(end)-1) '.set' ];
    end
    
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
            end
            save('-mat', tmpfile, fields{:});
        end
    else
        save('-v7.3' , tmpfile, '-struct', 'structure');
    end
