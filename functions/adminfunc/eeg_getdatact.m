% eeg_getdatact() - get EEG data from a specified dataset or
%                  component activity
%
% Usage:
%       >> signal = eeg_getdatact( EEG );
%       >> signal = eeg_getdatact( EEG, 'key', 'val');
%
% Inputs:
%   EEG       - Input dataset
%
% Optional input:
%   'channel'   - [integer array] read only specific channels.
%                 Default is to read all data channels.
%   'component' - [integer array] read only specific components
%
% Outputs:
%   signal      - EEG data or component activity
%
% Author: Arnaud Delorme, SCCN & CERCO, CNRS, 2008-
%
% See also: eeg_checkset()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 15 Feb 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function data = eeg_getdatact( EEG, varargin);
    
    data = [];
    if nargin < 1
        help eeg_getdatact;
        return;
    end;
    
    opt = finputcheck(varargin, { 'channel'   'integer' {} [1:EEG.nbchan];
                                  'component' 'integer' {} [] }, 'eeg_getdatact');
    if isstr(opt), error(opt); end;
    
    if strcmpi(EEG.data, 'in set file')
        filename = fullfile(EEG.filepath, EEG.filename);
        EEG = pop_loadset(filename);
    end;        
    
    if isnumeric(EEG.data)
        if ~isempty(opt.component)
            if isempty(EEG.icaact)
                data = (EEG.icaweights(opt.component,:)*EEG.icasphere)*EEG.data;
            else
                data = EEG.icaact(opt.component,:,:);
            end;
        else
            data = EEG.data(opt.channel,:,:);
        end;
    else
        % opening data file
        % -----------------
        if ~isempty(opt.component)
            filename = fullfile(EEG.filepath, [ EEG.data(1:end-3) 'icaact' ] );
            
            % reading ICA file
            % ----------------
            if exist(filename)
                data = repmat(single(0), [ length(opt.component) EEG.pnts EEG.trials ]);
                fid = fopen( filename, 'r', 'ieee-le'); %little endian (see also pop_saveset)
                if fid == -1, error( ['file ' filename ' could not be open' ]); end;
                for ind = 1:length(opt.component)
                    fseek(fid, (opt.component(ind)-1)*EEG.pnts*EEG.trials, -1);
                    data(ind,:) = fread(fid, [EEG.trials*EEG.pnts 1], 'float32')';
                end;
                fclose(fid);
                data = reshape(data, size(data,1), EEG.pnts, EEG.trials);
                return;
            end;
            
        end;
        
        filename = fullfile(EEG.filepath, EEG.data);
        fid = fopen( filename, 'r', 'ieee-le'); %little endian (see also pop_saveset)
        if fid == -1
            error( ['file ' filename ' not found. If you have renamed/moved' 10 ...
                    'the .set file, you must also rename/move the associated data file.' ]);
        else 
            fprintf('Reading float file ''%s''...\n', filename);
        end;
        
        % old format = .fdt; new format = .dat (transposed)
        % -------------------------------------------------
        datformat = 0;
        if length(filename) > 3
            if strcmpi(filename(end-2:end), 'dat')
                datformat = 1;
            end;
        end;
        EEG.datfile = EEG.data;
        
        % reading data file
        % -----------------
        if datformat
            if length(opt.channel) == EEG.nbchan
                data = fread(fid, [EEG.trials*EEG.pnts EEG.nbchan], 'float32')';
            else
                data = repmat(single(0), [ length(opt.channel) EEG.pnts EEG.trials ]);
                for ind = 1:length(opt.channel)
                    fseek(fid, (opt.channel(ind)-1)*EEG.pnts*EEG.trials, -1);
                    data(ind,:) = fread(fid, [EEG.trials*EEG.pnts 1], 'float32')';
                end;
            end;
        else
            data = fread(fid, [EEG.nbchan Inf], 'float32');
            data = data(opt.channel,:,:);
        end;
        fclose(fid);
        
        % recompute component activity
        % ----------------------------
        if ~isempty(opt.component)
            data = (EEG.icaweights(opt.component,:)*EEG.icasphere)*data;
        end;
        
        data = reshape(data, size(data,1), EEG.pnts, EEG.trials);
    end;
 
