% rmbase() - subtract basevector channel means from multi-epoch data matrix
%
% Usage:
%       >> [dataout] = rmbase(data); % remove whole-data channel means
%       >> [dataout datamean] = rmbase(data,frames,basevector);
%            % remove mean of basevector from each channel and epoch
% Inputs:
%   data       - data matrix (chans,frames*epochs) or (chans, frames, epochs);
%   frames     - data points per epoch {[]|0|default->data length}
%   basevector - vector of baseline frames per epoch
%                 Ex 1:128 {[]|0|default->whole epoch}
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 2-5-96 

% Copyright (C) 2-5-96 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 07-30-97 converted to rmbase() -sm
% 09-30-97 fixed! -sm
% 05-10-01 caught empty input data -sm
% 06-03-01 added test for length-1 basevector, added [] defaults -sm
% 01-05-02 added check for negative basevector indices -Luca Finelli
% 01-25-02 reformated help & license -ad 

function [dataout,datamean] = rmbase(data,frames,basevector)

	if nargin<3,
		basevector =0;
	end;
    if isempty(basevector)
		basevector =0;
	end;
    if length(basevector) == 1 & basevector(1) ~= 0
       fprintf('rmbase(): basevector should be a vector of frame indices.\n');
       return
    end

    if sum(basevector<0)~= 0
       fprintf('rmbase(): basevector should be 0 or a vector of positive frame indices.\n');
       return
    end

	if nargin < 2,
		frames = 0;
	end;
    if isempty(frames)
		frames =0;
	end;
	if nargin<1,
		help rmbase;
		fprintf('rmbase(): needs at least one argument.\n\n');
		return
	end
    if isempty(data)
		fprintf('rmbase(): input data is empty.\n\n');
		return
	end
    
    oridims = size(data);
	if ndims(data) == 3,
		data = reshape(data, size(data,1), size(data,2)*size(data,3));
	    reshape_flag=1;
	end	
	
	[chans framestot]= size(data);
	if frames ==0,
		frames = framestot;
	end;
    epochs = fix(framestot/frames);

	if length(basevector)>framestot,
		fprintf('rmbase(): length(basevector) > frames per epoch.\n\n');
		help rmbase;
		return
	end;

    datamean = zeros(chans,epochs);
    % fprintf('removing epoch means for %d epochs\n',epochs);

    dataout = data;
    for e=1:epochs
        for c=1:chans
            if basevector(1)~=0,
                rmeans = nan_mean(double(data(c,(e-1)*frames+basevector)'));
            else
                rmeans = nan_mean(double(data(c,(e-1)*frames+1:e*frames)'));
                   %if e==1
                   %    fprintf('rmbase(): whole-data channel means removed. \n\n');
                   %end
            end;
            datamean(c,e) = rmeans;
            dataout(c,(e-1)*frames+1:e*frames) = data(c,(e-1)*frames+1:e*frames) - rmeans;
        end;
    end;

    dataout = reshape(dataout, oridims);
    
    
% function out = nan_mean(in)
%     
%     nans = find(isnan(in));
%     in(nans) = 0;
%     sums = sum(in);
%     nonnans = ones(size(in));
%     nonnans(nans) = 0;
%     nonnans = sum(nonnans);
%     nononnans = find(nonnans==0);
%     nonnans(nononnans) = 1;
%     out = sum(in)./nonnans;
%     out(nononnans) = NaN;
% 
