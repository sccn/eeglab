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

% 07-30-97 converted to rmbase() -sm
% 09-30-97 fixed! -sm
% 05-10-01 caught empty input data -sm
% 06-03-01 added test for length-1 basevector, added [] defaults -sm
% 01-05-02 added check for negative basevector indices -Luca Finelli
% 01-25-02 reformated help & license -ad 

function [dataout,datamean] = rmbase(data,frames,basevector)

	if nargin<3,
		basevector =0;
	end
    if isempty(basevector)
		basevector =0;
	end
    if length(basevector) == 1 && basevector(1) ~= 0
       fprintf('rmbase(): basevector should be a vector of frame indices.\n');
       return
    end

    if sum(basevector<0)~= 0
       fprintf('rmbase(): basevector should be 0 or a vector of positive frame indices.\n');
       return
    end

	if nargin < 2,
		frames = 0;
	end
    if isempty(frames)
		frames =0;
	end
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
	end
    epochs = fix(framestot/frames);

	if length(basevector)>framestot,
		fprintf('rmbase(): length(basevector) > frames per epoch.\n\n');
		help rmbase;
		return
	end

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
            end
            datamean(c,e) = rmeans;
            dataout(c,(e-1)*frames+1:e*frames) = data(c,(e-1)*frames+1:e*frames) - rmeans;
        end
    end

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
