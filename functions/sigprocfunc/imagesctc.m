% imagesctc() - DEPRECATED. never completed or documented.
% imagesctc() - imagesc in true color. Can help plot different
%               colormap on the same window.
%
% Usage: same as imagesctc
%
% Exemple:
%         figure; 
%         colormap(jet); 
%         subplot(1,2,1); imagesctc(rand(10,10));
%         colormap(gray); 
%         subplot(1,2,2); imagesctc(rand(10,10));
% 
% Author: Arnaud Delorme, CNL / Salk Institute, 31 July 2002
%
% See also: image(), imagesc()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

if ~isempty(varargin)
	imagesc(a, varargin{:});
else
	imagesc(a);
end
c = caxis;

cm = colormap;
intervals = linspace(c(1), c(2), size(cm,1));
[tmp dest] = histc(a(:), intervals);

for wi = 1:size(a,2)
	for hi = 1:size(a,1)
		aa(hi,wi,1:3) = cm(dest(hi+(wi-1)*size(a,1)),1:3);
	end
end
if ~isempty(varargin)
	imagesc(aa, varargin{:});
else
	imagesc(aa);
end
