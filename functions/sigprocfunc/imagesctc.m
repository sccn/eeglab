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

if ~isempty(varargin)
	imagesc(a, varargin{:});
else
	imagesc(a);
end;
c = caxis;

cm = colormap;
intervals = linspace(c(1), c(2), size(cm,1));
[tmp dest] = histc(a(:), intervals);

for wi = 1:size(a,2)
	for hi = 1:size(a,1)
		aa(hi,wi,1:3) = cm(dest(hi+(wi-1)*size(a,1)),1:3);
	end;
end;
if ~isempty(varargin)
	imagesc(aa, varargin{:});
else
	imagesc(aa);
end;
