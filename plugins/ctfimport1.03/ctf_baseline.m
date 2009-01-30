function [ ctf ] = ctf_baseline(ctf)

% ctf_baseline - remove mean of pretrigger period from all ctf.data
%
% [ ctf ] = ctf_baseline(ctf)
%


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $

% Copyright (C) 2004  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Created: 01/2004, copyright 2004 Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for epoch = 1:size(ctf.data,3),
  
  pretrigger_mean = mean(ctf.data(1:ctf.setup.pretrigger_samples,:,epoch));
  
  pretrigger_mean = repmat(pretrigger_mean,ctf.setup.number_samples,1);
  
  ctf.data(:,:,epoch) = ctf.data(:,:,epoch) - pretrigger_mean;
  
end
