% varimax() - Perform orthogonal Varimax rotation on rows of a data 
%             matrix.
%
% Usage: >> V = varimax(data); 
%        >> [V,rotdata] = varimax(data,tol);
%        >> [V,rotdata] = varimax(data,tol,'noreorder') 
%
% Inputs:
%   data        - data matrix
%   tol         - set the termination tolerance to tol {default: 1e-4}
%   'noreorder' - Perform the rotation without component reorientation 
%                 or reordering by size. This suppression is desirable 
%                 when doing a q-mode analysis. {default|0|[] -> reorder}
% Outputs:
%   V           - orthogonal rotation matrix, hence
%   rotdata     - rotated matrix, rotdata = V*data;
%
% Author: Sigurd Enghoff - CNL / Salk Institute, La Jolla 6/18/98
%
% See also: runica(), pcasvd(), promax()

% Copyright (C) Sigurd Enghoff - CNL / Salk Institute, La Jolla 6/18/98
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

% Reference: % Henry F. Kaiser (1958) The Varimx criterion for 
% analytic rotation in factor analysis. Pychometrika 23:187-200.
%
% modified to return V alone by Scott Makeig, 6/23/98
% 01-25-02 reformated help & license, added link -ad 

function [V,data] = varimax(data,tol,reorder)

if nargin < 1
   help varimax
   return
end

DEFAULT_TOL = 1e-4;  % default tolerance, for use in stopping the iteration
DEFAULT_REORDER = 1; % default to reordering the output rows by size
                     % and adjusting their sign to be rms positive.
MAX_ITERATIONS = 50; % Default
qrtr = .25;          % fixed value

if nargin < 3
    reorder = DEFAULT_REORDER;
elseif isempty(reorder) || reorder(1) == 0
    reorder = 1; % set default
else
    reorder = strcmp('reorder',reorder);
end

if nargin < 2
    tol = 0;
end
if tol == 0
   tol = DEFAULT_TOL;
end
if ischar(tol)
   fprintf('varimax(): tol must be a number > 0\n');
   help varimax
   return
end

eps1 = tol; % varimax toler
eps2 = tol;

V = eye(size(data,1)); % do unto 'V' what is done to data
crit = [sum(sum(data'.^4)-sum(data'.^2).^2/size(data,2)) 0];
inoim = 0;
iflip = 1;
ict = 0;

fprintf(...
  'Finding the orthogonal Varimax rotation using delta tolerance %d...\n',...
                                                                eps1);
while inoim < 2 & ict < MAX_ITERATIONS & iflip,
    iflip = 0;
    for j = 1:size(data,1)-1,
        for k = j+1:size(data,1),
            u = data(j,:).^2-data(k,:).^2;
            v = 2*data(j,:).*data(k,:);
            a = sum(u);
            b = sum(v);
            c = sum(u.^2-v.^2);
            d = sum(u.*v);

            fden = size(data,2)*c + b^2 - a^2;
            fnum = 2 * (size(data,2)*d - a*b);

            if abs(fnum) > eps1*abs(fden)
                iflip = 1;
                angl = qrtr*atan2(fnum,fden);
                tmp    =  cos(angl)*V(j,:)+sin(angl)*V(k,:);
                V(k,:) = -sin(angl)*V(j,:)+cos(angl)*V(k,:);
                V(j,:) = tmp;

                tmp       =  cos(angl)*data(j,:)+sin(angl)*data(k,:);
                data(k,:) = -sin(angl)*data(j,:)+cos(angl)*data(k,:);
                data(j,:) = tmp;
            end
        end
    end

    crit = [sum(sum(data'.^4)-sum(data'.^2).^2/size(data,2)) crit(1)];
    inoim = inoim + 1;
    ict = ict + 1;

    fprintf('#%d - delta = %g\n',ict,(crit(1)-crit(2))/crit(1));

    if (crit(1) - crit(2)) / crit(1) > eps2
        inoim = 0;
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if reorder
    fprintf('Reordering rows...');
    [fnorm index] = sort(sum(data'.^2));
    V = V .* ((2 * (sum(data') > 0) - 1)' * ones(1, size(V,2)));
    data = data .* ((2 * (sum(data') > 0) - 1)' * ones(1, size(data,2)));
    V = V(fliplr(index),:);
    data = data(fliplr(index),:);
    fprintf('\n');
else
    fprintf('Not reordering rows.\n');
end
