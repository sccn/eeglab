% mi_pairs() - Calculate pairwise mutual information between rows of an input data matrix.
%
% Usage:
%              >> MI = mi_pairs (udata);
%              >> MI = mi_pairs (udata,returntype,verbose)
% Inputs:
%
%   udata      - (N,M) data matrix containing M samples for each of N components 
%
% Optional Inputs:
%
%   returntype - ('full' or 'diagonal') specify which portion of the elements 
%                of the udata matrix should be calculated {default: 'full'}
%   verbose    - (true or false) control how much information should be 
%                displayed during calculation {default: true}
% Output:
%
%   MI         - NxN mutual information matrix 
%
% Example:
%
%   % Compute mutual information matrix between ICA component activation pairs
%   >> MI = mi_pairs (EEG.icaact);
%
% See also: eeg_miclust(), getclusts(), showclusts(),

% Copyright (C) 2006 Jason Palmer and Nima Bigdely Shamlo, SCCN/INC/UCSD, nima@sccn.ucsd.edu
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

function MI = mi_pairs (u,returnType,verbose)

if nargin < 2 
   returnType = 'full'; 
end;

if nargin < 3
    verbose = true;
end;

if strcmp(returnType,'diagonal')
  fprintf('Only the diagonal elements will be calculated.\n');
elseif ~strcmp(returnType,'full')
  error('Unknown returntype argument');
end;

[m,N] = size(u);
ub = uint8(zeros(m,N));
nbins = 100;
H = zeros(m,m);
MI = zeros(m,m);
b = zeros(m,nbins);

%
% set up bins
%
for i = 1:m
    mx = max(u(i,:));
    mn = min(u(i,:));
    delta(i) = (mx-mn)/nbins;
    mx = mx + delta(i)/2; mn = mn - delta(i)/2;
    delta(i) = (mx-mn)/nbins;
    b(i,:) = mn + delta(i) * (1:nbins);
end

%
% Compute binned time series and marginal histograms
%
if verbose 
  disp('Binning the time series and computing marginal entropies ...'); pause(0.2);
end;
pdf = zeros(m,nbins);
for i = 1:m
    for r = 1:N
        if true% Nima
            for k = 1:nbins
                if u(i,r) <= b(i,k)
                    ub(i,r) = k;
                    pdf(i,k) = pdf(i,k) + 1;
                    break;
                end
            end
        end;
    end
    pdf(i,:) = pdf(i,:) / (sum(pdf(i,:))*delta(i));
end

%
% Compute marginal entropies
%
if verbose 
  disp('Computing marginal entropies ...'); pause(0.2);
end;
    for i = 1:m
        for k = 1:nbins
            if pdf(i,k) ~= 0    
                H(i,i) = H(i,i) - pdf(i,k)*log(pdf(i,k))*delta(i);
            end
        end
    end

if verbose 
  h = waitbar(0,'please wait...');
end;

%
% Compute pairwise histograms
%
if verbose 
  disp('getting pairwise entropies ...'); pause(0.2);
end;
for i = 1:m
    if verbose 
        waitbar(i/m);
    end;
    for j = (i+1):m
        if strcmp(returnType,'full') | strcmp(returnType,'diagonal')  & (i==j) 
            pdf2 = zeros(nbins,nbins);
            for r = 1:N
                bi = ub(i,r);
                bj = ub(j,r);
                
                pdf2(bi,bj) = pdf2(bi,bj) + 1;
            end
            di = delta(i);
            dj = delta(j);
            pdf2 = pdf2 / (N*di*dj);
            for k = 1:nbins
                for l = 1:nbins
                    if pdf2(k,l) ~= 0
                        H(i,j) = H(i,j) - pdf2(k,l)*log(pdf2(k,l))*di*dj;
                    end   
                end
            end
            MI(i,j) = H(i,i) + H(j,j) - H(i,j);
            MI(j,i) = MI(i,j);
         end;
    end        
    MI(i,i) = H(i,i);
end
if verbose 
  close(h);
end;

