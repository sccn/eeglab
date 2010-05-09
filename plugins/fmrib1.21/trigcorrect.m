% trigcorrect() -   correct for missing FMRI slice triggers
%
% Usage: >>trigdata=trigcorrect(trigdata,slices,volumes,trigch)
%
% Input:
%   trigdata:data containing trigger channel.  triggers in this
%   channel should have the minimum value anywhere in the channel.
%   slices: number of FMRI slices / vol.
%   volumes: number of FMRI vols.
%
% Output:
%   trigdata:data with trigger channel corrected.
%
% Author: Rami Niazy, FMRIB Centre, University of Oxford.
%
% 
%
% Copyright (c) 2004 University of Oxford.
%


% Copyright (C) 2004 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
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
%
%

%    DEC 23, 2004
%    Update (c)

%     01 OCT, 2004.


function T=trigcorrect(T,s,v,trigch)

trigval=min(T(trigch,:));
atrigs=find(T(trigch,:)==trigval);
atrig_l=length(atrigs);
dt=diff(atrigs);


ctrig_l=s*v;
trig_err=ctrig_l-atrig_l;

if trig_err==0
    return
end

vols_per_sec=5;
trigs_per_sec=s*vols_per_sec;
sections=floor((atrig_l-1)/trigs_per_sec);

dtrigmat=reshape(dt(1:sections*trigs_per_sec),trigs_per_sec,sections)';
maxdt=min(max(dtrigmat,[],2));
mindt=round(mean(min(dtrigmat,[],2)));
errdt=max(max(dtrigmat,[],2));
mxth=round((maxdt+mindt)/2);
errth=round((errdt+mindt)/2);
errIn=find(dt>errth);
errsz=length(errIn);
voltrigs=find(dt>=maxdt & dt<=(maxdt+2));


if trig_err~=errsz
    error('can not determine trigger error')
end

for er=1:errsz
    postV=find((voltrigs-errIn(er))>0);
    if (voltrigs(postV(1))-voltrigs(postV(1)-1))>s
        trigint=maxdt;
    else
        trigint=mindt;
    end
    T(1,atrigs(errIn(er))+trigint)=trigval;
end
