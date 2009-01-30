function [CC]=getclassifier(d,c, Mode)
% GETCLASSIFIER yields the classifier from labeled data
%   CC = getclassifier(d,c)
%   CC = getclassifier(d1,d2)
%
% d     DATA
% c     CLASSLABEL 
% d1    DATA of class 1 
% d2    DATA of class 0
% Mode	'LDA' and 'MDA' implemented. 
%
% number of rows in d and c must fit, and C must be a column vector, 
% OR number of columns in d1 and d2 must fit. 
% The functions COVM.M and SUMSKIPNAN.M from the NaN-toolbox are 
% required [1] for Mode 'LDA' and 'MDA'.
% 
% OUTPUT:
%   CC 	classifier
% 
% see also: LDBC, LLBC, MDBC, NaN/COVM
%
% Reference(s):
% [1] A. Schloegl, Missing values and NaN-toolbox for Matlab, 2000-2003.
% http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/NaN/

%	$Revision: 1.1 $
%	$Id: getclassifier.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%	Copyright (C) 1999-2003 by Alois Schloegl
%	a.schloegl@ieee.org	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if nargin<3,
        Mode = 'LDA';
end;

if exist('covm')~=2,
        fprintf(2,'Error GETCLASSIFIER: COVM.M is not in search path. You need to install the NaN-toolbox./n')
        return;
end;

if (size(d,1)==size(c,1)) & size(c,2)==1 & all(c==round(c));
        
elseif (size(d,2)==size(c,2)) ;
        d1 = d; d2 = c; 
        c  = [ones(size(d1,1),1); ones(size(d2,1),1)*2]; 
        d  = [d1; d2];
else
        fprintf(2,'Error GETCLASSIFIER: incorrect input arguments\n');
        return;
end;        

CL = sort(unique(c(~isnan(c)))); 

if strcmp(Mode,'LDA') | strcmp(Mode,'MDA'),
        for k = 1:length(CL);
                CC{k} = covm(d(c==CL(k),:),'E');        
        end;
else
        fprintf(2,'Error GETCLASSIFIER: classifier %s not implemented.\n');        
        return;
end;

