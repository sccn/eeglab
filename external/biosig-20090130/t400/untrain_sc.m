function [CC]=untrain_sc(CC,classlabel,x)
% UnTrain - decrementaal learning (untraining) of classifier 
% 
%  CC = untrain_sc(CC,classlabel,x)
%       
% CC is a statistical classifier, it contains basically the mean 
% and the covariance of the data of each class. This information 
% is incoded in the so-called "extended covariance matrices".  
%
% CC can be used for various statistical classifiers included
%  LDA, MDA, QDA, GRB, etc. 
%
% see also: TEST_SC, COVM, LDBC2, LDBC3, LDBC4, MDBC, GDBC

%	$Id: untrain_sc.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%	Copyright (C) 2005 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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


if (length(classlabel)~=size(x,1)) & any(size(classlabel)~=1),          
        error('length of classlabel must be a scalar or must fit size of data (number of rows)');
end;

if strcmp(CC.datatype,'classifier:statistical');
        if all(classlabel==classlabel(1)),
                [md,nn] = covm(x,'E');
                k = find(CC.Labels==classlabel(1));
                CC.MD(k,:,:) = CC.MD(k,:,:) - md;
                CC.NN(k,:,:) = CC.NN(k,:,:) - nn;
        else
                Labels = unique(classlabel(~isnan(classlabel)));
                for k = 1:length(Labels),
                        [md,nn] = covm(x(classlabel==Labels(k),:),'E');
                        
                        ix = find(CC.Labels==Labels(k));
                        CC.MD(ix,:,:) = CC.MD(ix,:,:) - md;
                        CC.NN(ix,:,:) = CC.NN(ix,:,:) - nn;
                end
        end;

elseif strcmp(CC.datatype,'classifier:SVM');
        error('decremental learning not implemented for SVM (yet)');         
        
end;