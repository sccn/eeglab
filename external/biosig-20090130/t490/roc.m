function [SEN,SPEC,d,ACC,AREA,YI,c]=roc(d,c,color);
% ROC receiver operator curve and derived statistics.
% [...] = roc(d,c);
% d     DATA
% c     CLASS, vector with 0 and 1 
% 
% [...]=roc(d1,d0);
% d1    DATA of class 1 
% d2    DATA of class 0
% 
% [SEN, SPEC, TH, ACC, AUC,Yi,idx]=roc(...);
% OUTPUT:
% SEN     sensitivity
% SPEC    specificity
% TH      Threshold
% ACC     accuracy
% AUC     area under ROC curve
% Yi 	  max(SEN+SPEC-1), Youden index 
% c	  TH(c) is the threshold that maximizes Yi 

%	$Id: roc.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%	Copyright (c) 1997-2003,2005,2007 by  Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.
%


% Problem : Wenn die Schwellwerte mehrfach vorkommen, kann es zu Ambiguiten kommen, welche sich auf die AUC auswirken.

MODE = all(size(d)==size(c)) & all(all((c==1) | (c==0)))
d=d(:);
c=c(:);
        
if ~MODE
        d2=c;
        c=[ones(size(d));zeros(size(d2))];
        d=[d;d2];
        fprintf(2,'Warning ROC: XXX\n')        
end;        

% handle (ignore) NaN's  
c = c(~isnan(d));
d = d(~isnan(d));

if nargin<3
        color='-';
end;

[D,I] = sort(d);
x = c(I);

% identify unique elements
if 0,
        fprintf(2,'Warning ROC: please take care\n');
        tmp= find(diff(D)>0);
        tmp=sort([1; tmp; tmp+1; length(d)]);%]',2*length(tmp+2),1);
        %D=d(tmp);
end;

FNR = cumsum(x==1)/sum(x==1);
TPR = 1-FNR;

TNR = cumsum(x==0)/sum(x==0);
FPR = 1-TNR;

FN = cumsum(x==1);
TP = sum(x==1)-FN;

TN = cumsum(x==0);
FP = sum(x==0)-TN;

SEN = TP./(TP+FN);
SPEC= TN./(TN+FP);
ACC = (TP+TN)./(TP+TN+FP+FN);

% SEN = [FN TP TN FP SEN SPEC ACC D];

%fill(TN/sum(x==0),FN/sum(x==1),'b');
%SEN=SEN(tmp,:);
%ACC=ACC(tmp);
%d=D(tmp);
d=D;

%plot(SEN(:,1)*100,SPEC*100,color);
plot(FPR*100,TPR*100,color);
%plot(FP*100,TP*100,color);
% fill([1; FP],[0;TP],'c');

%ylabel('Sensitivity (true positive ratio) [%]');
%xlabel('1-Specificity (false positive ratio) [%]');

% area under the ROC curve
AREA = -diff(FPR)' * (TPR(1:end-1)+TPR(2:end))/2;

% Youden index
[YI,c] = max(SEN+SPEC-1);
