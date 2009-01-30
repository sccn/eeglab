function [p] = glob(Input,n1,samp,varargin)
% The function glob performs different global tests. A global test provides one
% joint statement on all endpoints, i.e. the global hypotheses is testing. 
% This function returns the corresponding p-value.
% 
% 
% 
%----------
%INPUT
%
%These input arguments are required:
% Input: data matrix with the size [n,k]               
% n1:	number of patients in group one (0 < n1 <= n ), 
%	restricted by the kind of samp
% samp: kind of sample         
%                             
%           single sample       'single' (n1 = n)
%           paired sample       'paired' (n1 = n/2; n must be even)
%           independent sample  'indept' (n1 < n)
%-----
%
%   [...] = glob(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%   
%    Parameter        Value  
%
%      'tail'       The alternative hypothesis against which to compute
%                   p-values for testing the hypothesis of no differences.
%                   Choices are:
%
%		       tail		 Alternative Hypothesis			
%		'~=' (the default) "there is a significant difference" (two-sided test)
%               '>'                "the values of group 1 are higher than the values of group 2" (one-sided test)
%               '<'                "the values of group 1 are smaller than the values of group 2" (one-sided test)    
%
%
%       'B'            number of permutations
%		                default: 500 
%                       B must be in the intervall
%                       500 <= B <= 2^n1   for single and paired sample
%                       (for 2^n1 < 500 : B = min(B,2^n1)) 
%
%                       500 <= B <= n! / n1!*(n-1)  for independent sample
%                       (for n! / n1!*(n-n1)! < 500 : B = min(B,n!/n1!*(n-n1)!))
%
%---
%
%    'tstat'              teststatistic
%                         'tmax' => is sensitive to departures in only a few endpoints 
%                         'tsum' => is sensitive to departures of all endpoints in the same direction
%                         'tsumabs' (the default) => is sensitive to departures of all endpoints (independent of the
%		                                direction); only two-sided!!!
%                         'ta' => choose this test statistic if the relative number of false hypotheses is
%		                                small; only independent!
%-----------
%
% OUTPUT
%
% [p] = glob(Input,n1,samp) returns the p-value of the global test.
%
%-----------
%
% REFERENCES:
%
%  [1] Hemmelmann C, Horn M, Reiterer S, Schack B, Suesse T, Weiss S.
%	Multivariate tests for the evaluation of high-dimensional EEG data.
%	J Neurosci Methods. 2004 Oct 15;139(1):111-20. 
%
%
% Copyright (C) 2006,2007 Claudia Hemmelmann <claudia.hemmelmann@mti.uni-jena.de>
% Adapted by A Schloegl <a.schloegl@ieee.org> 2006,2007 
%
%***
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
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
%--------------------------------------------------------------------------

%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

if nargin < 3 %check the required input of completeness
	error('stats:fdrin:TooFewInputs', ...
	      'global requires at least three input arguments: X; n1; samp .');
end %if

[n,k]=size(Input); %Dimension of the data matrix Input

%Initialisation t
%t = zeros(1,k);

% check the required input argument n1
if n1 <= 0
    error('n1 > 0');
end%if

n1check = floor(n1);

if n1check - n1 == 0
else
    error('Use for n1 the data type interger in the interval [1,n]');
end%if

% required input argument samp; use 'single','paired','indept'
switch samp %choose the kind of sample; required input arument
    case 'single' %single sample
        sampout = 'single sample';
        
        if n1 == n
         else
         error('For single sample must be n1 = n');
        end%if
    case 'paired' %paired sample
        sampout = 'paired sample';
        
        ncheck1 = n/2;
        ncheck2 = floor(ncheck1);
        if ncheck2 - ncheck1 == 0
        else
            error('n, the number of rows of the data matrix I, must be even for paired sample');
        end%if
        
        if n1 == n/2
            else
        error('For paired sample must be n1 = n/2;');
        end%if
    case 'indept' %independent sample
        sampout = 'independent sample';
        
        if n1 < n
        else
            error('For independent sample must be n1 < n');
        end%if
    otherwise
        error('chose for the kind of sample the Values: ''single'' for single sample; ''paired'' for paired sample; ''indept'' for independent sample');
end %switch samp


%Set default values of optional arguments
%----------------------------------------
tail = '~=';

if strcmp(samp,'indept') 
        B = min(500,prod(1:n)/(prod(1:n1)*prod(1:n-n1)));
    else
        B = min(500,2^n1);
end%if

tstat = 'tsumabs';    
    
%Set values of submitted optinoal arguments
%------------------------------------------
for i=1:2:length(varargin)
    switch (varargin{i})
        case 'tail'
            tail = varargin{i+1};
        case 'B'
            B = varargin{i+1};
        case 'tstat'
            tstat = varargin{i+1};
        otherwise
            error(sprintf('Unknow option: %s',varargin{i}));
    end%switch(varargin{i})
end%for


switch tail  %choose tail
   case '>' %one-sided test with X > Y
    tailout = 'one-sided (>)';
   case '~=' %two-sided test with X unequal Y, the default value
    tailout = 'two-sided';
   case '<' %one-sided test with X < Y
    tailout = 'one-sided (<)';
   otherwise
        error('ts:global:Unknowntail', [...
           'The ''tail'' parameter value must be ''>'' for one-sided test with X > Y,', ...
           ' ''~='' for two-sided test with X unequal Y or', ...
           ' ''<'' for one-sided test with X < Y .']); 
end %switch tail


%check B
switch samp
    case 'single'
        
        Bcheck = floor(B);
        if Bcheck - B == 0
        else
        error('Use for B the data type interger in the interval [500,2^n1]');
        end%if Bcheck
       
        if  B > 2^n1
            B = 2^n1;
            fprintf('You calculate with B = 2^n1');
            B
        elseif B < 500
     B = min(500,2^n1);
     fprintf('You calculate with B = min(500,2^n1) ');
     B
        else
            B = B;
        end%if
               
    case 'paired'
        
          Bcheck = floor(B);
        if Bcheck - B == 0
        else
        error('Use for B the data type interger in the interval [500,2^n1]');
        end%if Bcheck
        
        if  B > 2^n1
            B = 2^n1;
            fprintf('You calculate with B = 2^n1');
            B
        elseif B < 500
            B = min(500,2^n1);
            fprintf('You calculate with B = min(500,2^n1)');
            B
        else
            B = B;
        end%if
      
    case 'indept'
        
        Bcheck = floor(B);
        if Bcheck - B == 0
        else
        error('Use for B the data type interger in the interval [500,n!/(n1!*(n-n1)!)]');
        end%if Bcheck
         
        if prod(1:n)/(prod(1:n1)*prod(1:n-n1)) < B
            B = prod(1:n)/(prod(1:n1)*prod(1:n-n1));
            fprintf('You calculate with B = n!/(n1!*(n-n1)!)');
            B
        elseif B < 500
            B = min(500,prod(1:n)/(prod(1:n1)*prod(1:n-n1)));
            fprintf('You calculate with B = min(500, n!/(n1!*(n-n1)!))');
            B
        else
            B = B;
        end%if
        
end%switch samp check B

switch tstat %choose teststatistic
case 'tmax'
    tstat = 'tmax';
    tstatout = 'tmax (is sensitive to departures in only a few endpoints)';
case 'tsum'
    tstat = 'tsum';
    tstatout = 'tsum (is sensitive to departures of all endpoints in the same direction)';

case 'tsumabs'
    if strcmp(tail, '~=')
    tstat = 'tsumabs';
    tstatout = 'tsumabs (is sensitive to departures of all endpoints (independent of the direction); only two-sided)';
    else
    error('stats:global:Unknowntstat', ...
          'The ''tstat'' parameter value can only be ''tsumabs'' for two-sided test ''~='' .');
    end%if

case 'ta'
    if strcmp(tail,'~=') & strcmp(samp,'indept')
    tstat = 'ta';
    tstatout = 'ta (test statistic if the relative number of false hypotheses is small; only independent!)';
    else
    error('stats:global:Unknowntstat', ...
        'The ''tstat'' parameter value can only be ''ta'' for two-sided test ''~='' and independent sample ''indept'':');
    end%if
otherwise
    error('stats:global:Unknowntstat', ...
        'The ''tstat'' parameter value for the test statistic must be ''tmax'' (is sensitive to departures in only a few endpoints), ', ...
        '''tsum'' (is sensitive to departures of all endpoints in the same direction),', ...
        '''tsumabs'' (is sensitive to departures of all endpoints (independent of the direction); only two-sided) or)', ...
        '''ta'' (choose this test statistic if the relative number of false hypotheses is small; only independent).');
end
%Calculation
%--------------------------------------------------------------------------

[p] = perm_gfwe(Input,n1,samp,B,tstat,tail);


%Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

disp(sprintf('\nThe following parameters have been used: '));
disp(sprintf('sample size: %g, (group 1: %g)',n,n1));
disp(sprintf('number of hypotheses: %g',k));
disp(sprintf('kind of sample: %s',sampout));
disp(sprintf('tail: %s',tailout));
disp(sprintf('number B of permutations: %g',B));
disp(sprintf('test statistic: %s',tstatout));

disp(sprintf('\nResult (p-value):'));



%end function[p] = glob(Input,n1,samp,varargin)
