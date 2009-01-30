function [O] = fdr(Input,n1,samp,varargin)
%
% The function FDR performs different multiple test procedures for 
% controlling the false discovery rate (FDR). 
% 
% function [O] = fdr(Input,n1,samp) returns the number of rejected hypotheses,
% the rank (O(:,1)), the indices of the rejected hypotheses (O(:,2)), the adjusted p-values
% (O(:,3)) and the unadjusted p-values (O(:,4)) for the procedure of Benjamini
% and Yekutieli (1995) with the significance level alpha=0.05.
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
%   [...] = fdr(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%   
%    Parameter        Value
%
%     'test'         Value for single sample
%    		          'ttest'               to compute the t-Test
%                                           assumption : normal(gaussian) distribution   
%                     'wilcox'              to compute the Wilcoxen signed rank test
%                                           assumption : symmetrical distribution
%                     'sign' (the default)  to compute the sign-test 
%                                           assumption : none  
%            
%                    Value for paired sample
%                     'ttest'                to compute the t-Test
%                                            assumption : normal(gaussian) distribution
%                     'wilcox'               to compute the Wilcoxen signed rank test
%                                            assumption : symmetrical distribution
%                     'sign' (the default)   to compute the sign-test
%                                            assumption : none  
%              
%                    Value for independent sample
%                     'ttest'                to compute the t-Test
%                                            assumption : normal(gaussian) distribution 
%                     'wilcox' (the default) to compute the Wilcoxen rank test (Wilcoxen-Man-Whitney-Test)
%                                            assumption : none
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
%---
%      'proc'        'BH' (the default)     chooses the procedure of Benjamini and Hochberg (1995)
%                    'BL'                   chooses the procedure of Benjamini and Liu (2001)
%                    'BKY'                  chooses the procedure of Benjamini, Krieger and Yekutieli (2001)
%---
%      'alpha'       0.05 (the default)    significance level
%                     for a other value: 0<alpha<=0.2
%-----------
%
% OUTPUT
%
% [O] = fdr(Input,n1,samp) returns the rank (O(:,1)), 
% the indices of the rejected hypotheses (O(:,2)), 
% the adjusted p-values (O(:,3)) and the unadjusted p-values (O(:,4)).
%
%-----------
%
% REFERENCES:
%
%  [1] Hemmelmann C, Horn M, Suesse T, Vollandt R, Weiss S.
%	New concepts of multiple tests and their use for evaluating high-dimensional EEG data.
%	J Neurosci Methods. 2005 Mar 30;142(2):209-17.
%  [2] Hemmelmann C, Horn M, Reiterer S, Schack B, Suesse T, Weiss S.
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
	      'fdr requires at least three input arguments: X; n1; samp .');
end %if

[n,k]=size(Input); %Dimension of the data matrix Input

%Initialisation t
t = zeros(1,k);

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
if  strcmp(samp,'indept')
    test = 'wilcox';
else
    test = 'sign';
end%if
tail = '~=';
proc = 'BH';
alpha = 0.05;

%Set values of submitted optinoal arguments
%------------------------------------------
for i=1:2:length(varargin)
    switch (varargin{i})
        case 'test'
            test = varargin{i+1};
        case 'tail'
            tail = varargin{i+1};
        case 'proc'
            proc = varargin{i+1};
        case 'alpha'
            alpha = varargin{i+1};
        otherwise
            error(sprintf('Unknow option: %s',varargin{i}));
    end%switch(varargin{i})
end%for


% Validate the test parameter
testChoices = {'ttest' 'wilcox' 'sign'};

if strcmp(samp ,'single')       %for single sample
	switch test           %choose test   
	  case 'ttest'        %t-Test
          testout = 't-test';
  	  case 'wilcox'       %Wilcoxen signed rank test
          testout = 'Wilcoxen signed rank test';
      case 'sign'         %sign test
          testout = 'sign test';
      otherwise
          error('stats:fdrin:Unknowntest', ...
          'The ''test'' parameter value must be ''ttest'' , ''wilcox'' or ''stest'' for single sample .');
   	end %switch test sing
elseif strcmp(samp,'paired')  % for paired sample
   	switch test          %choose test  
    	 case 'ttest'    %t-Test
             testout = 't-test';
    	 case 'wilcox'   %Wilcoxen signed rank test
             testout = 'Wilcoxen signed rank test';
    	 case 'sign'     %sign test
             testout = 'sign test';
         otherwise
          error('stats:fdrin:Unknowntest', ...
          'The ''test'' parameter value must be ''ttest'', '' wilcox'' or ''sign'' for paired sample .');
	end %switch test pair
else    % for independent sample
        switch test %choose test
	 case 'ttest'        %t-Test 
         testout = 't-test';
	 case 'wilcox'       %Wilcoxen rank sum test
         testout = 'Wilcoxen-Man-Whitney-Test';
     otherwise
          error('stats:fdrin:Unknowntest', ...
          'The ''test'' parameter value must be ''ttest'' or ''wilcox'' for independent sample .');
     end %switch test indept
end %if samp

switch tail  %choose tail
   case '>' %one-sided test with X > Y
    tailout = 'one-sided (>)';
   case '~=' %two-sided test with X unequal Y, the default value
    tailout = 'two-sided';
   case '<' %one-sided test with X < Y
    tailout = 'one-sided (<)';
   otherwise
        error('ts:fdrin:Unknowntail', [...
           'The ''tail'' parameter value must be ''>'' for one-sided test with X > Y,', ...
           ' ''~='' for two-sided test with X unequal Y or', ...
           ' ''<'' for one-sided test with X < Y .']); 
end %switch tail

switch proc       % choose procedure
   case 'BH'         % Benjamini Hochberg procedure
       procout = 'procedure of Benjamini and Hochberg (1995)';
   case 'BL'         % Benjamini Lui procedure
       procout = 'procedure of Benjamini and Liu (2001)';
   case 'BKY'        % Benjamini Krieger Yekutieli procedure
       procout = 'procedure of Benjamini, Krieger and Yekutieli (2001)';
   otherwise
        error('ts:fdrin:Unknownproc', [...
            'The ''proc'' parameter value must be ''BH'' for Benjamini and Hochberg,', ...
            ' ''BL'' for Benjamini and Liu or ''BKY'' for Benjamini, Krieger and Yekutieli.']);
end %switch proc

%check alpha 0<alpha<=0.5
if alpha<= 0
    error('alpha must be in the interval (0,0.2]');
elseif alpha>0.2
    error('apha must be in the interval (0,0.2]');
else
end%if

%Calculation
%--------------------------------------------------------------------------

%initialise matrices
if strcmp(samp,'paired')
    
        X = Input(1:n1,1:k);
        Y = Input(n1+1:n,1:k);
    
elseif strcmp(samp,'single')
    
        X = Input;
        Y = zeros(n,k);
    
elseif strcmp(samp,'indept')
    
        X = Input(1:n1,1:k);
        Y = Input(n1+1:n,1:k);
    
end%if

%choose test
if strcmp(samp,'single')       %for single sample
   switch test           %choose test   
	case 'ttest'        %t-Test
	  t=ttestC(X-Y); 
          FG=n1-1;
          p=tcdf(t,FG);
  	case 'wilcox'     %Wilcoxen signed rank test
   	  p=zeros(1,k);
          for j=1:k,
            [p(j),h]=wilcoxon_test(X(:,j),Y(:,j),alpha,tail); % Wilcoxentest für gepaarte Stichproben
	  end;%for
        case 'sign'       %sign test
          p=zeros(1,k);
          for j=1:k,
            [p(j),h]=signtest(X(:,j),Y(:,j),alpha,tail); % Vorzeichentest
          end;%for
   	end %switch test sing
elseif strcmp(samp,'paired')  %for paired sample
   switch test            
    	 case 'ttest'    %t-Test
	   t=ttestC(X-Y); 
           FG=n1-1;
           p=tcdf(t,FG); %einseitiger p-Wert des t-Tests  
    	 case 'wilcox'   %Wilcoxen signed rank test
    	   p=zeros(1,k);
           for j=1:k,
             [p(j),h]=wilcoxon_test(X(:,j),Y(:,j),alpha,tail); % Wilcoxentest für gepaarte Stichproben
           end;%for
    	 case 'sign'     %sign test
      	   p=zeros(1,k);
           for j=1:k,
             [p(j),h]=signtest(X(:,j),Y(:,j),alpha,tail); %Vorzeichentest
           end;%for
   end %switch test pair
else    % for independent sample
   switch test %choose test
	 case 'ttest'        %t-Test 
	   t=ttest3(Input,n1);
           FG=n-2;
           p=tcdf(t,FG);%einseitiger p-Wert des t-Tests 
	 case 'wilcox'       %Wilcoxen rank sum test
           p=zeros(1,k);
	   for j=1:k,
             [p(j),h] = u_test(X(:,j),Y(:,j),tail);
	   end;
   end %switch test indept
end %if samp

%choose tail
if strcmp(test, 'ttest')
   switch tail
	case '>' %one-sided test, with X > Y
	   p=1-p;
	case '~=' %two-sided test, with X unequal Y
	   p=2*min(p,1-p);;
	case '<' %one-sided test, with X < Y
    	   p=p;
   end %switch tail
end%if
%resort the p-values
[psd,indexsd]=sort(p);	        % sorted, start with the smallest
   h1(1,:)=psd;
   h1(2,:)=indexsd;
   h2=flipdim(h1,2);		% sorted, start with the greatest
   psu=h2(1,:);
   indexsu=h2(2,:);

% choose procedure
switch proc      
   case 'BH'         % Benjamini Hochberg procedure 
 	[AaH,pup,adpval] = bh95(psu,indexsu,k,alpha);
   case 'BL'         % Benjamini Liu procedure
 	[AaH,pup,adpval] = bl01(psd,indexsd,k,alpha);
   case 'BKY'   
	[AaH,pup,adpval]= zweistufen(psu,indexsu,k,alpha) ; % Benjamini Krieger Yekutieli procedure
end %switch proc


%Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

disp(sprintf('The following parameters have been used: '));
disp(sprintf('sample size: %g, (group 1: %g)',n,n1));
disp(sprintf('number of hypotheses: %g',k));
disp(sprintf('kind of sample: %s',sampout));
disp(sprintf('univariate test: %s',testout));
disp(sprintf('tail: %s',tailout));
disp(sprintf('multiple test procedure: %s',procout));
disp(sprintf('significance level alpha: %1.3f',alpha));

disp(sprintf('\nResult:'));
disp(sprintf('\nnumber of rejected hypotheses: %g',AaH));
disp(sprintf('\nOutputmatrix'));

rangzahl=[1:1:AaH]';
format short g
O = [rangzahl, pup', adpval', psd(1:AaH)'];

disp(sprintf('\nrank \t component \t p-values (adjusted) \t p-value (unadjusted)'));    
for ausg_O=1:AaH,
    disp(sprintf('%g\t\t',squeeze(O(ausg_O,:)))); 
     
end;


%end function[O] = fdr(Input,n1,samp,varargin)
