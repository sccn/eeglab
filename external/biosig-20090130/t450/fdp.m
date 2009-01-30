function [O] = fdp(Input,n1,samp,varargin)
%
% The function fdp performs different multiple test procedures for 
% controlling the false discovery proportion (FDP). 
% 
% function [O] = fdp(Input,n1,samp) returns the number of rejected hypotheses,
% the rank (O(:,1)), the indices of the rejected hypotheses (O(:,2)) and the unadjusted 
% p-values (O(:,3)) for the procedure B (conservative) of Korn et al. (2004)
% with the significance level alpha=0.05.
%
%
%----------
%INPUT
%
% These input arguments are required:
% Input: data matrix with the size [n,k]               
% n1:	number of patients in group one (0 < n1 <= n ), 
%	restricted by the kind of samp
% samp: kind of sample         
%                             
%           single sample       'single' (n1 = n)
%           paired sample       'paired' (n1 = n/2; n must be even)
%           independent sample  'indept' (n1 < n)
%
% [...] = fdp(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
% parameters and their values. Valid parameters are the following:
%   
%    Parameter        Value
%
%    'gamma'        gamma = 0.1; (the default)
%                     for an other value: 0< gamma <= 0.5
%
%     'B'            number of permutations
%		     default: 500 
%                    B must be in the intervall
%                       500 <= B <= 2^n1   for single and paired sample
%                           		   (for 2^n1 < 500 : B = min(B,2^n1)) 
%
%                       500 <= B <= n! / n1!*(n-1)  for independent sample
%                       (for n! / n1!*(n-n1)! < 500 : B = min(B,n!/n1!*(n-n1)!))
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
%			'~=' (the default)  "there is a significant difference" (two-sided test)
%            '>'                "the values of group 1 are higher than the values of group 2" (one-sided test)
%            '<'                "the values of group 1 are smaller than the values of group 2" (one-sided test)    
%
%---
%
%      'proc'        'ProcBv' (the default)  chooses the procedure B (conservative) of Korn et al. (2004)
%                    'ProcBe'                chooses the procedure B of Korn et al. (2004)
%                    'TL'                    chooses the procedure of Troendle (1995) and the extention of van der Laan et al.
%                    'LR1'                   chooses the procedure of Lehmann and Romano (2005) with some dependence
%                                                    assumptions or asymtotic control (see Romano and Wolf
%                                                    (2005) "Control of Generalized Error Rates in Multiple Tetsing")
%                    'LR2'                   chooses the procedure of Lehmann and Romano (2005) without any dependence assumptions (conservative!)             
%                    'HL'                    chooses the procedure of Holm and the extention of van der Laan et al.
%---
%
%      'alpha'       0.05 (the default)    significance level
%                    alpha must be a scalar and in the interval 0 < alpha <= 0.2
%
% OUTPUT
%
% [O] = fdp(Input,n1,samp) returns the rank (O(:,1)), 
% the indices of the rejected hypotheses (O(:,2)) and 
% the adjusted p-values (O(:,3)).
%
%
%-----------
%
% REFERENCES
%  [1]	Hemmelmann, C., Horn, M., Süße, T., Vollandt, R., Weiss, S. (2005):
%       New concepts of multiple tests and their use for evaluating
%       high-dimensional EEG data, Vol 142/2 pp 209-217.
%
%
%
%Copyright (C) 2006 by Claudia Hemmelmann <claudia.hemmelmann@mti.uni-jena.de>
%Institute of Medical Statistics, Computer Sciences and Documantation
%University of Jena
%This work was supported by DFG Project VO 683/2-1
%This is part of the BIOSIG-toolbox http://biosig.sf.net/
%---
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

%INPUT check
%--------------------------------------------------------------------------
%
%required argumets
%"""""""""""""""""
if nargin < 3 %check the required input of completeness
	error('stats:fdrin:TooFewInputs', ...
	      'fdri requires at least three input arguments: X; n1; samp .');
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

% choose the required input argument samp; use 'single','paired','indept'
switch samp
    case 'single'
        sampout = 'single sample'; %single sample
        if n1 == n
         else
         error('For single sample must be n1 = n');
        end%if
    case 'paired'
        sampout = 'paired sample'; %paired sample
        ncheck1 = n/2;
        ncheck2 = floor(ncheck1);
        if ncheck2 - ncheck1 == 0
        else
            error('n, the number of rows of the data matrix Input, must be even for paired sample');
        end%if
        
        if n1 == n/2
            else
        error('For paired sample must be n1 = n/2;');
    end%if
    case     'indept'
        sampout = 'independent sample'; %independent sample
        if n1 < n
        else
            error('For independent sample must be n1 < n');
        end%if
    otherwise
        error('chose for the kind of sample the Values: ''single'' for single sample; ''paired'' for paired sample; ''indept'' for independent sample');
end %switch samp

%
%optional arguments
%""""""""""""""""""


%Set default values of optional arguments

%gamma
gamma = 0.1;
%test
if  strcmp(samp,'indept')
    test = 'wilcox';
else
    test = 'sign';
end%if
%tail
tail = '~=';
%proc
proc = 'Bv';
%B
if strcmp(samp,'indept')
    B = min(500,prod(1:n)/(prod(1:n1)*prod(1:n-n1)));
else
B = min(500,2^n1);
end%if
%alpha    
alpha = 0.05;


%Set values of submitted optinoal arguments
for i=1:2:length(varargin)
    switch (varargin{i})
        case 'gamma'
            gamma = varargin{i+1};
        case 'test'
            test = varargin{i+1};
        case 'tail'
            tail = varargin{i+1};
        case 'proc'
            proc = varargin{i+1};
        case 'B'
            B = varargin{i+1};
        case 'alpha'
            alpha = varargin{i+1};
        otherwise
            error(sprintf('Unknow option: %s',varargin{i}));
    end%switch(varargin{i})
end%for


%check gamma 0<gamma<=0.5
if gamma < 0
    error('gamma must be in the interval (0,0.5]');
elseif gamma > 0.5
    error('gamma must be in the interval (0,0.5]');
else
end%if

% Validate the test parameter
testChoices = {'ttest' 'wilcox' 'sign'};

if strcmp(samp,'single')       %for single sample
	switch test           %choose test   
	  case 'ttest'        %t-Test
          testout = 't-test';
  	  case 'wilcox'       %Wilcoxen signed rank test
          testout = 'Wilcoxen signed rank test';
      case 'sign'         %sign test
          testout = 'sign test';
      otherwise
          error('stats:fdrin:Unknowntest', ...
          'The ''test'' parameter value must be ''ttest'' , ''wilcox'' or ''sign'' for single sample .');
   	end %switch test sing
elseif strcmp(samp,'paired')  % for paired sample
   	switch test          %chose test  
    	 case 'ttest'    %t-Test
             testout = 't-test';
    	 case 'wilcox'   %Wilcoxen signed rank test
             testout = 'Wicoxen signed rank test';
    	 case 'sign'     %sign test
             testout = 'sign test';
         otherwise
          error('stats:fdrin:Unknowntest', ...
          'The ''test'' parameter value must be ''ttest'', '' wilcox'' or ''sign'' for paired sample .');
	end %switch test pair
else    % for independent sample
        switch test %chose test
	 case 'ttest'        %t-Test 
         testout = 't-test';
	 case 'wilcox'       %Wilcoxen rank sum test
         testout = 'Wolcoxen-Man-Whitney-Test';
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
        error('ts:fdrin:Unknowntail',[...
           'The ''tail'' parameter value must be ''>'' for one-sided test with X > Y,', ...
           ' ''~='' for two-sided test with X unequal Y or', ...
           ' ''<'' for one-sided test with X < Y .']); 
end %switch tail

% Validate the procedure parameter

switch proc               
    case 'Bv'
        procout = 'procedure B (conservative) of Korn et al. (2004)';
    case 'Be'             
        procout = 'procedure B of Korn et al. (2004)';
    case 'TL'             
        procout = 'procedure of Troendle (1995) and the extention of van der Laan et al.';
    case 'LR1'            
        procout = 'procedure of Lehmann and Romano (2005) with some dependence assumptions';
    case 'LR2'            
        procout = 'procedure of Lehmann and Romano (2005) without any dependence assumptions';
    case 'HL'             
        procout = 'procedure of Holm and the extention of van der Laan et al.';
    otherwise
        error('ts:fdp:Unknownproc', [...
            'The ''proc'' parameter value must be ''Bv'' for procedure B (conservative),', ...
            ' ''Be'' for procedure B, ''TL'' for the procedure of Troendle + van der Laan,',...
            ' ''LR1'' for the procedure of Lehmann + Romano with some dependence assumptions, ''LR2'' for the procedure of Lehmann + Romano without any dependence assumptions',...
            ' or ''HL'' for the procedure of Holm + van der Laan.']);
end%switch 

if strcmp(proc,'Bv') | strcmp(proc,'Be') | strcmp(proc,'TL')  %only for procedures with permutation tests
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
end%if

%check alpha 0<alpha<=0.5
if alpha<= 0
    error('alpha must be in the interval (0,0.2]');
elseif alpha>0.2
    error('apha must be in the interval (0,0.2]');
else
end%if

%Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

if strcmp(proc,'LR1') | strcmp(proc,'LR2') | strcmp(proc,'HL') % Nicht-Permutationstests in pwerte

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
	switch test              
	  case 'ttest'        %t-Test
	  t=ttestC(X-Y); 
        FG=n1-1;
        p=tcdf(t,FG);
  	  case 'wilcox'     %Wilcoxen signed rank test
   	    p=zeros(1,k);
            for j=1:k,
                [p(j),h] = wilcoxon_test(X(:,j),Y(:,j),alpha,tail); % Wilcoxentest für gepaarte Stichproben
            end;%for
      case 'sign'       %sign test
          p=zeros(1,k);
          for j=1:k,
            [p(j),h] = signtest(X(:,j),Y(:,j),alpha,tail); %Vorzeichentest
          end;%for
   	end %switch test sing
elseif strcmp(samp,'paired')  %for paired sample
   	switch test          %chose test  
    	 case 'ttest'    %t-Test
	    t=ttestC(X-Y); 
        FG=n1-1;
    %einseitiger p-Wert des t-Tests  
        p=tcdf(t,FG);
    	 case 'wilcox'   %Wilcoxen signed rank test
    	  p=zeros(1,k);
            for j=1:k,
                [p(j),h] = wilcoxon_test(X(:,j),Y(:,j),alpha,tail); % Wilcoxentest für gepaarte Stichproben
            end;%for
    	 case 'sign'     %sign test
      	  p=zeros(1,k);
          for j=1:k,
            [p(j),h] = signtest(X(:,j),Y(:,j),alpha,tail); %Vorzeichentest
          end;%for
	end %switch test pair
else    % for independent sample
        switch test %chose test
	 case 'ttest'        %t-Test 
	    
        t=ttest3(Input,n1);
        FG=n-2;
        %einseitiger p-Wert des t-Tests 
        p=tcdf(t,FG);
	 case 'wilcox'       %Wilcoxen rank sum test
        p=zeros(1,k);	
        for j=1:k,
         [p(j),h] = u_test(X(:,j),Y(:,j),tail);
        end;
     end %switch test indept
end %if samp

if strcmp(test,'ttest')
    %choose tail
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

end%if proc
   
% choose procedure
switch proc      
case 'Bv'         % 
[AaH,pup,adpval,q]=vereinM_B(Input,n1,samp,gamma,B,test,tail,alpha,k);
case 'Be'         %
[AaH,pup,adpval,q]=exakteM_B(Input,n1,samp,gamma,B,test,tail,alpha,k);
case 'TL' %Troendle + van der Laan
[AaH,pup,adpval,q]=vereinM_A(Input,n1,samp,0,B,test,tail,alpha,k);

if AaH < k & AaH ~= 0,
   zahl=[0:1:k-AaH];
   suche=zahl./(zahl+ones(1,length(zahl))*AaH);
   vergleich=(suche<=gamma);
   maxi=sum(vergleich)-1;
   
   if AaH + maxi < k,
      AaH = AaH + maxi;
      %gFWE beachten
      [psd,indexsd]=sort(q);	        % sorted, start with the lowest
      pup=indexsd(1:AaH);
   else
      AaH =  k;
      [psd,indexsd]=sort(q);	        % sorted, start with the lowest
      pup=indexsd(1:AaH);
   end;%if
end;%if

case 'LR1' %
    p = psd;
[AaH,pwert]=lehrom(p,k,alpha,gamma,1); % mit Vorauss. oder nur asymptotisch
case 'LR2' %
p = psd;
[AaH,pwert]=lehrom(p,k,alpha,gamma,0); % ohne Voraussetzung
case 'HL' %
 p = psd;
[AaH,pwert]=homhof(p,k,alpha,0); %= Holm
if AaH < k & AaH ~= 0,
   zahl=[0:1:k-AaH];
   suche=zahl./(zahl+ones(1,length(zahl))*AaH);
   vergleich=(suche<=gamma);
   maxi=sum(vergleich)-1;
   
   if AaH + maxi < k,
      AaH = AaH + maxi;
   else
      AaH = k;
   end;%if
end;%if
end %switch proc






%Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%Output matrix O

disp(sprintf('\nThe following parameters have been used: '))
disp(sprintf('sample size: %g, (group 1: %g)',n,n1))
disp(sprintf('number of hypotheses: %g',k))
disp(sprintf('kind of sample: %s',sampout))
disp(sprintf('gamma: %1.3f',gamma))
if strcmp(proc,'Bv') | strcmp(proc,'Be') | strcmp(proc,'TL')
disp(sprintf('number B of permutations: %g',B));
end%if
disp(sprintf('univariate test: %s',testout));
disp(sprintf('tail: %s',tailout));
disp(sprintf('multiple test procedure: %s',procout));
disp(sprintf('significance level alpha: %1.3f\n',alpha));

disp(sprintf('Result:\n'));
disp(sprintf('number of rejected hypotheses: %g\n',AaH));

disp(sprintf('Outputmatrix'));

rangzahl=[1:1:AaH]';
if strcmp(proc,'Bv') | strcmp(proc,'Be') | strcmp(proc,'TL')
   if pup > 0 
    pw = q(pup)';
    pup = pup';
    format short g
    O = [rangzahl,pup,pw];
   else
       O = [];
   end
else
    if AaH > 0
    indexsd = indexsd(1:AaH)';	% index of the sorted p-values psd
    psd = psd(1:AaH)';	        % sorted p-values, start with lowest  
    format short g
    O = [rangzahl,indexsd,psd];
    else
        O = [];
    end
end;%if
    
disp(sprintf('\nrank \t  component \t p-value (unadjusted)'));    
for ausg_O=1:AaH,
    disp(sprintf('%g\t\t\t',squeeze(O(ausg_O,:)))); 
     
end;%for




%end function[O] = fdp(Input,n1,samp,varargin)
