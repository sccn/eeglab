% bootstat() - accumulate surrogate data to assess significance by bootstrap of 
%              some measure of two input variables. Fits the psd with a 4th-order 
%              polynomial using the distribution kurtosis. Reference: Ramberg, J.S., 
%              Tadikamalla, P.R., Dudewicz E.J., Mykkytka, E.F. "A probability 
%              distribution and its uses in fitting data." Technimetrics, 1979, 
%              21:201-214.
% Usage:
%            >> [rsignif,accarray] = bootstat(arg1, arg2, formula, varargin ...);
% Inputs:
%    arg1    - [array] 1-D, 2-D or 3-D array of values
%    arg2    - [array] 1-D, 2-D or 3-D array of values
%    formula - [string] formula to compute the given measure. Takes arguments
%                   'arg1', 'arg2' as inputs and 'res' (result, by default) as output.
%                   For data arrays of more than 1 dimension, the formula must be iterative 
%                   so that shuffling can occur at each step while scanning the last 
%                   array dimension.  Examples:
%                   'res = arg1 - arg2'              % difference of two 1-D data arrays
%                   'res = mean( arg1 .* arg2)'      % mean projection of two 1-D data arrays
%                   'res = res + arg1 .* conj(arg2)' % iterative, for use with 2|3-D arrays
% Optional inputs:
%   'boottype '   - ['first'|'second'|'both'|'both2']   {default: 'both'}
%                   'first' = accumulate surrogate data by shuffling the first dimension 
%                             only. Ex: for (trials,times) data, shuffle trials.
%                   'second'= shuffle the second dimension (Ex: times), keeping 
%                             the first dimension (Ex: trials) constant. 
%                   'both'  = shuffle both dimensions repetetively during accumulation.
%                   'both2' = same as 'both', but shuffle the second dimension only once 
%                             per cycle through the data.  
%   'alpha'       - [real] significance level (between 0 and 1) {default 0.05}.
%   'naccu'       - [integer] number of exemplars to accumulate {default 200}.
%   'bootside'    - ['both'|'upper'] side of the surrogate distribution to
%                   consider for significance. This parameter affects the size
%                   of the last dimension of the accumulation array ('accres') 
%                   (size is 2 for 'both' and 1 for 'upper') {default: 'both'}.
%   'basevect'    - [integer vector] time vector indices for baseline
%                   {default: all time points}.
%   'accarray'    - accumulation array (from a previous call). Allows faster 
%                   computation of the 'rsignif' output {default: none}.
%   'formulainit' - [string] for initializing the output variable. Ex: 'res=zeros(10,40);'
%                   {default: 'res = (size(arg1,3) x naccu)'}
%   'formulapost' - [string] transformation to apply after accumulation. 
%                   Ex: 'res = res /10;' {default: none}.
%   'formulaout'  - [string] name of the computed variable {default: 'res'}.
%   'distfit'     - ['on'|'off'] fit distribution with known function to compute more accurate 
%                   limits or exact p-value (see 'vals' option). the MATLAB statistical toolbox 
%                   is required. This option is currently implemented only for 1-D data.
%   'vals'        - [float array] significance values. 'alpha' is ignored and 
%                   rsignif returns the p-values. Requires 'distfit' (see above).
%                   This option currently implemented only for 1-D data.
%   'correctp'    - [phat pci zerofreq] parameters for correcting for a biased probability 
%                   distribution (requires 'distfit' above). See help of correctfit().
% Outputs: 
%    rsignif      - significance arrays. 2 values (low high) for each point (use
%                   'alpha' to change these limits).
%    accarray     - accumulated surrogate data values.
%
% Authors: Arnaud Delorme, Lars Kai Hansen & Scott Makeig
%          SCCN/INC, UCSD, La Jolla, 2002-
%
% See also: timef()

% NOTE: There is one hidden parameter 'savecoher', 0 or 1
% HELP TEXT REMOVED:          (Ex: Using option 'both', 
%                             coherence during baseline would be ignored since times
%                             are shuffled during each accumulation.

% Copyright (C) 9/2002  Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD
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

% $Log: not supported by cvs2svn $
% Revision 1.20  2003/12/09 22:29:07  arno
% typo in text
%
% Revision 1.19  2003/11/04 18:39:11  arno
% header
%
% Revision 1.18  2003/08/06 14:42:58  arno
% *** empty log message ***
%
% Revision 1.17  2003/07/09 22:04:00  arno
% rouding i to remove warning
%
% Revision 1.16  2003/07/09 21:37:39  arno
% typo
%
% Revision 1.15  2003/07/09 00:34:19  arno
% implementing correctp
%
% Revision 1.14  2003/07/08 16:36:25  arno
% no plot
%
% Revision 1.13  2003/07/08 15:37:02  arno
% *** empty log message ***
%
% Revision 1.12  2003/05/23 23:16:47  arno
% aded warning
%
% Revision 1.11  2003/05/23 00:58:17  arno
% larger range for 'vals' parameters
%
% Revision 1.10  2003/05/02 18:07:11  arno
% debuging 1-D bootstrap
%
% Revision 1.9  2003/04/23 00:31:20  arno
%  fitting with normal distribution for 1-D data
%
% Revision 1.8  2003/04/18 18:55:28  arno
% comment to fit with normal distribution
%
% Revision 1.7  2003/04/17 21:52:40  arno
% can now process 1-D arrays
%
% Revision 1.6  2003/01/06 19:46:21  arno
% debugging new bootstrap
%
% Revision 1.5  2003/01/06 19:26:56  arno
% implementing new bootstrap type
%
% Revision 1.4  2002/12/31 04:37:19  scott
% header edits
%
% Revision 1.3  2002/12/07 02:54:27  arno
% adding message for how to implement new statistics
%
% Revision 1.2  2002/10/02 00:36:04  arno
% update condstat, debug
%
% Revision 1.1  2002/10/01 16:07:33  arno
% Initial revision
%
% Revision 1.1  2002/09/24 23:28:02  arno
% Initial revision
%
% function for bootstrap initialisation
% -------------------------------------

% *************************************
% To fit the psd with as 4th order polynomial using the distribution kurtosis,
% Reference: Ramberg, J.S., Tadikamalla, P.R., Dudewicz E.J., Mykkytka, E.F. 
% "A probability distribution and its uses in fitting data." 
% Technimetrics, 1979, 21: 201-214.
% *************************************

function [accarrayout, Rbootout] = bootstat(oriarg1, oriarg2, formula, varargin)
%	nb_points, timesout, naccu, baselength, baseboot, boottype, alpha, rboot);
	
if nargin < 3
	help bootstat;
	return;
end;

if ~isstr(formula)
	error('The first argument must be a string formula');
end;

g = finputcheck(varargin, ...
                { 'dims'          'integer'  []                       []; ...
                  'naccu'         'integer'  [0 10000]                200; ...
                  'bootside'      {'cell' 'string'}   { [] {'both' 'upper'} }         'both'; ...
                  'basevect'      'integer'  []                       []; ...
                  'formulapost'   'string'   []                       ''; ...
                  'formulainit'   'string'   []                       'res = zeros(nb_points,g.naccu);'; ...
                  'formulaout'    {'cell' 'string'} []                       'res'; ...
				  'boottype'      'string'  {'first' 'second' 'both' 'times' 'trials' 'timestrials' 'timestrials2' 'both2'} 'both'; ...
				  'alpha'         'real'     [0 1]                    0.05; ...
				  'vals'          'real'     [0 Inf]                    []; ...
				  'distfit'       'string'   {'on' 'off' }            'off'; ...
				  'correctp'      'real'     []                       []; ...
				  'accarray'      'integer'  []                       NaN	});
if isstr(g)
	error(g);
end;
bootname = g.boottype;
unitname = '';
switch g.boottype
 case 'times',  g.boottype = 'first'; unitname = 'trials'; bootname = 'times';
 case 'trials', g.boottype = 'second'; unitname = 'times'; bootname = 'trials';
 case 'timestrials',  g.boottype = 'both'; unitname = 'trials'; bootname = 'timestrials';
 case 'timestrials2', g.boottype = 'both2'; unitname = 'trials'; bootname = 'timestrials';
end;
if ~iscell(g.formulaout)
    g.formulaout = { g.formulaout };
end;
if ~iscell(g.bootside)
    g.bootside = { g.bootside };
end;    
if isempty(g.accarray)
	g.accarray = NaN;
end;

% function for bootstrap computation
% ----------------------------------
[nb_points times trials]    = size(oriarg1);

%trials    = size(oriarg1, 1);
%times     = size(oriarg1, 2);
%nb_points = size(oriarg1, 3);

if isempty(g.basevect)
	g.basevect = 1:times;
    if length(g.basevect) == 1, disp('Warning 1 value only for shuffling dimension'); end;
end;

% vector of only one dimension
% -----------------------------
if (size(oriarg1,1) == 1 | size(oriarg1,2) == 1) & size(oriarg1,3) == 1
    if isnan(g.accarray)
        fprintf('Bootstrap type is 1-D\n');
        fprintf('Processing %d bootstrap accumulation\n', g.naccu);
        oriarg1 = oriarg1(g.basevect);
        oriarg2 = oriarg2(g.basevect);
        arg1 = oriarg1;
        for index = 1:g.naccu
            arg2 = shuffle(oriarg2);
            eval([ formula ';' ]);
            eval([ 'Rbootout(index) = ' g.formulaout{1} ';' ]);
        end
    else 
        Rbootout = g.accarray;
    end;
    tmpsort = sort(abs(Rbootout));
    i = round(g.alpha*g.naccu);
    sigval = [mean(tmpsort(1:i)) mean(tmpsort(g.naccu-i+1:g.naccu))];
    if strcmpi(g.bootside, 'upper'), sigval = sigval(2); end;
    accarrayout = sigval;
else

    % array of 2 or 3 dimensions
    % --------------------------
    fprintf('\nBootstrap baseline length is %d (out of %d) points\n', length(g.basevect), times);
    if isnan(g.accarray)
        eval( g.formulainit );
        if strcmpi( g.boottype, 'second') % get g.naccu bootstrap estimates for each time window
            fprintf('Bootstrap type is 3-D, shuffling oalong dimension 3 only.\n');
            fprintf('Processing %s (naccu=%d) bootstrap (of %s %d):\n',bootname, g.naccu, unitname, times(end));
            
            arg2 = zeros(nb_points, g.naccu);
            arg1 = zeros(nb_points, g.naccu);
            for index= 1:length(g.formulaout) 
                eval( ['accarray' int2str(index) '= zeros(nb_points, g.naccu, times);'] );
            end;
            
            for index=1:times % dim2
                if rem(index,10) == 0,  fprintf(' %d',index); end
                if rem(index,120) == 0, fprintf('\n'); end
                for allt=1:trials % dim1
                    j=1;
                    while j<=g.naccu
                        t = ceil(rand([1 2])*trials); % random ints [1,g.timesout]
                        arg1(:,j) = oriarg1(:, index, t(1));
                        arg2(:,j) = oriarg2(:, index, t(2));
                        %arg2(:,j) = squeeze(oriarg2(t(2),index,:));
                        j = j+1;
                    end
                    eval([ formula ';' ]);
                    %g.Coherboot = cohercomp(g.Coherboot, tmpsX, tmpsY, 1, 1:g.naccu);
                end;
                if ~isempty(g.formulapost)
                    eval([ g.formulapost ';' ]);
                end;
                %g.Coherboot = cohercomppost(g.Coherboot);  % CHECK IF NECSSARY FOR ALL BOOT TYPE
                for index2= 1:length(g.formulaout) 
                    eval([ 'accarray' int2str(index2) '(:,:,index) = ' g.formulaout{index2} ';' ]);
                end;
            end;
        elseif strcmpi(g.boottype, 'both') % handle timestrials bootstrap
            fprintf('Bootstrap type is 2-D, shuffling along dimension 2 and 3\n');
            fprintf('Processing %s (naccu=%d) bootstrap (of %s %d):\n',bootname, g.naccu,unitname,trials);
            arg1 = zeros(nb_points, g.naccu);
            arg2 = zeros(nb_points, g.naccu );
            for allt=1:trials
                if rem(allt,10) == 0,  fprintf(' %d',allt); end
                if rem(allt,120) == 0, fprintf('\n'); end
                
                j=1;
                while j<=g.naccu
                    t = ceil(rand([1 2])*trials); % random ints [1,trials]
                    
                    s = ceil(rand([1 2])* length(g.basevect)); % random ints [1,times]
                    s = g.basevect(s);
					
                    arg1(:,j) = oriarg1(:,s(1),t(1));
                    arg2(:,j) = oriarg2(:,s(2),t(2));
                    j = j+1;
                end
                eval([ formula ';' ]);
            end
            if ~isempty(g.formulapost)
                eval([ g.formulapost ';' ]);
            end;
            for index= 1:length(g.formulaout) 
                eval([ 'accarray' int2str(index) ' = ' g.formulaout{index} ';' ]);
            end;
        elseif strcmpi(g.boottype, 'both2') % handle timestrials bootstrap, shuffle time only once
            fprintf('Bootstrap type is 2-D, shuffling along dimension 2 and 3 (once per accumulation for 3)\n');
            fprintf('Processing %s (naccu=%d) bootstrap (of %s %d):\n',bootname, g.naccu,unitname,trials);
            arg1 = zeros(nb_points, g.naccu);
            arg2 = zeros(nb_points, g.naccu );	
            
            s = ceil(rand([1 2])* length(g.basevect)); % random ints [1,times]
            s = g.basevect(s);
            
            for allt=1:trials
                if rem(allt,10) == 0,  fprintf(' %d',allt); end
                if rem(allt,120) == 0, fprintf('\n'); end
                
                j=1;
                while j<=g.naccu
                    t = ceil(rand([1 2])*trials); % random ints [1,trials]
                    arg1(:,j) = oriarg1(:,s(1),t(1));
                    arg2(:,j) = oriarg2(:,s(2),t(2));
                    j = j+1;
                end
                eval([ formula ';' ]);
            end
            if ~isempty(g.formulapost)
                eval([ g.formulapost ';' ]);
            end;
            for index= 1:length(g.formulaout) 
                eval([ 'accarray' int2str(index) ' = ' g.formulaout{index} ';' ]);
            end;
        elseif strcmpi(g.boottype, 'first') % boottype is 'times'
            fprintf('Bootstrap type is 2-D, shuffling along dimension 2 only\n');
            fprintf('Processing %s (naccu=%d) 1-D bootstrap (of %s %d):\n',bootname, g.naccu,unitname,trials);
            arg1 = zeros(nb_points, g.naccu);
            arg2 = zeros(nb_points, g.naccu );
            for allt=1:trials
                if rem(allt,10) == 0,  fprintf(' %d',allt); end
                if rem(allt,120) == 0, fprintf('\n'); end
                j=1;
                while j<=g.naccu
                    goodbasewins = g.basevect; 
                    ngdbasewins = length(goodbasewins);
                    s = ceil(rand([1 2])*ngdbasewins); % random ints [1,times]
                    s=goodbasewins(s);
					
                    arg1(:,j) = oriarg1(:,s(1),allt);
                    arg2(:,j) = oriarg2(:,s(2),allt);
                    j = j+1;
                end
                eval([ formula ';' ]);
            end
            if ~isempty(g.formulapost)
                eval([ g.formulapost ';' ]);
            end;
            for index= 1:length(g.formulaout) 
                eval([ 'accarray' int2str(index) ' = ' g.formulaout{index} ';' ]);
            end;
        end;
    end;
    for index= 1:length(g.formulaout) 
        eval( [ 'accarray = accarray' int2str(index) ';' ]);
        Rbootout{index} = accarray;
        
        % 'boottype'='times' or 'timestrials', size(R)=nb_points*naccu
        % 'boottype'='trials',                 size(R)=nb_points*naccu*times
        if ~isreal(accarray)
            accarray = sqrt(accarray .* conj(accarray)); % faster than abs()
        end;
        accarray = sort(accarray,2); % always sort on naccu (when 3D, naccu is the second dim)
        
        % compute bootstrap significance level
        i = round(g.naccu*g.alpha);
        %rsignif = mean(accarray(:,g.naccu-i+1:g.naccu),2); % significance levels for Rraw
        if strcmpi(g.bootside{min(length(g.bootside), index)}, 'upper');
            accarray        = squeeze(mean(accarray(:,g.naccu-i+1:g.naccu,:),2));
        else 
            if strcmpi(g.boottype,'second') & ndims(accarray) ==3
                accarraytmp        = squeeze(mean(accarray(:,1:i,:),2));
                accarraytmp(:,:,2) = squeeze(mean(accarray(:,g.naccu-i+1:g.naccu,:),2));
                accarray = accarraytmp;
            else
                accarray = [mean(accarray(:,1:i),2) mean(accarray(:,g.naccu-i+1:g.naccu),2)];
            end;
        end;
        accarrayout{index} = accarray;
    end;
    if length(Rbootout) == 1, Rbootout = Rbootout{1}; end;
    if length(accarrayout) == 1, accarrayout = accarrayout{1}; end;
end; % 2-D or 3-D

% fit to distribution (currently only for 1-D data)
if strcmpi(g.distfit, 'on')
    if length(g.vals) ~= 1
        error('For fitting, vals must contain exactly one value');
    end;
    
    % fitting with Ramberg-Schmeiser distribution
    % -------------------------------------------
    accarrayout = 1 - rsfit(abs(Rbootout(:)), g.vals);
    if ~isempty(g.correctp)
        if length(g.correctp) == 2
            accarrayout = correctfit(accarrayout, 'gamparams', [g.correctp 0]); % no correction for p=0
        else
            accarrayout = correctfit(accarrayout, 'gamparams', g.correctp);
        end;
    end;
    return;
    
    % fitting with normal distribution (deprecated)
    % --------------------------------
    [mu sigma] = normfit(abs(Rbootout(:)));
    accarrayout = 1 - normcdf(g.vals, mu, sigma); % cumulative density distribution
                                        % formula of normal distribution
                                        % y = 1/sqrt(2) * exp( -(x-mu).^2/(sigma*sigma*2) ) / (sqrt(pi)*sigma);
    % plot results 
    % -------------------------------------
    % figure;
    % hist(abs(Rbootout)); tmpax = axis;
    % hold on; 
    % valcomp = linspace(min(abs(Rbootout(:))), max(abs(Rbootout(:))), 100);
    % normy = normpdf(valcomp, mu, sigma);
    % plot(valcomp, normy/max(normy)*tmpax(4), 'r');
    % return;
end;
    
    
return;

% % Gamma and Beta fits:
% elseif strcmpi(g.distfit, 'gamma')
%  [phatgam pcigam] = gamfit(abs(Rbootout(:)));
%  gamy = gampdf(valcomp, phatgam(1), pcigam(2))
%  p = 1 - gamcdf(g.vals, phatgam(1), pcigam(2)); % cumulative density distribution
% elseif strcmpi(g.distfit, 'beta')
%  [phatbeta pcibeta] = betafit(abs(Rbootout(:)));
%  betay = betapdf(valcomp, phatbeta(1), pcibeta(1));
%  p = 1 - betacdf(g.vals, phatbeta(1), pcibeta(1)); % cumulative density distribution
% end


