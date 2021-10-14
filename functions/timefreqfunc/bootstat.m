% bootstat() - accumulate surrogate data to assess significance by permutation of some 
%              measure of two input variables. 
%
%              If 'distfit','on', fits the psd with a 4th-order polynomial using the 
%              data kurtosis, as in Ramberg, J.S., Tadikamalla, P.R., Dudewicz E.J., 
%              Mykkytka, E.F. "A probability distribution and its uses in fitting data." 
%              Technometrics, 21:201-214, 1979.
% Usage:
%            >> [rsignif,rboot] = bootstat( { arg1 arg2 ...}, formula, varargin ...);
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
%   'boottype '   - ['rand'|'shuffle']
%                   'rand'  = do not shuffle data. Only flip polarity randomly (for real 
%                             number) or phase (for complex numbers).
%                   'shuffle' = shuffle values of first argument (see two options below). 
%                   Default.
%   'shuffledim'  - [integer] indices of dimensions to shuffle. For instance, [1 2] will
%                   shuffle the first two dimensions. Default is to shuffle along 
%                   dimension 2.
%   'shufflemode' - ['swap'|'regular'] shuffle mode. Either swap dimensions (for instance
%                   swap rows then columns if dimension [1 2] are selected) or shuffle 
%                   in each dimension independently (slower). If only one dimension is 
%                   selected for shuffling, this option does not change the result.
%   'randmode'    - ['opposite'|'inverse'] randomize sign (or phase for complex number,
%                   or randomly set half the value to reference.
%   'alpha'       - [real] significance level (between 0 and 1) {default 0.05}.
%   'naccu'       - [integer] number of exemplars to accumulate {default 200}.
%   'bootside'    - ['both'|'upper'] side of the surrogate distribution to
%                   consider for significance. This parameter affects the size
%                   of the last dimension of the accumulation array ('accres') 
%                   (size is 2 for 'both' and 1 for 'upper') {default: 'both'}.
%   'basevect'    - [integer vector] time vector indices for baseline in second dimension.
%                   {default: all time points}.
%   'rboot'       - accumulation array (from a previous call). Allows faster 
%                   computation of the 'rsignif' output {default: none}.
%   'formulaout'  - [string] name of the computed variable {default: 'res'}.
%   'dimaccu'     - [integer] use dimension in result to accumulate data.
%                   For instance if the result array is size [60x50] and this value is 2,
%                   the function will consider than 50 times 60 value have been accumulated.
%
% Fitting distribution:
%   'distfit'     - ['on'|'off'] fit distribution with known function to compute more accurate 
%                   limits or exact p-value (see 'vals' option). The MATLAB statistical toolbox 
%                   is required. This option is currently implemented only for 1-D data.
%   'vals'        - [float array] significance values. 'alpha' is ignored and 
%                   rsignif returns the p-values. Requires 'distfit' (see above).
%                   This option currently implemented only for 1-D data.
%   'correctp'    - [phat pci zerofreq] parameters for correcting for a biased probability 
%                   distribution (requires 'distfit' above). See help of correctfit().
% Outputs: 
%    rsignif      - significance arrays. 2 values (low high) for each point (use
%                   'alpha' to change these limits).
%    rboot        - accumulated surrogate data values.
%
% Authors: Arnaud Delorme, Bhaktivedcanta Institute, Mumbai, India, Nov 2004
%
% See also: timef()

% NOTE: There is an undocumented parameter, 'savecoher', [0|1]
% HELP TEXT REMOVED:  (Ex: Using option 'both', coherence during baseline would be 
%                      ignored since times are shuffled during each accumulation.

% Copyright (C) 9/2002  Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD
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

% *************************************
% To fit the psd with as 4th order polynomial using the distribution kurtosis,
% Reference: Ramberg, J.S., Tadikamalla, P.R., Dudewicz E.J., Mykkytka, E.F. 
% "A probability distribution and its uses in fitting data." 
% Technimetrics, 1979, 21: 201-214.
% *************************************

function [accarrayout, Rbootout, Rbootout2] = bootstat(oriargs, formula, varargin)
%	nb_points, timesout, naccu, baselength, baseboot, boottype, alpha, rboot);
	
if nargin < 2
	help bootstat;
	return;
end

if ~ischar(formula)
	error('The second argument must be a string formula');
end

g = finputcheck(varargin, ...
                { 'dims'          'integer'  []                       []; ...
                  'naccu'         'integer'  [0 10000]                200; ...
                  'bootside'      'string'   { 'both','upper' }       'both'; ...
                  'basevect'      'integer'  []                       []; ...
				  'boottype'      'string'  { 'rand','shuffle' }      'shuffle'; ...
				  'shufflemode'   'string'  { 'swap','regular' }      'swap'; ...
				  'randmode'      'string'  { 'opposite','inverse' }  'opposite'; ...
				  'shuffledim'    'integer'  [0 Inf]                  []; ...
				  'label'         'string'   []                       formula; ...
				  'alpha'         'real'     [0 1]                    0.05; ...
				  'vals'          'real'     []                       []; ...
				  'distfit'       'string'   {'on','off' }            'off'; ...
				  'dimaccu'       'integer'  [1 Inf]                  []; ...
				  'correctp'      'real'     []                       []; ...
				  'rboot'         'real'     []                       NaN	});
if ischar(g)
	error(g);
end
if isempty(g.shuffledim) && strcmpi(g.boottype, 'rand')
    g.shuffledim = []; 
elseif isempty(g.shuffledim)
    g.shuffledim = 2;
end; 
unitname = '';
if 2/g.alpha > g.naccu 
    if strcmpi(g.distfit, 'off') || ~((size(oriarg1,1) == 1 || size(oriarg1,2) == 1) && size(oriarg1,3) == 1)
        g.naccu = 2/g.alpha; 
        fprintf('Adjusting naccu to compute alpha value');
    end
end
if isempty(g.rboot)
	g.rboot = NaN;
end

% function for bootstrap computation
% ----------------------------------
if ~iscell(oriargs) || length(oriargs) == 1, 
    oriarg1 = oriargs;
    oriarg2 = []; 
else 
    oriarg1 = oriargs{1};
    oriarg2 = oriargs{2};
end
[nb_points times trials] = size(oriarg1);
if times == 1, disp('Warning 1 value only for shuffling dimension'); end

% only consider baseline
% ----------------------
if ~isempty(g.basevect)
    fprintf('\nPermutation statistics baseline length is %d (out of %d) points\n', length(g.basevect), times);
    arg1 = oriarg1(:,g.basevect,:);
    if ~isempty(oriarg2)
        arg2 = oriarg2(:,g.basevect,:);
    end
else
    arg1 = oriarg1;
    arg2 = oriarg2;
end

% formula for accumulation array
% ------------------------------
% if g.dimaccu is not empty, accumulate over that dimension
% of the resulting array to speed up computation
formula = [ 'res=' formula ];
g.formulapost = [ 'if index == 1, ' ...
                  '   if ~isempty(g.dimaccu), ' ...
                  '      Rbootout= zeros([ ceil(g.naccu/size(res,g.dimaccu)) size( res ) ]);' ...
                  '   else,' ...
                  '      Rbootout= zeros([ g.naccu size( res ) ]);' ...
                  '   end;' ...
                  'end,' ...
                  'Rbootout(count,:,:,:) = res;' ...
                  'count = count+1;' ...
                  'if ~isempty(g.dimaccu), ' ...
                  '   index = index + size(res,g.dimaccu);' ...
                  '   fprintf(''%d '', index-1);' ...
                  'else ' ...
                  '   index=index+1;' ...
                  '   if rem(index,10)  == 0, fprintf(''%d '', index); end;' ...
                  '   if rem(index,100) == 0, fprintf(''\n''); end;' ...
                  'end;' ];

% **************************
% case 1: precomputed values
% **************************
if ~isnan(g.rboot)
    Rbootout = g.rboot;
% ***********************************
% case 2: randomize polarity or phase
% ***********************************
elseif strcmpi(g.boottype, 'rand') && strcmpi(g.randmode, 'inverse')
    fprintf('Bootstat function: randomize inverse values\n');
    fprintf('Processing permutation statistics for %s (naccu=%d):', g.label, g.naccu);
    
    % compute random array
    % --------------------
    multarray = ones(size(arg1));
    totlen    = prod(size(arg1));
    if isreal(arg1), 
        multarray(1:round(totlen/2)) = 0;
    end
    for shuff = 1:ndims(multarray)
         multarray = supershuffle(multarray,shuff); % initial shuffling
    end
    if isempty(g.shuffledim), g.shuffledim = 1:ndims(multarray); end
    invarg1 = 1./arg1;
    
    % accumulate
    % ----------
    index = 1;
    count = 1;
    while index <= g.naccu
        for shuff = g.shuffledim
            multarray = supershuffle(multarray,shuff);
        end
        tmpinds = find(reshape(multarray, 1, prod(size(multarray))));
        arg1 = oriarg1;
        arg1(tmpinds) = invarg1(tmpinds);
        eval([ formula ';' ]);
        eval( g.formulapost ); % also contains index = index+1
    end
elseif strcmpi(g.boottype, 'rand') % opposite
    fprintf('Bootstat function: randomize polarity or phase\n');
    fprintf('Processing permutation statistics for %s (naccu=%d):', g.label, g.naccu);
    
    % compute random array
    % --------------------
    multarray = ones(size(arg1));
    totlen    = prod(size(arg1));
    if isreal(arg1), 
        multarray(1:round(totlen/2)) = -1;
    else
        tmparray            = exp(j*linspace(0,2*pi,totlen+1));
        multarray(1:totlen) = tmparray(1:end-1);
    end
    for shuff = 1:ndims(multarray)
         multarray = supershuffle(multarray,shuff); % initial shuffling
    end
    if isempty(g.shuffledim), g.shuffledim = 1:ndims(multarray); end
    
    % accumulate
    % ----------
    index = 1;
    count = 1;
    while index <= g.naccu
        for shuff = g.shuffledim
            multarray = supershuffle(multarray,shuff);
        end
        arg1 = arg1.*multarray;
        eval([ formula ';' ]);
        eval( g.formulapost ); % also contains index = index+1
    end
% ********************************************
% case 3: shuffle vector of only one dimension
% ********************************************
elseif length(g.shuffledim) == 1
    fprintf('Bootstat function: shuffling along dimension %d only\n', g.shuffledim);
    fprintf('Processing permutation statistics for %s (naccu=%d):', g.label, g.naccu);

    index = 1;
    count = 1;
    while index <= g.naccu
        arg1 = shuffleonedim(arg1,g.shuffledim);
        eval([ formula ';' ]);
        eval( g.formulapost );
    end
% ***********************************************
% case 5: shuffle vector along several dimensions
% ***********************************************
else 
    if strcmpi(g.shufflemode, 'swap') % swap mode
        fprintf('Bootstat function: shuffling along dimension %s (swap mode)\n', int2str(g.shuffledim));
        fprintf('Processing permutation statistics for %s (naccu=%d):', g.label, g.naccu);
        index = 1;
        count = 1;
        while index <= g.naccu
            for shuff = g.shuffledim
                arg1 = supershuffle(arg1,shuff);
            end
            eval([ formula ';' ]);
            eval( g.formulapost );
        end
    else  % regular shuffling
        fprintf('Bootstat function: shuffling along dimension %s (regular mode)\n', int2str(g.shuffledim));
        fprintf('Processing permutation statistics for %s (naccu=%d):', g.label, g.naccu);
        index = 1;
        count = 1;
        while index <= g.naccu
            for shuff = g.shuffledim
                arg1 = shuffleonedim(arg1,shuff);
            end
            eval([ formula ';' ]);
            eval( g.formulapost );
        end
    end
end
Rbootout(count:end,:,:,:) = [];

% **********************
% assessing significance
% **********************

% get accumulation array
% ----------------------
accarray = Rbootout;
if ~isreal(accarray)
    accarray = sqrt(accarray .* conj(accarray)); % faster than abs()
end
% reshape the output if necessary
% -------------------------------
if ~isempty(g.dimaccu)
    if g.dimaccu+1 == 3
        accarray = permute( accarray, [1 3 2]);
    end
    accarray = reshape( accarray, size(accarray,1)*size(accarray,2), size(accarray,3) );
end
if size(accarray,1) == 1, accarray = accarray'; end; % first dim contains g.naccu

% ******************************************************
% compute thresholds on array not fitting a distribution
% ******************************************************
if strcmpi(g.distfit, 'off')
  
    % compute bootstrap significance level
    % ------------------------------------
    accarray  = sort(accarray,1); % always sort on naccu
    Rbootout2 = accarray;
    i         = round(size(accarray,1)*g.alpha);
    accarray1 = squeeze(mean(accarray(size(accarray,1)-i+1:end,:,:),1));
    accarray2 = squeeze(mean(accarray(1:i                     ,:,:),1));
    if abs(accarray(1,1,1) - accarray(end,1,1)) < abs(accarray(1,1,1))*1e-15
        accarray1(:) = NaN;
        accarray2(:) = NaN;
    end

else
    % *******************
    % fit to distribution 
    % *******************
    sizerboot   = size (accarray);
    accarray1   = zeros(sizerboot(2:end));
    accarray2   = zeros(sizerboot(2:end));
    
    if ~isempty(g.vals{index})
        if ~all(size(g.vals{index}) == sizerboot(2:end) )
            error('For fitting, vals must have the same dimension as the output array (try transposing)');
        end
    end
    
    % fitting with Ramberg-Schmeiser distribution
    % -------------------------------------------
    if ~isempty(g.vals{index}) % compute significance for value
        for index1 = 1:size(accarrayout,1)
            for index2 = 1:size(accarrayout,2)
                accarray1(index1,index2) = 1 - rsfit(squeeze(accarray(:,index1,index2)), g.vals{index}(index1, index2));
                if length(g.correctp) == 2
                    accarray1(index1,index2) = correctfit(accarray1, 'gamparams', [g.correctp 0]); % no correction for p=0
                else
                    accarray1(index1,index2) = correctfit(accarray1, 'gamparams', g.correctp);
                end
            end
        end
    else % compute value for significance
        for index1 = 1:size(accarrayout,1)
            for index2 = 1:size(accarrayout,2)
                [p c l chi2] = rsfit(Rbootout(:),0);
                pval = g.alpha;   accarray1(index1,index2) = l(1) + (pval.^l(3) - (1-pval).^l(4))/l(2);
                pval = 1-g.alpha; accarray2(index1,index2) = l(1) + (pval.^l(3) - (1-pval).^l(4))/l(2);        
            end
        end
    end

    % plot results 
    % -------------------------------------
    % figure;
    % hist(abs(Rbootout)); tmpax = axis;
    % hold on; 
    % valcomp = linspace(min(abs(Rbootout(:))), max(abs(Rbootout(:))), 100);
    % normy = normpdf(valcomp, mu, sigma);
    % plot(valcomp, normy/max(normy)*tmpax(4), 'r');
    % return;
end

% set output array: backward compatible
% -------------------------------------
if strcmpi(g.bootside, 'upper'); % only upper significance
    accarrayout = accarray1;
else 
    if size(accarray1,1) ~= 1 && size(accarray1,2) ~= 1
        accarrayout        = accarray2;
        accarrayout(:,:,2) = accarray1;
    else
        accarrayout = [ accarray2(:) accarray1(:) ];
    end
end
accarrayout = squeeze(accarrayout);
if size(accarrayout,1) == 1 && size(accarrayout,3) == 1, accarrayout = accarrayout'; end

% better but not backward compatible
% ----------------------------------
% accarrayout = { accarray1 accarray2 }; 
    
return;

    % fitting with normal distribution (deprecated)
    % --------------------------------
    [mu sigma] = normfit(abs(Rbootout(:)));
    accarrayout = 1 - normcdf(g.vals, mu, sigma); % cumulative density distribution
                                        % formula of normal distribution
                                        % y = 1/sqrt(2) * exp( -(x-mu).^2/(sigma*sigma*2) ) / (sqrt(pi)*sigma);
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


    if strcmpi(g.distfit, 'off')
        tmpsort = sort(Rbootout);
        i = round(g.alpha*g.naccu);
        sigval = [mean(tmpsort(1:i)) mean(tmpsort(g.naccu-i+1:g.naccu))];
        if strcmpi(g.bootside, 'upper'), sigval = sigval(2); end
        accarrayout = sigval;
    end

% this shuffling preserve the number of -1 and 1
% for cloumns and rows (assuming matrix size is multiple of 2
% -----------------------------------------------------------
function array = supershuffle(array, dim)
    if size(array, 1) == 1 || size(array,2) == 1
        array = shuffle(array);
        return;
    end
    if size(array, dim) == 1, return; end

    if dim == 1
        indrows = shuffle(1:size(array,1));
        for index = 1:2:length(indrows)-rem(length(indrows),2) % shuffle rows
            tmparray                    = array(indrows(index),:,:);
            array(indrows(index),:,:)   = array(indrows(index+1),:,:);
            array(indrows(index+1),:,:) = tmparray;
        end
    elseif dim == 2
        indcols = shuffle(1:size(array,2));
        for index = 1:2:length(indcols)-rem(length(indcols),2) % shuffle columns
            tmparray                    = array(:,indcols(index),:);
            array(:,indcols(index),:)   = array(:,indcols(index+1),:);
            array(:,indcols(index+1),:) = tmparray;
        end
    else
        ind3d = shuffle(1:size(array,3));
        for index = 1:2:length(ind3d)-rem(length(ind3d),2) % shuffle columns
            tmparray                  = array(:,:,ind3d(index));
            array(:,:,ind3d(index))   = array(:,:,ind3d(index+1));
            array(:,:,ind3d(index+1)) = tmparray;
        end
    end

% shuffle one dimension, one row/columns at a time
% -----------------------------------------------
function array = shuffleonedim(array, dim)
    if size(array, 1) == 1 || size(array,2) == 1
        array = shuffle(array, dim);
    else
        if dim == 1
            for index1 = 1:size(array,3)
                for index2 = 1:size(array,2)
                    array(:,index2,index1) = shuffle(array(:,index2,index1));
                end
            end
        elseif dim == 2
            for index1 = 1:size(array,3)
                for index2 = 1:size(array,1)
                    array(index2,:,index1) = shuffle(array(index2,:,index1));
                end
            end
        else
            for index1 = 1:size(array,1)
                for index2 = 1:size(array,2)
                    array(index1,index2,:) = shuffle(array(index1,index2,:));
                end
            end
        end
    end
