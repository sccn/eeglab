% bootstat() - accumulate surrogate data to compare two conditions 
%
% Usage:
%     >> [accarray, rsignif] = bootstat(arg1, arg2, formula, varargin ...);
%
% Inputs:
%    arg1    - [array] 2-D or 3-D array of values
%    arg2    - [array] 2-D or 3-D array of values
%    formula - [string] formula to compute the given measure. Takes arguments
%              'arg1', 'arg2' as inputs and 'res' (result, by default) as 
%              output. i.e.    'res = res + arg1 .* conj(arg2)'
%
% Optional inputs:
%   'boottype '   - ['first'|'second'|'both'|'both2'] 
%                 'first'=  accumulate in the first dimension only (shuffling 
%                           the first dimension (e.g., trials)). 
%                 'second'= shuffle the second dimension (e.g., times), 
%                           keeping the first (e.g., trials) constant. 
%                 'both'=   shuffle both dimensions repetivelly during accumulation.
%                 'both2'=  same as 'both' but shuffle dimension 2 only once per 
%                           accumulation. In the context of EEG coherence, this 
%                           option allows to detect significant changes with respect
%                           to baseline synchronization (using option 'both' 
%                           synchronizations during baseline are ignored since time
%                           is shuffled during accumulation.
%                 Default is 'both'.
%   'alpha'       - [real] significance level (between 0 and 1). Default is 0.05.
%   'naccu'       - [integer] number of exemplars to accumulate. Default is 200.
%   'bootside'    - ['both'|'upper'] side of the surrogate distribution to
%                 consider for significance. This parameter affects the size
%                 of the last dimension of the accumulation array 'accres' (size
%                 is 2 for 'both' and 1 for 'upper'). Default is 'both'.
%   'basevect'    - [integer vector] time vector indices for baseline. Default 
%                 is all time points.
%   'accarray'    - accumulation array (from previous calls). Allows computing
%                 only the 'rsignif' output.
%   'formulainit' - [string] for initializing variable. i.e. 'res = zeros(10,40);'
%                 Default is initializing  'res = (size(arg1,3) x naccu)'
%   'formulapost' - [string] after the accumulation. i.e. 'res = res /10;'
%                 default is none.
%   'formulaout'  - [string] name of the computed variable. Default is 'res'.
%
% Outputs: 
%    accres  - result for shuffled data
%    res1    - result for first condition
%    res2    - result for second condition
%
% Authors: Arnaud Delorme, Lars & Scott Makeig
%          CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, La Jolla, 2002-
%
% See also: timef()

% NOTE: There is one hidden parameter 'savecoher', 0 or 1

% Copyright (C) 8/1/98  Arnaud Delorme, Sigurd Enghoff & Scott Makeig, SCCN/INC/UCSD
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
% to do: implement the fitting of the psd with as 4th order
% polynomial using the kurtosis of the distribution
% See: Ramberg, J.S., Tadikamalla, P.R., Dudewicz E.J., Mykkytka, E.F. 
% A probability distribution and its uses in fitting data. Technimetrics, 
% 1979, 21: 201-214, probably avalaible at the Script Librairy
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

Boot = finputcheck(varargin, ...
                { 'dims'          'integer'  []                       []; ...
                  'naccu'         'integer'  [0 10000]                200; ...
                  'bootside'      {'cell' 'string'}   {'both' 'upper'}         'both'; ...
                  'basevect'      'integer'  []                       []; ...
                  'formulapost'   'string'   []                       ''; ...
                  'formulainit'   'string'   []                       'res = zeros(nb_points,Boot.naccu);'; ...
                  'formulaout'    {'cell' 'string'} []                       'res'; ...
				  'boottype'      'string'  {'first' 'second' 'both' 'times' 'trials' 'timestrials' 'timestrials2' 'both2'} 'both'; ...
				  'alpha'         'real'     [0 1]                    0.05; ...
				  'accarray'      'integer'  []                       NaN	});
if isstr(Boot)
	error(Boot);
end;
bootname = Boot.boottype;
unitname = '';
switch Boot.boottype
 case 'times',  Boot.boottype = 'first'; unitname = 'trials'; bootname = 'times';
 case 'trials', Boot.boottype = 'second'; unitname = 'times'; bootname = 'trials';
 case 'timestrials',  Boot.boottype = 'both'; unitname = 'trials'; bootname = 'timestrials';
 case 'timestrials2', Boot.boottype = 'both2'; unitname = 'trials'; bootname = 'timestrials';
end;
if ~iscell(Boot.formulaout)
    Boot.formulaout = { Boot.formulaout };
end;
if ~iscell(Boot.bootside)
    Boot.bootside = { Boot.bootside };
end;    
if isempty(Boot.accarray)
	Boot.accarray = NaN;
end;

%   Boot.Rboot       = zeros(nb_points,naccu);  % summed bootstrap coher
%	Boot.boottype    = boottype;
%	Boot.baselength  = baselength;
%	Boot.baseboot    = baseboot;
%	Boot.Coherboot   = Coherboot;
%	Boot.naccu       = naccu;
%	Boot.alpha       = alpha;
%	Boot.rboot       = rboot;

% function for bootstrap computation
% ----------------------------------
[nb_points times trials]    = size(oriarg1);

%trials    = size(oriarg1, 1);
%times     = size(oriarg1, 2);
%nb_points = size(oriarg1, 3);

if isempty(Boot.basevect)
	Boot.basevect = 1:times;
end;

fprintf('\nBootstrap baseline length is %d (out of %d) points\n', length(Boot.basevect), times);

if isnan(Boot.accarray)
	eval( Boot.formulainit );
	if strcmpi( Boot.boottype, 'second') % get g.naccu bootstrap estimates for each time window
        fprintf('Bootstrap type is 3-D, shuffling oalong dimension 3 only.\n');
		fprintf('Processing %s (naccu=%d) bootstrap (of %s %d):\n',bootname, Boot.naccu, unitname, times(end));
		
		arg2 = zeros(nb_points, Boot.naccu);
		arg1 = zeros(nb_points, Boot.naccu);
        for index= 1:length(Boot.formulaout) 
            eval( ['accarray' int2str(index) '= zeros(nb_points, Boot.naccu, times);'] );
        end;
		
		for index=1:times % dim2
			if rem(index,10) == 0,  fprintf(' %d',index); end
			if rem(index,120) == 0, fprintf('\n'); end
			for allt=1:trials % dim1
				j=1;
				while j<=Boot.naccu
					t = ceil(rand([1 2])*trials); % random ints [1,g.timesout]
					arg1(:,j) = oriarg1(:, index, t(1));
					arg2(:,j) = oriarg2(:, index, t(2));
					%arg2(:,j) = squeeze(oriarg2(t(2),index,:));
					j = j+1;
				end
				eval([ formula ';' ]);
				%Boot.Coherboot = cohercomp(Boot.Coherboot, tmpsX, tmpsY, 1, 1:Boot.naccu);
			end;
			if ~isempty(Boot.formulapost)
				eval([ Boot.formulapost ';' ]);
			end;
			%Boot.Coherboot = cohercomppost(Boot.Coherboot);  % CHECK IF NECSSARY FOR ALL BOOT TYPE
            for index2= 1:length(Boot.formulaout) 
                eval([ 'accarray' int2str(index2) '(:,:,index) = ' Boot.formulaout{index2} ';' ]);
            end;
		end;
	elseif strcmpi(Boot.boottype, 'both') % handle timestrials bootstrap
        fprintf('Bootstrap type is 2-D, shuffling along dimension 2 and 3\n');
		fprintf('Processing %s (naccu=%d) bootstrap (of %s %d):\n',bootname, Boot.naccu,unitname,trials);
        arg1 = zeros(nb_points, Boot.naccu);
		arg2 = zeros(nb_points, Boot.naccu );
		for allt=1:trials
			if rem(allt,10) == 0,  fprintf(' %d',allt); end
			if rem(allt,120) == 0, fprintf('\n'); end
			
            j=1;
			while j<=Boot.naccu
				t = ceil(rand([1 2])*trials); % random ints [1,trials]
				
				s = ceil(rand([1 2])* length(Boot.basevect)); % random ints [1,times]
				s = Boot.basevect(s);
					
				arg1(:,j) = oriarg1(:,s(1),t(1));
				arg2(:,j) = oriarg2(:,s(2),t(2));
				j = j+1;
			end
			eval([ formula ';' ]);
		end
		if ~isempty(Boot.formulapost)
			eval([ Boot.formulapost ';' ]);
		end;
        for index= 1:length(Boot.formulaout) 
            eval([ 'accarray' int2str(index) ' = ' Boot.formulaout{index} ';' ]);
        end;
	elseif strcmpi(Boot.boottype, 'both2') % handle timestrials bootstrap, shuffle time only once
        fprintf('Bootstrap type is 2-D, shuffling along dimension 2 and 3 (once per accumulation for 3)\n');
		fprintf('Processing %s (naccu=%d) bootstrap (of %s %d):\n',bootname, Boot.naccu,unitname,trials);
        arg1 = zeros(nb_points, Boot.naccu);
		arg2 = zeros(nb_points, Boot.naccu );	
        
        s = ceil(rand([1 2])* length(Boot.basevect)); % random ints [1,times]
        s = Boot.basevect(s);
        
        for allt=1:trials
			if rem(allt,10) == 0,  fprintf(' %d',allt); end
			if rem(allt,120) == 0, fprintf('\n'); end
			
            j=1;
			while j<=Boot.naccu
				t = ceil(rand([1 2])*trials); % random ints [1,trials]
				arg1(:,j) = oriarg1(:,s(1),t(1));
				arg2(:,j) = oriarg2(:,s(2),t(2));
				j = j+1;
			end
			eval([ formula ';' ]);
		end
		if ~isempty(Boot.formulapost)
			eval([ Boot.formulapost ';' ]);
		end;
        for index= 1:length(Boot.formulaout) 
            eval([ 'accarray' int2str(index) ' = ' Boot.formulaout{index} ';' ]);
        end;
	elseif strcmpi(Boot.boottype, 'first') % boottype is 'times'
        fprintf('Bootstrap type is 2-D, shuffling along dimension 2 only\n');
		fprintf('Processing %s (naccu=%d) 1-D bootstrap (of %s%d):\n',bootname, Boot.naccu,unitname,trials);
		arg1 = zeros(nb_points, Boot.naccu);
		arg2 = zeros(nb_points, Boot.naccu );
		for allt=1:trials
			if rem(allt,10) == 0,  fprintf(' %d',allt); end
			if rem(allt,120) == 0, fprintf('\n'); end
			j=1;
			while j<=Boot.naccu
				goodbasewins = Boot.basevect; 
				ngdbasewins = length(goodbasewins);
				s = ceil(rand([1 2])*ngdbasewins); % random ints [1,times]
				s=goodbasewins(s);
					
				arg1(:,j) = oriarg1(:,s(1),allt);
				arg2(:,j) = oriarg2(:,s(2),allt);
				j = j+1;
			end
			eval([ formula ';' ]);
		end
		if ~isempty(Boot.formulapost)
			eval([ Boot.formulapost ';' ]);
		end;
        for index= 1:length(Boot.formulaout) 
            eval([ 'accarray' int2str(index) ' = ' Boot.formulaout{index} ';' ]);
        end;
	end;
end;
for index= 1:length(Boot.formulaout) 
    eval( [ 'accarray = accarray' int2str(index) ';' ]);
    Rbootout{index} = accarray;
	
    % 'boottype'='times' or 'timestrials', size(R)=nb_points*naccu
    % 'boottype'='trials',                 size(R)=nb_points*naccu*times
    if ~isreal(accarray)
		accarray = sqrt(accarray .* conj(accarray)); % faster than abs()
    end;
	accarray = sort(accarray,2); % always sort on naccu (when 3D, naccu is the second dim)
    
    % compute bootstrap significance level
    i = round(Boot.naccu*Boot.alpha);
    %rsignif = mean(accarray(:,Boot.naccu-i+1:Boot.naccu),2); % significance levels for Rraw
    if strcmpi(Boot.bootside{min(length(Boot.bootside), index)}, 'upper');
        accarray        = squeeze(mean(accarray(:,Boot.naccu-i+1:Boot.naccu,:),2));
	else 
        if strcmpi(Boot.boottype,'second') & ndims(accarray) ==3
            accarraytmp        = squeeze(mean(accarray(:,1:i,:),2));
            accarraytmp(:,:,2) = squeeze(mean(accarray(:,Boot.naccu-i+1:Boot.naccu,:),2));
            accarray = accarraytmp;
        else
            accarray = [mean(accarray(:,1:i),2) mean(accarray(:,Boot.naccu-i+1:Boot.naccu),2)];
        end;
    end;
    accarrayout{index} = accarray;
end;

if length(Rbootout) == 1, Rbootout = Rbootout{1}; end;
if length(accarrayout) == 1, accarrayout = accarrayout{1}; end;
