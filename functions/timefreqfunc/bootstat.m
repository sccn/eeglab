% bootstat() - accumulate surrogate data for comparing two conditions 
%
% Usage:
%     >> [accarray, rsignif] = bootstat(arg1, arg2, formula, varargin ...);
%
% Inputs:
%    arg1    - [array] 2-D or 3-D array of value
%    arg2    - [array] 2-D or 3-D array of value
%    formula - [string] formula to compute a given measure. Takes arguments
%              'arg1', 'arg2' as inputs and 'res' (by default) as output. i.e.
%              'res = res + arg1 .* conj(arg2)'
%
% Optional inputs:
%   'boottype '   - ['first'|'second'|'both'] 'fist' accumulate in the first 
%                 dimension only (shuffling the first dimension). 'second'
%                 shuffle the second dimension but keep the fist constant. 
%                 'both' shuffle the two dimension. Default is 'both'.
%   'alpha'       - [real] significance (between 0 and 1). Default is 0.05.
%   'naccu'       - [integer] number of exemplar to accumulate. Default is 200.
%   'bootside'    - ['both'|'upper'] side of the surrogate distribution to
%                 consider for significance. This parameter affect the size
%                 of the last dimension of accumulation array 'accres' (size
%                 is 2 for 'both' and 1 for 'upper'). Default is 'both'.
%   'basevect'    - [integer vector] time vector indices for baseline. Default 
%                 is all time points.
%   'accarray'    - accumulation array (from previous calls). Allow to compute
%                 only the 'rsignif' output.
%   'formulainit' - [string] for initializing variable. i.e. 'res = zeros(10,40);'
%                 Default is initialization of the res variable to
%                 (size(arg1,3) x naccu)
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

% NOTE: one hidden parameter 'savecoher', 0 or 1

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
% Revision 1.1  2002/10/01 16:07:33  arno
% Initial revision
%
% Revision 1.1  2002/09/24 23:28:02  arno
% Initial revision
%
% function for bootstrap initialisation
% -------------------------------------

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
				  'boottype'      'string'  {'first' 'second' 'both' 'times' 'trials' 'timestrials'} 'both'; ...
				  'alpha'         'real'     [0 1]                    0.05; ...
				  'accarray'      'integer'  []                       NaN	});
if isstr(Boot)
	error(Boot);
end;
bootname = Boot.boottype;
unitname = '';
switch Boot.boottype
 case 'times',  Boot.boottype = 'first'; unitname = 'trial '; bootname = 'times';
 case 'trials', Boot.boottype = 'second'; unitname = 'time '; bootname = 'trials';
 case 'timestrials',  Boot.boottype = 'both'; unitname = 'trial '; bootname = 'timestrials';
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

if isnan(Boot.accarray)
	eval( Boot.formulainit );
	if strcmpi( Boot.boottype, 'second') % get g.naccu bootstrap estimates for each trial
		fprintf('\nProcessing %s (naccu=%d) bootstrap (of %s%d):',bootname, Boot.naccu, unitname, times(end));
		
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
		fprintf('\nProcessing %s (naccu=%d) bootstrap (of %s%d):',bootname, Boot.naccu,unitname,trials);
        arg1 = zeros(nb_points, Boot.naccu);
		arg2 = zeros(nb_points, Boot.naccu );
		for allt=1:trials
			if rem(allt,10) == 0,  fprintf(' %d',allt); end
			if rem(allt,120) == 0, fprintf('\n'); end
			
            j=1;
			while j<=Boot.naccu
				t = ceil(rand([1 2])*trials); % random ints [1,trials]
				
				goodbasewins = Boot.basevect; 
				ngdbasewins = length(goodbasewins);
				s = ceil(rand([1 2])*ngdbasewins); % random ints [1,times]
				s=goodbasewins(s);
					
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
		fprintf('\nProcessing %s (naccu=%d) bootstrap (of %s%d):',bootname, Boot.naccu,unitname,trials);
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
