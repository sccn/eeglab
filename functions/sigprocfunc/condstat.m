% condstat() - accumulate surrogate data for comparing two conditions 
%
% Usage:
%     >> [diffres, accres, res1, res2] = condstat(formula, naccu, alpha, ...
%                                              bootside, condboot, arg1, arg2 ...);
%
% Inputs:
%    formula  - [string] formula to compute a given measure. Takes arguments
%               'arg1', 'arg2' ... as inputs. i.e.
%               'sum(arg1(:,:,X),3) ./ sqrt(sum(tfx(:,:,X))) ./ sqrt(sum(tfy(:,:,X)))'
%    naccu    - [integer] number of accumulation. i.e. 200
%    alpha    - [float] significance level (0<alpha<0.5)
%    bootside - ['both'|'upper'] side of the surrogate distribution to
%               consider for significance. This parameter affect the size
%               of the last dimension of accumulation array 'accres' (size
%               is 2 for 'both' and 1 for 'upper').
%    bootside - ['both'|'upper'] 
%    condboot - ['abs'|'angle'|'complex'|''] for comparing 2 conditions,
%               either absolute vales ('abs'), angles ('angles') or complex values 
%               ('complex'). '' and 'complex' let the formula unchanged.
%    arg1     - [cell_array] of 2 nD array of values to compare. The last dimensions
%               of the array is the dimention that will be shuffled to accumulate
%               data.
%    arg2...  - same as arg1
%
% Outputs: 
%    diffres  - differential array computed on the non-shuffled data
%    accdres  - result for shuffled data
%    res1     - result for first condition
%    res2     - result for second condition
%
% Authors: Arnaud Delorme, Lars & Scott Makeig
%          CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, La Jolla, 2002-
%
% See also: timef(), crossf()

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
% Revision 1.2  2002/10/01 16:06:16  arno
% compute statistics now
%
% Revision 1.1  2002/09/24 23:28:02  arno
% Initial revision
%

function [diffres, accres, res1, res2] = condstat(formula, naccu, alpha, bootside, condboot, varargin);

if nargin < 2
	help condstat;
	return;
end;

if ~isstr(formula) & ~iscell(formula)
	error('The first argument must be a string formula or cell array of string');
end;
if isstr(formula)
	formula = { formula };
end;
if isstr(bootside)
	bootside = { bootside };
end;
for index = 1:length(bootside)
	if ~strcmpi(bootside, 'both') & ~strcmpi(bootside, 'upper')
		error('Bootside must be either ''both'' or ''upper''');
	end;
end;	
if isstr(condboot)
	condboot = { condboot };
end;
for index = 1:length(condboot)
	if isempty(condboot{index}), condboot{index} = 'complex'; end;
end;	
		
for index = 1:length(varargin)
	if ~iscell(varargin) | length(varargin{index}) ~=2
		error('Except for the first arguments, all other arguments given to the function must be cell arrays of two numerical array');
	end;
end;

% accumulate coherence images (all arrays [nb_points x timesout x trials])
% ---------------------------
for index=1:length(varargin)
	tmpvar1 = varargin{index}{1};
	tmpvar2 = varargin{index}{2};
	if index == 1
		cond1trials = size(tmpvar1,ndims(tmpvar1));
		cond2trials = size(tmpvar2,ndims(tmpvar2));
		for tmpi = 1:length(formula)
			accres{tmpi} = zeros(size(tmpvar1,1), size(tmpvar1,2), naccu);
		end;
	end;
	
	if ndims(tmpvar1) == 2
		eval( [ 'arg' int2str(index) '=zeros(size(tmpvar1,1), cond1trials+cond2trials);' ] );
		eval( [ 'arg' int2str(index) '(:,1:cond1trials)=tmpvar1;' ] );
		eval( [ 'arg' int2str(index) '(:,cond1trials1+1:end)=tmpvar2;' ] );
	else
		eval( [ 'arg' int2str(index) '=zeros(size(tmpvar1,1), size(tmpvar1,2), cond1trials+cond2trials);' ] );
		eval( [ 'arg' int2str(index) '(:,:,1:cond1trials)=tmpvar1;' ] );
		eval( [ 'arg' int2str(index) '(:,:,cond1trials+1:end)=tmpvar2;' ] );
	end;
end;

fprintf('Accumulating bootstrap:');
alltrials = [1:cond1trials+cond2trials];

% processing formulas
% -------------------
formula1 = [];
formula2 = [];
for index = 1:length(formula)	
	% updating formula
	% ----------------
	switch lower(condboot{ min(length(condboot), index) })
	 case 'abs', formula{index} = [ 'abs(' formula{index} ')' ];
	 %case 'angle', formula{index} = [ 'exp(j*angle(' formula{index} '))' ];
	 case 'angle', formula{index} = [ 'angle(' formula{index} ')/(2*pi)' ];
	 case 'complex', ;
	 otherwise, condboot, error('condstat argument must be either ''abs'', ''angle'', ''complex'' or empty');
	end;

	% computing difference (non-shuffled)
	% -----------------------------------
	X = 1:cond1trials;
	eval( [ 'res1{index} = ' formula{index} ';'] );
	X = cond1trials+1:cond1trials+cond2trials;
	eval( [ 'res2{index} = ' formula{index}  ';'] );
	diffres{index} = res1{index}-res2{index};

	% build formula to execute
	% ------------------------
	arrayname = [ 'accres{' int2str(index) '}' ];
	if ndims(tmpvar1) == 2 % 2 dimensions
		formula1 = [ formula1 arrayname '(:,index) = ' formula{index} ';'];
		formala2 = [ formula2 arrayname '(:,index) = ' arrayname '(:,index) - ' formula{index}  ';'];	
	else % 3 dimensions
		formula1 = [ formula1 arrayname '(:,:,index) = ' formula{index}  ';'];
		formula2 = [ formula2 arrayname '(:,:,index) = ' arrayname '(:,:,index) - ' formula{index}  ';'];	
	end;
end;

% accumulating (shufling)
% -----------------------
for index=1:naccu
	if rem(index,10) == 0,  fprintf(' %d',index); end
	if rem(index,120) == 0, fprintf('\n'); end
	
	alltrials = shuffle(alltrials);
	X = alltrials(1:cond1trials);
	eval( formula1 );
	X = alltrials(cond1trials+1:end);
	eval( formula2 );
end;

% significance level
% ------------------
for index= 1:length(formula) 
	accarray = accres{index};
    
    % size = nb_points*naccu
    % size = nb_points*naccu*times
	if ~isreal(accarray)
		accarray = sqrt(accarray .* conj(accarray)); % faster than abs()
    end;
	
    % compute bootstrap significance level
    i = round(naccu*alpha);
    switch ndims(accarray)
	 case 3, 
	     accarray = sort(accarray,3); % always sort on naccu (when 3D, naccu is the second dim)
         if strcmpi(bootside{min(length(bootside), index)}, 'upper');
			 accarray = mean(accarray(:,:,naccu-i+1:naccu),3);
	     else
			 accarray = accarray(:,:,[end:-1:1]); 
			 accarraytmp(:,:,2) = mean(accarray(:,:,1:i),3);
			 accarraytmp(:,:,1) = mean(accarray(:,:,naccu-i+1:naccu),3);
			 accarray = accarraytmp;
		 end;
	 
	 case 2, 
	     accarray = sort(accarray,2); % always sort on naccu (when 3D, naccu is the second dim)
         if strcmpi(bootside{min(length(bootside), index)}, 'upper');
			 accarray = mean(accarray(:,naccu-i+1:naccu),2);
	     else
			 accarraytmp(:,2) = mean(accarray(:,1:i),2);
			 accarraytmp(:,1) = mean(accarray(:,naccu-i+1:naccu),2);
			 accarray = accarraytmp;
		 end;
	 case 3, 
	     accarray = sort(accarray,1); % always sort on naccu (when 3D, naccu is the second dim)
         if strcmpi(bootside{min(length(bootside), index)}, 'upper');
			 accarray = mean(accarray(naccu-i+1:naccu),1);
	     else
			 accarraytmp(2) = mean(accarray(1:i),1);
			 accarraytmp(1) = mean(accarray(naccu-i+1:naccu),1);
			 accarray = accarraytmp;
		 end;
    end;
    accres{index} = accarray;
end;

if length(res1) == 1
	res1 = res1{1};
	res2 = res2{1};
	diffres = diffres{1};
	accres = accres{1};
end;
fprintf('\n');
return;

% writing a function
% ------------------
% $$$ fid = fopen('tmpfunc.m', 'w');
% $$$ fprintf(fid, 'function [accres] = tmpfunc(alltrials, cond1trials, naccu,'); 
% $$$ for index=1:length(varargin)
% $$$ 	fprintf(fid, 'arg%d', index);
% $$$ 	if index ~=length(varargin), fprintf(fid,','); end;
% $$$ end;
% $$$ fprintf(fid, ')\n');
% $$$ commandstr = [ 'for index=1:naccu, ' ]
% $$$ % 			   'if rem(index,10) == 0,  disp(index); end;' ];
% $$$ % 			   'if rem(index,10) == 0,  fprintf('' %d'',index); end;' ...
% $$$ %	           'if rem(index,120) == 0, fprintf(''\n''); end;' ];
% $$$ commandstr = [ 	commandstr 'shuffle(alltrials);' ];
% $$$ commandstr = [ 	commandstr 'X = alltrials(1:cond1trials);' ];
% $$$ if ndims(tmpvar1) == 2 % 2 dimensions
% $$$ 	commandstr = [ 	commandstr 'accres(:,index) = ' formula ';'];
% $$$ 	commandstr = [ 	commandstr 'X = alltrials(cond1trials+1:end);'];
% $$$ 	commandstr = [ 	commandstr 'accres(:,index) = accres(:,index)-' formula ';end;'];	
% $$$ else
% $$$ 	commandstr = [ 	commandstr 'res1 = ' formula ';' 10];
% $$$ 	commandstr = [ 	commandstr 'X = alltrials(cond1trials+1:end);' 10];
% $$$ 	commandstr = [ 	commandstr 'res2 = ' formula ';' 10];	
% $$$ 	commandstr = [ 	commandstr 'accres(:,:,index) = res1 - res2; end;'];	
% $$$ end;	
% $$$ fprintf(fid, commandstr);
% $$$ fclose(fid);
% $$$ profile on;
% $$$ accres = tmpfunc(alltrials, cond1trials, naccu, arg1);
% $$$ profile report;
% $$$ profile off;
% $$$ return;

% evaluating a command
% --------------------
% $$$ commandstr = [ 'for index=1:naccu, ' ...
% $$$ 			   'if rem(index,10) == 0,  fprintf('' %d'',index); end;' ...
% $$$ 	           'if rem(index,120) == 0, fprintf(''\n''); end;' ];
% $$$ commandstr = [ 	commandstr 'shuffle(alltrials);' ];
% $$$ commandstr = [ 	commandstr 'X = alltrials(1:cond1trials);' ];
% $$$ if ndims(tmpvar1) == 2 % 2 dimensions
% $$$ 	commandstr = [ 	commandstr 'accres(:,index) = ' formula ';'];
% $$$ 	commandstr = [ 	commandstr 'X = alltrials(cond1trials+1:end);'];
% $$$ 	commandstr = [ 	commandstr 'accres(:,index) = accres(:,index)-' formula ';end;'];	
% $$$ else
% $$$ 	commandstr = [ 	commandstr 'accres(:,:,index) = ' formula ';'];
% $$$ 	commandstr = [ 	commandstr 'X = alltrials(cond1trials+1:end);'];
% $$$ 	commandstr = [ 	commandstr 'accres(:,:,index) = accres(:,:,index)-' formula ';end;'];	
% $$$ end;	
% $$$ eval(commandstr);
% $$$ return;
% $$$ 
