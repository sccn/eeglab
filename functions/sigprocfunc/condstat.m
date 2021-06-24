% condstat() - accumulate surrogate data for comparing two data conditions 
%
% Usage:
%     >> [diffres, accres, res1, res2] = condstat(formula, naccu, alpha, ...
%                                              bootside, condboot, arg1, arg2 ...);
%
% Inputs:
%    formula  - [string or cell array of strings] formula(s) to compute a given measure.
%               Takes arguments 'arg1', 'arg2' ... and X as inputs. e.g.,
%        'sum(arg1(:,:,X),3) ./ sqrt(sum(arg2(:,:,X))) ./ sqrt(sum(arg3(:,:,X)))'
%    naccu    - [integer] number of accumulations of surrogate data. e.g., 200
%    alpha    - [float] significance level (0<alpha<0.5)
%    bootside - ['both'|'upper'] side of the surrogate distribution to
%               consider for significance. This parameter affect the size
%               of the last dimension of accumulation array 'accres' (size
%               is 2 for 'both' and 1 for 'upper').
%    condboot - ['abs'|'angle'|'complex'|''] When comparing two conditions, compare
%               either absolute vales ('abs'), angles ('angles') or complex values 
%               ('complex'). Either '' or 'complex' leave the formula unchanged; 
%               'abs' takes its norm before subtraction, and 'angle' normalizes 
%               each value (to norm 1) before taking the difference.
%    arg1     - [cell_array] of two 1D,2D,or 3D matrices of values to compare. 
%               The last dimension of the array is shuffled to accumulate data, the
%               other dimensions must be the same size across matrices.
%               e.g. size(arg1{1})=[100 200 500], size(arg1{2})=[100 200 395]
%    arg2     - same as arg1, note that it is compared only to itself, and has
%               nothing to do with arg1 besides using the same formula, alpha, etc.
% ...argn     - may call n number of arguement pairs    
%
% Outputs: 
%    diffres  - difference array for the actual (non-shuffled) data, if more than one
%               arg pair is called, format is a cell array of matrices.
%    accres  -  [cell array of 3D numerical arrays] for shuffled data, one per formula. 
%    res1     - result for first condition
%    res2     - result for second condition
%
% Authors: Arnaud Delorme & Scott Makeig
%          CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, La Jolla, 2002-
%
% See also: timef(), crossf()

% Copyright (C) 2002  Arnaud Delorme, Lars Kai Hansen & Scott Makeig, SCCN/INC/UCSD
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

function [diffres, accres, res1, res2] = condstat(formula, naccu, alpha, bootside, condboot, varargin);

if nargin < 6
	help condstat;
	return;
end

if ~ischar(formula) && ~iscell(formula)
	error('The first argument must be a string formula or cell array of string');
end
if ischar(formula)
	formula = { formula };
end
if ischar(bootside)
	bootside = { bootside };
end
for index = 1:length(bootside)
	if ~strcmpi(bootside{index}, 'both') && ~strcmpi(bootside{index}, 'upper')
		error('Bootside must be either ''both'' or ''upper''');
	end
end
if ischar(condboot)
	condboot = { condboot };
end
for index = 1:length(condboot)
	if isempty(condboot{index}), condboot{index} = 'complex'; end
end
		
for index = 1:length(varargin)
	if ~iscell(varargin) || length(varargin{index}) ~=2
		error('Except for the first arguments, all other arguments given to the function must be cell arrays of two numerical array');
	end
end

% accumulate coherence images (all arrays [nb_points x timesout x trials])
% ---------------------------
for index=1:length(varargin)
	tmpvar1 = varargin{index}{1};
	tmpvar2 = varargin{index}{2};
	if index == 1   % Shouldn't this be recalculated for arg2, etc.? TF 2007.06.04
		cond1trials = size(tmpvar1,ndims(tmpvar1));
		cond2trials = size(tmpvar2,ndims(tmpvar2));
		for tmpi = 1:length(formula) 
			accres{tmpi} = zeros(size(tmpvar1,1), size(tmpvar1,2), naccu);
		end
	end
	
	if ndims(tmpvar1) == 2
		eval( [ 'arg' int2str(index) '=zeros(size(tmpvar1,1), cond1trials+cond2trials);' ] );
		eval( [ 'arg' int2str(index) '(:,1:cond1trials)=tmpvar1;' ] );
		eval( [ 'arg' int2str(index) '(:,cond1trials+1:end)=tmpvar2;' ] );
	else
		eval( [ 'arg' int2str(index) '=zeros(size(tmpvar1,1), size(tmpvar1,2), cond1trials+cond2trials);' ] );
		eval( [ 'arg' int2str(index) '(:,:,1:cond1trials)=tmpvar1;' ] );
		eval( [ 'arg' int2str(index) '(:,:,cond1trials+1:end)=tmpvar2;' ] );
	end
end

fprintf('Accumulating permutation statistics:');
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
	end

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
		formula2 = [ formula2 arrayname '(:,index) = ' arrayname '(:,index) - ' formula{index}  ';'];	
	else % 3 dimensions
		formula1 = [ formula1 arrayname '(:,:,index) = ' formula{index}  ';'];
		formula2 = [ formula2 arrayname '(:,:,index) = ' arrayname '(:,:,index) - ' formula{index}  ';'];	
	end
end

% accumulating (shuffling)
% -----------------------
for index=1:naccu
	if rem(index,10) == 0,  fprintf(' %d',index); end
	if rem(index,120) == 0, fprintf('\n'); end
	
	alltrials = shuffle(alltrials);
	X = alltrials(1:cond1trials);
	eval( formula1 );
	X = alltrials(cond1trials+1:end);
	eval( formula2 );
end

% significance level
% ------------------
for index= 1:length(formula) 
	accarray = accres{index};
    
    % size = nb_points*naccu
    % size = nb_points*naccu*times
	if ~isreal(accarray)    % might want to introduce a warning here: a complex 
                            % result may not be desirable, and a single complex value
                            % in accarray could turn this from a 2-tail to 1-tail
                            % bootstrap test, leading to false positives. Hard to
                            % think of a meaningful warning though...
                            % TF 2007.06.04
		accarray = sqrt(accarray .* conj(accarray)); % faster than abs()
    end
	
    % compute bootstrap significance level
    i = round(naccu*alpha);
    switch ndims(accarray)
	 case 3
	     accarray = sort(accarray,3); % always sort on naccu (when 3D, naccu is the second dim)
         if strcmpi(bootside{min(length(bootside), index)}, 'upper')
			 accarray = accarray(:,:,floor(naccu-i/2+1));
	     else
			 accarray = accarray(:,:,[end:-1:1]); 
			 accarraytmp(:,:,2) = accarray(:,:,ceil(i/2));
			 accarraytmp(:,:,1) = accarray(:,:,floor(naccu-i/2+1));
			 accarray = accarraytmp;
		 end
	 
	 case 2
	     accarray = sort(accarray,2); % always sort on naccu (when 3D, naccu is the second dim)
         if strcmpi(bootside{min(length(bootside), index)}, 'upper')
			 accarray = accarray(:,floor(naccu-i/2+1));
	     else
			 accarraytmp(:,2) = accarray(:,ceil(i/2));
			 accarraytmp(:,1) = accarray(:,floor(naccu-i/2+1));
			 accarray = accarraytmp;
		 end
	 case 1
	     accarray = sort(accarray,1); % always sort on naccu (when 3D, naccu is the second dim)
         if strcmpi(bootside{min(length(bootside), index)}, 'upper')
			 accarray = accarray(floor(naccu-i/2+1));
	     else
			 accarraytmp(2) = accarray(ceil(i/2));
			 accarraytmp(1) = accarray(floor(naccu-i/2+1));
			 accarray = accarraytmp;
		 end
    end
    accres{index} = accarray;
end

if length(res1) == 1
	res1 = res1{1};
	res2 = res2{1};
	diffres = diffres{1};
	accres = accres{1};
end
fprintf('\n');
return;

% writing a function
% ------------------
% $$$ fid = fopen('tmpfunc.m', 'w');
% $$$ fprintf(fid, 'function [accres] = tmpfunc(alltrials, cond1trials, naccu,'); 
% $$$ for index=1:length(varargin)
% $$$ 	fprintf(fid, 'arg%d', index);
% $$$ 	if index ~=length(varargin), fprintf(fid,','); end
% $$$ end
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
