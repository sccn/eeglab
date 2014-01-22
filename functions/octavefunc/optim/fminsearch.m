% Copyright (C) 2006   Sylvain Pelissier   <sylvain.pelissier@gmail.com>
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
% along with this program; If not, see <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {[@var{x}] =} fminsearch(@var{f},@var{X0},@var{options},@var{grad},@var{P1},@var{P2}, ...)
%
% Find the minimum of a funtion of several variables.
% By default the method used is the Nelder&Mead Simplex algorithm
% @seealso{fmin,fmins,nmsmax}
% @end deftypefn

function varargout = fminsearch(funfun, X0, varargin)
    if ismatlab
        p1 = fileparts(which('fminsearch'));
        rmpath(p1);
        p2 = fileparts(which('fminsearch'));
        if ~isequal(p1, p2)
            disp( [ 'Some Octave functions should not run on Matlab' 10 'removing path to Octave fminsearch and using Matlab fminsearch' ]);
            switch nargout
                case 1, varargout{1} = fminsearch(funfun, X0, varargin{:});
                case 2, [varargout{1} varargout{2}] = fminsearch(funfun, X0, varargin{:});
                case 3, [varargout{1} varargout{2} varargout{3}] = fminsearch(funfun, X0, varargin{:});
                case 4, [varargout{1} varargout{2} varargout{3} varargout{4}]= fminsearch(funfun, X0, varargin{:});
            end;
        else
            disp( [ 'Octave functions should not run on Matlab' 10 'remove path ' p1 ]);
        end;
        return;
    end;
	if (nargin == 0); usage('[x fval] = fminsearch(funfun, X0, options, grad, varargin)'); end
    if length(varargin) > 0, options = varargin{1}; varargin(1) = []; end;
    if length(varargin) > 0, grad = varargin{1}; varargin(1) = []; end;
	if (nargin < 3); options=[]; end
	if (nargin < 4); grad=[]; end
	if (nargin < 5); varargin={}; end
	varargout{1} = fmins(funfun, X0, options, grad, varargin{:});
	varargout{2} = feval(funfun, x, varargin{:});
