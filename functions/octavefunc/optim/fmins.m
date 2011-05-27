## Copyright (C) 2003 Andy Adler
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{x}] =} fmins(@var{f},@var{X0},@var{options},@var{grad},@var{P1},@var{P2}, ...)
## 
## Find the minimum of a funtion of several variables.
## By default the method used is the Nelder&Mead Simplex algorithm
##
## Example usage:
##   fmins(inline('(x(1)-5).^2+(x(2)-8).^4'),[0;0])
## 
## @strong{Inputs}
## @table @var 
## @item f 
## A string containing the name of the function to minimize
## @item X0
## A vector of initial parameters fo the function @var{f}.
## @item options
## Vector with control parameters (not all parameters are used)
## @verbatim
## options(1) - Show progress (if 1, default is 0, no progress)
## options(2) - Relative size of simplex (default 1e-3)
## options(6) - Optimization algorithm
##    if options(6)==0 - Nelder & Mead simplex (default)
##    if options(6)==1 - Multidirectional search Method
##    if options(6)==2 - Alternating Directions search
## options(5)
##    if options(6)==0 && options(5)==0 - regular simplex
##    if options(6)==0 && options(5)==1 - right-angled simplex
##       Comment: the default is set to "right-angled simplex".
##         this works better for me on a broad range of problems,
##         although the default in nmsmax is "regular simplex"
## options(10) - Maximum number of function evaluations
## @end verbatim
## @item grad
## Unused (For compatibility with Matlab)
## @item P1,P2, ...
## Optional parameters for function @var{f} 
##
## @end table
## @end deftypefn

function ret=fmins(funfun, X0, options, grad, varargin)
    stopit = [1e-3, inf, inf, 1, 0, -1];
    minfun = 'nmsmax'; 

    if nargin < 3; options=[]; end

    if length(options)>=1; stopit(5)= options(1); end
    if length(options)>=2; stopit(1)= options(2); end
    if length(options)>=5;
        if     options(6)==0; minfun= 'nmsmax'; 
            if     options(5)==0; stopit(4)= 0;
            elseif options(5)==1; stopit(4)= 1;
            else   error('options(5): no associated simple strategy');
            end
        elseif options(6)==1; minfun= 'mdsmax';
        elseif options(6)==2; minfun= 'adsmax';
        else   error('options(6) does not correspond to known algorithm');
        end
    end
    if length(options)>=10; stopit(2)= options(10); end

    ret = feval(minfun, funfun,  X0, stopit, [], varargin{:});
endfunction
