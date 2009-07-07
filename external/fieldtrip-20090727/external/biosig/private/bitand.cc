/* Copyright (C) 2000  Kai Habel
**
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 
*/

/*
INSTALLATION
- copy this file and the Makefile to directory of octave's LOADPATH
- compile this file:
	make
*/

#include <sys/types.h>
#include <climits>
#include <octave/oct.h>
#include <octave/lo-ieee.h>
#ifdef USE_OCTAVE_NAN
#define lo_ieee_nan_value() octave_NaN
#endif

//using namespace std;

typedef unsigned long bitop_int;
const unsigned int ULONG_SIZE=CHAR_BIT*sizeof(bitop_int);
const unsigned int BIT_AND = 1;
const unsigned int BIT_OR = 2;
const unsigned int BIT_XOR = 3;

inline unsigned int 
max(unsigned int x, unsigned int y) { 
	return x > y ? x : y; 
} 

double
scalar_bitop(double x,double y,unsigned int op) {
	double a=lo_ieee_nan_value();
	if (op ==BIT_AND) 
		if ((x>=0)&&(y>=0)&&((x<=ULONG_MAX)||(y<=ULONG_MAX))) {
			bitop_int xval=static_cast<bitop_int>( floor(x) ) & ULONG_MAX;
			bitop_int yval=static_cast<bitop_int>( floor(y) ) & ULONG_MAX;
			a = static_cast<double>(xval & yval);
		}	
	else if (op>=BIT_OR)
		if ((x>=0)&&(x<=ULONG_MAX)&&(y>=0)&&(y<=ULONG_MAX)) {
			bitop_int xval=static_cast<bitop_int>( floor(x) );
			bitop_int yval=static_cast<bitop_int>( floor(y) );
			if (op==BIT_OR)
				a = static_cast<double>(xval | yval);
			else if (op==BIT_XOR)
				a = static_cast<double>(xval ^ yval);
		}
	return(a);
}

octave_value_list
bitop(Matrix xmat,Matrix ymat,unsigned int op) {

	octave_value_list retval;

	bool is_scalar_op=false,is_matrix_op=false;
	unsigned int xr=xmat.rows();
	unsigned int yr=ymat.rows();
	unsigned int xc=xmat.columns();
	unsigned int yc=ymat.columns();

	if ( (xr*xc)==1 || (yr*yc)==1) 
		is_scalar_op=true;
	if ( (xr==yr)&&(xc==yc) )
		is_matrix_op=true;
	if (is_matrix_op || is_scalar_op) {
		unsigned int i,j,k,l;
		unsigned int r=max(xr,yr),c=max(xc,yc);
		Matrix a(r,c);

		for(i=0;i<xr;i++) {
			for(j=0;j<xc;j++) {
				if (is_scalar_op) {
					for(k=0;k<yr;k++) {
						for(l=0;l<yc;l++) {
							a(i+k,j+l)=scalar_bitop( xmat(i,j),ymat(k,l),op );
						}
					}
				}
				else {
					// is_matrix_op
					a(i,j)=scalar_bitop( xmat(i,j),ymat(i,j),op );
				}
			}
		}
		retval(0)=a;
	}
	else 
		error("size of x and y must match, or one operand must be a scalar");
	return(retval);
}

/*
%!assert(bitand(7,14),6);
*/
DEFUN_DLD (bitand, args, ,
	"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{A} =} bitand (@var{x}, @var{y})\n\
calculates the bitwise AND of nonnegative integers.\n\
@var{x},@var{y} must be in range [0..bitmax]\n\
@seealso{bitor,bitxor,bitset,bitget,bitcmp,bitshift,bitmax}\n\
@end deftypefn")
{
	octave_value_list retval;
	
	int nargin = args.length();
	if (!(nargin==2)) {
		print_usage ("bitand");
		return retval;
	}

	if (args(0).is_real_type()&&args(1).is_real_type()) {
		Matrix x = args(0).matrix_value();
		Matrix y = args(1).matrix_value();
		retval=bitop(x,y,BIT_AND);
	}
	else 
		error("both operands must be of real data type");
	return retval;
}

/*
%!assert(bitor(7,14),15);
*/
DEFUN_DLD (bitor, args, ,
	"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{A} =} bitor (@var{x}, @var{y})\n\
calculates the bitwise OR of nonnegative integers.\n\
@var{x},@var{y} must be in range [0..bitmax]\n\
@seealso{bitor,bitxor,bitset,bitget,bitcmp,bitshift,bitmax}\n\
@end deftypefn")
{
	octave_value_list retval;
	
	int nargin = args.length();
	if (!(nargin==2)) {
		print_usage ("bitor");
		return retval;
	}

	if (args(0).is_real_type()&&args(1).is_real_type()) {
		Matrix x = args(0).matrix_value();
		Matrix y = args(1).matrix_value();
		retval=bitop(x,y,BIT_OR);
	}
	else 
		error("both operands must be of real data type");
	return retval;
}

/*
%!assert(bitxor(7,14),9);
*/
DEFUN_DLD (bitxor, args, ,
	"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{A} =} bitxor (@var{x}, @var{y})\n\
calculates the bitwise XOR of nonnegative integers.\n\
@var{x},@var{y} must be in range [0..bitmax]\n\
@seealso{bitand,bitor,bitset,bitget,bitcmp,bitshift,bitmax}\n\
@end deftypefn")
{
	octave_value_list retval;
	
	int nargin = args.length();
	if (!(nargin==2)) {
		print_usage ("bitxor");
		return retval;
	}

	if (args(0).is_real_type()&&args(1).is_real_type()) {
		Matrix x = args(0).matrix_value();
		Matrix y = args(1).matrix_value();
		retval=bitop(x,y,BIT_XOR);
	}
	else 
		error("both operands must be of real data type");
	return retval;
}


/*
%!assert(bitmax != 0);
*/
DEFUN_DLD (bitmax, args, ,
	"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{A} =} bitmax\n\
returns the the maximum unsigned integer.\n\
@seealso{bitand,bitor,bitxor,bitset,bitget,bitcmp,bitshift}\n\
@end deftypefn")
{
	octave_value_list retval;
	if (args.length()!=0) 
		print_usage ("bitmax");
	else 
		retval(0)=octave_value(static_cast<double>(ULONG_MAX));
	return retval;
}
