% loadelec() - Load electrode names file for eegplot()
%
% Usage:  >> labelmatrix = loadelec('elec_file');
%
% Author: Colin Humprhies, CNL / Salk Institute, 1996
%
% See also: eegplot()

% Copyright (C) Colin Humphries, CNL / Salk Institute, Aug, 1996
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

% 01-25-02 reformated help & license, added links -ad 

function channames = loadelec(channamefile)

MAXCHANS = 256;
chans = MAXCHANS;
errorcode = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the channel names 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if channamefile ~= 0 & channamefile ~= '0'	% read file of channel names
   chid = fopen(channamefile,'r');
   if chid <3,
      fprintf('cannot open file %s.\n',channamefile);
      errorcode=2;
      channamefile = 0;
   else
      fprintf('Channamefile %s opened\n',channamefile);
   end;
   if errorcode==0,
      channames = fscanf(chid,'%s',[6 MAXCHANS]);
      channames = channames';
      [r c] = size(channames);
	 for i=1:r
	    for j=1:c
		if channames(i,j)=='.',
		   channames(i,j)=' ';
		end;
	    end;
	 end;
         % fprintf('%d channel names read from file.\n',r);
	 if (r>chans)
	    fprintf('Using first %d names.\n',chans);
	    channames = channames(1:chans,:);
	 end;
	 if (r<chans)
	    fprintf('Only %d channel names read.\n',r);
	 end;
   end;
end
if channamefile == 0 | channamefile == '0', % plot channel numbers
   channames = [];
   for c=1:chans
      if c<10,
         numeric = ['   ' int2str(c)];	% four-character fields
      else
	 numeric = ['  '  int2str(c)];
      end
      channames = [channames;numeric];
   end;
end; % setting channames

channames = str2mat(channames, ' ');	% add padding element to Y labels



