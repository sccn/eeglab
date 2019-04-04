% loadelec() - Load electrode names file for eegplot()
%
% Usage:  >> labelmatrix = loadelec('elec_file');
%
% Author: Colin Humprhies, CNL / Salk Institute, 1996
%
% See also: eegplot()

% Copyright (C) Colin Humphries, CNL / Salk Institute, Aug, 1996
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

% 01-25-02 reformated help & license, added links -ad 

function channames = loadelec(channamefile)

MAXCHANS = 256;
chans = MAXCHANS;
errorcode = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the channel names 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if channamefile ~= 0 && channamefile ~= '0'	% read file of channel names
   chid = fopen(channamefile,'r');
   if chid <3,
      fprintf('cannot open file %s.\n',channamefile);
      errorcode=2;
      channamefile = 0;
   else
      fprintf('Channamefile %s opened\n',channamefile);
   end
   if errorcode==0,
      channames = fscanf(chid,'%s',[6 MAXCHANS]);
      channames = channames';
      [r c] = size(channames);
	 for i=1:r
	    for j=1:c
		if channames(i,j)=='.',
		   channames(i,j)=' ';
		end
	    end
	 end
         % fprintf('%d channel names read from file.\n',r);
	 if (r>chans)
	    fprintf('Using first %d names.\n',chans);
	    channames = channames(1:chans,:);
	 end
	 if (r<chans)
	    fprintf('Only %d channel names read.\n',r);
	 end
   end
end
if channamefile == 0 || channamefile == '0', % plot channel numbers
   channames = [];
   for c=1:chans
      if c<10,
         numeric = ['   ' int2str(c)];	% four-character fields
      else
	 numeric = ['  '  int2str(c)];
      end
      channames = [channames;numeric];
   end
end; % setting channames

channames = str2mat(channames, ' ');	% add padding element to Y labels



