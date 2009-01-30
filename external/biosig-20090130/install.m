% BIOSIG runs on Matlab and Octave. 
% This is a script installing all components in an automatically.
%  
% 1) extract the files and
% 2) save the BIOSIG files in <your_directory>
% 3) start matlab
%    cd <your_directory>
%    biosig_installer 
% 4) For a permanent installation, save the default path with 
%     PATH2RC or
%     PATHTOOL and click on the "SAVE" button. 
% 5) For removing the toolbox 
%    remove the path to 
%       HOME/tsa
%       HOME/NaN
%       HOME/BIOSIG/ and all its subdirectories
% 
%  NOTE: by default also the NaN-toolbox is installed - 
%  - a statistical toolbox for handling missing values - which 
%  changes the behaviour of some standard functions. For more  
%  information see NaN/README.TXT . In case you do not want this, 
%  you can excluded the path to NaN/*. The BIOSIG tools will still 
%  work, but does not support the handling of NaN's.

%	$Id: install.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 2003-2005,2006,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

BIOSIG_HOME = pwd;	%
if exist('t200','dir')
	% install.m may reside in .../biosig/ or above (...)
        [BIOSIG_HOME,f,e] = fileparts(BIOSIG_HOME);
elseif exist('biosig','dir')
end; 

path([BIOSIG_HOME,'/biosig'],path);			% 
path([BIOSIG_HOME,'/biosig/demo'],path);		% demos
path([BIOSIG_HOME,'/biosig/doc'],path);		% docus, Eventtable etc. 
path([BIOSIG_HOME,'/biosig/t200'],path);		% dataformat
path([BIOSIG_HOME,'/biosig/t250'],path);		% trigger and quality control
path([BIOSIG_HOME,'/biosig/t300'],path);		% signal processing and feature extraction
path([BIOSIG_HOME,'/biosig/t400'],path);		% classification
path([BIOSIG_HOME,'/biosig/t450'],path);		% statistics, false discovery rates
path([BIOSIG_HOME,'/biosig/t490'],path);		% evaluation criteria
path([BIOSIG_HOME,'/biosig/t500'],path);		% display and presentation
path([BIOSIG_HOME,'/biosig/t501'],path);		% visualization ofcoupling analysis

if ~exist('OCTAVE_VERSION','builtin'),	
	path([BIOSIG_HOME,'/biosig/viewer'],path);		% viewer
	path([BIOSIG_HOME,'/biosig/viewer/utils'],path);	% viewer
	path([BIOSIG_HOME,'/biosig/viewer/help'],path);	% viewer
end;

path([BIOSIG_HOME,'/tsa'],path);		%  Time Series Analysis
%path([BIOSIG_HOME,'/tsa/inst'],path);		%  Time Series Analysis
% some users might get confused by this
path([BIOSIG_HOME,'/NaN'],path);		%  Statistics analysis for missing data
%path([BIOSIG_HOME,'/NaN/inst'],path);		%  Statistics analysis for missing data

%%% NONFREE %%%
if exist([BIOSIG_HOME,'/biosig/NONFREE/EEProbe'],'dir'),
	path(path,[BIOSIG_HOME,'/biosig/NONFREE/EEProbe']);	% Robert Oostenveld's MEX-files to access EEProbe data
end;
if exist([BIOSIG_HOME,'/biosig/NONFREE/meg-pd-1.2-4'],'dir'),
        path(path,[BIOSIG_HOME,'/biosig/NONFREE/meg-pd-1.2-4']);	% Kimmo Uutela's library to access FIF data
end;

ver = version; 
if (str2double(ver(1:3))<7.0)
%	path(path,[BIOSIG_HOME,'/biosig/maybe-missing']);
end

% test of installation 
fun = {'isdir','ischar','strtok','str2double','strcmpi','strmatch','bitand','bitshift','sparse','strfind'};
for k = 1:length(fun),
        x = which(fun{k});
        if isempty(x) | strcmp(x,'undefined'),
                fprintf(2,'Function %s is missing\n',upper(fun{k}));     
        end;
end;
        
if exist('OCTAVE_VERSION','builtin'),	% OCTAVE
        fun = {'bitand','xmldata'};
        for k = 1:length(fun),
                try,
                        xmlstruct('<xml>v<b>v</xml>');
                catch
                        mex([BIOSIG_HOME,'/biosig/maybe-missing/xmldata.c']);
                end;
                try,
                        bitand(5,7);
                catch
                        mkoctfile([BIOSIG_HOME,'/biosig/maybe-missing/bitand.cc']);
                end;
                try,
                        x = which(fun{k});
                catch
                        x = [];
                end;	

                if isempty(x) | strcmp(x,'undefined'),
                        fprintf(2,'Function %s is missing. \n',upper(fun{k}));     
                end;
        end;
	if any(size(sparse(5,4))<0)
        	fprintf(2,'Warning: Size of Sparse does not work\n')
	end;
else
        try,
                xmlstruct('<xml>v<b>v</xml>');
        catch
                unix(['mex ',BIOSIG_HOME,'/biosig/maybe-missing/xmldata.c']);
        end;
end;

% test of installation 
% 	Some users might get confused by it. 
% nantest;	
% naninsttest; 

disp('BIOSIG-toolbox activated');
if ~exist('OCTAVE_VERSION'),	% OCTAVE
    disp('	If you want BIOSIG permanently installed, use the command PATH2RC.')
    disp('	or use PATHTOOL to select and deselect certain components.')
end; 
