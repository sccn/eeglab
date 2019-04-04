% EEGLAB cross-platform compiling script
% should be run on a newly checked out EEGLAB version as 
% some folder are temporarily modified
%
% Arnaud Delorme - August 3rd, 2009

% Copyright (C) 2009 Arnaud Delorme
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

disp('This function will compile EEGLAB in the output folder');
disp('provided below. You may also enter a path relative to the EEGLAB');
disp('folder: ../compiled_EEGLAB for instance');
outputfolder = input('Enter output folder name:','s');

eeglab; close;
path_eeglab = fileparts(which('eeglab'));
cd(path_eeglab);

% deal with VisEd plugin (VisEd is both the name of the folder and the
% function inside and this creates a problem
path_vised = fileparts(which('VisEd'))
try, movefile( path_vised, [ path_vised '2' ]); catch, end
addpath([ path_vised '2' ]);

path_fileio = fileparts(which('chantype'));
try, movefile( fullfile(path_fileio, '@uint64'), fullfile(path_fileio, 'uint64') ); catch, end
try, movefile( fullfile(path_fileio, 'private', 'buffer.m')      ,  fullfile(path_fileio, 'private', 'bufferold.m') ); catch, end
try, movefile( fullfile(path_fileio, 'private', 'read_24bit.m')  ,  fullfile(path_fileio, 'private', 'read_24bitold.m')); catch, end
try, movefile( fullfile(path_fileio, 'private', 'read_ctf_shm.m'),  fullfile(path_fileio, 'private', 'read_ctf_shmold.m')); catch, end
try, movefile( fullfile(path_fileio, 'private', 'write_ctf_shm.m'), fullfile(path_fileio, 'private', 'write_ctf_shmold.m')); catch, end
path_fileio = path_fileio(length(path_eeglab)+2:end);
files_fileio         = fullfile(path_fileio, '*.m');
files_fileio_private = fullfile(path_fileio, 'private', '*.m');

path_fieldtrip = fileparts(which('electroderealign'));
addpath(fullfile(path_fieldtrip, 'public'));
try, movefile( fullfile(path_fieldtrip, 'fileio', '@uint64'), fullfile(path_fieldtrip, 'fileio', 'uint64') ); catch, end
try, movefile( fullfile(path_fieldtrip, '@uint64'), fullfile(path_fieldtrip, 'uint64') ); catch, end
try, movefile( fullfile(path_fieldtrip, 'topoplot.m'), fullfile(path_fieldtrip, 'topoplotold.m') ); catch, end
path_fieldtrip = path_fieldtrip(length(path_eeglab)+2:end);
files_fieldtrip         = fullfile(path_fieldtrip, '*.m');
files_public            = fullfile(path_fieldtrip, 'public', '*.m');
files_public_private    = fullfile(path_fieldtrip, 'public', 'private', '*.m');
files_fieldtrip_private = fullfile(path_fieldtrip, 'private', '*.m');
files_forwinv           = fullfile(path_fieldtrip, 'forward', '*.m');
files_forwinv_private   = fullfile(path_fieldtrip, 'forward', 'private', '*.m');
files_inverse           = fullfile(path_fieldtrip, 'inverse', '*.m');
files_inverse_private   = fullfile(path_fieldtrip, 'inverse', 'private', '*.m');

try
    rmpath('C:\Documents and Settings\delorme\My Documents\eeglab\plugins\editevents_arno');
catch, end
path_biosig = fileparts(which('install'));
path_biosig = path_biosig(length(path_eeglab)+2:end);
biosig  = ' sopen.m sclose.m sread.m ';
[allfiles3 fieldt]  = scanfold('external/fieldtrip-partial');
% note that the order is important if the two first folders are inverted,
% it does not work

% fieldt  = [ ' -a ' files_fieldtrip_private ...
%             ' -a ' files_fieldtrip ...
%             ' -a ' files_public ...
%             ' -a ' files_forwinv ...
%             ' -a ' files_forwinv_private ...
%             ' -a ' files_inverse ...
%             ' -a ' files_inverse_private ...
%             ' -a ' files_fileio ...
%             ' -a ' files_fileio_private ];
%fieldt  = [ ' -a external\fieldtrip-20090727\private\*.m -a external\fieldtrip-20090727\*.m  ' ...
%            ' -a external\fieldtrip-20090727\forwinv\*.m -a external\fieldtrip-20090727\forwinv\private\*.m ' ...
%            ' -a external\fileio-20090511\*.m -a external\fileio-20090511\private\*.m ' ];
% topoplot
% uint64 in fieldtrip and file-io
% other mex files in file-io private folder
[allfiles1 plugins]   = scanfold('plugins/');
[allfiles2 functions] = scanfold('functions/');

eval([ 'mcc -v -m eeglab' biosig plugins functions fieldt ]);
%eval([ 'mcc -v -C -m eeglab' biosig plugins functions fieldt ]);

mkdir(fullfile(outputfolder));
comp = computer;
if strcmpi(comp(1:2), 'PC')
    copyfile( 'eeglab.exe', fullfile(outputfolder, 'eeglab.exe'), 'f');
    copyfile( 'eeglab.ctf', fullfile(outputfolder, 'eeglab.ctf'), 'f');
else
    copyfile( 'eeglab', fullfile(outputfolder, 'eeglab'), 'f');
    copyfile( 'eeglab', fullfile(outputfolder, 'eeglab'), 'f');
    copyfile( 'eeglab.ctf', fullfile(outputfolder, 'eeglab.ctf'), 'f');
end

% copy BESA files etc
% -------------------
dipfitdefs;
mkdir(fullfile(outputfolder, 'help'));
tmpf = which('eeglablicense.txt');      copyfile(tmpf, fullfile(outputfolder, 'help', 'eeglablicense.txt'));
tmpf = which('eeg_optionsbackup.m');    copyfile(tmpf, fullfile(outputfolder, 'eeg_optionsbackup.txt'));
tmpf = which('eeg_options.m');          copyfile(tmpf, fullfile(outputfolder, 'eeg_options.txt'));
tmpf = which('mheadnew.xyz');           copyfile(tmpf, fullfile(outputfolder, 'mheadnew.xyz'));
tmpf = which('mheadnew.mat');           copyfile(tmpf, fullfile(outputfolder, 'mheadnew.mat'));
tmpf = which('mheadnew.elp');           copyfile(tmpf, fullfile(outputfolder, 'mheadnew.elp'));
tmpf = which('mheadnew.transform');     copyfile(tmpf, fullfile(outputfolder, 'mheadnew.transform'));
mkdir(fullfile(outputfolder, 'standard_BEM'));
mkdir(fullfile(outputfolder, 'standard_BEM', 'elec'));
copyfile(template_models(2).hdmfile , fullfile(outputfolder, 'standard_BEM', 'standard_vol.mat'));
copyfile(template_models(2).mrifile , fullfile(outputfolder, 'standard_BEM', 'standard_mri.mat'));
copyfile(template_models(2).chanfile, fullfile(outputfolder, 'standard_BEM',  'elec', 'standard_1005.elc'));
mkdir(fullfile(outputfolder, 'standard_BESA'));
copyfile(template_models(1).hdmfile , fullfile(outputfolder, 'standard_BESA', 'standard_BESA.mat'));
copyfile(template_models(1).mrifile , fullfile(outputfolder, 'standard_BESA', 'avg152t1.mat'));
copyfile(template_models(1).chanfile, fullfile(outputfolder, 'standard_BESA', 'standard-10-5-cap385.elp'));
copyfile(fullfile(path_biosig, 'doc', 'units.csv'),              fullfile(outputfolder, 'units.csv'));
copyfile(fullfile(path_biosig, 'doc', 'leadidtable_scpecg.txt'), fullfile(outputfolder, 'leadidtable_scpecg.txt'));
copyfile(fullfile(path_biosig, 'doc', 'elecpos.txt'),            fullfile(outputfolder, 'elecpos.txt'));
copyfile(fullfile(path_biosig, 'doc', 'DecimalFactors.txt'),     fullfile(outputfolder, 'DecimalFactors.txt'));
copyfile('sample_locs', fullfile(outputfolder, 'sample_locs'));
copyfile('sample_data', fullfile(outputfolder, 'sample_data'));

% copy all files for help
% -----------------------
disp('Copying help files');
allfiles = { allfiles1{:} allfiles2{:} };
for index = 1:length(allfiles)
    tmpp = which(allfiles{index});
    copyfile(tmpp, fullfile(outputfolder, 'help', allfiles{index}));
end

% copy MCR file and visual C++ librairies
% ---------------------------------------
if strcmpi(comp(1:2), 'PC')
    copyfile(fullfile(matlabroot, 'toolbox', 'compiler', 'deploy', 'win32', 'MCRInstaller.exe'), fullfile(outputfolder, 'MCRInstaller.exe'));
    copyfile(fullfile(matlabroot, 'bin', 'win32', 'vcredist_x86.exe'), fullfile(outputfolder, 'vcredist_x86.exe'));

    fid = fopen(fullfile(outputfolder, 'setup.bat'), 'w');
    if fid == -1, disp('Error: cannot create setup file');
    else
        fprintf(fid, 'echo off\r\n');
        fprintf(fid, 'echo Deploying EEGLAB project.\r\n');
        fprintf(fid, 'echo Running MCRInstaller\r\n');
        fprintf(fid, 'echo You may delete the MCRInstaller.exe and vcredist_x86.exe file after setup is complete\r\n');
        fprintf(fid, 'MCRInstaller.exe\r\n');
        fprintf(fid, 'echo Running Visual C++ deployment functions\r\n');
        fprintf(fid, 'vcredist_x86.exe\r\n');
        fprintf(fid, 'echo Installation complete.\r\n');
        fprintf(fid, 'echo Please refer to the documentation for any additional setup steps.\r\n');
        fprintf(fid, 'echo Now starting EEGLAB...\r\n');
        fprintf(fid, 'echo To start EEGLAB in the future, simply click on the EEGLAB.EXE file\r\n');
        fprintf(fid, 'eeglab.exe\r\n');
        fclose(fid);
    end
    fid = fopen(fullfile(outputfolder, 'eeglab.exe.manifest'), 'w');
    if fid == -1, disp('Error: cannot create manifest file');
    else
        fprintf(fid, '<?xml version=''1.0'' encoding=''UTF-8'' standalone=''yes''?>\r\n'); 
        fprintf(fid, '<assembly xmlns=''urn:schemas-microsoft-com:asm.v1'' manifestVersion=''1.0''>\r\n');
        fprintf(fid, '  <dependency>\r\n');
        fprintf(fid, '    <dependentAssembly>\r\n');
        fprintf(fid, '      <assemblyIdentity type=''win32'' name=''Microsoft.VC80.CRT'' version=''8.0.50727.762'' processorArchitecture=''x86'' publicKeyToken=''1fc8b3b9a1e18e3b'' />\r\n'); 
        fprintf(fid, '    </dependentAssembly>\r\n');
        fprintf(fid, '  </dependency>\r\n');
        fprintf(fid, '</assembly>\r\n');
        fclose(fid);
    end
end

% cleaning up
% -----------
try, movefile( [ path_vised '2' ], path_vised); catch, end
try, movefile( fullfile(path_fileio, 'uint64'), fullfile(path_fileio, '@uint64') ); catch, end
try, movefile( fullfile(path_fileio, 'private', 'bufferold.m')      ,  fullfile(path_fileio, 'private', 'buffer.m') ); catch, end
try, movefile( fullfile(path_fileio, 'private', 'read_24bitold.m')  ,  fullfile(path_fileio, 'private', 'read_24bit.m')); catch, end
try, movefile( fullfile(path_fileio, 'private', 'read_ctf_shmold.m'),  fullfile(path_fileio, 'private', 'read_ctf_shm.m')); catch, end
try, movefile( fullfile(path_fileio, 'private', 'write_ctf_shmold.m'), fullfile(path_fileio, 'private', 'write_ctf_shm.m')); catch, end
try, movefile( fullfile(path_fieldtrip, 'fileio', 'uint64'), fullfile(path_fieldtrip, 'fileio', '@uint64') ); catch, end
try, movefile( fullfile(path_fieldtrip, 'uint64'), fullfile(path_fieldtrip, '@uint64') ); catch, end
try, movefile( fullfile(path_fieldtrip, 'topoplotold.m'), fullfile(path_fieldtrip, 'topoplot.m') ); catch, end

return

%histforexe(allfiles1, 'help');

% help for lisence
% 

