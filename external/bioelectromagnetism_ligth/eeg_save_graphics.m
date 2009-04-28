function eeg_save_graphics(F,type,p,confirm)

% eeg_save_graphics - Save figure as graphics file
% 
% Usage: eeg_save_graphics(F,type,p,confirm)
% 
% 'F' is a figure handle, when empty uses 'gcf'
% 
% 'type' is the format of the graphic file, eg:
% 
%       'png'  - portable network graphics (default)
%       'jpg'  - JPEG image, quality level of 90
%       'tiff' - TIFF with packbits (lossless run-length encoding)
%       'eps'  - encapsulated postscript
% 
% Note that tiff and eps files can be very large.
% The default resolution is 300 DPI.
% 
% 'p' is the data struct of eeg_toolbox
% 
% 'confirm' is a boolean to control direct or GUI
% based saving.  0 = direct, 1 = GUI confirmation.
% 
% The output file is saved to p.volt.path with
% the filename of the voltage file and the timing
% of the graphic appended (assumes that the graphic
% is a topographic map).
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:53 $

% Licence:  GNU GPL, no implied or express warranties
% History:  04/2002, Darren.Weber_at_radiology.ucsf.edu
%           - extracted out of eeg_contours_engine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('F','var'),
    if isempty(F), F = gcf; end
else
    F = gcf;
end
if exist('type','var'),
    if isempty(type), type = 'png'; end
else
    type = 'png';
end
if exist('confirm','var'),
    if isempty(confirm), confirm = 0; end
else
    confirm = 0;
end


[path,file,ext] = fileparts(strcat(p.volt.path,filesep,p.volt.file));
if p.topoView,
    file = strcat(file, sprintf('_%s_%08.2f',p.topoView,p.volt.sampleTime),'.',type);
else
    file = strcat(file, sprintf('_%08.2f',p.volt.sampleTime),'.',type);
end
file = fullfile(path,file);

if confirm,
    [filename, filepath] = uiputfile(file, 'EEG Save Graphic');
    if ~isequal(filename,0),
        file = fullfile(filepath,filename);
        fprintf('EEG_SAVE_GRAPHIC: Saving to:...\n...%s\n',file);
        saveas(F,file);
    end
    
else
    type = lower(type);
    if strmatch(type,'png'),
        driver = '-dpng';
        option1 = '-noui';
        option2 = '-r300';
        %option3 = '-zbuffer';
        fprintf('EEG_SAVE_GRAPHIC: Saving to:...\n...%s\n',file);
        print(F,driver,option1,option2,file);
    elseif strmatch(type,'jpeg'),
        driver  = '-djpeg90'; % JPEG image, quality level of nn (90 here)
        option1 = '-r300';
        fprintf('EEG_SAVE_GRAPHIC: Saving to:...\n...%s\n',file);
        print(F,driver,option1,file);
    elseif strmatch(type,'eps'),
        driver = '-depsc2'; % Encapsulated Level 2 Color PostScript
        option1 = '-r300';
        option2 = '-tiff';
        fprintf('EEG_SAVE_GRAPHIC: Saving to:...\n...%s\n',file);
        print(F,driver,option1,option2,file);
    elseif strmatch(type,'tiff'),
        driver = '-dtiff'; % TIFF with packbits (lossless run-length encoding) compression
        option1 = '-r300';
        fprintf('EEG_SAVE_GRAPHIC: Saving to:...\n...%s\n',file);
        print(F,driver,option1,file);
    else
        fprintf('EEG_SAVE_GRAPHIC: Saving to:...\n...%s\n',file);
        saveas(F,file);
    end
end

return
