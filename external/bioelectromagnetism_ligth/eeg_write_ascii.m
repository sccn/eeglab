function eeg_write_ascii(file,data,format)

% eeg_write_ascii - Write matrix of EEG data into ascii data format
%
% Useage: eeg_write_ascii(file,data,format)
%
% where:    file is the path + filename, eg "c:\data.dat"
%           data is a data matrix (floating point)
%           format is an fprintf format string ('%12.6f' default)
%
% comment:  Uses the fprintf method, which is quicker than the 
%           dlmwrite() function. It remains to be seen whether it 
%           is quicker or provides more format options than the 
%           built-in SAVE function.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2001, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('format','var') format = '%12.6f'; end

    % transpose data for output, see note below
    data = data';
    
    [r,c] = size(data);         % size of matrix: r rows and c columns
    
    % create a format string long enough to provide a format for every 
    % row of the data - put a line return at the end of this string
    formats = '';
    for i = 1:r
        formats = strcat(formats, ' ', format);
    end
    formats = strcat(formats,'\n');
    
    % Output data
    fid = fopen(file,'w');  fprintf(fid,formats,data);  fclose(fid);
    
    fprintf('Completed ascii write to:\t%s\n', file);
    return

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NOTE on TRANSPOSE
%   this routine transposes in the input data matrix to
%   accommodate the matlab fprintf routine, which outputs in
%   column major format to effectively transpose the input.  For
%   example, 'help fprintf' provides this example:
%
%   x = 0:.1:1; y = [x; exp(x)];
%   fprintf('%6.2f  %12.8f\n',y);
%           
%   where y is a 2row, 11 column matrix and the output is:
%
%   0.00    1.00000000
%   0.10    1.10517092
%        ...
%   1.00    2.71828183
    
