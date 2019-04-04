% mmo() - create a memory-mapped data class
%
% Usage:
%   >> data_class = mmo(data);
%
% Inputs:
%   data         - input data or data file
%
% Outputs:
%   data_class    - output dataset class
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Nov. 2012-

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD
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

classdef mmo
    properties
        dataFile   = [];
        dimensions = [];
        writable   = true;  % duplicate of the writable field of memmapfile
        workspace  = [];
        type       = 'mmo';
        transposed = false;
        debug      = false;
    end
        
    methods
        function dataout = mmo(dataFileIn, datadims, writableVal, transposedVal, debugVal)
            if nargin < 3
                writableVal = true;
            end
            if nargin < 4
                transposedVal = false;
            end
            if nargin < 5
                debugVal = false;
            end
            
            % check that the file is not empty
            % --------------------------------
            if ~isempty(dataFileIn)
                if ~ischar(dataFileIn)
                    error('First input must be a file name');
                end

                dirContent = dir(dataFileIn);
                if isempty(dirContent)
                    error([ 'Data file ''' dataFileIn '''not found' ]);
                elseif dirContent(1).bytes == 0
                    error([ 'Empty data file ''' dataFileIn '''' ]);
                end
            else
                dataFileIn = mmo.getnewfilename;
                fid = fopen(dataFileIn, 'w');
                if fid == -1, error('Cannot open new file'); end
                if length(datadims) == 1, datadims(2) = 1; end
                tmpdata = zeros(1, prod(datadims(2:end)), 'single');
                for index = 1:datadims(1)
                    fwrite(fid, tmpdata, 'float');
                end
                fclose(fid);
            end
                
            % test memory map but discards it
            % -------------------------------
            test = memmapfile(dataFileIn, 'writable', writableVal, 'format', { 'single' datadims 'x' });
            clear test;
            
            % set fields
            % ----------
            while datadims(end) == 1 && length(datadims) > 1
                datadims(end) = [];
            end
            if transposedVal
                if length(datadims) == 1, datadims = [1 datadims];
                else                      datadims = [datadims(2:end) datadims(1)];
                end
            end
            
            dataout.dataFile   = dataFileIn;
            dataout.dimensions = datadims;
            dataout.writable   = writableVal;
            dataout.transposed = transposedVal;
            
            % set workspace
            % -------------
            dataout = updateWorkspace(dataout);
%             stack = dbstack;
%             stack(1) = [];
%             stack = rmfield(stack, 'line');
%             dataout.workspace = stack;
            dataout.debug     = debugVal;
        end
        
        % this function updates the function workspace
        % --------------------------------------------
        function obj = updateWorkspace(obj);
            stack = dbstack;
            stack(1:2) = [];
            stack      = rmfield(stack, 'line');
            obj.workspace = stack;
        end
        
        % function to check copies (only used for testing,
        % implemented in subasign as well)
        % --------------------------------
        function ncopies = checkcopies(obj);
            ncopies = checkworkspace(obj);
            if ncopies < 2
                s = evalin('caller', 'whos');
                for index = 1:length(s)
                    if strcmpi(s(index).class, 'struct') || strcmpi(s(index).class, 'cell')
                        tmpVar = evalin('caller', s(index).name);
                        ncopies = ncopies + checkcopies_local(obj, tmpVar);
                    elseif strcmpi(s(index).class, 'mmo')
                        if s(index).persistent || s(index).global
                            disp('Warning: mmo objects should not be made persistent or global. Behavior is unpredictable.');
                        end
                        tmpVar = evalin('caller', s(index).name);
                        if isequal(tmpVar, obj), ncopies = ncopies + 1; end
                        if ncopies > 1, break; end
                    end
                end
            end
        end
        
        % numerical implementations of basic functions
        % --------------------------------------------
        function obj2 = log(obj1); obj2 = unitaryopp(@log, obj1); end
        function val  = mean(obj,dim); if nargin <1, dim=1; end; val = sum(obj,dim)/size(obj,dim); end
        function val  = std(varargin); val = sqrt(var(varargin{:})); end
        function obj3 = minus(obj1, obj2); obj3 = binaryopp(@minus, obj1, obj2); end
        function obj3 = plus(obj1, obj2); obj3 = binaryopp(@plus, obj1, obj2); end
        function obj3 = time(obj1, obj2); obj3 = binaryopp(@time, obj1, obj2); end
    end
    
    methods (Static)
        function str  = getnewfilename; str = fullfile(gettempfolder(true), sprintf('memapdata_%.9d%.9d.fdt', round(rand(1)*10^9), round(rand(1)*10^9))); end
    end
end
