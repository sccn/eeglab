function [output,status] = plugin_urlread(urlChar,method,params)
%URLREAD Returns the contents of a URL as a string.
%   S = URLREAD('URL') reads the content at a URL into a string, S.  If the
%   server returns binary data, the string will contain garbage.
%
%   S = URLREAD('URL','method',PARAMS) passes information to the server as
%   part of the request.  The 'method' can be 'get', or 'post' and PARAMS is a 
%   cell array of param/value pairs.
%
%   [S,STATUS] = URLREAD(...) catches any errors and returns 1 if the file
%   downloaded successfully and 0 otherwise.
%
%   Examples:
%   s = urlread('http://www.mathworks.com')
%   s = urlread('ftp://ftp.mathworks.com/README')
%   s = urlread(['file:///' fullfile(prefdir,'history.m')])
% 
%   From behind a firewall, use the Preferences to set your proxy server.
%
%   See also URLWRITE.

%   Adapted by A. Delorme from 
%   Matthew J. Simoneau, 13-Nov-2001
%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.3.2.12 $ $Date: 2011/09/03 22:43:00 $

% This function requires Java.
if ~usejava('jvm')
   error(message('MATLAB:urlread:NoJvm'));
end

if exist('OCTAVE_VERSION', 'builtin') == 0
    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;

    % Be sure the proxy settings are set.
    com.mathworks.mlwidgets.html.HTMLPrefs.setProxySettings
end



% Check number of inputs and outputs.
if ~ischar(urlChar)
    error('MATLAB:urlread:InvalidInput','The first input, the URL, must be a character array.');
end
if (nargin > 1) && ~strcmpi(method,'get') && ~strcmpi(method,'post')
    error('MATLAB:urlread:InvalidInput','Second argument must be either "get" or "post".');
end

% Do we want to throw errors or catch them?
if nargout == 2
    catchErrors = true;
else
    catchErrors = false;
end

% Set default outputs.
output = '';
status = 0;

% GET method.  Tack param/value to end of URL.
if (nargin > 1) && strcmpi(method,'get')
    if mod(length(params),2) == 1
        error('MATLAB:urlread:InvalidInput','Invalid parameter/value pair arguments.');
    end
    for i=1:2:length(params)
        if (i == 1), separator = '?'; else separator = '&'; end
        param = char(java.net.URLEncoder.encode(params{i}));
        value = char(java.net.URLEncoder.encode(params{i+1}));
        urlChar = [urlChar separator param '=' value];
    end
end

% Create a urlConnection.
[urlConnection,errorid,errormsg] = plugin_urlreadwrite(mfilename,urlChar);
if isempty(urlConnection)
    if catchErrors, return
    else error(errorid,errormsg);
    end
end
urlConnection.setReadTimeout(5000); % timeout in 5 seconds

% POST method.  Write param/values to server.
if (nargin > 1) && strcmpi(method,'post')
    try
        urlConnection.setDoOutput(true);
        urlConnection.setRequestProperty( ...
            'Content-Type','application/x-www-form-urlencoded');
        printStream = java.io.PrintStream(urlConnection.getOutputStream);
        for i=1:2:length(params)
            if (i > 1), printStream.print('&'); end
            param = char(java.net.URLEncoder.encode(params{i}));
            value = char(java.net.URLEncoder.encode(params{i+1}));
            printStream.print([param '=' value]);
        end
        printStream.close;
    catch
        if catchErrors, return
        else error('MATLAB:urlread:ConnectionFailed','Could not POST to URL.');
        end
    end
end

% Read the data from the connection.
try
    inputStream = urlConnection.getInputStream;
    byteArrayOutputStream = java.io.ByteArrayOutputStream;
    % This StreamCopier is unsupported and may change at any time.
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
    isc.copyStream(inputStream,byteArrayOutputStream);
    inputStream.close;
    byteArrayOutputStream.close;
    output = native2unicode(typecast(byteArrayOutputStream.toByteArray','uint8'),'UTF-8');
catch
    if catchErrors, return
    else error('MATLAB:urlread:ConnectionFailed','Error downloading URL. Your network connection may be down or your proxy settings improperly configured.');
    end
end

status = 1;
