function [urlConnection,errorid,errormsg] = plugin_urlreadwrite(fcn,urlChar)
%URLREADWRITE A helper function for URLREAD and URLWRITE.

%   Adapted by A. Delorme from 
%   Matthew J. Simoneau, June 2005
%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $ $Date: 2011/10/22 22:05:21 $

% Default output arguments.
urlConnection = [];
errorid = '';
errormsg = '';

% Determine the protocol (before the ":").
protocol = urlChar(1:min(find(urlChar==':'))-1);

% Try to use the native handler, not the ice.* classes.
switch protocol
    case 'http'
        try
            handler = sun.net.www.protocol.http.Handler;
        catch exception %#ok
            handler = [];
        end
    case 'https'
        try
            handler = sun.net.www.protocol.https.Handler;
        catch exception %#ok
            handler = [];
        end
    otherwise
        handler = [];
end

%Try to fix proxy
useNewProxyInfo = 0;
try
    if exist('OCTAVE_VERSION', 'builtin') == 0
        [matlabversion matlabdate] = version;
        matlabversion2 = regexp(matlabversion,'R20(\d\d)([abcd])','match');
        if ~isempty(matlabversion2)
            matlabyear = regexp(matlabversion2,'20(\d\d)','match');
            if ~isempty(matlabyear)
                matlabyear = str2double(matlabyear{1,1});
                %matlabsubversion = regexp(matlabversion2,strcat(num2str(matlabyear), "([abcd])"),'match','once');
                useNewProxyInfo = matlabyear>=2018;
            end
        end
    end
catch
    warning('An error occurred when checking MATLAB/Octave version');
end

try
    if useNewProxyInfo == 1
        s = settings;
        proxyHost = s.matlab.web.ProxyHost.ActiveValue;
        proxyPort = s.matlab.web.ProxyPort.ActiveValue;
        switch protocol
            case 'http'
                java.lang.System.setProperty('http.proxyHost',proxyHost);
                java.lang.System.setProperty('http.proxyPort',proxyPort);
            case 'https'
                s = settings;
                proxyHost = s.matlab.web.ProxyHost.ActiveValue;
                proxyPort = s.matlab.web.ProxyPort.ActiveValue;
                java.lang.System.setProperty('https.proxyHost',proxyHost);
                java.lang.System.setProperty('https.proxyPort',proxyPort);
            otherwise
                warning('Unknown web protocol');
        end
    end
catch
    warning('An error occurred when checking proxy information using new format');
end

% Create the URL object.
try
    if isempty(handler)
        url = java.net.URL(urlChar);
        % url = javaObject('java.net.URL', urlChar); % Octave
    else
        url = java.net.URL([],urlChar,handler);
        % url = javaObject('java.net.URL', [],urlChar,handler); % Octave
    end
catch exception %#ok
    errorid = ['MATLAB:' fcn ':InvalidUrl'];
    errormsg = 'Either this URL could not be parsed or the protocol is not supported.';
    return
end

% Open a connection to the URL.
urlConnection = url.openConnection;

% build up the MATLAB User Agent
mlUserAgent = ['MATLAB R' version('-release') ' '  version('-description')];

% set User-Agent
urlConnection.setRequestProperty('User-Agent', mlUserAgent);
