function [polydata] = ERPLAB_polydetrend(data,window, orderpoly)

%function [polydata] = ERPLAB_polydetrend(data,window)
%
% Inputs:
%
%   EEG       - input dataset
%   window    - minimun time window (in sec) to represent the
%               mean of DC behaviour (into the window).
%               Row data will divided per this window, generating
%               a joint of points over which a polynomial will be
%               calculated.
%               Applying the same vector of time to polynomial we
%               will get the full DC behavior per channel. The
%               last one is subtracted from each channel finally.
%
%
% Outputs:
%
%   polydata   - output dataset. Detrending data
%
%     Author: ERPLAB Team, Center for Mind & Brain
%             Universidad de California, Davis. 2007

if nargin < 1
	help ERPLAB_polydetrend
	return
end

if exist('filtfilt','file') ~= 2
    disp('ERPLAB_polydetrend error: cannot find the Signal Processing Toolbox');
    return
end

if isempty(data)
    disp('ERPLAB_polydetrend error: cannot filter an empty dataset')
    return
end

if nargin < 3
    disp('ERPLAB_polydetrend error: please, enter all arguments!')
    return
end

[numchan frames] = size(data);
nwin=round(frames/window); % a priori

% calculas los puntos representativos por ventana
xf = linspace(1,nwin,frames);
polydata = 0*data;
ss = zeros(numchan,nwin);

for i = 1:numchan
    a=1;
    b=window;
    for j=1:nwin
        if b <= frames
        ss(i,j) = mean(data(i,a:b));
        a = b + 1;
        b = b + window;
        else
            ss(i,j) = ss(i,j-1);
        end
    end
    p = polyfit(1:nwin,ss(i,:),orderpoly);
    % Poner otros metodos aqui...
    polydata(i,:) = data(i,:) - polyval(p,xf);
    %polydata(i,:) = polyval(p,xf);

end


% Otros metodos...
    
    %polydata(i,:) = spline(1:nwin,ss(i,:),xf);
    %polydata(i,:) = interp1(1:nwin,ss(i,:),xf);
