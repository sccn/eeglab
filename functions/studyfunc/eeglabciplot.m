% eeglabciplot(lower,upper)       
% eeglabciplot(lower,upper,x)
% eeglabciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Author:   Raymond Reynolds 24/11/06
% Modified by Ramon Martinez Cancino
function eeglabciplot(lower,upper,x,colour, alphaval)

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<5
    alphaval = 1;
end
if nargin<4
    colour= 'r';
end

if nargin<3
    x=1:length(lower);
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

hdl_tmp = fill([x fliplr(x)],[upper fliplr(lower)],colour,'FaceAlpha', alphaval, 'EdgeColor', 'none');