% eeg_amplitudearea() - Resamples the averaged ERP with spline interpolation 
%                       at sample rate (resrate) in Hz to get the exact limits 
%                       of the window of integration. High samples the window 
%                       and adds together very narrow rectangles capped by 
%                       right-angled triangles under the window. Output is in uVms. 
%                       trade-off between speed and number of resamples and 
%                       channel selected occurs.
% Usage:
%      >> [channels, amplitude] = eeg_amplitudearea(EEG,channels, resrate, wstart, wend);
% Inputs:
%   channels    - vector of channel indices
%   resrate     - resampling rate for window of integration in Hz
%   wstart      - start of window of integration in ms post stimulus-onset
%   wend        - end of window of integration in ms post stimulus-onset
%   EEG         - EEGLAB data struct containing 3-D epoched data matrix 
%
% Outputs:
%   channels    - a vector of channel indices.
%   amplitude   - 1 dimensional array in uV for the channels
%
% Example
%    >> [channels, amplitude] = eeg_amplitudearea(EEG,[12 18 25 29], 2000, 90.52, 120.52);
%
% Author: Tom Campbell, Helsinki Collegium for Advanced Studies, Biomag Laboratory, 
%         Engineering Centre, Helsinki University Central Hospital Helsinki Brain 
%         Research Centre (tom.campbell@helsinki.fi) Spartam nanctus es: Hanc exorna. 
%         Combined with amplitudearea_msuV() by Darren Weber, UCSF 28/1/05

function [channels,overall_amplitude] = eeg_amplitudearea(EEG, channels, resrate, wstart, wend)

if wstart > wend
    error ('ERROR: wstart must be greater than wend')
else
    [channels, overall_amplitude1] = amplitudearea_msuV (EEG,channels, resrate, wstart, wend)
    overall_amplitude = overall_amplitude1/(wend -wstart)
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [channels, overall_area] = amplitudearea_msuV (EEG, channels, resrate, wstart, wend)

erp = mean(EEG.data,3);
restep = (100/resrate)
[tmp ind1] =min( abs( EEG.times - wstart ) ); % closest time index to wstart
[tmp ind2] =min( abs( EEG.times - wend ) ); % closest time index to wend

if EEG.times(ind1) > wstart
    ind1= ind1 -1;
end

if EEG.times(ind2) < wend
    ind2= ind2 +1;
end

for x= ind1:ind2
    t = (x -ind1)+1;
    tim(t) = EEG.times(x);
end

tr = 1
timr(tr) = wstart;
while timr(tr) < wend
    tr = tr + 1;
    timr(tr) = timr(tr-1)+ restep;
end

for x = 1:size(channels,2)
    channel = channels(x)
    %resamples
    rerp(x, 1:tr) = spline(tim(:),erp(channel, ind1:ind2), timr(1:tr))
    for y = 1:(tr -1)
          % identify which sample is the height of the rectangle under the
          % curve and which sample is the height of the triangle capping the
          % rectangle
        if (abs(rerp(x,y)) <= abs(rerp(x,y+1)))
            rhtsamp = y
            thtsamp = y + 1
        else
            rhtsamp = y + 1
            thtsamp = y
        end
        if y < (tr-1)
            if ((rerp(x,y) > 0)& (rerp(x,y+1) < 0))|((rerp(x,y) < 0)& (rerp(x,y+1) > 0))
              % handles samples containing a zerocrossing with two
              % triangles and trigonometry
                opposite = abs(rerp(x,y))+abs(rerp(x,y+1))
                adjacent = restep
                tantheta = opposite/adjacent
                zerocross = abs(rerp(x,y))/tantheta
                remainder = adjacent - zerocross
                area(x,y) = ((zerocross/2)* (rerp(x,y)))+ ((remainder/2)* (rerp(x,y+1)))
            else
                rectanglearea = rerp(x,(rhtsamp))*restep
                trianglearea = (rerp(x,(thtsamp))-rerp(x,(rhtsamp)))*restep/2
                area(x,y) =  rectanglearea  + trianglearea
            end
        elseif y == (tr-1)
              % special conditions for dealing with last and possibly curtailed
              % sample in window of integration
              % first calculate preconditions to see if zerocrossing occurs
              % prior to end of window of integration in the sample
            endstep =  wend - timr(tr-1)
            opposite = abs(rerp(x,y))+abs(rerp(x,y+1))
            adjacent = restep
            tantheta = opposite/adjacent
            zerocross = abs(rerp(x,y))/tantheta
            if zerocross < endstep
                %handles samples containing a zerocrossing with two
                %triangles and trigonometry
                remainder = endstep - zerocross
                opposite2 = tantheta * remainder * sign(rerp(x,y)) * -1
                area(x,y) = ((zerocross/2)* (rerp(x,y)))+ ((remainder/2)* opposite2)
            elseif rhtsamp > thtsamp
              % tests for amplitude decrease from left to right and fits area with
              % two rectangles capped by a triangle with trigonometry
                opposite = abs(rerp(x,(thtsamp))) - abs(rerp(x,(rhtsamp)))
                adjacent = restep
                tantheta = opposite/adjacent
                endstep = wend - timr(tr- 1)
                excessstep = restep - endstep
                extrarectangleht = tantheta* excessstep
                extrarectanglearea = extrarectangleht * endstep *sign(rerp(x,(thtsamp)))
                triht= opposite - extrarectangleht
                trianglearea = (endstep/2)*triht *sign(rerp(x,(thtsamp)))
                rectanglearea = rerp(x,(rhtsamp))*endstep
                area(x,y) =  rectanglearea + extrarectanglearea + trianglearea
            elseif thtsamp > rhtsamp
              % test for amplitude increase from left to right and fits
              % area with a rectangle capped by a triangle
                opposite = abs(rerp(x,(rhtsamp))) - abs(rerp(x,(rhtsamp)))
                adjacent = restep
                tantheta = opposite/adjacent
                endstep = timr(tr) - wend
                triht = rerp(x,thtsamp)- rerp(x,rhtsamp)
                trianglearea = (endstep/2)*triht
                rectanglearea = rerp(x,(rhtsamp))*endstep
                area(x,y) =  rectanglearea + trianglearea
            end
        end
    end
end


for x = 1:(size(channels,2))
    overall_area(x) = sum(area(x,:))
end

return
