function [channels,overall_amplitude] = eeg_amplitudearea(EEG, channels, resrate, wstart, wend)
% eeg_amplitudearea() - Resamples an ERP average using spline interpolation 
%                       at a new sample rate (resrate) in Hz to get the exact limits 
%                       of the window of integration. Finely samples the window 
%                       and adds together very narrow rectangles capped by 
%                       right-angled triangles under the window. Output is in uV. 
%                       Trade-off between speed and number of resamples and number of 
%                        channels selected occurs.
% Usage:
%      >> [channels, amplitude] = eeg_amplitudearea(EEG,channels, resrate, wstart, wend);
% Inputs:
%   EEG         - EEGLAB data struct containing a (3-D) epoched data matrix 
%   channels    - vector of channel indices
%   resrate     - resampling rate for window of integration in Hz
%   wstart      - start of window of integration in ms post stimulus-onset
%   wend        - end of window of integration in ms post stimulus-onset
%
% Outputs:
%   channels    - a vector of channel indices.
%   amplitude   - 1-dimensional array in uV for the channels
%
% Example
%    >> [channels, amplitude] = eeg_amplitudearea(EEG,[12 18 25 29], 2000, 90.52, 120.52);
%
% Author: Tom Campbell, Helsinki Collegium for Advanced Studies, Biomag Laboratory, 
%         Engineering Centre, Helsinki University Central Hospital Helsinki Brain 
%         Research Centre (tom.campbell@helsinki.fi) Spartam nanctus es: Hanc exorna. 
%         Combined with amplitudearea_msuV() by Darren Weber, UCSF 28/1/05
%         Retested and debugged Tom Campbell 2/2/05
%         Reconceived, factored somewhat, tested and debugged Tom Campbell 13:24 23.3.2005

if wstart > wend
    error ('ERROR: wstart must be greater than wend')
else
    [channels, overall_amplitude] = eeg_amplitudearea_msuV (EEG,channels, resrate, wstart, wend);
    overall_amplitude = overall_amplitude/(wend - wstart);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [channels, overall_area] = eeg_amplitudearea_msuV (EEG, channels, resrate, wstart, wend)

%if ndim(EEG.data) ~= 3
%  error('EEG.data must be 3-D data epochs');
%end
erp = mean(EEG.data,3);

[tmp ind1] =min( abs( EEG.times - wstart ) ); % closest time index to wstart
[tmp ind2] =min( abs( EEG.times - wend ) ); % closest time index to wend
restep = 1/resrate;

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

tr = 1;
timr(tr) = wstart;
while timr(tr) < wend
    tr = tr + 1;
    timr(tr) = timr(tr-1)+ restep;
end

for x = 1:size(channels,2)
    channel = channels(x);
    %resamples
    rerp(x, 1:tr) = spline(tim(:),erp(channel, ind1:ind2), timr(1:tr));
    pent = timr(tr - 1);
    overall_area(x) = 0;
    for y = 1:(tr -1)
        v1 =  rerp(x,(y));
        v2 =  rerp(x,(y+1));
        if ((v1 > 0) && (v2 < 0)) || ((v1 < 0) && (v2 > 0))
            if (y == (tr-1)) && (timr(y+1)> wend)
                area1 = zero_crossing_truncated(v1, v2, restep, wend, pent);
            else    
                area1 = zero_crossing(v1, v2, restep);
            end
        else
            if( y == (tr-1)) && (timr(y+1)> wend)
                area1 = rect_tri_truncated(v1, v2, restep,wend,pent);
            else
                area1 = rect_tri(v1, v2, restep);
            end
        end
        overall_area(x) = overall_area(x) + area1;
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [area] = zero_crossing(v1,v2,step)
    if (v1 > v2)
        T1 = v1;
        T2 = v2;
    else
        T1 = v2;
        T2 = v1;
    end
    tantheta = (abs(T1)+ abs(T2))/step;
    if (v1 > v2)
        %decline
        z = abs(T1)/tantheta;
        tr1= abs(T1)*(z/2);
        tr2= abs(T2)*((step-z)/2);
    else
        %incline
        z = abs(T2)/tantheta;
        tr2= abs(T2)*(z/2);
        tr1= abs(T1)*((step-z)/2);
    end
    [area] = (tr1 - tr2);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [area] = zero_crossing_truncated(v1,v2,step,wend,pent)
    if (v1 > v2)
        T1 = v1;
        T2 = v2;
    else
        T1 = v2;
        T2 = v1;
    end
    tantheta = (abs(T1)+ abs(T2))/step;
    s = wend - pent;
    if (v1 > v2)
        z = abs(T1)/tantheta
        if s < z   
            %decline,truncated before zerocrossing
            t1 = tantheta * s;
            r1 = abs(T1)-abs(t1);
            tr1= abs(t1)*(s/2);
            tr2= 0;
            rect1 = r1*s;
            rect2 = 0;
        else
            %decline,truncated after zerocrossing
            t2= tantheta*(s-z);
            tr1= abs(T1)*(z/2);
            tr2 = abs(t2)*((s-z)/2);
            rect1 = 0;
            rect2 = 0;
        end    
    else
        z = abs(T2)/tantheta;
        if s < z
            %incline,truncated before zerocrossing
            t2 = tantheta * s;
            r2 = abs(T2)-abs(t2);
            tr1= 0;
            tr2= abs(t2)*(s/2);
            rect1 = 0;
            rect2 = r2*s;
        else
            %incline,truncated after zerocrossing
            t1= tantheta*(s-z);
            tr1 = abs(t1)*((s-z)/2);
            tr2 = abs(T2) * (z/2);
            rect1 = 0;
            rect2 = 0;
        end
    end

[area] = ((rect1 + tr1) - (rect2 + tr2));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [area] = rect_tri(v1,v2,step)
    if (abs(v1) > abs(v2))
        T = abs(v1)-abs(v2);
        R = abs(v2);
    else
        T = abs(v2)-abs(v1);
        R = abs(v1);
    end
    rect = R*step;
    tri = T*(step/2);
    if v1 > 0
        area = 1* (rect+tri);
    else
        area = -1 * (rect+tri);
    end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [area] = rect_tri_truncated(v1,v2,step,wend,pent)
    if (abs(v1) > abs(v2))
        T = abs(v1)-abs(v2);
        R = abs(v2);
    else
        T = abs(v2)-abs(v1);
        R = abs(v1);
    end
    
    tantheta = abs(T)/step;
    s = wend -pent;
    
    if (v1>0)
        if v1 >v2
        %positive decline
            t = tantheta*s;
            e = abs(T)-abs(t);
            rect = s*R;
            exrect = s*e;
            tri = (s/2)*R;
        else
        %positive incline
            t = tantheta*s;
            rect = s*R;
            exrect = 0;
            tri = (s/2)*R;
        end
    else
       if v1 >v2
        %negative decline
            t = tantheta*s;
            rect = s*R;
            exrect = 0;
            tri = (s/2)*R;
        else
        %negative incline 
            t = tantheta*s;
            e = abs(T)-abs(t);
            rect = s*R;
            exrect = s*e;
            tri = (s/2)*R;
        end
    end
    tri = T*(step/2);
    if v1 > 0
        area = 1* (rect+exrect+tri);
    else
        area = -1 * (rect+exrect+tri);
    end
return
