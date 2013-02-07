% eegplot helper function to read key strokes

function eegplot_readkey(src,evnt)
    if strcmp(evnt.Key, 'rightarrow')==1
        eegplot('drawp',4);
    elseif strcmp(evnt.Key, 'leftarrow')==1
        eegplot('drawp',1);
    end
