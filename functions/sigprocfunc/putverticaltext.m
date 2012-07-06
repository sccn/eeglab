function putverticaltext(fig, tag, str);

    tmpobj = findobj(fig, 'tag', tag);
    tmppos = get(tmpobj, 'position');
    delete(tmpobj);
    axes('position', tmppos);
    axis('off');
    
    if ~iscell(str), str = { str }; end;
    maxX = 0;
    for index = 1:length(str)
        tmp = text(index, 0, str{index});
        maxX = index+1;
        set(tmp, 'rotation', 90, 'fontweight', 'bold');
    end;
    xlim([0 maxX]);
            