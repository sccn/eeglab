function putverticaltext(fig, tag, str);

alltmpobj = findobj(fig, 'tag', tag);
for tmpobj = alltmpobj'
    tmppos = get(tmpobj, 'position');
    delete(tmpobj);
    axes('position', tmppos);
    axis('off');
    
    if ~iscell(str), str = { str }; end;
    maxX = 0;
    for index = 1:length(str)
        tmp = text(index, 0, str{index});
        maxX = index+1;
        set(tmp, 'rotation', 90, 'fontweight', 'bold', 'tag', tag);
    end;
    xlim([0 maxX]);
end;
            