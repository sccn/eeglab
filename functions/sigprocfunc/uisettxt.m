function uisettxt(fig, tag, str, varargin);

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
        set(tmp, 'tag', tag, varargin{:});
    end;
    xlim([0 maxX]);
end;
            