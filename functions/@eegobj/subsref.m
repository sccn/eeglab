% support function for eegobj

function varargout = subsref(this,s)

if strcmpi(s(1).type, '()') && length(this) > 1
    %if length(s(1).subs{1})
    %    error('Unsuported object feature - a.field or a([x y]).field is not supported for object arrays');
    %end
    
%     if length(s) == 1
%         b = this(s(1).subs{1});
%     else
%         b = builtin('subsref', this(s(1).subs{1}).EEG, s(2:end));
%     end
    
    if length(s) == 1
        varargout{1} = this(s(1).subs{1});
    elseif length(s(1).subs{1}) > 1
        for index = 1:length(s(1).subs{1})
            varargout{1}{index} = builtin('subsref', this(s(1).subs{1}(index)).EEG, s(2:end));
        end
    else  
        varargout{1} = builtin('subsref', this(s(1).subs{1}).EEG, s(2:end));
    end
elseif strcmpi(s(1).type, '.') && length(this) > 1
    for index = 1:length(this)
        varargout{index} = builtin('subsref', this(index).EEG, s);
    end
    %if length(s) > 1 && iscell(s(2).subs) && length(s(2).subs{1})
    %    error('Unsuported object feature - a.field or a([x y]).field is not supported for object arrays');
    %end
else
    varargout{1} = builtin('subsref', this.EEG, s);
end


% return;
% 
% % SUBSREF 
% switch s.type
% case '()'
%     b = getfield(a.rawbit, s.subs);
% case '.'
%     b = getfield(struct(a), s.subs);
% otherwise
%    error('Wrong class argument')
% end
