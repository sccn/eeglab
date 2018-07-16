function ncopies = checkworkspace(obj);

stack = dbstack;
stack(1:2) = [];
stack      = rmfield(stack, 'line');
ncopies = 0;
if ~isempty(stack)
    % empty stack means the base workspace
    % check if we are in a different workspace
    if ~isequal(stack, obj.workspace)
        % is the current stack a subset of the obj workspace?
        % if yes, it means that the object was created in a subfunction
        % and that it is OK to modify it (because it does not exist in the
        % subfunction any more)
        subFlag = false;
        for index = 1:length(obj.workspace)
            if isequal(obj.workspace(index), stack(1))
                subFlag = true;
                if length(stack) > 1
                    for index2 = 1:length(obj.workspace)-index
                        if length(stack) < index2+1
                            subFlag = false;
                        elseif ~isequal(obj.workspace(index+index2), stack(index2+1))
                            subFlag = false;
                        end
                    end
                end
                if subFlag, return; end
            end
        end
                        
        % if subfunction, must be a copie
        if ~isempty(obj.workspace) && strcmpi(stack(end).file, obj.workspace(end).file) && ...
                ~strcmpi(stack(end).name, obj.workspace(end).name)
            % We are within a subfunction. The MMO must have
            % been passed as an argument (otherwise the current
            % workspace and the workspace variable would be
            % equal).
            ncopies = 2;
        else
            if ~isscript(stack(1).file)
                ncopies = 2;
                % we are within a function. The MMO must have
                % been passed as an argument (otherwise the current
                % workspace and the workspace variable would be
                % equal).
            else
                % we cannot be in a function with 0 argument
                % (otherwise the current workspace and the workspace
                % variable would be equal). We must assume that
                % we are in a script.
                while ~isempty(stack) && ~isequal(stack, obj.workspace)
                    stack(1) = [];
                end
                if ~isequal(stack, obj.workspace)
                    ncopies = 2;
                end
            end
        end
    end
end
