% pop_eventinfo() - GUI for BIDS event info editing, generated based on
%                   fields of EEG.event
%
% Usage:
%   >> [ALLEEG, eInfoDesc, eInfo] = pop_eventinfo( ALLEEG );
%                                              
% Inputs:
%   ALLEEG        - ALLEEG dataset structure. May only contain one dataset.
%
% Outputs:
%  'ALLEEG'       - [struct] Updated ALLEEG structure containing event BIDS information
%                in each EEG structure at EEG.BIDS.eInfoDesc and EEG.BIDS.eInfo
%
%  'eInfoDesc' - [struct] structure describing BIDS event fields as you specified.
%                See BIDS specification for all suggested fields.
%
%  'eInfo'     - [cell] BIDS event fields and their corresponding
%                event fields in the EEGLAB event structure. Note that
%                EEGLAB event latency, duration, and type are inserted
%                automatically as columns "onset" (latency in sec), "duration"
%                (duration in sec), "value" (EEGLAB event type)
%
% Author: Dung Truong, Arnaud Delorme
function [ALLEEG, eInfoDesc, eInfo] = pop_eventinfo(ALLEEG)
    % perform check to make sure EEG.event is consistent across ALLEEG
    try
       eventFields = fieldnames([ALLEEG.event]);
    catch ME
        if (strcmp(ME.identifier, 'MATLAB:catenate:structFieldBad'))
            numFields = cellfun(@(x) numel(fieldnames(x)), {ALLEEG.event});
            [num, idx] = max(numFields);
            eventFields = fieldnames(ALLEEG(idx).event);
            warning('There is mismatch in number of fields in EEG.event structures. Using fields of ALLEEG(%d) which has the highest number of fields (%d).', idx, num);
        end
    end
    % default settings
    appWidth = 800;
    appHeight = 500;
    bg = [0.65 0.76 1];
    fg = [0 0 0.4];
    levelThreshold = 20;
    fontSize = 12;
    eventBIDS = newEventBIDS();
    columnDefinition.LongName = 'Long, unabbreviated name of the field';
    columnDefinition.Description = 'Description of the field';
    columnDefinition.Levels = 'For categorical variables: possible values and their descriptions';
    columnDefinition.Units = 'Measurement units - format [<prefix>]<name>';
    columnDefinition.TermURL = 'URL pointing to a formal definition of this type of data in an ontology available on the web';
    
    eInfo = {};
    eInfoDesc = [];
    
    % create UI
    
    f = figure('MenuBar', 'None', 'ToolBar', 'None', 'Name', 'Edit BIDS event info - pop_eventinfo', 'Color', bg);
    f.Position(3) = appWidth;
    f.Position(4) = appHeight;
    uicontrol(f, 'Style', 'text', 'String', 'BIDS information for EEG.event fields', 'Units', 'normalized','FontWeight','bold','ForegroundColor', fg,'BackgroundColor', bg, 'Position', [0 0.9 1 0.1]);
    tbl = uitable(f, 'RowName', eventFields, 'ColumnName', { 'BIDS Field' 'Levels' 'LongName' 'Description' 'Unit Name' 'Unit Prefix' 'TermURL' }, 'Units', 'normalized', 'FontSize', fontSize, 'Tag', 'bidsTable');
    tbl.Position = [0.01 0.54 0.98 0.41];
    tbl.CellSelectionCallback = @cellSelectedCB;
    tbl.CellEditCallback = @cellEditCB;
    tbl.ColumnEditable = [true false true true true true];
    tbl.ColumnWidth = {appWidth/9,appWidth/9,appWidth*2/9,appWidth*2/9,appWidth/9,appWidth/9,appWidth/9};
    units = {' ','ampere','becquerel','candela','coulomb','degree Celsius','farad','gray','henry','hertz','joule','katal','kelvin','kilogram','lumen','lux','metre','mole','newton','ohm','pascal','radian','second','siemens','sievert','steradian','tesla','volt','watt','weber'};
    unitPrefixes = {' ','deci','centi','milli','micro','nano','pico','femto','atto','zepto','yocto','deca','hecto','kilo','mega','giga','tera','peta','exa','zetta','yotta'};
    tbl.ColumnFormat = {[] [] [] [] units unitPrefixes []};
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Ok', 'Units', 'normalized', 'Position', [0.85 0 0.1 0.05], 'Callback', @okCB); 
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Cancel', 'Units', 'normalized', 'Position', [0.7 0 0.1 0.05], 'Callback', @cancelCB); 
    
    % pre-populate table
    data = cell(length(eventFields),length(tbl.ColumnName));
    for i=1:length(eventFields)
        % pre-populate description
        field = eventFields{i};
        if isfield(eventBIDS, field)
            if isfield(eventBIDS.(field), 'BIDSField')
                data{i,1} = eventBIDS.(field).BIDSField;
            end
            if isfield(eventBIDS.(field), 'LongName')
                data{i,find(strcmp(tbl.ColumnName, 'LongName'))} = eventBIDS.(field).LongName;
            end
            if isfield(eventBIDS.(field), 'Description')
                data{i,find(strcmp(tbl.ColumnName, 'Description'))} = eventBIDS.(field).Description;
            end
            if isfield(eventBIDS.(field), 'Units')
                data{i,find(strcmp(tbl.ColumnName, 'Unit Name'))} = eventBIDS.(field).Units;
            end
            if isfield(eventBIDS.(field), 'TermURL')
                data{i,find(strcmp(tbl.ColumnName, 'TermURL'))} = eventBIDS.(field).TermURL;
            end
            if isfield(eventBIDS.(field), 'Levels') && ~isempty(eventBIDS.(field).Levels)
                data{i,find(strcmp(tbl.ColumnName, 'Levels'))} = strjoin(fieldnames(eventBIDS.(field).Levels),',');
            else
                if strcmp(field, 'latency') || strcmp(field, "usertags")
                    data{i,find(strcmp(tbl.ColumnName, 'Levels'))} = 'n/a';
                else
                    data{i,find(strcmp(tbl.ColumnName, 'Levels'))} = 'Click to specify below';
                end
            end
        end
        clear('field');
    end
    tbl.Data = data;
    waitfor(f);
  
    function cancelCB(src, event)
        clear('eventBIDS');
        close(f);
    end
    function okCB(src, event)
        % duration field is automatically calculated by EEGLAB
        eInfoDesc.duration.LongName = 'Event duration';
        eInfoDesc.duration.Description = 'Duration of the event (measured from onset) in seconds';
        eInfoDesc.duration.Units = 'second';
        
        % prepare return struct
        fields = fieldnames(eventBIDS);
        for idx=1:length(fields)
            eegField = fields{idx};
            bidsField = eventBIDS.(eegField).BIDSField;
            if ~isempty(bidsField)
                eInfo = [eInfo; {bidsField eegField}]; 
                
%                 eInfoDesc.(bidsField).EEGField = eegField;
                if ~isempty(eventBIDS.(eegField).LongName)
                    eInfoDesc.(bidsField).LongName = eventBIDS.(eegField).LongName;
                end
                if ~isempty(eventBIDS.(eegField).Description)
                    eInfoDesc.(bidsField).Description = eventBIDS.(eegField).Description;
                end
                if ~isempty(eventBIDS.(eegField).Units)
                    eInfoDesc.(bidsField).Units = eventBIDS.(eegField).Units;
                end
                if ~isempty(eventBIDS.(eegField).Levels)
                    eInfoDesc.(bidsField).Levels = eventBIDS.(eegField).Levels;
                end
                if ~isempty(eventBIDS.(eegField).TermURL)
                    eInfoDesc.(bidsField).TermURL = eventBIDS.(eegField).TermURL;
                end
            end
        end
        for e=1:numel(ALLEEG)
            ALLEEG(e).BIDS.eInfoDesc = eInfoDesc;
            ALLEEG(e).BIDS.eInfo = eInfo;
        end
        clear('eventBIDS');
        close(f);
    end

    function cellSelectedCB(arg1, obj) 
        if size(obj.Indices,1) == 1
            removeLevelUI();
            row = obj.Indices(1);
            col = obj.Indices(2);
            field = obj.Source.RowName{row};
            bfield = obj.Source.Data{row,1};
            columnName = obj.Source.ColumnName{col};
            
            if strcmp(columnName, 'BIDS Field')
                c6 = uicontrol(f, 'Style', 'text', 'String', sprintf('Enter BIDS field or choose one of the suggested below:'), 'Units', 'normalized', 'FontAngle','italic','ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'selectBIDSMsg');
                c6.Position = [0.01 0.44 1 0.05];
                c6.HorizontalAlignment = 'left';
                c = uicontrol(f,'Style','popupmenu', 'Units', 'normalized', 'Tag', 'selectBIDSDD');
                c.Position = [0.01 0.34 0.3 0.1];
                curBIDSFields = {obj.Source.Data{:,col}};
                predefined_bids_fields = setdiff({'onset', 'trial_type','value','stim_file','sample','HED','response_time'}, curBIDSFields(~cellfun(@isempty, curBIDSFields)));
                c.String = ['Choose BIDS predefined field' predefined_bids_fields];
                c.Callback = {@bidsFieldSelected, obj.Source, row, col};
            else % any other column selected
                if strcmp(columnName, 'Levels')
                    createLevelUI('','',obj,field);
                elseif strcmp(columnName, 'Description')
                    uicontrol(f, 'Style', 'text', 'String', sprintf('%s (%s):',columnName, columnDefinition.(columnName)), 'Units', 'normalized', 'Position',[0.01 0.44 0.98 0.05], 'HorizontalAlignment', 'left','FontAngle','italic','ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'cellContentHeader');
                    uicontrol(f, 'Style', 'edit', 'String', obj.Source.Data{row,col}, 'Units', 'normalized', 'Max',2,'Min',0,'Position',[0.01 0.24 0.7 0.2], 'HorizontalAlignment', 'left', 'Callback', {@descriptionCB, obj,field}, 'Tag', 'cellContentMsg');
                else
                    if strcmp(columnName, 'Unit Name') || strcmp(columnName, 'Unit Prefix')
                        columnName = 'Units';
                        content = [obj.Source.Data{row,find(strcmp(obj.Source.ColumnName, 'Unit Prefix'))} obj.Source.Data{row,find(strcmp(obj.Source.ColumnName, 'Unit Name'))}];
                    else
                        content = obj.Source.Data{row,col};
                    end
                    % display cell content in lower panel
                    uicontrol(f, 'Style', 'text', 'String', sprintf('%s (%s):',columnName, columnDefinition.(columnName)), 'Units', 'normalized', 'Position',[0.01 0.44 0.98 0.05], 'HorizontalAlignment', 'left','FontAngle','italic','ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'cellContentHeader');
                    uicontrol(f, 'Style', 'text', 'String', sprintf('%s',content), 'Units', 'normalized', 'Position',[0.01 0.39 0.98 0.05], 'HorizontalAlignment', 'left','ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'cellContentMsg');
                end
                if isempty(bfield)
                    c6 = uicontrol(f, 'Style', 'text', 'String', sprintf('Remember to specify BIDS name'), 'Units', 'normalized', 'FontWeight', 'bold', 'FontAngle','italic','ForegroundColor', [0.9 0 0],'BackgroundColor', bg, 'Tag', 'noBidsMsg');
                    c6.Position = [0.01 0.49 0.3 0.05];
                    c6.HorizontalAlignment = 'Left';
                end
            end
        end
    end

    function cellEditCB(arg1, obj)
        field = obj.Source.RowName{obj.Indices(1)};
        column = obj.Source.ColumnName{obj.Indices(2)};
        if ~strcmp(column, 'Levels')
            if strcmp(column, 'BIDS Field')
                otherBIDSFieldsIdx = setdiff(1:size(obj.Source.Data,1),obj.Indices(1));
                curBIDSFields = {obj.Source.Data{otherBIDSFieldsIdx,obj.Indices(2)}};
                if any(strcmp(obj.EditData, curBIDSFields)) % check for duplication of BIDS field
                    obj.Source.Data{obj.Indices(1),obj.Indices(2)} = obj.PreviousData; % reset if found
                else
                    eventBIDS.(field).BIDSField = obj.EditData;
                end
            elseif strcmp(column, 'Unit Name') || strcmp(column, 'Unit Prefix')
                unit = [obj.Source.Data{obj.Indices(1),find(strcmp(obj.Source.ColumnName, 'Unit Prefix'))} obj.Source.Data{obj.Indices(1),find(strcmp(obj.Source.ColumnName, 'Unit Name'))}];
                eventBIDS.(field).Units = unit;
            else
                eventBIDS.(field).(column) = obj.EditData;
            end
        end
    end
    
    function descriptionCB(src,event,obj,field) 
        obj.Source.Data{obj.Indices(1),obj.Indices(2)} = src.String; % reset if found
        eventBIDS.(field).Description = src.String;
    end

    function createLevelUI(src,event,table,field)
        removeLevelUI();
        levelCellText = table.Source.Data{find(strcmp(table.Source.RowName, field)), find(strcmp(table.Source.ColumnName, 'Levels'))}; % text (fieldName-Levels) cell. if 'n/a' then no action, 'Click to..' then conditional action, '<value>,...' then get levels
        if strcmp(field, 'usertags')
            uicontrol(f, 'Style', 'text', 'String', 'Levels editing not applied for HED. Use ''pop_tageeg(EEG)'' of HEDTools plug-in to edit event HED tags', 'Units', 'normalized', 'Position', [0.01 0.45 1 0.05],'ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'levelEditMsg');
        elseif strcmp(field, 'latency')
            uicontrol(f, 'Style', 'text', 'String', 'Levels editing not applied for EEG.event.latency field.', 'Units', 'normalized', 'Position', [0.01 0.45 1 0.05],'ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'levelEditMsg');
        else
            % retrieve all unique values from EEG.event.(field). 
            % Use ALLEEG(1) as representative EEG
            if isnumeric(ALLEEG(1).event(1).(field))
                values = arrayfun(@(x) num2str(x), [ALLEEG(1).event.(field)], 'UniformOutput', false);
                levels = unique(values)';
            else
                levels = unique({ALLEEG(1).event.(field)})';
            end
            if strcmp(levelCellText,'Click to specify below') && length(levels) > levelThreshold 
                msg = sprintf('\tThere are more than %d unique levels for field %s.\nAre you sure you want to specify levels for it?', levelThreshold, field);
                c4 = uicontrol(f, 'Style', 'text', 'String', msg, 'Units', 'normalized', 'FontWeight', 'bold', 'ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'confirmMsg');
                c4.Position = [0 0.38 1 0.1];
                c5 = uicontrol(f, 'Style', 'pushbutton', 'String', 'Yes', 'Units', 'normalized','Tag', 'confirmBtn', 'Callback', {@ignoreThresholdCB,table,field});
                c5.Position = [0.5-c5.Extent(3)/2 0.33 0.1 0.05];
            else
                % build table data
                t = cell(length(levels),2);
                for lvl=1:length(levels)
                    formattedLevel = checkFormat(levels{lvl}); % put level in the right format for indexing. Number is prepended by 'x'
                    t{lvl,1} = formattedLevel;
                    if ~isempty(eventBIDS.(field).Levels) && isfield(eventBIDS.(field).Levels, formattedLevel)
                        t{lvl,2} = eventBIDS.(field).Levels.(formattedLevel);
                    end
                end
                % create UI
                uicontrol(f, 'Style', 'text', 'String', ['Describe the categorical values of EEG.event.' field], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'Position', [0.31 0.45 0.7 0.05],'FontWeight', 'bold','ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'levelEditMsg');
                msg = 'BIDS allow you to describe the level for each of your categorical field. Describing levels help other researchers understand your experiment better';
                uicontrol(f, 'Style', 'text', 'String', msg, 'Units', 'normalized', 'HorizontalAlignment', 'Left','Position', [0.01 0 0.3 0.4],'ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'levelMsg');
                h = uitable(f, 'Data', t(:,2), 'ColumnName', {'Description'}, 'RowName', t(:,1), 'Units', 'normalized', 'Position', [0.31 0.07 0.68 0.38], 'FontSize', fontSize, 'Tag', 'levelEditTbl', 'CellEditCallback',{@levelEditCB,field},'ColumnEditable',true); 
                h.ColumnWidth = {appWidth*0.68*0.8};
            end
        end
    end
    
    function ignoreThresholdCB(src,event,table, field)
        table.Source.Data{find(strcmp(table.Source.RowName, field)), find(strcmp(table.Source.ColumnName, 'Levels'))} = 'Click to specify below (ignore max number of levels threshold)';
        createLevelUI('','',table,field);
    end
    
    % Event handler for level edit
    function levelEditCB(arg1, obj, field)
        level = checkFormat(obj.Source.RowName{obj.Indices(1)});
        description = obj.EditData;
        % update eventBIDS structure
        eventBIDS.(field).Levels.(level) = description;
        specified_levels = fieldnames(eventBIDS.(field).Levels);
        % Update main table
        mainTable = findobj('Tag','bidsTable');
        mainTable.Data{find(strcmp(field,mainTable.RowName)),find(strcmp('Levels',mainTable.ColumnName))} = strjoin(specified_levels, ',');
    end
    
    function bidsFieldSelected(src, event, table, row, col) 
        val = src.Value;
        str = src.String;
        selected = str{val};
        table.Data{row,col} = selected;
        field = table.RowName{row};
        eventBIDS.(field).BIDSField = selected;
    end
    function formatted = checkFormat(str)
        if ~isempty(str2num(str))
            formatted = ['x' str];
        else
            formatted = strrep(str,' ','_'); %replace space with _
        end
    end
    function removeLevelUI()
        % remove old ui items of level section if exist
        h = findobj('Tag', 'levelEditMsg');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'levelEditTbl');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'confirmMsg');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'confirmBtn');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'noLevelBtn');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'noBidsMsg');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'cellContentMsg');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'cellContentHeader');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'selectBIDSMsg');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'selectBIDSDD');
        if ~isempty(h)
            delete(h);
        end
        h = findobj('Tag', 'levelMsg');
        if ~isempty(h)
            delete(h);
        end
    end

    % generate 
    function event = newEventBIDS(eegIdx)
        event = [];
        bidsEEG = [];
        if isfield(ALLEEG,'BIDS') 
            bidsIdx = find(~cellfun(@isempty,{ALLEEG.BIDS}));
            if ~isempty(bidsIdx)
                if numel(ALLEEG) ~= numel(bidsIdx)
                    fprintf(['EEG.BIDS is found in ' num2str(numel(bidsIdx)) ' out of ' num2str(numel(ALLEEG)) ' EEG structure(s). ']);
                end
                fprintf(['Using BIDS info of ALLEEG(' num2str(bidsIdx(1)) ')...']);
                bidsEEG = ALLEEG(bidsIdx(1));
            end
        end
        
        % if resume editing
        if ~isempty(bidsEEG) && isfield(bidsEEG.BIDS,'eInfoDesc') && isfield(bidsEEG.BIDS,'eInfo')

            for idx=1:size(bidsEEG.BIDS.eInfo,1)
                field = bidsEEG.BIDS.eInfo{idx,2}; 
                bids_field = bidsEEG.BIDS.eInfo{idx,1};
                event.(field).BIDSField = bids_field;
                if isfield(bidsEEG.BIDS.eInfoDesc,bids_field) && isfield(bidsEEG.BIDS.eInfoDesc.(bids_field), 'LongName')
                    event.(field).LongName = bidsEEG.BIDS.eInfoDesc.(bids_field).LongName;
                else
                    event.(field).LongName = '';
                end
                if isfield(bidsEEG.BIDS.eInfoDesc,bids_field) && isfield(bidsEEG.BIDS.eInfoDesc.(bids_field), 'Description')
                    event.(field).Description = bidsEEG.BIDS.eInfoDesc.(bids_field).Description;
                else
                    event.(field).Description = '';
                end
                if isfield(bidsEEG.BIDS.eInfoDesc,bids_field) && isfield(bidsEEG.BIDS.eInfoDesc.(bids_field), 'Units')
                    event.(field).Units = bidsEEG.BIDS.eInfoDesc.(bids_field).Units;
                else
                    event.(field).Units = '';
                end
                if isfield(bidsEEG.BIDS.eInfoDesc,bids_field) && isfield(bidsEEG.BIDS.eInfoDesc.(bids_field), 'Levels')
                    event.(field).Levels = bidsEEG.BIDS.eInfoDesc.(bids_field).Levels;
                else
                    event.(field).Levels = [];
                end
                if isfield(bidsEEG.BIDS.eInfoDesc,bids_field) && isfield(bidsEEG.BIDS.eInfoDesc.(bids_field), 'TermURL')
                    event.(field).TermURL = bidsEEG.BIDS.eInfoDesc.(bids_field).TermURL;
                else
                    event.(field).TermURL = '';
                end
            end
            fields = setdiff(eventFields, {bidsEEG.BIDS.eInfo{:,2}}); % unset fields
            for idx=1:length(fields)
                event.(fields{idx}).BIDSField = '';
                event.(fields{idx}).LongName = '';
                event.(fields{idx}).Description = '';
                event.(fields{idx}).Units = '';
                event.(fields{idx}).Levels = [];
                event.(fields{idx}).TermURL = '';
            end

            clear('EEG');
        else % start fresh
            fields = eventFields; 
            for idx=1:length(fields)
                if strcmp(fields{idx}, 'type')
                    event.type.BIDSField = 'value';
                    event.type.LongName = 'Event marker';
                    event.type.Description = 'Marker value associated with the event';
                    event.type.Units = '';
                    event.type.Levels = [];
                    event.type.TermURL = '';
                elseif strcmp(fields{idx}, 'usertags')
                    event.usertags.BIDSField = 'HED';
                    event.usertags.LongName = 'Hierarchical Event Descriptor';
                    event.usertags.Description = 'Tags describing the nature of the event';      
                    event.usertags.Levels = [];
                    event.usertags.Units = '';
                    event.usertags.TermURL = '';
                elseif strcmp(fields{idx}, 'latency')
                    event.latency.BIDSField = 'onset';
                    event.latency.LongName = 'Event onset';
                    event.latency.Description = 'Onset (in seconds) of the event measured from the beginning of the acquisition of the first volume in the corresponding task imaging data file';
                    event.latency.Units = 'second';
                    event.latency.Levels = [];
                    event.latency.TermURL = '';
                else
                    event.(fields{idx}).BIDSField = '';
                    event.(fields{idx}).LongName = '';
                    event.(fields{idx}).Description = '';
                    event.(fields{idx}).Units = '';
                    event.(fields{idx}).Levels = [];
                    event.(fields{idx}).TermURL = '';
                end
            end
        end
    end
end