% pop_participantinfo() - GUI for BIDS participant info editing
%
% Usage:
%   >> [STUDY, pInfoDesc, pInfo] = pop_participantinfo( STUDY );
%                                              
% Inputs:
%   STUDY        - EEGLABE STUDY structure. May only contain one dataset.
%
% Outputs:
%  'STUDY'       - [struct] Updated STUDY structure containing event BIDS information
%                in STUDY.BIDS.pInfoDesc and STUDY.BIDS.pInfo
%
%  'pInfoDesc' - [struct] structure describing BIDS participant fields as you specified.
%                See BIDS specification for all suggested fields.
%
%  'pInfo'     - [cell] BIDS participant information.
%
% Author: Dung Truong, Arnaud Delorme
function [STUDY, pInfoDesc, pInfo] = pop_participantinfo(STUDY)
    %% default settings
    appWidth = 1000;
    appHeight = 500;
    bg = [0.65 0.76 1];
    fg = [0 0 0.4];
    levelThreshold = 20;
    fontSize = 12;
    pFields = { 'Participant_id' 'Sex' 'Age' 'Group' };
    pInfoBIDS = newpInfoBIDS();    
    pInfo = {};
    pInfoDesc = [];
    
    %% create UI
    f = figure('MenuBar', 'None', 'ToolBar', 'None', 'Name', 'Edit BIDS event info - pop_eventinfo', 'Color', bg);
    f.Position(3) = appWidth;
    f.Position(4) = appHeight;
    uicontrol(f, 'Style', 'text', 'String', 'Participant information', 'Units', 'normalized','FontWeight','bold','ForegroundColor', fg,'BackgroundColor', bg, 'Position', [0 0.9 0.3 0.1]);
    pInfoTbl = uitable(f, 'RowName', 'numbered', 'ColumnName', pFields, 'Units', 'normalized', 'FontSize', fontSize, 'Tag', 'pInfoTable', 'ColumnEditable', true);
    pInfoTbl.Position = [0.01 0.01 0.3 0.947];
%     pInfoTbl.Data = cell(500, length(pFields));
    % pre-populate pInfo table
    if isfield(STUDY, 'BIDS') && isfield(STUDY.BIDS, 'pInfo')
        pInfoTbl.Data = STUDY.BIDS.pInfo(2:end,:);
    else
        subjects = unique({STUDY.datasetinfo.subject});
        pInfoTbl.Data = cell(length(subjects), length(pFields));
        for i=1:length(subjects)
            pInfoTbl.Data{i,1} = subjects{i};
        end
    end
    
    uicontrol(f, 'Style', 'text', 'String', 'BIDS metadata for participants fields', 'Units', 'normalized','FontWeight','bold','ForegroundColor', fg,'BackgroundColor', bg, 'Position', [0.31 0.9 0.8 0.1]);
    tbl = uitable(f, 'RowName', pFields, 'ColumnName', {'Description' 'Levels' 'Units' }, 'Units', 'normalized', 'FontSize', fontSize, 'Tag', 'bidsTable');
    tbl.Position = [0.32 0.54 0.67 0.41];
    tbl.CellSelectionCallback = @cellSelectedCB;
    tbl.CellEditCallback = @cellEditCB;
    tbl.ColumnEditable = [true false true];
    tbl.ColumnWidth = {appWidth*0.67*2/5,appWidth*0.67/5,appWidth*0.67/5};
    units = {' ','ampere','becquerel','candela','coulomb','degree Celsius','farad','gray','henry','hertz','joule','katal','kelvin','kilogram','lumen','lux','metre','mole','newton','ohm','pascal','radian','second','siemens','sievert','steradian','tesla','volt','watt','weber'};
    unitPrefixes = {' ','deci','centi','milli','micro','nano','pico','femto','atto','zepto','yocto','deca','hecto','kilo','mega','giga','tera','peta','exa','zetta','yotta'};
    tbl.ColumnFormat = {[] [] [] [] units unitPrefixes []};
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Ok', 'Units', 'normalized', 'Position', [0.85 0 0.1 0.05], 'Callback', @okCB); 
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Cancel', 'Units', 'normalized', 'Position', [0.7 0 0.1 0.05], 'Callback', @cancelCB); 
    
    % pre-populate BIDS table
    data = cell(length(pFields),length(tbl.ColumnName));
    for i=1:length(pFields)
        % pre-populate description
        field = pFields{i};
        data{i,find(strcmp(tbl.ColumnName, 'Levels'))} = pInfoBIDS.(field).Levels;
        data{i,find(strcmp(tbl.ColumnName, 'Description'))} = pInfoBIDS.(field).Description;
        data{i,find(strcmp(tbl.ColumnName, 'Units'))} = pInfoBIDS.(field).Units;
        clear('field');
    end
    tbl.Data = data;
    %% wait
    waitfor(f);
  
    %% callback handle for cancel button
    function cancelCB(src, event)
        clear('eventBIDS');
        close(f);
    end

    %% callback handle for ok button
    function okCB(src, event)        
        % prepare return struct
        pTable = findobj('Tag', 'pInfoTable');
        idx = cellfun(@isempty,{pTable.Data{:,1}});
        pInfo = [pFields; pTable.Data(~idx,:)];
        
        fields = fieldnames(pInfoBIDS);
        for idx=1:length(fields)
            field = fields{idx};
            if ~isempty(pInfoBIDS.(field).Description)
                pInfoDesc.(field).Description = pInfoBIDS.(field).Description;
            end
            if ~isempty(pInfoBIDS.(field).Units)
                pInfoDesc.(field).Units = pInfoBIDS.(field).Units;
            end
            if ~isempty(pInfoBIDS.(field).Levels)
                pInfoDesc.(field).Levels = pInfoBIDS.(field).Levels;
            end
        end
        
        STUDY.BIDS.pInfoDesc = pInfoDesc;
        STUDY.BIDS.pInfo = pInfo;
        clear('pInfoBIDS');
        close(f);
    end

    %% callback handle for cell selection in the BIDS table
    function cellSelectedCB(arg1, obj) 
        if size(obj.Indices,1) == 1
            removeLevelUI();
            row = obj.Indices(1);
            col = obj.Indices(2);
            field = obj.Source.RowName{row};
            columnName = obj.Source.ColumnName{col};
            
            if strcmp(columnName, 'Levels')
                createLevelUI('','',field);
            elseif strcmp(columnName, 'Description')
                uicontrol(f, 'Style', 'text', 'String', sprintf('%s (%s):',columnName, 'Full description of the field'), 'Units', 'normalized', 'Position',[0.32 0.49 0.68 0.05], 'HorizontalAlignment', 'left','FontAngle','italic','ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'cellContentHeader');
                uicontrol(f, 'Style', 'edit', 'String', obj.Source.Data{row,col}, 'Units', 'normalized', 'Max',2,'Min',0,'Position',[0.32 0.29 0.5 0.2], 'HorizontalAlignment', 'left', 'Callback', {@descriptionCB, obj,field}, 'Tag', 'cellContentMsg');
            end
        end
    end
    
    %% callback handle for cell edit in BIDS table
    function cellEditCB(arg1, obj)
        field = obj.Source.RowName{obj.Indices(1)};
        column = obj.Source.ColumnName{obj.Indices(2)};
        if ~strcmp(column, 'Levels')
            pInfoBIDS.(field).(column) = obj.EditData;
        end
    end
    
    function descriptionCB(src,event,obj,field) 
        obj.Source.Data{obj.Indices(1),obj.Indices(2)} = src.String;
        pInfoBIDS.(field).Description = src.String;
    end

    function createLevelUI(src,event,field)
        removeLevelUI();
        
        if strcmp(field, 'Participant_id') || strcmp(field, 'Age')
            uicontrol(f, 'Style', 'text', 'String', 'Levels editing not applied.', 'Units', 'normalized', 'Position', [0.31 0.45 0.68 0.05],'ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'levelEditMsg');
        else
            pTable = findobj('Tag', 'pInfoTable');
            colIdx = find(strcmp(pTable.ColumnName, field));
            % retrieve all unique values from EEG.event.(field)
            if isnumeric(pTable.Data{1,colIdx})
                values = arrayfun(@(x) num2str(x), [pTable.Data{:,colIdx}], 'UniformOutput', false);
                levels = unique(values)';
            else
                values = {pTable.Data{:,colIdx}};
                idx = cellfun(@isempty, values);
                levels = unique(values(~idx))';
            end
            if isempty(levels)
                uicontrol(f, 'Style', 'text', 'String', 'No value found. Please specify values in Participant information table.', 'Units', 'normalized', 'Position', [0.31 0.45 0.68 0.05],'ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'levelEditMsg');
            elseif length(levels) > levelThreshold 
                msg = sprintf('\tThere are more than %d unique levels for field %s.\nAre you sure you want to specify levels for it?', levelThreshold, field);
                c4 = uicontrol(f, 'Style', 'text', 'String', msg, 'Units', 'normalized', 'FontWeight', 'bold', 'ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'confirmMsg');
                c4.Position = [0.31 0.38 0.68 0.1];
                c5 = uicontrol(f, 'Style', 'pushbutton', 'String', 'Yes', 'Units', 'normalized','Tag', 'confirmBtn', 'Callback', {@createLevelUI,field});
                c5.Position = [0.31+(1-c5.Extent(3))/2 0.33 0.68 0.05];
            else
                % build table data
                t = cell(length(levels),2);
                for lvl=1:length(levels)
                    formattedLevel = checkFormat(levels{lvl}); % put level in the right format for indexing. Number is prepended by 'x'
                    t{lvl,1} = formattedLevel;
                    if ~isempty(pInfoBIDS.(field).Levels) && isfield(pInfoBIDS.(field).Levels, formattedLevel)
                        t{lvl,2} = pInfoBIDS.(field).Levels.(formattedLevel);
                    end
                end
                % create UI
                uicontrol(f, 'Style', 'text', 'String', ['Describe the categorical values of participant field ' field], 'Units', 'normalized', 'HorizontalAlignment', 'left', 'Position', [0.48 0.45 0.53 0.05],'FontWeight', 'bold','ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'levelEditMsg');
                msg = 'BIDS allows you to describe the level for each of your categorical field. Describing levels helps other researchers to understand your experiment better';
                uicontrol(f, 'Style', 'text', 'String', msg, 'Units', 'normalized', 'HorizontalAlignment', 'Left','Position', [0.32 0 0.15 0.4],'ForegroundColor', fg,'BackgroundColor', bg, 'Tag', 'levelMsg');
                h = uitable(f, 'Data', t(:,2), 'ColumnName', {'Description'}, 'RowName', t(:,1), 'Units', 'normalized', 'Position', [0.48 0.07 0.5 0.38], 'FontSize', fontSize, 'Tag', 'levelEditTbl', 'CellEditCallback',{@levelEditCB,field},'ColumnEditable',true); 
                h.ColumnWidth = {appWidth*0.5*0.9};
            end
        end
    end

    function levelEditCB(arg1, obj, field)
        level = checkFormat(obj.Source.RowName{obj.Indices(1)});
        description = obj.EditData;
        pInfoBIDS.(field).Levels.(level) = description;
        specified_levels = fieldnames(pInfoBIDS.(field).Levels);
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
            formatted = str;
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
    function pBIDS = newpInfoBIDS()
        pBIDS = [];
%         if isfield(EEG,'BIDS') && isfield(EEG.BIDS,'pInfoDesc') && isfield(EEG.BIDS,'pInfo')
%             for idx=1:size(EEG.BIDS.pInfo,1)
%                 field = EEG.BIDS.pInfo{idx,2}; 
%                 bids_field = EEG.BIDS.pInfo{idx,1};
%                 event.(field).BIDSField = bids_field;
%                 if isfield(EEG.BIDS.pInfoDesc.(bids_field), 'LongName')
%                     event.(field).LongName = EEG.BIDS.pInfoDesc.(bids_field).LongName;
%                 else
%                     event.(field).LongName = '';
%                 end
%                 if isfield(EEG.BIDS.pInfoDesc.(bids_field), 'Description')
%                     event.(field).Description = EEG.BIDS.pInfoDesc.(bids_field).Description;
%                 else
%                     event.(field).Description = '';
%                 end
%                 if isfield(EEG.BIDS.pInfoDesc.(bids_field), 'Units')
%                     event.(field).Units = EEG.BIDS.pInfoDesc.(bids_field).Units;
%                 else
%                     event.(field).Units = '';
%                 end
%                 if isfield(EEG.BIDS.pInfoDesc.(bids_field), 'Levels')
%                     event.(field).Levels = EEG.BIDS.pInfoDesc.(bids_field).Levels;
%                 else
%                     event.(field).Levels = [];
%                 end
%                 if isfield(EEG.BIDS.pInfoDesc.(bids_field), 'TermURL')
%                     event.(field).TermURL = EEG.BIDS.pInfoDesc.(bids_field).TermURL;
%                 else
%                     event.(field).TermURL = '';
%                 end
%             end
%             fields = setdiff(fieldnames(EEG.event), {EEG.BIDS.pInfo{:,2}}); % unset fields
%             for idx=1:length(fields)
%                 event.(fields{idx}).BIDSField = '';
%                 event.(fields{idx}).LongName = '';
%                 event.(fields{idx}).Description = '';
%                 event.(fields{idx}).Units = '';
%                 event.(fields{idx}).Levels = [];
%                 event.(fields{idx}).TermURL = '';
%             end
%         else
            for idx=1:length(pFields)
                if strcmp(pFields{idx}, 'Participant_id')
                    pBIDS.Participant_id.Description = 'Unique subject identifiers';
                    pBIDS.Participant_id.Units = '';
                    pBIDS.Participant_id.Levels = 'n/a';
                elseif strcmp(pFields{idx}, 'Gender')
                    pBIDS.Gender.Description = 'Sex of the subject';      
                    pBIDS.Gender.Levels = [];
                    pBIDS.Gender.Units = '';
                elseif strcmp(pFields{idx}, 'Age')
                    pBIDS.Age.Description = 'Age of the subject';
                    pBIDS.Age.Units = 'years';
                    pBIDS.Age.Levels = 'n/a';
                elseif strcmp(pFields{idx}, 'Group')
                    pBIDS.Group.Description = 'Subject group';
                    pBIDS.Group.Units = '';
                    pBIDS.Group.Levels = [];
%                 else
%                     event.(fields{idx}).BIDSField = '';
%                     event.(fields{idx}).LongName = '';
%                     event.(fields{idx}).Description = '';
%                     event.(fields{idx}).Units = '';
%                     event.(fields{idx}).Levels = [];
%                     event.(fields{idx}).TermURL = '';
                end
            end
%         end
    end
end