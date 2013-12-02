function installRes = plugin_askinstall(pluginName, pluginFunc, forceInstall)

if nargin < 3, forceInstall = false; end;

if nargin < 2 || ~exist(pluginFunc)
    
    if ~forceInstall
        db = dbstack;
        if length(db) > 2 && ~strcmpi(db(end).name, 'checkouteeglab.m');
            error([ 'Cannot find ' pluginName ' extension, use EEGLAB Extension Manager to install it' ]);
        end;

        installRes = 0;

        % check is deactivated
        try, PLUGINLIST = evalin('base', 'PLUGINLIST'); catch, PLUGINLIST = []; end;
        if ~isempty(PLUGINLIST) && isfield(PLUGINLIST, 'plugin')
            indPlugin = strmatch(lower(pluginName), lower({ PLUGINLIST.plugin }), 'exact');
            if ~isempty(indPlugin) && strcmpi(PLUGINLIST(indPlugin(1)).status, 'deactivated')
                res = questdlg2( [ pluginName ' extension is de-activated. Do you want to reactivate it now?' ], [ pluginName ' extension installation' ], 'No', 'Yes', 'Yes' );
                if strcmpi(res, 'no'), return, end;
                plugin_reactivate(PLUGINLIST(indPlugin(1)).foldername);
                evalin('base', 'eeglab rebuild');
                installRes = 1;
                return;
            end;
        end;

        % check for installing
        res = questdlg2( [ pluginName ' extension is not installed. Do you want to download it now?' ], [ pluginName ' extension installation' ], 'No', 'Yes', 'Yes' );
    else
        res = 'yes';
    end;
    
    if strcmpi(res, 'no'), return, end;
    plugins = plugin_getweb('import', []);
    indPlugin = strmatch(lower(pluginName), lower({ plugins.name }), 'exact');
    if isempty(indPlugin), 
        plugins = plugin_getweb('process', []);
        indPlugin = strmatch(lower(pluginName), lower({ plugins.name }));
        if isempty(indPlugin), 
            error([ pluginName ' extension not found' ]); 
        end;
    end;
    result = plugin_install(plugins(indPlugin(1)).zip, plugins(indPlugin(1)).name, plugins(indPlugin(1)).version, forceInstall);
    if result == 1, installRes = 1; end;
    
    evalin('base', 'eeglab rebuild');
else
    installRes = 1;
end;
