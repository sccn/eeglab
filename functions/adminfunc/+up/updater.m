classdef updater < handle
    % Web-based software updater class with download and installation capabilities.
    % Copyright (C) 2012 by Nima Bigdely-Shamlo, Swartz Center for Computational Neuroscience, INC,
    % UCSD.
    % Usage is only allowed under the same (non-commercial) license as Measure Projection toolbox.
    
    properties
        currentVersionNumber
        currentVersionReleaseDate   % a string containing time and date of the current release
        latestVersionNumber
        releaseDate
        downloadUrl
        releaseNotes      % text about the release.
        releaseNotesUrl   % The URL in which more information about the release is provided.
        newMajorRevision
        xmlFileUrl
        lastTimeChecked % a date vector returned by clock().
        menuItemHandle = nan; % handle to the menu item that becomes visible if a new version is available.
        menuItemText = 'Upgrade to the Latest Version'; % text label for the item when a newer version is available.
        menuItemCallback = [];
        downloadedFileName % contains file name (and path) of the downloaded latest-version file.
        guiHandle
        softwareName       % the name of the software package to be upadted (e.g. EEGLAB)
    end
    
    methods (Static = true)

        function writeReleaseXml(xmlFileName, versionNumber, downloadUrl, releaseNotes, releaseNotesUrl, releaseDate)
            % writeReleaseXml(xmlFileName, versionNumber, downloadUrl, releaseNotes, releaseNotesUrl, releaseDate)
            % use this function to make the XML file to  be placed somewhere which is accessible from
            % the web.
            
            if nargin < 4
                releaseNotes = '';
            end;
            
            if nargin < 5
                releaseNotesUrl = '';
            end;
            
            if nargin < 6
                releaseDate = datestr(now); % e.g. '0-Jan-2012 11:54:10'
            end;
            
            % make sure version number is a string  as only strings can be set as attributes.
            if isnumeric(versionNumber)
                versionNumber = num2str(versionNumber);
            end;
            
            docNode = com.mathworks.xml.XMLUtils.createDocument('latestVersion');
            docRootNode = docNode.getDocumentElement;
            docRootNode.setAttribute('versionNumber', versionNumber);
            docRootNode.setAttribute('releaseDate', releaseDate);
            docRootNode.setAttribute('downloadUrl', downloadUrl);
            docRootNode.setAttribute('releaseNotes', releaseNotes);
            docRootNode.setAttribute('releaseNotesUrl', releaseNotesUrl);
            docRootNode.setAttribute('newMajorRevision', releaseNotesUrl);
            
            xmlwrite(xmlFileName, docNode);
        end
        
        function homeDir = getHomeDirectory
            if ispc
                homeDir= getenv('USERPROFILE');
            else
                homeDir = getenv('HOME');
            end
        end;
    end
    
    methods
        function obj = updater(currentVersionNumber, xmlFileUrl, softwareName, currentVersionReleaseDate)
            % obj = updater(currentVersionNumber, xmlFileUrl, softwareName, currentVersionReleaseDate)
                        
            obj.currentVersionNumber = currentVersionNumber;
            obj.xmlFileUrl = xmlFileUrl;
            
            if nargin > 2
                obj.softwareName = softwareName;
            end;
            
            if nargin > 3
                obj.currentVersionReleaseDate = currentVersionReleaseDate;
            end;
        end
        
        function obj = checkForNewVersion(obj, nameValuePairs)
            % nameValuePairs is a cell array of name1, value1, name 2, value 2,... send along with
            % the URL request.
            
            if nargin < 2
                nameValuePairs = {};
            end;
            
            % add current version number and some other basic setup info to request URL
            if ispc
                osType = 'pc';
            elseif isunix
                osType = 'unix';
            elseif ismac
                osType = 'max';
            else
                osType = 'unknown';
            end;
            
            nameValuePairs = [{'currentVersionNumber' num2str(obj.currentVersionNumber) 'matlabVersion' version 'OS' osType} nameValuePairs];
            
            try
                [xmlString successfullRead] = plugin_urlread([obj.xmlFileUrl], 'get', nameValuePairs);
            catch
                successfullRead = false;
            end;
            
            if successfullRead
                
                % write into a temp file
                temporaryFileName = tempname;
                fileId = fopen(temporaryFileName, 'w');
                fprintf(fileId, '%s', xmlString);
                fclose(fileId);
                
                % read as xml into a doc
                try
                    docNode = xmlread(temporaryFileName);
                    docRootNode = docNode.getDocumentElement;
                catch
                    successfullRead = false;
                    return;
                end;
                
                % delete temporary file afterwards
                delete(temporaryFileName);
                
                obj.latestVersionNumber = str2double(char(docRootNode.getAttribute('versionNumber')));
                obj.releaseDate = char(docRootNode.getAttribute('releaseDate'));
                obj.releaseNotes = char(docRootNode.getAttribute('releaseNotes'));
                obj.releaseNotesUrl = char(docRootNode.getAttribute('releaseNotesUrl'));
                obj.downloadUrl = char(docRootNode.getAttribute('downloadUrl'));
                obj.newMajorRevision = char(docRootNode.getAttribute('newMajorRevision'));
                
                obj.lastTimeChecked = clock;
                
                % if a menu item handle is assigned, make it visible and set its label
                % appropriately (show the menu if the version is up-to-date and hide it otherwise).
                if ishandle(obj.menuItemHandle)
                    if obj.latestVersionNumber > obj.currentVersionNumber
                        set(obj.menuItemHandle, 'label', [obj.menuItemText ' (' num2str(obj.latestVersionNumber) ')']);
                        set(obj.menuItemHandle, 'Callback', obj.menuItemCallback);
                        set(obj.menuItemHandle, 'visible', 'on');
                    else
                        set(obj.menuItemHandle, 'visible', 'off');
                    end;
                end;
            end;
        end;
        
        function answer = newerVersionIsAvailable(obj)
            answer = obj.latestVersionNumber > obj.currentVersionNumber;
        end;
        
        function timeAGo = howLongAgoChecked(obj) % how ago was the last version check, in seconds
            if isempty(obj.lastTimeChecked)
                timeAGo = Inf;
            else
                timeAGo = etime (clock, obj.lastTimeChecked);
            end;
        end;
        
        function obj = downloadLatestVersion(obj, folderToPlace)
            
            if usejava('jvm')
                obj.downloadedFileName = fullfile(tempdir, ['latest_version' num2str(obj.latestVersionNumber) '_' strrep(char(java.util.UUID.randomUUID),'-','_')]);
            else
                obj.downloadedFileName = fullfile(tempdir ,['latest_version' num2str(obj.latestVersionNumber) '_' num2str(feature('timing','cpucount'))]);
            end
            
            % make it a bit shoter
            obj.downloadedFileName = obj.downloadedFileName(1:end-30);
            
            % get the file part fromthe download URL and add it to the end of temp file.
            % One advantage of this is to preserve file type (e.g. zip)
            segmentFromUrl = obj.downloadUrl(find(obj.downloadUrl == '/', 1, 'last')+1:end);
            obj.downloadedFileName = [obj.downloadedFileName '_' segmentFromUrl];
            
            fprintf('Downloading the latest version (%s)....\n', num2str(obj.latestVersionNumber));
            [filestr downloadWasSucessfull] = urlwrite(obj.downloadUrl, obj.downloadedFileName);
            
            if downloadWasSucessfull
                fprintf(['Latest-version file successfully download and is located at: ' obj.downloadedFileName '\n']);
            else
                fprintf('Latest-version file could not be downloaded.\n');
            end;
        end;
        
        function obj = installDownloadedFile(obj, installFolder, goOneFolderLevelIn)
            
            if nargin < 2
                error(' You must specy a folder to install the latest version');
            end;
            
            if nargin < 3
                goOneFolderLevelIn = false;
            end;
            
            [pathstr, name, ext] = fileparts(obj.downloadedFileName);
            
            temporaryFolder = tempname;
            fprintf('Unpacking downloaded file...\n');
            if strcmpi(ext, '.zip')
                mkdir(temporaryFolder);
                unzip(obj.downloadedFileName, temporaryFolder);
            elseif strcmpi(ext, '.tar')
                mkdir(temporaryFolder);
                untar(obj.downloadedFileName, temporaryFolder);
            else
                fprintf('File type in unknown and cannot be unzipped.\n');
                return;
            end;
            
            % remove all files from installation folder
            if exist(installFolder)
                fprintf('Removing all files from installation folder...\n');
                rmdir(installFolder,'s');
            end;
            
            if goOneFolderLevelIn
                listing = dir(temporaryFolder);
                insideFolderName = '';
                for i=1:length(listing)
                    if ~(strcmp(listing(i).name, '.') || strcmp(listing(i).name, '..'))
                        insideFolderName = listing(i).name;
                        break;
                    end;
                end;
                
                temporaryFolder = [temporaryFolder filesep insideFolderName];
            end;
            
            % copy extracted files and folder into the install (destination) folder
            fprintf('Copying files from the latest version into the installation folder...\n');
            movefile(temporaryFolder, installFolder, 'f');
            
            fprintf('New version (%s) is now installed. You may need to clear functions and objects or exit Matlab in order to use it.\n', num2str(obj.latestVersionNumber));
        end;
        
        function numberOfDays = daysCurrentVersionIsOlder(obj)
            % numberOfDays = daysCurrentVersionIsOlder(obj)
            % return the number of days between current version release and the latest version
            % release. (empty is current release date in unknown).
            
            if isempty(obj.currentVersionReleaseDate)
                numberOfDays = []; % cannot be calculated
            else
                numberOfDays = round(etime(clock, datevec(datenum(obj.currentVersionReleaseDate))) / (3600*24) );
            end;
        end;
        
        function downloadButtonCallBack(handle, event, obj)
            % close the figure
            close(get(handle, 'parent'));
            drawnow;
            
            obj.downloadLatestVersion;
        end;
                
        function installButtonCallBack(handle, event, obj, installDirectory, goOneFolderLevelIn)
            
            % get the tag ifnromation from figure GUI.
            guiData = get(get(handle, 'parent'), 'userdata');
            
            % close the figure
            close(get(handle, 'parent'));
            drawnow;
            
            % download the file if it has not been downloaded yet.
            if isempty(obj.downloadedFileName)
                 obj.downloadLatestVersion;
            end;
            
            obj.installDownloadedFile(installDirectory, goOneFolderLevelIn);
            
            % evaulate post-install callback (text) after installation.
            if ~isempty(guiData.postInstallCallbackString)
                evalin('base', guiData.postInstallCallbackString);
            end;
        end;
        
        function launchGui(obj, installDirectory, goOneFolderLevelIn, backGroundColor, postInstallCallbackString)
            % launchGui(obj, installDirectory, goOneFolderLevelIn, backGroundColor, postInstallCallbackString)
            % installDirectory:     The directory into which the new version will nbe installed. 
            %                       All files in this directory will be deleted before installation.
            %
            % goOneFolderLevelIn:   If set to true, the folder one level inside the unpakced
            %                       (e.g. unzipped) file will be copied into the install directory.
            %                       This is useful if the downloaded file contains one folder which
            %                       is the actual software to be installed.
            % backGroundColor       GUI background color. {optional}
            %
            % postInstallCallbackString commands (callback as a string) to be executed after
            %                                    installation from the GUI. {optional}
            
            if nargin < 3
                goOneFolderLevelIn = false;
            end;
            
            if nargin < 5
                postInstallCallbackString = [];
            end;
            
            updaterFolder = fileparts(which('up.updater'));
            obj.guiHandle = open([updaterFolder filesep 'updater_gui.fig']);
            
            % save some information in tag property of the figure GUI.
            guiData = struct;
            guiData.postInstallCallbackString = postInstallCallbackString;
            set(obj.guiHandle, 'userData', guiData);
            
            if nargin > 3
                warningTextHandle = findobj(obj.guiHandle, 'tag', 'warningText');
                set(obj.guiHandle, 'color', backGroundColor);                
                set(warningTextHandle, 'backgroundColor', backGroundColor);
            end;
            
            % full the texbox
            textBoxHandle = findobj(obj.guiHandle, 'tag', 'textBox');
            
            if isempty(obj.softwareName)
                textBoxTextCell = {['A newer version is available and ready for download.']};
            else
                textBoxTextCell = {['A newer version of ' obj.softwareName ' is available and ready for download.']};
            end;
            
            currentVersionString = ['You are currently using version (' num2str(obj.currentVersionNumber) ')'];
            if isempty(obj.daysCurrentVersionIsOlder)
                currentVersionString = [currentVersionString '.'];
            else
                currentVersionString = [currentVersionString ' which is ' num2str(obj.daysCurrentVersionIsOlder) ' days older.'];
            end;
            
            textBoxTextCell  = cat(1, textBoxTextCell, currentVersionString);
            
            % release notes
            if ~isempty(obj.releaseNotes) && ~strcmp(obj.releaseNotes, '');
                textBoxTextCell  = cat(1, textBoxTextCell, '------------------------------');
                textBoxTextCell  = cat(1, textBoxTextCell, 'Release Notes:');
                textBoxTextCell  = cat(1, textBoxTextCell, obj.releaseNotes);
            end;
            
            set(textBoxHandle, 'string', textBoxTextCell);
            
            % set callbakcs for Install and Download buttonds
            downloadButtonHandle =  findobj(obj.guiHandle, 'tag', 'DownloadButton');
            set(downloadButtonHandle, 'callback', {@downloadButtonCallBack, obj});
            
            installButtonHandle =  findobj(obj.guiHandle, 'tag', 'installButton');
            if nargin<3
                % hide install button if install directory is not provided.
                set(installButtonHandle, 'visible', 'off');
            else
                set(installButtonHandle, 'callback', {@installButtonCallBack, obj, installDirectory, goOneFolderLevelIn});
            end;
            
        end;
    end
end