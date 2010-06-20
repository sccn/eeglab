Graphic Interface Toolkit
-------------------------
This package contains a set of function to easily create complex graphic
interface (see examples). This package is part of the EEGLAB public 
software for electro-encephalography analysis (which is the leading 
public software for processing electro-encephalography data with more than
70 000 download) and is distributed under the GNU GPL license.

Requirements
------------
Matlab version 7 is required.

Quickstart
----------
Complex graphic interface may be created with arbitrary uicontrols. A simple
example using the inputgui function is shown below

[res userdat strhalt restag] = inputgui('geometry', { 1 1 }, 'uilist', ...
                          { { 'style' 'text' 'string' 'Enter a value' } ...
                            { 'style' 'edit' 'string' '' 'tag' 'editstr'} });

res contains a cell array of output and restag contains a structure with
a field for each output (in the example above res.editstr will contain the
edited string).

Complex geometry for each control may be defined using the 'geom' input 
that works in a way similar to the subplot function (see function help 
message). For example:

res = inputgui('geom', { {2 1 [0 0] [1 1]} {2 1 [1 0] [1 1]} }, 'uilist', ...
                          { { 'style' 'text' 'string' 'Enter a value' } ...
                            { 'style' 'edit' 'string' '' } });

In addition, use the functions questdlg2, errordlg2, warndlg2, listdlg2, 
inputdlg2 to replace standard Matlab functions

ButtonName = questdlg2('What is your favorite color?', 'Color Question', ...
                         'Red', 'Green', 'Blue', 'Green');

GUI color may be customized for these graphic interfaces by defining a
file containing colors (see icadefs.m file in the EEGLAB software 
distribution).

For more information see the help message of each m function.

Content
-------
inputgui      - main function to create graphic interface
supergui      - support function creating the actual figure
finputcheck   - support function for checking function input validity
pophelp       - shows a function code in the Matlab browser
warmdlg2      - uses inputgui to replace the Matlab warndlg function
errordlg2     - uses inputgui to replace the Matlab errordlg function
listdlg2      - uses inputgui to replace the Matlab listdlg function
inputdlg2     - uses inputgui to replace the Matlab inputdlg function
questdlg2     - uses inputgui to replace the Matlab questdlg function

The latest revision of these function may be checked out using SVN at

http://sccn.ucsd.edu/repos/software/eeglab/functions/guifunc/

Arnaud Delorme, PhD
June 19th 2010
University of San Diego California, USA
CNRS, University of Toulouse Paul Sabatier, France
