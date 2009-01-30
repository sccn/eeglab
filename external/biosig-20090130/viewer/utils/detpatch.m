function result=detpatch(varargin)
% DETPATCH
% Select HELP in the Info menu 
%
% Version 1.0, November 2004
% Copyright by (C) Franz Einspieler <znarfi5@hotmail.com> and
%                  Alois Schloegl   <a.schloegl@ieee.org>
% University of Technology Graz, Austria
%
% This is part of the BIOSIG-toolbox http://biosig.sf.net/
% Comments or suggestions may be sent to the author.
% This Software is subject to the GNU public license.

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mousebutton=get(gcf,'SelectionType');
if isequal(mousebutton,'alt')
    if isequal(varargin{2},'move_line') | isequal(varargin{2},'scale') | isequal(varargin{2},'drawend')
        if isequal(mousebutton,'normal')
            varargin{2}='new_patch';
        end
    else
        try
            Data=get(findobj('Tag','sviewer'),'UserData');
            s = Data.Eventcodes_txt;
            akt_dettext_posEM = find(Data.Detection.EventMatrix(:,5) == str2num(varargin{6}));
            akt_dettext_dettype = Data.Detection.EventMatrix(akt_dettext_posEM,2);
            akt_dettext_posCI = find(s.CodeIndex == akt_dettext_dettype);
            akt_dettext = s.CodeDesc{akt_dettext_posCI};
            h=msgbox(['Selected Event: ' akt_dettext],'Selected Event:','help');
            return
        end
    end
end
try
	if isequal(mousebutton,'extend')
        Data = get(findobj('Tag','sviewer'),'UserData');
        tag = varargin{4};
        id = str2num(varargin{6});
        patch_pos = findobj(gcf,'Tag',tag); %gca
        leftline_pos= findobj(gcf,'Tag',[tag,'_left']); %gca
        rightline_pos= findobj(gcf,'Tag',[tag,'_right']);   %gca                
        slider_step = get(findobj('Tag','Slider1'),'SliderStep');
        pos = Data.Slider.Pos;
        X = pos/slider_step(1)*Data.HDR.SampleRate*Data.NoS / 2;
        del_line = find((Data.Detection.EventMatrix(:,5) == id));
        if Data.Detection.EventMatrix(del_line,3) == 0
            Data.Detection.EventMatrix(del_line,:)=[];
            set(findobj('Tag','sviewer'),'UserData',Data);
            delete(patch_pos);
            delete(leftline_pos);
            delete(rightline_pos);
            return;
        end
        Data.Detection.EventMatrix(del_line,:)=[];
        set(findobj('Tag','sviewer'),'UserData',Data); 
        delete(patch_pos);
        delete(leftline_pos);
        delete(rightline_pos); 
        return;
	end
end

if isstr(varargin{1})
  action = 'start';
  for i=1:2:length(varargin)
      switch lower(varargin{i})
          case 'action'
              action = varargin{i+1};
          case 'position'
              patchval.patch_pos = varargin{i+1};
          case 'facecolor'
              patchval.facecolor = varargin{i+1};
          case 'tag'
              tag = varargin{i+1};
          case 'id'
              id = num2str(varargin{i+1});
          case 'data'
              Data = varargin{i+1};
          otherwise
              disp(varargin{i});
      end
  end  
end 


switch action,
    % draw patch
    case 'start',
        x1 = patchval.patch_pos(1);
        x2 = patchval.patch_pos(1) + patchval.patch_pos(3);
        y1 = patchval.patch_pos(2);
        y2 = patchval.patch_pos(2) + patchval.patch_pos(4);
        
        try subplot(Data.SPlot(Data.Patch.actchannelpos));end
        patch_h = patch([x1 x2 x2 x1],[y1 y1 y2 y2],'k');
        set(patch_h,'facecolor',patchval.facecolor,...
                'erasemode','xor',...
                'linestyle','none',...
                'buttondownfcn',[mfilename,'(''action'',''new_patch'',''tag'',''',tag,''',''id'',''',id,''')'],...
                'Tag',tag);
        % draw line
        x = get(patch_h,'XData');
        y = get(patch_h,'YData');
        line_lr(1) = line([x(1) x(1)],[y(1) y(4)],'Tag',[tag,'_left']);
        line_lr(2) = line([x(2) x(2)],[y(2) y(3)],'Tag',[tag,'_right']);
        set(line_lr,'ButtonDownFcn',[mfilename,'(''action'',''move_line'',''tag'',''',tag,''',''id'',''',id,''')'],...
                'color',[.3 .3 .3], ...
                'linewidth',1.5 );
        if nargout > 0,
            result = patch_h;
        end
    % draw new line
    case 'new_patch',
        sviewer('Drawline_Callback',gcbo,[],guidata(gcbo))
    % move line
    case 'move_line',
        tag = varargin{4};
        leftline_pos= findobj(gcf, 'Tag', [tag,'_left'] );
        rightline_pos= findobj(gcf, 'Tag', [tag,'_right'] );
        startpos = get(leftline_pos,'XData');                 
        endpos = get(rightline_pos,'XData');                  
        
        if strcmp(get(gcbo,'Tag'),[tag,'_left']),
            set(gcf,'pointer','left');
        else
            set(gcf,'pointer','right');
        end
        set(gcf,'WindowButtonMotionFcn',...
            [mfilename,'(''action'',''scale'',''tag'',''',tag,''',''id'',''',id,''')']);
        set(gcf,'WindowButtonUpFcn', ...
            [mfilename,'(''action'',''drawend'',''tag'',''',tag,''',''id'',''',id,''')']);
    % activate change
    case 'drawend'
        Data=get(findobj('Tag','sviewer'),'UserData');
        slider_step = get(findobj('Tag','Slider1'),'SliderStep');
        pos = Data.Slider.Pos;
        X = pos/slider_step(1)*Data.HDR.SampleRate*Data.NoS / 2;
        sel_channel = get(gca, 'UserData');
        id = str2num(varargin{6});
        del_line = find((Data.Detection.EventMatrix(:,5) == id));
        
        if isequal(get(gcbo,'Pointer'),'right')
            rightline_pos= findobj(gca,'Tag', [tag,'_right'] );
            endpos = get(rightline_pos,'XData');
            set(gcf,'pointer','arrow',...
                'WindowButtonMotionFcn','',...
                'WindowButtonUpFcn','');
            startpos = Data.Detection.EventMatrix(del_line,1);
            Data.Detection.EventMatrix(del_line,4)=[round(X + endpos(1) - Data.Detection.EventMatrix(del_line,1))];
            Data.Detection.Update = Data.Detection.EventMatrix(del_line,:);
            set(findobj('Tag','sviewer'),'UserData',Data); 
            return
        end
        
        if isequal(get(gcbo,'Pointer'),'left')
            leftline_pos= findobj(gca,'Tag', [tag,'_left'] );
            startpos = get(leftline_pos,'XData'); 
            if startpos == 0
                startpos = 1;
            end
            set(gcf,'pointer','arrow',...
                'WindowButtonMotionFcn','',...
                'WindowButtonUpFcn','');
            startpos_old = Data.Detection.EventMatrix(del_line,1);
            endpos_old = Data.Detection.EventMatrix(del_line,4);
            startpos_new = X + startpos(1);
            Data.Detection.EventMatrix(del_line,1) = [round(startpos_new)];
            Data.Detection.EventMatrix(del_line,4) = [round(endpos_old + (startpos_old-startpos_new))];
            Data.Detection.Update = Data.Detection.EventMatrix(del_line,:);
            set(findobj('Tag','sviewer'),'UserData',Data); 
        end
    % scale
    case 'scale',
        patch_pos = findobj(gcf, 'Tag', tag );
        len_patch_pos = size(patch_pos,1);
        if len_patch_pos > 1
            patch_pos = findobj(gcf, 'Tag', tag );
            direction = get(gcf,'pointer');
            lh = findobj(gcf, 'Tag', [tag,'_',direction] );
            curr_point = get(gca,'currentpoint');
            pos_p = get(patch_pos(1),'XData');
            pos_min = max([min(xlim) -Inf]); 
            pos_max = min([max(xlim) Inf]); 
            scaling = 0;
            posnew_p = [1]*(fix(curr_point(1)/[1])+round(rem(curr_point(1),[1])/[1]));
            for i=1:len_patch_pos
                pos_p1 = get(patch_pos(i),'XData');
                switch direction
                    case 'left',
                        if ( posnew_p >= pos_min & posnew_p <= pos_p1(2) & posnew_p ~= pos_p1(1) )
                            pos_p1([1 4]) = [posnew_p posnew_p];
                            set(lh(i),'XData',[posnew_p posnew_p]);
                            set(patch_pos(i),'XData',pos_p1);
                            scaling = 1;
                        end
                    case 'right',
                        if ( posnew_p <= pos_max & posnew_p >= pos_p1(1) & posnew_p ~= pos_p1(2) )
                            pos_p1([2 3]) = [posnew_p posnew_p];
                            set(lh(i),'XData',[posnew_p posnew_p]);
                            set(patch_pos(i),'XData',pos_p1);
                            scaling = 1;
                        end
                    otherwise
                        disp([mfilename,' : Unknown direction ''',direction,''''])
                        return;
                end
            end
        else
            patch_pos = findobj(gca, 'Tag', tag );
            direction = get(gcf,'pointer');
            lh= findobj(gca, 'Tag', [tag,'_',direction] );
            curr_point = get(gca,'currentpoint');
            pos_p = get(patch_pos,'XData');
            pos_min = max([min(xlim) -Inf]); 
            pos_max = min([max(xlim) Inf]); 
            scaling = 0;
            posnew_p = [1]*(fix(curr_point(1)/[1])+round(rem(curr_point(1),[1])/[1]));
            switch direction,
                case 'left',
                    if ( posnew_p >= pos_min & posnew_p <= pos_p(2) & posnew_p ~= pos_p(1) ),
                        pos_p([1 4]) = [posnew_p posnew_p];
                        set(lh,'XData',[posnew_p posnew_p]);
                        set(patch_pos,'XData',pos_p);
                        scaling = 1;
                    end
                case 'right',
                    if ( posnew_p <= pos_max & posnew_p >= pos_p(1) & posnew_p ~= pos_p(2) ),
                        pos_p([2 3]) = [posnew_p posnew_p];
                        set(lh,'XData',[posnew_p posnew_p]);
                        set(patch_pos,'XData',pos_p);
                        scaling = 1;
                    end
                otherwise
                    disp([mfilename,' : Unknown direction ''',direction,''''])
                    return
            end
        end
    otherwise
        disp(['Unknown action : ',action]);  
end 