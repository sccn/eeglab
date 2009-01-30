function plot_coherence(data,locs,pars);
% usage plot_coherence(data,locs,pars);
% makes head-in-head plots
% data is an NxN matrix where N is the number of channels
% locs is (ideally) an Nx5 matrix:
%      1st column: a channel i gets a circle only if locs(i,1)>.5 
%                  if locs is an Mx5 matrix with M<N
%                  the absolute values of the first column 
%                  are interpreted as indices of scalp-electrodes.
%                  If the values themselves are smaller than 0 then
%                  these electrodes are only included within each 
%                  small circle but do not correspond to a small circle itself.
%      2nd and 3rd column: x,y coordinates of electrodes in 2D
%      4th and 5th column: x,y coordinates of centers of spheres
%                          slight deformation of electrode locations 
%                          to avoid overlapping spheres
%      if locs is Nx2, the 2 columns are interpreted as x,y coordinates 
%      for both electrodes and circle-centers. 
%      if locs is Nx3, the 2nd and 3rd column are interpreted as x,y coordinates 
%      for both electrodes and sphere-centers. The first as interpreted as for 
%      5xN case.
% pars sets parameters; 
%      pars.scale sets color-scale; it is a 1x2 vector denoting min 
%                 and max of colorbar. defaults is  pars.scale=[-max(max(abs(data))),max(max(abs(data)))];
%      pars.resolution sets resolution. Default is pars.resolution=25
%                      Increasing it makes better pictures but is slower
%      pars.global_size  sets global size factor of small circles and distancec between circles.                   
%                        default: global_size=1
%      pars.relative_size  sets size factor for distances leaving circle-size unchanged                   
%                        default: relative_size=1
%      pars.head_up    moves the big circle (the head) up. Default: head_up=0
%      pars.head_right    moves the big circle (the head) to the right. Default: head_right=0
%

if nargin>2
    if isfield(pars,'scale')
        scal=pars.scale;
        if length(scal)==1
            scal=[-scal scal];
        end
    else
        scal=[-max(max(abs(data))),max(max(abs(data)))];
    end

    if isfield(pars,'global_size');
        global_size=pars.global_size;
    else
        global_size=1;
    end
    if isfield(pars,'relative_size');
        rel_size=pars.relative_size;
    else
        rel_size=1;
    end
    if isfield(pars,'resolution')
        resolution=pars.resolution;
    else
        resolution=25;
    end

    if isfield(pars,'head_right')
        head_right=pars.head_right;
    else
        head_right=0;
    end
    if isfield(pars,'head_up')
        head_up=pars.head_up;
    else
        head_up=0;
    end

else
     scal=[-max(max(abs(data))),max(max(abs(data)))];
    resolution=25;
    global_size=1;
    rel_size=1;
    head_right=0;
    head_up=0;

end



[n,m]=size(locs);
locs_all=locs;
[nall,ndum]=size(locs_all);
ind2chan=(1:nall)';
[nc,nc]=size(data);

if m==5
  if nc>n
    ind2chan=abs(locs(:,1));
  end
  indd=ind2chan(locs(:,1)>.5);
  locs=locs(locs(:,1)>.5,:);
end


[n,m]=size(locs);
if m==2;
    locs=[(1:n)',locs,locs];
elseif m==3
    locs=[locs,locs(:,2:3)];
elseif m==4
    locs=[(1:n)',locs];
end
%ind2chan=locs(:,1);

locs(:,5)=detrend(locs(:,5),'constant');
%locs(:,4)=detrend(locs(:,4),'constant');


coor=locs(:,4:5);

ymin=min(coor(:,1));
ymax=max(coor(:,1));
xmin=max(coor(:,1));
xmax=max(coor(:,1));


%minmin=mindis(coor(:,1),coor(:,2));
[minmin,minmax,meanmin]=mindis(coor(:,1),coor(:,2));

tot_scale=max(abs([xmin,xmax,ymin,ymax]))+minmin/2;
tot_scale=tot_scale/1.1;
locs(:,2:5)=locs(:,2:5)*.5/tot_scale;
meanmin=meanmin*.5/tot_scale;


no=zeros(n,1);
for i=1:n
    no(i)=norm(coor(i,:));
end
faktor=1/max(no)/2/1.3;
%wfaktor=faktor*1.3
wfaktor=meanmin/1.6;

wfaktor=wfaktor*global_size;


meanx=mean(faktor*locs(:,4)*.85+.45);
meany=mean(faktor*locs(:,5)*.85+.45);

rad=max(sqrt( (faktor*locs(:,4)*.85+.45-meanx).^2+(faktor*locs(:,5)*.85+.45-meany).^2 ));
cirx=(rad*cos((1:1000)*2*pi/1000)+meanx)';ciry=(rad*sin((1:1000)*2*pi/1000)+meany)';
   faktor=.75;
faktor=faktor*global_size*rel_size;

nplot=n;
for chan=1:n
    zp=data(indd(chan),ind2chan)';
    %zp=z(ind2chan);
  
  
  
  
  subplot(nplot,nplot,n*chan);
%  set(gca,'Position',[faktor*locs(chan,4)*.85+.45 faktor*(locs(chan,5)*1.-.05)+.45 wfaktor*0.070 wfaktor*0.088]);
    set(gca,'Position',[faktor*locs(chan,4)*.85+.45 faktor*(locs(chan,5)*1.+.01)+.45 wfaktor wfaktor*.088/.07]);

  plot_elec_empty_lowres(zp,locs_all(:,2:3),scal,indd(chan),resolution);
  set(gca,'visible','off');
  set(gca,'xTick',[])
  set(gca,'yTick',[])
  caxis('manual');
  caxis(scal);
  if chan==1
      %set(gca,'Position',[faktor*locs(chan,4)*.85+.45 faktor*(locs(chan,5)*1.-.05)+.45 wfaktor*0.08 wfaktor*0.08]);
       set(gca,'Position',[faktor*locs(chan,4)*.85+.45 faktor*(locs(chan,5)*1.+.01)+.45 wfaktor*1.1 wfaktor*.088/.07]);

      h1=colorbar;
      set(h1,'position',[.85 .1 .05 .8]);
  end
end;
caxis('manual');
caxis(scal);
subplot(nplot,nplot,1);
set(gca,'Position',[0.12 0.01 .7/.88 .99]);
%set(gca,'Position',[0.12 0.01 1 1]);

scalx=1.0;




drawhead(meanx+wfaktor/4+head_right,.45+wfaktor/4+head_up,.43,scalx);
axis([ 0 1 0 1]);
%text(.1,.85,figtitle,'HorizontalAlignment','center');
%title('Titel')
set(gca,'visible','off');
%get(h1)
%h=colorbar; 
%set(h,'yticklabel',

return;

function plot_elec_empty_lowres(z,loc,skal,chan,resolution); 

x=loc(:,1);
y=loc(:,2);

xlin = linspace(1.4*min(x),1.4*max(x),resolution);
ylin = linspace(1.4*min(y),1.4*max(y),resolution);
[X,Y] = meshgrid(xlin,ylin);
Z = griddata(x,y,z,X,Y,'invdist');
%Z = griddata(x,y,z,X,Y,'nearest');


  % Take data within head
  rmax=1.1*max(sqrt(x.^2+y.^2));
  mask = (sqrt(X.^2+Y.^2) <= rmax);
  ii = find(mask == 0);
  Z(ii) = NaN;
  
  
surface(X,Y,zeros(size(Z)),Z,'edgecolor','none');shading interp;
%caxis([ - max(abs(z)) max(abs(z))]);
%disp([ - max(abs(z)) max(abs(z))]);
%caxis([ -skal  skal]);
%disp([ -skal  skal]);

hold on;
%plot(x,y,'.k');
axis([-rmax rmax -rmax rmax]);
%colorbar;
plot(loc(chan,1),loc(chan,2),'.k');
%set(gcf,'color','none');

set(gca,'xTick',[])
set(gca,'yTick',[])
plot(.985*rmax*sin((0:1000)/1000*2*pi), .985*rmax*cos((0:1000)/1000*2*pi),'linewidth',2,'color','k'); 
%set(gcf,'color','none');

return; 

function drawhead(x,y,size,scalx);

cirx=(x+scalx*size*cos((1:1000)*2*pi/1000) )';ciry=(y+size*sin((1:1000)*2*pi/1000))';

plot(cirx,ciry,'k','linewidth',1);
hold on;

ndiff=20;
plot( [x  cirx(250-ndiff) ],[y+1.1*size ciry(250-ndiff)],'k','linewidth',1);
plot( [x  cirx(250+ndiff) ],[y+1.1*size ciry(250+ndiff)],'k','linewidth',1);


return;

function [minmin,minmax,meanmin]=mindis(x,y);

[n,m]=size(x);

minall=zeros(n,1);



for i=1:n
    md=1.e8;
    for j=1:n
        if j~=i
          %dis=max( abs(x(i)-x(j)),abs(y(i)-y(j)));
          dis=norm([x(i)-x(j),y(i)-y(j)]);
          if dis<md
             md=dis;
          end
        end     
    end
    minall(i)=md;
end



minmin=min(minall);
meanmin=mean(minall);
[minmax,imax]=max(minall);

return;

