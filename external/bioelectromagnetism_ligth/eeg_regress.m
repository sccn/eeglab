function [H] = eeg_regress(data,robust,GP,labels)

% eeg_regress - plot multiple linear regression results
% 
% [H] = eeg_regress(data,[robust],[GP],labels)
% 
% data is NxP matrix, with N observations, P-1 predictors
% and the last column is the observed values to predict.
% So, if only 2 columns are given it does a simple linear
% regression.
% 
% robust = 1, use robust least-squares regression
% robust = 0, use least-squares regression
% note, current robust is forced to 1.
% 
% GP is Nx1 grouping variable
%
% labels = cell strings with variable names, keep to
% 5 characters or less for best format of figure text
% 
% H is an array of handles to figures
% 
% All the results are plotted to new figures, which can
% be saved/exported using print(H(i)...)
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:53 $

% Licence:  GNU GPL, no express or implied warranties
% History:  11/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('data','var'),
  fprintf('...no input data.\n\n');
  return
end

if ~exist('robust','var'),
  robust = 1;
else,
  warning('must use robust method, at present');
  robust = 1; % cheat for now, as some plotting stuff below requires it
end

if ~exist('labels','var'),
  for i = 1:size(data,2),
    if i == size(data,2),
      labels{i} = 'Y';
    else
      labels{i} = sprintf('X%d',i);
    end
  end
end

if ~exist('GP','var'),
  error('GP input variable is undefined\n');
end


plotmatrix(data); % to explore relationships
H(1) = gcf;


y = data(:,end);           % values to fit from last column
e = ones(length(data),1);  % for constant term in model
X = [e data(:,1:end-1)];   % predictor values from first columns


if ~robust,
  %beta = X\y;              % 3 x 1 matrix, const, beta1 & beta2
  [B,BCI,R,RCI,STATS] = regress(y,X,.05);
  P = X*B;      % fitted values
  
  % use B slope, with BCI intercepts, not strictly correct
  CI = [ [BCI(1,1); B(2:end)]  [BCI(1,2); B(2:end)] ];
  YCI = [X*CI(:,1),  X*CI(:,2) ];
  
  % Note if BCI(x,:) < or > 0, B(x) is significant
  
  % define figure rows & columns
  cols = 6;
  rows = size(X,2);
  
else
  [B,BCI,R,RCI,STATS] = regress(y,X,.05);
  clear B BCI R RCI;
  
  R2 = STATS(1);
  F  = STATS(2);
  Fp = STATS(3); % for total regression
  clear STATS;
  
  [B,STATS] = robustfit(X(:,2:end),y);
  
  P = X*B;           % predicted values
  R = STATS.resid;   % residuals
  
  BCI = [B - 2*STATS.se, B + 2*STATS.se]; % beta CI
  YCI = [P - 2*STATS.s,  P + 2*STATS.s ]; % predicted CI
  
  % define figure rows & columns
  cols = 5;
  rows = size(X,2);
  
end



for fig = 1:2,
  
  % Setup figure
  H(end+1) = figure('color',[0 0 0]);
  pos = get(gcf,'Position');
  set(gcf,'Position',[pos(1)-(pos(3)/2) pos(2)-pos(4) pos(3)*2 pos(4)*2]);  % double width/height
  
  colormap(prism);
  fontsize = 8;
  fontcolor = [1 1 1]; % white
  
  % plot observed, predicted and 95% CI
  if fig < 2,
    subplot(2,3,1);
  else
    subplot(rows,cols,1);
  end
  if max(GP),
    scatter(P,y,5,GP, 'filled'); hold on;
  else,
    scatter(P,y,5,'m','filled'); hold on;
  end
  %scatter(P,P,5,'w','filled');
  plot(P,P,'w-',...
    P,YCI(:,1),'w:',...
    P,YCI(:,2),'w:');
  %legend('predicted','95% CI',0);
  set(gca,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1]);
  
  title{1} = sprintf('y'' = %8.2f',B(1));
  for j = 2:length(B),
    title{1} = sprintf('%s + %6.2f*X%d',title{1},B(j),j-1);
  end
  Htitle = get(gca,'title');
  set(Htitle,'string',title{1}, 'fontsize',fontsize,'color',fontcolor);
  
  xlab = get(gca,'xlabel');
  set(xlab,'string','Y predicted','fontsize',fontsize,'color',fontcolor);
  ylab = get(gca,'ylabel');
  set(ylab,'string','Y observed', 'fontsize',fontsize,'color',fontcolor);
  
  
  
  % plot observed, predicted and 95% CI against Xmean
  if fig < 2,
    subplot(2,3,2);
  else
    subplot(rows,cols,2);
  end
  % calculate mean X values
  Xmean = mean(X(:,2:end),2);
  % plot the observed & predicted values
  if max(GP),
    scatter(Xmean,y,5,GP, 'filled'); hold on;
  else,
    scatter(Xmean,y,5,'m','filled'); hold on;
  end
  scatter(Xmean,P,5,'w','filled');
  set(gca,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1]);
  xlab = get(gca,'xlabel');
  set(xlab,'string','mean Xi', 'fontsize',fontsize,'color',fontcolor);
  ylab = get(gca,'ylabel');
  set(ylab,'string','Y observed','fontsize',fontsize,'color',fontcolor);
  
  
  % curve fitting, cubic polynomial
  [ Xsort, k ] = sort(Xmean);
  [p,S] = polyfit(Xsort,P(k),3);
  polfit = polyval(p,Xsort,S);
  % % plot observed, polyfit and CI
  plot(Xsort,polfit,'w-');
  % legend('data','polyfit',0);
  
  
  % plot residuals against predicted
  if fig < 2,
    subplot(2,3,4);
  else
    subplot(rows,cols,3);
  end
  if max(GP),
    scatter(P,R,5,GP,'filled'); hold on
  else,
    scatter(P,R,5,'m','filled'); hold on
  end
  line(P,zeros(size(R)),'linestyle','-','color','w');
  set(gca,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1]);
  xlab = get(gca,'xlabel');
  set(xlab,'string','Y predicted','fontsize',fontsize,'color',fontcolor);
  ylab = get(gca,'ylabel');
  set(ylab,'string','residual', 'fontsize',fontsize,'color',fontcolor);
  
  % plot normal QQ plot of residuals
  if fig < 2,
    subplot(2,3,5);
  else
    subplot(rows,cols,4);
  end
  qqplot(R);
  QQtitle = get(gca,'title');
  set(QQtitle,'string','QQ plot','fontsize',fontsize,'color',fontcolor);
  xlab = get(gca,'xlabel');
  set(xlab,'string','','fontsize',fontsize,'color',fontcolor);
  ylab = get(gca,'ylabel');
  set(ylab,'string','','fontsize',fontsize,'color',fontcolor);
  set(gca,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1]);
  
  
  
  
  
  % print model statistics
  if robust,
    
    %    [B,STATS] = ROBUSTFIT(...) also returns a STATS structure
    %     containing the following fields:
    %         stats.ols_s     sigma estimate (rmse) from least squares fit
    %         stats.robust_s  robust estimate of sigma
    %         stats.mad_s     MAD estimate of sigma; used for scaling
    %                         residuals during the iterative fitting
    %         stats.s         final estimate of sigma, the larger of robust_s
    %                         and a weighted average of ols_s and robust_s
    %         stats.se        standard error of coefficient estimates
    %         stats.t         ratio of b to stats.se
    %         stats.p         p-values for stats.t
    %         stats.coeffcorr estimated correlation of coefficient estimates
    %         stats.w         vector of weights for robust fit
    %         stats.h         vector of leverage values for least squares fit
    %         stats.dfe       degrees of freedom for error
    %         stats.R         R factor in QR decomposition of X matrix
    
    if fig < 2,
      subplot(2,3,3);
    else
      subplot(rows,cols,5);
    end
    axis off;
    txt{1} = sprintf('%6s %8s %8s','','Beta','SE');
    for i = 1:length(B),
      if i == 1,
        predictor = 'CONST';
      else
        if i-1 <= length(labels),
          predictor = labels{i-1};
        end
      end
      txt{i+1} = sprintf('%6s %8.4f %8.4f',predictor,B(i),STATS.se(i));
    end
    text(-.2,.5,txt,'Fontname','courier','fontsize',fontsize,'color',fontcolor)
    
    txt = '';
    txt{1} = sprintf('%6s %8s %6s %6s %5s','','F','DF','sig','r^2');
    txt{2} = sprintf('%6s %8.4f %6.2f %6.4f %4.2f','ALL',F,STATS.dfe,Fp,R2);
    for i = 1:length(B),
      if i == 1,
        predictor = 'CONST';
      else
        if i-1 <= length(labels),
          predictor = labels{i-1};
        end
      end
      txt{i+2} = sprintf('%6s %8.4f %6.2f %6.4f',predictor,STATS.t(i),STATS.dfe,STATS.p(i));
    end
    
    % output the correlation matrix
    datacorr = corrcoef(data);
    rtxt = num2str(datacorr,' %+5.2f');
    txt2{1} = '    Correlation Matrix';
    txt2{2} = '   ';
    for r = 1:length(datacorr),
      if r < length(datacorr),
        txt2{2}   = sprintf('%5s %5s',txt2{2},sprintf('X%d',r));
        txt2{r+2} = sprintf('%5s %s',labels{r},rtxt(r,:));
      else
        txt2{2}   = sprintf('%5s %5s',txt2{2},'Y');
        txt2{r+2} = sprintf('%5s %s',labels{r},rtxt(r,:));
      end
    end
    
    if fig < 2,
      text(-.2,.1,txt ,'Fontname','courier','fontsize',fontsize,'color',fontcolor);
      text(-.1,1 ,txt2,'Fontname','courier','fontsize',fontsize,'color',fontcolor);
    else,
      text(-.2,.1,txt ,'Fontname','courier','fontsize',fontsize,'color',fontcolor);
      text(-.1,1 ,txt2,'Fontname','courier','fontsize',fontsize,'color',fontcolor);
    end
    txt = '';
    txt2 = '';
  else
    
    % stats(1) = R^2, stats(2) = F, stats(3) = p; % for total regression
    if fig < 2,
      subplot(2,3,3);
    else
      subplot(rows,cols,5);
    end
    axis off;
    txt{1} = sprintf('%12s %10s  %10s  %10s\n','','R^2','F','sig');
    txt{2} = sprintf('%12s%10.4f  %10.4f  %10.4f','MODEL',STATS(1),STATS(2),STATS(3));
    text(-.2,.9,txt,'Fontname','courier','fontsize',fontsize,'color',fontcolor)
    
    if fig > 1,
      subplot(rows,cols,6); rcoplot(R, RCI);
    end
  end
  
end





% plot observed, predicted and 95% CI against individual predictors
for i = 2:size(X,2),
  subplot(rows,cols,cols*i-(cols-1));
  % plot the observed values
  if max(GP),
    scatter(X(:,i),y,5,GP, 'filled'); hold on;
  else
    scatter(X(:,i),y,5,'w','filled'); hold on;
  end
  set(gca,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1]);
  
  % Calculate predicted values for just this predictor
  beta = [B(1);B(i)];
  P = [X(:,1),X(:,i)]*beta;
  R = y - P;
  YCI = [P - 2*STATS.s,  P + 2*STATS.s ]; % predicted CI
  
  %scatter(X(:,i),P,5,'w','filled');
  plot(X(:,i),P,'w-',...
    X(:,i),YCI(:,1),'w:',...
    X(:,i),YCI(:,2),'w:');
  %legend('predicted','95% CI',0);
  
  title{1} = sprintf('y'' = %8.2f',B(1));
  title{1} = sprintf('%s + %6.2f*{\\bfX%d}',title{1},B(i),i-1);
  %Htitle = get(gca,'title');
  %set(Htitle,'string',title{1},'interpreter','tex');
  ylab = get(gca,'ylabel');
  set(ylab,'string',title{1},'fontsize',fontsize,'color',fontcolor);
  
  
  
  
  % plot residuals against predicted
  subplot(rows,cols,cols*i-(cols-3));
  if max(GP),
    scatter(X(:,i),R,5,GP, 'filled'); hold on;
  else
    scatter(X(:,i),R,5,'w','filled'); hold on;
  end
  line(X(:,i),zeros(size(R)),'linestyle','-','color','w');
  %xlab = get(gca,'xlabel');
  %set(xlab,'string','Xi','fontsize',fontsize,'color',fontcolor);
  ylab = get(gca,'ylabel');
  set(ylab,'string','residual','fontsize',fontsize,'color',fontcolor);
  set(gca,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1]);
  
  
  
  % plot normal QQ plot of residuals
  subplot(rows,cols,cols*i-(cols-4));
  qqplot(R);
  QQtitle = get(gca,'title');
  set(QQtitle,'string','QQ plot','fontsize',fontsize,'color',fontcolor);
  xlab = get(gca,'xlabel');
  set(xlab,'string','','fontsize',fontsize,'color',fontcolor);
  ylab = get(gca,'ylabel');
  set(ylab,'string','','fontsize',fontsize,'color',fontcolor);
  set(gca,'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1]);
end


return
