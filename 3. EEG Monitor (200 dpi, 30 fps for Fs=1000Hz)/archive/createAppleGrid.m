function createAppleGrid(outFile)
% Apple-style ECG grid, 10.2 : 1, 틈 없이

%% ❶  해상도(가로) = 255 × 정수 ------------------------------------------
major = 0.2;  minor = 0.04;          % grid spacing (s)
pxPerMinor = 20;                     % 한 셀당 px
nMinor     = 255;                    % 10.2 / 0.04
wPx        = nMinor * pxPerMinor;    % 255 minors horizontally
nMinorY    = round(1/minor);         % 1 s / 0.04 s = 25 minors vertically
hPx        = nMinorY * pxPerMinor;   % exact multiple of pxPerMinor

%%  Grid 색·간격 -----------------------------------------------------------
majorC = [173 172 176]/255;
minorC = [236 236 236]/255;
secC = [140 140 140]/255;           % Color for 1-second grid lines (every 5 major)

%% ❷  Figure ----------------------------------------------------------------
dpi = 150;
fig = figure('Pos',[100 100 wPx hPx],'Color','w', ...
             'InvertHardcopy','off', 'GraphicsSmoothing','off'); % ❸ alias 방지
ax  = axes(fig,'Pos',[0 0 1 1],'Color','w'); hold(ax,'on');
axis(ax,'off','equal')      % ❹ 1 data unit = 1 pixel 스케일 고정

ax.XLim = [0 10.2];
ax.YLim = [0 1];

%% ❺  Grid ------------------------------------------------------------------
xVals = 0:minor:10.2;

% Draw minor grid lines first (bottom layer)
for xm = xVals
    if abs(mod(xm,major)) >= 1e-12  % Only minor grid lines
        line([xm xm], ax.YLim, 'Color', minorC, 'LineWidth',0.6);
    end
end

yVals = 0:minor:1;
for yy = yVals
    if abs(mod(yy,major)) >= 1e-12  % Only minor grid lines
        line(ax.XLim, [yy yy], 'Color', minorC, 'LineWidth',0.6);
    end
end

% Draw major grid lines (middle layer)
for xm = 0:major:10.2
    % Skip the 1-second lines as they'll be drawn later with different style
    if abs(mod(xm,1)) >= 1e-12  % Not a 1-second mark
        line([xm xm], ax.YLim, 'Color', majorC, 'LineWidth',1.2);
    end
end

for yy = 0:major:1
    line(ax.XLim, [yy yy], 'Color', majorC, 'LineWidth',1.2);
end

% Draw 1-second grid lines (top layer, thickest)
for xm = 0:1:10.2
    line([xm xm], ax.YLim, 'Color', secC, 'LineWidth',1.8);
end

%%  Export -----------------------------------------------------------------
if nargin < 1, outFile = 'apple_grid.png'; end
exportgraphics(fig,outFile,'Resolution',dpi);
end
