%% Exploration analysis code: activity, area coverage, straightness over & w/in days
% Analysis code for "Ant colonies explore novel environments with increased
% activity and slower, curvier walks."
% Stefan Popp & Anna Dornhaus, Sep 2023
% For stats, see the R script

% Requires Mapping toolbox for gaze analysis (wrapTo180 function)
% violin function:
% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de

% 1. Creates 'dataClean' & 'dataWithWalls' file of all tracks from all
% individual trial track files OR: loads them in
% 2. Calculates local straightness values
% 3. Creates figures

% Glossary
%{
ants = dataClean file. Contains all track data. Each row is 1 point
col = colony. Either letter name (T-Y) or number 1-5, respectively
idUni = unique track ID across trials
s = step length (distance between 2 points)
v = speed
nestDispk = Euclidean distance from that point to the nest [mm]
st = local straightness
nu = angle relative to nest
%}

% === User input === %
% You can make the analysis data files from the raw track data provided
% (takes a few minutes), OR load in the provided analysis data files made,
% which should be identical to what you get running the data analysis maker
newData = false; % true: make data from raw files; flase: load it in
% === End of user input === %

% The rest of the script needs no further input. Just run each section
% (ctrl+enter) to produce data or figures reported in the paper.

ftsz = {'fontsize',12}; % Of plot labels
lw = {'LineWidth', 1}; % Line plot line widths

if ispc; slash = '\'; else; slash = '/'; end % OS path compatibility
topFolderDir = 'AntExploration'; % Top-level folder (downloaded)
inputFolderDir = [topFolderDir slash 'dataRaw']; % Input ant tracks
plotPath = [topFolderDir slash 'plots']; % Where plots will be saved

% Adds folders to path for MATLAB to access them
if ~isfolder(plotPath); mkdir(plotPath); end % Makes plot folder if missing
addpath(genpath(topFolderDir))

%% Make or load raw data file (depending on newData)
% Collating all raw data + adding straightness

cols = ["T" "U" "V" "W" "Y"]; % List of colony names to be included
colsCh = cell2mat(convertStringsToChars(cols)); % Char instead of string

if newData
    w = 21; % Window for local straightness
    wBar = waitbar(0,'Loading & prepping','Name','Making analysis data');
    % Load all 3 days of each colony into the same table
    antsC = cell(numel(cols)*3,1); % Preallocating to save computing time
    for col = 1:numel(cols)
        for day = 1:3
            waitbar((col*3-3+day)/(numel(cols)*3), wBar, ['Col ' colsCh(col) num2str(day)])
            antsCurr = readtable(['HRM_' colsCh(col) num2str(day) '_ants']); % Loads in individual raw tracks
            load(['HRM_' colsCh(col) num2str(day) '_procParams']) % Loads parameter file for nestXY position
            
            % We sweep all track gaps (wall, tracker not picking up,
            % artefacts) into one row, to count them in activity but not
            % speed & straightness calculations
            antsCurr.gap(antsCurr.s>2.3) = 1; % gap points not included in st calculation
            antsCurr.gap(abs(antsCurr.a)>100 | antsCurr.v>25) = 1; % If there are still artefacts left
            antsCurr.gap(antsCurr.a<-40 & antsCurr.v>17) = 1; % Deceleration & still being fast is unlikely
            antsCurr{end+1,:} = nan; % Prepping for next 2 lines
            antsCurr.s(end) = nansum(antsCurr.s(antsCurr.gap==1)); % Last row is s of all wall parts for activity
            antsCurr.s(antsCurr.gap==1) = NaN; % Would be a double count; AND for ST calculation below
            antsCurr.gap(end) = 1;
            
            % Adding colony & day columns
            antsCurr.col = repmat(col,height(antsCurr),1);
            antsCurr.day = repmat(day,height(antsCurr),1);
            antsCurr.st = localSTfun(antsCurr,w); % Local straightness: euclid. distance between 2 pts/walked length
            
            % === For Gaze analysis === %
            % Angle ant heading to nest
            nestAngle = atan2d(antsCurr.x-params.nestXY(1), antsCurr.y-params.nestXY(2));
            antsCurr.nu = wrapTo180(nestAngle + antsCurr.theta);
            
            antsCurr = array2table(single(antsCurr{:,:}),'variablenames',antsCurr.Properties.VariableNames); % Saving storage
            antsCurr = sortrows(antsCurr,{'col' 'day' 'id' 't'});
            antsC{col*3-3+day,1} = antsCurr;
        end
    end
    waitbar(1, wBar, 'Saving files')
    tracks = vertcat(antsC{:}); % All 3 days of all 5 colonies
    [~,~,tracks.idUni] = unique(tracks(:,{'col' 'day' 'id'}),'rows');
    
    tracks(:,{'a' 'disp' 'theta'}) = []; % Unnecessary columns

    % Saving data
    writetable(tracks,[topFolderDir slash 'dataWithWalls']) % For fig 1 activity/dispersivity: Includes tracking artefacts & walls
    tracks(tracks.gap==1,:) = []; % For rest: kick out tracking artefacts & walls
    writetable(tracks,[topFolderDir slash 'dataClean'])
    close(wBar)
    clearvars -except ftsz inputFolderDir lw plotPath slash topFolderDir
    newData = false;
end

% Loading in data
cols = ["T" "U" "V" "W" "Y"]; % List of colony names to be included
antsWalls = readtable([topFolderDir slash 'dataWithWalls']);
tracks = readtable([topFolderDir slash 'dataClean']);
wallDim = readmatrix([inputFolderDir slash 'wallDim.txt']);

%% General results (walls, #of ants)
% Points near walls
disp([num2str(nansum(antsWalls.wall)/height(antsWalls)*100,3) '% of points are near walls (= excluded)'])

% #of simultaneous ants
for col = 1:5; tmax(col) = max(antsWalls.t(antsWalls.col==col)); end
clear antCountStr
i = 1;
for day = 1:3
    for col = 1:5
        idx = antsWalls.day == day & antsWalls.col == col;
        antCountStr.nrAnts(i,1) = numel(unique(antsWalls.id(idx & antsWalls.t>min(tmax)-10)));
        antCountStr.day(i,1) = day;
        antCountStr.col(i,1) = col;
        i = i+1;
    end
end
antCount = struct2table(antCountStr);
antCountMuSD = groupsummary(antCount,'day',{'mean' 'std'},'nrAnts')

% boxplot(antCount.nrAnts, antCount.day)
% xlabel('Days')
% ylabel('Number of ants')

% Repeated bouts per ant
disp([num2str(numel(unique(tracks.idUni(tracks.nestDisp<10))) / 750,3) ' tracks per ant'])

%% Data for Fig 1: Pixels of #of points
antsWalls.xBin = round(antsWalls.x,-1);
antsWalls.yBin = round(antsWalls.y,-1);
cHeat = groupsummary(antsWalls,{'xBin' 'yBin' 'day'}); % For a)

%%       Figure 1: Heatmaps of point densities

% c heatmap
f1 = figure('units','centimeters','Outerposition',[.2 1.5 12 20]);
ti = tiledlayout(3,1,"TileSpacing","compact",Padding="compact");
ax = gobjects(3,1);
for day = 1:3
    ax(day) = nexttile;
    idx = cHeat.day == day;
    scatter(cHeat.xBin(idx),cHeat.yBin(idx),3,log(cHeat.GroupCount(idx)),'s')
    rectangle('Position',wallDim) % Rectangle indicating included data
    % hold on % Nest circle
    % scatter(nestXY.x(2),nestXY.y(2),20,'w',LineWidth=1)
    % hold off 
    xlim([0 2680])
    ylim([0 1800])
    title(['Day ' num2str(day)])
end
ylabel(ti,'y [mm]',ftsz{:})
xlabel(ti,'x [mm]',ftsz{:})
ax(1).XTick = [];
ax(2).XTick = [];

set(ax, 'Colormap', jet, 'CLim', [0 max(log(cHeat.GroupCount))])
cb = colorbar(ax(2));
cb.Label.String = 'Log number of points per cm^2';
cb.Label.FontSize = ftsz{2};
cb.Layout.Tile = 'east';

annotation('textbox',[.18 .94 0 0],'String','A',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.18 .62 0 0],'String','B',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.18 .31 0 0],'String','C',"LineStyle","none",ftsz{:},'FontWeight',"bold")

exportgraphics(f1,[plotPath slash 'f1_heat.png'],'Resolution',300) % For initial subs only
% exportgraphics(f1,[plotPath slash 'f1_heat.pdf'],'Resolution',300,'contenttype','vector')

%% Data for Fig 2: Distance walked, grouped by colony over days; Dispersivity
% Estimated number of exploration bouts (multiple bouts per ants)

disp(['Total ' num2str(nansum(antsWalls.s)/1000000,4) ' km walked'])
disp(['Of that, ' num2str(nansum(antsWalls.s(antsWalls.wall==1))/1000000,3) 'm near the walls.'])

activ = groupsummary(antsWalls,{'col' 'day'},'sum');
disp([num2str(nanmean(activ.sum_s)/1000000,4) ' ± ' num2str(nanstd(activ.sum_s)/1000000,3),...
        ' (µ ± SD) km walked per colony and day'])

nrOfAnts = numel(unique(antsWalls.idUni(antsWalls.nestDisp<30 & [diff(antsWalls.idUni)~=0; false])));
disp([num2str(nrOfAnts) ' tracks start from the nest'])
disp(['That is on average ' num2str(nrOfAnts/75,4) ' per hour'])

dispTab = groupsummary(antsWalls,{'col' 'day'}, 'mean','nestDisp'); % For dispersivity col means (2left)

% For plots
activhc = zeros(200,3);
activhcEdges = zeros(201,3);
for day = 1:3
    [activhc(:,day), activhcEdges(:,day)] = histcounts(antsWalls.nestDisp(antsWalls.day==day),200);
end
activhc(:,4) = activhc(:,1) - activhc(:,2);
activhc(:,5) = activhc(:,2) - activhc(:,3);
activhc(:,6) = activhc(:,1) - activhc(:,3);
actComp = array2table(activhc, 'variablenames',{'Day1' 'Day2' 'Day3' 'Day1-2' 'Day2-3' 'Day1-3'});

%%       Figure 2: Dispersivity over days & differences thereof (hist)
f2 = figure('units','centimeters','Outerposition',[.2 1 17.7 18]);
ti = tiledlayout(3,2,"TileSpacing","compact",'padding','compact');

ax = gobjects(3,2);
% Left column
for day = 1:3
    ax(day*2-1) = nexttile(day*2-1);
    bar(activhcEdges(1:end-1,day),actComp{:,day}, 'facecolor',[.5 .5 .5])
    hold on % Mean nestDisp per colony as vertical colored bars
    for col = 1:numel(cols)
        pl(col) = plot(dispTab.mean_nestDisp(dispTab.day==day & dispTab.col==col).*ones(2,1), [0 120000],LineWidth=1);
    end
    hold off
    title(['Day ' num2str(day)])
    xlim([0 max(antsWalls.nestDisp)])
end
% Legend
leg = legend(pl,cols);
title(leg,'Colonies')
leg.ItemTokenSize = [10;18]; % Shortens lines in legend

% Right column
for day = 1:3
    ax(day*2) = nexttile(day*2);
    bar(activhcEdges(1:end-1,day),actComp{:,day+3}, 'facecolor',[.5 .5 .5])
    ax(day*2).YAxisLocation = 'right';
    title(actComp.Properties.VariableNames{day+3}) % E.g. Day2-1
    xlim([0 max(antsWalls.nestDisp)])
end

% Left column (absolute)
ax(1).XTickLabel = '';
ax(3).XTickLabel = '';
ax(5).XLabel.FontSize = ftsz{2};
ax(3).YLabel.String = 'Number of points';
ax(3).YLabel.FontSize = ftsz{2};
% Right column (differences between days)
ax(2).XTickLabel = [];
ax(4).XTickLabel = [];
ax(4).YAxis.Exponent = 4;
ax(4).YLim = [0 100000];
ax(6).XLabel.FontSize = ftsz{2};

linkaxes([ax(:)],'y')
xlabel(ti,'Distance from the nest [mm]');

% Panel letters
annotation('textbox',[.07 .94 0 0],'String','A',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.07 .63 0 0],'String','B',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.07 .34 0 0],'String','C',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.55 .94 0 0],'String','D',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.55 .63 0 0],'String','E',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.55 .34 0 0],'String','F',"LineStyle","none",ftsz{:},'FontWeight',"bold")

exportgraphics(f2,[plotPath slash 'f2_disp.png'],'Resolution',300) % For initial subs only
% exportgraphics(f2,[plotPath slash 'f2_activ_disp.pdf'],'Resolution',300,'contenttype','vector')

%%       >>Legacy<< Figure 2: Activity & dispersivity (requires running of 'distance walked' section above)
% f1 = figure('units','centimeters','Outerposition',[.2 1 17.7 15]);
% ti = tiledlayout(3,4,"TileSpacing","tight",'padding','tight');
% nexttile([2 2])
% hold on
% for col = 1:5; plot(activ.sum_s(activ.col==col)/1000000,lw{:}); end
% hold off
% xlim([0.9 3.1])
% xticks([1 2 3])
% xlabel('day',ftsz{:})
% ylabel('Total walked distance [km]',ftsz{:})
% leg = legend(cols);
% title(leg,'Colonies')
% leg.ItemTokenSize = [10;18]; % Shortens lines in legend
% 
% ax = cell(3,1);
% for day = 1:3
%     ax{day} = nexttile(4*day-1,[1 2]);
%     histogram(antsWalls.nestDisp(antsWalls.day==day),200,"FaceColor",[.5 .5 .5], 'HandleVisibility','off')
%     hold on % Mean nestDisp per colony as vertical colored bars
%     for col = 1:numel(cols)
%         plot(dispTab.mean_nestDisp(dispTab.day==day & dispTab.col==col).*ones(2,1), [0 120000],LineWidth=1)
%     end
%     hold off
%     title(['Day ' num2str(day)])
%     xlim([0 max(antsWalls.nestDisp)])
% end
% ax{1}.XTickLabel = '';
% ax{2}.XTickLabel = '';
% linkaxes([ax{1} ax{2} ax{3}],'y')
% ax{3}.XLabel.String = 'Distance from the nest [mm]';
% ax{3}.XLabel.FontSize = ftsz{2};
% ax{2}.YLabel.String = 'No. of points';
% ax{2}.YLabel.FontSize = ftsz{2};
% % xlabel(ti,'Distance from the nest [mm]',ftsz{:})
% % ylabel(ti,'Nr of points',ftsz{:})
% 
% annotation('textbox',[.09 .95 0 0],'String','A',"LineStyle","none",ftsz{:},'FontWeight',"bold")
% annotation('textbox',[.93 .94 0 0],'String','B',"LineStyle","none",ftsz{:},'FontWeight',"bold")
% annotation('textbox',[.93 .61 0 0],'String','C',"LineStyle","none",ftsz{:},'FontWeight',"bold")
% annotation('textbox',[.93 .28 0 0],'String','D',"LineStyle","none",ftsz{:},'FontWeight',"bold")

%% Legacy Data for Fig 3 intra-day: binning
% For within days: Binning into 50 bins of equal width
[~,~,tracks.tBin] = histcounts(tracks.t,30);
lineSum = groupsummary(tracks,{'day','tBin'},'median');
lineSum = array2table(double(lineSum{:,:}),'variablenames',lineSum.Properties.VariableNames); % For regressions

%% Data for fig 3
tracks.nearFar(tracks.nestDisp>=1000) = 1; % 0 for <1m, 1 for >1m
vMedians = groupsummary(tracks,{'day' 'nearFar'},'median',{'v' 'st'});
% Wide format for bar3
wide = zeros(3,2,2);
for vari = 1:2
    for day = 1:3
        for nearFar = 1:2
            wide(day,nearFar,vari) = vMedians.(['median_' char(vars(vari))])(vMedians.day==day & vMedians.nearFar==nearFar-1);
        end
    end
end

%%      Figure 3 A-L: Violins between days & near/far
vars = ["v";"st"];
med = zeros(12,1);
f3_1 = figure('units','centimeters','Outerposition',[.2 1 17.4 16.4]);
ti = tiledlayout(3,4, 'Tilespacing','tight');
i = 1;
for day = 1:3
    for vari = 1:numel(vars)
        for nearFar = 1:2
            tn = tilenum(ti,day,vari*2-2+nearFar);
            ax(tn) = nexttile(tn);
            trCurr = tracks(tracks.nearFar==nearFar-1 & tracks.day==day,vars(vari));
            [~,~,~,med(i,:)] = violin(trCurr.(vars(vari)));
            ax(tn).YLabel.String = '';
            ax(tn).XLabel.String = '';
            ax(tn).XTick = [];
            i = i+1;
        end
    end
end

xlabel(ax([9 11]),'nest\_dist<1m',ftsz{:})
xlabel(ax([10 12]),'nest\_dist>=1m',ftsz{:})
ylabel(ax(5),'Speed [mm/s]',ftsz{:})
ylabel(ax(8),'Straightness',ftsz{:})
set(ax([2 3 6 7 10 11]),'ytick','')
set(ax([4 8 12]),'yaxislocation','right')
ax(8).YAxisLocation = 'right';

% Vertical divider speed from straightness
annotation('line',[0.515 0.515],[0.11 0.93],'LineWidth',1);

% Median value annotations
medPos = [repmat([0.23; 0.42; 0.62; 0.81],3,1),... % x positions
          repelem([0.79; 0.90; 0.51; 0.63; 0.25; 0.37],2,1),... % y posiitons
          zeros(12,2)]; % Extent of the (nonexistent) box
for i = 1:12
    annotation('textbox',medPos(i,:),'String',num2str(med(i),2),"LineStyle","none",ftsz{:})
    % Panel letters
    annotation('textbox',ax(i).Position,'String',char(i+64),"LineStyle","none",ftsz{:},'FontWeight',"bold")
end
annotation('textbox',[.47 .8 .1 .05],'String','Day 1',"LineStyle","none",ftsz{:},'FontWeight',"bold",'backgroundcolor','w')
annotation('textbox',[.47 .51 .1 .05],'String','Day 2',"LineStyle","none",ftsz{:},'FontWeight',"bold",'backgroundcolor','w')
annotation('textbox',[.47 .27 .1 .05],'String','Day 3',"LineStyle","none",ftsz{:},'FontWeight',"bold",'backgroundcolor','w')

exportgraphics(f3_1,[plotPath slash 'f3A_violin.png'],'Resolution',300) % For initial subs only
% exportgraphics(f3_1,[plotPath slash 'f3A_violin.pdf'],'Resolution',300,'contenttype','vector')

%%      Figure 3M,N: 3D Bar plot
f3_2 = figure('units','centimeters','Outerposition',[.2 1 17.4 8.2]);
ti = tiledlayout(1,2, 'Tilespacing','tight');
for vari = 1:2
    % ax(12+vari) = nexttile(11+vari*2,[1 2]);
    ax(vari) = nexttile(vari);
    wdCurr = wide(:,:,vari); % Data of current variable (st or v)
    b = bar3(wdCurr, 1);
    rnge = range(wdCurr,'all');
    zlim([min(wdCurr,[],'all')-rnge, max(wdCurr,[],'all')+rnge]) % Zooming in
    set(gca, 'YDir','reverse','XTickLabel',{'<1m' '>=1m'}) % To make days increase towards the bottom
    b(1).FaceColor = [.5 .5 .5];
    b(2).FaceColor = [0.8500 0.3250 0.0980]; % Red
    b(1).FaceAlpha = 0.5; % Front bars transparent
    % xlabel('Distance to the nest',ftsz{:})
    ylabel('Day',ftsz{:})
    zlabel(vars(vari),ftsz{:})
end
zlabel(ax(1),'Speed [mm/s]',ftsz{:})
zlabel(ax(2),'Straightness',ftsz{:})
annotation('textbox',ax(1).Position,'String',char(13+64),"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',ax(2).Position,'String',char(14+64),"LineStyle","none",ftsz{:},'FontWeight',"bold")

exportgraphics(f3_2,[plotPath slash 'f3M_3D.png'],'Resolution',300) % For initial subs only
% exportgraphics(f3_2,[plotPath slash 'f3M_3D.pdf'],'Resolution',300,'contenttype','vector')

%%       Figure 3: v & ST between days (violin) & within days (line)
vars = ["v";"st"];
ylabelList = ["Speed [mm/s]";"Local straightness"];

f3 = figure('units','centimeters','Outerposition',[.2 1 17.7 23.8]);
ti = tiledlayout(4,3, 'Tilespacing','tight');
ax = cell(12,2);

% Between days
for vari = 1:numel(vars)
    ax{vari} = nexttile([1 3]);
    plotC = splitapply( @(x){x}, tracks{:,vars(vari)}, tracks.day);
    [~,~,~,med(vari,:)] = violin(plotC');
    xticks([1 2 3])
    ylabel(ylabelList(vari),ftsz{:})
end
ax{1}.XTickLabel = [];
ax{2}.XTickLabel = [];
% xlabel('Day',ftsz{:})
annotation('textbox',[.15 .9 0 0],'String','A',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.15 .7 0 0],'String','B',"LineStyle","none",ftsz{:},'FontWeight',"bold")
% Significance letters
annotation('textbox',[.25 .92 0 0],'String','b',"LineStyle","none",ftsz{:})
annotation('textbox',[.52 .92 0 0],'String','a',"LineStyle","none",ftsz{:})
annotation('textbox',[.77 .92 0 0],'String','c',"LineStyle","none",ftsz{:})
annotation('textbox',[.24 .72 0 0],'String','a',"LineStyle","none",ftsz{:})
annotation('textbox',[.5 .72 0 0],'String','b',"LineStyle","none",ftsz{:})
annotation('textbox',[.76 .72 0 0],'String','c',"LineStyle","none",ftsz{:})
% Median value annotation
annotation('textbox',[.26 .83 0 0],'String',num2str(med(1,1),2),"LineStyle","none",ftsz{:})
annotation('textbox',[.52 .83 0 0],'String',num2str(med(1,2),2),"LineStyle","none",ftsz{:})
annotation('textbox',[.78 .83 0 0],'String',num2str(med(1,3),2),"LineStyle","none",ftsz{:})
annotation('textbox',[.26 .69 0 0],'String',num2str(med(2,1),2),"LineStyle","none",ftsz{:})
annotation('textbox',[.52 .69 0 0],'String',num2str(med(2,2),2),"LineStyle","none",ftsz{:})
annotation('textbox',[.78 .69 0 0],'String',num2str(med(2,3),2),"LineStyle","none",ftsz{:})

% Within days
for vari = 1:numel(vars)
    for day = 1:3
        ax{6+ vari*3-3 + day} = nexttile;
        idxDay = tracks.day == day;
        if vars(vari) == "st"
            idxDay(isnan(tracks.st)) = false;
        end
        [N,C] = hist3([tracks{idxDay,"t"},tracks{idxDay,vars(vari)}],[100,100]);
        maxN = max(max(N));
        contourf(C{1},C{2},N',1:maxN/10:maxN,'LineStyle','none')
        colormap(flipud(pink)); % clim([-200 maxN+100]) % Lowest level still shaded
        hold on
        plot(lineSum{lineSum.day==day,"median_t"},lineSum{lineSum.day==day,append("median_",vars(vari))},'k',lw{:})
        hold off
    end
end
ax{8}.YTick = []; ax{11}.YTick = [];
ax{9}.YAxisLocation = 'right'; ax{12}.YAxisLocation = 'right';
ax{7}.YLabel.String = ylabelList(1);
ax{7}.YLabel.FontSize = ftsz{2};
ax{10}.YLabel.String = ylabelList(2);
ax{10}.YLabel.FontSize = ftsz{2};
ax{8}.YLabel.String = ''; ax{3}.YLabel.String = '';
ax{11}.YLabel.String = ''; ax{6}.YLabel.String = '';
ax{1}.Title.String = 'Day 1                       Day 2                       Day 3';
ax{1}.Title.FontSize = ftsz{2};
ax{7}.XTick = [];
ax{8}.XTick = [];
ax{9}.XTick = [];
% ax{2}.Title.String = 'Day 2'; ax{2}.Title.FontSize = ftsz{2};
% ax{3}.Title.String = 'Day 3'; ax{3}.Title.FontSize = ftsz{2};
xlabel(ti,'Time [s]',ftsz{:})
linkaxes([ax{7} ax{8} ax{9}],'y')
linkaxes([ax{10} ax{11} ax{12}],'y')
annotation('textbox',[.15 .5 0 0],'String','C',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.15 .25 0 0],'String','D',"LineStyle","none",ftsz{:},'FontWeight',"bold")

% exportgraphics(f3,[plotPath slash 'f3_vST_days.png'],'Resolution',300) % For initial subs only
% exportgraphics(f3,[plotPath slash 'f3_vST_days.pdf'],'Resolution',300,'contenttype','vector')

%% Data for Fig 4: Gaze analysis
% Do gazes towards the nest contain more stops than gazes in other directions?
vThresh = 2; % [s] stop = ant takes >vThresh to move the 2mm between 2 pts
tracks.stop = tracks.dt > vThresh;

% Only first ones
% Time of first ant entering nest for each colony
nestStartIdC = cell(size(cols'));
for col = 1:numel(cols)
    tMax = min(tracks.t(tracks.nestDisp<50 & [diff(tracks.idUni)~=0; false] & tracks.col==col & tracks.day==1));
    nestStartIdC{col,1} = tracks.idUni(tracks.col == col & tracks.nestDisp<50 & tracks.t<tMax & [false; diff(tracks.idUni)~=0]);
end
nestStartId = vertcat(nestStartIdC{:});
antsFirst = tracks(ismember(tracks.idUni, nestStartId),:);
antsFirst.nuBin = round(antsFirst.nu./30).*30;
disp(['n = ' num2str(sum(antsFirst.stop)) ' stops and '...
             num2str(sum(~antsFirst.stop)) ' non-stops in '...
             num2str(numel(nestStartId)) ' tracks.'])
         
% Is the angular distribution of stops different from that of the #of
% non-stop points?
% Watson's U^2 test from Pierre Mégevand (2023).
% https://github.com/pierremegevand/watsons_u2, GitHub. Retrieved January 31, 2023.
% Could be done in R package 'circular' to get test stats
pGaze = watsons_U2(antsFirst.nu(antsFirst.stop==1),...
                   antsFirst.nu(antsFirst.stop==0))

%%       Figure 4A: No conspicuous patterns
% ToDo: merge w/ 3b&c; add letters; y-axis 0
gazePropFirst = groupsummary(antsFirst,'nuBin','mean');

% A: Proportions of stops by nest-gaze angle
f4 = figure('units','centimeters','Outerposition',[.2 1 17.7 13]);
ti = tiledlayout(4,1,'tilespacing','none');
ax{1} = nexttile([3 1]);
scatter(gazePropFirst.nuBin, gazePropFirst.mean_stop.*100,10,'filled','k')
ylabel({'Proportion of'; 'points being stops'},ftsz{:})
ytickformat('percentage')
ylim([0 0.65])
xlim([-180 180])

% Tile below: sample size histogram + stops as vertical lines
ax{2} = nexttile;
histogram(antsFirst.nu,'facecolor',[.5 .5 .5],'edgecolor','none')
hold on
xline(antsFirst.nu(antsFirst.stop==1),'b')
hold off
xlim([-180 180])
xlabel('Gaze direction (0°=towards nest)',ftsz{:})
ylabel({'Nr of steps into';'that direction'}, ftsz{:})

ax{1}.YTick(1) = ''; % Removes 0, which overlaps w/ bottom histogram y-tick label
annotation('textbox',[.15 .9 0 0],'String','A',"LineStyle","none",ftsz{:},'FontWeight',"bold")

exportgraphics(f4,[plotPath slash 'f4A_gaze.png'],'Resolution',300) % For initial subs only
% exportgraphics(f4,[plotPath slash 'f4A_gaze.tif'],'Resolution',300,'contenttype','vector')

%%       Figure 4B,C Plotting 5 first tracks in col W, day 1, starting @ the nest, first 30cm
[~,stInd] = findgroups(tracks.idUni); % Start index
[~,stInd] = ismember(stInd,tracks.idUni);
stIdx = false(size(tracks.x)); for i = 1:numel(stInd); stIdx(stInd(i)) = true; end % Logical index vector
trackStarts = tracks(stIdx & tracks.nestDisp < 50 & tracks.col==4,:); % First points of tracks starting @ nest
startSum = groupsummary(tracks(ismember(tracks.idUni,trackStarts.idUni),:),'idUni','sum','t');
trackStarts(ismember(trackStarts.idUni,startSum.idUni(startSum.GroupCount < 150)),:) = []; % Sorting out too short tracks

trackStartsFirst = trackStarts(trackStarts.day==1,:);
trackStartsLast = trackStarts(trackStarts.day==3,:);
[~,frstInd] = mink(trackStartsFirst.t,5);
[~,lastInd] = maxk(trackStartsLast.t,5);

load('HRM_W1_procParams'); W1nestXY = params.nestXY./params.pxpermm; % nest
load('HRM_W3_procParams'); W3nestXY = params.nestXY./params.pxpermm;

f4b = figure('units','centimeters','Outerposition',[.2 1 17.7 10]);
ti = tiledlayout(1,2, 'Tilespacing','tight');
nexttile
hold on
for id = 1:numel(frstInd)
    idx = tracks.idUni == trackStartsFirst.idUni(frstInd(id));
    plot(tracks.x(idx),tracks.y(idx),lw{:})
end
axis([1200 2700 0 1200])
pbaspect([1.25 1 1])
viscircles(W1nestXY,3);
hold off

nexttile
hold on
for id = 1:numel(lastInd)
    idx = tracks.idUni == trackStartsLast.idUni(lastInd(id));
    plot(tracks.x(idx),tracks.y(idx),lw{:})
end
axis([1200 2700 0 1200])
pbaspect([1.25 1 1])
viscircles(W1nestXY,3);
hold off
set(gca,'YTickLabel',[]);
xlabel(ti,'x [mm]',ftsz{:}); ylabel(ti,'y [mm]',ftsz{:})

annotation('textbox',[.15 .9 0 0],'String','B',"LineStyle","none",ftsz{:},'FontWeight',"bold")
annotation('textbox',[.55 .9 0 0],'String','C',"LineStyle","none",ftsz{:},'FontWeight',"bold")

% exportgraphics(f4b,[plotPath slash 'f4b_tracks.png'],'Resolution',300) % For initial subs only
