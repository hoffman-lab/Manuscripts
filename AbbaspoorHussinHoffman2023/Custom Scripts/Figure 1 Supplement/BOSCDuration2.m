clear all, clc

fh = figure()
fh.WindowState = 'maximized';
ax1 = axes('Position',[0.1 0.55 0.35 0.35]);
ax2 = axes('Position',[0.55 0.55 0.35 0.35]);
tprob = [];

for Condition = 1:2
    if Condition == 1
        load('D:\HiTGun\Data\BOSC_LE_RestSearch_forDurDist')
        load('D:\HiTGun\Data\BOSC_LU_RestSearch_forDurDist')
        
        BOSC = [LEBOSC; LUBOSC];
        RestThetaDur = cellfun(@(x,y) cat(2, x,y), BOSC.PreRestThetaDur, BOSC.PostRestThetaDur, 'UniformOutput', false);
        BOSC.RestThetaDur = RestThetaDur;
        RestGammaDur = cellfun(@(x,y) cat(2, x,y), BOSC.PreRestGammaDur, BOSC.PostRestGammaDur, 'UniformOutput', false);
        BOSC.RestGammaDur = RestGammaDur;
        
        [pvalue_theta,h] = ranksum(cell2mat([BOSC.RestThetaDur]'), cell2mat([BOSC.SearchThetaDur]'));
        [pvalue_gamma,h] = ranksum(cell2mat([BOSC.SearchGammaDur]'), cell2mat([BOSC.RestGammaDur]'));
        
    elseif Condition == 2
        load('D:\HiTGun\Wi TH Data\Analyses\BOSC\BOSC_Wi_RestSearch_forDurDist.mat')
        load('D:\HiTGun\LU TH Data\Analyses\BOSC\BOSC_LU_RestSearch_forDurDist.mat')
        
        BOSC = [WIBOSC; LUBOSC];
        
        [pvalue_theta,h] = ranksum(cell2mat([BOSC.RestThetaDur]'), cell2mat([BOSC.SearchThetaDur]'));
        [pvalue_gamma,h] = ranksum(cell2mat([BOSC.SearchGammaDur]'), cell2mat([BOSC.RestGammaDur]'));
    end
    
    for PLOT = 1:4
        switch PLOT
            case 1
                RESTBOUTS = cell2mat([BOSC.RestThetaDur]');
                BOUTS = cell2mat([BOSC.RestThetaDur]');
                pd = fitdist(RESTBOUTS', 'Kernel', 'kernel', 'normal', 'Width', 50);
                x_values = 0:0.1:3000; y = pdf(pd,x_values); tprob = [tprob, 1 - cdf(pd, 800)];
                COLOR = [56 138 191]/255; ALPHA = 1;
                axt = ax1;
                if Condition==1, LINESTYLE = '--', else LINESTYLE= '-'; end
                XLIM = [0 3000];
            case 2
                SEARCHBOUTS = cell2mat([BOSC.SearchThetaDur]');
                BOUTS = cell2mat([BOSC.SearchThetaDur]');
                [pvalue,h] = ranksum(RESTBOUTS,SEARCHBOUTS);
                pd = fitdist(SEARCHBOUTS', 'Kernel', 'kernel', 'normal', 'Width', 50);
                x_values = 0:0.1:3000; y = pdf(pd,x_values); tprob = [tprob, 1 - cdf(pd, 800)];
                COLOR = [255 8 32]/255; ALPHA = 0.7;
                axt = ax1;
                if Condition==1, LINESTYLE = '--', else LINESTYLE= '-'; end
                XLIM = [0 3000];
            case 3
                RESTBOUTS = cell2mat([BOSC.RestGammaDur]');
                BOUTS = cell2mat([BOSC.RestGammaDur]');
                pd = fitdist(RESTBOUTS', 'Kernel', 'kernel', 'normal', 'Width', 10);
                x_values = 0:0.1:1000; y = pdf(pd,x_values); tprob = [tprob, 1 - cdf(pd, 300)];
                COLOR = [56 138 191]/255; ALPHA = 1;
                axt = ax2;
                if Condition==1, LINESTYLE = '--', else LINESTYLE= '-'; end
                XLIM = [0 1000];
            case 4
                SEARCHBOUTS = cell2mat([BOSC.SearchGammaDur]');
                BOUTS = cell2mat([BOSC.SearchGammaDur]');
                pd = fitdist(SEARCHBOUTS', 'Kernel', 'kernel', 'normal', 'Width', 10);
                x_values = 0:0.1:1000; y = pdf(pd,x_values); tprob = [tprob, 1 - cdf(pd, 300)];
                COLOR = [255 8 32]/255; ALPHA = 0.7;
                axt = ax2;
                if Condition==1, LINESTYLE = '--', else LINESTYLE= '-'; end
                XLIM = [0 1000];
        end
        
        hold(axt, 'on')
        plot(axt, x_values, y, 'Color', COLOR, 'LineWidth', 2, 'LineStyle', LINESTYLE)
        if PLOT == 2
            plot(axt, [800 800], [0 2e-3], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1)
        elseif PLOT == 4
            plot(axt, [300 300], [0 5e-3], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1)
        end
        set(axt, 'LineWidth', 1, 'Box', 'off', 'Tickdir', 'out', 'FontSize', 15)
        ylabel(axt, 'Probability')
        xlabel(axt, 'Time(ms)')
        xlim(axt, XLIM)
%         if PLOT == 1
%             title(axt, pvalue_theta)
%         elseif PLOT == 3
%             title(axt, pvalue_gamma)
%         end
    end
    
end

FigureDirectory = 'D:\HiTGun\0 Main\Figures - refined'
FigureName = 'BoutDurationDistribution'
print(fh, fullfile(FigureDirectory, FigureName),'-dsvg','-r300')
print(fh, fullfile(FigureDirectory, FigureName),'-dpng','-r300')
savefig(fh, fullfile(FigureDirectory, FigureName))
% legend(h)


