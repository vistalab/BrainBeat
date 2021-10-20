% Implements equations from Bianciardi et al., 2015
% Play with parameters: 
% 1) test our Stanford scan parameters with different flip angles

% since whole brain --> Ernst angle equations would be better

%%

% MRI acquisition settings:
FA_deg = [5 20 36 48 90]; % flip angle in degrees to plot
FA_all = (FA_deg*2*pi)/360; % flip angle in radians
TE = 0.0116;
% [WM GM CSF arterial_blood veinous_blood]
% T2s_all = [0.021 0.015 2 .049 0.015]; % T2s from Weiskopf 2013 Frontiers in Neuroscience, with Viesmann
T2s_all = [0.021 0.015 .300 .049 0.015]; % T2s from Weiskopf 2013 Frontiers in Neuroscience, with Viesmann
TR = 0.250;

V = 0:0.01:11; % cm/s
Vc = 0.4/TR; % cm/s ST/TR (ST = slice thickness) speed at which there is complete inflow of new spins

% Tissue parameters [WM GM CSF blood] from Bianciardi at 7T
% T1 = [1.2 2.1 4.4 1.7]; % seconds longitudinal recovery rate

% % Tissue parameters [WM GM CSF blood] from Bianciardi at 3T
T1 = [0.8 1.3 4.4 1.7]; % seconds longitudinal recovery rate
% T1=[0.9  1.67  4]; % old: previous used longitudinal recovery rate 

PD = [0.7 0.85 1 1]; % this is not used in Bianciardi, from Mezer et al, blood is a guess

%% plot for V==0, V<Vc, V>=Vc

FA_deg = [1:90]; % flip angle in degrees to plot
FA_all = (FA_deg*2*pi)/360; % flip angle in radians

% tt_colors = {[0 .3 1],[.5 .5 .5],'c'};
% tt_colors = {'r','b','g'};

Mxy1 = NaN(3,length(FA_all));
Mxy2 = NaN(3,length(FA_all));
Mxy3 = NaN(3,length(FA_all));

for FA_ind = 1:length(FA_all)
    FA = FA_all(FA_ind);
    for tt = 1:4 % [1:WM 2:GM 3:CSF 4:blood]

        T2s = T2s_all(tt);

        G = 1; % a gain factor
        M0 = 1 * PD(tt);

        q = exp(-TR/T1(tt))*cos(FA);

        % regime 1 v==0:
        kk = find(V==0);
        MzV0 = (M0 * (1-exp(-TR/T1(tt))))/(1-q);
        Mxy1(tt,FA_ind) = G * sin(FA) * MzV0 *...
            exp(-TE/T2s);
         
        % regime 2 v<vc:
        kk = find(V>.4,1);% find(V>0 & V<Vc) 
        MzV0 = (M0 * (1-exp(-TR/T1(tt))))/(1-q);
        Mxy2(tt,FA_ind) = G * sin(FA) * (MzV0 + (M0-MzV0)*(1-q.^(Vc/V(kk)))/((Vc/V(kk))*(1-q)) ) *...
            exp(-TE/T2s);

        % regime 3 v>vc:
        kk = find(V>=Vc,1);
        Mxy3(tt,FA_ind) = G * sin(FA) * M0 *...
            exp(-TE/T2s);

    end
end

figure('Position',[0 300 800 250])
subplot(1,3,1), plot(Mxy1')
ylim([0 1])
legend({'WM','GM','CSF','blood'})
xlabel('flip angle'),ylabel('Mxy')
title('v==0')


subplot(1,3,2), plot(Mxy2')
ylim([0 1])
legend({'WM','GM','CSF'})
xlabel('flip angle'),ylabel('Mxy')
title('v<vc')

subplot(1,3,3), plot(Mxy3')
ylim([0 1])
legend({'WM','GM','CSF'})
xlabel('flip angle'),ylabel('Mxy')
title('v>=vc')

%%
FA_deg = [5 20 36 48 90]; % flip angle in degrees to plot
FA_all = (FA_deg*2*pi)/360; % flip angle in radians

figure('Position',[0 300 800 400])
tt_colors = {[0 .3 1],[.5 .5 .5],'c','r'};
% tt_colors = {'r','b','g'};
for FA_ind = 1:length(FA_all)
    FA = FA_all(FA_ind);
    subplot(2,length(FA_all),FA_ind),hold on
    title(['FA = ' int2str(FA_deg(FA_ind))])
    for tt = 1:4 % [1:WM 2:GM 3:CSF 4:blood]

        T2s = T2s_all(tt);

        G = 1; % a gain factor
        M0 = 1 * PD(tt);

        q = exp(-TR/T1(tt))*cos(FA);

        Mxy = NaN(size(V));

        % regime 1 v==0:
        for kk = find(V==0)   
            MzV0 = (M0 * (1-exp(-TR/T1(tt))))/(1-q);
            Mxy(kk) = G * sin(FA) * MzV0 *...
                exp(-TE/T2s);
        end

        % regime 2 v<vc:
        for kk = find(V>0 & V<Vc) 
            MzV0 = (M0 * (1-exp(-TR/T1(tt))))/(1-q);
            Mxy(kk) = G * sin(FA) * (MzV0 + (M0-MzV0)*(1-q.^(Vc/V(kk)))/((Vc/V(kk))*(1-q)) ) *...
                exp(-TE/T2s);
        end

        % regime 3 v>vc:
        for kk = find(V>=Vc) 
            Mxy(kk) = G * sin(FA) * M0 *...
                exp(-TE/T2s);
        end

        subplot(2,length(FA_all),FA_ind),hold on
        plot(V,Mxy,'LineWidth',2,'Color',tt_colors{tt})
        xlabel('v')
        ylabel('Mxy')
        xlim([V(1) V(end)]),ylim([0 1])

        subplot(2,length(FA_all),length(FA_all)+FA_ind),hold on
        plot(V,Mxy,'LineWidth',2,'Color',tt_colors{tt})
        xlabel('v')
        ylabel('Mxy')
        xlim([V(1) .8]),ylim([0 1])
    end
end

subplot(2,length(FA_all),1)
legend({'WM','GM','CSF','blood'})
