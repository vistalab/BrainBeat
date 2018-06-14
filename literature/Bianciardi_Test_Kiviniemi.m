% Implements equations from Bianciardi et al., 2015
% Play with parameters: 
% 1) test Kiviniemi et al., 2016 parameters

% MRI acquisition settings:
FA_deg = [25]; % flip angle in degrees to plot
FA = (FA_deg*2*pi)/360; % flip angle in radians
TE = 0.036;
T2s_all = [0.021 0.015  2]; % T2s from Weiskopf 2013 Frontiers in Neuroscience
TR = 0.100;

V = 0:0.01:11; % cm/s
Vc = 0.3/TR; % ST/TR (ST = slice thickness) speed at which there is complete inflow of new spins

% Tissue parameters [WM GM CSF] from Bianciardi at 7T
% T1 = [1.2 2.1 4.4]; % seconds longitudinal recovery rate

% % Tissue parameters [WM GM CSF] from Bianciardi at 3T
T1 = [0.8 1.3 4.4]; % seconds longitudinal recovery rate

%%
figure,hold on
tt_colors = {'y',[.5 .5 .5],'c'};
for tt = 1:3 % [1:WM 2:GM 3:CSF]

T2s = T2s_all(tt);

G = 1; % a gain factor
M0 = 1;

q = exp(-TR/T1(tt))*cos(FA);

Mxy = NaN(size(V));

% regime 1 v==0:
for kk = find(V==0)   
    MzV0 = (M0 * (1-exp(-TR/T1(tt))))/(1-q);
    Mxy(kk) = G * sin(FA) * MzV0 * exp(-TE/T2s);
end

% regime 2 v<vc:
for kk = find(V>0 & V<Vc) 
    MzV0 = (M0 * (1-exp(-TR/T1(tt))))/(1-q);
    Mxy(kk) = G * sin(FA) * (MzV0 + (M0-MzV0)*(1-q.^(Vc/V(kk)))/((Vc/V(kk))*(1-q)) )* exp(-TE/T2s);
end

% regime 3 v>vc:
for kk = find(V>=Vc) 
    Mxy(kk) = G * sin(FA) * M0 * exp(-TE/T2s);
end

plot(V,Mxy,'LineWidth',2,'Color',tt_colors{tt})
xlabel('v (cm/s)')
ylabel('Mxy')
end