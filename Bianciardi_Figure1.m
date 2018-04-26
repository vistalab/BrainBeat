
% Implements equations from Bianciardi et al., 2015
% Reproduces figure 1

% MRI acquisition settings:
FA_deg = [10 30 45 65 90]; % flip angle in degrees to plot
FA_all = (FA_deg*2*pi)/360; % flip angle in radians
TE = 0;
TR = 0.033;

V = 0:0.01:11; % cm/s
Vc = 3.6; % ST/TR (ST = slice thickness) speed at which there is complete inflow of new spins

% 7T tissue parameters [blood WM GM CSF] from Bianciardi
T1 = [2.5 1.2 2.1 4.4]; % seconds longitudinal recovery rate
T2s = 1; % displays are for TE = 0, T2s can be anything as equations have TE/T2s;

% 3T tissue parameters [blood WM GM CSF] from Bianciardi
% T1 = [1.7 0.8 1.3 4.4]; % seconds longitudinal recovery rate

tt = 1; % tissue type: 1=blood 2=WM 3=GM 4=CSF

%%

% figure,hold on
for FA = FA_all
    
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
    plot(V,Mxy)
end
legend(int2str(FA_deg'))
