function subs = bb_subs(subj)

switch subj
    case 1
        subs.subj='20140911_1411';
        subs.anat = '8_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '7944_8_1';
        subs.scan{1}='3_1_mux8fov4_r1_35s_3mm';
        subs.scanName{1}='7944_3_1';
        subs.scan{2}='6_1_mux8fov4_r1_25s_4mm';
        subs.scanName{2}='7944_6_1';
    case 2
        subs.freesurferDir = '/sni-storage/wandell/data/anatomy/wandell/20141017/';    % freesurfer data in anatomy folder
        subs.subj = '20141017_1242';    % Date _ Time out of NIMS
        subs.anat = '9_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '8202_9_1';
        subs.veno = '3_1_2D_MRV';   % Venogram data
        subs.venoName = '8202_3_1';
        subs.scan{1} ='5_1_mux8fov4_r1_25s_4mm';  
        subs.scanName{1}='8202_5_1';
        subs.scanFA{1}=34;
        subs.scan{2} ='6_1_mux8fov4_r1_25s_4mmFA25';  
        subs.scanName{2}='8202_6_1';
        subs.scanFA{2}=25;
        subs.scan{3} ='7_1_mux8fov4_r1_25s_4mmFA48';  
        subs.scanName{3}='8202_7_1';
        subs.scanFA{3}=48;
    case 3
        subs.freesurferDir = '/sni-storage/wandell/data/anatomy/rosemary2015/';
        subs.subj = '20150106_1621';    % Date _ Time out of NIMS
        subs.anat = '3_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '8603_3_1';
        subs.veno = '14_1_2D_MRV';   % Venogram data
        subs.venoName = '8603_14_1';
        subs.scan{1} ='10_1_mux8_25s_4mmFA36';  % subject awake, watching movie 
        subs.scanName{1}='8603_10_1';
        subs.scanFA{1}=36;
        subs.scan{2} ='11_1_mux8_25s_4mmFA20';  % subject awake, watching movie 
        subs.scanName{2}='8603_11_1';
        subs.scanFA{2}=20;
        subs.scan{3} ='12_1_mux8_25s_4mmFA48';  % subject awake, watching movie 
        subs.scanName{3}='8603_12_1';
        subs.scanFA{3}=48;
        subs.scan{4} ='16_1_mux8_25s_4mmFA20';  % subject asleep (subjective report)
        subs.scanName{4}='8603_16_1'; 
        subs.scanFA{4}=20;
        subs.scan{5} ='23_1_mux8_25s_4mmFA20';  % subject awake (subjective report)
        subs.scanName{5}='8603_23_1'; 
        subs.scanFA{5}=20;
end