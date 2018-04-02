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
        subs.freesurferDir = '/sni-storage/wandell/data/anatomy/';
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
    case 4
        subs.freesurferDir = '/sni-storage/wandell/data/anatomy/';
        subs.subj = '20180319_1232';    % Date _ Time out of NIMS
        subs.anat = '22_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '17178_22_1';
        subs.veno = '21_1_2D_MRV';   % Venogram data
        % mux 4 mm scans
        subs.venoName = '17178_21_1';
        subs.scan{1} ='5_1_mux8_25s_4mmFA20';  % subject awake, eyes closed
        subs.scanName{1}='17178_5_1';
        subs.scanFA{1}=20;
        subs.scan{2} ='4_1_mux8_25s_4mmFA34';  % subject awake, eyes closed
        subs.scanName{2}='17178_4_1';
        subs.scanFA{2}=34;
        subs.scan{3} ='6_1_mux8_25s_4mmFA48';  % subject awake, eyes closed
        subs.scanName{3}='17178_6_1';
        subs.scanFA{3}=48;
        
        % reverse order 
        subs.scan{4} ='8_1_mux8_25s_4mmFA20';  % subject awake, eyes closed
        subs.scanName{4}='17178_8_1';
        subs.scanFA{4}=20;
        subs.scan{5} ='7_1_mux8_25s_4mmFA34';  % subject awake, eyes closed
        subs.scanName{5}='17178_7_1';
        subs.scanFA{5}=34;
        subs.scan{6} ='9_1_mux8_25s_4mmFA48';  % subject awake, eyes closed
        subs.scanName{6}='17178_9_1';
        subs.scanFA{6}=48;
        
        % 3 mm voxel
        subs.scan{7} ='17_1_mux8_25s_3mmFA20';  % subject awake, eyes closed
        subs.scanName{7}='17178_17_1';
        subs.scanFA{7}=20;
        subs.scan{8} ='19_1_mux8_25s_3mmFA34';  % subject awake, eyes closed
        subs.scanName{8}='17178_19_1';
        subs.scanFA{8}=48;
        subs.scan{9} ='18_1_mux8_25s_3mmFA48';  % subject awake, eyes closed
        subs.scanName{9}='17178_18_1';
        subs.scanFA{9}=34;
        
        % GE hyperband, no physio (order is top: 2 bottom: 40)
        subs.scan{10} ='15_1_HB8_025s_4mm_FA20';  % subject awake, eyes closed
        subs.scanName{10}='17178_15_1';
        subs.scanFA{10}=20;
        subs.scan{11} ='13_1_HB8_025s_4mm_FA34';  % subject awake, eyes closed
        subs.scanName{11}='17178_13_1';
        subs.scanFA{11}=48;
        subs.scan{12} ='11_1_HB8_025s_4mm_FA48';  % subject awake, eyes closed
        subs.scanName{12}='17178_11_1';
        subs.scanFA{12}=34;
        
end