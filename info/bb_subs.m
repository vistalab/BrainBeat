function subs = bb_subs(subj)

switch subj
    case 1
        subs.subj='20140911_1411';
        subs.age = 24;
        subs.gender = 'f';
        subs.anat = '8_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '7944_8_1';
        subs.scan{1}='3_1_mux8fov4_r1_35s_3mm';
        subs.scanName{1}='7944_3_1';
        subs.scan{2}='6_1_mux8fov4_r1_25s_4mm';
        subs.scanName{2}='7944_6_1';
    case 2
        subs.subj = '20141017_1242';    % Date _ Time out of NIMS
        subs.age = 63;
        subs.gender = 'm';
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
        subs.subj = '20150107_1621';    % Date _ Time out of NIMS
        subs.age = 24;
        subs.gender = 'f';
        subs.anat = '3_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '8603_3_1';
        subs.veno = '14_1_2D_MRV';   % Venogram data
        subs.venoName = '8603_14_1';
        subs.scan{1} = '10_1_mux8_25s_4mmFA36';  % subject awake, watching movie 
        subs.scanName{1} = '8603_10_1';
        subs.scanFA{1} = 36;
        subs.scan{2} = '11_1_mux8_25s_4mmFA20';  % subject awake, watching movie 
        subs.scanName{2} = '8603_11_1';
        subs.scanFA{2} = 20;
        subs.scan{3} = '12_1_mux8_25s_4mmFA48';  % subject awake, watching movie 
        subs.scanName{3} = '8603_12_1';
        subs.scanFA{3} = 48;
        subs.scan{4} = '16_1_mux8_25s_4mmFA20';  % subject asleep (subjective report)
        subs.scanName{4} = '8603_16_1'; 
        subs.scanFA{4} = 20;
        subs.scan{5} ='23_1_mux8_25s_4mmFA20';  % subject awake (subjective report)
        subs.scanName{5} = '8603_23_1'; 
        subs.scanFA{5} = 20;
    case 4
        subs.subj = '20180319_1232';    % Date _ Time out of NIMS
        subs.age = 25;
        subs.gender = 'f';
        subs.anat = '22_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '17178_22_1';
        subs.veno = '21_1_2D_MRV';   % Venogram data
        % mux 4 mm scans
        subs.venoName = '17178_21_1';
        subs.scan{1} = '5_1_mux8_25s_4mmFA20';  % subject awake, eyes closed
        subs.scanName{1} = '17178_5_1';
        subs.scanFA{1} = 20;
        subs.scan{2} ='4_1_mux8_25s_4mmFA34';  % subject awake, eyes closed
        subs.scanName{2} = '17178_4_1';
        subs.scanFA{2} = 34;
        subs.scan{3} = '6_1_mux8_25s_4mmFA48';  % subject awake, eyes closed
        subs.scanName{3} = '17178_6_1';
        subs.scanFA{3} = 48;
        
        % reverse order 
        subs.scan{4} = '8_1_mux8_25s_4mmFA20';  % subject awake, eyes closed
        subs.scanName{4} = '17178_8_1';
        subs.scanFA{4} = 20;
        subs.scan{5} = '7_1_mux8_25s_4mmFA34';  % subject awake, eyes closed
        subs.scanName{5} = '17178_7_1';
        subs.scanFA{5} = 34;
        subs.scan{6} = '9_1_mux8_25s_4mmFA48';  % subject awake, eyes closed
        subs.scanName{6} = '17178_9_1';
        subs.scanFA{6} = 48;
        
        % 3 mm voxel
        subs.scan{7} ='17_1_mux8_25s_3mmFA20';  % subject awake, eyes closed
        subs.scanName{7}='17178_17_1';
        subs.scanFA{7}=20;
        subs.scan{8} ='19_1_mux8_25s_3mmFA34';  % subject awake, eyes closed
        subs.scanName{8}='17178_19_1';
        subs.scanFA{8}=34;
        subs.scan{9} ='18_1_mux8_25s_3mmFA48';  % subject awake, eyes closed
        subs.scanName{9}='17178_18_1';
        subs.scanFA{9}=48;
        
        % GE hyperband, no physio (order is top: 2 bottom: 40)
        subs.scan{10} ='15_1_HB8_025s_4mm_FA20';  % subject awake, eyes closed
        subs.scanName{10}='17178_15_1';
        subs.scanFA{10}=20;
        subs.scan{11} ='13_1_HB8_025s_4mm_FA34';  % subject awake, eyes closed
        subs.scanName{11}='17178_13_1';
        subs.scanFA{11}=34;
        subs.scan{12} ='11_1_HB8_025s_4mm_FA48';  % subject awake, eyes closed
        subs.scanName{12}='17178_11_1';
        subs.scanFA{12}=48;
    case 5
        % people at scan: Hua Wu and Dora Hermes
        subs.age = 30;
        subs.gender = 'f';
        subs.subj = '20180529_17677';    % Date _ Time out of NIMS
        subs.anat = '11_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '17677_11_1';
        subs.veno = '12_1_2D_MRV';   % Venogram data
        subs.venoName = '17677_12_1';
        
        % mux 4 mm scans
        subs.scan{1} = '4_1_mux8_25s_4mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{1} = '17677_4_1';
        subs.scanFA{1} = 48;
        
        % mux 4 mm scans reverse order 
        subs.scan{2} = '5_1_mux8_25s_4mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{2} = '17677_5_1';
        subs.scanFA{2} = 48;
        
        % 3 mm voxel
        subs.scan{3} = '9_1_mux8_25s_3mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{3} = '17677_9_1';
        subs.scanFA{3} = 48;
        
        % mux 4 mm scans, double echo, echo 1
        subs.scan{4} = '6_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{4} = '17677_6_1_te0008';
        subs.scanFA{4} = 48;
        
        % mux 4 mm scans, double echo, echo 2
        subs.scan{5} = '6_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{5} = '17677_6_1_te0028';
        subs.scanFA{5} = 48;

        % mux 4 mm scans, double echo, echo 1
        subs.scan{6} = '7_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{6} = '17677_7_1_te0008';
        subs.scanFA{6} = 48;
        
        % mux 4 mm scans, double echo, echo 2
        subs.scan{7} = '7_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{7} = '17677_7_1_te0028';
        subs.scanFA{7} = 48;

        % GE hyperband, did not run, don't know why 
                
    case 6
        subs.age = 30;
        subs.gender = 'm';
        subs.subj = '20180607_17734';    % Date _ Time out of NIMS
        subs.anat = '5_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '17734_5_1';
        subs.veno = '6_1_2D_MRV';   % Venogram data
        subs.venoName = '17734_6_1';
        
        % mux 4 mm scans
        subs.scan{1} = '4_1_mux8_25s_4mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{1} = '17734_4_1';
        subs.scanFA{1} = 48;
        
        % mux 4 mm scans ? reverse order ?
        subs.scan{2} = '10_1_mux8_25s_4mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{2} = '17734_10_1';
        subs.scanFA{2} = 48;
        
        % 3 mm voxel
        subs.scan{3} = '11_1_mux8_25s_3mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{3} = '17734_11_1';
        subs.scanFA{3} = 48;
    
        % mux 4 mm scans, double echo, echo 1
        subs.scan{4} = '8_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{4} = '17734_8_1_te0008';
        subs.scanFA{4} = 48;
        
        % mux 4 mm scans, double echo, echo 2
        subs.scan{5} = '8_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{5} = '17734_8_1_te0028';
        subs.scanFA{5} = 48;

        % mux 4 mm scans, double echo, echo 1
        subs.scan{6} = '9_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{6} = '17734_9_1_te0008';
        subs.scanFA{6} = 48;
        
        % mux 4 mm scans, double echo, echo 2
        subs.scan{7} = '9_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{7} = '17734_9_1_te0028';
        subs.scanFA{7} = 48;
        
    case 7
        subs.age = 66;
        subs.gender = 'm';
        subs.subj = '20180619_17796';    % Date _ Time out of NIMS
        subs.anat = '11_1_T1w_1mm_sag';   % Anatomical data
        subs.anatName = '17796_11_1';
        subs.veno = '12_1_2D_MRV';   % Venogram data
        subs.venoName = '17796_12_1';
        
        % mux 4 mm scans
        subs.scan{1} = '5_1_mux8_25s_4mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{1} = '17796_5_1';
        subs.scanFA{1} = 48;
        
        % mux 4 mm scans ? reverse order ?
        subs.scan{2} = '6_1_mux8_25s_4mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{2} = '17796_6_1';
        subs.scanFA{2} = 48;
        
        % 3 mm voxel
        subs.scan{3} = '9_1_mux8_25s_3mmFA48';  % subject awake, eyes open, watching a movie
        subs.scanName{3} = '17796_9_1';
        subs.scanFA{3} = 48;
    
        % mux 4 mm scans, double echo, echo 1
        subs.scan{4} = '7_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{4} = '17796_7_1_te0008';
        subs.scanFA{4} = 48;
        
        % mux 4 mm scans, double echo, echo 2
        subs.scan{5} = '7_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{5} = '17796_7_1_te0028';
        subs.scanFA{5} = 48;

        % mux 4 mm scans, double echo, echo 1
        subs.scan{6} = '8_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{6} = '17796_8_1_te0008';
        subs.scanFA{6} = 48;
        
        % mux 4 mm scans, double echo, echo 2
        subs.scan{7} = '8_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{7} = '17796_8_1_te0028';
        subs.scanFA{7} = 48;

        % mux 4 mm scans
        subs.scan{8} = '16_1_mux8_25s_4mmFA48'; % subject awake/asleep?
        subs.scanName{8} = '17796_16_1';
        subs.scanFA{8} = 48;
        
        % mux 4 mm scans, double echo, echo 1
        subs.scan{9} = '17_1_fMRI_mux8_me2_4mm'; % subject awake, eyes open, watching a movie
        subs.scanName{9} = '17796_17_1_te0008';
        subs.scanFA{9} = 48;

        % mux 4 mm scans, double echo, echo 2
        subs.scan{10} = '17_1_fMRI_mux8_me2_4mm'; % subject awake/asleep?
        subs.scanName{10} = '17796_17_1_te0028';
        subs.scanFA{10} = 48;
end


