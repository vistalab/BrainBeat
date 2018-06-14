 function bb_roi = bb_subs_rois(subj)

switch subj
    case 1

    case 2
    % Figure traces medial regions:
    bb_roi(1).curPos = [-1 21 -86]; 
    bb_roi(1).voxelLabel = 'Basilar1';
    bb_roi(2).curPos = [-1 75 -30]; 
    bb_roi(2).voxelLabel = 'ACA';
    bb_roi(3).curPos = [-1 -46 -46];%another site is [-1 -50 -49]; 
    bb_roi(3).voxelLabel = 'MOA';
    bb_roi(4).curPos = [-1 -63 -17];
    bb_roi(4).voxelLabel = 'SagitalSinus1';
    bb_roi(5).curPos = [-1 -41 -65];
    bb_roi(5).voxelLabel = 'StraightSinus';
    bb_roi(6).curPos = [-1 -8 21];
    bb_roi(6).voxelLabel = 'Subarachnoid';
    bb_roi(7).curPos = [-1 43 -22];%[-1 38 -22] - Dual peaked?
    bb_roi(7).voxelLabel = 'LateralVentricle1';
    bb_roi(8).curPos = [-1 -30 -15];%[-1 38 -22] - Dual peaked?
    bb_roi(8).voxelLabel = 'PosteriorGM';

    case 3
    % Figure traces medial regions:
    bb_roi(1).curPos = [1 4 -23];
    bb_roi(1).voxelLabel = 'Basilar1';

    case 4
    % Figure traces medial regions:
    bb_roi(1).curPos = [0 19 39];
    bb_roi(1).voxelLabel = 'LateralVentricle1';

end