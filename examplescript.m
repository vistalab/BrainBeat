clear all
close all

data_path='/biac4/wandell/data/BrainBeat/data/';
subj_name='20140911_1411';
scan_name='3_1_mux8fov4_r1_35s_3mm';
fmri_data=fullfile(data_path,subj_name,scan_name,'7944_3_1.nii.gz');
scan_name='8_1_T1w_1mm_sag';
anat_data=fullfile(data_path,subj_name,scan_name,'/7944_8_1.nii.gz');

ni=niftiRead(fmri_data);

ni=niftiRead(anat_data);

%%

% get slice times:

nslices = BB_get(ni,'mux');
 

mux = 8
nslices = 5
tr = 
mux_slice_acq_order = range(0,nslices,2) + range(1,nslices,2)
mux_slice_acq_time = [float(s)/nslices*tr for s in xrange(nslices)]
unmux_slice_acq_order = [nslices*m+s for m in xrange(mux) for s in mux_slice_acq_order]
unmux_slice_acq_time = mux_slice_acq_time * 3
for acq,slice in zip(unmux_slice_acq_time,unmux_slice_acq_order):
%    print "    slice %02d acquired at time %.3f sec" % (slice+1,acq)
   
end

%%

figure
for k=1:size(ni.data,4)
imagesc(double(ni.data(:,:,30,k))); axis image;
colormap(hot)
pause(.35)

end
