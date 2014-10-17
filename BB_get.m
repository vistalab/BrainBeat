function val = BB_get(ni,param,varargin)
% Get value from BB structure
%
% val = BB_get(ni, param, varargin)
%
% param list:              - arguments:
%
% 'tr'
% 'sliceduration'
        % duration of 1 slice within a tr
% nslices
        % number of slices within 1 tr
% mux
        % the mux factor
% 'timing'
        % timing for every slice in ms
% 'physio'
        % physio data: ppg (pulse) and resp
% 'ppg_peaks' 
        % returns the peaks in the ppg signal in seconds
% 'ppg_response_function'
        % gets the heart-rate locked brain response for every voxel
        % and gets the time (t)
        % varargin{1}: adding a slice number is optional
%
%
% example: BB_get(ni,'physio')
%
% Wandell Copyright Vistasoft Team, 2013
% Written by Aviv and Dora 2014


% remove spaces and upper case
param = mrvParamFormat(param);

%% Get the requested parameter
switch(param)
    case{'tr'}
        val = ni.pixdim(4);
    case('sliceduration')
        val=ni.slice_duration;
    case{'nslices'}
        val = round(ni.pixdim(4)./ni.slice_duration);
    case{'mux'}
        val = round(ni.dim(3)./BB_get(ni,'nslices'));
    case{'timing'}
        nslices=BB_get(ni,'nslices');
        mux=BB_get(ni,'mux');
        tr=ni.pixdim(4);
        mux_slice_acq_order = [[1:2:round(nslices)] [2:2:round(nslices)]];
        
        clear mux_slice_acq_time
        for s=0:nslices-1
            mux_slice_acq_time(s+1) = s/nslices*tr;
        end
        
        clear unmux_slice_acq_order
        ii=0;
        for m=0:mux-1
        for s=mux_slice_acq_order
            ii=ii+1;
            unmux_slice_acq_order(ii)=nslices*m+s;
        end
        end
        
        unmux_slice_acq_time=[];
        for k=1:mux
            unmux_slice_acq_time = [unmux_slice_acq_time mux_slice_acq_time];
        end
        % we have timing per volume (unmux_slice_acq_order unmux_slice_acq_time)
        
        % now we create a matrix with columns for nr of slices and rows for
        % nr of tr
        timing=zeros(ni.dim(3),ni.dim(4));
        for k=1:size(timing,2)
            timing(:,k)=unmux_slice_acq_time+((k-1)*tr);
        end
        val=timing;
    case{'physio'}
        [file_path,file_name]=fileparts(ni.fname);
        [~,file_name]=fileparts(file_name);
        physio_name=fullfile(file_path,[file_name '_physio.tgz']);
        physio_dir=fullfile(file_path,[file_name '_physio']);
        UnTarpPhysio_dir=fullfile(file_path,[file_name '_physio_Unzip']);
        if ~exist(physio_dir,'dir')
            physio=gunzip(physio_name);
            physio=untar(physio{1},UnTarpPhysio_dir);
        else
            physioTmp1=dir(UnTarpPhysio_dir);
            physioTmp2=dir(fullfile(UnTarpPhysio_dir,physioTmp1(3).name));
            physio={physioTmp2(3:end).name};
            for k=1:length(physio)
                physio{k}=fullfile(UnTarpPhysio_dir,physioTmp1(3).name,physio{k});
            end
        end
        physio_output=[];
        for k=1:length(physio)
            [~, d]=fileparts(physio{k});
            physio_output(k).name=d(1:6);
            physio_output(k).rawdata=load(physio{k});
            if isequal(physio_output(k).name,'PPGDat');
                physio_output(k).CNIsrate=100;
                % here, we are hardcoding the CNI sampling rate for
                % ECG
            elseif isequal(physio_output(k).name,'RESPDa')
                physio_output(k).CNIsrate=25;
                % here, we are hardcoding the CNI sampling rate for
                % respiration
            end
        end
        % now we are going to time-lock it to the scan, and include a
        % parameter for scan-onset
        scan_duration=ni.dim(4)*ni.pixdim(4);% in sec
        for k=1:length(physio_output)
            if isfield(physio_output(1),'CNIsrate') % if there is a sampling rate
                % chop of the beginning
                physio_output(k).data=physio_output(k).rawdata(end-round(scan_duration*physio_output(k).CNIsrate)+1:end);
                % add scan onset:
                physio_output(k).scan_onset=zeros(size(physio_output(k).data));
                for m=1:round(ni.dim(4))
                    physio_output(k).scan_onset(round(ni.pixdim(4)*(m-1)*physio_output(k).CNIsrate+1))=1;
                end
            end
        end
        
        val=physio_output;
        
    case{'ppg_peaks'}
        % returns the peaks in the ppg signal in seconds
        
        physio=BB_get(ni,'physio');
        ppg_ind=0;
        for k=1:length(physio)
            if isequal(physio(k).name,'PPGDat')
                ppg_ind=k;
            end
        end
        % first get a rough estimate of onsets for the heartbeat
        signal=physio(ppg_ind).data;
        phys_srate=physio(ppg_ind).CNIsrate;
        % set minimum inter-heartbeat interval
        PPG_int=0.7;
        
        % we already got our value, but keep this code if necessary for
        % later:
        % low-pass filter % not necessary, keep in for lower quality data?
        band=5;
        Rp=3; Rs=60; % third order Butterworth
        high_p=band(1)*2/phys_srate;
        high_s=(band(1)+20)*2/phys_srate;
        [n_band,wn_band]=buttord(high_p,high_s,Rp,Rs);
        [bf_b,bf_a]=butter(n_band,wn_band,'low');
        band_sig=filtfilt(bf_b,bf_a,signal);
        % detect peaks:
        [pks1,locs1]=findpeaks(band_sig,'minpeakdistance',PPG_int*phys_srate);
        % now move back to seconds and return value
        val=locs1/phys_srate; % ppg onsets

        % add some plots so we can check stuff:
        % make epochs and plot:
        epoch_pre=.5*phys_srate;
        epoch_post=4*phys_srate;
        hb_epochs1=NaN(length(pks1),length(-epoch_pre+1:epoch_post));

        for k=1:length(pks1)
            if locs1(k)+epoch_post<length(signal) && locs1(k)-epoch_pre>0
                hb_epochs1(k,:)=signal(locs1(k)-epoch_pre+1:locs1(k)+epoch_post);
            end
        end

        % plot averaged PPG response 
        figure('Position',[0 0 800 700])
        subplot(2,2,1),hold on, plot(signal,'k'), plot(locs1,pks1,'r.')
        title('PPG signal and detected peaks')
        subplot(2,2,2),hold on, plot(signal,'k'), plot(locs1,pks1,'r.')
        xlim([0 phys_srate*30])
        title('PPG signal and detected peaks - zoom')
        subplot(2,1,2),hold on
        t=[-epoch_pre+1:epoch_post]/phys_srate;
        plot(t,hb_epochs1,'k'), plot(t,mean(hb_epochs1),'r')
        title('PPG epochs')

    case{'ppg_response_function'}
        % gets the heart-rate locked brain response for every voxel
        % adding a slice number is optional in varargin{1}
        [response_matrix,t] = BB_response2physio(ni,varargin{1});
        val.response_matrix=response_matrix;
        val.t=t;
    otherwise
        error('Uknown afq parameter');
        
end