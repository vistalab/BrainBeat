function val = BB_get(ni,param,varargin)
% Get value from BB structure
%
%    val = BB_get(ni, param, varargin)
%
% param list:              - arguments:
%
%  'tr'
%  'slice duration'  - duration of 1 slice within a tr
%   'n slices        - number of slices within 1 tr
%   'mux'            - the number of simultaneous slices (mux factor)
%   'timing'         - timing for every slice in ms
%   'physio'         - physio data: ppg (pulse) and resp (resp)
%   'ppg_peaks'      - returns the peaks in the ppg signal in seconds
%                      (ppg = photoplethethysmogram)
%   'ppg_response_function' - heart-rate locked brain impulse response for
%                             every voxel and gets the time (t) 
%                           - varargin{1} = slice number is optional
%                             (default all slices, but that takes a couple of minutes)
%   'ppg_correlation - correlation with ppg signal size voxels X voxels X
%                      slices
%                    - varargin{1} = slice number is optional
%                      (default all slices, but that takes a couple of minutes)
% Examples: 
%    BB_get(ni,'physio')
%     
%
% Written by Aviv and Dora, Copyright Vistasoft Team 2014


%% remove spaces and upper case
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
        % p = bbGet(ni,'physio')
        % Read the physio and respiratory data specified in the NI file

        % We first check to see if there is an unziped/tar file with the
        % name
        
        % Should add physio file handling to niftiGet(ni,'physio ....');
        [fPath,fName] = fileparts(ni.fname);
        [~, fName]    = fileparts(fName);       
        UnTarpPhysio_dir = fullfile(fPath,[fName '_physio_Unzip']);
        
        if ~exist(UnTarpPhysio_dir,'dir')
            % Make the untar'd gunzip'd physio directory
            physio_name = fullfile(fPath,[fName '_physio.tgz']);
            if ~exist(physio_name,'file')
                error('No physio tgz file %s\n',physio_name);
            end
            physioTar = gunzip(physio_name);    % This is a cell of the tar file name
            untar(physioTar{1},UnTarpPhysio_dir);
            
            % Eliminate the unwanted tar file
            delete(physioTar{1});
        end
        
        % physio_dir  = fullfile(file_path,[file_name '_physio']);
        
        % So the directory now exists.  Get the physio data.  One of them
        % is the PPG and the other is RESP
        physioData = dir(fullfile(UnTarpPhysio_dir,[fName,'_physio'],'*Data*'));
       
        % This the Matlab structure we return.  The first one is PPG and
        % the second one is RESP.  Maybe we should structure it as
        %  physio.ppg and physio.resp
        %
        physio_output = [];
        for k=1:length(physioData)
            if isequal(physioData(k).name(1:6),'PPGDat');
                physio_output(1).name = 'PPG';
                physio_output(1).CNIsrate = 100; % here, we are hardcoding the CNI sampling rate for ECG
                physio_output(1).rawdata = load(fullfile(UnTarpPhysio_dir,[fName,'_physio'],physioData(k).name));
            elseif isequal(physioData(k).name(1:6),'RESPDa')
                physio_output(2).name = 'RESP';
                physio_output(2).CNIsrate = 25; % here, we are hardcoding the CNI sampling rate for respiration
                physio_output(2).rawdata = load(fullfile(UnTarpPhysio_dir,[fName,'_physio'],physioData(k).name));
            end
        end
        % now we are going to time-lock it to the scan, and include a
        % parameter for scan-onset
        scan_duration=ni.dim(4)*ni.pixdim(4);% in sec
        for k=1:length(physio_output)
            if isfield(physio_output(k),'CNIsrate') % if there is a sampling rate
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
        ppg_ind=1; % 1 for ppg, 2 for resp

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
        
    case{'ppg_correlation'}
        % gets the correlation with PPG for every voxel
        % adding a slice number is optional in varargin{1}
        val = BB_correlate2physio(ni,varargin{1});
        
    otherwise
        error('unknown BB parameter');
        
end