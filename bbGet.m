function val = bbGet(ni,param,varargin)
% Get value from BB (Brain beat) structure
%
%    val = bbGet(ni, param, varargin)
%
% The brain beat project tries to measure fMRI signals at a high rate and
% then understand which aspects of the response vary with the heart beat
% and sometimes with the respiratory cycle.
%
% The heart beat is measured with a Plethysmograph based on light and worn
% by the subject on their index finger.  So GE calls it a Photo
% Plethysmograph (PPG).
% (http://en.wikipedia.org/wiki/Plethysmograph) 
%
% We also have subjects wear a breathing belt to record their respiration.
%
%
% Inputs:
%   ni    - a nifti 4D file with the imaging data over time (niftiRead)
%   param - See choices below
%   
% param list:              - arguments:
%
%   'tr'             - Repetition time
%   'slice duration' - duration of one slice within a tr
%   'n slices        - number of slices within one tr
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
%   bbGet(ni,'physio')
%   bbGet(ni,'tr')
%   bbGet(ni,'slice duration')
%   bbGet(ni,'super slices')
%   bbGet(ni,'total slices')
%   bbGet(ni,'simultaneous slices')
%   bbGet(ni,'slice acquisition order')
%   bbGet(ni,'n volumes')
%   T = bbGet(ni,'timing'); 
%   mrvNewGraphWin; plot(T(:,1),'-o'); hold on; plot(T(:,2),'-x')
%
%   p = bbGet(ni,'physio');
%   
%
% AM/DH/BW Vistasoft Team , 2014

%% Programming todo
% Let's move the physio part out and re-write a set of functions called
% physioCreate/Get/Set
%    p = physioCreate(filename)
% if filename is empty then we return a default physio structure
% The physio structure should have physio.resp and physio.ppg as two main
% slots.
%
% The physio stuff in here can be moved out.  This routine is good as
% taking a nifti structure in and returning the SMS parameters.  We should
% probably try to put a label into the NIFTI structure somewhere that
% identifies it as a SMS (brain beat) type of structure.

%% remove spaces and upper case from para
param = mrvParamFormat(param);

%% Get the requested parameter
switch(param)
    case{'tr'}
        % The fourth dimension of pixdim is a timing value, the TR
        % The first three are spatial dimensions
        val = ni.pixdim(4);
    case('sliceduration')
        % Read time for a single slice in seconds
        % Consider adding varargin of time unit (e.g., 'sec','ms', so
        % forth)
        val = ni.slice_duration;
    case {'nvolumes'}
        % This is the number of times the whole brain is measured.
        % A super slice is measured every slice duration, and a volume is
        % acquired every tr.
        val = ni.dim(4);
    case{'nslices','totalslices'}
        % Total number of slices
        val = ni.dim(3);
    case {'superslices'}
        % The slices are acquired continuously, so the TR divided by the
        % slice duration is the number of super slices (or blocks of
        % slices, or something)
        %
        % We could check that the value is very close to an integer
        val = round(bbGet(ni,'tr')./ bbGet(ni,'slice duration'));
    case{'mux','simultaneousslices'}
        % Number of slices acquired simultaneously
        val = round(ni.dim(3)./bbGet(ni,'super slices'));
    case {'sliceacquisitionorder'}
        % This is currently hard coded for the CNI.  In a better brighter
        % world this information will be stored in the NIFTI
        nslices= bbGet(ni,'super slices');
        val = [(1:2:round(nslices)) (2:2:round(nslices))];

    case{'timing'}
        % This is the moment in time (secs) that each slice is acquired,
        % for all of the slices throughout the data set.  So the returned
        % size is (total slices)*(n 
        % bbGet(ni,'timing')
        % 
        % The moment in time that we start to measure each slice
        nslices= bbGet(ni,'super slices');
        mux    = bbGet(ni,'simultaneous slices');
         
        tr     = bbGet(ni,'tr');
        sDuration = bbGet(ni,'slice duration');
        mux_slice_acq_time = 0:sDuration:tr;
        
        % This is the timing for the super slices
        mux_slice_acq_order   = bbGet(ni,'slice acquisition order');
        [~,idx] = sort(mux_slice_acq_order);
        sliceTime = mux_slice_acq_time(idx);
        
        sAcquisitionTiming = zeros(mux,nslices);
        for ii=1:mux
            sAcquisitionTiming(ii,:) = sliceTime;
        end
        sAcquisitionTiming = sAcquisitionTiming';
        sAcquisitionTiming = sAcquisitionTiming(:);
               
        % The timing matrix
        nVolumes = bbGet(ni,'n volumes');
        timing = zeros(bbGet(ni,'total slices'),nVolumes);
        for k=1:nVolumes
            dt = (k-1)*tr;
            timing(:,k) = sAcquisitionTiming + dt;
        end
        val = timing;
        
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
       
        % This the Matlab structure we return. 
        physio_output = [];
        for k=1:length(physioData)
            if isequal(physioData(k).name(1:6),'PPGDat');
                physio_output.ppg.name = 'PPG';
                physio_output.ppg.srate = 100; % here, we are hardcoding the CNI sampling rate for ECG
                physio_output.ppg.rawdata = load(fullfile(UnTarpPhysio_dir,[fName,'_physio'],physioData(k).name));
            elseif isequal(physioData(k).name(1:6),'RESPDa')
                physio_output.resp.name = 'RESP';
                physio_output.resp.srate = 25; % here, we are hardcoding the CNI sampling rate for respiration
                physio_output.resp.rawdata = load(fullfile(UnTarpPhysio_dir,[fName,'_physio'],physioData(k).name));
            end
        end
        
        % now we are going to time-lock it to the scan, and include a
        % parameter for scan-onset
        scan_duration=ni.dim(4)*ni.pixdim(4);% in sec
        %%%% PPG:
        % chop of the beginning
        physio_output.ppg.data=physio_output.ppg.rawdata(end-round(scan_duration*physio_output.ppg.srate)+1:end);
        % add scan onset:
        physio_output.ppg.scan_onset=zeros(size(physio_output.ppg.data));
        for m=1:round(ni.dim(4))
            physio_output.ppg.scan_onset(round(ni.pixdim(4)*(m-1)*physio_output.ppg.srate+1))=1;
        end
        %%%% RESP:
        % chop of the beginning
        physio_output.resp.data=physio_output.resp.rawdata(end-round(scan_duration*physio_output.resp.srate)+1:end);
        % add scan onset:
        physio_output.resp.scan_onset=zeros(size(physio_output.resp.data));
        for m=1:round(ni.dim(4))
            physio_output.resp.scan_onset(round(ni.pixdim(4)*(m-1)*physio_output.resp.srate+1))=1;
        end
        
        val=physio_output;
        
    case{'ppg_peaks'}
        % returns the peaks in the ppg signal in seconds
        
        physio=bbGet(ni,'physio');

        % first get a rough estimate of onsets for the heartbeat
        signal=physio.ppg.data;
        phys_srate=physio.ppg.srate;
        % set minimum inter-heartbeat interval in seconds
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
        [response_matrix,t] = bbResponse2physio(ni,varargin{1});
        val.response_matrix=response_matrix;
        val.t=t;
        
    case{'ppg_correlation'}
        % gets the correlation with PPG for every voxel
        % adding a slice number is optional in varargin{1}
        val = bbCorrelate2physio(ni,varargin{1});
        
    otherwise
        error('unknown BB parameter');
        
end