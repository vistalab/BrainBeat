function val = BB_get(ni,param)
% Get value from BB structure
%
% val = BB_get(afq, param, varargin)
%
% param list:              - arguments:
%
% 'tr'
% 'te'
% 'fa'
% 'nr_points'
% 'resolution'
% 'dimentions'
% 'physio'
%
% BB_get(ni,'tr')
%
% Wandell Copyright Vistasoft Team, 2013
% Written by Aviv an Dora 2014


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
        for k=1:length(physio_output(k))
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
    otherwise
        error('Uknown afq parameter');
        
end