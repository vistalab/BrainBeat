function p = physioCreate(varargin)
% Create a physiological data structure for heart and respiration
%
%  p = physioCreate(varargin)
%
% The heart beat is estimated with a photoplethysmogram
% The respiration is estimated with a breathing belt.
%
% Example:
%    p = physioCreate; %returns an empty structure
% 
%    p = physioCreate('filename',fname); 
% 
%    WORK WITH THE NIFTI STRUCTURE
%       enter the nifti structure loaded with
%       niftiRead, this assumes that you have the physio data as they come
%       from the CNI, and are named ni.filename_physio.tgz
%    p = physioCreate('nifti',ni); 
%       optional figure flag:
%    p = physioCreate('nifti',ni,'figure',f_flag) f_flag = 1 to add a figure
%
%    WORK WITH SEPARATE FILENAMES
%    p = physioCreate('niftiFilename',fname,'ppg',ppgname,'resp',respname); 
%     
%    p = physioCreate('ppg srate',50);
%
%    p = physioCreate('resp srate',33)
%
% See also:  bbGet
%
% DH/BW Vistasoft team, 2014

ppg = ppgInit;
p.ppg = ppg;

resp = respInit;
p.resp = resp;

p.filename = 'filename';

if ~isempty(varargin)
    if ~mod(length(varargin),2)  % Must be param/val pairs
        if strcmp(varargin{1},'filename')
            p.filename = varargin{2};
            % Do that bbGet thing in here to read the file data
            % Figure that the file is right.
            
            return
            
            
            
        elseif strcmp(varargin{1},'nifti')
            ni=varargin{2};
            p.filename=ni.fname;
            
            % Should add physio file handling to niftiGet(ni,'physio ....');
            [fPath,fName] = fileparts(ni.fname);
            [~, fName]    = fileparts(fName);       
            UnTarpPhysio_dir = fullfile(fPath,[fName '_physio']);

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
            physioData = dir(fullfile(UnTarpPhysio_dir,'*Data*'));

            % This the Matlab structure we return. 
            physio_output = [];
            for k=1:length(physioData)
                if isequal(physioData(k).name(1:6),'PPGDat');
                    physio_output.ppg.name = 'PPG';
                    physio_output.ppg.srate = 100; % here, we are hardcoding the CNI sampling rate for ECG
                    physio_output.ppg.rawdata = load(fullfile(UnTarpPhysio_dir,physioData(k).name));
                elseif isequal(physioData(k).name(1:6),'RESPDa')
                    physio_output.resp.name = 'RESP';
                    physio_output.resp.srate = 25; % here, we are hardcoding the CNI sampling rate for respiration
                    physio_output.resp.rawdata = load(fullfile(UnTarpPhysio_dir,physioData(k).name));
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

            p=physio_output;
            
            % make a figure of output as a check
            if length(varargin)>2
                if isequal(varargin{3},'figure') && varargin{4}==1
                    figure
                    subplot(2,1,1),hold on
                    plot([1:length(p.ppg.data)]/p.ppg.srate,p.ppg.data)
                    plot(find(p.ppg.scan_onset>0)/p.ppg.srate,mean(p.ppg.data),'r.')
                    legend('ppg data','scan onset')
                    xlabel('sec'),ylabel('ppg data')
                    subplot(2,1,2),hold on
                    plot([1:length(p.resp.data)]/p.resp.srate,p.resp.data)
                    plot(find(p.resp.scan_onset>0)/p.resp.srate,mean(p.resp.data),'r.')
                    legend('resp data','scan onset')
                    xlabel('sec'),ylabel('resp data')
                end
            end
            
            return
            
            
            
            
        elseif strcmp(varargin{1},'niftiFilename')
            ni=niftiRead(varargin{2});
            p.filename=ni.fname;
          
            % This the Matlab structure we return. 
            physio_output = [];
            physio_output.ppg.name = 'PPG';
            physio_output.ppg.srate = 100; % here, we are hardcoding the CNI sampling rate for ECG

            physio_output.resp.name = 'RESP';
            physio_output.resp.srate = 25; % here, we are hardcoding the CNI sampling rate for respiration

            % now load the data
            % quick check on ppg and resp order
            ppgFilename = [];
            respFilename = [];
            for k=1:2:length(varargin)
                if isequal(varargin{k},'ppg')
                    ppgFilename = varargin{k+1};
                elseif isequal(varargin{k},'resp')
                    respFilename = varargin{k+1};
                end
            end
            
            % if there is a PPG/RESP filename, put the data in the ppg structure
            % and time-lock it to the scan, include a parameter for scan-onset
            scan_duration=ni.dim(4)*ni.pixdim(4);% in sec
            
            if ~isempty(ppgFilename) 
                physio_output.ppg.rawdata = load(ppgFilename);
                %%%% PPG:
                % chop of the beginning
                physio_output.ppg.data=physio_output.ppg.rawdata(end-round(scan_duration*physio_output.ppg.srate)+1:end);
                % add scan onset:
                physio_output.ppg.scan_onset=zeros(size(physio_output.ppg.data));
                for m=1:round(ni.dim(4))
                    physio_output.ppg.scan_onset(round(ni.pixdim(4)*(m-1)*physio_output.ppg.srate+1))=1;
                end
            else
                physio_output.ppg.data=[];
                warning('no ppg filename entered')
            end
            
            if ~isempty(respFilename)
                physio_output.resp.rawdata = load(respFilename);
                %%%% RESP:
                % chop of the beginning
                physio_output.resp.data=physio_output.resp.rawdata(end-round(scan_duration*physio_output.resp.srate)+1:end);
                % add scan onset:
                physio_output.resp.scan_onset=zeros(size(physio_output.resp.data));
                for m=1:round(ni.dim(4))
                    physio_output.resp.scan_onset(round(ni.pixdim(4)*(m-1)*physio_output.resp.srate+1))=1;
                end
            else
                physio_output.resp.data=[];
                warning('no resp filename entered')
            end


            p=physio_output;
            
            % make a figure of output as a check
            if length(varargin)>6
                if isequal(varargin{7},'figure') && varargin{8}==1
                    figure
                    if ~isempty(physio_output.ppg.data)
                        subplot(2,1,1),hold on
                        plot([1:length(p.ppg.data)]/p.ppg.srate,p.ppg.data)
                        plot(find(p.ppg.scan_onset>0)/p.ppg.srate,mean(p.ppg.data),'r.')
                        legend('ppg data','scan onset')
                        xlabel('sec'),ylabel('ppg data')
                    end
                    if ~isempty(physio_output.resp.data)
                        subplot(2,1,2),hold on
                        plot([1:length(p.resp.data)]/p.resp.srate,p.resp.data)
                        plot(find(p.resp.scan_onset>0)/p.resp.srate,mean(p.resp.data),'r.')
                        legend('resp data','scan onset')
                        xlabel('sec'),ylabel('resp data')
                    end
                end
            end

            return

            
            

        else
            % Overwrite default arguments with the input arguments
            for ii=1:2:length(varargin)
                p = physioSet(p,varargin{ii},varargin{ii+1});
            end
        end
    end
end

end


%%%%
function ppg = ppgInit
% Photoplethysmogram data structure
%
%
% Time scale is always at the CNI sampling rate
%
%
ppg.name = 'PPG';
ppg.srate = 100;        % Number of samples per second

ppg.rawdata    = [];    % Values in the file
ppg.timeseries = [];
ppg.scan_onset = [];    % Times when the scan acq begins

end

%%%%
function resp = respInit

resp.name = 'RESP';
resp.srate = 25;

resp.rawdata    = [];    % Values in the file
resp.timeseries = [];
resp.scan_onset = [];    % Times when the scan acq begins

end
