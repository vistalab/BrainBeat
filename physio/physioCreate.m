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

% resp = respInit;
% p.resp = resp;

p.filename = 'filename';

if ~isempty(varargin)
    if ~mod(length(varargin),2)  % Must be param/val pairs
        if strcmp(varargin{1},'nifti')
            ni = varargin{2};
            p.filename = ni.fname;
        elseif strcmp(varargin{1},'niftiFilename')
            ni = niftiRead(varargin{2});
            p.filename = ni.fname;
        else
            % Overwrite default arguments with the input arguments
            for ii = 1:2:length(varargin)
                p = physioSet(p,varargin{ii},varargin{ii+1});
            end
        end
        % get scan duration to now we are going to time-lock it to the scan
        scan_duration = ni.dim(4)*ni.pixdim(4);% in sec

        % check whether PPG physio data exist
        ppg_exist = 0;
        physio_name = [extractBefore(p.filename,'_bold.nii') '_recording-PPG_physio.tsv.gz'];
        physio_name_json = [extractBefore(p.filename,'_bold.nii') '_recording-PPG_physio.json'];
        if exist(physio_name,'file') && exist(physio_name_json,'file') 
            ppg_exist = 1;
        else
            error('ERROR: PPG physio file %s or its json sidecar does not exist',physio_name)
        end
        
        % check whether RESP physio data exist
        resp_exist = 0;
        resp_name = [extractBefore(p.filename,'_bold.nii') '_recording-RESP_physio.tsv.gz'];
        resp_name_json = [extractBefore(p.filename,'_bold.nii') '_recording-RESP_physio.json'];
        if exist(resp_name,'file') && exist(resp_name_json,'file') 
            resp_exist = 1;
        else
            warning('WARNING: RESP physio file %s or its json sidecar does not exist',resp_name)
        end
        
        
        if ppg_exist==1
            % add PPG data in the physio_output
            physio_output.ppg.rawdata = bids.util.tsvread(physio_name);

            % for the sampling frequency srate we have to read the json sidecar
            ppg_jsondata = bids.util.jsondecode(physio_name_json);
            physio_output.ppg.srate = ppg_jsondata.SamplingFrequency;
            physio_output.ppg.name = 'PPG';

            % time-lock PPG data to the scan, and include a parameter for scan-onset
            % chop of the beginning
            physio_output.ppg.data = physio_output.ppg.rawdata(end-round(scan_duration*physio_output.ppg.srate)+1:end);
            % add scan onset:
            physio_output.ppg.scan_onset = zeros(size(physio_output.ppg.data));
            for mm = 1:round(ni.dim(4))
                physio_output.ppg.scan_onset(round(ni.pixdim(4)*(mm-1)*physio_output.ppg.srate+1)) = 1;
            end
        end
        
        
        if resp_exist==1
            % add PPG data in the physio_output
            physio_output.resp.rawdata = bids.util.tsvread(resp_name);

            % for the sampling frequency srate we have to read the json sidecar
            resp_jsondata = bids.util.jsondecode(resp_name_json);
            physio_output.resp.srate = resp_jsondata.SamplingFrequency;
            physio_output.resp.name = 'RESP';
            
            % time-lock RESP data to the scan, and include a parameter for scan-onset
            % chop of the beginning
            physio_output.resp.data = physio_output.resp.rawdata(end-round(scan_duration*physio_output.resp.srate)+1:end);
            % add scan onset:
            physio_output.resp.scan_onset = zeros(size(physio_output.resp.data));
            for mm = 1:round(ni.dim(4))
                physio_output.resp.scan_onset(round(ni.pixdim(4)*(mm-1)*physio_output.resp.srate+1))=1;
            end
        end
        
        p = physio_output;

        % make a figure of output as a check
        if length(varargin)>2
            if isequal(varargin{3},'figure') && varargin{4}==1
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
