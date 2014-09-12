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
        val = ni.pixdim(4)./ni.slice_duration;
    case{'mux'}
        val = ni.dim(3)./BB_get(ni,'nslices');
    case{'timing'}
        nslices=BB_get(ni,'nslices');
        mux=BB_get(ni,'mux');
        tr=ni.pixdim(4);
        mux_slice_acq_order = [[1:mux-1:nslices] [2:mux-1:nslices]];
%         for s=1:nslices
%             mux_slice_acq_time = s/nslices*tr;
%         end
%         mux_slice_acq_time=ni.slice_duration;
        ii=0;
        for s=mux_slice_acq_order
        for m=1:mux
            ii=ii+1;
            unmux_slice_acq_order(ii)=nslices*m+s;
        end
        end
        unmux_slice_acq_time = mux_slice_acq_time * mux;
        
        
%         unmux_slice_acq_order = [nslices*m+s for m in xrange(mux) for s in mux_slice_acq_order]
%         unmux_slice_acq_time = mux_slice_acq_time * 3
% for acq,slice in zip(unmux_slice_acq_time,unmux_slice_acq_order):
% %    print "    slice %02d acquired at time %.3f sec" % (slice+1,acq)
% end
    otherwise
        error('Uknown afq parameter');
        
end