function [output] = meg_cfc_pipeline_demo(data,vertices,amps,phases,wind,blwind,sr,entfrq,surrs)
% meg_cfc_pipeline_demo
% Phase-amplitude coupling routine for Murphy et al, 2019
%
% References
%   Seymour, R. A., Rippon, G., & Kessler, K. (2017). The Detection of Phase Amplitude Coupling during Sensory Processing. Frontiers in neuroscience, 11, 487. doi:10.3389/fnins.2017.00487
%   Canolty, R. T., & Knight, R. T. (2010). The functional role of cross-frequency coupling. Trends in cognitive sciences, 14(11), 506–515. doi:10.1016/j.tics.2010.09.001
%
% Inputs
% - data = timeseries data from MEG/EEG (chans x time x trials) --- this
%   example uses an roi approach to average over channels in a region. this
%   can be edited to consider individual channels
% - vertices = logical index of channels/vertices to average over to create
%   an ROI
% - amps = amplitude frequencies (Hz)
% - phases = phase frequencies (Hz)
% - wind = window to denote the none baseline portion of the data epochs
%   (samples)
% - blwind = window to denote the baseline portion of the data epochs
%   (samples)
% - sr = sampling rate (Hz)
% - entfrq = entrainment frequency (Hz) --- this is the amplitude frequency
%   of interest
% - surrs = number of surrogates
% Outputs
% - output.precfc.ampdata = ampdata;
% - output.precfc.phasedata = phasedata;
% - output.precfc.data = data;
% - output.precfc.Hilb_ampdata = Hamp; -- hilbert transformed amp
% - output.precfc.Hilb_phasedata = Hphase; -- hilbert transformed phase
% - output.ratio.ratiosS = ratiosS; -- non-sinusoidal to sinusoidal ratio
% - output.ratio.ratiosB = ratiosB;
% - output.PAC.stim = Xs;
% - output.PAC.baseline = Xb;
% - output.PAC.stim_surr = Xss; -- surrogates for stim period
% - output.PAC.bl_surr = Xsb; -- surrogates for baseline period
% - output.PAC.ZXs = Zhypcfcs; -- z transformed against surrogates
% - output.PAC.ZXb = Zhypcfcb;
% ---- for the z transformation the value that is considered is the
% hypothesised phase.

%% Organise variables for CFC
dataclone = data; % backup
data = nanmean(data(vertices,:,:),1); % must be in format chan/vert * time * trial
% get dimension sizes
[~,tp,trr] = size(data);
% set indices for cfc normalisation
ioi = false(tp,1);
ioi(wind) = true;
ioib = false(tp,1);
ioib(blwind) = true;
% create window to normalise time series amplitudes to
normwind = false(tp,1);
normwind(200:length(normwind)-200) = true; % exclude edges


%% Filter the timeseries
% see subroutines below
[ampdata,phasedata] = bpfilterSUB(amps,phases,data,sr,tp,trr);
[Hamp,Hphase] = sigtransSUB(ampdata,phasedata);

%% Check for non-sinusoidal waveforms
% routine described in Seymour et al, 2017
[ratiosS] = ns_rise_decaySUB(phasedata(:,wind,:));
[ratiosB] = ns_rise_decaySUB(phasedata(:,blwind,:));


%% Phase-amplitude coupling

[Xs]=cfcSUB(Hphase,Hamp,phases,sr,ioi,normwind); % stimulation cfc
[Xb]=cfcSUB(Hphase,Hamp,phases,sr,ioib,normwind); % baseline cfc
% prepare surrogate cfc storage
%%% surrogate options are included as an optional extra
Xss = zeros(numel(amps),numel(phases),surrs);
Xsb = zeros(numel(amps),numel(phases),surrs);
[~,surrtemplates] = randpermsurr(phasedata,surrs,trr,2,wind);
[~,surrtemplateb] = randpermsurr(phasedata,surrs,trr,2,blwind);
parfor ii = 1:surrs
    % routine described in Canolty et al, 2010
    [Xss(:,:,ii)]=cfcSUB(surrtemplates(:,:,:,ii),Hamp,phases,sr,ioi,normwind);
    [Xsb(:,:,ii)]=cfcSUB(surrtemplateb(:,:,:,ii),Hamp,phases,sr,ioib,normwind);
    disp(num2str(ii))
end

%% Z-transform
% hypothesis driven approach

hdxss = zeros(1,surrs);
hdxsb = zeros(1,surrs);
% find indices to use
test=phases.*amps';
hypphase = 6;
hypamp = entfrq;
f = find(phases==hypphase);
ff = find(test(:,f)==hypphase*hypamp);
% z test val = matrix(ff,f);
for ii = 1:surrs
    hdxss(:,ii) = squeeze(Xss(ff,f,ii));
    hdxsb(:,ii) = squeeze(Xsb(ff,f,ii));
end
Zhypcfcs = zeros(size(Xs));
Zhypcfcb = zeros(size(Xb));
Zhypcfcs= (Xs-mean(hdxss))./std(hdxss);
Zhypcfcb = (Xb-mean(hdxsb))./std(hdxsb);

%% Output
output.precfc.ampdata = ampdata;
output.precfc.phasedata = phasedata;
output.precfc.data = data;
output.precfc.Hilb_ampdata = Hamp;
output.precfc.Hilb_phasedata = Hphase;
output.ratio.ratiosS = ratiosS;
output.ratio.ratiosB = ratiosB;
output.PAC.stim = Xs;
output.PAC.baseline = Xb;
output.PAC.stim_surr = Xss;
output.PAC.bl_surr = Xsb;
output.PAC.ZXs = Zhypcfcs;
output.PAC.ZXb = Zhypcfcb;

end

%% Subroutines



function [ampdata,phasedata] = bpfilterSUB(amps,phases,TS,sr,tp,tr)
% High Frequency
ampband = zeros(length(amps),2);
ampband(:,1) = amps-max(phases);
ampband(:,2) = amps+max(phases);
ampdata = zeros(length(amps),tp,tr);
parfor iter = 1:length(amps)
    tempA = eegfilt(TS,sr,ampband(iter,1),ampband(iter,2),0,[],0,'fir1',0);
    tempA = squeeze(reshape(tempA,size(TS,1),tp,tr));
    ampdata(iter,:,:) = tempA;
end
clear iter tempA
% Low Frequency
phaseband = zeros(length(phases),2);
phaseband(:,1) = phases-1;
phaseband(:,2) = phases+1;
phasedata = zeros(length(phases),tp,tr);
parfor iter = 1:length(phases)
    tempP= eegfilt(TS,sr,phaseband(iter,1),phaseband(iter,2),0,[],0,'fir1',0);
    tempP = squeeze(reshape(tempP,1,tp,tr));
    phasedata(iter,:,:) = tempP;
end
clear iter tempP
end

function [Hamp,Hphase] = sigtransSUB(ampdata,phasedata)
Hamp = abs(hilbert(permute(ampdata,[2,1,3]))).^2;
Hphase = angle(hilbert(permute(phasedata,[2,1,3])));
Hamp = permute(Hamp,[2,1,3]);
Hphase = permute(Hphase,[2,1,3]);
end

function [ratios] = ns_rise_decaySUB(phasedata)

for ii = 1:size(phasedata,1)
    test = phasedata(ii,:,:);
    p=1;
    for iter = 1:size(phasedata,3)
        trl = test(:,:,iter);
        trlfl = trl.*-1;
        [~,peak_locations] = findpeaks(trl);
        [~,trough_locations] = findpeaks(trlfl);
        % Equalise the number of peak and trough events
        if ~isempty(peak_locations) && ~isempty(trough_locations)
            if length(peak_locations) > length(trough_locations)
                peak_locations(1) = [];
            elseif length(peak_locations) < length(trough_locations)
                trough_locations(1) = [];
            end
            
            % Calculate time to peak and time to decay
            time_to_decay = [];
            time_to_peak = [];
            if peak_locations(1)<trough_locations(1) %if peak first
                for i = 1:length(peak_locations)-1
                    time_to_decay(i) = trough_locations(i)-peak_locations(i);
                    time_to_decay_all(p) = trough_locations(i)-peak_locations(i);
                    time_to_peak(i) = abs(peak_locations(i+1)-trough_locations(i));
                    time_to_peak_all(p) = abs(peak_locations(i+1)-trough_locations(i));
                    p = p+1;
                end
            elseif peak_locations(1)>trough_locations(1) %if trough first
                for i = 1:length(peak_locations)-1
                    time_to_decay(i) = peak_locations(i)-trough_locations(i);
                    time_to_decay_all(p) = peak_locations(i)-trough_locations(i);
                    time_to_peak(i) = abs(trough_locations(i+1)-peak_locations(i));
                    time_to_peak_all(p) = abs(trough_locations(i+1)-peak_locations(i));
                    p = p+1;
                end
            end
            ratios(ii,iter) = mean(time_to_decay)./mean(time_to_peak);
        elseif ~isempty(peak_locations) || ~isempty(trough_locations)
            ratios(ii,iter) = 0;
            
        end
    end
end
f = isnan(ratios);
ratios(f) = 0;
end

function [surrdata,surrtemplate] = randpermsurr(data,surrs,tr,method,wind)
xV = zeros(1,length(wind)*tr);
% go through filtered data and shuffle
surrdata = zeros(size(data,1),length(wind),tr,surrs);
temp = data(:,wind,:);
[xx,yy,zz] = size(temp);
temp = temp(:,:);
% get samples
randtime = randsample(round(length(xV)*.8),surrs)+round(length(xV)*.1);
parfor ii = 1:surrs
    temp2 = [temp(:,randtime(ii):end) temp(:,1:randtime(ii)-1)];
    temp2 = reshape(temp2,xx,yy,zz);
    if method == 1
        temp2 = abs(hilbert(permute(temp2,[2,1,3]))).^2;
        temp2 = permute(temp2,[2,1,3]);
    elseif method == 2
        temp2 = angle(hilbert(permute(temp2,[2,1,3])));
        temp2 = permute(temp2,[2,1,3]);
    end
    surrdata(:,:,:,ii) = temp2;
    %     clear temp*
end
surrtemplate = zeros(size(data,1),size(data,2), size(data,3),size(surrdata,4));
surrtemplate(:,wind,:,:) = surrdata;
end

function [X]=cfcSUB(ph,amp,fPH,sr,ioi,normwind)

[phoi,ntp,ntr]=size(ph);
ampoi=size(amp,1);
if ~exist('ioi','var')||isempty(ioi)
    ioi=true(1,ntp);
end
% amplitude normalisation
mn=repmat(min(amp(:,normwind,:),[],2),[1 ntp 1]);
rng=repmat(range(amp(:,normwind,:),2),[1 ntp 1]);
amp=eps+(amp-mn)./rng;
amp_old=amp;
amp=resh(amp(:,ioi,:));
ph=resh(ph(:,ioi,:));
ioi=repmat(ioi,[1 ntr]);

if exist('fPH','var')
    fPH=fPH(:);
end
if ~exist('sr','var')||isempty(sr)
    sr=250;
end

% cfc - based on Canolty et al, 2010
for nph=1:phoi
    X(:,nph)=abs(mean(amp.*exp(1i*repmat(ph(nph,:),[ampoi 1])),2)) ;
end

end