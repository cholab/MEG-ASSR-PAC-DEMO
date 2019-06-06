function [ampdata,phasedata] = bpfilterDEMO(amps,phases,TS,sr,tp,tr)
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