%% Simulate PAC Data

% In this example we create a signal with a baseline period that contains
% 6/40 hz pac at a low amplitude with white noise, and an active period
% with greater amplitude PAC and white noise.


%simulate signal for pac
timer = 2; % seconds
sr = 1000;
trls = 500;
low = 6;
high = 40;
lowfrq = randi([-2 2], trls, 1);
hifrq = randi([-2 2], trls, 1);
phaseshift = pi/2.*rand(trls,1);
noiseratio = 0.1;
ampratio = 1;
pacsig = zeros(1,sr*timer,trls);
sigx = 1:(sr*timer);
lvars = low+lowfrq;
hvars = high+hifrq;
for iter = 1:trls
    noisesignal = noiseratio.*randn(max(sigx),1);
    carriersig = sin((lvars(iter)/sr)*2*pi*sigx+phaseshift); % modulating signal
    messagesig = (sin((hvars(iter)/sr)*2*pi*sigx)).*ampratio; % signal to be modulated
    modulatedsig = (messagesig.*(carriersig+1)); % signal with PAC
    simsig = modulatedsig'+noisesignal; % add noise to PAC signal
    simsig = (simsig(:,1) + ampratio*noisesignal(:,1)); % scale amplitude of noise
    pacsig(1,:,iter) = simsig;
    clear noisesignal
end

% create baseline and event distinction
noisesignal = noiseratio.*randn(max(sigx),1);
bl = 1:1000;
for iter = 1:500
    pacsig(1,bl,iter) = (pacsig(1,bl,iter)/5).*rand; % scale amplitude of baseline
    pacsig(1,:,iter) = pacsig(1,:,iter)+noisesignal'; % add additional noise
    pacsig(1,:,iter) = circshift(pacsig(1,:,iter),-100);
    % we are shifting the end of the trial back by 100 ms to leave breathing
    % room at the end of the trial
end
newbl = 1:900; % this will be used to window the baseline in the pac routine
cfc_data = pacsig;











