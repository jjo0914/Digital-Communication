mcs = 4;  
psduLen = 512;
wlantx = WLANTx('psduLen', psduLen, 'mcs', mcs); % 'mcs' 기본값3 설정4(QPSK)
x = wlantx();
% TODO:  Compltee WLANPktAnalyzer constructor and extractSym method, then
% run the following.
pktAnal = WLANPktAnalyzer();
pktAnal.extractLTFSym(x); % 원본데이터에서 LTF뽑아내는거니깐..

% TODO:  Plot the constellation in pktAnal.symLTF.
% TODO:  Complete the code in WLANPktAnalyzer.extractPayloadSym.  Then, run
% the following code
pktAnal = WLANPktAnalyzer();
pktAnal(x); %  PAYLOAD만 뽑아내기


% TODO:  Plot the pilot and data constellation points that should be stored
% in pktAnal.symPilot and pktAnal.symData

if mcs<= 1
    bitsPerSym = 1;      % BPSK
elseif mcs<=3
    bitsPerSym = 2;      % QPSK
elseif mcs<=5
    bitsPerSym = 4;      % 16-QAM
else
    bitsPerSym = 6;      % 64-QAM
end

% TODO:  아래변수들을 출력 및 계산
% nscData/ nsymData:  Number of OFDM symbols for the data/nREData =
% nsymData*nscData:  Total number of modulation symbols or resource elements (REs) for the data /
%/ncodedBits:  Total number of coded data bits (use the bitsPerSym above)/codeRate:  Use the psduLen and ncodeBits.  This may be a little < the code rate in the table since there is some zero padding.

%  Compute the upsampled data rate
ovRatio = 2;
fsampUp = fsamp*ovRatio;

% Create the channel
chan = RandMPChan("fsamp", fsampUp); % 기본20개경로

% Generate random multi-paths
chan.genPath()
% TODO:   
%  t = ... t(i) = chan.dly(i)-min(chan.dly) 
%"가장 빠른 경로보다 얼마나 더 늦게 오는지(초과 지연 t)를 구해서, 지연(ns) vs **이득(dB)**
%  stem(...) Plot chan.gain(i) vs. t(i) 

% Create the TX and packet analyzer
wlantx = WLANTx('psduLen', psduLen, 'mcs', mcs);
pktAnal = WLANPktAnalyzer();

% Compute the signal parameters
info = wlanNonHTOFDMInfo('NonHT-Data');
fsamp = wlantx.fsamp;
sigBW = (max(info.ActiveFFTIndices)-min(info.ActiveFFTIndices)+1)/info.FFTLength*fsamp;

% Create the upsampling and downsampling filters
ovRatio = 2;
fsampUp = fsamp*ovRatio;
txFilt = TxFilt('ovRatio', ovRatio, 'rateIn', fsamp, 'sigBW', sigBW);
rxFilt = RxFilt('ovRatio', ovRatio, 'rateIn', fsampUp, 'sigBW', sigBW);

% Create the packet detector 
wlanDet = WLANDetect(); % 전수업에서완성? 

% Create a RX signal
snr = 10;

% Generate the TX packet, analyze and up-sample
x = wlantx();
pktAnal(x);
xup = txFilt(x);

% Get the noise variance to match the SNR
Ex = mean(abs(xup).^2);
wvar = db2pow(-snr)*Ex;
chan.wvar = wvar;

% Pass data through channel
rup = chan(xup);

% Downsample
r = rxFilt(rup);

% Perform the STF and LTF detection
[pktFound, indltf] = wlanDet(r);

wlanrx = WLANRx();
wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData', nsymData);
% TODO complete the code in WLANRx.compChanEstRaw.  Then, run the following
% code
wlanrx = WLANRx();
wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData', nsymData);
wlanrx.compChanEstRaw(r, indltf);

% TODO:   Plot the real components of these two channel estimates.

% TODO complete the code in WLANRx.compChanEst.  Then, run the following
% code to compute the 
wlanrx = WLANRx();
wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData', nsymData);
% -> 원본정답을 WLANRx에넘김
wlanrx.compChanEst(r, indltf);

% TODO:  Plot the real component of the raw channel estimate on both symbols with `o` marker. 
% On the same plot, plot the real component of the smoothed channel
% estimate

% TODO:  Complete the code in WLANRx.eqSym, then run the following code
% which will perform the equalization
eqMethod = 'inversion';  % Set to either 'inversion' or 'MMSE'
wlanrx = WLANRx('eqMethod', eqMethod);
wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData', nsymData);
wlanrx(r, indltf);

% TODO:  The equalized data symbols are stored in wlanrx.symDataEq.
% Plot the equalized constellation points on the complex plane.
    
