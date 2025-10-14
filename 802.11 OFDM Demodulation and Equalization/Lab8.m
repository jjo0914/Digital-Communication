% 기초부터다시 패킷이 어떻게만들어지는지부터
% 1. 지정된길이의 비트들을생성 -> 프리엠블합치기
%  
% 2. -> 일종방법으로 변조
%   -> CP를제외히고 행(Nfft=64)은 서브캐리어수 열은 심볼갯수(비트길이)
%   -> 
% 3. -> nfft=X
% 802.11g = 20MHZ
% 4096비트-> mcs1 = 172심볼 ,mcs4=43 ofdm심볼  
%   STF(8us)+LTF 2개(8us) + L_sig1개 (4us) + 43=48개


% LTF만드는법

clc; clear all;
mcs = 4;  
psduLen = 512; % 바이트
wlantx = WLANTx('psduLen', psduLen, 'mcs', mcs); % MCS: 1(BPSK-3/4) 2(QPSK-1/2) 3(QPSK3/4) 4(16QAM 1/2)
                                                 % MCS가달라도 FFT크기 LTF STF이런건 안바뀜
x = wlantx(); % mcs변조설정에따라 802.11 a/g 심볼 패킷만들고

pktAnal = WLANPktAnalyzer(); % 수신기입장에서 순수 802.11 a/g의 패킷분석
pktAnal.extractLTFSym(x);  % 그냥 DATA CP길이써도된다
                           % 4096비트= MCS 4에서는 43개심볼
                           % 8us는 2심볼
                           % 복조하면 총 48개심볼
                           % stf 2개 +ltf 2개 +l-sig 1개 + 43개심볼 =48

% TODO:  Plot the constellation in pktAnal.symLTF.
scatter(real(pktAnal.symLTF(:,1)), imag(pktAnal.symLTF(:,1)), 12, 'filled'); 
axis equal; grid on; 
xlabel('In-Phase'); ylabel('Quadrature');
title('LTF-1 constellations');
figure;
scatter(real(pktAnal.symLTF(:,2)), imag(pktAnal.symLTF(:,2)), 12, 'filled');
axis equal; grid on; 
xlabel('In-Phase'); ylabel('Quadrature');
title('LTF-2 constellations');

% TODO:  Complete the code in WLANPktAnalyzer.extractPayloadSym.  Then, run
% the following code
pktAnal = WLANPktAnalyzer();
pktAnal(x); %  PAYLOAD만 뽑아내기


% TODO:  Plot the pilot and data constellation points that should be stored
% in pktAnal.symPilot and pktAnal.symData
% 파일럿
figure;
zPil  = pktAnal.symPilot(:);           % 원래 4×Nsym → (4Nsym)×1
scatter(real(zPil), imag(zPil), 14, 'filled'); 
axis equal; grid on; title('Pilots (all symbols)');
xlabel('I'); ylabel('Q');

% 데이터
figure;
zDat  = pktAnal.symData(:);            % 48×Nsym → (48Nsym)×1
scatter(real(zDat), imag(zDat), 10, 'filled');
axis equal; grid on; title('Data (all symbols)');
xlabel('I'); ylabel('Q');

% mcs가 뭐냐에따라 , Mb/s, 총비트가 몇 data 심볼이되는지 결정됨
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

% 총 ofdm symols.
R=mod(mcs,2);% 1/2= 1비트를 보내려 2비트전송 3/4 3비트보내려고 4비트전송?
if R==0
    R=1/2;
elseif R==1
    R=3/4;
end
% 총비트를 R만큼배수한뒤 비트를 늘리고 1변조 심볼당 22비트 붙히고 ->변조 ?
Ncbps=48*bitsPerSym; % 1ofdm심볼에 코딩된비트수
Ndbps=48*bitsPerSym*R;    % 1ofdm심볼에 포함된 databit(48은 데이터인덱스 48개)
Nofdmsym= ceil((psduLen*8 +16+6)/Ndbps)% data구간의 총 ofdm심볼
Ncodeddatabit=Nofdmsym*Ncbps
%  Compute the upsampled data rate
%  802.11a/g 기본 fsamp 20mhz긴해
fsamp=wlantx.fsamp;
ovRatio = 2;
fsampUp = fsamp*ovRatio; %-> 수신샘플링은 은 표준의 2배

% Create the channel
% RandMPChan <- 이미 만들어져서 제공함
% 기능 : 20개의 경로지연 & gain만듬
chan = RandMPChan("fsamp", fsampUp); % 기본20개경로

% Generate random multi-paths
chan.genPath();
% chan.dly 1초동안 20개 경로의 지연
% chan.gain < 각경로의  gain  
%  t = ... t(i) = chan.dly(i)-min(chan.dly) 
%"가장 빠른 경로보다 얼마나 더 늦게 오는지(초과 지연 t)를 구해서,
%  stem(...) Plot chan.gain(i) vs. t(i) 
t=chan.dly-min(chan.dly);
figure;
stem(-1*chan.gain,t*1e9);
xlabel('dB 손실');ylabel('ns 지연'); % 축바뀔수잇음



% Create the TX and packet analyzer
% 왜또 Tx랑 PktAnalyzer하는데?
wlantx = WLANTx('psduLen', psduLen, 'mcs', mcs);
pktAnal = WLANPktAnalyzer();

% Compute the signal parameters
info = wlanNonHTOFDMInfo('NonHT-Data');

% Active= data+pilot 나머지는 가드 밴드
sigBW = (max(info.ActiveFFTIndices)-min(info.ActiveFFTIndices)+1)/info.FFTLength*fsamp;

% Create the upsampling and downsampling filters
% Txfilt Rxfilt는 LAB7에서 만듬
% 같은필터를 송수신에 쓰면 매치드 필터??! 실수 FIR필터=
% Txfilt: 2배 오버샘플링 0패딩(isi억제?기억안남) -> gain 낮춤 (true일경우)-> ofdm신호주파수필터링
% Rxfilt: 오버샘플링만큼 낮추고 -> 필터링 
ovRatio = 2;
fsampUp = fsamp*ovRatio;
fp=sigBW/2; % 시그널이존재하는 주파수
%fst= 스탑밴드를 어케잡아야하는지.. 그냥디폴트값써?
fst=fsamp/2;

% %두개 필터는같음=매치드필터=실수필터일경우(만든게실수필터래)
% 필터: 1. fdesign.lowpass: 지정한 주파수만 통과시키는 저역필터
% 필터원리,탭3 : y[n]=h[0]x[n]+h[1]x[n-1]+h[2]x[n-3] -> 인접샘플끼리 평균=저주파통과이런느낌
% 필터메소드(원리몰루)들이 탭을 결정하고 이런탭들이 필터의 특징을 결정함
%  FIR필터(equip)는 분자(탭) /사용법:conv(계수,신호)
%  IIR필터(butter)는 분자와 분모가있음/ 사용법: filter(분자,분모,신호)
txFilt = TxFilt('ovRatio', ovRatio, 'rateIn', fsamp, 'Fp', fp,'Fst',fst); % equripple필터
% 원래신호20mhz-> 40mhz -> 20mhz제한(필터) -> fsampup으로샘플링 -> 20mhz제한(필터)
rxFilt = RxFilt('ovRatio', ovRatio, 'rateIn', fsampUp, 'Fp', fp,'Fst',fst);

% Create the packet detector 
wlanDet = WLANDetect('psduLen',psduLen,'mcs',mcs); % WLANDetect 없는데 -> Unit7의 WLANRx
                        % Unit7, WLANRx의설정바꿔저야함
%Create a RX signal
snr = 10; % ->채널의잡음추가

% Generate the TX packet, analyze and up-sample
x = wlantx(); % -> x=802.11 ag신호생성
pktAnal(x); % -> LTF랑 Payload(data+pilot)추출
xup = txFilt(x); % -> x 송신필터링 

% Get the noise variance to match the SNR
Ex = mean(abs(xup).^2);
wvar = db2pow(-snr)*Ex;
chan.wvar = wvar;  % ->채널의잡음추가

% Pass data through channel
rup = chan(xup); % -> 채널적용 잡음+멀티패스20개 ㄷ

% Downsample
 r = rxFilt(rup); % -> 이제 패팃 검증

% % Perform the STF and LTF detection
%wlanDet.set('pktLen',length(r));
wlanDet.pktDetect(r); % 못찾는데?ㅋㅋ. 식바꾸니찾노?
%   wlanDet.rhoSTF -> 시작인덱스에서 윈도144 일때 상관값 
%   즉 istf 에서최대값이나온다면 그위치가맞는듯?

 pktFound= wlanDet.pktFound; %[pktFound, indltf] = wlanDet(r);
 indltf=wlanDet.istf ;  % iltf해야하는데.. istf로일단넘겨
 indltf=wlanDet.istf + 160; % 20mhz에서 160샘플아님?

% % TODO complete the code in WLANRx.compChanEstRaw.  Then, run the following
% % code
% wlanrx = WLANRx();
% wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData', nsymData);
% wlanrx.compChanEstRaw(r, indltf);
  wlanrx = WLANRx(); % TODO: 채울것
  wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData', Nofdmsym); % nsymData=시간축의 data의 ofdm심볼
  wlanrx.compChanEstRaw(r, indltf); % 원리? 그냥구간추출해서 ofdmdemod
                                    % Raw채널추정 2가지경우로 나옴

% % TODO:   Plot the real components of these two channel estimates.
 figure;
 title('Raw채널추정(주파수관점)')

 plot(real(wlanrx.chanEstRaw))
 xlabel('통과하는주파수 index  (1~64=20mhz)'); 


% TODO complete the code in WLANRx.compChanEst.  Then, run the following
% code to compute the 
wlanrx = WLANRx();
wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData', Nofdmsym);
% -> 원본정답을 WLANRx에넘김
wlanrx.compChanEst(r, indltf);

% TODO:  Plot the real component of the raw channel estimate on both symbols with `o` marker. 
% On the same plot, plot the real component of the smoothed channel
% estimate
hold on

 plot(real(wlanrx.chanEst),'o-')
 legend('Raw추정1','Raw추정2','스무딩(ㄹㅇ채널가정)')



% TODO:  Complete the code in WLANRx.eqSym, then run the following code
% which will perform the equalization

% Inversion 방법: data구간에서 ofdmdemod하고 추정채널로 나눔=복조data
% MMSE 방법 : data구간에서 ofdmdemod하고 ?
 eqMethod = 'inversion';  % Set to either 'inversion' or 'MMSE'\
                     % MMSE 산점도 잘 안모여있듬..
  wlanrx = WLANRx('eqMethod', eqMethod);
 wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData',Nofdmsym);
  wlanrx(r, indltf);

% % TODO:  The equalized data symbols are stored in wlanrx.symDataEq.
% % Plot the equalized constellation points on the complex plane.
% 
% 데이터
figure;
zDat2  = wlanrx.symDataEq;           % 48×Nsym → (48Nsym)×1
scatter(real(zDat2), imag(zDat2), 10, 'filled');
axis equal; grid on; title('채널 추정 복조 Data (all symbols)');
xlabel('I'); ylabel('Q');

%-> ------------------SDR전송부분  ----------------------

loopback = true;  

% Select to run TX and RX
runTx = true;
runRx = true;
% Create the TX and packet analyzer
psduLen = 512;
mcs = 4;
wlantx = WLANTx('psduLen', psduLen, 'mcs', mcs);
pktAnal = WLANPktAnalyzer();

% Compute the signal parameters
info = wlanNonHTOFDMInfo('NonHT-Data');
fsamp = wlantx.fsamp;
sigBW = (max(info.ActiveFFTIndices)-min(info.ActiveFFTIndices)+1)/info.FFTLength*fsamp;

% Create the upsampling and downsampling filters
ovRatio = 2;
fsampUp = fsamp*ovRatio;
txFilt = TxFilt('ovRatio', ovRatio, 'rateIn', fsamp, 'Fp', fp,'Fst',fst); % equripple필터
% 원래신호20mhz-> 40mhz -> 20mhz제한(필터) -> fsampup으로샘플링 -> 20mhz제한(필터)
rxFilt = RxFilt('ovRatio', ovRatio, 'rateIn', fsampUp, 'Fp', fp,'Fst',fst);
% Create the packet detector 
wlanDet = WLANDetect('psduLen',psduLen,'mcs',mcs); 

% TODO:  Generate the TX packet, analyze and up-sample 
%   x = wlantx();
%   pktAnal(...);
x = wlantx(); % -> x=802.11 ag신호생성
pktAnal(x); % -> LTF랑 Payload(data+pilot) OFDM추출



% TODO:  Get the number of symbols of data from the 
% OFDM 심볼 갯수
% size of pktAnal.symData;
nsymData=size(pktAnal.symData,2); % 43개 이응이응

% Set the packet size for the detector to ensure that it detects 
wlanDet.pktLen = length(x); % 실제 패킷길이를넘겨서
                            % 받은신호에서 데이터부분의 샘플길이를 잘라
                            % iltf or istf 더찾기쉽게 할려고


% TODO:  Create the up-sampled TX data
%   xup = ...   
xup = txFilt(x); % -> x 송신필터링 

% clear previous instances
clear sdrtx sdrrx

% add path to the common directory where the function is found


% Run the creation function.  
nsampsFrame = length(xup);
[sdrtx, sdrrx] = plutoCreateTxRx(createTx = runTx, createRx = runRx, loopback = loopback, ...
    nsampsFrameTx= nsampsFrame, nsampsFrameRx=nsampsFrame, sampleRate = fsampUp);

if runTx
    % TODO:  Use the sdrtx.release() and sdrtx.transmitRepeat() commands to
    % continuously send xup
    sdrtx.release();
    sdrtx.transmitRepeat(xup);
end
% TODO:  Capture data at the upsampled rate.  Capture 2*nsampsFrame since
% we want to make sure that we capture a full packet
%   rup = rx.capture(...)
rup=sdrrx.capture(2*nsampsFrame); 

% TODO:  Scale to floating point
%   rup = ... 
rup=single(rup)./2048;
% TODO:  Down-sample
%   r = rxFilt(...);
r=rxFilt(rup);

% TODO:  Run detection with the wlanDet class
%   [pktFound, indltf] = ...;
wlanDet.pktDetect(r); 

 pktFound= wlanDet.pktFound; %[pktFound, indltf] = wlanDet(r);
 indltf=wlanDet.istf ;  % iltf해야하는데.. istf로일단넘겨
 indltf=wlanDet.istf + 160; % 20mhz에서 160샘플


% Create the RX and set the packet info
eqMethod = 'inversion';  % Set to either 'inversion' or 'MMSE'
wlanrx = WLANRx('eqMethod', eqMethod);
wlanrx.setPktData("symLTFTx", pktAnal.symLTF, 'nsymData', nsymData);

% TODO:  Run the receiver
 wlanrx(r, indltf);

% TODO: Plot the equalized constellation points on the complex plane.
figure;
zDat3  = wlanrx.symDataEq;           % 48×Nsym → (48Nsym)×1
scatter(real(zDat3), imag(zDat3), 10, 'filled');
axis equal; grid on; title('SDR Loop back 채널 추정 복조 Data (all symbols)');
xlabel('I'); ylabel('Q');
% 복원잘된다

% % -------------------------지속모니터링 ->  readme github작성 (귀찮..)
% Number of iterations
nit = 100;

% Initialize arrays to store data
snrEst = zeros(nit,1);                  
rxTime = zeros(nit,1);

% Initialize plots
f=figure;
subplot(1,2,1);
symEqReal = [0];
symEqImag = [0];
peq = scatter(symEqReal, symEqImag, 'o');
grid on;
xlim([-15,15]);
ylim([-15,15]);
title('Equalized symbols');

subplot(1,2,2);
snrToNow = [0,0];
timeToNow = [0,3];
psnr = scatter(timeToNow,snrToNow, 'o');
xlabel('Time [sec]');
ylabel('SNR [dB]');
grid on;
ylim([0,40]);
title('SNR Estimate');

% Set pointers to X and Y data for both plots
peq.XDataSource = 'symEqReal';
peq.YDataSource = 'symEqImag';
psnr.XDataSource = 'timeToNow';
psnr.YDataSource = 'snrToNow';

for it = 1:nit



    % TODO:  Get the data and convert to floating point
    %   rup = ...
  rup=sdrrx.capture(2*nsampsFrame);
rup=single(rup)./2048;
    % Get the time
    if (it==1)
        tic;
    end    
    rxTime(it) = toc();

    % TODO:  Down-sample
    %   r = rxFilt(...);
r=rxFilt(rup);    
    % TODO:  Run detection with the rx.pktDetect method.
    %   [pktFound, indltf] = ...
wlanDet.pktDetect(r); 

 pktFound= wlanDet.pktFound; %[pktFound, indltf] = wlanDet(r);
 indltf=wlanDet.istf ;  % iltf해야하는데.. istf로일단넘겨
 indltf=wlanDet.istf + 160; % 20mhz에서 160샘플

    if ~pktFound
        % If packet is not found set SNR to a low value
        % and empty the RX constellation
        snrEst(it) = 0;
        symEqReal = [0];
        symEqImag = [0];

    else

        % TODO:  Run the receiver
        %   wlanrx(...);
wlanrx(r, indltf);
        % Set the estimated SNR
        snrEst(it) = wlanrx.snrEst;
zDat4  = wlanrx.symDataEq;

        % TODO: Get the equalized data symbols 
        %   symEqReal = ...
        %   symEqImag = ...
            symEqReal=real(zDat4(:)); % 왜 애만 *1 벡터로바꿔야지?
            symEqImag=imag(zDat4(:));
    end


    % Update the plots   

    timeToNow = rxTime(1:it);
    snrToNow = snrEst(1:it);    
    refreshdata;
    drawnow;
    exportgraphics(f, 'Lab8result.gif', 'Append', true);  % 프레임 추가 exportgraphics(f exportgraphics(f, 'anim.gif', 'Append', true);  % 프레임 추가, 'anim.gif', 'Append', true);  % 프레임 추가
end
