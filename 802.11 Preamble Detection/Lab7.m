clc; clf; clear all;
% ------Unit07_Lab-------------
% 문제: 802.11g Wi-fi(가장 간단한 모델) 패킷안에있는 프리엠블탐지하기 (WLAN toolbox 설치)
% 802.11 패킷앞에는 STF,LTF 동기화용신호가 있어서 잡아야함
% 1. auto-corrlation 으로 STF감지
% 2. matched-filter detector로 LTF감지(뭐였지?..)
%------------------------------

%-------802.11g 패킷구조-------------
% STF: 패킷잠지+gain control용 (8us=t1+t2+...t10)
% LTF: 정확한 패킷시작점 추정+ 채널추정(8us=GARD+LT1+LT2)
%-----------------------------------

%------함수써서 802.11g 패킷만들기!-----------
% 패킷을 만드는 WLANTx.m 함수를 만드시오
%-------------------------

wlantx = WLANTx('psduLen', 512); % 512바이트생성

x = wlantx();  % 일단 802.11패킷을만들어 (소스뜯기x)
               %

% TODO:  Get the sample rate from wlantx.fsamp.
fsamp=wlantx.fsamp;

% TODO:  Plot the absolute values of the samples in x vs. time in us.  If you did everything
% correctly, you should see that the packet is about 250 us long. 
t=0+1/fsamp:1/fsamp:1/fsamp*length(x); % 252마이크로초인데?
% 프리엠블 부분이 되게 특이하게 생겼타르
plot(t*1e+6,abs(x));


xlabel('time(us)');
ylabel('value if WLAN samples')
legend('real','abs')

%---프리엠블을 감지하기위해 프리앰블을 생성하기---
wlanrx = WLANRx(); % 프리엠블은똑같으니 나머지설정은 필요x
tstf=wlanrx.lstf; % 160샘플 = 8us
tslf=wlanrx.lltf; % 160샘플 = 8us

%----SDR로 보낼때는 2배로 오버샘플링 -> 스펙트럼 복제생기니 주파수 필터통과------- (원래샘플도 -f/2~f/2  ,f/2~f 쭉복제생김)
%OFDM Bandwidth 계산법 Data+pilot 서브캐리어 갯수 /총 FFTLength * fsamp(=한개의 서브캐리어가 차지하는 주파수)
% Data+Pilot+GUARD=FFTLength  //FFTLength+ CP=OFDM 길이
% fsamp= 유효서브캐리어/총FFT
info=wlanNonHTOFDMInfo('NonHT-Data') % 대역폭 기본값20MHZ
                                     %  NonHT=802.11a/g  OFDM 파라미터들을 알려줌 
                                     % 'NonHT-Data'=NonHT의 Data구간
sigBW=(max(info.ActiveFFTIndices)-min(info.ActiveFFTIndices)+1)/info.FFTLength *fsamp

% Oversampling ratio
ovRatio = 2;

% Set the filter backoff level as in the previous lab
backoffLev = 9; 
% TODO:  Set the filter parameters
 Fp = sigBW/2    % 외워 sigBandwidth/2 = Passband
 Fst = fsamp/2   % fsamp/2 > Stopband
 % Construct the TX filter object
txFilt = TxFilt('rateIn', fsamp, 'Fp', Fp, 'Fst', Fst, 'ovRatio', ovRatio, 'backoffAutoScale', true, 'backoffLev', backoffLev);
%등리플 필터?
%필터공부 ㄱ
% Upsample the output
xup = txFilt(x);

% Get the upsampled rate 
fsampUp = txFilt.rateOut;

% Get the TX filter coefficients
b = txFilt.bfilt;

% TODO:  Use freqz to plot the frequency response of the filter in dB



% Get the average energy per sample in the packet
Etx = mean(abs(xup).^2);

% Channel parameters
dlyRange = [50,100]*1e-6;   % Min and max of the delay in us
gain = 0;                   % channel gain in dB
snr = 20;                   % channel SNR in dB
wvar = Etx*db2pow(-snr);    % channel noise variance
nsampOut = 2^14;            % number of samples of output

% Create the channel 
chan = RandDelayChan('dlyRange', dlyRange, 'fsamp', fsampUp,...
    'gain', gain, 'wvar', wvar, 'nsampOut', nsampOut);

% Run the channel
[rup,dly] = chan(xup);

% Get the true delay in us
dlyus = dly*1e6;

% TODO:   plot abs(rup) vs. time.  Also, plot the delay, dlyus with  a vertical line.
% TODO:  Set the passband and stopband frequencies as before
%   Fp = ...
%   Fst = ...

% Construct the RX filter
rxFilt = RxFilt('ovRatio', ovRatio, 'rateIn', fsampUp, 'Fp', Fp, 'Fst', Fst, 'rxScale',false);

% Down-sample and scale the receive signal
r = rxFilt(rup); % -> 어떤구조지?

% TODO:  Complete the code in WLANRx.detectSTF() as above
% Create an instance of the RX
wlanrx = WLANRx();

% Run the STF detector
wlanrx.detectSTF(r); % ->원리가뭐지?

% TODO:  Plot the correlation wlanrx.rhoSTF vs. time in us.  Also, plot the
% location of the detected maximum
% TODO:  Complete the in WLANRx.pktDetect() to compute the correlation and peak detection with the L-LTF.  
% Test it with the following code.

% Create simulation objects
wlanrx = WLANRx();


% Perform the LTF correlation
wlanrx.pktDetect(r); % ->원리가뭐지?

% TODO:  Plot the LTF correlation in rx.rhoLTF.  

% Parameters
ntest = 100;              % number of points per SNR
snrTest = [-6:2:4]';      % SNR values to test
                    % 일정 SNR구간에서  패킷검출율급락
nsnr = length(snrTest);

% Initialize arrays
dlyEst = zeros(ntest,nsnr);
snrEst = zeros(ntest,nsnr);
dlyTrue = zeros(ntest,nsnr);
pmd = zeros(nsnr,1);
nfound = zeros(nsnr,1);

% Create system objects
wlantx = WLANTx('psduLen', 512);
wlanrx = WLANRx('psduLen', 512);
txFilt = TxFilt('rateIn', fsamp, 'Fp', Fp, 'Fst', Fst, 'ovRatio', ovRatio, 'backoffAutoScale', true, 'backoffLev', backoffLev);
rxFilt = RxFilt('ovRatio', ovRatio, 'rateIn', fsampUp, 'Fp', Fp, 'Fst', Fst, 'rxScale',false);

for isnr = 1:nsnr
    % Construct the simulation objects
    snr = snrTest(isnr);   
    chan = RandDelayChan('dlyRange', dlyRange, 'fsamp', fsampUp,...
        'gain', gain, 'wvar', wvar, 'nsampOut', nsampOut);
    
    i = 0;
    for it = 1:ntest
        
        % TODO:  Generate a random packet using wlantx()
        %    x = ... 
        

        % TODO:   Upsample with the txFilt   
        %    xup = ...
        

        % Compute the noise variance
        Etx = mean(abs(xup).^2);
        wvar = Etx*db2pow(-snr);
        chan.wvar = wvar;


        % TODO:  Pass through a random single path channel using the
        % chan.randChan method
        %     [rup, dly] = chan(...)        

        % TODO:  Down-sample
        %    r = rxFilt(...);
        

        % Run detection with the rx.pktDetect method.  
        wlanrx.pktDetect(r); % 1 or 0

        % If packet is found:
        if wlanrx.pktFound
            i = i + 1;  % Increment counter of number of packets found

            % TODO:  Store:
            %     dlyTrue(i,isnr) = true delay 
            %     dlyEst(i,isnr) = estimated delay based on rx.iltf
              
        end
    end
    
    % Set number of packets found
    nfound(isnr) = i;

    % TODO:  Set     
    %   pmd(isnr) = fraction of packets that were missed.  
    
    
    % Print results
    fprintf(1,"SNR=%7.2f MD=%12.4e \n", snr, pmd(isnr));
end

% TODO:  Plot the missed detection vs. SNR.  You should see that the missed
% detection is low once the SNR is about 2 dB.

% TODO:  Compute the delay offset
%    dlyOff = ...


% TODO:  Set dlyEstAdj = dlyEst + dlyOff, which is the adjusted delay
% estimate.

% SNR values to plot are snr(Iplot(j))
Iplot = [3,5];
nplot = length(Iplot);
for iplot = 1:nplot
    subplot(1,nplot,iplot);
    isnr = Iplot(iplot);     
    n = nfound(isnr);

    % TODO:  Create a scatter plot of dlyTrue(:,isnr) vs. dlyEstAdj(:,isnr)
    % Label the axes
    
end
