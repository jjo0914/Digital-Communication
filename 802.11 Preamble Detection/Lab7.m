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
% Scatter plot of the absolute values of the samples in x

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
t=0:1/fsamp:1/fsamp*(length(x)-1); % 252마이크로초인데?
% 프리엠블 부분이 되게 특이하게 생겼타르
plot(t*1e+6,abs(x));


xlabel('time(us)');
ylabel('value if WLAN samples')
legend('abs')

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
sigBW=(max(info.ActiveFFTIndices)-min(info.ActiveFFTIndices)+1)/info.FFTLength *fsamp;

% Oversampling ratio
ovRatio = 2;

% Set the filter backoff level as in the previous lab
backoffLev = 9; 
% TODO:  Set the filter parameters
 Fp = sigBW/2  ;  % 외워 sigBandwidth/2 = Passband
 Fst = fsamp/2;   % fsamp/2 > Stopband

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
[H, F] = freqz(b, 1, 1024,fsampUp);
figure;
plot(F*1e-6, 20*log10(abs(H)));
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Frequency Response of the TX Filter');
grid on;


% 필터 통과한 패킷의 평균에너지 or  샘플당 평균전력(이쪽이맞다)
Etx = mean(abs(xup).^2);

% 이건 채널 h 설정
dlyRange = [50,100]*1e-6;   % Min and max of the delay in us
gain = 0;                   % channel gain in dB
snr = 20;                   % channel SNR in dB -> 조정가능
wvar = Etx*db2pow(-snr);    % channel noise variance
nsampOut = 2^14;            % number of samples of output

% 채널은 64탭? 많으면연산량증가 등..
% 채널 sinc도 아무런 지연이없으면 채널 임펄스함수랑 결과가 같다
% 하지만 분수지연은 임펄스가 구현x sinc함수를쓴
chan = RandDelayChan('dlyRange', dlyRange, 'fsamp', fsampUp,...
    'gain', gain, 'wvar', wvar, 'nsampOut', nsampOut);

% Run the channel
[rup,dly] = chan(xup);

% Get the true delay in us/
dlyus = dly*1e6;

% TODO:   plot abs(rup) vs. time.  Also, plot the delay, dlyus with  a vertical line.
t=0:1:length(rup)-1; t=t./fsampUp;
figure;
plot(t*1e+6,abs(rup))
hold on;
yLimits = ylim; % Get current y-axis limits
xline(dlyus, 'r--', 'Delay');
ylim(yLimits); % Restore y-axis limits
ylabel('Time(us)');
hold off;
% TODO:  Set the passband and stopband frequencies as before
%  업샘플링된 신호니깐 패스밴스스탑밴드 다시설정
Fp = sigBW/2    ;% 외워 sigBandwidth/2 = Passband
Fst = fsamp/2 ;  % fsamp/2 > Stopband

% Construct the RX filter
rxFilt = RxFilt('ovRatio', ovRatio, 'rateIn', fsampUp, 'Fp', Fp, 'Fst', Fst, 'rxScale',false);

% Down-sample and scale the receive signal
r = rxFilt(rup); % 디자인된 필터의 분모랑 신호랑 conv
                 % 원래 임펄스응답함수는 impz
             

% TODO:  Complete the code in WLANRx.detectSTF() as above
% ----------- WLAN Preamble 감지는 두파트로 진행됩니다
% 원래 시작지점 감지는 원래신호와 수신신호를  상관값 계산해서 최대위치를 찾으면되는데
% 여기서 STF 길이는 160샘플이고 주기는16  -> 그럼 160샘플 구간에 16샘플의반복이있나? 보는법(계산량적음)
% 그래서 받은신호를 16샘플지연과 곱해서 크기를구하고 -> 1111배열과 컨볼루션 (이유: 상관값을 누적하는역할) 하고 ->상관값계산
% 원리는 모름 ㅇㅇ
wlanrx = WLANRx();

% % Run the STF detector
% wlanrx.detectSTF(r); % ->원리가뭐지?
wlanrx.detectSTF(r); % wlanrx.pktFound 1찾음 0못찾음

% TODO:  Plot the correlation wlanrx.rhoSTF vs. time in us.  Also, plot the
% location of the detected maximum
figure;
t=0:1/fsamp:1/fsamp*(length(wlanrx.rhoSTF)-1); 
plot(t*1e6,wlanrx.rhoSTF) ;
hold on;
xline(wlanrx.istf/fsamp * 1e6,'r--','detection location');
hold off;

% TODO:  Complete the in WLANRx.pktDetect() to compute the correlation and peak detection with the L-LTF.  
% Test it with the following code.
% STF 로 찾은 위치는 대략적이다 이제 LTF로검증해보자
% LTF는 알려진 LTF샘플로 매치드필터 ->정규화 매치드필터=XCORR(A,B)=CONV(A,CONJ(FLIPUD(B)) 길이가
% XCORR이 길이가 훨씨인길다


% % Create simulation objects
% wlanrx = WLANRx();
wlanrx=WLANRx(); 
% 
% % Perform the LTF correlation
% wlanrx.pktDetect(r); % -  
wlanrx.pktDetect(r);
% % TODO:  Plot the LTF correlation in rx.rhoLTF.  
t=0:1/fsamp:1/fsamp*(length(wlanrx.rhoLTF)-1); 
figure;
plot(t*1e6,wlanrx.rhoLTF) ; % 피크 양옆에 작은피크가 있는데 그 이유는 ltf내에 주기반복되는신호가잇대
hold on;
dlyEst=(wlanrx.iltf -length(wlanrx.lstf))/fsamp * 1e6  ; % 정밀한 ltf 위치 -8마이크로초빼기=dly추정
fprintf("설정된 딜레이 : %f ,STF로추정한 딜레이:%f, LTF로추정한 딜레이 %f",dlyus,wlanrx.istf/fsamp * 1e6,dlyEst); % <?
%---------------------여기서 부터 SNR 에 기반한 반복문-------------------------------------------
%Parameters
ntest = 100;              % SNR하나에 100개테스트
snrTest = [-6:2:4]';      % SNR values to test
                          % 일정 SNR구간에서  패킷검출율급락
nsnr = length(snrTest);

% Initialize arrays
dlyEst = zeros(ntest,nsnr);
snrEst = zeros(ntest,nsnr);
dlyTrue = zeros(ntest,nsnr);
pmd = zeros(nsnr,1);  % 100번 시도중에 패킷을 놓친횟수
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
    for it = 1:ntest %100번테스트

        % TODO:  Generate a random packet using wlantx()
        %    x = ... 
          x=wlantx();

        % TODO:   Upsample with the txFilt   
        %    xup = ...
         xup=txFilt(x);

        % Compute the noise variance
        Etx = mean(abs(xup).^2);
        wvar = Etx*db2pow(-snr);
        chan.wvar = wvar; %chan은 snr-> wvar만 바뀌는걸로


        % TODO:  Pass through a random single path channel using the
        % chan.randChan method
        %     [rup, dly] = chan(...)        
        % chan

      % Run the channel
      [rup,dly] = chan(xup);

        % TODO:  Down-sample
        %    r = rxFilt(...);
        r=rxFilt(rup);

        % Run detection with the rx.pktDetect method.  
        wlanrx.pktDetect(r); % lstf추정

        
        if wlanrx.pktFound % 패킷을 찾을때마다 dlyTrue와 딜레이추정
            i = i + 1;  % Increment counter of number of packets found

            % TODO:  Store:
            %     dlyTrue(i,isnr) = true delay 
            %     dlyEst(i,isnr) = estimated delay based on rx.iltf
            dlyTrue(i,isnr)=dly*1e6;
            dlyEst(i,isnr) = (wlanrx.iltf -length(wlanrx.lstf))/fsamp * 1e6; % Store estimated delay based on rx.iltf
        end
    end

    % Set number of packets found
    nfound(isnr) = i; % snr1개당 패킷을 찾은횟수

    % TODO:  Set     
    %   pmd(isnr) = fraction of packets that were missed.  
    pmd(isnr)=nfound(isnr)/ntest;

    % Print results
    fprintf(1,"SNR=%7.2f MD=%12.4e  \n", snr, pmd(isnr));
end

% TODO:  Plot the missed detection vs. SNR.  You should see that the missed
% detection is low once the SNR is about 2 dB.
figure;
plot(snrTest,pmd);
xlabel('SNR(dB)');ylabel('Packet Detection 확률');
ylim([0 1.2]);
% TODO:  Compute the delay offset
%    dlyOff = ...


% TODO:  Set dlyEstAdj = dlyEst + dlyOff, which is the adjusted delay
% estimate.

% SNR values to plot are snr(Iplot(j))
figure;
Iplot = [3,5];
nplot = length(Iplot);
for iplot = 1:nplot
    subplot(1,nplot,iplot);
    isnr = Iplot(iplot);     % 3번째 snr0 5번째 snr 2
    n = nfound(isnr);

    % TODO:  Create a scatter plot of dlyTrue(:,isnr) vs. dlyEstAdj(:,isnr)
    % Label the axes
    scatter(dlyTrue(:,isnr),dlyEst(:,isnr));
    axis tight; axis equal;       % 축 비율 맞추면 45도선 비교가 쉬움
    grid on;
     hold on;
   % y=x선 그리면 비교 쉬울거같은데
  lims = [min([xlim ylim]) max([xlim ylim])];
    plot(lims, lims, 'r--'); 
  xlabel('True delay [\mus]');
    ylabel('Estimated delay [\mus]');
end


% %----------------SDR 전송파트----------------------------------
% % TODO:  Set parameters
% % Set to true for loopback, else set to false
% loopback = true; % true 1개 false 2개  
% 
% % Select to run TX and RX
% runTx = true;
% runRx = true;
% 
% % Flag indicating if pre-recorded data is to be used
% usePreRecorded = false;
% saveData = false;
% 
% % Compute the signal bandwidth
% info = wlanNonHTOFDMInfo('NonHT-Data');
% sigBW = (max(info.ActiveFFTIndices)-min(info.ActiveFFTIndices)+1)/info.FFTLength*fsamp;
% 
% 
% % Intialize the communication objects
% wlantx = WLANTx('psduLen', 512);
% wlanrx = WLANRx('psduLen', 512);
% fsamp = wlantx.fsamp;
% ovRatio = 2;
% fsampUp = fsamp*ovRatio;
% txFilt = TxFilt('rateIn', fsamp, 'Fp', Fp, 'Fst', Fst, 'ovRatio', ovRatio, 'backoffAutoScale', true, ...
%     'backoffLev', backoffLev);
% rxFilt = RxFilt('ovRatio', ovRatio, 'rateIn', fsampUp, 'Fp', Fp, 'Fst', Fst, 'rxScale', true); % -> true로바꼇넼
% 
% 
% % TODO:  Generate a random packet using wlantx()
% %    x = ... 
% x=wlantx();
% xup=txFilt(x);
% % Set the packet length
% wlanrx.set('pktLen', length(x)); %? 원래패킷의 길이?
% % clear previous instances
% clear devtx devrx
% % Set the parameters.  Note that we capture twice the number of samples
% % to ensure that we always have a full packet in the window
% nsampsFrameTx = length(xup); % 만샘플= 보내는  샘플갯수 =1프레임?
% nsampsFrameRx = 2*length(xup); % 받는샘플갯수는 더많이
% centerFrequency = 2.4e9; % 보내는 주파수 2.4GHZ 
% 
% % Run the creation function.  Skip this if we are using pre-recorded
% % samples
% if ~usePreRecorded
%     [devtx, devrx] = plutoCreateTxRx(createTx = runTx, createRx = runRx, loopback = loopback, ...
%         nsampsFrameTx = nsampsFrameTx, nsampsFrameRx=nsampsFrameRx, sampleRate = fsampUp, ...
%         centerFrequency=centerFrequency);
% 
% end
% 
% 
% if runTx
%     % TODO:  Use the devtx.release() and devtx.transmitRepeat() commands to
%     % continuously send xup
%     devtx.release();
%     devtx.transmitRepeat(xup);
% 
% end
% 
% 
% if usePreRecorded  
%     load singData;
% else
%     % 12비트 fullscale로 조정안하는 이유 -? RxFilt 에서 함
% 
%     % TODO:  Capture data at the upsampled rate.  Capture nsampsFrameRx.
%     %   rup = rx.capture(...)    
% 
%     rup=devrx.capture(nsampsFrameRx);% 프레임당 샘플 명시!
% 
% end
% 
% 
% 
% % TODO:  Down-sample
% %   r = rxFilt(...);
% r=rxFilt(rup);
% % TODO:  Run detection with the rx.pktDetect method.  
% wlanrx.detectSTF(r);
% 
% % TODO:  Plot wlanrx.rhoSTF vs. time in us
% t=0:1/fsamp:1/fsamp*(length(wlanrx.rhoSTF)-1); 
% figure;
% plot(t*1e6,wlanrx.rhoSTF) ; 
% 
% %--------------------Continous Monitoring-----------------------
% % Number of iterations
% nit = 100;
% 
% % Initialize arrays to store data
% rhoMax = zeros(nit,1);                  
% rxTime = zeros(nit,1);
% 
% % Initialize plots
% f=figure;
% title('채널의 딜레이 추정')
% subplot(1,2,1);
% rhoMaxToNow = [0,0];
% timeToNow = [0,3];
% prhoMax = plot(timeToNow,rhoMaxToNow, 'o-', 'LineWidth', 2);
% ylim([0,1]);
% xlabel('Time [sec]');
% ylabel('rhoMax');
% grid on;
% 
% subplot(1,2,2);
% 
% n = length(wlanrx.rhoSTF);
% trho = (0:n-1)/fsamp*1e6;
% rhoSTF = wlanrx.rhoSTF;
% prho = plot(trho, rhoSTF);
% ylim([0,1]);
% xlabel('Time [us]');
% ylabel('rho');
% grid on;
% 
% 
% % Set pointers to X and Y data for both plots
% prhoMax.XDataSource = 'timeToNow';
% prhoMax.YDataSource = 'rhoMaxToNow';
% prho.XDataSource = 'trho';
% prho.YDataSource = 'rhoSTF';
% 
% % Initialize array to save data
% measData = zeros(nsampsFrameRx, nit);
% 
% % Load pre-recorded data if need
% if usePreRecorded
%     load multData;
% 
% end
% 
% for it = 1:nit
% 
% 
%     if usePreRecorded
%         rup = measData(:,it);
% 
%     else
%         % TODO:  Get the data
%         %  rup = devrx.capture(...)
% 
%          rup=devrx.capture(nsampsFrameRx);
%         % Get the time
%         if (it==1)
%             tic; % 시작시간
%         end    
%         rxTime(it) = toc(); %끝난시간
%     end
%     if saveData
%         measData(:,it) = rup;
%     end
% 
% 
% 
%     % TODO:  Down-sample
%     %   r = rxFilt(...);
% 
%     r=rxFilt(rup);
%     % TODO:  Run detection with the wlanrx.pktDetect method.  
%     wlanrx.detectSTF(r);
% 
% 
%     % Get the STF correlation and create the vectors to plot
%     n = length(wlanrx.rhoSTF);
%     trho = (0:n-1)/fsamp*1e6;
%     rhoSTF = wlanrx.rhoSTF;
% 
%     % TODO:  Get the maximum correlation
%     %   rhoMax(it) = ...    
%    rhoMax(it)=max(rhoSTF);wlanrx.detectSTF(r);
% 
%     % Update the plots   
%     timeToNow = rxTime(1:it);
%     rhoMaxToNow = rhoMax(1:it);    
%     refreshdata;
%     drawnow;
%   exportgraphics(f, 'Lab7result.gif', 'Append', true);
% end
% 
