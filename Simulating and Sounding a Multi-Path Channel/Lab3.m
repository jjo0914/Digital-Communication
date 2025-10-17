clc;
clear all;
clf;
% 채널모델 시뮬레이션 코드
% TDL 모델
% 3GPP라는 기관에서 사용하는  신호가 여러경로를 통해오는것을 모방하는 5G모델을 공뷰
% 수식은 h(t)=각채녈 gain * 시간지연델타 들의 합

% ----------------섹션1. 채널모델 비주얼라이징 하기----------------------
%
%일단 다음채널들을 사용하라
tdl=nrTDLChannel('DelayProfile','TDL-A','DelaySpread',1e-7);
%ntTDLChaneel : 5g채널 TDL 모델 객체생성함수
%인자 DelayProfile은 지연모델을 정함.. 예시로 TDL-A,B,C,D 등이잇다.. <이게모임?
%  TDL-A는 23개의 경로
%인자 DelaySpread는  지연rms값(계산수식이있음 지연값과 Gain값이용)을 받음
chaninfo=info(tdl); % info는 객체오브젝트의값들만 받는다는데...
gainPath=chaninfo.AveragePathGains';%각경로당 Gain(dB) (23경로) 
dlyPath=chaninfo.PathDelays'; % 각 경로당 delay (23경로)
% 
%To visualize the channel, plot gainPath vs. the delay.  
%Set the BaseValue to -40 so that x-axis is well below the maximum path
%Put the delay in nanoseconds.
dlyns=dlyPath*10e9;

stem(dlyns,gainPath,'BaseValue',-40); %base가 0부터시작해서 -40dB으로 수정해야함
xlabel('delay(ns)')
ylabel('gainPath(dB)')


%----------섹션2.멀티-패스 채널 시뮬레이팅하기--------------
%지연을 적용하기위해  지연값이 단순히 샘플링시간에 정수배라면..
%배열을 몇칸씩 땡기는걸로 쉽게 지연을 구연할 수 있지만, 1.3T이럴경우
%1.3칸씩 이동하는 인덱스는 존재하지않으니 샘플 중간값을 추정해야한다,
%그러기 위해서 원래 델타함수지만 밴드제한된영역에서는
% 일반적으로 sinc함수를 필터로 사용해서 구현한다고한다.(샘플링타이밍에 0임으로 지연구현가능..?)
%보간 수식들이있긴한데 Matlab에서는 이과정을 쉽개하기 위해 함수를제공함.dsp.VariableFractionalDelay
% Create a fractional delay object 
fracDly = dsp.VariableFractionalDelay(...
    'InterpolationMethod', 'Farrow','FilterLength',8,...
    'FarrowSmallDelayAction','Use off-centered kernel',...
    'MaximumDelay', 1024); %보간함수farrow ,필터함수길이는8=입력샘플도8
           %Use off-centered kernel은 지연이0임에도 샘플이 가중합되는걸 방지하는 거라나머라나
           %헷갈

%Create a ramp signal
nt = 100; %100개샘플
x = (0:nt-1)'/nt; % 값이 0부터 0.01씩늘어나는 그냥 직선함수.
dlySamp = 10.7; %10.7샘플 지연

y = fracDly(x, dlySamp); %원래함수가 100개면 지연함수도 그냥100개네 100+지연샘플이아니라..
figure(2)
plot(x)
hold on
plot(y)
legend('Oringinal 함수','dlySamp만큼 샘플 지연된 함수')

% SISOchan.m 함수 채우기.. dlyPath, fsamp, gainPath 을 값으로 보냄
% 이제 신호를 TDL-A모델링하고 + 시간지연을 함수로 재현하는거 (SISOchan)
fsamp=30.72e6*2; %60MHz 수신기 샘플링속도
%1샘플당 16.2나노초
chan=SISOChan(dlyPath=dlyPath,gainPath=gainPath,fsamp=fsamp);  % 클래스내에 변수값입룍 obj.으로 dlyPath,gainPath,fsamp를 입력으로받을수잇음
                                                               % 클래스내변수는 dlyPath,gainPath,fsamp,phasePath가있음
% Create a vector of times at the sample rate.  We start a little before
% zero since the channel filters take some initial zeros.
% 필터초기에는 입력앞에 0을 패딩하는게 일반적이래 왜인지 알거같았는데
tmin= -0.5e-6;
tmax= 2e-6;
t=(tmin:1/fsamp:tmax)'; %수신기측에서 시간
nt=length(t);
% TODO:  Create a signal ximp of length nt representing an impulse at zero.
% Specifically, we want 
%      ximp(i0) = 1 at the sample i0 where t(i0) ~= 0 
%      ximp(i) = 0 for all other i ~= i0
%
%    ximp = ...
%  델타함수=임펄스함수? 
ximp=zeros(nt,1); % 샘플링시간만큼 배열만들고 시간0에서 1
[~,i0]=min(abs(t-0));
ximp(i0)=1;
figure(3);
plot(t,ximp);
hold on;
%으에..
% TODO:  Get the impulse response by sending ximp through the channel
%   yimp = chan(...);
yimp=chan(ximp);
chan.release;

% Compute the energy in each output sample, we limit the smallest energy to
% avoid taking log of zero
% SISOChan에서 db를 진폭으로 바꿔서 y를 계산하긴햇어
ypow = pow2db(max(abs(yimp).^2, 1e-8)); %abs: 복소수 크기변환
                                        %max( ,1e-8) 1e-8 이하값은 강제로
                                        %1e-8로바꿔버림 이거안하면 이상함값이
%ypow는 "채널이 입력 신호를 어느 시간 지점에서 얼마나 강하게 통과시키는가"를 dB로 나타낸 시간 에너지 분포입니다.

% TODO:  Plot ypow vs. time using the stem plot as before.
% On the sample plot, plot the channel impulse response.
plot(t,abs(yimp))
legend('지연X 채널h','경로지연모델 채널 h')
xlabel('time(s)')
ylabel('신호 크기')
figure(4);
stem(t,ypow,'BaseValue',-80)
ylabel('dB');
title('채널이 입력 신호를 어느 시간 지점에서 얼마나 강하게 통과시키는가"');
% 위에서 ximp를 보내서 응답을 보고 채널을 관측할 수 있었지만,
% % 실제로 ximp신호를 보내는 것은 비효율적이다 (짧은 시간에 너무 많은 에너지를 내보내야 하기 때문입니다.)
% % 그래서 시간영역대신 주파수 영역을 관찰하는데.(imp 함수 푸리에변환은 1)
% % 현실에서는 임펄스 대신 **넓은 주파수 대역을 커버하는 신호(예: OFDM, PN 시퀀스)**를 사용하여 채널 특성을 추정합니다.
% 
% %이 방법에서는 신호를 먼저 주파수도메인으로 만듬(심볼이 주파수도메인) 
% %그다음에 ifft를 통해신호를 만들고 전송?..
% 
% Parameters
nfft = 1024;
fsamp = 30.72e6*2;

%난수고정을 한다는데왜?
rng(1,'twister') % twister:난수 생성알고리즘 중 하나

% TODO:  Create a vector, xfd, or nfft random QPSK symbols.
%   xfd = ...
bit=randi([0 1],nfft*2,1); %몇개를 생성해야하지?..몰루,,,
xfd=qammod(bit,4,'InputType','bit','UnitAveragePower',false); %nfft개의 주파수도메인심볼
%save txData xfd;
x=ifft(xfd);
% %이게뭐지 그냥심볼생성하고 ifft하는게?,.
% 
% %이제 다음은 Unit02에서 생성한 txfilt사용하는데
% %디지털심볼을 오버샘플링하고 그ㅡ심볼을를 원하는 주파수 & 원하는 채널대역폭으로
% %가진 필터계수들과 콘볼루션?해서 신호를 만듬
% %그리고 backoff를받아서 전력을 좀낮춤..
% 
% % Create the TX filter object
backoffLev = 12;
txFilt = TxFilt(ovRatio=1, rateIn=fsamp, backoffLev = backoffLev, backoffAutoScale = true);

% Filter the signal
x_txFilt = txFilt(x);
% 
%이제신호를 4번반복하겠다는데?.
%  배열크기가 내가맞는거같은데..
nrep=4;
xrep=repmat(x_txFilt,4,1); 

% 이제 이신호를 TDL-A ,다중경로지연채널모델에넣음
chan = SISOChan(dlyPath = dlyPath, fsamp=fsamp, gainPath=gainPath);
r0=chan(xrep);

%그다음 다중경로 지연채널에 30dB낮은 신호를 넣는데 이유는나중에
snr=30; % -30db=10log(목표신호^2/원신호^2)
wvar=mean(abs(r0).^2)*db2pow(-snr) %(목표신호^2 전력이나옴)
n=length(r0);
r=r0+sqrt(wvar/2)*(randn(n,1) + 1i*randn(n,1));
%randn(n,1) + 1i*randn(n,1) 이식은 각각 평균0 분산1 가우시안잡음을이며
%전력은 2래(각각 전력 1 (=분산))
%여기다  원하는 전력으로 맞추려면@
%가우시안잡음이 실수부 복소수로 나눠있어서
% x+yj 에서 각 x y앞에 있는 계수^2 을 더하면 전력이될거임 그걸수식으로맞춰야함

% %이제 어찌되었든 4번 반복된 신호 xrep와 채널h  그걸수신한신호r 이있음
% % 우리는 이제 채널h의 응답함수가 궁금하다
% % 주파수 도메인에서 계산하면 채널h의 응답함수는 x와r을 단순히나눗셈하면됨
% % 이때 순환 콘볼루션 개념을 사용하는데.. 신호x가 반복전송되니
% % ???? 코드해보자
% % Extract one FFT period of the data
r1 = r(nfft+1:2*nfft);%수신신호r을 nfft 2번째주기부터 nfft갯수만큼 추출하는데..
                     %이유몰루.. 1024가 총 심볼갯수인데..
                     %추출위치 인덱스는 심볼시작지점 이랑일치해야한다..!
%           %인덱스가 졍렬되지않으면 hfd=rfd/xfd가 성립x래... 무엥?..
% % TODO:  Complete the TODO sections in estChanResp up to hfd = ...
% % Then run:
% %     hfd = estChanResp(r1,xfd);
[hfd,h]=estChanResp(r1,xfd,'normToMean',false);
% %hfd는 주파스응답함수 h는 시간영역응답함수
% %hfd를 이용해서 주파수 대역에서 dB구하라
% % TODO:  Compute the channel power gain and shift it using fftshift and
% % plot against the frequency in MHz.  Label your aees.
% %   hfddB = ...
% %   fMHz = ...
% %   plot(...)
hfddB=pow2db(max(abs(hfd).^2,1e-8));% power gain의 표준 표현 db아님!!
fMHz=(-nfft/2:nfft/2-1)*fsamp/nfft;
figure(5);
plot(fMHz/10e6,fftshift(hfddB))
title('H=주파수영역에서 채널추정(넓은주파수대역신호)')
xlabel('MHz')
hpow=pow2db(max(abs(h).^2,1e-8));
tns=(1:length(hpow))/fsamp;
figure(6);
plot(tns,hpow);
title('시간영역 채널 추정')
ylabel('dB')
xlabel('seconds')

%그래서 위 주파수도매인 심볼을 통해서
%채널추정을 하면 어느 시간대에서 gain이 큰지 추정할수있지만
%실제로 과연 맞았을까하는? 코드
[Pm,im]=max(gainPath); %estChanResp함수에서 h를 가장큰 gain에 위치를 맞추었지, 왼쪽에서8번째..
realDlyTrue=dlyPath-dlyPath(im); % 시간값들 배열을 가장큰 gain이 나온 시간으로 빼면 결국 시간이 왼쪽으로 shift되는거자
gainTrue=gainPath-Pm; % gain을 왜 최대값 기준으로 다빼지?
%얻은 채널h 임펄스 함수도 똑같이 스케일링
[Pm,im]=max(hpow);
realDlyMeasure=tns-tns(im);
gainTrueMeasure=hpow-Pm;

figure(7);
plot(realDlyTrue,gainTrue,'DisplayName','True');
hold on;
plot(realDlyMeasure,gainTrueMeasure,'DisplayName','True');
legend('실제 채널 gain','채널 추정 gain')
title('gain값 비교')
%와 이게 비슷하게 나오네?.


%%---ADAM PLUTO 실습 부분!
%%%-------------이제 sdr로 real channel을 추정해보자?
% 채널임펄스응답합수를 추정하기위한 주파수도메인신호를만듬!
% Parameters
nfft = 1024;
fsamp = 30.72e6*2;

%난수고정을 한다는데왜?
rng(1,'twister') % twister:난수 생성알고리즘 중 하나

% TODO:  Create a vector, xfd, or nfft random QPSK symbols.
%   xfd = ...
bit=randi([0 1],nfft*2,1); %몇개를 생성해야하지?..몰루,,,
xfd=qammod(bit,4,'InputType','bit','UnitAveragePower',false); %nfft개의 주파수도메인심볼
save txData xfd;
x=ifft(xfd);
%이게뭐지 그냥심볼생성하고 ifft하는게?,.

%이제 다음은 Unit02에서 생성한 txfilt사용하는데
%디지털심볼을 오버샘플링하고 그ㅡ심볼을를 원하는 주파수 & 원하는 채널대역폭으로
%가진 필터계수들과 콘볼루션?해서 신호를 만듬
%그리고 backoff를받아서 전력을 좀낮춤..

% Create the TX filter object
backoffLev = 12;
txFilt = TxFilt(ovRatio=1, rateIn=fsamp, backoffLev = backoffLev, backoffAutoScale = true);

% Filter the signal
x_txFilt = txFilt(x);

%이제신호를 4번반복하겠다는데?.
%  배열크기가 내가맞는거같은데..
nrep=4;
xrep=repmat(x_txFilt,4,1); 
% TODO:  Set parameters
% Set to true for loopback, else set to false
loopback=true;

%select to run Tx and Rx;
runTx=true;
runRx=true;
saveData=true;
usePreRecorded=false;
fsamp=30.72e6*2;

%clear previous instance;
clear tx rx

if ~usePreRecorded
    [tx,rx]=plutoCreateTxRx(createTx=runTx,createRx=runRx,loopback=loopback,...
        nsampsFrameRx=nfft,nsampsFrameTx=nfft,sampleRate=fsamp);
end

if runTx && ~usePreRecorded
     % TODO:  Use the tx.release and tx.transmitRepeat commands to
    % continuously send xrep
    tx.release();
    tx.transmitRepeat(xrep); %이러면 채널응답함수를 추정할때 심볼시작위치를 어떻게 추정하지.?
end
if ~runRx
    return;
end

%수신할때 single왜했지 fullscale왜했지?
nbits=12;
fullScale=2^(nbits-1);


if ~usePreRecorded
    
    %rx.set('GainSource','Manual');
    %rx.set('Gain',0);  % 추정채널h가 +20dB나오길래 원래 0dB나와야정상아닌가해서 수신게인0으로절
    % TODO:  Capture data for nfft samples
    %   r = rx.capture(...)
    r=rx.capture(nfft); %값이 -1024~1024사이 11비트
    % TODO:  Scale to floating point
    %   r = ...
    r=single(r)/fullScale; % 원래int16값인데 소수점을 위해 32bit single바꿈


    % Save the data
     if saveData
            save rxDataSing r;
     end

else 
    % Load pre-recorded data
    load rxDataSing;

end
[hfd,h]=estChanResp(r,xfd,'normToMean',true);  % h를 노말라이즈안하면..?

hpow=pow2db(abs(h).^2);
% TODO:  Plot hpow vs. time in ns
%   tns = ...
%   plot(tns, ...);
   
tns=(1:1024)/fsamp;
figure;
plot(tns*1e6,hpow);
xlabel('microSeconds');
ylabel('dB')
title(' SDR 루프백 채널 h 추정') %왜 dB값이 그냥int값이나 single값이나 정규화값이나 차이가없지?..

%%-------이제 실기간 모니터링 하는방법----%%%
f=figure;
%보내는거는 그대로
if runTx && ~usePreRecorded
     % TODO:  Use the tx.release and tx.transmitRepeat commands to
    % continuously send xrepj
    tx.release
    tx.transmitRepeat(xrep); %이러면 채널응답함수를 추정할때 심볼시작위치를 어떻게 추정하지.?
end
if ~runRx
    return;
end
ncaptures=100;
%실시간으로 수신하기
if usePreRecorded
    load rxDatMult; %변수이름이 rdat임
    ncaptures=size(rdat,2) % 1:행 2:열크기
else
    rdat=zeros(nfft,ncaptures);
end
% 초기값 캡처
if usePreRecorded
    r=rdat(:,1)
    
else
    r=rx.capture(nfft);
     r=single(r)/fullScale;
end
% TODO:  Get the channel response 
[hfd,h]=estChanResp(r,xfd,'normToMean',true);% normtomean true하니깐 좀더작아짐
hpow1=pow2db(abs(h).^2);


tus1=tns*1e6;
p=plot(tus1,hpow1);
xlim([min(tus1),max(tus1)]);
ylim([-60, 100]);
xlabel('Time [us]');
ylabel('Path gain [dB]');
title(' SDR 루프백 채널 h 추정');
grid on;
% y값 변하면 새로플로팅
p.YDataSource='hpow1'; %단순변수명?

snr=zeros(ncaptures,1); % snr은또뭐임?
for t=1:ncaptures
    if ~usePreRecorded
        r=rx.capture(nfft);
        r=single(r)/fullScale;
        if saveData
            rdat(:,t)=r;
        end
    else
        r=rdat(:,t);
    end
     % TODO:  Re-compute the impulse response    
     [~,h]=estChanResp(r,xfd,'normToMean',true); %norm to mean true하니간 좀더작아짐.
     % TODO:  Re-compute hpow1 = impulse response on the first nplot samples.
     hpow1=pow2db(abs(h).^2);
         % Redraw plot
    refreshdata(p); % p의 xdatasource 혹은 ydatasource가 변하면 p를 업데이트
    drawnow %그린다.
    exportgraphics(f, 'lab3result3.gif', 'Append', true);
end

% Save data
if saveData
    save rxDatMult rdat;
end





