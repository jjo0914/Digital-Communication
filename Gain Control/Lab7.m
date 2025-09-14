% LAB 5. GAin control
% Gain 제어  수학모델  SNR=G/(G/SNR1 + 1/SNR2)

snr1 = 30;  % Ex/wvar1 in firststage  in dB
snr2 = 10;  % Ex/wvar2 in second stage in dB
rxGainTest = linspace(0,40,100)';  % gain in dB

% TODO: 
%   snr = ...
%   plot(...);
snr=rxGainTest./(rxGainTest/snr1 + 1/snr2);  % Ex/(gain*wvar+wvar2)
plot(rxGainTest,snr);
title('이론 Gain값 선형 구간 비 선형구간')
legend('snr1=30dB , snr2=10dB')
% Gain을 무한대로늘려면 snr은 snr1으로 가까이간다. (비선형영역)
% 그래서 gain을 선형영역에서 조절하는게중요하다

% 그리고 또한 아날로그 수신기는 입력 한계가있다. 신호가 잘려 왜곡될 수있어 gain을 잘조절해라
% 위 모델을 위해 함수sat = @(x,a)  max(-a,min(a,real(x))) + 1i*max(-a,min(a,imag(x))); 를사용하라

% 이하 파라미터들 복붙
% Parameters
saveData = true;
usePreRecorded = false;

% Signal parameters
nfft = 1024;   % signal length
backoff = 9;   % backoff in dB  (포화구간 지점에서 얼마나 낮게운용할것인지)
fsamp = 30.72e6*2;  % Sample rate in Hz
nsampsFrame = nfft; % Signal length or 1프레임당 1024신호


% Set the random generator to ensure that the TX and RX use the same data
rng(0,'twister');  % 난수발생 고정

if usePreRecorded
    load txData;
else
    % TODO:  Create a vector, xfd, of nfft random QPSK symbols.
    %   xfd = ...
    bits=randi([0 1],nsampsFrame*2,1);
    xfd=qammod(bits,4,'InputType','bit','UnitAveragePower',true);

    % Save frequency-domain TX data
    if saveData
        save txData xfd;
    end
end

% TODO: Take IFFT
%   x = ifft(...);
% 심볼을 fft도아니고 ifft? (OFDM이래)
x=ifft(xfd,nfft);
figure(2)
plot(abs(xfd));
legend('이 샘플범위내에서(주파수영역)을 값을 차지하는 심볼들모임?')
figure(3);
plot(abs(x));
legend('ofdm 1심볼 시간영역')

% TODO:  Rescale and clip values for the backoff.  You can use the
% saturation function sat above
%   x = sat(x,1); 백오프만큼 신호를 낮추고 플리핑해라
x=db2mag(-backoff).*x;
sat = @(x,a)  max(-a,min(a,real(x))) + 1i*max(-a,min(a,imag(x)));
x=sat(x,1);

% 이제 만든 4-QPSK ofdm신호를 Gain control하겠다
wvar1 = 1;   % power평균= E|w1|^2 in linear scale -> 복소가우시안신호의 분산=N0
snr1 = 30;   % Erx / E|w1|^2 in dB
snr2 = 10;   % Erx / E|w2|^2 in dB
gainRx = 10; % RX gain in dB

% wvar2를 게산하라
% TODO:   
%   Erx = ...
%   wvar2 = ...
Erx=db2pow(snr1);
wvar2=Erx/db2pow(snr2);

% r=h*x + 잡음
% r이 잡  음이 없다고하면 Erx/Ex=channelgain (Erx=는 신호제곱해서 더한거라가정)
% TODO:  
%   chanGain = ...
%   r = chanGain * x;  
Ex=sum(abs(x).^2); % -> 이거 /N인것같은데 수정: 원래sum(abs(x).^2);:
chanGain=Erx/Ex;  % 에너지 단위네?
r=sqrt(chanGain).*x;
% TODO:
%    w1 = ... 
%    w2 = ...
%    y = db2mag(gainRx)*(r + w1) + w2; 
w1=1/sqrt(2).*(randn(nsampsFrame,1)+1i.*(randn(nsampsFrame,1)));
w2=sqrt(wvar2)/sqrt(2).*(randn(nsampsFrame,1)+1i.*(randn(nsampsFrame,1)));
y = db2mag(gainRx)*(r + w1) + w2; 

% 함수 estChanResp를 쓰면 h를 (r=h*x + w ) 얻을 수 있는데
% h를 이용해서 snr 계산하는게가능함 (Ex/Ew)
% TODO:
%   [~,~,snr] =  estChanResp(y,xfd,'normToNoise',true);
%   snrTheory = ...
%   fprintf(...);  % Print snr and snrTheory

% 여기서는 y=root(G)*(h*x+w1) + w2 를 보내고 원신호x를보내니
% SNR= Gh Ex/ GEw1 + Ew2 가 나올거임 
% 이론 SNR= GEr/GEw1 + Ew2 이랑 비슷한가? x->r 갈때 잡음이없어서 똑같네?ㅇㅅㅇ
[~,h,snr]=estChanResp(y,xfd,'normToNoise',true);
snrtheory=snr1/(1+snr1/gainRx/snr2);
fprintf("snrTheory: %f snr: %f \n",snrtheory,snr); % 오 값비슷하노무현

%%%=---- 자이제 여러 게인을 테스트 해야겠지?..------------
% Gain이낮을수록 이론값이랑 안맞을거래
% Range of gain values to test
gainRxTest = linspace(-10,40,25)'; % -10 40을 25분할... (40 - (-10) /25 =2 칸)
ngain = length(gainRxTest);
ntrials = 100;      % number of trials at each gain level 이건뭐지

% Initialize vectors
snr = zeros(ngain,1);        % Median SNR measured at each gain level
snrTheory = zeros(ngain,1);  % Theoretical SNR with no non-linearity

% Loop over gain levels
for i = 1:ngain

    % Get gain to test
    gainRx = gainRxTest(i);
    snrIt = zeros(ntrials,1);

    for it = 1:ntrials % snr값을 100번계산해서 중간값추출?
        
        % TODO:  Generate random noises and output y
        %    w1 = ...
        %    w2 = ...
        %    y = ...
       w1=1/sqrt(2).*(randn(nsampsFrame,1)+1i.*(randn(nsampsFrame,1)));
       w2=sqrt(wvar2)/sqrt(2).*(randn(nsampsFrame,1)+1i.*(randn(nsampsFrame,1)));
       y = db2mag(gainRx)*(r + w1) + w2; 
        % TODO:  Estimate snr and store it in snrIt
        %   [~,~,snrIt(it)] =  estChanResp(...);
        [~,~,snrIt(it)]=estChanResp(y,xfd,'normToNoise',true);
    end

    % Save the median SNR over the trials
    %    snr(i) = ...   
    snr(i)=median(snrIt);
    
    % TODO:  Measure the theoretical SNR
    %    snrTheory(i) = ...
    snrTheory(i)=snr1/(1+snr1/gainRx/snr2);
end

% TODO:  Plot the measured and theoretical SNR
figure(5);
plot(gainRxTest,snr);
hold on
plot(gainRxTest,snrTheory);
xlabel('gain(dB)')
ylabel('SNR')
legend('추정한 SNR값','이론 SNR값 ');
% Gain이 낮을수록 지랄나는 경우를 볼수 있는데 그 이유는
% estChanResp에서 피크를 오검출하기때문이다.(SNR과대추정..)
% (잡음 구간을 신호라고 판단)

% Gain커질수록 비선형 구간이 생기는걸 없애기 위해
% Gain이 커질수록 잡음의 영향을 줄이겟다(비선형 요인은 잡음의 영향이커져서 ㅇㅇ)
% ------------------------ 잡음의 Gain값 자르기 ---------------
% Range of gain values to test
gainRxTest = linspace(-10,40,25)';
ngain = length(gainRxTest);
ntrials = 100;      % number of trials at each gain level

% Nonlinear parameters
satLev = 15;  %  15dB Gain넘으면 잡읍을 컷 Saturation level in dB above noise
satLevLin = sqrt( db2pow(satLev)*wvar1 ); % 15dB이 넘을때의 잡음의 크기를 포화레벨로설정

% Initialize vectors
snrNL = zeros(ngain,1);      % Median SNR with non-linearity

% Loop over gain levels
for i = 1:ngain

    % Get gain to test
    gainRx = gainRxTest(i);
    snrIt = zeros(ntrials,1);

    for it = 1:ntrials
        
        % TODO:  Generate random noises and output with the non-linearity
        %    w1 = ...
        %    w2 = ...
        %    u1 = ...
        %    u2 = sat(u1,satLevLin);     
        %    y = ...
       w1=1/sqrt(2).*(randn(nsampsFrame,1)+1i.*(randn(nsampsFrame,1)));
       w2=sqrt(wvar2)/sqrt(2).*(randn(nsampsFrame,1)+1i.*(randn(nsampsFrame,1)));
        u1=db2mag(gainRx).*w1; 
        u2=sat(u1,satLevLin); % 잡읍만컷트?
       y = db2mag(gainRx).*r+u2+ w2; %   y = db2mag(gainRx)*(r + w1) + w2;
        % TODO:  Estimate snr and store it in snrIt
        %   [~,~,snrIt(it)] =  estChanResp(...);
         [~,~,snrIt(it)]=estChanResp(y,xfd,'normToNoise',true);
    end

    % Save the median SNR over the trials
    %    snrNL(i) = ...
   snrNL(i) =median(snrIt);

end

% TODO:  Plot the measured SNR with and without the non-linearity and
% theoretical SNR.  Add a legend to your plot so that you can see the
% different curves
figure(6);
plot(gainRxTest,snrNL);
hold on
plot(gainRxTest,snrTheory);
xlabel('gain(dB)')
ylabel('SNR')
legend('잡음만 자르고 추정한 SNR값','이론 SNR값 ');
% 오진짜 잡음만 자르니깐 점점 선형이되네?

%-------------이제 SDR로 실습-0--------------------------
% 복붙
% TODO:  Set parameters
% Set to true for loopback, else set to false
loopback = false;  

% Select to run TX and RX
runTx = true;
runRx = true;
% clear previous instances
clear tx rx
% addpath('..\common'); % 혹시 createTxRx가 다른 디렉토리에 있으면 추가하라고
% Run the creation function.  Skip this if we are using pre-recorded
% samples
if ~usePreRecorded
    [tx, rx] = plutoCreateTxRx(createTx = runTx, createRx = runRx, loopback = loopback, ...
        nsampsFrame = nsampsFrame, sampleRate = fsamp);
end

if runTx && ~usePreRecorded
    % TODO:  Use the tx.release() and tx.transmitRepeat() commands to
    % continuously send x
    tx.release();
    tx.transmitRepeat(x);% OFDM=4QPSK심볼만들고 IFFT
end
% If not running the RX, stop the live script now
if ~runRx
    return;
end


% PLUTO SDR에는 훌륭한 AGC가있지만 여기서는 끄고 수동으로 게인을 컨트롤을 해보겠다
% Gain range in dB
gainRxMin = -3;
gainRxMax = 71;

% Place the Pluto in manual gain mode.  When you do this the first time,
% you will need to release it
if ~usePreRecorded
    rx.release();
    rx.set('GainSource', 'Manual');
end

% Manual gain setting.  You may need to change this value to get a good response in the next step
if loopback
    gainRx = 20;      % 루프백이면 gain작게
else
    gainRx = 50;
end
if ~usePreRecorded
    rx.set('Gain', gainRx);
end

nbits = 12;               % number of ADC bits(2^12= 01111111111(12개)=2047,10000000000=-2048) 
fullScale = 2^(nbits-1);  % full scale 2048

if ~usePreRecorded

    % TODO:  Capture data
    %   y = rx.capture(...)
    y=rx.capture(nfft);
    % TODO:  Scale to floating point
    %   y = ...
    y=single(y)./fullscale;   % double형태가 아니라 single형태로 해야할지도?
    % Save the data 
    if saveData
        save rxDataSing y;
    end
  
else 
    % Load pre-recorded data
    load rxDataSing;

end

% TODO:  Measure the channel response and SNR with the estChanResp function
[channelresponse,SNR]=estChanResp(y,xfd,'normToNoise',true);
fprintf("수신Gain20일때 snr : %f",SNR);
% TODO:  Plot the channel response
figure(7);
plot(channelresponse);
legend("SDR 숫신 gain20일때 채널 response")

%%-----------수신 게인  테스트------------
% Rx gain levels to test
gainRxTest = linspace(0,60,7)'; % Gain 0부터 60 (7포인트=6구간,0 10 20 30),,
ngain = length(gainRxTest);

% Number of trials per gain step
ntrials = 10;

% SNR and RX signal energy per sample in trial and gain
snr = zeros(ngain,ntrials);
ypow = zeros(ngain,ntrials);

% Get pre-recorded data if requested
if usePreRecorded
    load rxDatGain;
else
    ydat = zeros(nsampsFrame,ngain,ntrials);
end

snrMed=zeros(ngain,1);
powMed=zeros(ngain,1);
rxbackoff=zeros(ngain,1);
% Loop over gain levels
for i = 1:ngain

    % Get gain to test
    gainRx = gainRxTest(i);

    % Set the gain
    if ~usePreRecorded        
        rx.set('Gain', gainRx);
    end
    

    % Loop over trials
    snrIt = zeros(ntrials,1);   
    for it = 1:ntrials
                       
        if usePreRecorded
            y = ydat(:,i,it);
        else
            % TODO:  Get the data and convert to floating point
            y = rx.capture(nfft);
            y = single(y)./fullScale;
        end

        % Save the data
        ydat(:,i,it) = y; %캡처할때마 새로운 잡음이 섞이니 SNR다다르지

        % TODO:  Estimate snr and store it snr(i,it)
        %   [~,~,snr(igain,it)] =  estChanResp(...);
        [~,~,snr(i,it)]=estChanResp(y,xfd,'normToNoise',true); % snr ngain*ntrials
           
        % TODO: Measure the energy per sample in y and store it in ypow
        %   ypow(i,it) = ...
        ypow(i,it)=mean(abs(y).^2); % |y|.^2 순간전력 / 다더한걸 평균내면 평균전력=샘플당 에너지?
                                    % 다더하면 총에너지  ㅝ지..?( 총전력이란 말은 없어?!=전력은=가속도,에너지는속도)
    end

    % Print the gain and median SNR
    fprintf(1,'gainRx = %.2f snr=%.2f\n', gainRx, median(snr(i,:)));
    snrMed(i)=median(snr(i,:));
    powMed(i)=median(ypow(i,:));
    rxbackoff(i)=-pow2db(powMed(i)/2); % sdr최대 평균 전력은 2 (real=1,imag=1)
                                    % 근데 Gain수집된신호가 최대전력에비해 몇dB 낮은지?
                                    % 어느 Gain부터 포화되는 db값이 있을거임
end
% Save the data
if saveData
    save rxDatGain ydat;
end
% 
figure(8);

subplot(1,2,1);
% TODO:  Plot snrMed vs. gainRxTest.  Label your axes
plot(gainRxTest,snrMed)
xlabel('sdr Gain 설정')
ylabel('측정된 snr값')
subplot(1,2,2);
% TODO:  Plot rxBackoffMed vs. gainRxTest.  Label your axes
plot(gainRxText,rxbackoff)
xlabel('sdr Gain 설정')
ylabel('Gain적용된신호가 최대전력(2)에비해 몇db낮은지?');


%% d이제 rxbackoff가 12db이되게 Gain설정을(AGC) 코드로 짜라는거
% 알빠노?싶
% Number of iterations
nit = 100;

% RX backoff target in dB
rxBackoffTgt = 12;

% Initialize arrays to store data
snr = zeros(nit,1);             
rxBackoff = zeros(nit,1);       
gainRx = zeros(nit,1);
rxTime = zeros(nit,1);
rxPow = zeros(nit,1);

% Get pre-recorded data, if requested
if usePreRecorded
    load rxDatContinuous;
else
    ydat = zeros(nsampsFrame,nit);

    rx.release();
    rx.set('GainSource', 'Manual');
end

% Initialize RX gain
gainInit = 0;
gainRx(it) = gainInit;

% Initialize plots
figure(9);
subplot(1,2,1);
snrToNow = [0,0];
timeToNow = [0,3];
psnr = plot(timeToNow,snrToNow, 'o-', 'LineWidth', 2); %플롯생성
ylim([0,60]);
xlabel('Time [sec]');
ylabel('SNR[dB]');
grid on;

subplot(1,2,2);
rxPowToNow = [-50,50];
ppow = plot(timeToNow,rxPowToNow, 'o-', 'LineWidth', 2); %플롯생성
%ylim([-80,0]);
xlabel('Time [sec]');
ylabel('RX power [dB]');
grid on;


% Set pointers to X and Y data for both plots
psnr.XDataSource = 'timeToNow';
psnr.YDataSource = 'snrToNow';
ppow.XDataSource = 'timeToNow';
ppow.YDataSource = 'rxPowToNow';

for it = 1:nit
    
    if usePreRecorded
        y = ydat(:,it);
    else
        % TODO:  Set the gain to gainRx(it)
        %   rx.set(...); 'Gain',숫자

        % TODO:  Get the data and convert to floating point
        %   y = ...
    end    
    ydat(:,it) = y;

    % Get the time
    if (it==1)
        tic;  %시간측정함수
    end    
    rxTime(it) = toc(); %시간측정종료

    % TODO:  Estimate the snr 
    %   [~,~,snr(it)] =  estChanResp(...);

    % TODO:  Estimate the energy per sample, RX backoff and RX power
    %    ypow = ...
    %    rxBackoff(it) = ...
    %    rxPow(it) = ...
    
    if (it < nit-1)
    
        % TODO: Compute the gain step.  
        %   step = ...

        % Limit the step size
        stepMax = 10;
        step = max(min(stepMax, step), -stepMax);

        % TODO:  Update the gain for the next iteration
        %    gainRx(it+1) = ...
        gainRx(it+1) = gainRx(it) + step;
        
        % Limit the gain to the max and minimum values and ensure the value
        % an integer
        gainRx(it+1) = max(gainRxMin, gainRx(it+1));
        gainRx(it+1) = min(gainRxMax, gainRx(it+1));
        gainRx(it+1) = round(gainRx(it+1));
    end
     
    % Update the plots
    snrToNow = snr(1:it);
    timeToNow = rxTime(1:it);
    rxPowToNow = rxPow(1:it);
    refreshdata;
    drawnow;

end

if saveData
    save rxDatContinuous ydat;
end