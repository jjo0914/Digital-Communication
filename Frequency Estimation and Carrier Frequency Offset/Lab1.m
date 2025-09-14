clc; clear all;
%%%%%%%%%% UG버전 굳이 안해도될듯?..%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1. Creating a Complex Exponential in Time-Domain
% 
% nsamp = 4096;       % number of samples
% fsamp = 30.72e6;    % 30.72MHz sample frequency in Hz 
% f0 = 35e3;          % 35kHZ frequency in Hz
% A = 1;              % complex gain
% 
% % TODO:  Create discrete-time signal
% %   t = ...  시간도메인
% %   x = ...  디지털샘플링
% %t=A
% time=linspace(0,1.3333333333333333333333333333333e-4,4096);
% x=A*exp(1i*2*pi*(f0/fsamp)*linspace(0,nsamp-1,nsamp));  % 
% % 1/fs =t 잖아 1/fs*samp=t , 이면 t=1.3333333333333333333333333333333e-4 초 동안
% t=A*exp(1i*2*pi*f0*linspace(0,1.3333333333333333333333333333333e-4,4096));
% figure(1);
% plot(linspace(0,1.3333333333333333333333333333333e-4,4096),real(x))
% hold on
% plot(linspace(0,1.3333333333333333333333333333333e-4,4096),imag(x))
% legend('real','imag')
% xlabel('time(s)')
% 
% %Visualizing the Signal in Frequency-Domain
% Xf = fft(x)/nsamp;  % 샘플로 나누노
% Xf = fftshift(Xf);
% f = (-nsamp/2:nsamp/2-1)/nsamp*fsamp; % HZ단위
% %PSD : fft값을 dBM으로
% Xfpow = mag2db( abs(Xf) ); % dB/sample
% figure(2);
% plot(f/1e6,Xfpow);
% xlabel('MHz');ylabel('dB');
% 
% %Estimating the Frequency via a Peak in the FFT
% %스펙트럼에서 피크주파수를 찾아서 캐리어주파수를 추정하라 
% [~,idx]=max(Xfpow);  
% f0est=f(idx);
% fprintf("추정 캐리어 주파수 %.2f kHZ /",f0est/1e3)
% fprintf("실제 캐리어 주파수 %.2f kHZ \n",f0/1e3)
% fres=fsamp/nsamp/2; % 음..어떻게계산하는? fs/(N)은 간격, 2로 나눠서 +- 구역
% fprintf("에러가 날수있는 허용오차 +- %.2f kHZ \n",fres/1e3)
%% 


%%%%%%%%%%%%%%%%%%%SDR로 전송하기.. %%%%%%%%%%%%%%%%%%%%5
% TODO:  Set parameters
% Set to true for loopback, else set to false
loopback = true;   %그냥논리값 배졍 true:기기1개 ,false:기기2개

% Select to run TX and RX
runTx = true;  %그냥 논리값 배정
runRx = true;  %그냥 논리값 배정
usePreRecorded = false;  % SDR없는 사람들은 TRUE

% clear previous instances
clear devtx devrx

% add path to the common directory where the function is found
addpath('C:\Users\juwon\Desktop\Digital Communication\Unit01_Lab');  %나중에 estFreq함수쓸건데 함수 위치 추가

% Parameters

sampleRate = 30.72e6; %샘플링 주파수 30.72MHz
nsampsFrame = 2^12;    %프레임당 샘플
fc = 2.4e9;  %  캐리어주파수or 중심주파수2.4GHZ

% Run the creation function.  Skip this if we are using pre-recorded
% samples
if ~usePreRecorded %FALSE면
   % [devtx, devrx] = plutoCreateTxRx(createTx = runTx, createRx = runRx, loopback = loopback, ...
   %     nsampsFrameRx = nsampsFrame, nsampsFrameTx = nsampsFrame, sampleRate = sampleRate, centerFrequency = fc);
   % 위의 함수 않되는데?
   devtx=sdrtx('Pluto','SamplesPerFrame',nsampsFrame,'BasebandSampleRate',sampleRate,'CenterFrequency',fc);
   devrx=sdrrx('Pluto','SamplesPerFrame',nsampsFrame,'BasebandSampleRate',sampleRate,'CenterFrequency',fc);
   %장치여러개일경우 'RadioId'로 지정가능
end

%복소수 신호 만들기
% Set the normalized digital frequency and gain
%X= exp(2pift)신호만들건데 디지털 정규화주파수에서 프레임마다 반복되게만들려면 nu=정수/nsampsFrame
nu0 = 4/nsampsFrame; %nu= fc/fsamp?= 78.125
A = 1;

% TODO:  Create a digital signal the length of one frame with a digital frequency
% of nu0
%    x = ...
x=A*exp(1i*2*pi*nu0*linspace(0,nsampsFrame-1,nsampsFrame)); %1프레임마다 4번반복하는 일종의신호
% fc fsamp 일단생각x




% TODO:  Plot the real and imaginary components of x.
% Use subplot to plot the two components on different plots
% side-by-side.  Also plot the x vs. the time in micro-seconds
clf;
figure(1)
time=linspace(0,nsampsFrame-1,nsampsFrame)/sampleRate;
subplot(2,1,1)
plot(time.*1e+06,real(x))
xlabel('microSecond')
legend('real')
subplot(2,1,2)
plot(time.*1e+6,imag(x))
legend('imag')
xlabel('microSecond')

%%%%%%%%%%전송하기%%

if runTx && ~usePreRecorded
    devtx.release();
    devtx.transmitRepeat(x.');  %디지털신호를 보내면 알아서 passband로변환해주는구나
    
    %devtx.transmitRepeat(complex(zeros(4096,1)));
end
% If not running the RX, stop the live script now
if ~runRx
    return;  %% 스크립트 즉시종료.
end

% If using pre-recorded samples, load them now
% This will return:  
% r = receive samples, 
% centerFreq = center frequency in Hz
if usePreRecorded   %%% 미리저장된샘플 사용할때,
    load freqSamples;
end

fprintf("The normalized digital frequency (nu0): %.2f \n",nu0);
f0=nu0*sampleRate;
fprintf("The baseband frequency f0 = %.2f \n",f0);  % f0/fsamp = nu0 30kHz
fprintf("The passband frequency f0pb= %.2f \n",fc+f0)  % 알아서 전송해주는구나..

%%수신하기%%%
nbits = 12;               % number of ADC bits, 2^12=4096
fullScale = 2^(nbits-1);  % 12비트 adc는 -2047~2048의값을가짐
                          % full scale  그냥 정규화값으로 나누는걸 matlab으로함



if ~usePreRecorded

    % Capture data
    r = devrx.capture(nsampsFrame); % 지정된샘플만 수집하나봐
    
    % TODO:  Scale to floating point
    %   r = ...
    r = single(r)/fullScale;  % 2048로 나누어 정규화
end


% TODO:  Plot the real and imaginary RX samples vs. time in micro-seconds
figure(2);
subplot(2,1,1)
plot(time.*1e+06,real(r))
xlabel('microSecond')
legend('SDR real')
%ylim([-1 1])

subplot(2,1,2)
plot(time.*1e+6,imag(r))
legend('SDR imag')
xlabel('microSecond')
%ylim([-1 1])

if ~usePreRecorded
   save freqSamples r  %신호 r을 freqSamples로 저장하라
end


% Correlation 방법을 이용한 캐리어 주파수 추정하기..
% 가정: r[n]=  A*exp(2*pi*1i*nu*n) + noise 일때,
% z = \sum_{n=1 ~ r} r[n+1]*conj( r[n] ) -> 크기제곱*exp(f
%  nu = angle(z)/2/pi
% gpt물어보면 잘알려줌
% 어그냥함수실행할거야
nuestTx = estFreq(x.', method='correlation'); %정규화주파수 nu를 구하는거
nuestRx = estFreq(r, method='correlation');

fprintf("신호x 추정주파수 : %d Hz \n",nuestTx*sampleRate); %% ㄷㄷ 잘되노..
fprintf("수신 신호r 추정주파수 : %d Hz \n",nuestRx*sampleRate);

%CFO 와 PPM계산하기
% cfo:  The carrer frequency offset in Hz (the difference between the TX and RX frequencies)
% ppm:  The error in PPM defined as:  ppm = cfo / fc * 1e6 .  The devices are rated for about ~10 ppm error.  So, the total error should be |ppm| < 20.  
cfo=nuestTx*sampleRate-nuestRx*sampleRate
ppm=cfo/fc*1e6

%이제 원신호x의 크기 추정하기, Least-square method, 수신신호r을통해
% 
% r[n]=A*exp(2 pi nu n) +noise 로 이미 nu를 추정함
% noise를 무시하고 r[n]과 A*exp를 빼면 크기가 0이겠지? 그래서 제일작은값이 나오는 A를찾는다
% 식을 미분해서 0이 되는 값을찾으면 최소 A를 찾을수있다. (GPT ㄱㄱ)
% 기준 신호 u[n] 만들기
u = exp(2*pi*1i * nuestRx * linspace(0,nsampsFrame-1,nsampsFrame)).';           % u[n] = e^(j*2pi*nu*n)
% 최소제곱으로 c 추정
c_hat =  sum(r .* conj(u)) / sum(abs(u).^2);   % Least-squares 공식 / 잘않되노 뭐지..?
fprintf("추정되는 진폭 : %.2f",abs(c_hat))

%% 진폭추정 잘않되네...

