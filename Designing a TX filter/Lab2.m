clc;
clear all;
clf;
%Unit2 Lab1
% QAM심볼만들기 -> 업샘플링 ->디지털필터(펄스 Shaping) -> 전송 -> PSD측정

% Wifi 시스템 Parameters
fsamp = 16e6;   % Signal sample rate in Hz (before upconversion)
fchan = 20e6;   % Channel bandwidth


%QAM심볼 만들기
nsym = 2^14;    % Number of symbols 
Rmod = 4;       % Modulation rate 4비트=1심볼 ?
M = 2^Rmod;     % QAM order M=16? 16심볼?

% TODO
% nbits = ...
% bits = ...
nbits=nsym*Rmod;
bits=randi([0 1],nbits,1);

% You should set 'InputType' to 'bit'. Also, set 'UnitAveragePower' to true so the symbols have unit average power.
sym=qammod(bits,M,'InputType','bit','UnitAveragePower',true);
%플롯팅
% scatterplot(sym); 또는 cd = comm.ConstellationDiagram; cd(y);   

% 이제 업샘플링과 Pulse shaping을 할건데 , 만들어진 함수를 수정해서 쓰라고 합니다.

% Set parameters
Astop = 40;    % dB 감쇠/Stopband rejection in dB 
ovRatio = 2;   % 2배 업샘플링/Oversampling ratio

% TODO:  Set remaining parameters
%   Fp = passband frequency
%   Fst = stopband frequency
%   fsampUp = sampling freq after upsampling
Fp= fsamp/2; % 필터는 절반값만 입력받아서 그런가 ? Fp should be set to >= sampling rate / 2 before upsampling
Fst=  fchan/2; %필터는 절반만 그려줘서 ? Fst should be set to the channel bandwidth / 2
fsampUp=ovRatio * fsamp;

% Construct the filter class
txFilt = TxFilt("ovRatio", 2, "rateIn", fsamp, "Fp", Fp, "Fst", Fst, "Astop", Astop);

% Run the design filter method 
txFilt.designFilter();

% Get the filter taps
bfilt = txFilt.bfilt; % 그 필터 분모계수?라는데 아 다른필터는 계수 지정을해줘야하는데!
npts = 512;

% TODO
%   [H,w] = freqz(bfilt, 1, npts);
%   plot(...);
[H,w]=freqz(bfilt,1,npts);
figure(1);
plot(w/2/pi,abs(H));

% Create the upsampled and filtered symbols
xup = txFilt(sym); % 필터링을 안하면 그냥심볼은 샘플링주파수 내에 전부존재하는 사각펄스가되버린다->구현못함x
                   % sym:시간영역에서 거의 랜덤신호 or 직사각형신호 -> 주파수영역에서 평탄한스펙트럼
                   % -> 
figure(2)
plot(real(sym))
hold on;
plot(real(xup))
legend('Original Symbol','Upsampling&Pulse Shaping')

%%PSD그리기 함수딸깍ㅋㅋㅋ
[Px,fx] = pwelch(xup,hamming(512),[],[],fsampUp,'centered');
figure(3)
plot(fx/10e6,pow2db(Px))
hold on;

xline(Fp/10e6,'Color','b');
hold on;

xline(Fst/10e6,'color','r')

xlabel('MHz')

% % 플루토로 신호를 보내기전에 전력크기를조정해야한다.
% 즉 신호의 크기는 -1 1 사이(full-scale일때,)로 전력은 -9dB낮게.
% TODO:  Complete the code in `applyBackoff` section of the TxFilt.stepImpl().

% 함수 수정했고 backoff수치에따른 PSD그리기
backoffLevTest = [0,6,12];
figure(4);
ntest = length(backoffLevTest);
for i = 1:ntest

    % Enable auto-scaling for backoff and set the backoff level
    backoffLev = backoffLevTest(i);
    txFilt.set('backoffAutoScale', true, 'backoffLev', backoffLev);

    % Create the upsampled and filtered symbols 
    xup = txFilt(sym);

    % TODO:  Compute the PSD of xup
    %   [Px,fx] = pwelch(...)
    [Px,fx] = pwelch(xup,hamming(512),[],[],fsampUp,'centered');
    % Create a sub-plot
    subplot(1,ntest,i);

    % TODO:  Plot the PSD.  Make sure all sub-plots have the same ylimit,
    % so you can compare the levels.
    %    plot(...)
    %    ylim([...]) 
    %    title(...)  % Set title for the backoff
    plot(fx/10e6,pow2db(Px));
    ylim([-140 -65])
    title('backoff %d PSD값 비교',backoffLev);
    % 0일때다른 이유는 이제 그 최대 최소값 컷 했기때문인듯..?
    
end
% o ㅅ ㅇ ㅋㅋㅋㅋㅋㅋ
%플루토 설정방법
% TODO:  Set parameters
% Set to true for loopback, else set to false
loopback = true;  

% Select to run TX and RX
runTx = true;
runRx = true;

% Center frequency
centerFrequency = 2.4e9;  

% Flag indicating if pre-recorded samples are to be used.  Use this if you
% have no Pluto devices.
usePreRecorded = false;
saveData = true;

% clear previous instances
clear devtx devrx

% add path to the common directory where the function is found
addpath('C:\Users\juwon\Desktop\Digital Communication\Unit02\Lab');

% Parameters
nsampsFrame = length(xup);

% Run the creation function.  Skip this if we are using pre-recorded
% samples
if ~usePreRecorded
    [devtx, devrx] = plutoCreateTxRx(createTx = runTx, createRx = runRx, loopback = loopback, ...
        nsampsFrameTx = nsampsFrame, nsampsFrameRx=nsampsFrame, sampleRate = fsampUp, ...
        centerFrequency=centerFrequency);
  
end
% 플루토 발싸~ Now TX the clipped data xup repreatedly through the Pluto device.
if runTx && ~usePreRecorded
    % TODO:  Use the devtx.release and devtx.transmitRelease commands to
    % continuously send up 
    devtx.release(); % 한번 초기화>?
    devtx.transmitRepeat(xup);
else% If not running the RX, stop the live script now if ~runRx
    return;
end
nbits=12;
fullScale=2^(nbits-1) % ADC 가 2047~2048 값을 가지는데 -1 ~1로 맞추기 위해서?

if ~usePreRecorded
    r=devrx.capture(nsampsFrame); %값이 -2048~2047로 저장되는데...
    r=single(r)/fullScale; % 만약single안하면?.. 필터계산을 할때 float32형태가 필요하기 때문
    if saveData
        save symModSing r; %r을 symModSing로저장
    end
    % Load pre-recorded data
    load symModSing;
end
%--Measuring PSD
% TODO:  Create an initla plot as before
%    [Prx,fx] = pwelch(...);
%    plot(...
[Prx,fx]=pwelch(r,hamming(512),[],[],fsampUp,'centered');
figure(5);
plot(fx/10e6,pow2db(Prx));
title('SDR 수신값 PSD backoff는 12(마지막 전송값)');

%Create spectrum analyzer
%반복문 사용해서 실시간 캡처하기
if runTx && ~usePreRecorded
    % TODO:  Use the devtx.release and devtx.transmitRelease commands to
    % continuously send xup
    devtx.release();     
    devtx.transmitRepeat(xup);
end
figure(6);
% TODO:  Compute the initial PSD in dB and frequency in MHz
%    [Prx,fx] = pwelch(...);
%    fxMHz = fx/1e6;
%    PrxdB = ...;


% Clear plot and plot data initial time
fxMHz=fx/10e6
Prxdb=pow2db(Prx)
p = plot(fxMHz, Prxdb, 'LineWidth', 2); %플롯을 변수로생성할수잇구나

% TODO:  Fix the limits to some reasonable values
%   xlim([...]);
%   ylim([...]);
%   xlabel(...);
%   ylabel(...);
xlim([-2 2]);
ylim([-130 -70]);
xlabel('MHz');
ylabel('dB');
% Set pointers to X and Y data
p.XDataSource = 'fxMHz'; %변수명입력
p.YDataSource = 'Prxdb';

% Number of loops before the spectrum analyzer stops
ncaptures = 100; %캡처할 프레임수

% Load pre-recorded data if required
if usePreRecorded
    load symModMult;
    ncaptures = size(rdat,2);
else
    rdat = zeros(nsampsFrame, ncaptures);
end


for t = 1:ncaptures
    if ~usePreRecorded

        % TODO:  Capture data
        %    r = devrx.capture(...);
         r = devrx.capture(nsampsFrame  );
          r=single(r)/fullScale;
        % Save data if required
        if saveData
            rdat(:,t) = r;
        end
    else
        % Load pre-recorded
        r = rdat(:,t);        
    end


    % TODO:  Re-compute the PSD
    %    [Prx,fx] = pwelch(...);
    %    fxMHz = ...
    %    PrxdB = ...    
    [Prx,fx] = pwelch(r,hamming(512),[],[],fsampUp,'centered');
    fxMHz=fx/10e6;
    Prxdb=pow2db(Prx);

    
    % Redraw plot
    refreshdata(p); %p가 바뀌면
    drawnow           %새로그린다.
    
end

% Redraw plot

% Save data
if saveData
    save symModMult rdat;
end
