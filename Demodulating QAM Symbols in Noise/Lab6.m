clc; clf; clear all;
% Unit 6 LAB..
% 1단계.16 QAM변복조하고 BER계산하기 ㅇㅇ
modRate = 4;    % num bits per symbol (4=16 QAM)
M = 2^modRate;  % QAM order
nsym = 1e5;     % num symbols

% TODO
% nbits = ...
% bits = ...
% sym = ...
nbits=modRate*nsym;
bits=randi([0 1],nsym*modRate,1);
sym=qammod(bits,M,'InputType','bit','UnitAveragePower',true); % -> 평균 POWER =1

% NOISE 추가하기 Eb/N0=5 dB. ->ES로 바꿔야하는거아님?
% Remmeber to scale the noise variance correctly!
% SNR이 총 두종류거든  ber&ser에 쓰이는 snr은 Eb(비트1개의 평균에너지)/N0 -> (Watt)/(W/HZ) ->심볼1개당 SNR
%                     -> Pa/Pb
% Eb(비트1개에너지) 랑 P(평균전력) 랑 구하는식은 같은데 의미가 다르다
EbN0=5; % dB
EsN0=pow2db(db2pow(EbN0)*modRate); %   W -> dB
symbol_noise=db2mag(-EsN0); %            아데시벨에 마이너스 해줘야하네
                                         % dB -> W -> 비트 수 곱
                                         % bit per symbol * E_b= E_s
w=(1/sqrt(2)).*(randn(nsym,1)+1i.*randn(nsym,1));
r=sym+symbol_noise.*w;   % 심볼노이즈가 0.2812 인데 작으면 작을수록 경계가 뚜렷해지네
scatterplot(sym);
scatterplot(w);
scatterplot(r);

% 디모듈레이트 하기.
% TODO
% bithat = qamdemod(...)
% ber = ...
bithat=qamdemod(r,M,'OutputType','bit','UnitAveragePower',true);
ber=mean(bits~=bithat)
bertheory=0.75*qfunc(sqrt(0.2.*db2pow(EsN0))) % 오잉 이론값이 조금 다른데
                                              % 예제에서 준거보다 gpt가 준 이론값이
                                              % 더정확해,,,!
%%---------- 위 코드를 이제 다양한 EbN0에서 심볼을 디모듈레이션하고 
%%  BER 구하고 이론값과 비교하는건데 Exercise 4 에서 이미 함
%% 이제 Phase error 를 시뮬레이션 할껀데.
%% 심볼 모듈레이션 하면 별자리가 생기지? 이걸 임의의각도 theta로 그대로 회전시킬거임

% Phase errors to test
thetaStdTest = [0,0.01,0.02,0.04]*2*pi; % 테스트할 심볼 회전각도
ntest = length(thetaStdTest);

% SNRs to test 여기다가 SNR까지추가 싯팔..
EbN0Test = (-3:20)';  
nsnr = length(EbN0Test);

% Measure the BER on the phase corrupted symbols
ber1 = zeros(nsnr,ntest); % SNR + 위상각도 때마다 생기는 BER오류측정..
legStr = cell(ntest,1);
for itheta = 1:ntest

    % Get phase error
    thetaStd = thetaStdTest(itheta);

    for it = 1:nsnr
        % Get EbN0 to test
        EbN0 = EbN0Test(it);

        % TODO:  Create the phase corrupted symbols
        %    theta = ...
        %    sym1 = ...        
        theta=thetaStd*randn(1,1); % 분산조절해서 가우시안분포에서 한개뽑음
        sym1=sym.*exp(1i*theta);
       
        % TODO:  Add random noise
        %    w = ...
        %    r = sym1 + w;  ...
        EsN0=pow2db(db2pow(EbN0)*modRate); %   W -> dB
        symbol_noise=db2mag(-EsN0); %            아데시벨에 마이너스 해줘야하네
                                         % dB -> W -> 비트 수 곱
                                         % bit per symbol * E_b= E_s
        w=(1/sqrt(2)).*(randn(nsym,1)+1i.*randn(nsym,1));
        r=sym1+symbol_noise.*w;
        % TODO:  Demodulate bits and compute BER
        %     bithat = qamdemod(...);
        %     ber1(it,itheta) = ...
        bithat=qamdemod(r,M,'OutputType','bit','UnitAveragePower',true);
        ber1(it,itheta)=mean(bits~=bithat);
                 
    end
     scatterplot(sym1);
    % Display progress
    fprintf(1,'phase error =%5.3f*2pi completed\n', thetaStd/2/pi);    
    legStr{itheta} = sprintf('phase error = %5.3f *2pi', thetaStd/2/pi);
end

% TODO:  Plot the results on one plot with one curve per phase error
% Add the legend and adjust the axes

plot(EbN0Test,ber1(:,1));
hold on
plot(EbN0Test,ber1(:,2));
hold on
plot(EbN0Test,ber1(:,3));
hold on
plot(EbN0Test,ber1(:,4));
hold on
legend(legStr(1:4));
xlabel('Eb/N0')
ylabel('BER')



