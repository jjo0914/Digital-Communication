
% 레일리 페이딩 모델은 신호가 직접적으로 수신기에 도달하지않고
% 여러 확산경로를 통해 도달하는 상황을 모델랑한것 (ex: 도시) (신호가 반사 회절 간섭 된 경로로 도달함,,)
% 여기서 진폭이 레일리분포에서 추출된값들이다..
% (왜? : 실수 허수 부분이 평균0인 가우시안분포를따른다고 하면 그 합들의 크기는 레일리분포라서,,,)



%% Rayleigh Fading Theory
% Rayleigh fading occurs when a receiver sees the same signal from multiple
% paths that arrive at different delays and angles of arrivals.  
% When the receiver moves in such a multi-path environment,
% it will experience a time-varying complex channel gain of the form:
%
%    h(t) = sqrt(H0/npath)* \sum_k exp(i*(phi0(k) + 2*pi*fdHz(k)*t) <이식은
%    그냥이런게있다..?
% 
%    위상차이 : 도달하는 시간이 다르니깐 주파수 모양이다름=위상차이
%    도플러 주파수: 여기서 수신기가 움직인깐 위상차이가 시간에대한 함수로나타니고
%    위상을 시간이대헤 미분해서 각주파수로 만들고 그걸 주파수로 표현한게 도플러 주파수

% where |npath| is the number of paths, and |phi0(k)| and |fdHz(k)| 
% are the phase (in radians) and Doppler shift (in Hz) 
% of the |k|-th path and |H0| is the path power gain in linear scale.  
% The Doppler shifts are given by:
%
%    fdHz(k) = fdmaxHz*cos(theta(k))  
%    왜 코사인이냐? 하면, 서로 방향이 달라서
%    도플러 효과를 파의 이동방향과 수신기의 이동방향과의 내적
%    fdmaxHz = v/vc*fc
% 
% where |theta(k)| is the angle of arrival of the |k|-th path relative
% to the direction of motion, |v| is the receiver velocity in m/s;
% |vc=3e8| is the speed of light and |fc| is the carrier frequency in Hz.
% The fluctuation of the gain due to constructive and destructive
% interference between paths is one of the central challenges in wireless
% communication.

%% Simulation parameters
% In this lab, we will simulate what is called a single *path cluster*, 
% where a large number of paths arrive with slightly different angles of
% arrival (AoAs) but similar time delays.  Typically, a path cluster occurs
% whenever there is a diffuse reflection causing scattering of the paths.  
% We will simulate the path cluster with following parameters:

% 단일 path cluster란  많은 다중경로신호들이 거의 같은 시간 지연으로 도착하지만,
% 도달방햑(입사각) 이 다른 경우 (ex: 건물 벽안 에서 신호가 확산 반사가 일어난다.)

% AoA : 벽에 튕기고난 각도
clear all;
npath = 100;           % Number of sub-paths in the cluster
thetaCenDeg = 0;       % Center AoA in degrees
thetaSpreadDeg = 10;   % AoA spread in degrees
pathGaindB = -10;      % Path gain in dB
vkmph = 30;            % RX velocity (km/h)
fcGHz = 28.0;          % Carrier freq in GHz


%% Creating a PathCluster object
% When performing complex simulations, it is useful to write the code in an
% object-oriented manner.  For this lab, a skeleton file, 
% |PathCluster_partial.m| has been included.  Copy the file and rename it as
% |PathCluster.m|.  The |properties| and constructor are already completed.
% You will fill in the remaining parts.  Note that if you edit the 
% |PathCluster| object, you will need to clear the environment to 
% ensure that MATLAB reloads the new version:
%
%     clear PathCluster
%
% For this section, you do not need to modify |PathCluster.m|.  
% You just need to create an instance of the class with a command of the
% form:
%
%     cluster = PathCluster('prop1', val1, 'prop2', val2, ...)

% 복잡한 시뮬레이션을 위해 객체지향 코드방식을 사용한 PatchCluster.m을제공할테니
% 너가 수정해라. 지금 단계예서는 class instance를 불러와라

% TODO: 
%
% * Create a |cluster| object and set the parameters of the cluster
% * Print the properties using |disp(cluster)| to ensure the properties are
% set
cluster=PathCluster('pathGaindB',pathGaindB,'npath',npath,'thetaCenDeg',thetaCenDeg,'thetaSpreadDeg',thetaSpreadDeg,'vkmph',vkmph,'fcGHz',fcGHz);
disp(cluster);

%% Compute path cluster parameters
% Now, modify the |computeFd| method in the |PathCluster| class to generate
% random AoA and Doppler shifts for the paths.
% 
% * Set |theta(k) =| Gaussian with mean thetaCenDeg and 
%   std. dev thetaSpreadDeg.  Convert to radians
% * Set |phi0(k) =| uniform between 0 and 2*pi
% * set |fdHz(k) =| Doppler shift with the formula in the Theory section.
% % Set |H0 =| path gain in linear scale from |pathGaindB|

% TODO:
% 
% * Modify the |computeFd()| method
% * Run |cluster.computeFd()| 
% * Plot the locations of Doppler spreads 
%  경로별 도플러 주파수값들을 시각화하라
cluster.computeFd();
[x,idx]=sort(cluster.theta);
y=cluster.fdHz(idx);
plot(x/pi*180,y);
ylabel('도플러 주파수');
xlabel('Angle of arrivals(degree)')


%% Generate a random fading trajectory
% Next, complete the code in 
% TODO:
%
% * Complete the code in the method |genFading()|.
% * Run the code for 1024 time points over 0.1 sec.
% * Plot the channel gain in dB.
% h(t)를 만들어라
h=cluster.genFading(0.1); % 0.1 초동안 샘플만들기
figure(2);
plot(linspace(0,0.1,1024).',mag2db(abs(h)))
xlabel('시간(s)')
ylabel('Channel Gain(dB)')


%% Measure the autocorrelation
% Now, we will measure the autocorrelation.  Accurate estimation of the
% autocorrelation requires large numbers of samples.  
%
% TODO:  Generate a random trajectory of |h(t)| with 
% |nt=2^16| samples at a sampling period of |tsamp = 0.1| ms.
% 똑같이하면되는데 자기상관은 더많은 샘플이필요하다 R=E[x(t)*x(t + tau)]를  못구해서
% 샘플이많아야 기대값에 가까워진다..

% nt = 2^16;
% tsamp = 0.1;
% h = ...
nt=2^16;
tsamp=0.1;
t=linspace(0,tsamp,nt).';
h=zeros(length(t),1);
for k=1:length(t)
    h(k)=sqrt(cluster.H0/cluster.npath) * sum(exp(i*(cluster.phi0 + 2*pi*cluster.fdHz*t(k))));
end


% Use the |xcorr| function to estimate the autocorrelation at |nlags=1000|
% lags.  Use the |'unbiased'| scaling option.  
%
% TODO
% nlags = 1000;
% [rh, lags] = xcorr(...);
nlags=1000; % +- tsamp/nt * 1000초
[rh,lags]=xcorr(h,nlags,'unbiased');
% -nlags ~ lags 지연값계산..
% TODO:  Plot the magnitude of the estimated autocorrelation as a 
% function of the time lag.
figure(3)
plot(lags*(t(2)-t(1)),rh);
xlabel('time lags');
ylabel('자기상관값 (클 수록 신호 변화없음)')

%% 
% The time it takes the process to be uncorrelated is called the
% *coherence* time of the channel.  Within the coherence time, the channel
% does not vary significantly.  This time is important for various
% procedures in wireless communication such as channel estimation.

%신호가 자기 자신과 상관성이 사라질 때까지 걸리는 시간을
%**채널의 코히런스 타임 (coherence time)이라고 부른다.
%코히런스 타임 이내에서는 채널이 크게 변하지 않는다.
%이 시간(코히런스 타임)은 채널 추정과 같은
% 무선 통신의 다양한 기술 절차에 있어 매우 중요한 기준이 된다.

%% Compute the theoretical autocorrelation
% Finally, we will compare the measured auto-correlation with the
% theoretical value.  When the number of paths is large, the
% auto-correlation should approximately be:
%
%     R(tau) = H0*\int exp(2*pi*1i*fdmax*tau*cos(theta))*p(theta)dtheta,
%
% where |p(theta)| is the PDF on the AoA's theta.  This integral has no 
% closed form expression.  So, it needs to be evaluated numerically.  You
% can do this in MATLAB with the function |integral| where you supply it a
% function pointer to the integrand.
%

% 레일리 페이딩의 자기상관값 이론값이 R(tau)일때 이론값을 구해라
% p(theta) = 입사각의 확률밀도 함수 ..를 어떻게구하지?.. -> (0,10/180 * pi)가우시안 분포
%   ->> normpdf(x범위, 평균,sigma)  로 구할수있다

% TODO:
fdmax=cluster.fdmax
p=@(theta) normpdf(theta,thetaCenDeg,sqrt(thetaSpreadDeg/180*pi));

t=linspace(-0.1,0.1,200);
for tau=1:length(t)
    R(tau)=cluster.H0*integral(@(theta) exp(2*pi*1i*fdmax*t(tau)*cos(theta)).*p(theta),-pi,pi);
end

hold on
plot(t,abs(R.'));

legend('자기상관계수값','자기상관계수 이론값');
% 
% * Using numerical integration compute the theoretical autocorrelation
% |R(tau)| above for 200 values of |tau| in [-100,100] ms.
% * Plot the theoretical autocorrelation vs. tau on the same graph as the
% estimated autocorrelation.  There may be a significant difference since
% the estimate autocorrelation was computed with only 100 paths.


