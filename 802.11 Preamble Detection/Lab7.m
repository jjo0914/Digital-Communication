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

x = wlantx();  % This calls the stepImpl() method.

% TODO:  Get the sample rate from wlantx.fsamp
%    fsamp = ...


% TODO:  Plot the absolute values of the samples in x vs. time in us.  If you did everything
% correctly, you should see that the packet is about 250 us long. 
