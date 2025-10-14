function [hfd,h,snr] = estChanResp(r,xfd,opt)
% estChanResp:  Estimates the channel response using frequency-domain
% correlation
% 
% The model is that the TX repeatedly sends x=ifft(xfd) and the RX 
% receives samples
%     r = h*x + w
% where h is the channel impulse response and w is noise.  
%
% Parameters
% ----------
% r:  RX complex baseband samples of length nfft
% xfd:  TX samples in frequency-domain
% nleft, nright:  Number of samples to the left and right of peak
%     used for signal energy computation
% normToNoise:  Boolean flag indicating that the energy / tap
%     should be normalized to the noise energy estimate.
%     In this case abs(h(k))^2 represents the SNR per tap.
%
% Returns
% -------
% hfd:   Frequency domain channel estimate.  hfd(k) is the estimate
%    of the channel frequency response at f(k) = k/nfft*fsamp
% h:  Time-domain channel impulse response estimate.
% snr:  Total snr in dB.
arguments
    r (:,1) double;
    xfd (:,1) double;
    
    opt.nleft (1,1) {mustBeInteger} = 8; 
    opt.nright (1,1) {mustBeInteger} = 8;
    opt.normToNoise (1,1) = false;

end

% Create empty outputs until they are set.
% You can delete this code when you have set the variables
hfd = [];
h = [];
snr = [];

% TODO:  Take the FFT of r 
%    rfd = ...
rfd=fft(r,length(xfd)); % rfd=hfd xfd + W   

% TODO:  Estimate the channel frequency response by dividing by the 
% FFT of x, xfd
%    hfd = ...
hfd=rfd./xfd;


% TODO:  Compute the time-domain response
%    h = ...
h=ifft(hfd,length(xfd));  % h탭들이 결정되고 가장강한값을 첫번째인덱스에 맞추면(주경로)
                           %
           
% TODO:  Find the peak location in h
% Then, circularly shift h so that the peak is at position opt.nleft
%    h = circshift(...): 배열을 지정한 숫자만큼 왼쪽으로 이동 

[~,idx]=max(abs(h));  % 그냥 max(h) 가아니라 max(abs(h))ㄷㄷ
h=circshift(h,opt.nleft+1-idx);  % 어쨋든 피크를 경로에맞춤


% TODO:  Compute the SNR
%   Enoise = ...
%   Esig = ...
%   snr = ... (in dB)
 
% 여기서 Rfd=Hfd*Xfd+Wfd  (원래는)
% Hfd=Rfd/XFD (채널추정은 이렇게) 
% Hfd'(추정)=Hfd(진짜) + Wfd/Xfd
% (ifft) h'(추정)=h(진짜)+ w' ((h(진짜)=h'일정길이))


% 여기서 h=h'(일정길이)    % nleft+1 ~ nright까지 구간이 신호라 가정
                           % 왜? 물리적으로 채널임펄스응답이  어느 구간 이상을 넘을수없어
                          % 어찌보면 탭수를 제한해서 진짜 h를 가정하는거지
    % (2) RFD/XFD=Hfd+Wfd/XFD
    %     r'=h+w'
  
nsig=opt.nright+opt.nleft+1;
Enoise=mean(abs(h(opt.nright+opt.nleft+1:end)).^2); %  =N0 단위는 W/HZ (N0는 일정 ㅇ)  h' - h =w'    (h=h'(일정길이)) 
Etot=sum(abs(h(1:opt.nright+opt.nleft+1)).^2)+nsig*Enoise; % r'=h+w' 애는 Watt*sample

Esig=Etot-nsig*Enoise; % 수정 . 원래:Etot-Enoise [J] -[W] 샘플곱하면 W->J
%Esig=mean(abs(h(opt.nleft+1:opt.nright+opt.nleft+1)).^2)
%snr=pow2db(Esig/Enoise) % 단위 안맞는이유? 이산신호는 원래 ?(Watt*[sample]/ Watt)=[sample] 다이렇다고?)
                        % sample이라는게 원래 무차원 단위래
                        % [총신호에너지/평균잡음전력]= SNR !! or [신호전력/N0]=SNR
snr  = pow2db( Esig / Enoise );
% TODO:  Normalize to the noise level
%    h = ...
if opt.normToNoise  % 그냥 플롯팅용이래
  h=h./sqrt(Enoise);
end

end