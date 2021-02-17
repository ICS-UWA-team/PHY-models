function y = replayfilter_bb(x, fs_x, h, fs_tau)
%     persistent mem
%     if isempty(mem)
%         mem = 0;
%     end

%     % bring signal to baseband
%     t = (0:length(x)-1)'/fs_x;
%     x=x.*exp(-2*pi*fc*1i*t);

%     h=zeros(1, 100);
%     h(78) = 1+1j;

    % resample signal to sampling frequency of channel estimate
    [N, D]=rat(fs_tau / fs_x);
    x=resample(x, N, D);

    y=0*x; % initiate output signal
    h=flipud(h.'); % rotate and flip
    [K, dummy]=size(h);
    x=[zeros(K-1,1); x]; % zero padding
    Ly=length(y);

%     Direct-replay channel simulation
    for k=0:Ly-1
        f = k/K - floor(k/K); 
        n = floor(k/K) + 1; % impulse response counter
        ir = (1 - f)*h(:, n) + f*h(:, n+1); % linear interpolation
        y(k+1) = x(k+1:k+K).'*ir; % filtering
    end
    
% %         Direct-replay channel simulation
%     for k=0:Ly-1
%         f = (k + mem)/K - floor((k + mem)/K); 
%         n = floor((k + mem)/K) + 1; % impulse response counter
%         ir = (1 - f)*h(:, n) + f*h(:, n+1); % linear interpolation
% %         H=fft(ir);
% %         H1 = H(1023:2044,1);
% %         H2 = H(1:1022,1);
% %         H = [H1; H2];
% %         ir=ifft(H);
%         y(k+1) = x(k+1:k+K).'*ir; % filtering
%         mem = mem + 1;
%         if n == dummy - 1
%             mem = 0;
%         end
%     end

%     H=fft(h(:,1));
% %     PlotDatShit(H, length(H) - 1, fs_x);
%     H1 = H(1023:2044,1);
%     H2 = H(1:1022,1);
%     H = [H1; H2];
%     h=ifft(H);
% 
%     for k=0:Ly-1
%         y(k + 1) = x(k+1:k+K).' * h(:,100); % filtering
%     end

%     y = x;
    
    % resample waveform to sampling frequency of input signal
    [N,D]=rat(fs_x/fs_tau);
    y=resample(y,N,D);

%     % shift to passband
%     t=(0:length(y)-1)'/fs_x;
%     y=real(y.*exp(2*pi*fc*1i*t));
return
