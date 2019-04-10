function [psdx, freq] = lfp_tfa_powerspectrum(x, Fs, foi)
    N = length(x);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(x):Fs/2;
    psdx = psdx(freq >= foi(1) & freq <= foi(end));
end