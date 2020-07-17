function g=Bluestein_czt(x,k,w,a)
[m,n]=size(x);
nfft=2^nextpow2(m+k-1);
kk=((-m+1):max(k-1,m-1)).';
ww=w.^((kk.^2)./2);
nn=(0:(m-1))';
aa=a.^(-nn).*ww(m+nn);
y=x.*aa(:,ones(1,n));
fy=fft(y,nfft);
fv=fft(1./ww(1:(k-1+m)),nfft);
fy=fy.*fv(:,ones(1,n));
g=ifft(fy);
g=g(m:(m+k-1),:).*ww(m:(m+k-1),ones(1,n));
end