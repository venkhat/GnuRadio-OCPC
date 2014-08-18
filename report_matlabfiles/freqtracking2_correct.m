%with unscneted kalman with genreal L samples and as given in oscillator phase
%noise paper 
clc;
close all;
clear all;
fs=1;
delf=.01*fs;
var_pn=4*pi*(delf/fs);
var_noise=1;
alphasq=10^-3;
beta=2;
L=1; %number of samples of unscnted is 2L+1
lambda=L*(alphasq-1);
%weights of unscented klamna filter
wm(1)=lambda/(L+lambda);
wc(1)=wm(1)+(1-alphasq+beta);
wm(2:2*L+1)=1/(2*(L+lambda));
wc(2:2*L+1)=wm(2:2*L+1);
%let phi be x
iter=100;
pn=sqrt(var_pn)*randn(1,iter);
w=randn(1,iter)+1j*randn(1,iter); %complex noise
w=sqrt(var_noise/2)*w;
x(1)=0+pn(1);%initial phase=0;
y(1)=exp(1j*x(1))+w(1);

xcap(1)=0;
statevar(1)=var_pn;

for n=2:iter
    x(n)=x(n-1)+pn(n);
    y(n)=exp(1j*x(n))+w(n);
    
    predictx(n)=xcap(n-1);
    predictvar(n)=statevar(n-1)+var_pn;
    %unscented transform
    u(1)=predictx(n);
    m=1:L;
    u(m+1)=predictx(n)+((L+lambda)*predictvar(n))^0.5;
    u(L+m+1)=predictx(n)-((L+lambda)*predictvar(n))^0.5;
    %non linear observation at unscented values
    v=exp(j*u);
    %observstion y
    predicty(n)=sum(wm.*v);
    var_yy=sum(wc.*(v-predicty(n)).*conj(v-predicty(n)));
    var_xy=sum(wc.*(u-predictx(n)).*conj(v-predicty(n)));
    K(n)=var_xy/(var_yy+var_noise);
    xcap(n)=predictx(n)+K(n)*(y(n)-predicty(n));%trial
    %xcap(n)=predictx(n)+K(n)*(y(n)-predicty(n)); %actual
   % xcap(n)=real(xcap(n)); %trial
    statevar(n)=predictvar(n)-K(n)*conj(var_xy);
    
end

figure,plot(1:iter,x,'bo--');

hold on;
plot(1:iter,xcap,'r+-');

legend('original phase','kalman estimnated phase');

%figure,plot(1:iter,(x-xcap).*conj(x-xcap));

q=real(x-xcap).*real(x-xcap);
meanerror=sum(q)/numel(q)
c1=1-(sqrt(q)./real(x));
c2=mean(c1)