%actual algorithm given in paper directly implemented.
clc;

clear all;
fs=1;%symbol frequency
delf=.01*fs;
var_pn=4*pi*(delf/fs);
var_noise=1;
alphasq=10^-3;
beta=2;
lambda=alphasq-1;
wm(1)=lambda/(1+lambda);
wc(1)=wm(1)+(1-alphasq+beta);
wc(2)=1/(2+2*lambda);
wc(3)=wc(2);
wm(2)=wc(2);
wm(3)=wc(2);

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
    u(2)=predictx(n)+((1+lambda)*predictvar(n))^0.5;
    %u(3)=u(2)
    %trial
    u(3)=predictx(n)-((1+lambda)*predictvar(n))^0.5;
    %non linear observation at unscented values
    v=exp(j*u);
    %observstion y
    predicty(n)=sum(wm.*v);
    var_yy=sum(wc.*(v-predicty(n)).*conj(v-predicty(n)));
    var_xy=sum(wc.*(u-predictx(n)).*conj(v-predicty(n)));
    K(n)=var_xy/(var_yy+var_noise);
    xcap(n)=predictx(n)+K(n)*((y(n)-predicty(n)));%trial
    %xcap(n)=predictx(n)+K(n)*(y(n)-predicty(n)); %actual
    %xcap(n)=real(xcap(n)); %trial
    statevar(n)=predictvar(n)-K(n)*conj(var_xy);
    
end

figure,plot(1:iter,x,'bo-');

hold on;
plot(1:iter,xcap,'r+-');

legend('original phase','kalman estimnated phase');
title('with real value of xcap')
q=(x-xcap).*conj(x-xcap);
%figure,plot(1:iter,q);
meanerror=sum(q)/numel(q)
