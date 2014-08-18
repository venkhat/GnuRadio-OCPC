%  Symbol error calculation for qpsk with comparison of Pe with phase n
%  without phase.
clc;
close all; clear all;
iter=10000;
snrcount=0;
global xcap;
global statevar;
global delf;
delf=0.001;
var_pn=4*pi*delf;
Pe=zeros(1,10);
Pe1=zeros(1,10);
x=zeros(1,iter);
Predictxx=zeros(1,iter);
for snrdb=10:20
    snr=10^(snrdb/10);
    snrcount=snrcount+1;
    %%Transmitter section 
    a=randsrc(1,1)+1j*randsrc(1,1);
    h=sqrt(1/2)*randn(1,1)+1j*randn(1,1);
    sigpow=abs(a*h)^2;
    var_noise=sigpow/snr;
    w=randn(1,1)+1j*randn(1,1); %complex noise
    w=sqrt(var_noise/2)*w;
    
    pn=sqrt(var_pn)*randn(1,1);
    x(1)=0+pn;%initial phase=0;
    
    rx=a*h*exp(1j*x(1))+w;
    %%receiver
    rx1=rx*conj(h)/abs(h);
    data=decoderv(rx1);
    error=0;
    if data==a
        error=error+0;
    else error=error+1;
    end
    %% without phase part
    rx=a*h+w;
    rx1=rx*conj(h)/abs(h);
    data=decoderv(rx1);
    error1=0;
    if data==a
        error1=error1+0;
    else error1=error1+1;
    end
    %%
    xcap=0;
    statevar=var_pn;
    predictxx(1)=xcap;
    for n=2:iter
        
        a=randsrc(1,1)+1j*randsrc(1,1);
        h=sqrt(1/2)*randn(1,1)+1j*randn(1,1);
        sigpow=abs(a*h)^2;
        var_noise=sigpow/snr;
        w=randn(1,1)+1j*randn(1,1); %complex noise
        w=sqrt(var_noise/2)*w;
        pn=sqrt(var_pn)*randn(1,1);
        x(n)=x(n-1)+pn;
        rx=a*h*exp(1j*x(n))+w;
        rx1=rx*conj(h)*exp(-1j*predictxx(n-1))/abs(h);
        data=decoderv(rx1);
        if data==a
            error=error+0;
        else error=error+1;
        end
        y=rx/(data*h);
        var_noisek=var_noise/(abs(a*h)^2);%noise variance when it goes to kalman
        predictxx(n)=kalmanv(y,var_noisek,1);%output
        %% without phase
        rx=a*h+w;
        rx1=rx*conj(h)/abs(h);    
        data=decoderv(rx1);
        if data==a
            error1=error1+0;
        else error1=error1+1;
        end
        
    end
%     figure,plot(1:iter,x,'bo-');
%     hold on;
%     plot(1:iter,predictxx,'r+-');
    Pe(snrcount)=error/iter;
    Pe1(snrcount)=error1/iter;
end
% figure,plot(1:iter,x,'bo-');
% hold on;
% plot(1:iter,predictxx,'r+-');

figure, plot(10:20,Pe,'ro--');
hold on;
plot(10:20,Pe1,'b+-');
legend('With phase', 'without phase');