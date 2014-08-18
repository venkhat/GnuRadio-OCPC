function predictphi=kalamnv(y,var_noise,l)

%% Initial assumptiona  for kalman independent of iteraions
global delf;
var_pn=4*pi*delf;
%var_noise=1;
alphasq=10^-3;
beta=2;
L=l; %number of samples of unscnted is 2L+1
lambda=L*(alphasq-1);
%weights of unscented klamna filter
wm(1)=lambda/(L+lambda);
wc(1)=wm(1)+(1-alphasq+beta);
wm(2:2*L+1)=1/(2*(L+lambda));
wc(2:2*L+1)=wm(2:2*L+1);
global xcap;
global statevar;
%%
    predictx=xcap;
    predictvar=statevar+var_pn;
    %unscented transform
    u(1)=predictx;
    m=1:L;
    u(m+1)=predictx+((L+lambda)*predictvar)^0.5;
    u(L+m+1)=predictx-((L+lambda)*predictvar)^0.5;
    %non linear observation at unscented values
    v=exp(j*u);
    %observstion y
    predicty=sum(wm.*v);
    var_yy=sum(wc.*(v-predicty).*conj(v-predicty));
    var_xy=sum(wc.*(u-predictx).*conj(v-predicty));
    K=var_xy/(var_yy+var_noise);
    xcap=predictx+K*(y-predicty);
    predictphi=real(xcap);%output
    statevar=predictvar-K*conj(var_xy);
end
