function data=decoderv(y)
%%% for 4QPSk symbols
if real(y)>0
    data= 1+0j;
else
    data=-1+0j;
    
end
if imag(y)>0
    data=data+1j;
else
    data=data-1j;
end
end

