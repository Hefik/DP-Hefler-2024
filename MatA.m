function [A,Aw,Av]= MatA(T,L,N)
%% Parametry modelu

F =[1 T
    0 1]; %% dinamics matrix
H=[1 0];

O_L=zeros(L,2);
for i=0:L-1
    O_L(i+1,:)=H*F^i;
end
Ov=O_L;

F_L=F^N;
Xi=F^(N-1);
for i=N-2:-1:0
    Xi=[Xi,F^i];
end

% pro N=0
if N==0
    Xi=zeros(2,N*2);
end

Gamma=zeros(L,2*(L-1));
k=0;
for i=1:2:2*(L-1)+1
    for j=(i+3)/2:L
        Gamma(j,i:i+1)=H*F^(j-2-k);
    end
    k=k+1;
end

A1=[Ov*Xi, Gamma];
A2= [-Ov*F_L*pinv(Ov)*Gamma,zeros(L,N*2)];

A3=[zeros(L,N), eye(L)];
A4=[-Ov*F_L*pinv(Ov),zeros(L,N)];

Aw=[eye(L) eye(L)]*[A1;A2];
Av=[eye(L) eye(L)]*[A3;A4];
A=[Aw Av];

end

