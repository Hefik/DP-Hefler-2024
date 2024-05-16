clear all; close all; clc
%% Parametry modelu
t=1e-2;
F =[1 t
    0 1]; %% dynamics matrix
H=[1 0];
c=299792458;
q1=0.04/c^2;
q2=0.01/c^2;

 q1=0.02/c^2;
 q2=0.003/c^2;

Qd=[q1*t+q2*t^3/3 q2*t*t/2
    q2*t*t/2 q2*t]; % discrete noise covariance
Qdd=chol(Qd)';
R=2.1e-19;
R=1e-19;
Rd=sqrt(R);
Theta_ex=[q1;q2;R];
%% Matice MDM
L=5;
N=1;
P=L+N;
S_index=1+N;
E_index=S_index+L;
r=P+2*(P-1); % Velikost Psi
Psi=zeros(r*r,4);

Q11N=(1:2*r+2:r*2*(P-1))';
Q12N=[(2:2*r+2:r*2*(P-1))'; (r+1:2*r+2:r*2*(P-1))'];
Q22N=(r+2:2*r+2:r*2*(P-1))';
RN=(r*2*(P-1)+r-P+1:r+1:r*r)';

Psi(Q11N,1)=1;
Psi(Q12N,2)=1;
Psi(Q22N,3)=1;
Psi(RN,4)=1;


% Upsilon=@(T)[T T^3/3 0
%          0 T^2/2 0
%          0 T 0
%          0 0 1];
%% Monte Carlo
par=3;
cykly=10000;
n=1e5;%1,1574 dne ---- %43200000; % 5 dní ---- T=0,01s
rng(1)
lts=100;
ts=linspace(t,3,lts);

Theta_MC=nan(par, cykly,lts);
Theta_MCap=nan(par, cykly,lts);
Theta_MCk=nan(par, cykly,lts);
mW=[0;0];
tic
%%

parfor k=1:cykly
    % System
    ts=linspace(t,3,lts);

for p=1:lts
    tt=round(ts(p)/t);
    nt=length(1:tt:n);
    F =[1 ts(p)
    0 1];
    Qd=[q1*ts(p)+q2*ts(p)^3/3 q2*ts(p)*ts(p)/2
    q2*ts(p)*ts(p)/2 q2*ts(p)]; % discrete noise covariance
    Qdd=chol(Qd)';
    X=zeros(2,nt);
    W=Qdd*randn(2,nt);
    V=Rd*randn(1,nt);
    for i = 1:nt
        X(:,i+1)= F*X(:,i)+W(:,i);
    end
    Z=H*X(:,2:end)+V;
    
    [A,Aw,Av]=MatA(ts(p),L,N);
    Ups=[ts(p) ts(p)^3/3 0
         0 ts(p)^2/2 0
         0 ts(p) 0
         0 0 1];
    Zk=zeros(L+N,nt-L);
    for j=N+1:nt-L+1
        Zk(:,j)=Z(:,j-N:j-N+P-1)';
    end
    Zd=Av*Zk(:,N+1:end);
    Z_tilda=Zd*Zd';
    Zdk=Z_tilda(:)./length(Zd(1,:));
    Theta_MC(:,k,p)=pinv(kron(A,A)*Psi*Ups)*Zdk;
    Theta_MCap(:,k,p)=pinv(kron(A,A)*Psi(:,[1,3,4]))*Zdk;
    Theta_MCk(:,k,p)=pinv(kron(A,A)*Psi(:,[1,3,4]))*(kron(A,A)*Psi*Ups)*[q1;q2;R];
    
end

    k
    

end

toc
%% Vykresleni - novy odhad 
%close all
%Ts=[t,t*5,t*20,t*50,t*75,t*100,t*150,t*200,t*250,t*300,t*350,t*400,t*450,t*500,t*600,t*700];
Ts=ts';
Means1=squeeze(mean(Theta_MC(1,:,:)));
Means2=squeeze(mean(Theta_MC(2,:,:)));
Means3=squeeze(mean(Theta_MC(3,:,:)));
% Var1=1./sqrt(Ts).*[var(Theta_MC(1,:)), var(Theta_MC5(1,:)),var(Theta_MC20(1,:)),var(Theta_MC50(1,:)),var(Theta_MC75(1,:)),var(Theta_MC100(1,:)),var(Theta_MC150(1,:)),var(Theta_MC200(1,:)),var(Theta_MC250(1,:)),var(Theta_MC300(1,:)),var(Theta_MC350(1,:)),var(Theta_MC400(1,:)),var(Theta_MC450(1,:)),var(Theta_MC500(1,:)),var(Theta_MC600(1,:)),var(Theta_MC700(1,:))];
% Var2=1./sqrt(Ts).*[var(Theta_MC(2,:)), var(Theta_MC5(2,:)),var(Theta_MC20(2,:)),var(Theta_MC50(2,:)),var(Theta_MC75(2,:)),var(Theta_MC100(2,:)),var(Theta_MC150(2,:)),var(Theta_MC200(2,:)),var(Theta_MC250(2,:)),var(Theta_MC300(2,:)),var(Theta_MC350(2,:)),var(Theta_MC400(2,:)),var(Theta_MC450(2,:)),var(Theta_MC500(2,:)),var(Theta_MC600(2,:)),var(Theta_MC700(2,:))];
% Var3=[var(Theta_MC(3,:)), var(Theta_MC5(3,:)),var(Theta_MC20(3,:)),var(Theta_MC50(3,:)),var(Theta_MC75(3,:)),var(Theta_MC100(3,:)),var(Theta_MC150(3,:)),var(Theta_MC200(3,:)),var(Theta_MC250(3,:)),var(Theta_MC300(3,:)),var(Theta_MC350(3,:)),var(Theta_MC400(3,:)),var(Theta_MC450(3,:)),var(Theta_MC500(3,:)),var(Theta_MC600(3,:)),var(Theta_MC700(3,:))];
Var1=nan(lts,1);
Var2=nan(lts,1);
Var3=nan(lts,1);
for o=1:lts
    Var1(o)=var(Theta_MC(1,:,o));
    Var2(o)=var(Theta_MC(2,:,o));
    Var3(o)=var(Theta_MC(3,:,o));
end
figure
plot(Ts,Means1)
hold on
grid on
yline(q1,'k-.')
plot(Ts,Means1+1*sqrt(Var1),'g--')
plot(Ts,Means1-1*sqrt(Var1),'g--')
%ylim([0.035,0.048])
xlim([t,3])
title('Vliv velikosti T_s na kvalitu odhadu parametru q_1')
xlabel('T_s')
ylabel('q_1')
legend('Střední hodnota odhadu', 'Parametr q_1', '1-\sigma okolí')

figure
plot(Ts,Means2)
hold on
grid on
yline(q2,'k-.')
plot(Ts,Means2+1*sqrt(Var2),'g--')
plot(Ts,Means2-1*sqrt(Var2),'g--')
%ylim([0.005,0.015])
xlim([t,3])
title('Vliv velikosti T_s na kvalitu odhadu parametru q_2')
xlabel('T_s')
ylabel('q_2')
legend('Střední hodnota odhadu', 'Parametr q_2', '1-\sigma okolí')

figure
plot(Ts,Means3)
hold on
grid on
yline(R,'k-.')
plot(Ts,Means3+1*sqrt(Var3),'g--')
plot(Ts,Means3-1*sqrt(Var3),'g--')
%ylim([-0.003,0.004])
xlim([t,3])
title('Vliv velikosti T_s na kvalitu odhadu parametru R')
xlabel('T_s')
ylabel('R')
legend('Střední hodnota odhadu', 'Parametr R', '1-\sigma okolí')

%% Vykresleni odhad s aproximaci
Means1ap=1./Ts.*squeeze(mean(Theta_MCap(1,:,:)));
Means2ap=1./Ts.*squeeze(mean(Theta_MCap(2,:,:)));
Means3ap=squeeze(mean(Theta_MCap(3,:,:)));

Means1k=1./Ts.*squeeze(mean(Theta_MCk(1,:,:)));
Means2k=1./Ts.*squeeze(mean(Theta_MCk(2,:,:)));
Means3k=squeeze(mean(Theta_MCk(3,:,:)));

% Var1=1./sqrt(Ts).*[var(Theta_MC(1,:)), var(Theta_MC5(1,:)),var(Theta_MC20(1,:)),var(Theta_MC50(1,:)),var(Theta_MC75(1,:)),var(Theta_MC100(1,:)),var(Theta_MC150(1,:)),var(Theta_MC200(1,:)),var(Theta_MC250(1,:)),var(Theta_MC300(1,:)),var(Theta_MC350(1,:)),var(Theta_MC400(1,:)),var(Theta_MC450(1,:)),var(Theta_MC500(1,:)),var(Theta_MC600(1,:)),var(Theta_MC700(1,:))];
% Var2=1./sqrt(Ts).*[var(Theta_MC(2,:)), var(Theta_MC5(2,:)),var(Theta_MC20(2,:)),var(Theta_MC50(2,:)),var(Theta_MC75(2,:)),var(Theta_MC100(2,:)),var(Theta_MC150(2,:)),var(Theta_MC200(2,:)),var(Theta_MC250(2,:)),var(Theta_MC300(2,:)),var(Theta_MC350(2,:)),var(Theta_MC400(2,:)),var(Theta_MC450(2,:)),var(Theta_MC500(2,:)),var(Theta_MC600(2,:)),var(Theta_MC700(2,:))];
% Var3=[var(Theta_MC(3,:)), var(Theta_MC5(3,:)),var(Theta_MC20(3,:)),var(Theta_MC50(3,:)),var(Theta_MC75(3,:)),var(Theta_MC100(3,:)),var(Theta_MC150(3,:)),var(Theta_MC200(3,:)),var(Theta_MC250(3,:)),var(Theta_MC300(3,:)),var(Theta_MC350(3,:)),var(Theta_MC400(3,:)),var(Theta_MC450(3,:)),var(Theta_MC500(3,:)),var(Theta_MC600(3,:)),var(Theta_MC700(3,:))];
Var1ap=nan(lts,1);
Var2ap=nan(lts,1);
Var3ap=nan(lts,1);
for o=1:lts
    Var1ap(o)=1./sqrt(Ts(o)).*var(Theta_MCap(1,:,o));
    Var2ap(o)=1./sqrt(Ts(o)).*var(Theta_MCap(2,:,o));
    Var3ap(o)=var(Theta_MCap(3,:,o));
end
figure
plot(Ts,Means1ap)
hold on
grid on
yline(q1,'k-.')
%plot(Ts,Means1k,'g--');
plot(Ts,Means1ap+1*sqrt(Var1ap),'m--')
plot(Ts,Means1ap-1*sqrt(Var1ap),'m--')
%ylim([-0.145,0.16])
xlim([t,3])
title('Vliv velikosti T_s na kvalitu odhadu parametru q_1')
xlabel('T_s')
ylabel('q_1')
legend('Střední hodnota odhadu', 'Parametr q_1', '1-\sigma okolí')

figure
plot(Ts,Means2ap)
hold on
grid on
yline(q2,'k-.')
%plot(Ts,Means2k,'g.');
plot(Ts,Means2ap+1*sqrt(Var2ap),'m--')
plot(Ts,Means2ap-1*sqrt(Var2ap),'m--')

%ylim([-0.018,0.026])
xlim([t,3])
title('Vliv velikosti T_s na kvalitu odhadu parametru q_2')
xlabel('T_s')
ylabel('q_2')
legend('Střední hodnota odhadu', 'Parametr q_2', '1-\sigma okolí')

figure
plot(Ts,Means3ap)
hold on
grid on
%plot(Ts,Means3k,'g--');
yline(R,'k-.')
plot(Ts,Means3ap+1*sqrt(Var3ap),'m--')
plot(Ts,Means3ap-1*sqrt(Var3ap),'m--')

%ylim([-0.055,0.06])
xlim([t,3])
title('Vliv velikosti T_s na kvalitu odhadu parametru R')
xlabel('T_s')
ylabel('R')
legend('Střední hodnota odhadu', 'Parametr R', '1-\sigma okolí')

%% Porovnani

figure
plot(Ts,Means1)
hold on
grid on
plot(Ts,Means1ap)
yline(q1,'k-.')
plot(Ts,Means1+1*sqrt(Var1),'g--')

plot(Ts,Means1ap+1*sqrt(Var1ap),'m--')
plot(Ts,Means1-1*sqrt(Var1),'g--')
plot(Ts,Means1ap-1*sqrt(Var1ap),'m--')
%ylim([0.035,0.048])
xlim([t,3])
title('Porovnání nového odhadu a odhadu s aproximací parametru q_1')
xlabel('T_s')
ylabel('q_1')
legend('Střední hodnota odhadu','Střední hodnota aprox. odhadu' , 'Parametr q_1', '1-\sigma okolí odhadu', '1-\sigma okolí aprox. odhadu')

figure
plot(Ts,Means2)
hold on
grid on
plot(Ts,Means2ap,'-.')
yline(q2,'k-.')
plot(Ts,Means2+1*sqrt(Var2),'g--')

plot(Ts,Means2ap+1*sqrt(Var2ap),'m--')
plot(Ts,Means2-1*sqrt(Var2),'g--')
plot(Ts,Means2ap-1*sqrt(Var2ap),'m--')
%ylim([0.005,0.015])
xlim([t,3])
title('Porovnání nového odhadu a odhadu s aproximací parametru q_2')
xlabel('T_s')
ylabel('q_2')
legend('Střední hodnota odhadu','Střední hodnota aprox. odhadu' , 'Parametr q_2', '1-\sigma okolí odhadu', '1-\sigma okolí aprox. odhadu')

figure
plot(Ts,Means3)
hold on
grid on
plot(Ts,Means3ap,'--')
yline(R,'k-.')
plot(Ts,Means3+1*sqrt(Var3),'g')

plot(Ts,Means3ap+1*sqrt(Var3ap),'m--')
plot(Ts,Means3-1*sqrt(Var3),'g')
plot(Ts,Means3ap-1*sqrt(Var3ap),'m--')
%ylim([-0.002,0.0022])
xlim([t,3])
title('Porovnání nového odhadu a odhadu s aproximací parametru R')
xlabel('T_s')
ylabel('R')
legend('Střední hodnota odhadu','Střední hodnota aprox. odhadu' , 'Parametr R', '1-\sigma okolí odhadu', '1-\sigma okolí aprox. odhadu')


%%
%save('Simulace_nova2')
