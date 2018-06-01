clc
clear all
%PID-PID
global E c_p c_pc dH k0 Mc U V rho t_m t_c w Caf T_w T_f T_m w_c A
 
E=1422; %K
c_p=3.8; c_pc=4.2; %kJ/(kg*K)
dH=-2017; %kJ/kg
k0=0.09; %m^3/(kg*min)
Mc=2000; %kg
U=1471; %kJ/(min*m^3)
V=250; %m^3
rho=1281; %kg/m^3
t_c=1; t_m=0.2; %min
w=454; w_c=476; %kg/min
Caf=144; %kg/m^3
T_w=308; T_f=339; T_m=339; %K
A=2000; %m^2

C=6500; % C=Cv*sqrt(dPv/Gf)


X0=[144;339;320;339];
tspan=0:0.02:200;

 Kc_TC=-0.0001; Ti_TC=0.001; Td_TC=(1/4)*Ti_TC;
Kc_FC=0.0000001; Ti_FC=0.3; Td_FC=(1/4)*Ti_FC;
wc_sp=0; w_c=476; r=339; l=0.5; 
t_st=[]; x_st=[]; r_st=[]; w_st=[]; wc_st=[]; wc_sp_st=[];

for i=1:length(tspan)-1
      if i==round(length(tspan)*0.5)
        r=350;
      elseif i==round(length(tspan)*0.7)
              w=500;
      end  
    %Temperature controller
    if i==1
        epp=0; ep=0; e=0;  
    elseif i==2
        epp=0; ep=e; e=r-X(end,4); 
    else
        epp=ep; ep=e; e=r-X(end,4);
    end  
    
    dwc_sp=Kc_TC*(e-ep + (1/Ti_TC)*e + Td_TC*(e-2*ep+epp));
    wc_sp=wc_sp+dwc_sp;
    
    if wc_sp<0
        wc_sp=0;
    end


   %Flow Controller
   if i==1
       epp_FC=0; ep_FC=0; e_FC=0;  
    elseif i==2
       epp_FC=0; ep_FC=e_FC; e_FC=wc_sp-w_c;
    else
       epp_FC=ep_FC; ep_FC=e_FC; e_FC=wc_sp-w_c;
    end  
    
   dl=Kc_FC*(e_FC-ep_FC + (1/Ti_FC)*e_FC + Td_FC*(e_FC-2*ep_FC+epp_FC));
    l=l+dl;
    
    if l>1
        l=1;
    elseif l<0
       l=0;
    end
   
   w_c=C*l;
 
%     process
    if i==1
        x0_a=X0;
    else
        x0_a=X(end,:)';
    end
    
    tspan_a=[tspan(i) tspan(i+1)];
    [t,X]=ode15s(@func_cstr,tspan_a,x0_a);
    
    %data save
    t_st=[t_st;t]; x_st=[x_st;X]; r_st=[r_st;r*ones(length(t),1)]; w_st=[w_st;w*ones(length(t),1)];
    wc_st=[wc_st;w_c*ones(length(t),1)]; wc_sp_st=[wc_sp_st;wc_sp*ones(length(t),1)];
end

e_sum=sum(abs(r-x_st(:,4)))/length(x_st(:,4));
figure(1),plot(t_st,x_st(:,4),t_st,r_st,t_st,w_st);  % hold on;

xlabel('time [s]')
ylabel('Temperature of material in the reactor [K]')
legend('CV','SP','DV');

figure(2),plot(t_st,wc_st,t_st,wc_sp_st); legend('CV-FC','SP-FC');
