function dX=func_cstr(t,X)
        

global E c_p c_pc dH k0 Mc U V rho t_m t_c wc_sp w Caf T_w T_f w_c A



dCadt=(w*(Caf-X(1,1))-V*rho*k0*exp(-E/X(2,1))*X(1,1)^2)/(V*rho);
dX(1,1)=dCadt;
dTdt=(w*c_p*(T_f-X(2,1))-U*A*(X(2,1)-X(3,1))+(-dH)*V*k0*exp(-E/X(2,1))*X(1,1)^2)/(V*rho*c_p);
dX(2,1)=dTdt;
dT_cdt=(U*A*(X(2,1)-X(3,1))+w_c*c_pc*(T_w-X(3,1)))/(Mc*c_pc);
dX(3,1)=dT_cdt;
dT_mdt=(X(2,1)-X(4,1))/t_m;
dX(4,1)=dT_mdt;




