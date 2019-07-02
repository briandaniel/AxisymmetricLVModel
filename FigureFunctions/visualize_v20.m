%%
close all

if strcmp( param.odeSystem, 'full') 
    Ppa = Y(:,1);
    Ppp = Y(:,2);
    Psp = Y(:,5);
    Pra = Y(:,7);
    Ppv = Y(:,8);
    Psa = Y(:,10);
    Vla = Y(:,11);
    Pao = Y(:,12);
    a1 = Y(:,13);
    a2 = Y(:,14);
    a3 = Y(:,15);

    Plv = otherVara_vec(:,1);
    Rmv = otherVara_vec(:,2);
    Raov = otherVara_vec(:,3);
    Prv = otherVara_vec(:,4);
    Vlv = otherVara_vec(:,5);
    Vrv = otherVara_vec(:,6);
    Pla = otherVara_vec(:,7);
    phi = otherVara_vec(:,8);
    phi_a = otherVara_vec(:,9);

    Tshow = .01;
    dt = param.dt0;

    ic = floor(Tshow/dt);

    figure('outerposition', [700,1000,1300,400])
    plot(T(ic:end), Plv(ic:end), T(ic:end), Prv(ic:end), T(ic:end), ...
            Pao(ic:end), T(ic:end), Psa(ic:end), T(ic:end), Psp(ic:end))    
    legend('P_l_v', 'P_r_v','P_a_o','P_s_a','P_s_p')

    figure('outerposition', [700,600,1300,400])
    plot(T(ic:end), Pla(ic:end), T(ic:end), Pra(ic:end), T(ic:end), Plv(ic:end),...
        T(ic:end), Prv(ic:end), T(ic:end), Ppv(ic:end), T(ic:end), Pra(ic:end),...
        T(ic:end), Psp(ic:end))
    legend('P_l_a', 'P_r_a','P_l_v','P_r_v','P_p_v', 'P_r_a', 'P_s_p')
    % axis([Tshow, max(T), -.2, 2 ])

    figure('outerposition', [700,280,1300,400])
    plot(T(ic:end),Vlv(ic:end),T(ic:end),Vrv(ic:end),T(ic:end),Vla(ic:end))
    legend('V_l_v', 'V_r_v', 'V_l_a')

    figure('outerposition', [700,100,1300,200])
    plot(T(ic:end),Rmv(ic:end),T(ic:end),Raov(ic:end))
    legend('R_m_v', 'R_a_o_v')

    figure('outerposition', [200,300,600,800])
    plot(Vlv(ic:end), Plv(ic:end), Vrv(ic:end), Prv(ic:end))

    figure('outerposition', [700,100,1300,200])
    plot(T(ic:end),phi(ic:end),T(ic:end),phi_a(ic:end))
    legend('phi', 'phi_a')

    figure('outerposition', [1400,100,600,600])
    plot(Vla(ic:end), Pla(ic:end))
    legend('P_l_a', 'V_l_a')

elseif strcmp( odeSystem, 'circOnly')
      
    Ppa = Y(:,1);
    Ppp = Y(:,2);
    Psp = Y(:,5);
    Pra = Y(:,7);
    Ppv = Y(:,8);
    Psa = Y(:,10);
    Pla = Y(:,11);
    Pao = Y(:,12);
    a1 = Y(:,13);
    a2 = Y(:,14);
    a3 = Y(:,15);
    
    Plv = otherVara_vec(:,1);
    Rmv = otherVara_vec(:,2);
    Raov = otherVara_vec(:,3);
    Prv = otherVara_vec(:,4);
    Vlv = otherVara_vec(:,5);
    Vrv = otherVara_vec(:,6);
    Pla = otherVara_vec(:,7);
    phi = otherVara_vec(:,8);
    phi_a = otherVara_vec(:,9);

    
    Tshow = .01;
    dt = param.dt0;
    ic = floor(Tshow/dt);

    figure('outerposition', [700,1000,1300,400])
    plot(T(ic:end), Plv(ic:end), T(ic:end), Prv(ic:end), T(ic:end), ...
            Pao(ic:end), T(ic:end), Psa(ic:end), T(ic:end), Psp(ic:end))    
    legend('P_l_v', 'P_r_v','P_a_o','P_s_a','P_s_p')

    figure('outerposition', [700,600,1300,400])
    plot(T(ic:end), Pla(ic:end), T(ic:end), Pra(ic:end), T(ic:end), Plv(ic:end),...
        T(ic:end), Prv(ic:end), T(ic:end), Ppv(ic:end), T(ic:end), Pra(ic:end),...
        T(ic:end), Psp(ic:end))
    legend('P_l_a', 'P_r_a','P_l_v','P_r_v','P_p_v', 'P_r_a', 'P_s_p')
    % axis([Tshow, max(T), -.2, 2 ])

    figure('outerposition', [700,280,1300,400])
    plot(T(ic:end),Vlv(ic:end),T(ic:end),Vrv(ic:end))
    legend('V_l_v', 'V_r_v')

    figure('outerposition', [700,100,1300,200])
    plot(T(ic:end),Rmv(ic:end),T(ic:end),Raov(ic:end))
    legend('R_m_v', 'R_a_o_v')

    figure('outerposition', [200,300,600,800])
    plot(Vlv(ic:end), Plv(ic:end), Vrv(ic:end), Prv(ic:end))

    figure('outerposition', [700,100,1300,200])
    plot(T(ic:end),phi(ic:end),T(ic:end),phi_a(ic:end))
    legend('phi', 'phi_a')
     
elseif strcmp( odeSystem, 'lvOnly')
    
    Pla = Y(:,1);
    Pao = Y(:,2);
    a1 = Y(:,3);
    a2 = Y(:,4);
    a3 = Y(:,5);

    Plv = otherVara_vec(:,1);
    Rmv = otherVara_vec(:,2);
    Raov = otherVara_vec(:,3);
    Vlv = otherVara_vec(:,4);
    At = otherVara_vec(:,5);
    
    Tshow = .01;
    dt = param.dt0;
    ic = floor(Tshow/dt);

    figure('outerposition', [700,1000,1300,400])
    plot(T(ic:end), Plv(ic:end), T(ic:end), Pao(ic:end), T(ic:end), Pla(ic:end))  
    legend('P_l_v', 'P_a_o', 'P_l_a')
    
    figure('outerposition', [700,600,1300,400])
    plot(T(ic:end), a1(ic:end), T(ic:end), a2(ic:end), T(ic:end), a3(ic:end))    
    legend('a_1', 'a_2','a_3')
    
    figure('outerposition', [700,280,1300,400])
    plot(T(ic:end),Vlv(ic:end))
    legend('V_l_v')

    figure('outerposition', [700,100,1300,200])
    plot(T(ic:end),At(ic:end))
    legend('A(t)')
    
end































