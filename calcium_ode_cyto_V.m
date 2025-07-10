% Equation 6 in Li&Rinzel paper 
% including cyto and ER volumes

function dydt = calcium_ode_cyto_V(~, y, p)
    c_cyto = y(1);
    h = y(2);
    u_er = y(3);
    

    % IP3 receptor term
    m_inf = (p.IP3/(p.IP3 + p.d_ip3)) * (c_cyto/(c_cyto + p.d_act));
    %J_chan = p.c1* (p.v1 * m_inf^3 * h^3 + p.v2);
    J_chan = (p.v1 * m_inf^3 * h^3 + p.v2);

    % SERCA pump (cyto to ER)
    J_pump = p.v3 * c_cyto^2 / (p.k_er^2 + c_cyto^2);

    % Clear pump (cyto to extracelullar)
    J_clear = p.v4 * c_cyto^2 / (p.k_pl^2 + c_cyto^2);

    % SOCE flux
    J_soce = p.k_soce * (p.c_cs - u_er);

    % h dynamics
    h_inf = p.d_inh / (c_cyto + p.d_inh);
    tau_h = 1 / (p.a * (c_cyto + p.d_inh));
    dhdt = (h_inf - h) / tau_h;

    % du_er/dt
    dcerdt = (-J_chan * (u_er - c_cyto) ...
             + J_pump ...
             + J_soce)/p.VER;

        % d(c_er)/dt = (1 + kD*S/(u_er + kD)^2) * du_er/dt
    gamma = 1 + p.k_D * p.S / (u_er + p.k_D)^2;

    % du_er/dt (unbound calcium change in ER)
    duerdt = dcerdt / gamma;

    % dc_cyto/dt
    dc_cytodt = (J_chan * (u_er - c_cyto) ...
            - J_pump - J_clear)/p.Vcyto; ...

    dydt = [dc_cytodt; dhdt; duerdt; dcerdt];
end
