

function val = weightFunction_Riemann(s, ti, ti1, alpha, beta,proxy,velocityShift)
    eps_val = 1e-10;
    val = zeros(size(s));

    for k = 1:length(s)
        sk = s(k);

        % Beta verschoben um "velocityShift"
        arg1 = (1 / max(sk - ti1, eps_val) - velocityShift(1))/(velocityShift(2)-velocityShift(1));
        arg1 = min(max(arg1, 0), 1);

        arg2 = (1 / max(sk - ti, eps_val) - velocityShift(1))/(velocityShift(2)-velocityShift(1));
        arg2 = min(max(arg2, 0), 1);
        
        if proxy
     % Dirac -> Heaviside
            F1 = double(arg1 >= alpha/(alpha+beta));
            F2 = double(arg2 >= alpha/(alpha+beta));
        else
            % Beta Verteilung
            F1 = betacdf(arg1, alpha, beta);
            F2 = betacdf(arg2, alpha, beta);

            F1 = unifcdf(max(0,(1 / max(sk - ti1, eps_val))),velocityShift(1),velocityShift(2));
            F2 = unifcdf(max(0,(1 / max(sk - ti, eps_val))),velocityShift(1),velocityShift(2));
        end

        
        val(k) = max(F1 - F2, 0);
    end
end