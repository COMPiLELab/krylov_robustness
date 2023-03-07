function betha = get_step_size(beta_init, rho, x_curr, a_curr, f, grad_f)
    c = 1e-4;
    betha = beta_init;
    fk = f(x_curr, a_curr);
    [gk, Gk] = grad_f(x_curr,a_curr);
    
    while f(x_curr - betha * Gk, a_curr - betha * gk) > fk + c * betha * (norm(Gk,'fro')^2 + gk^2)
        betha = rho * betha;
    end
end
