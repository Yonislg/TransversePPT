% Function calculates secondary electron emissions
function ge = SEE(gi, E_i, W)
    if E_i-W>W
        gamma = 0.016 *(E_i - 2* W);
        ge = gi .* gamma;
    else
        ge = 0;
    end
end