function dCdp = figure6(CloverClmax,segment)

    if segment == 1
        dCdp = 0.1107*CloverClmax^3 - 0.0437*CloverClmax^2 - 0.0483*CloverClmax + 0.0321;
    else
        dCdp = 0.0821*CloverClmax^3 + 0.003*CloverClmax^2 - 0.066*CloverClmax + 0.0408;
    end

end