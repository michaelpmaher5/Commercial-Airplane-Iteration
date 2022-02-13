function K = getFuselageFormFactor(LoD)

    K = -0.0019*LoD^3 + 0.0477*LoD^2 - 0.4157*LoD + 2.4156;
end