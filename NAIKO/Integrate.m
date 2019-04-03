function [integral] = Integrate(Forces)
    integral = trapz(Forces);
end

