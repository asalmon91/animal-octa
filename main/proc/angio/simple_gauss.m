function sg = simple_gauss(x, m, bw)
%simple_gauss is a simple Gaussian, see Jia 2012 Opt Express Eq.2

sig = bw*max(x)/(2*sqrt(2*log(2)));
sg = exp(-((x-m).^2 / (2*sig^2)));

end

