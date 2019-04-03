%% Isp Model Thrust
m_prop = 0.12; %kg
m_i = 0.22; %kg
m_f = 0.1; %kg
g = 9.8; %m/s^2
Forces = ReadIn("LA8am_test1");
Integral = Integrate(Forces);

Isp = Integral./(m_prop*g);
del_V = Isp*g*ln(m_i/m_f);

%% Interpolate Data
