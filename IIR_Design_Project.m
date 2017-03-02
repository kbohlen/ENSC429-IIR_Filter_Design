%%%
% IIR Impulse Invariance Design Project
% Summer 2016
%
% Author(s): Kurtis Bohlen
%%%

%%%
% DESIGN SPECS for DISCRETE TIME
% 
% 1 - delta_l <= |H_d(e^(jw))| <= 1 + delta_l ,          |w| <= 0.25pi 
% |H_d(e^(jw))| <= delta_h                    , 0.4pi <= |w| <= pi
%%%

%%%
% DISCRETE TIME FILTER SAMPLING
%
% In the impulse invariance design procedure for transforming a continuous-
% time filter h_c(t) into a discrete-time filter, the discrete-time 
% filter is chosen:
%     h[n] = T * h_c(n*T)
% For symplicity we will use a sampling period (T) of 1
%%%

%%%
% DESIGN SPECS for CONTINUOUS TIME (Qusetion 1)
% 0.7071 <= |H(e^(jW))| <= 1.2929 ,            |w| <= 0.25pi/T 
% |H(e^(jW))| <= 0.3162           , 0.4pi/T <= |w| <= infinity
%%%

%%%
% CONTINUOUS TIME BUTTERWORTH FILTER DESIGN (Question 2)
%                          1
% |H(e^(jW))|^2 = --------------------
%                 1 + ( W / W_c )^(2N)
%%%
delta_l = [0.2929 0.1 0.01 0.001];
delta_h = [sqrt(0.1) sqrt(0.1) sqrt(0.01) sqrt(0.001)];

for a = 1:numel(delta_l)

H_pass = 1 - delta_l(a)
H_stop = delta_h(a)

H_pass2 = H_pass^2
H_stop2 = H_stop^2

H_pass2Inv = 1/H_pass2
H_stop2Inv = 1/H_stop2

% 
% |H(e^(jW))|^-2 = 1 + ( W / W_c )^(2N)
%

H_passWFrac = H_pass2Inv - 1
H_stopWFrac = H_stop2Inv - 1

H_stop_pass = H_stopWFrac / H_passWFrac

%
%                ( 0.4*pi / (T*Wc) )^(2N)
% H_stop_pass =  ------------------------
%                ( 0.25pi / (T*Wc) )^(2N)
%
% H_stop_pass = ( 0.4pi / 0.25pi )^(2N)
%
% log10(H_stop_pass) = 2N * log10( 0.4pi / 0.25pi )
%

N(a) = 1/2 * log10(H_stop_pass) / log10( (0.4*pi) / (0.25*pi) )
% N must be an integer number so round up
N(a) = ceil(N(a))

% 
% H_passWFrac = ( 0.25*pi / W_c )^(2N)
%
W_c(a) = ( (0.25*pi)^(2*N(a)) / H_passWFrac )^(1/(2*N(a)))
% and get it in terms of pi
W_c_pi(a) = W_c(a) / pi

%%%
% Transfer Function for Nth Order Analog Butterworth filter
% normalized to W_c = 1
%%%
% calcualte the zeros, poles, and gain
[r,p,k] = butter(N(a), 1, 's')

%%%
% creation of the transfer function for vaiable number of poles
%              1
% H0 = ------------------
%      (s - p_1)(s - p_2)
%%%
s = tf('s')
H0 = 1
for b = 1:numel(p)
   H0 = H0 / (s - p(b))
end

% shift the cutoff frequency to W_c/T calcuated above
T = 1
s_shift = s / (W_c(a) / T)
% recreate the transfer function using shifted s
H1 = 1
for c = 1:numel(p)
   H1 = H1 / (s_shift - p(c))
end

% recalcuate the roots and poles for new transfer function
[r,p,k] = residue(H1.num{1}, H1.den{1})

%%%
% take the inverse laplace transform
%%%
syms x
f_t{a} = 0
for d = 1:numel(r)
    f_t{a} = f_t{a} + ilaplace( r(d) / (x - p(d)) )
end

%%%
% take the Z transfrom
%%%
syms t z
h_z{a} = ztrans(f_t{a}, t, z)

%%%
% Transfrom into a highpass filter
%%%
alpha(a) = -(cos((W_c(a) + 0.7*pi)/2))/(cos((W_c(a) - 0.7*pi)/2))

%%%
% substitute the z with the highpass version
%%%
h_z_hp{a} = subs(h_z{a}, z, -((1 + alpha(a)*(1/z))/(alpha(a) + (1/z))))

%%%
% plotting the magnitude resonses of the Z transform
%%%
plot = 1
if plot
frac{a} = collect(simplifyFraction(vpa(h_z{a}, 12)))
[num_rev{a} den_rev{a}] = numden(frac{a})

num_coeffs{a} = coeffs(vpa(num_rev{a}))
den_coeffs{a} = coeffs(vpa(den_rev{a}))

num{a} = double(fliplr(num_coeffs{a}))
den{a} = double(fliplr(den_coeffs{a}))

figure(2*a - 1)
freqz(num{a}, den{a})
title(['Discrete Low-pass Butterworth Filter of Order ', num2str(N(a))])
ax1 = subplot(2,1,1);
limsy=get(gca,'YLim');
set(gca,'Ylim',[limsy(1) 10])
ax2 = subplot(2,1,2);
limsy=get(gca,'YLim');
set(gca,'Ylim',[limsy(1) 0])
print(['Discrete Low-pass Butterworth Filter of Order ', num2str(N(a))], '-dpng' )

frac_hp{a} = collect(simplifyFraction(vpa(h_z_hp{a}, 12)))
[num_rev_hp{a} den_rev_hp{a}] = numden(frac_hp{a})

num_coeffs_hp{a} = coeffs(vpa(num_rev_hp{a}))
den_coeffs_hp{a} = coeffs(vpa(den_rev_hp{a}))

num_hp{a} = double(fliplr(num_coeffs_hp{a}))
den_hp{a} = double(fliplr(den_coeffs_hp{a}))

figure(2*a)
freqz(num_hp{a}, den_hp{a})
title(['Discrete High-pass Butterworth Filter of Order ', num2str(N(a))])
ax1 = subplot(2,1,1);
limsy=get(gca,'YLim');
set(gca,'Ylim',[limsy(1) 10])
ax2 = subplot(2,1,2);
%%% shift the highpass phase by n2Pi
axesObjs = gca;  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
ydata = get(dataObjs, 'YData');
set(dataObjs, 'YData', ydata + (a-1)*360)
%%%
print(['Discrete High-pass Butterworth Filter of Order ', num2str(N(a))], '-dpng' )
end

end

%%%
% output the transfroms to a file
%%%
fileID = fopen('output.txt','wt');

for d = 1:numel(delta_l)
fprintf(fileID,'delta_l = %f \n', delta_l(d));
fprintf(fileID,'f_t = %s \n', char(vpa(f_t{d}, 6)));
fprintf(fileID,'h_z = %s \n', char(vpa(h_z{d}, 6)));
fprintf(fileID,'h_z_highpass = %s \n\n', char(vpa(h_z_hp{d}, 6)));
end

fclose(fileID);
