% write_line_profile.m

% %--- Example of file format for a 2D line profile with N values:
% ((profilename line N)
% (x
% 4.00000E+00 4.00000E+00 4.00000E+00 4.00000E+00
% 4.00000E+00 4.00000E+00 4.00000E+00 4.00000E+00 )
% (y
% 1.06443E-03 3.19485E-03 5.33020E-03 7.47418E-03
% 2.90494E-01 3.31222E-01 3.84519E-01 4.57471E-01 )
% (dataname
% 5.47866E+00 6.59870E+00 7.05731E+00 7.40079E+00
% 1.01674E+01 1.01656E+01 1.01637E+01 1.01616E+01 )
% )
% %----------------------------------------------------

clear all
close all

%--- Parameters 
D    = 0.03;  
Umax = 4;       

%--- Define profile
N = 21;   
y = linspace(0, D, N);
x = zeros(size(y));

u1 = -Umax*(y-D).*(y+D)/D^2;
u2 = u1 + (0.2*Umax)*(rand(size(u1))-0.5)*2;

figure, hold on, box on, grid on, xlabel ux, ylabel y,
plot(u1, y, '--'),
plot(u2, y, '--'),

%--- Write file 
fid = fopen(['line_profile.txt'], 'w');    

fprintf(fid, '((velocity-uy line %u) \n', length(u2));

fprintf(fid, '(x \n', []);
fprintf(fid, '%6.12g ', x);
fprintf(fid, ')\n', []);

fprintf(fid, '(y \n', []);
fprintf(fid, '%6.12g ', y);
fprintf(fid, ')\n', []);
 
fprintf(fid, '(u \n', []);
fprintf(fid, '%6.12g ', u2);
fprintf(fid, ')\n', []);

fprintf(fid, ')\n', []);
fclose(fid);
