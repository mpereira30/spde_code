fig_name = 'open_loop_surf.fig';
fig = openfig(fig_name);
hold on;

x1 = [0.1719, 0.2109];
x1 = [x1;x1];

x2 = [0.4688, 0.5156];
x2 = [x2; x2];

x3 = [0.7734, 0.8125];
x3 = [x3; x3];

y = [0; 100];
y = [y, y];

z1 = 5 * ones(2, 2);
z2 = 2.5 * ones(2, 2);

surf(x1, y, z1, 'FaceColor', 'r');
surf(x2, y, z2, 'FaceColor', 'r');
surf(x3, y, z1, 'FaceColor', 'r');