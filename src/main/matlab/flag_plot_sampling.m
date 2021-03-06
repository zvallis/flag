function flag_plot_sampling(L, P, R)

if P > 1
    [rs, thetas, phis] = flag_sampling(L, P, R);
else
    [thetas, phis] = ssht_sampling(L);
    rs = R;
end

x = [];
y = [];
z = [];
s = [];
for i = 1:P
    for j = 1:L
        for k = 1:2*L-1
            x = [x rs(i)*cos(phis(k))*sin(thetas(j))];
            y = [y rs(i)*sin(phis(k))*sin(thetas(j))];
            z = [z rs(i)*cos(thetas(j))];
            s = [s i];
        end
    end   
end
c = numel(x):-1:1;

figure('Position',[1 1 600 600],'Color',[1 1 1]);
title('Fourier-Laguerre transform: 3D sampling theorem', 'FontSize', 20);

if P > 1    
    h = scatter3(x,y,z,0.5*s,c,'filled');
else
    f = zeros(L,2*L-1);
    ssht_plot_sphere(f,L,'PlotSamples', true);
    h = scatter3(x,y,z,30,'filled');
end
set(gca, 'visible', 'off'); 
view(60,30);
zoom(1.5);

end