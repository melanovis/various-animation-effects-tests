format compact
clear
clc
%close all
clf reset


img_full = imread("test_image.jpg");
img_full = flip(img_full);
[img_x,img_z] = meshgrid(1:width(img_full), 1:height(img_full));
img_y = ones(height(img_x),width(img_x));

panel_cellwidth = 80; %pixels
panel_xedges = -panel_cellwidth/10 : panel_cellwidth : width(img_full) + panel_cellwidth;
panel_zedges = -panel_cellwidth/10 : panel_cellwidth : height(img_full) + panel_cellwidth;

for n=1:length(panel_zedges)-1
    for m=1:length(panel_xedges)-1
        panel_centercoords(n,m,:) = [mean([panel_zedges(n),panel_zedges(n+1)]), mean([panel_xedges(m),panel_xedges(m+1)])];
    end
end

x_lin = img_x(1,1:end);
z_lin = img_z(1:end,1).';


theta_desired = 20; %degrees
normal_vec = [0.05,0,1];

hold on
grid on
%find a specific panel set and draw it
for n=1:length(panel_zedges)-1
    for m=1:length(panel_xedges)-1

        %find bounding indexes for panel
        [~,ind_xlower] = min(abs(x_lin-panel_xedges(m)));
        [~,ind_xupper] = min(abs(x_lin-panel_xedges(m+1)));
        [~,ind_zlower] = min(abs(z_lin-panel_xedges(n)));
        [~,ind_zupper] = min(abs(z_lin-panel_xedges(n+1)));

        panel_x = img_x(ind_zlower:ind_zupper,ind_xlower:ind_xupper) - panel_centercoords(n,m,2);
        panel_y = ones(length(ind_zlower:ind_zupper), length(ind_xlower:ind_xupper));
        panel_z = img_z(ind_zlower:ind_zupper,ind_xlower:ind_xupper) - panel_centercoords(n,m,1);
        panel_img = img_full(ind_zlower:ind_zupper, ind_xlower:ind_xupper,:);

        % panel_x = img_x(ind_zlower:ind_zupper,ind_xlower:ind_xupper);
        % panel_y = zeros(length(ind_zlower:ind_zupper), length(ind_xlower:ind_xupper));
        % panel_z = img_z(ind_zlower:ind_zupper,ind_xlower:ind_xupper);
        % panel_img = img_full(ind_zlower:ind_zupper, ind_xlower:ind_xupper,:);

        % warp(panel_x, panel_y, panel_z, panel_img)
        % view(0,0)
        drawnow()

        panel_x_rot=[];
        panel_y_rot=[];
        panel_z_rot=[];
        for ind_z = 1:height(panel_z)
            for ind_x = 1:width(panel_x)
                coord_vec = [panel_x(ind_z, ind_x), panel_y(ind_z, ind_x), panel_z(ind_z, ind_x)];
                coord_new = coord_vec*cosd(theta_desired) + cross(normal_vec, coord_vec)*sind(theta_desired) + normal_vec*dot(normal_vec, coord_vec)*(1 - cosd(theta_desired));

                panel_x_rot(ind_z, ind_x) = coord_new(1);
                panel_y_rot(ind_z, ind_x) = coord_new(2);
                panel_z_rot(ind_z, ind_x) = coord_new(3);
            end
        end
        panel_x_rot = panel_x_rot + panel_centercoords(n,m,2);
        panel_z_rot = panel_z_rot + panel_centercoords(n,m,1);
        warp(panel_x_rot, panel_y_rot, panel_z_rot, panel_img)

    end
end



scatter(nan,nan)
hold on
grid on
%warp(img_x,img_y,img_z,img_full)
%warp(new_x,new_y,new_z,img_pixcut)
for n=1:length(panel_zedges)-1
    for m=1:length(panel_xedges)-1
        scatter3(panel_centercoords(n,m,2), 1, panel_centercoords(n,m,1), "k+")
    end
end
axis equal tight
%view(30,30)
view(0,0)
xlim([min(panel_xedges)-50,max(panel_xedges)+50])
zlim([min(panel_zedges)-50,max(panel_zedges)+50])

% xlim([-100,1e3])
% ylim([-0.5e3,0.5e3])
% zlim([-100,1e3])