format compact
clear
clc
%close all
clf reset


% textentry = [
% "$\Delta_{\sigma} \frac{\gamma}{\alpha} \quad a$"
% ];
% 
% 
% cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
% 
% textobj = text(0,0,textentry, 'VerticalAlignment', 'middle', interpreter="latex", FontSize=50);
% 
% axis tight equal
% 
% text_textent= get(textobj, 'Extent');
% xmargin = 0.1;
% ymargin = 0.5;
% xlim([text_textent(1) - xmargin*text_textent(3), text_textent(1) + text_textent(3) + xmargin*text_textent(3)]);
% ylim([text_textent(2) - ymargin*text_textent(4), text_textent(2) + text_textent(4) + ymargin*text_textent(4)]);
% 
% frame_plot = getframe(gca);
% img = flip(frame_plot.cdata);
% img = single(imresize(img, 0.25));
% img = img./max(max(img));
% img_processed = squeeze(img(:,:,3));
% img_processed = 1-img_processed;
% 
% imwrite(img_processed, "image.png")

img_full = imread("test_image.jpg");
img_full = flip(img_full);
[img_x,img_z] = meshgrid(-width(img_full)/2:width(img_full)/2, 1:height(img_full));
img_y = ones(height(img_x),width(img_x))+height(img_full);

theta_desired = 0; %degrees
normal_vec = [0,0,1];

for n=1:height(img_full)
    for m=1:width(img_full)
        coord_vec = [img_x(n,m), img_y(n,m), img_z(n,m)];
        coord_new = coord_vec*cosd(theta_desired) + cross(normal_vec, coord_vec)*sind(theta_desired) + normal_vec*dot(normal_vec, coord_vec)*(1 - cosd(theta_desired));
        img_x_rot(n,m) = coord_new(1);
        img_y_rot(n,m) = coord_new(2);
        img_z_rot(n,m) = coord_new(3);
    end
end

scatter(nan,nan)
hold on
grid on
warp(img_x_rot,img_y_rot,img_z_rot,img_full)
axis equal
%view(30,30)
view(0,0)

xlim([-1e3,1e3])
ylim([-1e3,1e3])
zlim([-1e3,1e3])