format compact
clear
clc
%close all
clf reset

%warning('off','all')

record_video = false;

flip_frame_length = 22;

if record_video == true
    v = VideoWriter("panel transition test", 'MPEG-4');
    v.FrameRate = 30;
    open(v);
end

img_full = imread("test_image.jpg");
img_full = flip(img_full);
[img_x,img_z] = meshgrid(1:width(img_full), 1:height(img_full));
img_y = ones(height(img_x),width(img_x));

img_background = uint8(zeros(height(img_z),width(img_z),3));
img_green = img_background;
img_green(:,:,2) = uint8(255*ones(height(img_z),width(img_z)));

panel_cellwidth = 70; %pixels
panel_xedges = -panel_cellwidth+1 : panel_cellwidth : width(img_full) + panel_cellwidth*0.1;
panel_zedges = -panel_cellwidth+1 : panel_cellwidth : height(img_full) + panel_cellwidth*0.1;

for n=1:length(panel_zedges)-1
    for m=1:length(panel_xedges)-1
        panel_centercoords(n,m,:) = [mean([panel_zedges(n),panel_zedges(n+1)]), mean([panel_xedges(m),panel_xedges(m+1)])];
    end
end

x_lin = img_x(1,1:end);
z_lin = img_z(1:end,1).';

maxwait = max([length(panel_zedges)-1, length(panel_xedges)-1]);

%panel rotation PID control
PID_time = 1:flip_frame_length;
k = 0.5;
numerator = [1,k];
denominator = [3,1.5,k];
sys = tf(numerator,denominator);
[stepsys,time]=step(sys, flip_frame_length);
stepsys = stepsys.*180;
time = time+2;
PID_theta = interp1(time,stepsys,PID_time,"linear");
PID_theta = PID_theta;
PID_theta(end) = 180;
PID_theta(1) = 0;
PID_theta(2) = mean([PID_theta(1),PID_theta(3)]);

%rotation_theta_seq = linspace(0,180,flip_frame_length);
rotation_theta_seq = PID_theta;
rotation_theta_seq = [rotation_theta_seq, repelem(180,maxwait+flip_frame_length) , 180-rotation_theta_seq];
animation_frames = round(length(rotation_theta_seq)*1.4);

%generate rotation sequence per frame
panel_rotation_sequence = zeros(length(panel_zedges)-1, length(panel_xedges)-1, animation_frames);
ind_xlim=0;
ind_zlim=0;
diag_prev=[];
for ind_frame = 1:animation_frames
    %diag_pattern = zeros(ind_zlim, ind_xlim);
    for n = 1:length(panel_zedges)-1
        for m = 1:length(panel_xedges)-1
            diag_pattern(n, m) = max(0, ind_zlim - (n + m - 2));
            if diag_pattern(n, m) > length(rotation_theta_seq)
                diag_pattern(n, m) = length(rotation_theta_seq);
            end
        end
    end
    
    diag_pattern = fliplr(flipud(diag_pattern));

    for n = 1:length(panel_zedges)-1
        for m = 1:length(panel_xedges)-1
            if diag_pattern(n,m) > 0
                panel_rotation_sequence(n,m,ind_frame) = rotation_theta_seq(diag_pattern(n,m));
            end
        end
    end

    ind_zlim = ind_zlim+1;
    ind_xlim = ind_xlim+1;
end



normal_vec = [0,0,1];
max_x_prev = 0;
max_z_prev = 0;
for ind_frame = 1:animation_frames
    
    
    scatter(nan,nan) %clear screen

    hold on
    grid on

    %find a specific panel set and draw it
    for n=1:length(panel_zedges)-1
        for m=1:length(panel_xedges)-1
    
            theta_desired = panel_rotation_sequence(n,m,ind_frame); %degrees
            if theta_desired ~=0 && theta_desired ~= 180
                theta_desired = theta_desired + rand()*5;
            end
    
            %find bounding indexes for panel
            [~,ind_xlower] = min(abs(x_lin-panel_xedges(m)));
            [~,ind_xupper] = min(abs(x_lin-panel_xedges(m+1)));
            [~,ind_zlower] = min(abs(z_lin-panel_xedges(n)));
            [~,ind_zupper] = min(abs(z_lin-panel_xedges(n+1)));
    
            panel_x = img_x(ind_zlower:ind_zupper,ind_xlower:ind_xupper) - panel_centercoords(n,m,2);
            panel_y = -1 .* ones(length(ind_zlower:ind_zupper), length(ind_xlower:ind_xupper));
            panel_z = img_z(ind_zlower:ind_zupper,ind_xlower:ind_xupper) - panel_centercoords(n,m,1);
            
            panel_img = img_full(ind_zlower:ind_zupper, ind_xlower:ind_xupper,:);
            panel_img = fliplr(panel_img);
    
            panel_y_green = zeros(length(ind_zlower:ind_zupper), length(ind_xlower:ind_xupper));
            panel_img_green = img_green(ind_zlower:ind_zupper, ind_xlower:ind_xupper,:);
            panel_img_background = img_background(ind_zlower:ind_zupper, ind_xlower:ind_xupper,:);

            if ind_frame == 1
                %find image limits

                if ind_zupper > max_z_prev
                    max_z_prev = ind_zupper;
                end
                if ind_xupper > max_x_prev
                    max_x_prev = ind_xupper;
                end
            end

            warp(panel_x+panel_centercoords(n,m,2), panel_y-100, panel_z+panel_centercoords(n,m,1), panel_img_background)

            %for regular
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

            %for regular
            panel_x_rot=[];
            panel_y_rot=[];
            panel_z_rot=[];
            for ind_z = 1:height(panel_z)
                for ind_x = 1:width(panel_x)
                    coord_vec = [panel_x(ind_z, ind_x), panel_y_green(ind_z, ind_x), panel_z(ind_z, ind_x)];
                    coord_new = coord_vec*cosd(theta_desired) + cross(normal_vec, coord_vec)*sind(theta_desired) + normal_vec*dot(normal_vec, coord_vec)*(1 - cosd(theta_desired));
    
                    panel_x_rot(ind_z, ind_x) = coord_new(1);
                    panel_y_rot(ind_z, ind_x) = coord_new(2);
                    panel_z_rot(ind_z, ind_x) = coord_new(3);
                end
            end
            panel_x_rot = panel_x_rot + panel_centercoords(n,m,2);
            panel_z_rot = panel_z_rot + panel_centercoords(n,m,1);
            warp(panel_x_rot, panel_y_rot, panel_z_rot, panel_img_green)
    
        end
    end

    for n=1:length(panel_zedges)-1
        for m=1:length(panel_xedges)-1
            if panel_rotation_sequence(n,m,ind_frame) ~= 0 && panel_rotation_sequence(n,m,ind_frame) ~= 180
                scatter3(panel_centercoords(n,m,2), 1, panel_centercoords(n,m,1), "k+")
            end
        end
    end
    
    axis equal tight
    view(30,30)
    xlim([0,max_x_prev])
    ylim([-1e3,1e3])
    zlim([0,max_z_prev])

    set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
    hold off

    drawnow()
    

    if record_video == true
        frame = getframe(gcf);
        % imwrite(frame.cdata, 'test.jpg')

        writeVideo(v,frame.cdata)
        fprintf("-----------\n")
        fprintf("written frame %i.\n", ind_frame)
    end

end

if record_video == true
    close(v);
end

sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);
