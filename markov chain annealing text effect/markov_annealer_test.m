format long
clear
clc
%close all
clf reset

record_video = true;

if record_video == true
    v = VideoWriter("traq_markov_fadein", 'MPEG-4');
    v.FrameRate = 30;
    open(v);
end

image_raw = imread("image.png");
image_raw = flip(single(image_raw));
image_raw = squeeze(image_raw(:,:,1));
image_raw = image_raw./max(max(image_raw));
image_raw = 1-image_raw;
image_raw = imresize(image_raw,[1080/3,nan]);
image_raw(:,1:2) = 0;
image_raw(end-1:end,:) = 0;
image_raw(image_raw>1) = 1;
image_raw(image_raw<0) = 0;
target_image = round(image_raw);

target_comp = bwconncomp(target_image,4);
num_cells = target_comp.NumObjects;
img_dim = target_comp.ImageSize;

target_celladj = ones(img_dim);
for n=1:num_cells
    target_celladj(target_comp.PixelIdxList{n}) = num_cells-n+2;
end
% target_celladj = flipud(target_celladj);

inital_landscape_bandind = round(linspace(1,img_dim(2),num_cells+3)); %one extra cell for background region
intial_landscape = ones(img_dim);
for n=3:length(inital_landscape_bandind)-1
    %intial_landscape(inital_landscape_bandind(n-1):inital_landscape_bandind(n),1:img_dim(2)) = n-2;
    intial_landscape(:,inital_landscape_bandind(n-1):inital_landscape_bandind(n)) = length(inital_landscape_bandind)-n + 1;
    %intial_landscape(:,inital_landscape_bandind(n)-3:inital_landscape_bandind(n)) = 1;
end
intial_landscape(1:end*0.05,:) = 1;
intial_landscape(end*0.95:end,:) = 1;



%landscape_current = imdilate(target_celladj, strel('disk', 10));
landscape_current = intial_landscape;
landscape_prev = landscape_current;
fitness_prev = -inf;

% delete(gcp('nocreate'));
% parpool('local',8);

%cell flipping control
ind_1 = 1;
ind_sf = 1;
fitness_current = double(0);
f_countdown = 1e3;
landscape_backup = intial_landscape;
rand_modifier = 0;

if record_video == true
    frame = getframe(gcf);
    frame.cdata = imresize(frame.cdata, [1895,3840]);
    writeVideo(v,frame)
end

while fitness_current < 1

    for n=1:num_cells+1

        landscape_squeezed = reshape(landscape_current,[],1);
        cell_spec = zeros(img_dim);
        cell_spec = landscape_squeezed == n;
        cell_spec(cell_spec==1) = n;
        cell_spec = reshape(cell_spec,img_dim);
    
        %exterior region we can flip to this specific colour
        exterior_periregion = bwperim(cell_spec,4);
        exterior_periregion = imdilate(exterior_periregion, strel('disk', 1)); %adds a 1 pixel border
        exterior_periregion = ~cell_spec & exterior_periregion;
        exterior_indexregions = landscape_current.*exterior_periregion;
    
        %interior region we can flip to the bordering colour
        interior_periregion = bwperim(cell_spec,4);
        interior_periregion = cell_spec & interior_periregion;
        interior_indexregions = landscape_current.*interior_periregion;
        interior_nearbymask = interior_periregion .* imdilate(exterior_indexregions, strel('disk', 1));

        check_good = false;
        while ~check_good

            % if rand() < 0.5
            %     %flip interior cell to exterior index
            %     peri_spec_inds = find(interior_periregion);
            %     rand_spec_cellind = peri_spec_inds(randi(numel(peri_spec_inds)));
            %     [ind_r, ind_c] = ind2sub(size(interior_periregion), rand_spec_cellind);
            %     %flip to whatever cell index is in nearby_mask
            %     landscape_current(ind_r, ind_c) = interior_nearbymask(ind_r, ind_c);
            % else
            %     %flip exterior cell to interior index
            %     peri_spec_inds = find(exterior_periregion);
            %     rand_spec_cellind = peri_spec_inds(randi(numel(peri_spec_inds)));
            %     [ind_r, ind_c] = ind2sub(size(exterior_periregion), rand_spec_cellind);
            %     %flip to n
            %     landscape_current(ind_r, ind_c) = n;
            % end
            
            peri_spec_inds = find(interior_periregion);
            rand_spec_cellind = peri_spec_inds(randi(numel(peri_spec_inds)));
            [ind_r, ind_c] = ind2sub(size(interior_periregion), rand_spec_cellind);
            %flip to whatever cell index is in nearby_mask
            landscape_current(ind_r, ind_c) = interior_nearbymask(ind_r, ind_c);

            check_good = true;
            if n==1
                if ind_r == 1 || ind_r == img_dim(1) || ind_c == 1 || ind_c == img_dim(2)
                    check_good = false;
                end
            end
        end
        
        %checking if any cells have been broken (if so immediately rejected)
        %cell_broken_list = zeros(num_cells,1);
        cell_broken = false;
        landscape_squeezed_new = reshape(landscape_current,[],1);
        for m = 1:num_cells+1 
            cell_spec = landscape_squeezed_new == m;
            cell_spec(cell_spec==1) = m;
            cell_spec = reshape(cell_spec,img_dim);
            cell_spec_conn = bwconncomp(cell_spec,4);
            if cell_spec_conn.NumObjects > 1 || cell_spec_conn.NumObjects == 0
                %cell has been broken
                %cell_broken_list(m) = 1;
                cell_broken = true;
                break
            end
        end
        
        %descision making
        progress_plot = false;
        if cell_broken
            landscape_current = landscape_prev;
        else
            fitness_current = 1 / immse(target_celladj.^2 - 1, landscape_current.^2 - 1);
            if fitness_current >= fitness_prev
                landscape_prev = landscape_current;
                fitness_prev = fitness_current;
                progress_plot = true;
                
                %reset countdown
                f_countdown = 1e3;
                landscape_backup = landscape_current;
            else
                f_countdown = f_countdown - 1;
                if rand() < 0.02 + rand_modifier
                    landscape_prev = landscape_current; %accept new map anyway
                    %progress_plot = true;
                else
                    landscape_current = landscape_prev; %reject map
                end
            end
        end
   
        if f_countdown < 0
            %stops us getting stuck
            f_countdown = rand()*50 + 10;

            for m = 1:2
                if rand()<0.05
                    f_countdown = f_countdown*10;
                end
            end

            landscape_current = landscape_backup;
            rand_modifier = rand()/2 - 0.1;
            if rand_modifier > 0.3
                rand_modifier = 0;
            end
            fprintf("rolling back. \n")
        end
    
        if progress_plot || ind_1 < 10
            ind_sf = ind_sf+1;
        end
        ind_1 = ind_1 + 1;

        if rem(ind_sf,200) == 0 || rem(ind_sf,1e5) == 0 || ind_1 < 10
            cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
            colormap(cmap)
            surf(landscape_current,EdgeColor="none")
            view([0,90])
            axis equal tight
            set(gca,'XTickLabel',[]);
            set(gca,'YTickLabel',[]);
            clim([1,12])
            fprintf("\n step forward taken (%i, %i, %3.9f) " + string(datetime("now"))+ ". \n",ind_1, ind_sf, fitness_current)
            drawnow()
            ind_sf = ind_sf+1;

            if record_video == true
                frame = getframe(gcf);
                frame.cdata = imresize(frame.cdata, [1895,3840]);
                writeVideo(v,frame)
            end

        end
    end
end

if record_video == true
    close(v);
end



% cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
% colormap(cmap)
% hold on
% axis equal tight
% surf(target_celladj,EdgeColor="none")
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);
