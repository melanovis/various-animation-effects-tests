format compact
clear
clc
%close all
clf reset

image_raw = imread("image.png");
image_raw = flip(single(image_raw));
image_raw = squeeze(image_raw(:,:,1));
image_raw = image_raw./max(max(image_raw));
image_raw = 1-image_raw;
image_raw = imresize(image_raw,[1080/4.75,nan]);
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

inital_landscape_bandind = round(linspace(1,img_dim(2),num_cells+3)); %one extra cell for background region
intial_landscape = ones(img_dim);
for n=3:length(inital_landscape_bandind)-1
    intial_landscape(:,inital_landscape_bandind(n-1):inital_landscape_bandind(n)) = length(inital_landscape_bandind)-n + 1;
end
intial_landscape(1:end/4,:)=1;
intial_landscape(end*3/4:end,:)=1;


landscape_current = intial_landscape;
landscape_prev = intial_landscape;

while true

    %cell flipping control
    pix_changelist = zeros([num_cells,3]);
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


        if rand() > 0.5
            %flip interior cell to exterior index

            %find a random exterior index
            peri_spec_inds = find(interior_periregion);
            rand_spec_cellind = peri_spec_inds(randi(numel(peri_spec_inds)));
            [ind_r, ind_c] = ind2sub(size(interior_periregion), rand_spec_cellind);

            %flip to whatever cell index is in nearby_mask
            pix_changelist(n,:) = [ind_r, ind_c, interior_nearbymask(ind_r, ind_c)];
        else
            %flip exterior cell to interior index

            %find a random interior index
            peri_spec_inds = find(exterior_periregion);
            rand_spec_cellind = peri_spec_inds(randi(numel(peri_spec_inds)));
            [ind_r, ind_c] = ind2sub(size(exterior_periregion), rand_spec_cellind);

            %flip to n
            pix_changelist(n,:) = [ind_r, ind_c, n];
        end
    end
    
    for n=1:height(pix_changelist)
        landscape_current(pix_changelist(n,1), pix_changelist(n,2)) = pix_changelist(n,3);
    end

    cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
    colormap(cmap)
    subplot(2,1,1)
    imagesc(exterior_periregion)
    axis equal tight
    subplot(2,1,2)
    imagesc(landscape_current)
    axis equal tight
    clim([1,12])
    fprintf("----\n")
    drawnow()

    %checking if any cells have been broken (if so immediately rejected)
    some_cell_broken = false;
    landscape_squeezed = reshape(landscape_current,[],1);
    for n=1:num_cells+1 
        cell_spec = landscape_squeezed == n;
        cell_spec(cell_spec==1) = n;
        cell_spec = reshape(cell_spec,img_dim);
        cell_spec_conn = bwconncomp(cell_spec,4);
        if cell_spec_conn.NumObjects > 1 || cell_spec_conn.NumObjects == 0
            %cell has been broken
            some_cell_broken = true;
            break
        end
    end

    %descision making
    if some_cell_broken
        landscape_current = landscape_prev;
    else
        landscape_prev = landscape_current;
    end

end

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
colormap(cmap)

subplot(2,1,1)
hold on
grid on
axis equal tight  
surf(target_celladj,EdgeColor="none")
clim([1,12])

subplot(2,1,2)
hold on
grid on
axis equal tight
surf(landscape_current,EdgeColor="none")
clim([1,12])