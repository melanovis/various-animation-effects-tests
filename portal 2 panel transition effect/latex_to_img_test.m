format compact
clear
clc
%close all
clf reset


textentry = [
"$line_1$"
"$line_2 \quad \Delta_{\sigma}$"
"$line_3 \quad \vec{r}^* = \vec{r} \cos (\psi) + (\hat{n} \times \vec{r}) \sin (\psi) + \hat{n}(\hat{n} \cdot \vec{r})(1-\cos (\psi))$"
];


cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

textobj = text(0,0,textentry, 'VerticalAlignment', 'middle', interpreter="latex", FontSize=50);
axis tight equal

text_textent= get(textobj, 'Extent');
xmargin = 0.1;
ymargin = 0.5;
xlim([text_textent(1) - xmargin*text_textent(3), text_textent(1) + text_textent(3) + xmargin*text_textent(3)]);
ylim([text_textent(2) - ymargin*text_textent(4), text_textent(2) + text_textent(4) + ymargin*text_textent(4)]);

frame_plot = getframe(gca);
img = flip(frame_plot.cdata);
img = single(imresize(img, 0.5));
img = img./max(max(img));
img_processed = squeeze(img(:,:,3));

scatter(nan,nan)
hold on
grid on
imagesc(img_processed)
axis equal tight
colormap("gray")