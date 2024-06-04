## Run MouseActivity in MATLAB
1. Download MouseActivity5.m, xlwrite.m and the folder poi_library into a folder on your computer;
2. Open MATLAB R2018b;
3. Go to the folder containing MouseActivity5.m in MATLAB;
4. Double-click the file MouseActivity5.m;
5. Click "Run", and follow the onscreen instructions to analyze video files. 
6. Please note, certain video file format may not work. It has been tested that .mp4 and .mov files are good to go. 

## Run MouseActivity GUI .mlappinstall for MatLab

# Installation
To use run this gui app, either MouseActivityGui.mlapp or MouseActivity.mlappinstall could be used. 
For MouseActivityGui.mlapp:
```bash
Open MouseActivityGui.mlapp in MatLab and tap 'Run' in the toolbar.
```
For MouseActivity.mlappinstall:
```bash
Open MatLab install the app by using 'Install App' under the 'APPS' menu.
```

## Tips 
A 'draw next' pushbutton is added to the window, to draw next rectangle area for mouse activity. 

For white mouse in black background, please 
1) Comment out the following lines:

imshow(bwareaopen(im2bw(255 - frame100, kk(k)),50) * 255);

frame2 = bwareaopen(im2bw(255-read(video_obj,x(jj)), threshold),objpixels) * 255; % remove any object smaller than 50 pixels

3) Uncomment the following lines:

%imshow(bwareaopen(im2bw(frame100, kk(k)),50) * 255); % use this line for white mouse

%frame2 = bwareaopen(im2bw(read(video_obj,x(jj)), threshold),objpixels) * 255; % use this for white mouse
