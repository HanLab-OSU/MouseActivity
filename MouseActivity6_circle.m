 function MouseActivity6_circle()
%% Mouse Activity Analyzer version 6
% Author: Renzhi Han
% Affiliation: Indiana University School of Medicine
% This code analyzes mice in a circular open field, each in its own sub-arena.
% It stores the following information:
% Position
% Pathlength
% Accumulative travel distance
% Mouse size (area)
% Major axis length assuming an oval shape for a mouse
% Minor axis length assuming an oval shape for a mouse
% Orientation angle (theta)
% Eccentricity (the ratio of foci distance to major axis length, between 0-1)
%% code starts here

clc; % clear the comman window
close all; % close all figures
clear; % erase all existing variables

    [filename, dirpath] = uigetfile('*.mov;*.wmv;*.mp4;*.avi','Open video');
    
    if isequal(filename, 0) || isequal(dirpath, 0)
        logf('Cancel opening video');
        return;
    end
    
    filepath = [dirpath filename];
    
    [u1, video_name, u2] = fileparts(filepath);
                
    result_name = [video_name '_result.mat'];
    resultpath = [dirpath video_name '/' result_name];
            
if exist(resultpath, 'file') ~= 2
 
    logf('Opening the video: %s', filepath);
                
    try
        video_obj = VideoReader(filepath);  
    catch exception
        close_handle(dialog);
        msgbox(['Error opening video file. Message: ' exception.message], ...
               'Open video', 'error');
        
        set_status('Error when opening video file', 'button_cancel_red', true);
        return;
    end    
        nframes = video_obj.NumberOfFrames;
        frame100 = read(video_obj, round(nframes/2));
        videosize = [video_obj.Width,video_obj.Height];
        frameRate = video_obj.FrameRate;
        duration = video_obj.Duration;
        
        if(~exist([dirpath video_name],'dir'))
            mkdir([dirpath video_name]);
        end
        imgfilename1 = [dirpath video_name '/' video_name '_threshold'];
        
        %figure('Renderer', 'painters', 'Position', [10 10 900 600]) 
        subplot(2,3,1)
        imshow(frame100); 
%         text(10, -200,'Please click on the image to choose the threshold','Color','r','FontSize',20,'FontWeight','bold');
        
        
        kk = [0.5 0.6 0.7, 0.75, 0.8];
        
        for k = 1:length(kk)
            subplot(2,3,k+1)
            imshow(bwareaopen(im2bw(255 - frame100, kk(k)),50) * 255);
            set(gca,'tag',num2str(kk(k)))
            title(['Threshold: ' num2str(kk(k))]);            
        end
        set(gcf, 'Position', [0, 0.04, 900, 600]);
        w = waitforbuttonpress;
        if w == 0
            th = get(gca,'tag');
        end        
        %saveas(gcf,[imgfilename1 '.tif']);
        print(gcf,[imgfilename1 '.tif'],'-dtiff','-r300');

        prompt = {'Enter the threshold:','Enter the minimal pixels of object:','Mouse number:', 'StartFrame to be analyzed:', 'LastFrame to be analyzed:', 'Frame steps:'};
        dlg_title = 'Collect movie processing parameters';
        num_lines = 1;
        defaultans = {th,'100','1','1000',num2str(nframes),'1'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        threshold = str2double(answer{1}); % default 0.70; % the threshold is a key parameter for successful segmentation
        objpixels = floor(str2double(answer{2})); % the minimal pixel size of the mouse, default is 50
        MouseN = floor(str2double(answer{3})); % the number of mouse in the movie, default is 2
        StartFrame = floor(str2double(answer{4})); % the first frame to be analyzed, default is 1
        LastFrame = floor(str2double(answer{5})); % the last frame to be analyzed, default is the maximum frame number
        if LastFrame > nframes
            LastFrame = nframes;
        end
        Step = floor(str2double(answer{6})); % the frame steps to be analyzed (e.g. analyze every 'x' frame), default is 1
        close all;
        
        %bg = read(video_obj,300);
        figure('Position', [10 10 900 900]) 
        imshow(read(video_obj,StartFrame));   
        
        arenaD = cell(MouseN,1);
        xa = cell(MouseN,1);
        ya = cell(MouseN,1);
        arena_r = cell(MouseN,1);

        x0 = 662;
        y0 = 363;
        r0 = 150;
        for n = 1:MouseN
            logf(['please draw the arena for mouse' num2str(n)]);
            %arena = drawrectangle('Position',[20 20 w0 h0], 'LineWidth',1);;
            arena = drawcircle('Center', [x0, y0], 'Radius', r0, 'StripeColor','red','LineWidth', 1);
            h = customWait(arena);
            logf(num2str(h));
            xa{n} = floor(arena.Center(1));
            ya{n} = floor(arena.Center(2));
            logf(['the center is: [' num2str(xa{n}) ' ' num2str(ya{n}) ']']);
            arena_r{n} = floor(arena.Radius);
            logf(['the radius is: ' num2str(arena_r{n}) ]);
            arenaD{n} = [xa{n} ya{n} arena_r{n}];

        end

        f = waitbar(0,'1','Name','Tracking mouse...');
                        
        close all;
        
        result.filepath = filename;
        result.frameRate = frameRate;
        result.duration = duration;
        result.nframes = nframes;
        result.videosize = videosize;
        result.StartFrame = StartFrame;
        result.LastFrame = LastFrame;
        result.arenaD = arenaD;
        result.mouseN = MouseN;
        result.step = Step;
        
        result.positions = [];        

     tic;
     tStart = tic;
     x = StartFrame:Step:floor((LastFrame-StartFrame)/Step)*Step+StartFrame; % floor((LastFrame-StartFrame)/Step)*Step+StartFrame
     
     c1 = cell(MouseN, length(x), 2);
     m1_majl = cell(MouseN, length(x), 1);
     m1_minl = cell(MouseN, length(x), 1);
     m1_ori = cell(MouseN, length(x), 1);
     m1_ecc = cell(MouseN, length(x), 1);
     m1_area = cell(MouseN, length(x), 1);
     
     for jj = 1:length(x)
         
            frame2 = bwareaopen(im2bw(255-read(video_obj,x(jj)), threshold),objpixels) * 255; % remove any object smaller than 50 pixels
            for m = 1:MouseN 
                %mouse = regionprops(frame2,'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
                frame3 = frame2(ya{m}-arena_r{m}:ya{m}+arena_r{m}, xa{m}-arena_r{m}:xa{m}+arena_r{m});
                [X,Y]=ndgrid(1:size(frame3,1),1:size(frame3,2));
                X = X - arena_r{m};
                Y = Y - arena_r{m};
                L = sqrt(X.^2 + Y.^2)>arena_r{m};
                maskedImage = frame3;
                maskedImage(L) = 0; 
                %imshow(maskedImage);
                
                mouse = regionprops(maskedImage, 'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
               
                if ~isempty([mouse.Area])                    
                    areaArray = [mouse.Area];
                    [~,idx] = max(areaArray);
                    c1{m}(jj,:) = floor(mouse(idx).Centroid +[xa{m}-arena_r{m} ya{m}-arena_r{m}]);
                    
                    m1_majl{m}(jj,:) = mouse(idx).MajorAxisLength;
                    m1_minl{m}(jj,:) = mouse(idx).MinorAxisLength;
                    m1_ori{m}(jj,:) = mouse(idx).Orientation;
                    m1_ecc{m}(jj,:) = mouse(idx).Eccentricity;
                    m1_area{m}(jj,:) = mouse(idx).Area;                    
                else
                    %the following is needed because mice sometimes jump
                    %up at the edge (out of the ROI).
                    c1{m}(jj,:) = c1{m}(jj-1,:);
                    m1_majl{m}(jj,:) = m1_majl{m}(jj-1);
                    m1_minl{m}(jj,:) = m1_minl{m}(jj-1);
                    m1_ori{m}(jj,:) = m1_ori{m}(jj-1);
                    m1_ecc{m}(jj,:) = m1_ecc{m}(jj-1);
                    m1_area{m}(jj,:) = m1_area{m}(jj-1);
                 end
            end
            %the following three lines greatly accelerate the read
            %videoframe process
            if Step>1 && x(jj)+1<LastFrame
                read(video_obj,[x(jj)+1 min(LastFrame,x(jj)+Step-1)]);            
            end

         waitbar(jj/length(x),f,sprintf(['Tracking Progress: ' num2str(round(jj*10000/length(x))/100) '%%']))
     end      
     delete(f)   
     result.positions = {};
     result.area = {};
     result.orientation = {};
     
     for p = 1:MouseN         
         result.positions = [result.positions c1{p}];
         result.area = [result.area m1_area{p}];
         result.orientation = [result.orientation [m1_majl{p} m1_minl{p} m1_ori{p} m1_ecc{p}]];     
     end
     
     save(resultpath, 'result');
     tElapsed = toc(tStart);   
     logf(['tracking completed and data saved! used: ' num2str(tElapsed) 'seconds']);          
end
               
    if exist(resultpath, 'file') == 2 && exist([dirpath video_name '/' [video_name '_result.xls']], 'file') ~= 2
        logf('Tracking has been done, running analysis...');
        analysis(filename, dirpath);
    end
    
        logf('Tracking and analysis have been done');

         prompt = sprintf('Tracking and anlysis completed, check the tracking results?');
         button = questdlg(prompt,'Check the results','Yes','No','Yes');

         trackmovie_name = [video_name '_result.avi'];
         trackmoviepath = [dirpath video_name '/' trackmovie_name];
         
         if strcmp(button,'Yes')
            if exist(trackmoviepath, 'file') ~= 2
                close all;
                video_obj = VideoReader(filepath);
                R = load(resultpath);
                MouseN = R.result.mouseN;
                nframes = R.result.nframes;
                vidSize = R.result.videosize;
                startframe = R.result.StartFrame;
                lastframe = R.result.LastFrame;
                c1 = R.result.positions;
                ori = R.result.orientation;
                area = R.result.area;
                step = R.result.step;
                x = startframe:step:floor((lastframe-startframe)/step)*step+startframe;
                
                F(min(length(x),900)) = struct('cdata',[],'colormap',[]);

                colors = ['b', 'r', 'g', 'm', 'y'];

                for i=1:min(length(x),900) %here you can change to for i=1:length(x)

                    vframe = read(video_obj,x(i));
                    
                    %image(vframe);
                    %figure('Position', [10 10 vidSize(1) vidSize(2)]) 
                    image(vframe);
                    %set(gcf, [10, 10, vidSize(1), vidSize(2)]);
                    truesize([vidSize(2)./2 vidSize(1)./2]); 
                    hold on

                    for mn=1:MouseN
                        m1_majl{mn} = ori{mn}(:,1);
                        m1_minl{mn} = ori{mn}(:,2);
                        m1_ori{mn} = ori{mn}(:,3);
                        m1_area{mn} = area{mn}(:,1);
                        m1_ecc{mn} = ori{mn}(:,4);
                        p1 = calculateEllipse(c1{mn}(i,1),c1{mn}(i,2),m1_majl{mn}(i,1)./2,m1_minl{mn}(i,1)./2,m1_ori{mn}(i,1));
                        plot(p1(:,1), p1(:,2), ['.-' colors(mn)])
                        text(c1{mn}(i,1),c1{mn}(i,2),num2str(mn),'Color',colors(mn),'FontSize',20,'FontWeight','bold')
    %                     text(40+270*(mn-1),400,['mouse area: ' num2str(m1_area{mn}(i,1)) ' ecc: ' num2str(m1_ecc{mn}(i,1))],'Color','r','FontSize',12,'FontWeight','bold')
                    end

                    title(['Frame ' num2str(x(i))])
                    hold off
                    F(i) = getframe;
                end

                v = VideoWriter(trackmoviepath);
                open(v);
                writeVideo(v,F);
                close(v);
                logf('Video saved.');
            else
                h=implay(trackmoviepath);
                play(h.DataSource.Controls);
            end
        end
        
end        

function analysis(filename1, dirpath1)
   
    filepath = [dirpath1 filename1];
    
    [u1 video_name u2] = fileparts(filepath);
                
    result_name = [video_name '_result.mat'];
    resultpath = [dirpath1 video_name '/' result_name];

    R = load(resultpath);
    startframe = R.result.StartFrame;
    lastframe = R.result.LastFrame;
    step = R.result.step;
    totalframe = floor((lastframe-startframe)/step)+1;
    frameRate = R.result.frameRate;
    videosize = R.result.videosize;
    arena = R.result.arenaD;
    pos = R.result.positions;
    ori = R.result.orientation;
    area = R.result.area;
    mouseN = R.result.mouseN;
    colors = ['b', 'r', 'g', 'm'];

    for i=1:mouseN
        m1_majl{i} = ori{i}(:,1);
        m1_ecc{i} = ori{i}(:,4);
        m_area{i} = area{i}(:,1);
    end
% subplot(3,2,1)
% hold on
% plot(m1_ecc,area(:,1),'r',m2_ecc,area(:,2),'b');
% xlabel('Eccentricity','FontSize',16);
% ylabel('Mouse Area','FontSize',16);
% set(gca,'fontsize',16,'linewidth',2,'box','on')
% hold off
% 
% subplot(3,2,2)
% hold on
% plot(m1_ecc,m1_majl,'r',m2_ecc,m2_majl,'b');
% xlabel('Eccentricity','FontSize',16);
% ylabel('Major Axis Length','FontSize',16);
% set(gca,'fontsize',16,'linewidth',2,'box','on')
% hold off

    subplot(3,1,1)
    for i=1:mouseN
        ecc1{i} = histfit(m1_ecc{i},20,'gamma');
        ecc1{i}(1).FaceColor = [i./mouseN 0.8 1./i];
        ecc1{i}(2).Color = [1./i 0.2 i./mouseN];
        set(ecc1{i}(1),'facealpha',0.5)
        hold on
    end    
        xlabel('Eccentricity','FontSize',16);
        ylabel('Distribution','FontSize',16);
        yt = get(gca, 'YTick');
        set(gca, 'YTick', yt, 'YTickLabel', round(100*yt/numel(m1_ecc{1}))/100)
        set(gca,'fontsize',16,'linewidth',2,'box','on')
        hold off

    subplot(3,1,2)
    for i=1:mouseN
        h1{i} = histfit(area{i}(:,1),20,'Normal');
        h1{i}(1).FaceColor = [i./mouseN 0.8 1./i];
        h1{i}(2).Color = [1./i 0.2 i./mouseN];
        set(h1{i}(1),'facealpha',0.5)
        %n1{i} = sum(area{i}(:,1)<600);
        %text(200,400-i*50,['mouse' num2str(i) ' standing up: ' num2str(n1{i})]);    
        hold on
    end
    xlabel('Mouse Area','FontSize',16);
    ylabel('Distribution','FontSize',16);
    set(gca,'fontsize',16,'linewidth',2,'box','on')
    hold off

    subplot(3,1,3)
    for i=1:mouseN
        majl1{i} = histfit(m1_majl{i},40,'Normal');
        majl1{i}(1).FaceColor = [i./mouseN 0.8 1./i];
        majl1{i}(2).Color = [1./i 0.2 i./mouseN];        
        set(majl1{i}(1),'facealpha',0.5)
        hold on
    end
    xlabel('Major Axis Length','FontSize',16);
    ylabel('Distribution','FontSize',16);
    set(gca,'fontsize',16,'linewidth',2,'box','on')

    imgfilename1 = [dirpath1 video_name '/' video_name '_orientation'];

    % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    set(gcf, 'Position', [0, 0.04, 450, 900]);
    print(gcf,[imgfilename1 '.tif'],'-dtiff','-r300');

    hold off

    x = startframe:step:lastframe;
    ctimearray=zeros(length(x),1);
    for j=1:length(x)
        ctimearray(j) = (j-1).*step./frameRate;
    end

    xlsfilename = [video_name '_result.xls'];
    logf('Saving Results into Excel and TIFF files, please wait...', [dirpath1 video_name '/' xlsfilename]);

    for i=1:mouseN

        arenaR = arena{i}(:,3);
        pixels = 6.*25.4./arenaR; %pixel size in mm; the arena is 6 inche in radius
        peripheral = floor(arenaR./6.*2); % 2 inches; Check the mouse length distribution (major axis length), at least half of that can be set as peripheral. 
        x0 = arena{i}(:,1);
        y0 = arena{i}(:,2);
        r0 = arenaR - peripheral;
        th = 0:pi/1000:2*pi;
        
        xa = arenaR * cos(th) + x0;
        ya = arenaR * sin(th) + y0;
        xv = r0 * cos(th) + x0;
        yv = r0 * sin(th) + y0;
        
        mouse = ['mouse' num2str(i)]; 
        xq1 = pos{i}(:,1);
        yq1 = pos{i}(:,2);
        in1 = inpolygon(xq1,yq1,xv, yv);  
        timesin1 = numel(xq1(in1));
        pathl1 =[0];
        dd1 = 0;
        distance1 = [0];

        for k = 2:length(pos{i})
              D1 = sqrt((pos{i}(k,1)-pos{i}(k-1,1))^2+(pos{i}(k,2)-pos{i}(k-1,2))^2);
              D1 = D1.* pixels; % path length in mm
              pathl1 = [pathl1; D1];
              dd1 = dd1 + D1./1000; % travel distance in m
              distance1 = [distance1; dd1];
        end

        titlerow = horzcat({'Time (s)'},{[mouse '_X']},{[mouse '_Y']},{[mouse '_PathL (mm)']},{[mouse '_Distance(m)']});
        finaldata = horzcat(ctimearray,xq1,yq1,pathl1,distance1);

        TimeInner1 = timesin1.*step./ frameRate;
        TimeOuter1 = totalframe.*step./frameRate - TimeInner1;
        Thigmotaxis1 = 1-timesin1./length(x);
        summarytitle = [{'FileName'},{[mouse '_Outer (s)']},{[mouse '_Inner (s)']},{[mouse '_Thigmotaxis']}];
        summaryresults = [{[dirpath1 video_name '.mov']},TimeOuter1,TimeInner1,Thigmotaxis1];

        xlwrite([dirpath1 video_name '/' xlsfilename],titlerow,mouse,'A1');
        xlwrite([dirpath1 video_name '/' xlsfilename],finaldata,mouse,'A2');
        xlwrite([dirpath1 video_name '/' xlsfilename],summarytitle,mouse,'G1');
        xlwrite([dirpath1 video_name '/' xlsfilename],summaryresults,mouse,'G2');

        xv1{i} = xv;
        yv1{i} = yv;
        xa1{i} = xa;
        ya1{i} = ya;
        finaldata1{i} = finaldata;

    end


    hsum=figure('Visible','on'); 
    Bkg = 230 * ones(videosize(2), videosize(1), 3, 'uint8');

    subplot(1,2,1); image(Bkg); axis image
    %freezeColors;
    title([video_name '.mov'],'Interpreter','none');
    hold all
    %traveldistanceplot = [];
    for i=1:mouseN
        plot(pos{i}(:,1),pos{i}(:,2),'Color',colors(i));
        plot(xv1{i},yv1{i},'--k','LineWidth',2); % plot inner arena for each mouse
        plot(xa1{i},ya1{i},'k','LineWidth',2); % plot the entire arena
        %traveldistanceplot = [traveldistanceplot finaldata1{i}(:,5)];
    end
    set(gca,'fontsize',20)
    ymax = 0;

    subplot(1,2,2);
    for i=1:mouseN 
        ymax = max(ymax, max(finaldata1{i}(:,5)));
        plot(finaldata1{1}(:,1),finaldata1{i}(:,5),'Color',colors(i),'LineWidth',2);
        hold on
    end
    hold off
%     plot(finaldata1{1}(:,1),traveldistanceplot,'LineWidth',2);
    xlim([0 length(x).*step./frameRate]);
    ylim([0 ymax]);
    xlabel('Times (s)');
    ylabel('Travel Distance (m)');
    set(gca,'linewidth',2,'fontsize',20,'box', 'off');
    set(gcf,'Units','Normalized','Position',[0 0 1 0.5],'PaperPositionMode','auto','PaperSize',[14 14]);
    % title(sprintf(video_name));
    title([video_name '.mov'],'Interpreter','none');

    imgfilename = [dirpath1 video_name '/' video_name '_summary'];
    print(gcf,[imgfilename '.tif'],'-dtiff','-r300');
    hold off
    %delete(hsum);

    logf('Results saved to folder %s', [dirpath1 video_name '/']);

end    

function pos = customWait(hROI)

    % Listen for mouse clicks on the ROI
    l = addlistener(hROI,'ROIClicked',@clickCallback);

    % Block program execution
    uiwait;

    % Remove listener
    delete(l);

    % Return the current position
    pos = hROI.Position;

end

function clickCallback(~,evt)

    if strcmp(evt.SelectionType,'double')
        uiresume;
    end

end

function  th=clicksubplot  
    w = waitforbuttonpress;
    if w == 0
        th = get(gca,'tag');
    end
end    

function logf(varargin)
    message = sprintf(varargin{1}, varargin{2:end});
    str = ['[' datestr(now(), 'HH:MM:SS') '] ' message];
    disp(str);
end

function [X,Y]=calculateEllipse(x, y, a, b, angle, steps)
%     %# This functions returns points to draw an ellipse
%     %#  @param x     X coordinate
%     %#  @param y     Y coordinate
%     %#  @param a     Semimajor axis
%     %#  @param b     Semiminor axis
%     %#  @param angle Angle of the ellipse (in degrees)

    narginchk(5, 6);
    if nargin<6, steps = 36; end

    beta = -angle * (pi / 180);
    sinbeta = sin(beta);
    cosbeta = cos(beta);

    alpha = linspace(0, 360, steps)' .* (pi / 180);
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);

    X = x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta);
    Y = y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta);

    if nargout==1, X = [X Y]; end

end