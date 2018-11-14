

function MyoCount(varargin)


inputs.SmallestMyotubePixelCount = 500;
inputs.SmallestNucleusPixelCount = 300;
inputs.MyotubeChannel = 2;
inputs.NucChannel = 3;
inputs.FillSize = 10;
inputs.NucFillSize = 5;
inputs.MinCircleRad = 15;
inputs.MaxCircleRad = 40;
inputs.TubeThresh = 1.0;
inputs.MinNuclei = 3;
inputs.Output = "";
inputs.FilePath = "";

known_vars = fieldnames(inputs);

numerics = ["SmallestMyotubePixelCount","SmallestNucleusPixelCount","MyotubeChannel","NucChannel","FillSize","NucFillSize","MinCircleRad","MaxCircleRad","TubeThresh","MinNuclei"];

if mod(nargin,2) ~= 0
    error('Expecting an even number of arguments');
end


for idx = 1 : 2 : nargin-1
    if ~ismember(varargin{idx}, known_vars)
        error('Argument "%s" is not a known parameter name', varargin{idx});
    end
    if ismember(varargin{idx},numerics)
        val = str2double(varargin{idx+1});
        if isnan(val)
            error('Value "%s" for variable "%s" is not a numeric scalar', varargin{idx+1}, varargin{idx});
        end
        inputs.(varargin{idx}) = val;
    else
        inputs.(varargin{idx}) = varargin{idx+1};
    end
end

MyocountVersion = 'Myocount beta version 1.2';
disp(MyocountVersion);

%if no input image specified prompt to select one.
if(strcmp(inputs.FilePath,"")~=0)
    [FileName,PathName] = uigetfile('*.tif','Please select tif');
    inputs.FilePath = fullfile(PathName, FileName);
end

%if there's a wildcard in the path given, find all files matching it and
%run tool against all sequentially.
fileList = dir(inputs.FilePath);
filecount = size(fileList);
disp(fileList);

csvdata = [];
for i=1:filecount
    if any(strcmp(fileList(i).name, {'.', '..'}))
        continue
    end
         
    File = fullfile(fileList(i).folder, fileList(i).name);
    disp(File);
    
    [a,b,c,d] = processimagefile( ...
        File, ...
        inputs.SmallestMyotubePixelCount, ...
        inputs.SmallestNucleusPixelCount, ...
        inputs.MyotubeChannel, ...
        inputs.NucChannel, ...
        inputs.FillSize, ...
        inputs.NucFillSize, ...
        inputs.MinCircleRad, ...
        inputs.MaxCircleRad, ...
        inputs.TubeThresh, ...
        inputs.MinNuclei, ...
        inputs.Output);
    
    name = fileList(i).name;
    if(name=="")
        error("ERROR:Empty file name returned, this may be a bug");
        exit(0);
    end
    if(a=="")
        error("ERROR:Invalid coverage returned, this may be a bug");
        exit(0);
    end
    if(b=="")
        error("ERROR:Invalid count returned, this may be a bug");
        exit(0);
    end
    if(c=="")
        error("ERROR:Invalid fused count returned, this may be a bug");
        exit(0);
    end
    if(d=="")
        error("ERROR:Invalid final count returned, this may be a bug");
        exit(0);
    end
    
    

    x = {name,a,b,c,d};
    csvdata = [csvdata; x ];
end

if(strcmp(inputs.Output,"")==0)
    
    disp("Generating CSV file");
    OutputCSV = fullfile(inputs.Output, 'Summary.csv');
    
    fid = fopen(OutputCSV,'wt');
    if fid>0
        for k=1:size(csvdata,1)
            fprintf(fid,'%s,%s,%s,%s,%s\n',csvdata{k,:});
        end
        fclose(fid);
    end
    disp("Completed generating CSV summary file");
    exit(0);
end


end



function [a] = validateoutputfolder(Outputfolder,name)

Output = "";
Outputfolder = strtrim(Outputfolder);
if(  strcmp( Outputfolder(end),'\')==0 && strcmp( Outputfolder(end),'/')==0 )
    Outputfolder = strcat(Outputfolder,'/');
end


if(strcmp(Outputfolder ,'')==0)
    if (7==exist(fullfile(Outputfolder),'dir'))
        disp(Outputfolder);
        Output = fullfile(Outputfolder, name);
    else
        disp(strcat(Outputfolder,' is not a valid existing folder'));
        disp(Outputfolder);
        disp(exist(Outputfolder,'dir'));
        exit(1)
    end
    
end
a = Output;
end


function [a,b,c,d] = processimagefile(imagefilelocation, ...
    SmallestMyotubePixelCount, ...
    SmallestNucleusPixelCount, ...
    MyotubeChannel, ...
    NucChannel, ...
    FillSize, ...
    NucFillSize, ...
    MinCircleRad, ...
    MaxCircleRad, ...
    TubeThresh, ...
    MinNuclei, ...
    Outputfolder)

[~,name,~] = fileparts(imagefilelocation);

Output = validateoutputfolder(Outputfolder,name);

img = imread(imagefilelocation); % Read image
imshow(img);

tubes = getmyotubes(img,MyotubeChannel,SmallestMyotubePixelCount,FillSize,TubeThresh);
nuclei = getnuclei(img,NucChannel,SmallestNucleusPixelCount,NucFillSize);

[x,y]=size(img(:,:,MyotubeChannel));
pixelcount=x*y;


[~,~,tubesCount]=size(tubes);
[x,y,nucleiCount]=size(nuclei);

bw = zeros(x, y);
bw2 = zeros(x, y);
bw3 = zeros(x, y);
bw4 = zeros(x, y);
bw5 = zeros(x, y);
bw6 = zeros(x, y);

for i=1:tubesCount
    %imshow( imoverlay(img, tubes(:,:,i) , [.0 1 0]))
    bw = bw+tubes(:,:,i);
end

for i=1:nucleiCount
    %imshow( imoverlay(img, nuclei(:,:,i) , [.0 0 1]))
    bw2 = bw2+nuclei(:,:,i);
end



figure1 = imshow( imoverlay( imoverlay(img, bw , [.0 1 .0]),bw2 , [.0 0 1]));
if(strcmp(Output,"")==0)
    saveas(figure1,strcat(Output,'.file1.png'));
end

mastercenters = [];
masterradii = [];
[columnsInImage, rowsInImage] = meshgrid(1:y, 1:x);

for i=1:nucleiCount
    [centres, radii, ~] = imfindcircles(nuclei(:,:,i), [MinCircleRad, MaxCircleRad]);
    mastercenters = [mastercenters;centres];
    masterradii = [masterradii;radii];
end

figure1 = imshow(img);
hold on;
plot(mastercenters(:,1), mastercenters(:,2), 'r*'); %// Plot centres
viscircles(mastercenters, masterradii, 'EdgeColor', 'b'); %// Plot circles - Make edge blue
hold off;

if(strcmp(Output,"")==0)
    saveas(figure1,strcat(Output ,'.file2.png'))
end

circlePixels = zeros(x, y);

n = size(mastercenters,1);

for i = 1:n
    centerX = mastercenters(i,1);
    centerY = mastercenters(i,2);
    %cut out detected circles and a margin
    radius = masterradii(i)*1.3;
    circlePixels = circlePixels | (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
end

%remove the detected circles.
bw3 = bw2-circlePixels;

%further split connected shapes.
bw4 = splitconnected(bw3);

%remove anything too small and get the centroids
[indeximgarray,count] = bwlabel(bwareafilt(imbinarize(bw4,0.5),[SmallestNucleusPixelCount , pixelcount/100]));
points = mastercenters;
for i=1:count
    objseparate(:,:,i) = indeximgarray == i;
    s = regionprops(objseparate(:,:,i),'centroid');
    centroid = cat(1, s.Centroid);
    points = [centroid;points];
end

%overlay
figure1 = imshow(img(:,:,NucChannel));
hold on
plot(points(:,1),points(:,2), 'b*')
hold off

if(strcmp(Output,"")==0)
    saveas(figure1,strcat(Output , '.file3.png'))
end


n = size(points,1);
for i = 1:n
    bw5(floor(points(i,2)),floor(points(i,1)))=1;
end

fused = 0;
fused2 = 0;

for i=1:tubesCount
    %imshow( imoverlay(img, tubes(:,:,i) , [.0 1 0]))
    tube = tubes(:,:,i);
    overlappingnuclicount = sum(sum((tube & bw5)==1));
    if overlappingnuclicount >= MinNuclei
        bw6 = bw6+tubes(:,:,i);
        fused2 = fused2 + overlappingnuclicount;
    end
    fused = fused + overlappingnuclicount;
end

nucTubePixelcount = sum(sum(bw6));




percent = (nucTubePixelcount/pixelcount)*100;

figure1 = imshow( imoverlay(img, bw6, [.0 1 .0]));
hold on
plot(points(:,1),points(:,2), 'b*')
legend ('Location','south','Orientation','horizontal',num2str(percent))
hold off
PercentOccupied = strcat('Percent of image occupied by myotubes containing at least ',num2str(MinNuclei),' nuclei:',num2str(percent));
NucleiCount = strcat('Total nuclei Count:',num2str(size(points,1)));
InternalNuclei = strcat('Total nuclei within myotubes:',num2str(fused));
FinalNuclei = strcat('Total nuclei within myotubes containing at least ',num2str(MinNuclei),' nuclei:',num2str(fused2));

disp(PercentOccupied);
disp(NucleiCount);
disp(InternalNuclei);
disp(FinalNuclei);

if(strcmp(Output,"")==0)
    saveas(figure1,strcat(Output , '.file4.png'))
    
    fid=fopen( strcat(Output , '.txt'),'w');
    
    fprintf(fid,  PercentOccupied );
    fprintf(fid,  NucleiCount );
    fprintf(fid,  InternalNuclei );
    fprintf(fid,  FinalNuclei );
    
    fclose(fid);
    
end



a = num2str(percent);
b = num2str(size(points,1));
c = num2str(fused);
d = num2str(fused2);
end




function y = getmyotubes(image,MyotubeChannel,SmallestMyotubePixelCount,FillSize,thresh)
singlecolour = image(:,:,MyotubeChannel);

I_cropped = singlecolour;
I_eq = adapthisteq(I_cropped);

% mark myotubes
marked = imbinarize(I_eq, graythresh(I_eq)*thresh);
filled = imfill(marked,'holes');
opened = imopen(filled, ones(FillSize,FillSize));
%ignore anything smaller than 500 pixels.
openmask = bwareaopen(opened, SmallestMyotubePixelCount);

[x,y]=size(I_eq);
pixels=x*y;


Finalmask = bwareafilt(openmask,[SmallestMyotubePixelCount pixels]);

%Split into parts
[indeximgarray,count] = bwlabel(Finalmask);

for i=1:count
    objseparate(:,:,i) = indeximgarray == i;
end

y = objseparate;
end

function y = getnuclei(image,NucChannel,SmallestNucleusPixelCount,NucFillSize)
%get the correct channel
singlecolour = image(:,:,NucChannel);
I_eq = adapthisteq(singlecolour);
I_eq = wiener2(I_eq,[NucFillSize NucFillSize]);

% use multiple thresholds, extract an image of the scene for each threshold

thresh = multithresh(I_eq,2);
seg_I = imquantize(I_eq,thresh);
map = [1, 0, 0
    0, 1, 0
    0, 0, 1];
RGB = label2rgb(seg_I,map);


singlecolour1 = ~logical(RGB(:,:,1));
singlecolour2 = ~logical(RGB(:,:,2));
singlecolour3 = logical(RGB(:,:,3));

%extract only elements which are small enough to be nuclei

[x,y]=size(I_eq);
pixels=x*y;
%throw out anything larger than 1% of the total image size to discard
%washed out regions
upper = pixels/100;
%throw out anything less than this many pixels
lower = SmallestNucleusPixelCount;

BW2 = bwareafilt(singlecolour1,[lower upper]);
BW3 = bwareafilt(singlecolour2,[lower upper]);
BW4 = bwareafilt(singlecolour3,[lower upper]);

[indeximgarray,count] = bwlabel(BW2 + BW3 + BW4);

for i=1:count
    objseparate(:,:,i) = indeximgarray == i;
end

y = objseparate;
end


function y = splitconnected(shapes)
%credit https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
D = -bwdist(~shapes);
mask = imextendedmin(D,2);
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
test2 = shapes;
test2(Ld2 == 0) = 0;
y = test2;
end

