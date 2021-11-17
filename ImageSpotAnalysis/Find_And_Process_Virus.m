function [BindingDataToSave, OtherDataToSave] =...
    Find_And_Process_Virus(StackFilePath,ThresholdInput,StackFilename, ...
    StackNum, DefaultPathname,Options)
   

    %The first image in the stack is read and displayed. A 3D array
    %(ImageStackMatrix) is created which will contain the data for all images in
    %the stack.  ImageStackMatrix is pre-allocated with zeros to make the next
    %for loop faster.
    StackInfo = imfinfo(StackFilePath);
    if isnan(Options.FrameNumberLimit)
        NumFrames = length(StackInfo);
    else
        NumFrames = Options.FrameNumberLimit;
    end
        ImageWidth = StackInfo.Width; %in pixels % PUT DEBUG STOP HERE!!!!!!!!
        ImageHeight = StackInfo.Height;
        BitDepth = StackInfo.BitDepth;
      ImageStackMatrix = zeros(ImageWidth, ImageHeight, NumFrames, 'uint16');
        BWStackMatrix = ImageStackMatrix > 0; %Create a logical stack the same size as the image stack.
        
        %Preallocate threshold vectors as well
        ThresholdToFindViruses = zeros(NumFrames,1);
        RoughBackground = zeros(NumFrames,1);
        %Back_MedMed = zeros(NumFrames,1);
        
        %Set up figures
        [FigureHandles] = Setup_Figures(Options);
        
    %This for loop populates the ImageStackMatrix with the data from each image
    %in the stack.  The 1st two dimensions are the x,y of the image plane and the 3rd 
    %dimension is the frame number.
    for b = 1:NumFrames
        CurrentFrameImage = imread(StackFilePath,b);
        if  strcmp(Options.DualColor,'y')
            if  strcmp(Options.ReverseFind,'y')
                 if rem(b,2) == 0
                     ImageStackMatrix(:,:,b-1) = CurrentFrameImage;
                 else 
                     ImageStackMatrix(:,:,b+1) = CurrentFrameImage;
                 end

            else
                ImageStackMatrix(:,:,b) = CurrentFrameImage;
            end
        else
            ImageStackMatrix(:,:,b) = CurrentFrameImage;
        end 
    end
    
    for b = 1:NumFrames
        CurrentFrameImage = ImageStackMatrix(:,:,b);
        
        if b == 1
%            Display first image
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            imshow(ImageStackMatrix(:,:,1), [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
            drawnow
        end
        
        %We define the moving threshold which will be used to find Viruses
        %throughout the stream.
        RoughBackground(b) = mean(median(CurrentFrameImage));
        ThresholdToFindViruses(b) = (RoughBackground(b) + ThresholdInput)/2^BitDepth;
        
        %We apply the threshold to create a big logical matrix
        CurrThresh = ThresholdToFindViruses(b);
%         BWStackMatrix(:,:,b) =  CurrentFrameImage >(RoughBackground(b) + ThresholdInput);
        BWStackMatrix(:,:,b) = im2bw(CurrentFrameImage, CurrThresh);
        BWStackMatrix(:,:,b) = bwareaopen(BWStackMatrix(:,:,b), Options.MinParticleSize, 8);
        
        if rem(b,20)==0
            set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow);
            title(strcat('Loading Frame :', num2str(b),'/', num2str(NumFrames)));
            drawnow
        end
        
    end

    
    
    set(0,'CurrentFigure', FigureHandles.BackgroundTraceWindow);
    hold on
    plot(RoughBackground,'r-')
    title('Background intensity versus frame number');
    
    NumFramesAnalyzed = 0;
    
    %Now we find all the viruses In each image
    if  strcmp(Options.DualColor,'y')
        FrameNumberList = 1:2:NumFrames;
    else
        FrameNumberList = 1:NumFrames;
    end
    
    for CurrFrameNum = FrameNumberList
        NumFramesAnalyzed = NumFramesAnalyzed + 1;
        NumberColocalizedGood = 0;
        NumberColocalizedTotal = 0;
        
        if strcmp(Options.DualColor,'y')
            CurrentImage2 = ImageStackMatrix(:,:,CurrFrameNum+1);
            CurrentImage = ImageStackMatrix(:,:,CurrFrameNum);
            BinaryCurrentImage = BWStackMatrix(:,:,CurrFrameNum);
            
            set(0,'CurrentFigure',FigureHandles.Image2Window);
            hold off
            imshow(CurrentImage2, [Options.MinImage2Show, Options.MaxImage2Show], 'InitialMagnification', 'fit','Border','tight');
            title('Second Color Image')
            hold on
        else
            CurrentImage = ImageStackMatrix(:,:,CurrFrameNum);
            BinaryCurrentImage = BWStackMatrix(:,:,CurrFrameNum);
        end
        
        if strcmp(Options.IgnoreAreaNearBigParticles,'y')
                [CurrentImage,BinaryCurrentImage] = Remove_Area_Around_Big_Particles(Options,...
                    ImageWidth,ImageHeight,CurrentImage,BinaryCurrentImage,FigureHandles);
        end
        
        if strcmp(Options.RemoveBleedthroughAreas,'y')
            [CurrentImage,BinaryCurrentImage] = Remove_Bleedthrough_Areas(CurrFrameNum,CurrentImage,BinaryCurrentImage,...
                ImageStackMatrix,Options,FigureHandles,ImageHeight,ImageWidth);
        end
            NonzeroPixels = CurrentImage(CurrentImage ~= 0);
            NumberPixelsNotBlackedOut = length(NonzeroPixels);
            NumberMicronsNotBlackedOut = NumberPixelsNotBlackedOut*(0.16)^2;
        
        %All of the isolated regions left behind are "virus regions" and will
        %be analyzed.
            VirusComponentArray = bwconncomp(BinaryCurrentImage,8);

        %The properties associated with each virus in the binary image are
        %extracted.
            VirusProperties = regionprops(VirusComponentArray, CurrentImage, 'Centroid',...
                'Eccentricity', 'PixelValues', 'Area','PixelIdxList');
            NumberOfVirusesFound = length(VirusProperties);
            NumberGoodViruses = 0;
            NumberBadViruses = 0;
            
        %Plot the image
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            hold off
            imshow(CurrentImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
            hold on
            
            set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
            imshow(BinaryCurrentImage, 'InitialMagnification', 'fit','Border','tight');
            drawnow
        
        %Analyze each region
        % FIX THIS
        for n = 1:NumberOfVirusesFound
            CurrentVirusProperties = VirusProperties(n);
            CurrVesX = round(VirusProperties(n).Centroid(1)); 
            CurrVesY = round(VirusProperties(n).Centroid(2));
                CurrentVirusProperties.Centroid = [CurrVesX, CurrVesY];
                
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            hold off
            title(strcat('Virus :', num2str(n),'/', num2str(NumberOfVirusesFound)));
            drawnow
            hold on
            
            %Apply many tests to see if Virus is good
            [IsVirusGood, CroppedImageProps, CurrentVirusBox, OffsetFrom2DFit, Noise,...
                DidFitWork, CroppedVesImageThresholded, SizeOfSquareAroundCurrVesicle,...
                CurrVesicleEccentricity, NewArea, ReasonVirusFailed] =...
                Simplified_Test_Goodness(CurrentImage,CurrentVirusProperties,BitDepth,...
                CurrThresh, Options.MinParticleSize, Options.MaxEccentricity,ImageWidth, ImageHeight,Options.MaxParticleSize,BinaryCurrentImage);
            
            if strcmp(IsVirusGood,'y')
                LineColor = 'g-';
                NumberGoodViruses = NumberGoodViruses + 1;
                
            elseif strcmp(IsVirusGood,'n')
                LineColor = 'r-';
                NumberBadViruses = NumberBadViruses + 1;
                if strcmp(Options.DispReasonFailed,'y')
                    disp(ReasonVirusFailed)    
                end
            end
                                               
            %Plot a box around the Virus
                CVB = CurrentVirusBox;
                BoxToPlot = [CVB.Bottom,CVB.Left;CVB.Bottom,CVB.Right;CVB.Top,CVB.Right;CVB.Top,CVB.Left;CVB.Bottom,CVB.Left];

                set(0,'CurrentFigure',FigureHandles.ImageWindow);
                plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor)
                hold on
                drawnow
            
            %Now we grab the intensity of the current virus particle

                CurrentVirusArea = ImageStackMatrix(...
                    CurrentVirusBox.Top:CurrentVirusBox.Bottom,...
                    CurrentVirusBox.Left:CurrentVirusBox.Right,...
                    CurrFrameNum);

                CurrentRawIntensity = sum(sum((CurrentVirusArea)));

                CurrentIntensityBackSub = CurrentRawIntensity -...
                    RoughBackground(CurrFrameNum).*(CurrentVirusBox.Bottom - CurrentVirusBox.Top + 1)^2;
            
                % Now we calculate the intensity in the same region of 
                % interest in the second color image.
                if strcmp(Options.DualColor,'y')
                    if strcmp(Options.SimplifyDualColor, 'n')
                        [IsVirusGood,LineColor,RoughIntensity2,GaussianIntensity2] = Calculate_Intensity_Second_Color(...
                            CurrentImage2,CurrentVirusBox,RoughBackground,...
                            FigureHandles,BoxToPlot,LineColor,CurrFrameNum,IsVirusGood);

                        if GaussianIntensity2 > Options.GaussianColocCutoff
                                if strcmp(Options.ShowColoc,'y')
                                    LineColor2 = 'y-';
                                end
                                NumberColocalizedTotal = NumberColocalizedTotal +1;
                                if strcmp(IsVirusGood,'y')
                                    NumberColocalizedGood = NumberColocalizedGood +1;
                                end
                        else
                            if strcmp(Options.ShowColoc,'y')
                                LineColor2 = 'm-';
                            end
                        end

                        if strcmp(Options.ShowColoc,'y')
                            
                            if strcmp(Options.ReverseFind,'y')
                                set(0,'CurrentFigure',FigureHandles.Image2Window);
                            else
                                set(0,'CurrentFigure',FigureHandles.ImageWindow);
                            end
                            
                            plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor2)
                            hold on
                            drawnow
                            
                        end
                    else
                        RoughIntensity2 = NaN;
                        GaussianIntensity2 = NaN;
                    end
                end


            %Save the data
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).RawIntensity = CurrentRawIntensity;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).IntensityBackSub = CurrentIntensityBackSub;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).Coordinates = VirusProperties(n).Centroid;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).Area = VirusProperties(n).Area;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).Eccentricity = VirusProperties(n).Eccentricity;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).FullFilePath = StackFilePath;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).StreamFilename = StackFilename;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).BoxAroundVirus = CurrentVirusBox;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).IsVirusGood = IsVirusGood;
                BindingDataToSave(NumFramesAnalyzed).VirusData(n).ReasonVirusFailed = ReasonVirusFailed;
                
                
                if strcmp(Options.DualColor,'y')
                    BindingDataToSave(NumFramesAnalyzed).VirusData(n).RoughIntensity2 = RoughIntensity2;
                    BindingDataToSave(NumFramesAnalyzed).VirusData(n).GaussianIntensity2 = GaussianIntensity2;
                end
        end
        BindingDataToSave(NumFramesAnalyzed).TotalVirusesBound = NumberOfVirusesFound;
        BindingDataToSave(NumFramesAnalyzed).NumberGoodViruses = NumberGoodViruses;
        BindingDataToSave(NumFramesAnalyzed).NumberBadViruses = NumberBadViruses;
        BindingDataToSave(NumFramesAnalyzed).NumberPixelsNotBlackedOut = NumberPixelsNotBlackedOut;
        BindingDataToSave(NumFramesAnalyzed).NumberMicronsNotBlackOut = NumberMicronsNotBlackedOut;
        
        
        if strcmp(Options.DualColor,'y')
            if strcmp(Options.SimplifyDualColor, 'n')
                BindingDataToSave(NumFramesAnalyzed).NumberColocalizedGood = NumberColocalizedGood;
                BindingDataToSave(NumFramesAnalyzed).NumberColocalizedTotal = NumberColocalizedTotal;
            else
                BindingDataToSave(NumFramesAnalyzed).NumberColocalizedGood = NaN;
                BindingDataToSave(NumFramesAnalyzed).NumberColocalizedTotal = NaN;
            end
        end
        
        BindingDataToPrint = BindingDataToSave(NumFramesAnalyzed)
        
    end
    
    OtherDataToSave.ThresholdsUsed = ThresholdToFindViruses;
    OtherDataToSave.RoughBackground = RoughBackground;
    
end