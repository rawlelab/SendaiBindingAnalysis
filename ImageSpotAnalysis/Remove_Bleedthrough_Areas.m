function [CurrentImage,BinaryCurrentImage] = Remove_Bleedthrough_Areas(CurrFrameNum,CurrentImage,BinaryCurrentImage,...
    ImageStackMatrix,Options,FigureHandles,ImageHeight,ImageWidth)

MinimumAreaSize = 8;
BilayerThreshold = 450;
DisplayFactor = 0.3;
MinimumBlackoutArea = 50;

BilayerImage = ImageStackMatrix(:,:,CurrFrameNum+1);
    BilayerMinimumShow = median(min(BilayerImage));
    BilayerMaximumShow = max(max(BilayerImage)) - DisplayFactor*max(max(BilayerImage));
    BilayerIntensityCutoff = mean(median(BilayerImage)) + BilayerThreshold;
BinaryBilayerImage = BilayerImage >BilayerIntensityCutoff;
BinaryBilayerImage = bwareaopen(BinaryBilayerImage,MinimumAreaSize,8);


%Plot the bilayer image
        set(0,'CurrentFigure',FigureHandles.BilayerWindow);
        hold off
        imshow(BilayerImage, [BilayerMinimumShow, BilayerMaximumShow], 'InitialMagnification', 'fit','Border','tight');
        hold on
        drawnow

        VirusComponentArray = bwconncomp(BinaryBilayerImage,8);

    %The properties associated with each virus in the binary image are
    %extracted.
        VirusProperties = regionprops(VirusComponentArray, BilayerImage, 'Centroid',...
            'Eccentricity', 'PixelValues', 'Area','PixelIdxList');
        NumberOfVirusesFound = length(VirusProperties);

    for j = 1:NumberOfVirusesFound
        if VirusProperties(j).Area > MinimumAreaSize
            CurrentVirusArea = max(VirusProperties(j).Area,MinimumBlackoutArea);
            CurrentVirusCoordinates = VirusProperties(j).Centroid;
                CurrentVirusX = CurrentVirusCoordinates(1);
                CurrentVirusY = CurrentVirusCoordinates(2);

            SizeOfSquareAroundCurrVesicle = (sqrt(CurrentVirusArea)*2);

                %Coordinates are set up to define the area around the vesicle
                %of interest (i.e. the ROI)  
                CurrentBoxLeft = max([round(CurrentVirusX) - round(SizeOfSquareAroundCurrVesicle/2),...
                    1]);
                CurrentBoxRight = min([round(CurrentVirusX) + round(SizeOfSquareAroundCurrVesicle/2),...
                    ImageWidth]);
                CurrentBoxTop = max([round(CurrentVirusY) - round(SizeOfSquareAroundCurrVesicle/2),...
                    1]);
                CurrentBoxBottom = min([round(CurrentVirusY) + round(SizeOfSquareAroundCurrVesicle/2),...
                    ImageHeight]);

 
                BilayerImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;

                BinaryBilayerImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;
                
                
                CurrentImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;

                BinaryCurrentImage(CurrentBoxTop:CurrentBoxBottom,...
                    CurrentBoxLeft:CurrentBoxRight) = 0;
        end

    end

   %Plot the corrected image
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            hold off
            imshow(CurrentImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
            hold on
            
            set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
            imshow(BinaryCurrentImage, 'InitialMagnification', 'fit','Border','tight');
            drawnow
        
    %Plot the corrected image
            set(0,'CurrentFigure',FigureHandles.BilayerWindow);
            hold off
            imshow(BilayerImage, [BilayerMinimumShow, BilayerMaximumShow], 'InitialMagnification', 'fit','Border','tight');
            hold on
        
end