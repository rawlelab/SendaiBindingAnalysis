function [Options] = Setup_Options()

% NOTE: This current version only counts "good" viruses before checking the
% 2nd color (typically the SLB). But it will still plot on the 2nd image a
% "bad" virus even if the original designation is "good".

    Options.Threshold = 60; %60 150 300
        % This is the number of counts above background which will be used
        % to detect virus particles. You will need to optimize this number 
        % for each set of imaging conditions and/or each data set. An optimal 
        % threshold value has been reached when you are able to see that the
        % program is accurately finding the particles in each image. To avoid
        % bias introduced by the optimization process, you should make sure 
        % that changing the optimal threshold value somewhat doesn't 
        % significantly affect your data. Assuming similar particle 
        % densities/intensities between data sets, and that the same imaging 
        % conditions are used, you shouldn't need to change the threshold 
        % value much if at all between data sets.
        
    Options.DataFileLabel = strcat('Thresh=',num2str(Options.Threshold));
        % This is the label for the output .mat file.
        
    Options.SimplifyDualColor = 'n';
    Options.DispReasonFailed = 'n';
    
    Options.FrameNumberLimit = NaN; 
        % Determines the number of frames which will be loaded and analyzed. 
        % Use 'NaN' to indicate no limit (i.e. all frames will be included).        
        % You would likely only use this option if you wanted to quickly 
        % assess how the program is running without having 
        % to wait for an entire set of images to load. Alternatively, you could 
        % use it to exclude frames at the end of the image stack from your analysis.

    Options.GaussianColocCutoff = 300; %optimize, look at green intensities look by eye
        % For a given region in which a particle was found in the first color 
        % or finding image, this is the number of counts above background in 
        % the second color (background determined by Gaussian fit) above 
        % which the program will consider the region to contain a real 
        % particle in the second color image. This is only for the purposes 
        % of determining the percentage of co-localized particles, calculated 
        % in the Find_And_Process_Virus.m script.
        
    Options.MinImageShow = 350;
    Options.MaxImageShow = 8000; %8000;  
    Options.MinImage2Show = 350;
    Options.MaxImage2Show = 700;
        % These determine the minimum and maximum intensity counts that will 
        % be used to set the contrast for the grayscale images that are displayed.
        % The minimum value will be displayed as black and the maximum value 
        % will be displayed as white. 2 refers to the second color image.

% ---------Parameters Used To Find Particles/Assess Particle 'Goodness'---------
    Options.MinParticleSize = 4;
        % This is the minimum particle size (defined as the number of connected pixels 
        % above the threshold) in order for a particle to be found. 
        % Particles smaller than this size will not be found (i.e. the program 
        % will assume that they are noise).
    Options.MaxParticleSize = 300; %250
        % This is the maximum particle size (defined as the number of connected pixels 
        % above the threshold) in order for a particle to be considered "good". 
        % Particles larger than this size will be designated as "bad".
    Options.MaxEccentricity = 0.8;
        % This is the maximum eccentricity that a particle can have in order to still 
        % be considered "good". 0 = perfect circle, 1 = straight line. If the 
        % eccentricity is too high, that may indicate that the particle being 
        % analyzed is actually two diffraction limited particles close together.
         
    Options.IgnoreAreaNearBigParticles = 'y';
        % 'y' OR 'n'
        % Choose 'y' if you want to ignore the region around bright particles 
        % because the noise nearby from those particles is being incorrectly 
        % identified as other particles. If you choose 'y', then you need 
        % to specify what size is considered "big" below. You will also need 
        % to scale your data to account for the regions which have been ignored, 
        % as this can artificially skew your results lower than they should be.
    Options.MinAreaBig = 300; %250
        % The number of pixels above threshold for a particle to be considered 
        % "too big", in which case the region around that particle will be 
        % ignored if that option is chosen above.
        
        
    Options.RemoveBleedthroughAreas = 'n';
        % 'y' OR 'n'
        % Choose 'y' if you'd like to remove the contribution from bright 
        % debris in the bilayer (un-ruptured vesicles, etc.) which bleeds through 
        % into the viral particle channel and is incorrectly identified as 
        % viral particle. The contribution from these spots should be 
        % minimal, but this is a way to remove them entirely. If you choose 'y', 
        % then your data should be formatted as viral image followed by bilayer image.
        
    Options.DualColor = 'y'; %n?
        % 'y' OR 'n'
        % Choose 'y' to find the particles in the odd-numbered images and 
        % to also grab the intensity in the corresponding region of interest 
        % in the  even-numbered images (i.e. the second color). 
    Options.ReverseFind = 'y'; 
        % 'y' OR 'n'
        % Choose 'y' to switch the order of which images are used to 
        % find (i.e. the even numbered images will be used to find 
        % the particles instead of the odd-numbered images).
        
    Options.ShowColoc = 'n'; %y for coloc
        % 'y' OR 'n'
        % If 'y', noncolocalized spots will show in magenta; colocalized spots will show in yellow.
        % If 'n', good virus will show in green; bad virus will show in red.
        
        if strcmp(Options.ReverseFind,'y')
            OldMax2 = Options.MaxImage2Show;
            OldMin2 = Options.MinImage2Show;
            Options.MaxImage2Show = Options.MaxImageShow;
            Options.MinImage2Show = Options.MinImageShow;
            Options.MaxImageShow = OldMax2;
            Options.MinImageShow = OldMin2;
        end
end