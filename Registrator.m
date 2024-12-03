classdef Registrator < handle
   
% C:\Users\r114m\Desktop\HotzMat\cubes 

    properties
        files
        ct = []
        cst = []
        contours
        cube_for_mask
        bothContours
        structures
        Xshift
        Yshift
        bothStructureSets
        cubeHU
        HUCube
        cubes
        StructuresName
        parentDir1
        FullSubFolderName
        TwoFolders
        parentDir2
        inputFolder
        pathToCubesFolder
        mergedCube
        deformStructures
        maskOfRegImg
        SetOfMasks
        minHU
        maxHU
        meanHUs
        MasksSet
        OtherStructures
        bothMaskSets
    end

    methods

        function obj = Registrator(path)
            %loading files
            fileList = dir(path);
            fileList = fileList(~[fileList.isdir]);
            for i = 1:2
                firstFileName = fileList(i).name;
                fullFilePath = fullfile(path, firstFileName);
                obj.files{i} = load(fullFilePath);
            end

        end


        function obj = mainFunction(obj)

            for i = 1:2
                if ~isempty(obj.files)
                    obj.ct = obj.files{1, i}.ct;
                    obj.cst = obj.files{1, i}.cst;
                end
                
                Xshifts = [251, 251];
                Yshifts = [417, 220];

                if ~isempty(obj.cst)
                    obj.Xshift = Xshifts(i);
                    obj.Yshift = Yshifts(i);
                    obj = createMask(obj); 
                    obj.bothContours{i} = obj.contours; 
                end

                if ~isempty(obj.bothContours{i})
                    obj = fillContours(obj);
                    obj.bothMaskSets{i} = obj.MasksSet;
                    obj.bothStructureSets{i} = obj.structures;
                    obj.bothStructureSets{i} = rmfield(obj.bothStructureSets{i}, 'Body');
                end

                if ~isempty(obj.bothStructureSets{i})

                    obj.parentDir1 = 'C:\Users\r114m\Documents\registration\RegistrationFolder'; % creating of the folder with all the structures
                    folderName = sprintf('Cubes_%01d',i);
                    FullFolderName = fullfile(obj.parentDir1, folderName);
                    mkdir(FullFolderName);
                    obj.parentDir2 = FullFolderName;

                    fieldNames = fieldnames(obj.bothStructureSets{i});

                    for j = 1:size(fieldNames)
                        subFolderName = sprintf('Cube_%01d',j); % creating of the folder with only one structure 
                        obj.FullSubFolderName = fullfile(obj.parentDir2, subFolderName);
                        obj.TwoFolders{i} = obj.parentDir2; % Folders with ref png data and synth png data
                        mkdir(obj.FullSubFolderName);
                        
                        obj.StructuresName = fieldNames{j};
                        obj.cubeHU = obj.bothStructureSets{i}.(obj.StructuresName);
                        obj = DcmToPng(obj);
                        
                    end
                end
            end

            obj = registration(obj);
            obj.mergedCube = zeros(obj.ct.cubeDim);
      
            for i = 1:length(fieldNames)
                InputFolder = 'C:\Users\r114m\Documents\registration\RegistrationFolder\PngRegStructures';
                NameOfFolder = sprintf('Cube_%01d', i);
                obj.inputFolder = fullfile(InputFolder, NameOfFolder);
                Fields = fieldnames(obj.bothStructureSets{1, 1});
                
                obj = PngToDicom(obj);

                % FURTHER PROCESSING %
                mask = zeros(obj.ct.cubeDim);
                maskWithoutContours1 = zeros(obj.ct.cubeDim);
                maskWithoutContours2 = zeros(obj.ct.cubeDim);

                % return filled voxels outside of the structure
                mask(~(obj.bothStructureSets{1, 1}.(Fields{i})(:) == min(obj.bothStructureSets{1, 1}.(Fields{i})(:)))) = 1;
                maskOfregImg = mask;
                obj.HUCube(~(mask == 1)) = obj.minHU;
                mask(obj.HUCube(:) ~= obj.minHU) = 1;
                obj.HUCube = obj.HUCube.*mask;
                obj.HUCube(obj.HUCube == 0) = obj.minHU;

                
                % border value issue
                mask(obj.HUCube(:) == obj.minHU) = 0;
                 
                for k = 1:obj.ct.cubeDim(3)

                    mask2d = mask(:,:,k);

                    boundaries = bwboundaries(mask2d);

                    for b = 1:length(boundaries)

                        boundary = boundaries{b};

                        for p = 1:size(boundary, 1)
                            x = boundary(p, 1); 
                            y = boundary(p, 2); 
                            mask2d(x, y) = 0; 
                        end
                    end

                    maskWithoutContours1(:,:,k) = mask2d;
                end

                for k = 1:obj.ct.cubeDim(3)

                    mask2d = maskWithoutContours1(:,:,k);

                    boundaries = bwboundaries(mask2d);

                    for b = 1:length(boundaries)

                        boundary = boundaries{b};

                        for p = 1:size(boundary, 1)
                            x = boundary(p, 1); 
                            y = boundary(p, 2); 
                            mask2d(x, y) = 0; 
                        end
                    end

                    maskWithoutContours2(:,:,k) = mask2d;
                end

                meanHU = mean(obj.HUCube(maskWithoutContours2(:) == 1)); % mean for synth structure after deformable registration
                
                for k = 1:obj.ct.cubeDim(3)

                    mask2d1 = mask(:,:,k);
                    mask2d2 = mask(:,:,k);

                    boundaries1 = bwboundaries(mask2d1);
                    boundaries2 = bwboundaries(mask2d2);

                    for b = 1:length(boundaries1)

                        boundary = boundaries1{b};

                        for p = 1:size(boundary, 1)
                            x = boundary(p, 1); 
                            y = boundary(p, 2); 
                            obj.HUCube(x, y, k) = meanHU; 
                        end
                    end
                    for b = 1:length(boundaries2)

                        boundary = boundaries2{b};

                        for p = 1:size(boundary, 1)
                            x = boundary(p, 1); 
                            y = boundary(p, 2); 
                            obj.HUCube(x, y, k) = meanHU; 
                        end
                    end
                end              


                % quantities to remove the difference

                obj.meanHUs.(NameOfFolder) = mean(obj.ct.cubeHU{1}(obj.MasksSet.(fieldNames{i})(:) == 1)); % mean for synth structure on the previos picture
                Difference =  obj.meanHUs.(NameOfFolder) - meanHU;
                obj.HUCube(obj.HUCube ~= obj.minHU) = obj.HUCube(obj.HUCube ~= obj.minHU) + Difference;
                
                minHUref = min(obj.files{1, 1}.ct.cubeHU{1}(:));
                % A = sum(~(obj.bothStructureSets{1, 1}.(Fields{i})(:) == minHUref)); % Number of all the voxels inside structure
                % B = sum(~(obj.HUCube(:) == obj.minHU)); 

                % fill 'empty' (min) value inside of structure
                % DifferenceHU =  (obj.meanHUs.(NameOfFolder)*A - meanHU*B)/(A - B);
                obj.HUCube(~(obj.bothStructureSets{1, 1}.(Fields{i})(:) == minHUref) & obj.HUCube(:) <= -110) = obj.meanHUs.(NameOfFolder);
               
                obj.deformStructures.(NameOfFolder) = obj.HUCube;
                
                obj.SetOfMasks.(NameOfFolder) = maskOfregImg;

            end
            fieldNames = fieldnames(obj.SetOfMasks);
            for i = 1:length(fieldNames)
                fieldName1 = fieldNames{i};
                for j = i+1:length(fieldNames)
                    fieldName2 = fieldNames{j};
                    intersection = obj.SetOfMasks.(fieldName1) & obj.SetOfMasks.(fieldName2);
                    obj.deformStructures.(fieldName2)(intersection) = 0;
                end
            end
            obj.mergedCube = zeros(obj.ct.cubeDim);
            mergedMasks = zeros(obj.ct.cubeDim);
            for i = 1:length(fieldNames)
                fieldName = fieldNames{i};
                obj.SetOfMasks.(fieldName)(obj.SetOfMasks.(fieldName) ~= 0) = 1;
                obj.mergedCube = obj.mergedCube + obj.SetOfMasks.(fieldName).*obj.deformStructures.(fieldName);
                mergedMasks = mergedMasks + obj.SetOfMasks.(fieldName);
            end
            mergedMasks(mergedMasks == 0) = min(obj.ct.cubeHU{1}(:));
            obj = AllOtherStructures(obj);
            obj.mergedCube = mergedMasks + obj.mergedCube + obj.OtherStructures; 
        end
        
        function obj = registration(obj)
            OutputFolder = 'PngRegStructures';
            FullOutputFolderName = fullfile(obj.parentDir1, OutputFolder);
            mkdir(FullOutputFolderName)
            for i = 1:numel(fieldnames(obj.contours)) - 1
                subFolderName = sprintf('Cube_%01d',i);
                fullSubFolderName = fullfile(FullOutputFolderName, subFolderName);
                mkdir(fullSubFolderName);
                for j = 1:obj.ct.cubeDim(3)
                    NameOfSubFolder = sprintf('Cube_%01d', i);
                    NameOfFile = sprintf('slice_%03d.png', j);% if length(unique(RefersenceImage)) ~= 1
                    
                    ReferenceImage = imread(fullfile(obj.TwoFolders{1}, NameOfSubFolder, NameOfFile)); 
                    SynthImage = imread(fullfile(obj.TwoFolders{2}, NameOfSubFolder, NameOfFile));

                    if length(unique(ReferenceImage)) ~= 1
                        if size(ReferenceImage, 3) == 3
                            ReferenceImage = rgb2gray(ReferenceImage);
                        end
                        if size(SynthImage, 3) == 3
                            SynthImage = rgb2gray(SynthImage);
                        end
                        [~, registeredImage] = imregdemons(SynthImage, ReferenceImage);
                    else
                        registeredImage = ReferenceImage;
                    end

                    filename = fullfile(fullSubFolderName, sprintf('slice_%03d.png', j));
                    imwrite(registeredImage, filename);
                end
            end
        end
        
        function obj = createMask(obj)

            NeededStructures = {'Femur head right', 'Femur head left', 'Sigma', 'Rectum', 'Bladder', 'CTV', 'Bones', 'BODY'};

            for k = 1:size(obj.cst, 1)
                obj.cube_for_mask = zeros(obj.ct.cubeDim);
                NameOfStruct = obj.cst{k, 2};
                Bool = any(strcmp(NeededStructures, NameOfStruct));
            
                if Bool
                    for i = 1:obj.ct.cubeDim(3)
                
                        if isempty(obj.cst{k, 7}{1, 1}) || length(obj.cst{k, 7}{1, 1}) < i || isempty(obj.cst{k, 7}{1, 1}{i, 3})
                            continue;
                        end
                        
                        contour_coords = obj.cst{k, 7}{1, 1}{i, 3};
                
                
                        if isempty(contour_coords) || size(contour_coords, 2) < 2
                            continue;
                        end
                
                
                        for j = 2:size(contour_coords, 2)
                            currentCoordX = contour_coords(1, j)*2 - obj.Xshift; %251 for synth
                            currentCoordY = contour_coords(2, j)*2 - obj.Yshift; %220 for synth
                            [~, idxY] = min(abs(obj.ct.x - currentCoordX));
                            [~, idxX] = min(abs(obj.ct.y - currentCoordY));
                
                            obj.cube_for_mask(idxX, idxY, i) = 1;
                        end
                        obj = removeIsolatedOnes(obj);
                    end
                    if strcmp(NameOfStruct, NeededStructures{1})
                            obj.contours.FHR = obj.cube_for_mask;
                        elseif strcmp(NameOfStruct, NeededStructures{2})
                            obj.contours.FHL = obj.cube_for_mask;
                        elseif strcmp(NameOfStruct, NeededStructures{3})
                            obj.contours.Sigma = obj.cube_for_mask;
                        elseif strcmp(NameOfStruct, NeededStructures{4})
                            obj.contours.Rectum = obj.cube_for_mask;
                        elseif strcmp(NameOfStruct, NeededStructures{5})
                            obj.contours.Bladder = obj.cube_for_mask;
                        elseif strcmp(NameOfStruct, NeededStructures{6})
                            obj.contours.CTV = obj.cube_for_mask;
                        elseif strcmp(NameOfStruct, NeededStructures{7})
                            obj.contours.Bones = obj.cube_for_mask;
                        elseif strcmp(NameOfStruct, NeededStructures{8})
                            obj.contours.Body = obj.cube_for_mask;
                    end
                end
            end
        end

        
        function obj = removeIsolatedOnes(obj)
        
            padded_matrix = padarray(obj.cube_for_mask, [1, 1], 0);
            cleaned_matrix = obj.cube_for_mask;
        
        
            for i = 1:size(obj.cube_for_mask, 1)
                for j = 1:size(obj.cube_for_mask, 2)
                    if obj.cube_for_mask(i, j) == 1
                        neighborhood = padded_matrix(i:i+2, j:j+2);
                        if sum(neighborhood(:)) == 1
                            cleaned_matrix(i, j) = 0;
                        end
                    end
                end
            end
        end


        function obj = fillContours(obj)
            obj.minHU = min(obj.ct.cubeHU{1}(:));
            mask = zeros(obj.ct.cubeDim);
            fieldNames = fieldnames(obj.contours);
            
            for i = 1:size(fieldNames)
                fieldName = fieldNames{i};
                for j = 1:obj.ct.cubeDim(3)
                    if length(unique(obj.contours.(fieldName)(:,:,j))) == 2
                        mask(:,:,j) = imfill(obj.contours.(fieldName)(:,:,j), "holes");
                    else
                        mask(:,:,j) = obj.contours.(fieldName)(:,:,j);
                    end
                end
                zeros_and_struct = mask.*obj.ct.cubeHU{1};
                % obj.meanHUs.mean
                mask(mask == 0) = obj.minHU; %minHU
                mask(mask == 1) = 0;
                obj.structures.(fieldName) = zeros_and_struct + mask;
                mask(mask == 0) = 1;
                mask(mask == -1000) = 0;
                mask(mask == -1024) = 0;
                obj.MasksSet.(fieldName) = mask;
                
            end
            
        end


        function obj = DcmToPng(obj)
            obj.minHU = min(obj.ct.cubeHU{1}(:));
            obj.maxHU = max(obj.ct.cubeHU{1}(:));
            
            
            cube_normalized = uint16(65535 * (obj.cubeHU - obj.minHU) / (obj.maxHU - obj.minHU));
            cube_normalized(cube_normalized < 0) = 0;
            cube_normalized(cube_normalized > 65535) = 65535;
            
            for i = 1:size(cube_normalized, 3)
                slice_2d = cube_normalized(:, :, i);
                filename = fullfile(obj.FullSubFolderName, sprintf('slice_%03d.png', i));
                imwrite(slice_2d, filename, 'BitDepth', 16);  
            end
        end

        function obj = PngToDicom(obj)
            obj.minHU = min(obj.ct.cubeHU{1}(:));
            obj.maxHU = max(obj.ct.cubeHU{1}(:));

            pngFiles = dir(fullfile(obj.inputFolder, '*.png'));

            firstSlice = imread(fullfile(obj.inputFolder, pngFiles(1).name));
            [rows, cols] = size(firstSlice);
            
            numSlices = length(pngFiles);
            obj.HUCube = zeros(rows, cols, numSlices);
      
            obj.maskOfRegImg = zeros(obj.ct.cubeDim);
            for i = 1:numSlices
                slicePath = fullfile(obj.inputFolder, pngFiles(i).name);
                slice_2d = imread(slicePath);
                doubleSlice2d = double(slice_2d);
                obj.maskOfRegImg(:,:,i) = doubleSlice2d;
                normalizedSlice2d = doubleSlice2d*(obj.maxHU - obj.minHU)/65535 + obj.minHU;
                obj.HUCube(:, :, i) = normalizedSlice2d;
            end
        end
        
        function obj = AllOtherStructures(obj)
            
            fieldNames = fieldnames(obj.bothMaskSets{1});
            SumOfMasks = zeros(obj.ct.cubeDim);
            obj.OtherStructures = zeros(obj.ct.cubeDim);

            for i = 1:length(fieldNames)
                if ~(strcmp(fieldNames{i}, 'Body'))
                    SumOfMasks = SumOfMasks + obj.bothMaskSets{1}.(fieldNames{i});
                end
            end

            BodyMask = obj.bothMaskSets{1}.Body;

            obj.OtherStructures(SumOfMasks(:) == 0 & BodyMask(:) == 1) = 1;
            obj.OtherStructures = obj.OtherStructures.*(971);
        end

    end

end