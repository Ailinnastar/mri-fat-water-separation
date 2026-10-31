r = [0.087, 0.694, 0.128, 0.004, 0.039, 0.048];
fF_m = [-0.498, -0.447, -0.345, -0.261, -0.063, 0.064]* 1e3;
subjectFolder = getenv('MRI_DICOM_DIR');
if isempty(subjectFolder) || ~isfolder(subjectFolder)
    error('Set MRI_DICOM_DIR to an authorised local DICOM directory.');
end
[Data4D, rdataEchoTimes, rdataMetadata, rdataSlicePositions] = loadDicomSeries4D(subjectFolder);
fprintf('Loaded %d slices with %d echo times\n', size(Data4D, 3), size(Data4D, 4));

rows = size(Data4D, 1);
cols = size(Data4D, 2);
sliceIndex = str2double(getenv('MRI_SLICE_INDEX'));
if isnan(sliceIndex)
    sliceIndex = min(10, size(Data4D, 3));
end
if ~isscalar(sliceIndex) || ~isfinite(sliceIndex) || sliceIndex < 1 || ...
        sliceIndex ~= floor(sliceIndex) || sliceIndex > size(Data4D, 3)
    error('MRI_SLICE_INDEX must identify an available slice.');
end
rdataFF_map = zeros(rows, cols);
rdataR2s_map = zeros(rows, cols);

no_sig_area = abs(Data4D(1:min(30, rows), 1:min(30, cols), sliceIndex, 1));
totalVoxels = rows * cols;
D = parallel.pool.DataQueue;
afterEach(D, @(~) updateProgress(totalVoxels));


% Use magnitude of raw data for fitting (set true only to simulate low-SNR for testing)
ADD_SYNTHETIC_NOISE = false;
noise_floor = max(std(no_sig_area(:)), 1e-6);  % avoid division by zero in SNR

parfor row = 1:rows
    for col = 1:cols
        raw_cpx = squeeze(Data4D(row, col, sliceIndex, :));
        if ADD_SYNTHETIC_NOISE
            snr_target = 70;
            signal = abs(raw_cpx + (mean(abs(raw_cpx)) / snr_target)*(randn(size(raw_cpx)) + 1i*randn(size(raw_cpx))));
        else
            signal = abs(raw_cpx);
        end
        snr = mean(signal) / noise_floor;
        R2s = 0;
        pw = 0; pf = 0;

        if mean(signal) > 8 && snr > 1
            [pw, pf, R2s] = fitMultiPeakDixonModel(signal, rdataEchoTimes, r, fF_m);

            if (pw + pf) > 0 && snr > 5
                FF = pf / (pw + pf);
                FF = min(1, max(0, FF));
                rdataFF_map(row, col) = FF * 100;
            end

            if snr > 2 && R2s >= 50 && R2s <= 90
                rdataR2s_map(row, col) = R2s;
            end
        end

        send(D, 1);
    end
end

rdataFF_map = medfilt2(rdataFF_map, [3 3]);
rdataFF_map = imgaussfilt(rdataFF_map, 0.5);
rdataR2s_map = medfilt2(rdataR2s_map, [5 5]);
rdataR2s_map = imgaussfilt(rdataR2s_map, 1.2);

figure('Name', sprintf('Dixon MRI Results - Slice %d', sliceIndex), 'Position', [50, 50, 1200, 500]);


subplot(1,2,1);
imagesc(rdataFF_map, [0 100]);
colormap(gca, jet);
colorbar;
title(sprintf('Fat Fraction (FF) Map - Slice %d', sliceIndex));
axis image; axis off;


subplot(1,2,2);
imagesc(rdataR2s_map, [0 450]);
colormap(gca, jet);
colorbar;
title(sprintf('R2* Map - Slice %d', sliceIndex));
axis image; axis off;

saveas(gcf, sprintf('rdata_Slice%d.png', sliceIndex));

function updateProgress(total)
    persistent cnt t0
    if isempty(cnt)
        cnt = 0; t0 = tic;
    end
    cnt = cnt + 1;
    if mod(cnt, max(1,floor(total/100)))==0
        fprintf('Progress: %5.1f%% (%d/%d) \n', 100*cnt/total, cnt, total);
    end
end

function [data4D, echoTimes, metadata, slicePositions] = loadDicomSeries4D(folderPath)
    allFiles = dir(fullfile(folderPath, '*.dcm'));
    if isempty(allFiles)
        error('No .dcm files were found in MRI_DICOM_DIR.');
    end
    fileInfos = cell(length(allFiles), 1);
    for i = 1:length(allFiles)
        fileInfos{i} = dicominfo(fullfile(folderPath, allFiles(i).name));
    end

    magnitudeIndices = [];
    for i = 1:length(fileInfos)
        if isfield(fileInfos{i}, 'ImageType')
            imgType = fileInfos{i}.ImageType;
            isMagnitude = false;

            if ischar(imgType)
                isMagnitude = contains(imgType, 'M_FFE\M') || (contains(imgType, 'M_FFE') && ~contains(imgType, 'I_FFE') && ~contains(imgType, 'R_FFE'));
            end

            if isMagnitude
                magnitudeIndices = [magnitudeIndices, i];
            end
        end
    end

    if ~isempty(magnitudeIndices)
        fprintf('Found %d magnitude images using M_FFE pattern in ImageType.\n', length(magnitudeIndices));
    end

    selectedFiles = allFiles(magnitudeIndices);
    selectedInfos = fileInfos(magnitudeIndices);
    if isempty(selectedFiles)
        error('No supported M_FFE magnitude images were found.');
    end
    numFiles = length(selectedFiles);
    allSliceLocations = zeros(numFiles, 1);
    allEchoTimes = zeros(numFiles, 1);

    for i = 1:numFiles
        allEchoTimes(i) = selectedInfos{i}.EchoTime / 1000;
        allSliceLocations(i) = selectedInfos{i}.SliceLocation;
    end
    uniqflie = unique(allSliceLocations);
    uniqet = unique(allEchoTimes);

    if length(uniqflie) > length(allFiles)/3

        [sortedPositions, sortOrder] = sort(uniqflie);
        startIdx = floor(length(sortedPositions) / 3) + 1;
        endIdx = 2 * (startIdx-1);

        mid = sortedPositions(startIdx:endIdx);
        uniqflie = mid;

        validIndices = [];
        for i = 1:numFiles
            if ismember(allSliceLocations(i), mid)
                validIndices = [validIndices, i];
            end
        end

        selectedFiles = selectedFiles(validIndices);
        selectedInfos = selectedInfos(validIndices);
        numFiles = length(selectedFiles);
        allSliceLocations = allSliceLocations(validIndices);
        allEchoTimes = allEchoTimes(validIndices);

    end

    [uniqflie, ~, sliceIndices] = unique(allSliceLocations);
    [uniqflie, sortIdx] = sort(uniqflie);
    idxMap = zeros(size(uniqflie));
    for i = 1:length(sortIdx)
        idxMap(sortIdx(i)) = i;
    end
    sliceIndices = idxMap(sliceIndices);

    sampleFilePath = fullfile(folderPath, selectedFiles(1).name);
    sampleData = dicomread(sampleFilePath);
    [rows, cols] = size(sampleData);
    numSlices = length(uniqflie);
    numEchoTimes = length(uniqet);
    data4D = zeros(rows, cols, numSlices, numEchoTimes);
    metadata = cell(numEchoTimes, numSlices);

    for i = 1:numFiles
        echoIdx = find(uniqet == allEchoTimes(i));
        sliceIdx = sliceIndices(i);
        imgData = double(dicomread(fullfile(folderPath, selectedFiles(i).name)));

        if isfield(selectedInfos{i}, 'RescaleSlope') && isfield(selectedInfos{i}, 'RescaleIntercept')
            rescaleSlope = selectedInfos{i}.RescaleSlope;
            rescaleIntercept = selectedInfos{i}.RescaleIntercept;
            imgData = imgData * rescaleSlope + rescaleIntercept;
        end

        data4D(:, :, sliceIdx, echoIdx) = imgData;
        metadata{echoIdx, sliceIdx} = selectedInfos{i};


    end

    echoTimes = uniqet;
    slicePositions = uniqflie;

    fprintf('\nSuccessfully loaded 4D dataset with %d slices and %d echo times\n', numSlices, numEchoTimes);
end

function [pw, pf, R2s] = fitMultiPeakDixonModel(signal, echoTimes, r, fF_m)
    signal_max = max(abs(signal));
    signal = abs(signal) / max(signal_max, eps);
    signal_decay = abs(signal(end)) / abs(signal(1));
    est_R2 = -log(signal_decay) / (echoTimes(end) - echoTimes(1));
    est_R2 = min(90, max(50, est_R2));
    dom1 = [1, 0, est_R2];
    dom2 = [0, 1, est_R2];

    starts = [[1, 0, est_R2]; [0, 1, est_R2]; repmat(dom1, 2, 1) + 0.01 * (2*rand(2,3) - 1); repmat(dom2, 2, 1) + 0.01 * (2*rand(2,3) - 1)];
    starts = max([0, 0, 50], min(starts, [1, 1, 90]));

    lb = [0, 0, 50];
    ub = [1, 1, 90];

    opts = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 1000, 'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8);

    try
        bval = inf;
        bp = [0, 0, 0];

        for i = 1:size(starts, 1)
            objectiveFunc = @(params) sum((abs(signal) - modelObjective(params, echoTimes, r, fF_m)).^2);

            [params, fval] = fmincon(objectiveFunc, starts(i,:), [], [], [], [], lb, ub, [], opts);

            if fval < bval
                bval = fval;
                bp = params;
            end
        end

        pw = bp(1);
        pf = bp(2);
        R2s = bp(3);

        fit_quality = bval / (signal_max * length(signal));
        confidence = exp(-fit_quality * 5);
        R2s = min(90, max(50, R2s));

    catch
        pw = NaN;
        pf = NaN;
        R2s = NaN;
    end

    function signalMag = modelObjective(p, t, r, fF_m)
        pw = p(1);
        pf = p(2);
        R2s = p(3);

        fsig = zeros(size(t));
        for m = 1:length(r)
            fsig = fsig + r(m) * exp(1i*2*pi*fF_m(m)*t);
        end

        signalMag = abs((pw + pf * fsig) .* exp(-t*R2s));
    end
end
