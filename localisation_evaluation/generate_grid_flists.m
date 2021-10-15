function generate_grid_flists

% number of speakers
nspeakers = 34;

% Wav file path
[~, twoearsRoot] = get_data_root;
wav_path = fullfile(twoearsRoot, 'twoears-data/sound_databases/grid_subset');

% Define angular resolution
azRes = 5;

% All possible azimuth angles
angles = 0:azRes:(360-azRes);
nangles = length(angles);

% number of files per azimuth angle
nfiles_train = 2;
nfiles_test = 2;
nfiles = nfiles_train + nfiles_test;

fnames = cell(nfiles, nspeakers, nangles);

% Loop over all speakers
for spkid = 1:nspeakers
    
    fprintf('speaker %d\n', spkid);
    
    spkwav_path = sprintf('%s/s%d', wav_path, spkid);
    
    % Retrieve all file names per speaker
    all_files = dir(fullfile(spkwav_path, '*.wav'));
    
    % Randomly select nfiles
    file_indices = randperm(length(all_files));
    
    % Loop over all azimuths
    for n = 1:nangles
        % Loop over all files
        for k = 1:nfiles
            wavfn = all_files(file_indices((n-1)*nfiles+k)).name;
            fnames{k, spkid, n} = wavfn(1:end-4);
        end
    end
end

% Write file lists
fid_train = fopen(sprintf('flist_train_%ddeg.txt',azRes), 'w');
fid_test = fopen('flist_test.txt', 'w');

for n = 1:nangles
    % Loop over all speakers
    for spkid = 1:nspeakers
        % Loop over all files
        for k = 1:nfiles
            if k <= nfiles_train
                fprintf(fid_train, 'az%d/s%d_%s\n', angles(n), spkid, fnames{k, spkid, n});
            else
                fprintf(fid_test, 's%d_%s\n', spkid, fnames{k, spkid, n});
            end
        end
    end
end
fclose(fid_train);
fclose(fid_test);
