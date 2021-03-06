classdef drf_channel
    % class drf_channel is a private class to describe a single
    % DigitalRFReader channel - do not create directly
    %
    % $Id: drf_channel.m 756 2015-04-03 15:20:19Z brideout $
    
    properties
        channel % channel name (string)
        samples_per_file % samples per file in this drf_channel
        sample_rate % sample rate in Hz in this drf_channel
        is_complex % 1 if channel has real and imag data, 0 if real only
        num_subchannels % number of subchannels - 1 or greater
        subdirectory_array % an ordered struct array of subdirectory metadata
            % fields - fullpath, first_sample, last_sample, is_continuous,
            % first_sample_list, rf_file_list
            % first_sample_list is an ordered list of first sample in each rf
            % file as determined by file name convention only
            % rf_file_list is an ordered list of rf_file basenames
    end
    
    methods
        function channel = drf_channel(channel, subdirectory_list)
            % drf_channel constructor
            % Inputs:
            %   channel -the channel name (string)
            %   subdirectory_array - a sorted char array of full paths to 
            %       each subdirectory.
            channel.channel = channel;
            channel.subdirectory_array = {};
            for i = 1:length(subdirectory_list)
                [sample_rate, samples_per_file, is_complex, num_subchannels, ...
                    start_index, end_index, is_continuous] = ...
                    channel.get_subdirectory_metadata(char(subdirectory_list(1,i)));
                first_sample_list = channel.get_first_sample_list(char(subdirectory_list(1,i)), sample_rate);
                rf_file_list = channel.get_rf_file_list(char(subdirectory_list(1,i)));
                channel.subdirectory_array{end+1} = struct('fullpath', char(subdirectory_list(1,i)), ...
                    'first_sample', start_index, 'last_sample', end_index, ...
                    'is_continuous', is_continuous, 'first_sample_list', first_sample_list, ...
                    'rf_file_list', rf_file_list);
            end
            channel.samples_per_file = samples_per_file;
            channel.sample_rate = sample_rate;
            channel.is_complex = is_complex;
            channel.num_subchannels = num_subchannels;
        end
        
        
        function [sample_rate, samples_per_file, is_complex, num_subchannels, ...
                start_index, end_index, is_continuous] ...
            = get_subdirectory_metadata(obj, subdirectory)
            % get_subdirectory_metadata returns metadata about one
            %   subdirectory
            rf_file_glob = 'rf@[0-9]*.[0-9][0-9][0-9].h5';
            rfResult = glob(fullfile(subdirectory, rf_file_glob));
            rfSize = size(rfResult);
            if rfSize(1) == 0
                ME = MException('drf_channel:invalidArg', ...
                    'No rf files found in %s', ...
                    subdirectory);
                throw(ME)
            end
            rfResult = sort(rfResult);
            % get info from first and last file
            [sample_rate, samples_per_file, is_complex, num_subchannels, ...
                start_index, end_index_1] = ...
                    obj.get_rf_file_metadata(char(rfResult(1,:)));
            [sample_rate, samples_per_file, is_complex, num_subchannels, ...
                start_index_2, end_index] = ...
                    obj.get_rf_file_metadata(char(rfResult(end,:)));
            if samples_per_file * length(rfResult) == (end_index - start_index) + 1
                is_continuous = 1;
            else
                is_continuous = 0;
            end
        end
        
        
        function [sample_rate, samples_per_file, is_complex, num_subchannels, ...
                start_index, end_index] ...
            = get_rf_file_metadata(obj, rf_file)
            % get_rf_file_metadata returns metadata about a single rf file
            %   by using hdf5 read methods
            sample_rate = h5readatt(rf_file, '/rf_data', 'sample_rate');
            samples_per_file = h5readatt(rf_file, '/rf_data', 'samples_per_file');
            is_complex = h5readatt(rf_file, '/rf_data', 'is_complex');
            num_subchannels = h5readatt(rf_file, '/rf_data', 'num_subchannels');
            index_data = h5read(rf_file, '/rf_data_index');
            start_index = index_data(1,1);
            end_index = index_data(1,end) + ((samples_per_file - index_data(2,end))-1);
        end
        
        
        function gap_array = get_subdirectory_gaps(obj, subdirectory, ...
                first_sample, last_sample, samples_per_file, first_sample_list, ...
                rfBasenameList)
            % get_subdirectory_metadata returns a gap array for one
            % subdirectory. Returns a N x 2 
            % array of int64, where each row represents a data gap, with 
            % column 1 being the first missing sample, and the 2 column 
            % being the last missing sample.  Gaps outside of range of
            % first_sample and last sample are excluded. first_sample_list
            % is an ordered list
            % of first samples for all rf files in fullpath (used to speed
            % data access) rfBasenameList is a char array of ordered rf
            % file basnames
            gap_array = int64.empty(0,2);
            start_file_index = obj.get_start_file_index(first_sample, first_sample_list);
            % get gap info from each file
            rfSize = size(rfBasenameList)
            for i = start_file_index:rfSize(1)
                rfFileName = fullfile(subdirectory, char(rfBasenameList(i,:)));
                [sample_rate, samples_per_file, is_complex, num_subchannels, ...
                    start_index, end_index] = ...
                        obj.get_rf_file_metadata(rfFileName);
                if end_index < first_sample
                    last_end_sample = end_index;
                    continue;
                end
                if i > start_file_index
                    if (last_end_sample + 1 < start_index) && (first_sample < start_index)
                        % gap found between samples
                        gap_array = cat(1, gap_array, ...
                            [max([last_end_sample+1, first_sample]), ...
                            min([start_index-1, last_sample])]);
                    end
                end
                last_end_sample = end_index;
                index_data = h5read(rfFileName, '/rf_data_index');
                % add gaps from index_data
                this_last_end_sample = 0;
                index_size = size(index_data);
                for j = 1:index_size(2)+ 1
                    if j < index_size(2)+ 1
                        this_start_sample = index_data(1,j);
                        this_start_index = index_data(2,j);
                    end
                    if j == 1
                        last_start_sample = this_start_sample;
                        last_start_index = this_start_index;
                        continue
                    elseif j < index_size(2)+ 1
                        active_start_sample = last_start_sample;
                        active_end_sample = last_start_sample + ((this_start_index - last_start_index)-1);
                        last_start_sample = this_start_sample;
                        last_start_index = this_start_index;
                    else
                        active_start_sample = last_start_sample;
                        active_end_sample = last_start_sample + (samples_per_file - last_start_index);
                    end
                    if active_end_sample < first_sample
                        this_last_end_sample = active_end_sample;
                        continue;
                    end
                    if this_last_end_sample >= last_sample
                        break;
                    end
                    if (first_sample < active_start_sample) && (j>2)
                        % gap found between samples
                        gap_array = cat(1, gap_array, ...
                            [max([this_last_end_sample+1, first_sample]), ...
                            min([active_start_sample-1, last_sample])]);
                    end
                    this_last_end_sample = active_end_sample;
                    % see if we are done
                    if active_end_sample > last_sample
                        break;
                    end
                end % done loop through /rf_data_index
                % see if we are done
                if end_index > last_sample
                    break;
                end
            end % end rf file loop
        end
        
        
        function first_sample_list = get_first_sample_list(obj, subdirectory, sample_rate)
            % get_first_sample_list returns an array of first sample values
            % based only on file naming convention
            rf_file_glob = 'rf@[0-9]*.[0-9][0-9][0-9].h5';
            rfResult = glob(fullfile(subdirectory, rf_file_glob));
            rfSize = size(rfResult);
            if rfSize(1) == 0
                ME = MException('drf_channel:invalidArg', ...
                    'No rf files found in %s', ...
                    subdirectory);
                throw(ME)
            end
            rfResult = sort(rfResult);
            first_sample_list = zeros(length(rfResult), 1,'int64');
            for i = 1:length(rfResult)
                [pathstr,name,ext] = fileparts(char(rfResult(i,:)));
                timeStamp = str2double(name(4:end));
                samples = floor(timeStamp * sample_rate);
                first_sample_list(i,1) = samples;
            end
        end % end get_first_sample_list
        
        
        function rf_file_list = get_rf_file_list(obj, subdirectory)
            % get_rf_file_list returns an array of sorted rf_file basenames
            rf_file_glob = 'rf@[0-9]*.[0-9][0-9][0-9].h5';
            rfResult = glob(fullfile(subdirectory, rf_file_glob));
            rfSize = size(rfResult);
            if rfSize(1) == 0
                ME = MException('drf_channel:invalidArg', ...
                    'No rf files found in %s', ...
                    subdirectory);
                throw(ME)
            end
            rfResult = sort(rfResult);
            rf_file_list = [];
            for i = 1:length(rfResult)
                [pathstr,name,ext] = fileparts(char(rfResult(i,:)));
                rfFile = strcat(name, ext);
                rf_file_list = cat(1, rf_file_list, rfFile);
            end
        end % end get_rf_file_list
        
        
        function start_file_index = get_start_file_index(obj, start_sample, first_sample_list)
            % get_start_file_index returns the index before the index in
            % first_sample_list where start_sample >= first_sample_list(i)
            % for now simple loop
            start_file_index = 0; % determine if value found
            for i = 1:length(first_sample_list)
                if first_sample_list(i,1) >= start_sample
                    if i > 1
                        start_file_index = i-1;
                    else
                        start_file_index = 1;
                    end
                    break
                end
            end
            if start_file_index == 0
                start_file_index = length(first_sample_list);
            end
        end
        
        
        function vector = read_subdirectory(obj, fullpath, start_sample, last_sample, ...
                first_sample_list, rfBasenameList)
            % read_subdirectory returns a vector with all data between
            % start_sample and last_sample for the data in fullpath. Data 
            % type will be complex if data was complex, otherwise data type
            % of that stored in Hdf5. first_sample_list is an ordered list
            % of first samples for all rf files in fullpath (used to speed
            % data access) rfBasenameList is a cell array of ordered rf
            % file basenames
            vector = []; 
            start_file_index = obj.get_start_file_index(start_sample, first_sample_list);
            % get data from each file
            rfSize = size(rfBasenameList);
            for i = start_file_index:rfSize(1)
                rfFileName = fullfile(fullpath, char(rfBasenameList(i,:)));
                [sample_rate, samples_per_file, is_complex, num_subchannels, ...
                    start_index, end_index] = ...
                        obj.get_rf_file_metadata(rfFileName);
                if end_index < start_sample
                    continue;
                end
                if start_index > last_sample
                    break;
                end
                % we need data from this file
                data = h5read(rfFileName, '/rf_data');
                index_data = h5read(rfFileName, '/rf_data_index');
                if obj.is_complex
                    this_vector = complex(data.r, data.i);
                else
                    this_vector = data;
                end
                this_vector = transpose(this_vector);
                % see if we need to slice
                if start_sample > start_index || last_sample < end_index
                    % we only need a slice - find slice_start and slice_end
                    % default
                    slice_start = 1;
                    slice_end = samples_per_file;
                    index_size = size(index_data);
                    for j = 1:index_size(2)+1
                        if j < index_size(2)+ 1
                            this_start_sample = index_data(1,j);
                            this_start_index = index_data(2,j);
                        end
                        if j == 1
                            last_start_sample = this_start_sample;
                            last_start_index = this_start_index;
                            continue
                        elseif j < index_size(2)+ 1
                            active_start_sample = last_start_sample;
                            active_end_sample = last_start_sample + ((this_start_index - last_start_index)-1);
                            active_start_index = last_start_index;
                            last_start_sample = this_start_sample;
                            last_start_index = this_start_index;
                        else
                            active_start_sample = last_start_sample;
                            active_end_sample = last_start_sample + (samples_per_file - last_start_index);
                            active_start_index = last_start_index;
                        end
                        if start_sample >= active_start_sample && start_sample <= active_end_sample
                            slice_start = (active_start_index + (start_sample - active_start_sample)) + 1;
                        end
                        if last_sample >= active_start_sample && last_sample <= active_end_sample
                            slice_end = (active_start_index + (last_sample - active_start_sample)) + 1;
                        end
                    end % end loop through index_data
                    this_vector = this_vector(slice_start:slice_end, :);
                    
                end % end slice if logic
                if length(vector) == 0
                    vector = this_vector;
                else
                    vector = cat(1, vector, this_vector);
                end
            end % end rf file loop
            
        end % end read_subdirectory
   
        
    end % end methods
end % end class

