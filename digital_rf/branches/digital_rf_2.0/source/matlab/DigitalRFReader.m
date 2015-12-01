classdef DigitalRFReader
    % class DigitalRFReader allows easy read access to static Digital RF data
    %   See testDigitalRFReader.m for usage, or run <doc DigitalRFReader>
    %
    % $Id: DigitalRFReader.m 791 2015-07-07 17:43:03Z brideout $
    
    properties
        topLevelDirectories % a char array of one or more top level directories
        channel_map % a Map object with key=channel_name, value = drf_channel object
            
    end
    
    methods
        function reader = DigitalRFReader(topLevelDirectories)
            % DigitalRFReader is the contructor for this class.  
            % Inputs - topLevelDirectories - a char array of one or more 
            % top level directories, where a top level directory holds
            % channel directories
            
            
            % constants
            % define glob string for sub_directories in form YYYY-MM-DDTHH-MM-SS
            sub_directory_glob = '[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]-[0-9][0-9]-[0-9][0-9]';
            rf_file_glob = 'rf@[0-9]*.[0-9][0-9][0-9].h5';
            
            % topLevelDirectories - a char array of one or more top level 
            %   directories.
            if ~(ischar(topLevelDirectories))
                ME = MException('DigitalRFReader:invalidArg', ...
                  'topLevelDirectories arg not a string or char array');
                throw(ME)
            end
            dims = size(topLevelDirectories);
            if length(dims) == 1
                reader.topLevelDirectories = char(topLevelDirectories);
            else
                reader.topLevelDirectories = topLevelDirectories;
            end
            
            % make sure all exist
            dims = size(reader.topLevelDirectories);
            for i = 1:dims(1)
                if exist(strtrim(reader.topLevelDirectories(i,:)), 'dir') == 0
                    ME = MException('DigitalRFReader:invalidArg', ...
                        'topLevelDirectory %s not found', ...
                        reader.topLevelDirectories(i,:));
                    throw(ME)
                end
            end
            
            % the rest of this constructor fills out this map
            reader.channel_map = containers.Map();
            
            % fill out temp structure dirArr with fields: 1) top_level_dir, 
            %   2) channel, 3) subdirectory.  Ignore empty subdirectories
            dirFlag = 0;
            for i = 1:dims(1)
                topGlobPath = fullfile(strtrim(reader.topLevelDirectories(i,:)),'*', sub_directory_glob);
                top_level_dir_len = length(strtrim(reader.topLevelDirectories(i,:)));
                result = glob(topGlobPath);
                resultDims = size(result);
                for j = 1:resultDims(1)
                    data = char(result(j));
                    % make sure there are rf files
                    rfResult = glob(fullfile(data, rf_file_glob));
                    rfSize = size(rfResult);
                    if rfSize(1) == 0
                        % skip empty subdirectories
                        continue
                    end
                    remainder = data(top_level_dir_len + 2:end-1);
                    [pathstr,name,ext] = fileparts(remainder);
                    if dirFlag == 0
                        dirArr = struct('top_level_dir', strtrim(reader.topLevelDirectories(i,:)), ...
                            'channel', pathstr, 'subdirectory', name);
                        dirFlag = 1;
                    else
                        newArr = struct('top_level_dir', strtrim(reader.topLevelDirectories(i,:)), ...
                            'channel', pathstr, 'subdirectory', name);
                        dirArr(end+1) = newArr;
                    end
                end
                
            end
            
            % the next task is to order dirArr by subdirectory using
            % sortDirArr
            sortedDirArr = reader.sortDirArr(dirArr);
            
            % now loop through each unique channel name and get all ordered
            % subdirectories
            channels = char(dirArr.channel);
            chanDims = size(dirArr);
            unique_channels = unique(channels, 'rows');
            uniqueDims = size(unique_channels);
            top_levels = char(dirArr.top_level_dir);
            subdirs = char(dirArr.subdirectory);
            
            for i = 1:uniqueDims(1)
                thisChannel = strtrim(char(unique_channels(i,:)));
                sub_list = {}; % cell array of full paths in this channel
                for j = 1:chanDims(2)
                    thisItem = strtrim(char(dirArr(j).channel));
                    if strcmp(thisItem, thisChannel)
                        this_top_level = strtrim(top_levels(j,:));
                        this_subdir = strtrim(subdirs(j,:));
                        thisFullpath = fullfile(this_top_level, thisChannel, this_subdir);
                        sub_list{end+1} = thisFullpath;
                    end
                end
                % found all subdirectories for this channel - add it
                new_drf_channel = drf_channel(thisChannel, sub_list);
                reader.channel_map(thisChannel) = new_drf_channel;
            end % end loop through unique channel names
            
        end % end DigitalRFReader constructor
        
        
        
        function channels = get_channels(obj)
            % get_channels returns a cell array of channel names found
            % Inputs: None
            channels = keys(obj.channel_map);
        end % end get_channels
        
        
        
        function [lower_sample, upper_sample] = get_bounds(obj, channel)
            % get_bounds returns the first and last sample in channel.
            % sample bounds are in samples since 0 seconds unix time
            % (that is, unix time * sample_rate)
            drf_chan = obj.channel_map(channel);
            first_sub_struct = drf_chan.subdirectory_array{1};
            last_sub_struct = drf_chan.subdirectory_array{end};
            lower_sample = first_sub_struct.first_sample;
            upper_sample = last_sub_struct.last_sample;
        end
        
        
        
        function samples_per_file = get_samples_per_file(obj, channel)
            % get_samples_per_file returns samples_per_file for given channel
            drf_chan = obj.channel_map(channel);
            samples_per_file = drf_chan.samples_per_file;
        end % end get_channels
        
        
        
        function sample_rate = get_sample_rate(obj, channel)
            % get_sample_rate returns sample_rate in Hz for given channel
            drf_chan = obj.channel_map(channel);
            sample_rate = drf_chan.sample_rate;
        end % end get_sample_rate
        
        
        
        function is_complex = get_is_complex(obj, channel)
            % get_is_complex returns is_complex (1 or 0) for given channel
            drf_chan = obj.channel_map(channel);
            is_complex = drf_chan.is_complex;
        end % end get_is_complex
        
        
        
        function num_subchannels = get_num_subchannels(obj, channel)
            % get_num_subchannels returns num_subchannels (1 or greater) for given channel
            drf_chan = obj.channel_map(channel);
            num_subchannels = drf_chan.num_subchannels;
        end % end get_num_subchannels
        
        
        
        function gap_array = get_gap_array(obj, channel, first_sample, last_sample)
            % get_gap_array returns a N x 2 array of int64, where each row represents a data gap
            % with column 1 being the first missing
            % sample, and the 2 column being the last missing sample. May
            % be an empty array if no gaps found. Error raised if range of
            % first_sample, last_sample outside of available data as 
            % returned by get_bounds
            gap_array = int64.empty(0,2);
            [chan_first_sample, chan_last_sample] = obj.get_bounds(channel);
            if first_sample < chan_first_sample
                ME = MException('DigitalRFReader:invalidArg', ...
                  'first_sample %i before beginning of data at %i', ...
                  first_sample, chan_first_sample);
                throw(ME)
            end
            if last_sample > chan_last_sample
                ME = MException('DigitalRFReader:invalidArg', ...
                  'last_sample %i after end of data at %i', ...
                  last_sample, chan_last_sample);
                throw(ME)
            end
            % walk through ordered subdirectories to find all gaps
            drf_chan = obj.channel_map(channel);
            last_dir_end_sample = 0;
            for i = 1:length(drf_chan.subdirectory_array)
                this_struct = drf_chan.subdirectory_array{i};
                this_first = this_struct.first_sample;
                this_last = this_struct.last_sample;
                this_first_sample_list = this_struct.first_sample_list;
                this_rf_file_list = this_struct.rf_file_list;
                if this_last < first_sample
                    last_dir_end_sample = this_last;
                    continue;
                end
                % first see if there is a gap between subdirectories
                if i > 1
                    if (this_first - 1 > last_dir_end_sample) && (first_sample < this_first)
                        gap_array = cat(1, gap_array, ...
                            [max([last_dir_end_sample+1, first_sample]), ...
                            min([this_first-1, last_sample])]);
                    end
                end
                last_dir_end_sample = this_last;
                % analyze this subdirectory if  not continuous
                if this_struct.is_continuous == 0
                    this_sub_gaps = drf_chan.get_subdirectory_gaps(this_struct.fullpath, ...
                        first_sample, last_sample, obj.get_samples_per_file(channel), ...
                        this_first_sample_list, this_rf_file_list);
                    gap_array = cat(1, gap_array, this_sub_gaps);
                end
                % see if we are done
                if last_sample < this_last
                    break;
                end
            end
            
        end
        
        
        
        function continuous_array = get_continuous_array(obj, channel, first_sample, last_sample)
            % get_continuous_array returns a N x 2 array of int64, where each row represents a block of continuous data
            % with column 1 being
            % the first continuous
            % sample, and the 2 column being the last continuous sample. May
            % be an empty array if no data found. Error raised if range of
            % first_sample, last_sample outside of available data as 
            % returned by get_bounds
            gap_array = obj.get_gap_array(channel, first_sample, last_sample);
            continuous_array = int64.empty(0,2);
            gap_size = size(gap_array);
            if gap_size(1) == 0
                continuous_array = cat(1, continuous_array, ...
                            [first_sample, last_sample]);
            else
                for i = 1:gap_size(1)+1
                    if i == 1
                        if gap_array(1,1) > first_sample
                            continuous_array = cat(1, continuous_array, ...
                                [first_sample, gap_array(1,1)-1]);
                        end
                    elseif i < gap_size(1)+1
                        continuous_array = cat(1, continuous_array, ...
                            [gap_array(i-1,2)+1, gap_array(i,1)-1]);
                    else
                        if gap_array(i-1,2)+1 <= last_sample
                            continuous_array = cat(1, continuous_array, ...
                                [gap_array(i-1,2)+1, last_sample]);
                        end
                    end
                end % end loop through gaps
            end % end if
        end % end get_continuous_array
        
        
        
        function vector = read_vector(obj, channel, start_sample, sample_length)
            % read_vector returns a data vector sample_length x num_subchannels.  
            % Data type will be complex if data was complex, otherwise data type
            % as stored in same format as in Hdf5 file
            
            % verify no gaps
            last_sample = start_sample + (sample_length-1);
            data_gaps = obj.get_gap_array(channel, start_sample, last_sample);
            gap_size = size(data_gaps);
            if gap_size(1) > 0
                ME = MException('DigitalRFReader:invalidArg', ...
                  'data not continuous between %i and %i', ...
                  start_sample, last_sample);
                throw(ME)
            end
            % walk through all subdirectories to get data
            vector = [];
            drf_chan = obj.channel_map(channel);
            for i = 1:length(drf_chan.subdirectory_array)
                this_struct = drf_chan.subdirectory_array{i};
                this_first = this_struct.first_sample;
                this_last = this_struct.last_sample;
                this_fullpath = this_struct.fullpath;
                this_first_sample_list = this_struct.first_sample_list;
                this_rf_file_list = this_struct.rf_file_list;
                if this_last < start_sample
                    continue;
                end
                if this_first > last_sample
                    break;
                end
                this_vector = drf_chan.read_subdirectory(this_fullpath, ...
                    start_sample, last_sample, this_first_sample_list, ...
                    this_rf_file_list);
                if length(vector) > 0
                    vector = cat(1, vector, this_vector);
                else
                    vector = this_vector;
                end
            end % end subdirectory loop
        end
        
        
        
        function sortedDirArr = sortDirArr(obj, A)
            % sortDirArr is a private method that sorts struct array dirArr 
            % by field sundirectory, which puts them in time order
            % - input A is struct array with fields top_level_dir, channel,
            %  and subdirectory
            % uses algorithm in http://blogs.mathworks.com/pick/2010/09/17/sorting-structure-arrays-based-on-fields/
            Afields = fieldnames(A);
            Acell = struct2cell(A);
            sz = size(Acell);
            Acell = reshape(Acell, sz(1), []);
            Acell = Acell';
            Acell = sortrows(Acell, 3); % 3 is third field, subdirectory
            Acell = reshape(Acell', sz);
            sortedDirArr = cell2struct(Acell, Afields, 1);
        end % end sortedDirArr

    end % end methods
    
end % end DigitalRFReader calss


