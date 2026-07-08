function chunk_analysis(timestamps, keypress_id, ktt_matrix)

if ~isempty(timestamps) && length(timestamps) > 1
    chunk = [];
    i_cluster = 1;
    for i = 1:length(kp)-1
        chunk(i) = i_cluster;
        if (diff(timestamps(i:i+1)) > ktt_matrix(keypress_id(i), keypress_id(i+1)))
            i_cluster = i_cluster + 1;
        end
    end
    chunk(end+1) = i_cluster;

elseif isempty(timestamps)
    chunk = [];
elseif length(timestamps) == 1
    chunk = 1;
end

end
