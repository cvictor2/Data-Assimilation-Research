function result = resolution_change(data, new_size)
    % Get the size of the input data
    [original_size, ~] = size(data);

    % Create the new array of zeros
    result = zeros(new_size, new_size);

    % Calculate centers
    center_old = original_size / 2;
    center_new = new_size / 2;

    if new_size > original_size
        % Upsampling
        pos_indices = 1:original_size/2;
        neg_indices = original_size:-1:original_size/2+2;

        result(pos_indices,pos_indices) = data(pos_indices,pos_indices);
        result(pos_indices,new_size/2+neg_indices) = data(pos_indices,neg_indices);
        result(new_size/2+neg_indices,pos_indices) = data(neg_indices,pos_indices);
        result(new_size/2+neg_indices,new_size/2+neg_indices) = data(neg_indices,neg_indices);

        nyquist_freq = original_size/2+1;

        new_nyquist_x = zeros(1,new_size);
        new_nyquist_x(1:1:original_size/2+1) = data(nyquist_freq,1:1:original_size/2+1);
        new_nyquist_x(new_size - original_size/2+2:1:end) = data(nyquist_freq,original_size/2+2:1:end);

        new_nyquist_y = zeros(new_size,1);
        new_nyquist_y(1:1:original_size/2+1) = data(1:1:original_size/2+1,nyquist_freq);
        new_nyquist_y(new_size - original_size/2+2:1:end) = data(original_size/2+2:1:end,nyquist_freq);

        % result(nyquist_freq,1:1:old_size) = nyquist_data_x;
        % result(1:1:old_size,) = nyquist_data_x;

        result(nyquist_freq,:) = new_nyquist_x;
        result(:,nyquist_freq) = new_nyquist_y;
        result(new_size - nyquist_freq+2,:) = conj(new_nyquist_x);
        result(:, new_size - nyquist_freq+2) = conj(new_nyquist_y);
        result(new_size-nyquist_freq+2, new_size-nyquist_freq+2) = conj(data(nyquist_freq,nyquist_freq));

        result = result.*2*(new_size / original_size);

        % result(pos_indices, nyquist_freq) = data(pos_indices, nyquist_freq);
        % result(nyquist_freq, pos_indices) = data(nyquist_freq,pos_indices);
        % result(nyquist_freq, nyquist_freq) = data(nyquist_freq,nyquist_freq);

        % result(new_size/2+neg_indices,nyquist_freq) = data(neg_indices,pos_indices);
        % result(new_size/2+neg_indices,new_size/2+neg_indices) = data(neg_indices,neg_indices);


      
    else
        % Downsampling
        % Copy the original data to the center of the new array
        % result = data(center_old - center_new + 1:center_old + center_new, ...
                      % center_old - center_new + 1:center_old + center_new);
    end
end