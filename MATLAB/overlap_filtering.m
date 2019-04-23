function convolution_struct = overlap_filtering(input, polyphase_filter, convolution_struct)

    %Getting the size of the convolution between the filter's coefficients and the input buffer
    convsize = length(input) + length(polyphase_filter) - 1;
    
    %Overlap filtering
    y = fftfilt([input zeros(1,length(polyphase_filter) - 1)],[polyphase_filter zeros(1, length(input) - 1)]);

    %Output 
    convolution_struct.buffer = zeros(1, length(input));
    for i = 1:length(input)
        convolution_struct.buffer(i) = y(i);
    end
    
    if length(convolution_struct.tail) == 1 %when tail is empty
        convolution_struct.tail = zeros(1,length(polyphase_filter) - 1);
    end    
    
    for i = 1:(length(polyphase_filter) - 1)
        convolution_struct.buffer(i) = convolution_struct.buffer(i) + convolution_struct.tail(i);
    end
    
    for i = (length(input) + 1):(convsize)
        convolution_struct.tail(i - length(input)) = y(i);
    end
    
end