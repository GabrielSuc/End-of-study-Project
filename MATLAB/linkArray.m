function array3 = linkArray(array1,array2)
    % Use this function to combine two arrays made of groude delays and
    % related frequencies. The goal is to merge the two frequencies columns
    % by ascending order while keeping the groupe delay values associated to
    % each frequency.
    
    % The arrays should have this shape
    
    % | f1 gd1 |
    % | f2 gd2 |
    % |   .    |
    % |   .    |
    % |   .    |
    % | fn gdn |
    
   
    % One of the two arrays might already have duplicates. Remove them
    if length(array1) ~= length(unique(array1))
        C = unique(array1(:,1));
        n = C(histc(array1(:,1),C)>=2);
        indexes = find(ismember(array1(:,1),n)); % find the indics of same frequency values

        for l = 1:2:length(indexes)
            array1(indexes(l),2) = array1(indexes(l),2) + array1(indexes(l)+1,2);
            array1(indexes(l)+1,:) = 0;
        end

        % Removing the 0 
        array1(~all(array1(:,2),2),:) = [];

    end    

    if length(array2) ~= length(unique(array2))
        C = unique(array2(:,1));
        n = C(histc(array2(:,1),C)>=2);
        indexes = find(ismember(array2(:,1),n)); % find the indics of same frequency values

        for l = 1:2:length(indexes)
            array2(indexes(l),2) = array2(indexes(l),2) + array2(indexes(l)+1,2);
            array2(indexes(l)+1,:) = 0;
        end

        % Removing the 0 
        array2(~all(array2(:,2),2),:) = [];

    end
   
    % Verify that two array have same dimension. If they have not, zero-pad
    % the smallest. It will add a lot of 0 Hz frequency with 0 s delays.
    % Summing this delays will not affect the output then.  
    if size(array1,1) ~= size(array2,1)
        if size(array1,1) > size(array2,1)
            array2 = [array2; zeros(size(array1,1) - size(array2,1),2)];
            which_array = 2;
        else
            array1 = [array1; zeros(size(array2,1) - size(array1,1),2)];
            which_array = 1;
        end
        
    else
        which_array = 0;
    end
    

    
    % Number of lines in each table
    N1 = size(array1,1);
    N2 = size(array2,1);
    
    % Create a thrid array which will be the final one
    array3 = zeros(N1+N2,2);
    
    % Sart sorting the arrays
    i = 1;
    j = 1;
    k = 1;
    while i <= N1 && j <= N2
        
    % Check if current element of first 
    % array is smaller than current element 
    % of second array. If yes, store first  
    % array element and increment first array 
    % index. Otherwise do same with second array 

        if array1(i,1) < array2(j,1)
            array3(k,:) = array1(i,:); 
            k = k + 1;
            i = i + 1;
        else
            array3(k,:) = array2(j,:); 
            k = k + 1;
            j = j + 1;
        end
    end
    
    % Store remaining elements of first array
    while i <= N1
        array3(k,:) = array1(i,:);
        k = k + 1;
        i = i + 1;
    end    
        
    % Store remaining elements of second array
    while j <= N2
        array3(k,:) = array2(j,:);
        k = k + 1;
        j = j + 1;
    end 
    
    % Now that we have obtained the final array, we can remove the 0
    % introduced by the zero-padding
    % Removing the 0 
    array3(~all(array3(:,2),2),:) = [];
    
    if which_array == 1 
        array1(~all(array1(:,2),2),:) = [];
    elseif which_array == 2 
        array2(~all(array2(:,2),2),:) = [];
    end
    
    % Make sure there is not the same frequency in the array
    % Works this way because frequencies are sorted in final array
    if length(array3) ~= length(unique(array3))
        C = unique(array3(:,1));
        n = C(histc(array3(:,1),C)>=2);
        indexes = find(ismember(array3(:,1),n)); % find the indics of same frequency values
    end    
    
    % Sum group delay values for same frequencies. The frequencies in each
    % array are unique, thus the final array cannot contain more then 2
    % duplicates of the same frequency value. Therefore, two consecutive
    % values in indexes point to the same frequency.
    % We can still make sure it is the case 
    
    if (length(array1(:,1)) ~= length(unique(array1(:,1)))) || (length(array2(:,1)) ~= length(unique(array2(:,1))))
        error('The frequncies in one (or the two) arrays are not unique')
    end
    
    for l = 1:2:length(indexes)
        array3(indexes(l),2) = array3(indexes(l),2) + array3(indexes(l)+1,2);
        array3(indexes(l)+1,:) = 0;
    end
    
    % Removing the 0 
    array3(~all(array3(:,2),2),:) = [];

end