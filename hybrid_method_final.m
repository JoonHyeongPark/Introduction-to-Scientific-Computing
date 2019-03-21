function hybrid_method(square_dimension, matrix, b_vector)
  
  matrix
  b_vector
  augmented_matrix = [matrix, b_vector]
  
  for i = 1 : square_dimension, # initialize row pointer and row max.
    row_pointer(i) = i;
    row_max(i) = max(abs(matrix(i,:)));
    
    if(row_max(i) == 0), # if for some row max = 0, then there is no unique solution.
      printf("There is no unique solution.\n");
      return;
    end;
    
  end;
  
  row_max
  
  for i = 1 : square_dimension - 1, # elimination process to make upper triangular matrix.
    
    M = -1;
    pointer = -1;
    for j = i : square_dimension, 
      if abs(augmented_matrix(row_pointer(j), i)) / row_max(row_pointer(j)) > M, # scale pivotting, find optimal pivot.
        M = abs(augmented_matrix(row_pointer(j), i)) / row_max(row_pointer(j));
        pointer = j;
      end;
    end;
    
    if augmented_matrix(row_pointer(pointer), i) == 0, # if for some pivot = 0, then there is no unique solution.
      printf("There is no unique solution.\n");
      return;
    end;
    
    if row_pointer(pointer) != row_pointer(i), # interchange row by scaled pivotting.
      temp = row_pointer(i);
      row_pointer(i) = row_pointer(pointer);
      row_pointer(pointer) = temp;
    end;
    
    for j = i + 1 : square_dimension, # row elimination using diagonal pivot.
      multiplier = augmented_matrix(row_pointer(j), i) / augmented_matrix(row_pointer(i), i);
      augmented_matrix(row_pointer(j), i) = 0; # A[NROW[j]][i] = A[NROW[j]][i] - A[NROW[j]][i] / A[NROW[i]][i] * A[NROW[i]][i] is obviously zero.
      augmented_matrix(row_pointer(j), i + 1 : square_dimension + 1) -= multiplier * augmented_matrix(row_pointer(i), i + 1 : square_dimension + 1);
    end;

    # if you want more accuracy, use complete pivotting not partial pivotting.
    # just calculate row max again in here.
    # if you add this code below, you can get more accuracy but less speed using complete pivotting.
    
    # for j = i + 1 : square_dimension, row_max(row_pointer(j)) = max(abs(augmented_matrix(row_pointer(j), i + 1 : square_dimension))); end;
    
  end;
  
  if augmented_matrix(row_pointer(square_dimension), square_dimension) == 0, # if for some pivot = 0, then there is no unique solution.
    printf("There is no unique solution.\n")
    return;
  end;
  
  printf("Gaussian Elimination with Scaled Partial Pivotting has been finished.\n");
  
  printf("\n[Upper-trianglar Matrix, b vector] : \n");
  disp(augmented_matrix(row_pointer,:));
  
  iterator = square_dimension - 1 : -1 : 1;
  
  for i = iterator, # reducing process to make diagonal matrix.
    
    inverse_row_pointer(row_pointer(i)) = i;
    pivot = augmented_matrix(row_pointer(i + 1), i + 1);
    
    for j = 1 : i, # reduce column using diagonal pivot.
        multiplier = augmented_matrix(row_pointer(j), i + 1) / pivot;
        augmented_matrix(row_pointer(j), i + 1) = 0; # A[NROW[j]][i + 1] = A[NROW[j]][i + 1] - A[NROW[j]][i + 1] / A[NROW[i + 1]][i + 1] * A[NROW[i + 1]][i + 1] is obviously zero.
        augmented_matrix(row_pointer(j), square_dimension + 1) -= augmented_matrix(row_pointer(i + 1), square_dimension + 1) * multiplier;
    end;
    
  end;
  inverse_row_pointer(row_pointer(square_dimension)) = square_dimension;
  
  printf("\nReducing has been finished.\n");
  printf("\n[Diagonal Matrix, b vector] : \n");
  disp(augmented_matrix(row_pointer,:));

  for i = 1 : square_dimension, x_vector(i, 1) = augmented_matrix(i, square_dimension + 1) / augmented_matrix(i, inverse_row_pointer(i)); end; # calculate answer.
  
  printf("\nX vector has been calculated.\n");
  printf("\nX : \n");
  disp(x_vector(row_pointer));
  
end;