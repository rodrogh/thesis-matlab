clear all;
clc;
%Create sudoku - expert level (17] numbers)

board_size = 9;
square_size = 3;
init_numbers = [1 3];
total_squares = board_size/square_size;
sudoku = struct('board',zeros(board_size,board_size));

%Fill initial cells
num_appear = zeros(3,1);
total_inputs = 0;
square_numbers = zeros(3,3);
while sum(sum(square_numbers,1),2) < 17 %Keep filling board until desired number
    i = randi(3); 
    j = randi(3);
    while square_numbers(i,j) == 3
        i = randi(3); 
        j = randi(3);    
    end
    %Make sure the combination of different numbers do not exceed the
    %desired total number for desired level
    
    if num_appear(3) == 4
        if num_appear(2) == 2
            num_size = 1; %number of inputs in a given square
        elseif num_appear(1) == 1
            num_size = 2;
        else
            num_size = randi([1 2]);
        end
    else
        if num_appear(2) == 2
            num_size = randi([1 3]);
            while num_size == 2
                num_size = randi([1 3]);
            end
        elseif num_appear(1) == 1
            num_size = randi([2 3]);
        else
            num_size = randi([1 3]);
        end
    end
    %Count number appearences
    num_appear(1) = length(find(square_numbers==1));
    num_appear(2) = length(find(square_numbers==2));
    num_appear(3) = length(find(square_numbers==3));
    if total_inputs + num_size > 17
       num_size = num_size - (total_inputs + num_size - 17);
    end
    num_size = num_size - square_numbers(i,j);    
    p=0;
    while p < num_size
        board_i = randi(square_size) + (i-1)*square_size;
        board_j = randi(square_size) + (j-1)*square_size;
        %Check for repeated coordinates
        if sudoku.board(board_i,board_j) == 0
            sudoku.board(board_i,board_j) = randi(9);
            p = p + 1;
            square_numbers(i,j) = square_numbers(i,j) + 1; %counter for nmbers in square            
        end
    end
    total_inputs = sum(sum(square_numbers,1),2);
end

sudoku.board
total = find(sudoku.board);
