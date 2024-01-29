
function states = generate_all_states(N, K)
    % Generate all K^N states for the cellular automata model.
    % Each state is represented as a row in the matrix 'states'.

    % Total number of states
    KhatN = K^N;
    
    % Initialize the matrix to store states
    states = zeros(KhatN, N);
    
    for i = 1:KhatN
        % Compute the ith state
        state = i-1;
        for j = 1:N
            % Convert the state to base K and store it
            states(i, N-j+1) = mod(state, K);
            state = floor(state / K);
        end
    end
end
