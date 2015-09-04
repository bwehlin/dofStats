% Copyright (c) 2015, UC San Diego CURE Program
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [parking_functions, dof_dist] = getParkingFunctions( m, n )
    % Generates all Park^m_n functions and checks degrees of freedom
    %
    % Parameters:
    %   m                 - "maximum" generalized Shi plane
    %   n                 - dimension
    %
    % Output:
    %
    %   parking_functions - matrix of parking functions
    %   dof_dist          - distribution of degrees of freedom

    parking_functions = zeros(0,n);
    
    next_pf = ones(1, n);

    done = 0;
    
    % Generate parking functions
    while done == 0
        if next_pf ~= 0
            parking_functions(end+1, :) = next_pf;
            next_pf = nextParkingFunction(next_pf, m, n);
        else
            done = 1;
        end
    end
    
    % Discard permuted duplicates
    for i=1:size(parking_functions, 1)
        parking_functions(i, :) = sort(parking_functions(i, :));
    end
    
    parking_functions = unique(parking_functions, 'rows');
    
    
    % Check degrees of freedom
    
    dof_dist = zeros(n, 2);
    dof_dist(:,1) = transpose(1:n);
    
    for i=1:size(parking_functions, 1)
        parking_function = parking_functions(i,:);
        
        prime_components = countPrimeComponents(parking_function, m, n);
        
        dof_dist(prime_components, 2) = dof_dist(prime_components, 2) + ...
            countCombinations(parking_function, n);
    end
end

function next = nextParkingFunction(current, m, n)
    % The next parking function, ordered in increasing order
    %
    % Parameters:
    %   current - current parking function
    %   m       - "maximum" generalized Shi plane
    %   n       - dimension
    %
    % Output:
    %
    %   next    - next parking function

    next = current;

    % Last position i at which b_i < 1+m(i-1)
    last_pos_below = 1;
    
    for i=n:-1:2
        if current(1,i) < 1+m*(i-1)
            last_pos_below = i;
            break;
        end
    end
    
    if last_pos_below == 1
        next = 0;
        return;
    end
    
    % Reset everything after last_pos_below
    
    for i=last_pos_below+1:n
        next(1,i) = 1;
    end
    
    next(1,last_pos_below) = next(1,last_pos_below) + 1;
end

function prime_components = countPrimeComponents(parking_function, m, n)
    % Counts the number of prime components in the Dyck path corresponding
    % to a parking function
    %
    % Parameters:
    %   parking_function    - parking function
    %   m                   - "maximum" generalized Shi plane
    %   n                   - dimension
    %
    % Output:
    %
    %   prime_components    - the # of prime components in the
    %                         corresponding Dyck path
    
    prime_components = 1;
    
    for i=2:n
        if parking_function(1, i) == 1+m*(i-1)
            prime_components = prime_components + 1;
        end
    end
end

function combinations = countCombinations(parking_function, n)
    % Counts the number of unique combinations of the letters in a
    % parking function
    %
    % Parameters:
    %   parking_function - parking function
    %   n                - dimension
    %
    % Output:
    %
    %  combinations      - the number of possible combinations
    
    parking_histo = histc(parking_function, unique(parking_function));
    
    divisor = prod(arrayfun(@(x) factorial(x), parking_histo));
    
    combinations = factorial(n) / divisor;
end
