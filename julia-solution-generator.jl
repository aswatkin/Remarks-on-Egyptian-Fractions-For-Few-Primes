using Base

# Function to calculate the least common multiple (LCM)

# lcm(a,b) - calcluates the least common multiple of a and b  
# Input : a,b - integers (BigInt)
# Output : lcm (BigInt) - least common multiple of a and b
function lcm(a::BigInt, b::BigInt)::BigInt
    return abs(a * b) ÷ gcd(a, b)
end

# lcm_list(list) - calculates the least common multiple of a list of integers
# Input : list - list of integers (BigInt)
# Output : lcm (BigInt) - least common multiple of the list
function lcm_list(list::Vector{BigInt})::BigInt
    
    result = list[1]
    for i in 2:length(list)
        result = lcm(result, list[i])
    end
    return result
end
# get_proposed_bound (p,q,rank) - calculates the greedy algotithm bound for given parameters
# Input : p,q - integers (BigInt) - the two primes
#         rank - integer (Int) - the rank of the solution
# Output : proposed_bound (Int) - the proposed bound for the solution
function get_proposed_bound(p::BigInt, q::BigInt, rank::Int)::BigInt
    highest_power = 25
    first_denoms = BigInt[]

    for i in 0:highest_power-1
        for j in 0:highest_power-1
            if i + j > 0             
                push!(first_denoms, BigInt(p)^i * BigInt(q)^j)
            end
        end
    end

    first_denoms = sort(first_denoms)
    possible_denoms = unique(first_denoms)
    partial_solution = BigInt[]
    total_sum = Rational{BigInt}(0, 1)

    for denominator in possible_denoms
        if total_sum + Rational{BigInt}(1, denominator) < 1 && length(partial_solution) < rank - 1
            total_sum += Rational{BigInt}(1, denominator)
            push!(partial_solution, denominator)
        end
    end

    denom_bound = Rational{BigInt}(1, 1) / (1 - total_sum)

    counter = 1
    while possible_denoms[counter] < denom_bound
        counter += 1
    end
    proposed_bound = possible_denoms[counter]

    return proposed_bound
end

# sum_list(denom_list) - calculates the sum of the reciprocals of a list of integers
# Input : denom_list - list of integers (BigInt)
# Output : sum (Rational{BigInt}) - the sum the reciprocals of the list

function sum_list(denom_list::Vector{BigInt})::Rational{BigInt}
    running_sum = Rational{BigInt}(0, 1)
  
    if length(denom_list) == 0
        return Rational{BigInt}(0, 1)
    end
    
    denominator = BigInt(lcm_list(denom_list))  # Use BigInt for the denominator
    
    numerator = BigInt(0)
    for denom in denom_list
        
        numerator += denominator ÷ denom
    end
    return Rational{BigInt}(numerator, denominator)
end

#get_solutions(possible_denoms, rank, running_sol, list_of_solutions) - recursive function to find all solutions
# Input : possible_denoms - list of integers (BigInt) - the possible denominators
#         rank - integer (Int) - the rank of the solution
#         running_sol - list of integers (BigInt) - the current solution being built
#         list_of_solutions - list of lists of integers (BigInt) - the list of all solutions found
# Output : None - the function modifies list_of_solutions in place

function get_solutions(possible_denoms::Vector{Any}, rank::Int, running_sol::Vector{BigInt}, list_of_solutions::Vector{Any})
    running_sum = sum_list(running_sol)  # Returns Rational{BigInt}
    if running_sol == [2, 3, 8, 27, 243, 2048, 41472, 524288, 10616832]
       println("Found the solution: ", running_sol)
        print("Running sum is: ", running_sum, "\n")       
    end

    if running_sum == 1 && length(running_sol) < rank
        return
    end
    if running_sum > 1
        return
    elseif running_sum == 1 && length(running_sol) == rank
        push!(list_of_solutions, copy(running_sol))
        return
    elseif length(running_sol) >= rank
        return
    else
        m_n_frac = Rational{BigInt}(1, 1) - running_sum
        m, n = numerator(m_n_frac), denominator(m_n_frac)
        if m == 0
            println("Error: m is 0 during biggest_denom calculation.")
            return
        end
        s = length(running_sol)
        index = length(running_sol) > 0 ? findfirst(x -> x == running_sol[end], possible_denoms) : 1
        biggest_denom = (n * (rank - s)) ÷ m
        while index <= length(possible_denoms) && possible_denoms[index] <= biggest_denom
            push!(running_sol, possible_denoms[index])
            get_solutions(possible_denoms, rank, running_sol, list_of_solutions)
            pop!(running_sol)
            index += 1
        end
    end
end

function main()

    #Initialize  parameters
    p = 2
    q = 13
    rank = 7
    highest_power = 25


    # Generate list of possible denominators
    first_denoms = []
    println("Starting with p = ", p, " q = ", q, " rank = ", rank)
    for i in 0:highest_power-1
        for j in 0:highest_power-1
            if i + j > 0
                #print("2^", i, " * 3^", j, " = ", BigInt(p^i * q^j), " should be, ", BigInt(2^i * 3^j), "\n")
                push!(first_denoms, BigInt(BigInt(p)^i * BigInt(q)^j))
            end
        end
    end

    possible_denoms = sort(unique(first_denoms))


    list_of_solutions = []
  
    #Generate all possible solutions
    get_solutions(possible_denoms, rank, BigInt[], list_of_solutions)

    
    max_denominator = 0
    for solution in list_of_solutions
        for item in solution
            if item > max_denominator
                max_denominator = item
            end
        end
    end



    # Find the prime factorization of the maximum denominator

    exp_p, exp_q = 0, 0
    to_factor = max_denominator

    if max_denominator == 0
        println("Error: max_denominator is 0. No valid solutions found.")
        return
    end
    while to_factor % p == 0
        to_factor ÷= p
        exp_p += 1
    end

    while to_factor % q == 0
        to_factor ÷= q
        exp_q += 1
    end

    #write the list of solution to a file 
    open("solutions.txt", "w") do file
        for solution in list_of_solutions
            write(file, join(solution, ", "), "\n")  # Write each solution as a comma-separated line
        end
    end
    println("Solutions written to solutions.txt")

    print("There are ", length(list_of_solutions), " solutions\n")

    # UNCOMMENT TO PRINT ALL SOLUTIONS
    #for solution in list_of_solutions
        #print( solution, "\n")
   # end

    println("The maximum denominator that appears is: ", max_denominator, " which is 2^", exp_p, " * ", q, "^", exp_q)

    greedy_bound = get_proposed_bound(BigInt(p), BigInt(q), rank)
    println("The greedy bound is: ", greedy_bound)
    if max_denominator > greedy_bound
        println("WARNING: The greedy bound fails in this case")
    else
        println("The maximum denominator is less than the greedy bound")
    end

end 

main()
