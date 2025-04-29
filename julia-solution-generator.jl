using Base

# Function to calculate the least common multiple (LCM)
function lcm(a::BigInt, b::BigInt)::BigInt
    return abs(a * b) ÷ gcd(a, b)
end
function lcm_list(list::Vector{BigInt})::BigInt
    
    result = list[1]
    for i in 2:length(list)
        result = lcm(result, list[i])
    end
    return result
end
# This does the greedy algorithm
function get_proposed_bound(p::Int, q::Int, rank::Int)::Int
    highest_power = 25
    first_denoms = []

    for i in 0:highest_power-1
        for j in 0:highest_power-1
            if i + j > 0             
                push!(first_denoms, p^i * q^j)
            end
        end
    end

    first_denoms = sort(first_denoms)
    possible_denoms = unique(first_denoms)
    partial_solution = []
    total_sum = 0

    for denominator in possible_denoms
        if total_sum + 1 // denominator < 1 && length(partial_solution) < rank - 1
            total_sum += 1 // denominator
            push!(partial_solution, denominator)
        end
    end

    proposed_bound = reduce(lcm, partial_solution)
    return proposed_bound
end

function sum_list(denom_list::Vector{Int})::Float64
    running_sum = 0 // 1
    for denom in denom_list
        running_sum += 1 // denom
    end
    return running_sum
end

function modified_sum_list(denom_list::Vector{BigInt})::Rational{BigInt}
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

function get_solutions(possible_denoms::Vector{Any}, rank::Int, running_sol::Vector{BigInt}, list_of_solutions::Vector{Any})
    running_sum = modified_sum_list(running_sol)  # Returns Rational{BigInt}
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
        #Uncommment this if you want to look at a small portion of the solutions
        # -- honestly, you can modify this however you want and can look for specific denominator or solution also
        #=if length(list_of_solutions) < 100
            println(running_sol[end])
        end=#
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

    #CHANGE THESE TO USE DIFFERENT PRIMES AND RANKS
    p = 2
    q = 5
    rank = 11
    highest_power = 25
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
  

    get_solutions(possible_denoms, rank, BigInt[], list_of_solutions)

    
    print("There are ", length(list_of_solutions), " solutions\n")
    #for solution in list_of_solutions
        #print( solution, "\n")
   # end

   max_denominator = 0
    for solution in list_of_solutions
       for item in solution
            if item > max_denominator
                max_denominator = item
            end
        end
    end
    print("Max denominator is: ", max_denominator, "\n")
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

    open("solutions.txt", "w") do file
        for solution in list_of_solutions
            write(file, join(solution, ", "), "\n")  # Write each solution as a comma-separated line
        end
    end
    println("Solutions written to solutions.txt")

    println("The maximum denominator that appears is: ", max_denominator, " which is 2^", exp_p, " * ", q, "^", exp_q)
end 

main()
