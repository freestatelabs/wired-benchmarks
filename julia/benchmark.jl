""" Wired.jl - Julia Benchmark Run File 
    (c) 2024 ryan@freestatelabs.com 
"""

include("wired.jl")
using BenchmarkTools, PrettyTables

N = (100, 1_000, 10_000) 
results = Matrix{Any}(undef, 5, 4)
header = ["Opt Level", "100x100", "1000x1000", "10000x10000"]
results[:,1] = ["Simple Loops", "Vectorization", "Parallel, 2 cores", "Parallel, 4 cores", "Parallel, 8 cores"]

for i in eachindex(N)
    print("Benchmarking with N = ")
    print(N[i]) 
    print("\n")
    nodes, sources = createproblem(N[i],N[i])
    results[1,i+1] = @belapsed solve($nodes, $sources)
    results[2,i+1] = @belapsed solve2($nodes, $sources)
    results[3,i+1] = @belapsed solve3($nodes, $sources, Nt=2)
    results[4,i+1] = @belapsed solve3($nodes, $sources, Nt=4)
    results[5,i+1] = @belapsed solve3($nodes, $sources, Nt=8)
end
println("Benchmarking complete.")
pretty_table(results, header=header)


# Display the answer to the problem using each function to ensure accuracy of 
#   the solution
nodes, sources = createproblem(100,100)
display(solve(nodes, sources))
display(solve2(nodes, sources))
display(solve3(nodes, sources))