""" Wired.jl - Julia Benchmark 
    (c) 2024 ryan@freestatelabs.com 
"""

using LinearAlgebra: norm, dot
const mu0 = 4pi * (1e-7)

# Create the benchmark problem with `N` nodes and `M` sources in a row-major format
#  (slower). `nodes` is 3xN, `sources` is 7xM. 
function createproblem_rm(N::Integer, M::Integer; I=5000, Px=3, Py=3)
    nodes = zeros(3,N) 
    sources = zeros(7,M) 

    zstart = 0; zend = 1
    dz = (zend - zstart) / (N-1)
    z = zstart

    for i = 1:N
        nodes[:,i] = [Px, Py, z]
        z += dz 
    end

    zstart = -1000; zend = 1000 
    dz = (zend - zstart) / M
    z = zstart 
    for i = 1:M 
        sources[:,i] = [0, 0, z, 0, 0, z+dz, I] 
        z += dz
    end

    return nodes, sources
end


# Create the benchmark problem with `N` nodes and `M` sources in a column-major 
#   format (faster). Nodes is Nx3, sources is 7xM.
function createproblem(N::Integer, M::Integer; I=5000, Px=3, Py=3)
    nodes = zeros(N,3) 
    sources = zeros(7,M) 

    zstart = 0; zend = 1
    dz = (zend - zstart) / (N-1)
    z = zstart

    for i = 1:N
        nodes[i,:] = [Px, Py, z]
        z += dz 
    end

    zstart = -1000; zend = 1000 
    dz = (zend - zstart) / M
    z = zstart 
    for i = 1:M 
        sources[:,i] = [0, 0, z, 0, 0, z+dz, I] 
        z += dz
    end

    return nodes, sources
end


# In-place cross-product to reduce memory allocations, c = a x b 
# This version is defined for 3-length vectors
# See: https://discourse.julialang.org/t/excessive-memory-allocation/108841/4 
function cross!(c::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64})
    if !(length(a) == length(b) == 3)
        throw(DimensionMismatch("Vectors must be of length=3!"))
    end

    c .= (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
end


# In-place cross-product for row-major format; cross the rows of A with b and place in C
@views function cross!(C::AbstractArray{Float64}, A::AbstractArray{Float64}, b::Vector{Float64})

    C[1,:] .= A[2,:] .* b[3] .- A[3,:] .* b[2]
    C[2,:] .= A[3,:] .* b[1] .- A[1,:] .* b[3]
    C[3,:] .= A[1,:] .* b[2] .- A[2,:] .* b[1]

end


# In-place cross-product for column-major format; cross the columns of A with b and place in C
@views function crosscols!(C::AbstractArray{Float64}, A::AbstractArray{Float64}, b::Vector{Float64})

    C[:,1] .= A[:,2] .* b[3] .- A[:,3] .* b[2]
    C[:,2] .= A[:,3] .* b[1] .- A[:,1] .* b[3]
    C[:,3] .= A[:,1] .* b[2] .- A[:,2] .* b[1]

end


# In-place dot-product of two arrays; note that a and b must both be 3xN
@views function dotcols!(c::Vector{Float64}, A::AbstractArray{Float64}, B::AbstractArray{Float64})

    c .= (A[1,:] .* B[1,:] .+ A[2,:] .* B[2,:] .+ A[3,:] .* B[3,:])
end


# In-place dot-product of the rows of A by vector b; place in C
@views function dotrows!(C::Vector{Float64}, A::AbstractArray{Float64}, b::Vector{Float64})

    C .= (A[:,1] .* b[1] .+ A[:,2] .* b[2] .+ A[:,3] .* b[3])
end


# In-place vector norm of the 3-length columns of an array 
#   This is for row-major format (slower)
@views function normcols!(c::Vector{Float64}, A::AbstractArray{Float64})

    c .= (sqrt.(A[1,:].^2 .+ A[2,:].^2 .+ A[3,:].^2))

end


# In-place vector norm of the 3-length rows of an array 
#   This is for column-major format (faster)
@views function normrows!(c::Vector{Float64}, A::AbstractArray{Float64})

    c .= (sqrt.(A[:,1].^2 .+ A[:,2].^2 .+ A[:,3].^2))

end


# Multiply the rows of a matrix `A` by a vector b 
@views function mult!(A::AbstractArray{Float64}, b::Vector{Float64})

    A[1,:] .*= b 
    A[2,:] .*= b 
    A[3,:] .*= b

end


# Multiply the rows of a matrix `A` by a vector b 
@views function multcols!(A::AbstractArray{Float64}, b::Vector{Float64})

    A[:,1] .*= b 
    A[:,2] .*= b 
    A[:,3] .*= b

end


# Solve the Biot-Savart problem with optimization level 1: simple loops
#   This is for row-major implementation; it turns out this is about as fast 
#   as the column-major implementation
@views function solve_rm(nodes::AbstractArray{Float64}, sources::AbstractArray{Float64})

    B = zeros(size(nodes))
    a = Vector{Float64}(undef,3)
    b = similar(a)
    c = similar(a) 
    d = 0.0
    cxa = similar(a)

    for j in axes(sources)[2]
        d = mu0 * sources[7,j] / (4pi)
        a .= sources[4:6,j] .- sources[1:3,j] 
        for i in axes(nodes)[2]
            b .= sources[1:3,j] .- nodes[:,i] 
            c .= sources[4:6,j] .- nodes[:,i]
            cross!(cxa, c, a)

            B[:,i] .+= cxa .* d .* (1/(norm(cxa)^2)) .* (dot(a,c)/norm(c) .- dot(a,b)/norm(b))
        end
    end

    return B
end


# Solve the Biot-Savart problem with optimization level 1: simple loops
@views function solve(nodes::AbstractArray{Float64}, sources::AbstractArray{Float64})

    B = zeros(size(nodes))
    a = Vector{Float64}(undef,3)
    b = similar(a)
    c = similar(a) 
    d = 0.0
    cxa = similar(a)

    for j in axes(sources)[2]
        d = mu0 * sources[7,j] / (4pi)
        a .= sources[4:6,j] .- sources[1:3,j] 
        for i in axes(nodes)[1]
            b .= sources[1:3,j] .- nodes[i,:] 
            c .= sources[4:6,j] .- nodes[i,:]
            cross!(cxa, c, a)

            B[i,:] .+= cxa .* d .* (1/(norm(cxa)^2)) .* (dot(a,c)/norm(c) .- dot(a,b)/norm(b))
        end
    end

    return B
end


# Solve the Biot-Savart problem with optimization level 1: vectorization
# Replace the inner loop with matrix/vector operations to allow Julia to make 
#   low-level optimizations not otherwise visible to the programmer
# This is the original row-major implementation (slower)
@views function solve2_rm(nodes::AbstractArray{Float64}, sources::AbstractArray{Float64})

    Nn = size(nodes)[2]
    B = zeros(3, Nn)
    a = zeros(3) 
    b = zeros(3, Nn) 
    c = zeros(3,Nn) 
    cxa = zeros(3,Nn) 
    norm_cxa = zeros(Nn) 
    dot_ac = zeros(Nn) 
    norm_c = zeros(Nn) 
    dot_ab = zeros(Nn) 
    norm_b = zeros(Nn) 

    d = 0.0

    for j in axes(sources)[2]
        d = mu0 * sources[7,j] / (4pi)
        a .= sources[4:6,j] .- sources[1:3,j] 
        b .= sources[1:3,j] .- nodes 
        c .= sources[4:6,j] .- nodes

        cross!(cxa, c, a)
        normcols!(norm_cxa, cxa)
        dotcols!(dot_ac, a, c) 
        dotcols!(dot_ab, a, b) 
        normcols!(norm_c, c) 
        normcols!(norm_b, b)

        mult!(cxa, (d .* (norm_cxa.^(-2)) .* (dot_ac./norm_c .- dot_ab./norm_b))) 
        B .+= cxa
    end

    return B
end


# Solve the Biot-Savart problem with optimization level 2: vectorization
# Replace the inner loop with matrix/vector operations to allow Julia to make 
#   low-level optimizations not otherwise visible to the programmer
# This is the column-major version (faster)
@views function solve2(nodes::AbstractArray{Float64}, sources::AbstractArray{Float64})

    Nn = size(nodes)[1]
    B = zeros(Nn,3)
    a = zeros(3) 
    b = zeros(Nn,3) 
    c = zeros(Nn,3) 
    cxa = zeros(Nn,3) 
    norm_cxa = zeros(Nn) 
    dot_ac = zeros(Nn) 
    norm_c = zeros(Nn) 
    dot_ab = zeros(Nn) 
    norm_b = zeros(Nn) 

    d = 0.0

    for j in axes(sources)[2]

        # The use of "dot assignment" i.e. `.=` prevents re-allocation of memory 
        #  within the loop
        d = mu0 * sources[7,j] / (4pi)
        a .= sources[4:6,j] .- sources[1:3,j] 
        b .= sources[1:3,j]' .- nodes 
        c .= sources[4:6,j]' .- nodes

        # In-place operators also prevent re-allocation of memory
        crosscols!(cxa, c, a)
        normrows!(norm_cxa, cxa)
        dotrows!(dot_ac, c, a) 
        dotrows!(dot_ab, b, a) 
        normrows!(norm_c, c) 
        normrows!(norm_b, b)

        multcols!(cxa, (d .* (norm_cxa.^(-2)) .* (dot_ac./norm_c .- dot_ab./norm_b))) 
        B .+= cxa
    end

    return B
end


# Split a problem up for multi-threaded operation 
# For thread number `it`, total number of threads `Nt`, and total number of 
#   tasks `N`, determine the start/stop index for that particular thread number
function threadindices(it::Integer, Nt::Integer, N::Integer) 

    Nperthread = div(N, Nt)
    remainder = rem(N, Nt) 

    if it == 1
        i1 = 1 
        i2 = i1 + Nperthread - 1 
    elseif it == Nt 
        i2 = N 
        i1 = i2 - Nperthread - remainder + 1 
    else
        i1 = (it-1)*Nperthread + 1
        i2 = i1 + Nperthread - 1 
    end 

    return i1:i2
end


# Solve the Biot-Savart problem with optimization level 3: multi-threading
# Default is to use all threads available to Julia
@views function solve3(nodes::AbstractArray{Float64}, sources::AbstractArray{Float64}; 
                        Nt::Integer = 0)

    Ns = size(sources)[2]
    if Nt == 0 
        Nt = Threads.nthreads()
    end

    # Spawn a new task for each thread by splitting up the source array
    tasks = Vector{Task}(undef, Nt)
    for it = 1:Nt 
        @views tasks[it] = Threads.@spawn solve2(nodes, sources[:,threadindices(it, Nt, Ns)])
    end
    
    # Get the result from each calculation and add it to the output array 
    B = zeros(size(nodes))
    for it = 1:Nt 
        B .+= fetch(tasks[it]) 
    end 

    return B
end