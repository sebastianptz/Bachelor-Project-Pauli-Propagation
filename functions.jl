using LinearAlgebra
# using Combinatorics
using StatsBase
# using Plots
using CSV
using DataFrames
# using TimerOutputs
# const to = TimerOutput()
#include("Sampling.jl")
# function sampling(rho,n)
#     state = Array{Int64}(undef,n)
#     D = 1
#     for i = 1:n
#         state[i] = sample(1:4, Weights(abs.(rho[i])))
#         D = D * sum(abs.(rho[i])) #nicht das D aus dem paper D(paper)=D*1/2^n
#         global sgn = sign(rho[i][state[i]])
#     end
#     c = 1/2^n * D * sgn
#     return (state,c)
# end
################
function samplingzero(n,N)
    states = zeros(Int8,N)
    states = bits.(states,Ref(n))
    # for i = 1:N
    #     push!(states, 3 .* bits(n) .+ 1)
    # end
    return states
end

function bits(state,n)
    return 3 .* digits(rand(0:2^n-1),base=2,pad=n) .+ 1
end
#include("Propagation.jl")
function propagate(state,c,indices,gates,H)
    #1=Hadamard,2=Phase,3=CNOT
    clifford = [[1;4;3;2],[1;3;2;4],[[1;1],[1;2],[4;3],[4;4],[2;2],[2;1],[3;4],[3;3],[3;2],[3;1],[2;4],[2;3],[4;1],[4;2],[1;3],[1;4]]]
    clsign = [[1;1;-1;1],[1;1;-1;1],[1;1;1;1;1;1;1;-1;1;1;-1;1;1;1;1;1]]
    #sigma = [[[1 0]; [0 1]], [[0 1]; [1 0]], [[0 -im]; [im 0]], [[1 0]; [0 -1]]]
    ############
    # H = [[1;0;0;0],[0;0;0;1],[0;0;-1;0],[0;1;0;0]]
    # S = [[1;0;0;0],[0;0;1;0],[0;-1;0;0],[0;0;0;1]]
    # T = [[1;0;0;0],[0;1/sqrt(2);1/sqrt(2);0],[0;-1/sqrt(2);1/sqrt(2);0],[0;0;0;1]]
    #CNOT = [[1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0],[0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;0],[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1],[0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0],[0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0],[0;0;0;0;0;0;0;0;0;0;-1;0;0;0;0;0],[0;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;-1;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;0;0;0;0;1;0;0;0],[0;0;0;0;0;0;0;0;0;0;0;0;0;1;0;0],[0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0],[0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0]]
    #cnot = [1;2;15;16;6;5;12;-11;10;9;-8;7;13;14;3;4]
    #[1;2;15;16;6;5;12;-11;10;9;-8;7;13;14;3;4]
    #[[1;1],[1;2],[4;3],[4;4],[2;2],[2;1],[3;4],[?],[3;2],[3;1],[?],[2;3],[4;1],[4;2],[1;3],[1;4]]
    ############
    for j = 1:length(gates)
        if gates[j] <= 2
            state[indices[j]]=clifford[gates[j]][state[indices[j]]]
            c = c * clsign[gates[j]][state[indices[j]]]
        elseif gates[j] == 3
            cnotstate = ((state[indices[j][1]]-1)*4+(state[indices[j][2]]-1))+1
            c = c * clsign[3][cnotstate]
            state[indices[j][1]]=clifford[3][cnotstate][1]
            state[indices[j][2]]=clifford[3][cnotstate][2]
        else
            if state[indices[j]] == 3
                state[indices[j]] = rand(2:3)
                c = c * sqrt(2)
            elseif state[indices[j]] == 2
                rndm = rand(2:3)
                if rndm == 3
                    state[indices[j]] = 3
                    c = -c
                end
                c = c * sqrt(2)
            end
        end
    end
    clsign[2] = [1;-1;1;1]
    for j = length(gates):-1:1
        if gates[j] <= 2
            state[indices[j]]=clifford[gates[j]][state[indices[j]]]
            c = c * clsign[gates[j]][state[indices[j]]]
        elseif gates[j] == 3
            cnotstate = ((state[indices[j][1]]-1)*4+(state[indices[j][2]]-1))+1
            c = c * clsign[3][cnotstate]
            state[indices[j][1]]=clifford[3][cnotstate][1]
            state[indices[j][2]]=clifford[3][cnotstate][2]
        else
            if state[indices[j]] == 2
                state[indices[j]] = rand(2:3)
                c = c * sqrt(2)
            elseif state[indices[j]] == 3
                rndm = rand(2:3)
                if rndm == 2
                    state[indices[j]] = 2
                    c = -c
                end
                c = c * sqrt(2)
            end
        end
    end
    for i = 1:length(H)
        state[H[i]] = clifford[1][state[H[i]]]
        c = c * clsign[1][state[H[i]]]
    end
    return (state,c)
end
####################
function randomgates(ngates,T,n)
    indices = []
    gates = Int64[]
    for i = 1:(ngates-T)
        P = [0.25, 0.25, 0.5]
        push!(gates,sample(1:3, Weights(P)))
        if gates[i] != 3
            push!(indices,rand(1:n))
        else
            ix = rand(1:n)
            a = collect(1:n)
            deleteat!(a,ix)
            push!(indices,[ix,rand(a)])
        end
    end
    for i = 1:T
        r = rand(1:length(gates)+1)
        insert!(gates,r,4)
        insert!(indices,r,rand(1:n))
    end
    return(gates,indices)
end
##################
# function notsorandomgates(ngates,T,n)
#     indices = []
#     gates = []
#     for i = 1:ngates
#         P = [0.25 * (1-p), 0.25 * (1-p), 0.5 * (1-p), p]
#         push!(gates,sample(1:4, Weights(P)))
#         if gates[i] != 3
#             push!(indices,rand(1:n))
#         else
#             r = Int64(rand(1:n/2)*2)
#             push!(indices,[r-1,r])
#         end
#     end
#     return(gates,indices)
# end
##################
# function exact(gates,indices,n)
#     g = [1/sqrt(2) * [[1 1];[1 -1]],[[1 0];[0 im]],[[1 0 0 0];[0 1 0 0];[0 0 0 1];[0 0 1 0]],[[1 0];[0 1/sqrt(2)*(1+1im)]],[[1 0 0 0];[0 0 0 1];[0 0 1 0];[0 1 0 0]]]
#     I = [[1 0];[0 1]]
#     null = []
#     for i = 1:Int64(n/2)
#         push!(null,[1;0;0;0])
#     end
#     for i = 1:length(gates)
#         if gates[i] != 3
#             if isodd(indices[i])==true
#                 null[Int64((indices[i]+1)/2)]=kron(g[gates[i]],I)*null[Int64((indices[i]+1)/2)]
#             else
#                 null[Int64(indices[i]/2)]=kron(I,g[gates[i]])*null[Int64(indices[i]/2)]
#             end
#         else
#             if isodd(indices[i][1])==true
#                 null[Int64(indices[i][2]/2)]=g[gates[i]]*null[Int64(indices[i][2]/2)]
#             else
#                 null[Int64(indices[i][1]/2)]=g[5]*null[Int64(indices[i][1]/2)]
#             end
#         end
#     end
#     for i = 1:length(null)
#         null[i]=null[i][1]
#     end
#     expv = (abs(prod(null)))^2
#     return expv
# end
##################
function samplingzero1(n,N)
    states = Array{Int8}[]
    tindex = Array{Int8}[]
    # a = 3 .* digits(rand(0:2^n*N-1),base=2,pad=n*N) .+1
    a = sample([1,4],n*N)
    b = sample([2,3],2*T*N)
    for i = 1:N
        push!(states, a[n*i-(n-1):n*i])
        push!(tindex, b[2*T*i-(2*T-1):2*T*i])
    end
    return states,tindex
end
#include("Propagation.jl")
function propagate1(state,c,tindex,indices,gates,H)
    #1=Hadamard,2=Phase,3=CNOT
    clifford = [[1;4;3;2],[1;3;2;4],[[1;1],[1;2],[4;3],[4;4],[2;2],[2;1],[3;4],[3;3],[3;2],[3;1],[2;4],[2;3],[4;1],[4;2],[1;3],[1;4]]]
    clsign = [[1;1;-1;1],[1;1;-1;1],[1;1;1;1;1;1;1;-1;1;1;-1;1;1;1;1;1]]
    #sigma = [[[1 0]; [0 1]], [[0 1]; [1 0]], [[0 -im]; [im 0]], [[1 0]; [0 -1]]]
    ############
    # H = [[1;0;0;0],[0;0;0;1],[0;0;-1;0],[0;1;0;0]]
    # S = [[1;0;0;0],[0;0;1;0],[0;-1;0;0],[0;0;0;1]]
    # T = [[1;0;0;0],[0;1/sqrt(2);1/sqrt(2);0],[0;-1/sqrt(2);1/sqrt(2);0],[0;0;0;1]]
    #CNOT = [[1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0],[0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;0],[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1],[0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0],[0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0],[0;0;0;0;0;0;0;0;0;0;-1;0;0;0;0;0],[0;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;-1;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;0;0;0;0;1;0;0;0],[0;0;0;0;0;0;0;0;0;0;0;0;0;1;0;0],[0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0],[0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0]]
    #cnot = [1;2;15;16;6;5;12;-11;10;9;-8;7;13;14;3;4]
    #[1;2;15;16;6;5;12;-11;10;9;-8;7;13;14;3;4]
    #[[1;1],[1;2],[4;3],[4;4],[2;2],[2;1],[3;4],[?],[3;2],[3;1],[?],[2;3],[4;1],[4;2],[1;3],[1;4]]
    ############
    k = 1
    for j = 1:length(gates)
        if gates[j] <= 2
            state[indices[j]] = clifford[gates[j]][state[indices[j]]]
            c = c * clsign[gates[j]][state[indices[j]]]
        elseif gates[j] == 3
            cnotstate = ((state[indices[j][1]]-1)*4+(state[indices[j][2]]-1))+1
            c = c * clsign[3][cnotstate]
            state[indices[j][1]] = clifford[3][cnotstate][1]
            state[indices[j][2]] = clifford[3][cnotstate][2]
        else
            if state[indices[j]] == 3
                state[indices[j]] = tindex[k]
                c = c * sqrt(2)
                k += 1
            elseif state[indices[j]] == 2
                rndm = tindex[k]
                k += 1
                if rndm == 3
                    state[indices[j]] = 3
                    c = -c
                end
                c = c * sqrt(2)
            end
        end
    end
    clsign[2] = [1;-1;1;1]
    for j = length(gates):-1:1
        if gates[j] <= 2
            state[indices[j]] = clifford[gates[j]][state[indices[j]]]
            c = c * clsign[gates[j]][state[indices[j]]]
        elseif gates[j] == 3
            cnotstate = ((state[indices[j][1]]-1)*4+(state[indices[j][2]]-1))+1
            c = c * clsign[3][cnotstate]
            state[indices[j][1]] = clifford[3][cnotstate][1]
            state[indices[j][2]] = clifford[3][cnotstate][2]
        else
            if state[indices[j]] == 2
                state[indices[j]] = tindex[k]
                c = c * sqrt(2)
                k += 1
            elseif state[indices[j]] == 3
                rndm = tindex[k]
                k += 1
                if rndm == 2
                    state[indices[j]] = 2
                    c = -c
                end
                c = c * sqrt(2)
            end
        end
    end
    x = sample(1:n,H,replace=false)
    for i = 1:H
        state[x[i]] = clifford[1][state[x[i]]]
        c = c * clsign[1][state[x[i]]]
    end
    return (state,c)
end

function evaluatezero(statec)
    for i = 1:length(statec[1])
        if statec[1][i] in [2,3]
            return 0
            break
        end
    end
    return statec[2]
end
##################
# function randomrho(n)
#     rho = Array{Array{Float64,1},1}(undef, n)
#     for i = 1:n
#         rho[i] = rand(Float64,3)
#         for j = 1:3
#             rho[i][j] = rho[i][j] * rand([-1,1])
#         end
#         pushfirst!(rho[i],1)
#     end
#     return(rho)
# end
##################
n = 30
N = 100000
NN = collect(1:0.3:4.9)
T = 5
Hn = 3
ngates = 150
set_zero_subnormals(true)
(gates,indices) = randomgates(ngates,T,n)
H = sample(1:n,Hn,replace=false)
C = 0
Threads.@threads for q = 1:1
    states = samplingzero(n,N)
    c = ones(Int64,N)
    global C += sum(evaluatezero.(propagate.(states,c,Ref(indices),Ref(gates),Ref(H))))
end
C = 0
