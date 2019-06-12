# Ntotal=10000
# n = 30
# N = 100000
# NN = pushfirst!(collect(1:0.25:6),0)
# T = 5
# Hn = 3
# ngates = 150
# set_zero_subnormals(true)

# (gates,indices) = randomgates(ngates,T,n)
# H = sample(1:n,Hn,replace=false)
    # CSV.write("results.csv",convert(DataFrame,[n exct]);append = true)
# k = 1

for j = 1:3
    (gates,indices) = randomgates(ngates,T,n)
    H = sample(1:n,Hn,replace=false)
    global C = 0
    start = time()
    M0 = 0
    for k = 1:length(NN)
        M1=round(Int,10^NN[k])
        Threads.@threads for i = M0:M1-1
            states = samplingzero(n,N)
            #println(exact(gates,indices,n))
            #expv = zeros(Float64,N)
            c = ones(Int64,N)
            global C += sum(evaluatezero.(propagate.(states,c,Ref(indices),Ref(gates),Ref(H))))
            # expvout = C/(N)
            # CSV.write("firstresults.csv",convert(DataFrame,[1 expvout]);append = true)
            #delta = expvout - exct
        end
        t = time()-start
        expvout = C/M1*10^-5
        CSV.write("Run$j.csv",convert(DataFrame,[M1*10^5 t expvout]);append = true)
        M0 = M1
        #CSV.write("results.csv",convert(DataFrame,[1 expvout]);append = true)
    end
end
