
function CG(A,b,ϵ,maxit)
    #Lets begin with the vector zero as starting point
    x=zeros(n)
    #The starting gradient
    g=A*x-b
    #The starting descent direction
    d=-g
    #iterations
    it=0
    #starting the residuals
    res=[norm(b-A*x)]
    # we begin the while loop
    while it<maxit && norm(g)>ϵ
        #We set everything as the CG method
        α=(norm(g)^2)/(transpose(d)*A*d)[1]
        x+=α*d
        β=(norm(g+α*A*d)^2)/(norm(g)^2)
        g=g+α*A*d
        d=-g+β*d
        #Append 
        push!(res,norm(b-A*x))
        it+=1
    end
    #We return the outputs
    Any[x,it,res]
end

#Define a constructor    
function Hilbert(n)
    #Define a zero matrix to fill it with the values
    H=zeros(n,n)
    for i in 1:n
        for j in 1:n
            H[i,j]=1/(i+j-1)
        end
    end
    H
end

ϵ=1e-5
maxit=10000

n=5
A=Hilbert(n);
# b vector with ones
b=ones(n);

CG_result=CG(A,b,ϵ,maxit)

CG_result[2]

n=8
A=Hilbert(n);
# b vector with ones
b=ones(n);

CG_result=CG(A,b,ϵ,maxit)

CG_result[2]

n=12
A=Hilbert(n);
# b vector with ones
b=ones(n);

CG_result=CG(A,b,ϵ,maxit)

CG_result[2]

n=20
A=Hilbert(n);
# b vector with ones
b=ones(n);

CG_result=CG(A,b,ϵ,maxit)

CG_result[2]
