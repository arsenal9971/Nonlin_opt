
# CG method
function CG(A,b,ϵ,maxit)
    #Lets begin with the vector zero as starting point
    x=zeros(size(A)[1])
    #The starting gradient
    g=A*x-b
    #The starting descent direction
    d=-g
    #iterations
    it=0
    # we begin the while loop
    while it<maxit && norm(g)>ϵ
        #We set everything as the CG method
        # Using the modification in the Remark 6.2
        α=max(0,(norm(g)^2)/(transpose(d)*A*d)[1])
        x+=α*d
        β=max(0,transpose(α*A*d)*(g+α*A*d)/norm(g)^2)[1]
        g=g+α*A*d
        d=-g+β*d
        it+=1
    end
    #We return the outputs
    x
end

# Newton point
function NewtonPoint(xa,f,∇f,∇2f,ϵ,maxit)
    xa+CG(∇2f(xa),-∇f(xa),ϵ,maxit)
end

function CauchyPoint(xa,f,∇f,∇2f,Δa)
    if (transpose(∇f(xa))*∇2f(xa)*∇f(xa))[1] <=0
        αa=Δa/norm(∇f(xa))
    else
        αa=min(Δa/norm(∇f(xa)),
        norm(∇f(xa))^2/((transpose(∇f(xa))*∇2f(xa)*∇f(xa))[1]))
    end
    xa-αa*∇f(xa)
end

function dogleg(xa,f,∇f,∇2f,Δ)
    # Lets compute the Cauchy and newton points in this Benchmark
    # Using as tolerance for the Newton-CG method 1e-5 and maxiteration for the
    # same 10000
    xN=NewtonPoint(xa,f,∇f,∇2f,1e-5,10000)
    xCP=CauchyPoint(xa,f,∇f,∇2f,Δ)
    # Lets check the options
    Cond1=(transpose(xN-xCP)*(xCP-xa))[1]>=0
    Cond2=(transpose(∇f(xa))*∇2f(xa)*∇f(xa))[1] >0
    # With a tolerance in the equality of 1e-15
    if (Δ-1e-15<=norm(xa-xCP) && norm(xa-xCP)<=Δ+1e-15) || (Cond2 && ~Cond1)
        xv=xCP
    elseif (norm(xa-xCP)<norm(xa-xN) && norm(xa-xN)<=Δ) && Cond1
        xv=xN
    else
        a=norm(xCP)^2+norm(xN)^2-2*dot(xCP,xN)
        b=-2*norm(xN)^2+2*dot(xCP,xN)-2*dot(xa,xCP)+2*dot(xa,xN)
        c=norm(xN)^2+norm(xa)^2-2*dot(xa,xN)-Δ^2
        t=(-b+sqrt(b^2-4*a*c))/2*a
        xv=t*xCP+(1-t)*xN
    end
    # We return dv=xv-xa
    xv-xa
end

function trustregion_step(xa,xv,f,∇f,∇2f,Δ,maxiter)
    #First lets define the parameters
    # As just three different steps are needed to do the calculations
    # we are just going to save three l's
    μ0=0.01
    μ_=0.1
    μ̄=0.9
    ω_=0.4
    ω̄=2.0
    C=1.0
    # In order to have the least memory usage possible we are going to 
    # arrays with 3 entries that will represent the steps [l-1,l,l+1]
    z=Any[xa]
    zv=Any[xv]
    Δ=[Δ]
    ared=Any[]
    pred=Any[]
    dv=Any[]
    l=1
    while z[l]==xa && l<maxiter
        # First calculate the actual reduction and the predicted reduction 
        push!(ared,f(xa)-f(zv[l]))
        push!(dv,zv[l]-xa)
        push!(pred,-(transpose(∇f(xa))*dv[l])[1]-(1/2)*(transpose(dv[l])*∇2f(xa)*dv[l])[1])
        if ared[l]/pred[l] < μ0
            push!(z,xa)
            push!(Δ,ω_*Δ[l])
            if l>=2 && Δ[l]>Δ[l-1]
                z[l+1]=z[l-1]
                Δ[l+1]=Δ[l-1]
            else
                push!(dv,dogleg(xa,f,∇f,∇2f,Δ[l+1]))
                push!(zv,xa+dv[l+1])
            end
        elseif μ0<=ared[l]/pred[l] && ared[l]/pred[l]<=μ_
            push!(z,zv[l])
            push!(Δ,ω_*Δ[l])
        elseif μ_<=ared[l]/pred[l] && ared[l]/pred[l]<=μ̄
            push!(z,zv[l])
        elseif μ̄<=ared[l]/pred[l] 
            if (Δ[l]-1e-15<=norm(dv[l]) && norm(dv[l])<=Δ[l]+1e-15) && Δ[l]<=C*norm(∇f(xa))
                push!(z,xa)
                push!(Δ,ω̄*Δ[l])
                push!(dv,dogleg(xa,f,∇f,∇2f,Δ[l+1]))
                push!(zv,xa+dv[l+1])
            else 
                push!(z,zv[l])
            end
        end
        l+=1
    end
    Any[z[length(z)],Δ[length(Δ)]]
end 

function trustregion(f,∇f,∇2f,x0,r0,ϵ,niter)
    #Initialize the values
    X=Any[x0]
    r=r0
    k=1
    while norm(∇f(X[k]))>ϵ && k<=niter
        # Calculate with the dogleg method the trial point
        xv=X[k]+dogleg(X[k],f,∇f,∇2f,r)
        # Get the new point and the new trust region radio with the trust
        # region algorithm 7.3
        step=trustregion_step(X[k],xv,f,∇f,∇2f,r,1000)
        push!(X,step[1])
        # We change to the new radio
        r=step[2]
        k=k+1
    end
    X
end

# We define the Rosenbrock function 
function f(x)
    100*(x[2]-x[1]^2)^2+(1-x[1])^2
end

# We define the gradient of the Rosenbrock function
function ∇f(x)
    [400*x[1]^3+2*x[1]-400*x[1]*x[2]-2,-200*x[1]^2+200*x[2]]
end

# Finally, we define the laplacian of the Rosenbrock function
function ∇2f(x)
    reshape([1200*x[1]^2-400*x[2]+2,200,200,-400*x[1]],2,2)
end

# Lets define the initial point x0, r0 the initial trust region radius, ϵ is the termination tolerance, and niter the
# maximal number of iteration steps
x0=[-5.0,5.0]
r0=1.0
ϵ=1.0e-5
niter=1000;

# Lets calculate the solution in this case
X=trustregion(f,∇f,∇2f,x0,r0,ϵ,niter)
