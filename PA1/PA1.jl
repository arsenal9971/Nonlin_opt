
function armijo(f,∇f,x,σ,d,β)
    #Lets obtain l until armijo rule holds
    l=0
    α=β^l
    while f(x+α*d)>f(x)+σ*α*dot(∇f(x),d)
        α=α*β
    end
    α
end

function steepestdec(f,∇f,x0,ϵ,maxit)
    #We create an array of the points in the iterations, with first entry, the starting point
    x=x0
    y=x
    i=1 
    #Start the iterations
    while i<=maxit && norm(∇f(x)) > ϵ
        #negative gradient descent direction
        d=-∇f(x)
        y=[y,x+armijo(f,∇f,x,σ,d,β)*d]
        x=x+armijo(f,∇f,x,σ,d,β)*d
        i+=1
    end
    #reshape the big bector
    reshape(y,2,int(length(y)/2))
end

f(x)=100*(x[2]-x[1]^2)^2+(1-x[1])^2

∇f(x)=[-400*(x[2]-x[1]^2)*x[1]-2*(1-x[1]),200*(x[2]-x[1]^2)]

ϵ=1e-3
maxit=10000
σ=1e-3
β=1/2

x0=[1,-0.5]
y1=steepestdec(f,∇f,x0,ϵ,maxit)

x0=[-1.2,1]
y2=steepestdec(f,∇f,x0,ϵ,maxit)

using PyPlot
using Distributions

#Define the points to plot
N=100
X=linspace(-2,2,N)
Y=linspace(-2,2,N)
Xgrid=repmat(X',N,1)
Ygrid=repmat(Y,1,N)
Z=zeros(N,N)
for i in 1:N
    for j in 1:N
        Z[i:i,j:j] = f([X[i],Y[j]])
    end
end

X1=transpose(y1[1,:])
Y1=transpose(y1[2,:])
M=length(X1)
Z1=[f([X1[i],Y1[i]]) for i in 1:M];

plot_wireframe(Xgrid,Ygrid,Z)
plot(X1, Y1, Z1, "ro")
plot(X1, Y1, Z1, "r")

#Define the points to plot
N=100
X=linspace(-2,2,N)
Y=linspace(-2,2,N)
Xgrid=repmat(X',N,1)
Ygrid=repmat(Y,1,N)
Z=zeros(N,N)
for i in 1:N
    for j in 1:N
        Z[i:i,j:j] = f([X[i],Y[j]])
    end
end

X2=transpose(y2[1,:])
Y2=transpose(y2[2,:])
M=length(X2)
Z2=[f([X2[i],Y2[i]]) for i in 1:M];

plot_wireframe(Xgrid,Ygrid,Z)
plot(X1, Y2, Z2, "ro")
plot(X1, Y2, Z2, "r")
