using JuMP
using GLPKMathProgInterface
using Ipopt
using ECOS
using SCS
#using PyPlot

delta=0.2
###
cost_bins=[1,2,3]
cap_bins=[5,6,7]
items=[2,3,2,4,5]
number_bins=3
number_items=5
lambda=0.3*ones(number_bins)
x_temp=zeros(number_items,number_bins)
total_cost=[]
for t=1:100
    println("allocation",x_temp)
    for k=1:number_items
        myModel1 = Model(solver=GLPKSolverLP())
        #myModel1= Model(solver=ECOSSolver())
        #myModel1= Model(solver=SCSSolver())
        #define activation interface of devices
        @variable(myModel1,x[1:number_bins]>=0)
        @constraint(myModel1,sum(x[i] for i=1:number_bins)==1)
        for p=1:number_bins
            #print("cost",sum(x_temp[:,p])-x_temp[k,p]-cap_bins[p])
            @constraint(myModel1,x[p]*items[k]+sum(x_temp[i,p]*items[i] for i=1:number_items)-x_temp[k,p]*items[k]<=cap_bins[p])
        end
        obj=sum(x[i]*cost_bins[i]*items[k] for i=1:number_bins) +sum(x[p]*lambda[p]*items[k] for p=1:number_bins)
        @objective(myModel1, Min, obj) # Sets the objective to be minimized. For maximization use Max
        #print(myModel1)
        status = solve(myModel1) # solves the model
        x_temp[k,:]=getValue(x)
        #println(x_temp)
        #println(getObjectiveValue(myModel1))
    end

    for j=1:number_bins
         lambda[j]=lambda[j]+delta*(sum(x_temp[i,j]*items[i] for i=1:number_items)-cap_bins[j])

        if lambda[j]<0
            lambda[j]=abs(delta*(sum(x_temp[i,j]*items[i] for i=1:number_items)-cap_bins[j]))
            #lambda[j]=0
            #global delta=delta/2

        end
        #println(lambda[j])
    end
    cost=0
    for i=1:number_items
        cost=cost+sum(x_temp[i,j]*cost_bins[j] for j=1:number_bins)
    end
    append!(total_cost,cost)


end
#plot(total_cost)
println(total_cost)
