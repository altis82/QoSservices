using JuMP
using GLPKMathProgInterface
using Ipopt
#using ECOS
using SCS

#Number of nodes
N=5
#Number of flow
F=5
#link cost matrix
E=[[0,      1,   2,-1,-1],
   [1,      0,   2,3   ,-1],
   [2,      1,   0,-1,1],
   [-1,   3,-1,0   ,2],
   [-1,-1,   1,   2,0]]
 #finding paths from source to destination
paths=[]
flag=0
label=zeros(N,N)

numberof_path=0
list_path=[]
#######################################
function find_paths(source_node, destination_node)
    flag=0
    visited_node=[]
    for i=1:N
        if E[source_node][i]>0 && i==destination_node

            temp=[]
            for j=1:length(path)
                append!(temp,path[j])
            end
            append!(temp,i)
            append!(list_path,temp)

            #pop!(path)
            println("path found:",temp)

        elseif E[source_node][i]>0  &&in(i,path)==false

            append!(path,i)
            #println("source",source_node,"index",i)
            find_paths(i,destination_node)

            pop!(path)
        end


    end

end
######################

path=[]

append!(path,1)
#find_paths(1,5)
#####################

#1==>2==>4==>5
#|  |      |
#|=3=======

#bandwidth matrix
B=[[0,      10,   20,0,0],
   [10,      0,   20,30   ,0],
   [20,      10,   0,0,10],
   [10,   30,0,0   ,20],
   [10,10,   10,   20,0]]
#assume that all the flows can go by 4 paths
P=4
path_cost  =[6,4,9,3]
path_delay =[1,2,2,3]
path_bw    =[20,10,20,10] #calculate based on the minimum BW of the link on that path
packet_lost=[0.1,0.2,0.1,0.3]
#flow requirement
Flow_bw    =[1,3,3,2,2]
Flow_delay =[4,6,5,7,9]

#
#Initialize previous variable
rate=[10.0 0.0 0.0 0.0; 0.0 10.0 0.0 0.0; 10.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 10.0 0.0 0.0 0.0]
#
#dual variables
lambda=[0.1,0.1,0.1,0.1,0.1]
beta=0.1*ones(P)

function Optimize_flow_f(rate_f,flow_index) #a vectoc 1xp
    delta=1000
    myModel = Model(solver=GLPKSolverLP())
    #myModel= Model(solver=ECOSSolver())
    #myModel= Model(solver=SCSSolver())
    #define activation interface of devices

    @variable(myModel,r[1:P]>=0)
    @variable(myModel, temp_r[1:P]>=0)
    for p=1:P
        @constraint(myModel, temp_r[p]==r[p]*1/(rate_f[p]+delta))
        #relaxation variable should be less than 1
        @constraint(myModel, temp_r[p]<=1)
        #delay constraint
        @constraint(myModel,path_delay[p]*(temp_r[p])<=Flow_delay[flow_index])
        #requirement throughput

    end
    @constraint(myModel, sum(r[p] for p=1:P) >=Flow_bw[flow_index])

    #calculate the cost
    term1=0
    #first term
    for p=1:P
        term1=term1+path_cost[p]*temp_r[p]+ packet_lost[p]*r[p]
    end

    term2=lambda[flow_index]*(sum(temp_r[p] for p=1:P)-1)

    term3=0
    for p=1:P
        term3=term3+beta[p]*(sum(r[j] for j=1:P))
    end
    #cost function

    obj=term1+term2+term3
    @objective(myModel, Min, obj) # Sets the objective to be minimized. For maximization use Max
    #print(myModel)
    status = solve(myModel) # solves the model
    r_value=getValue(r)
    print(r_value)
    println(status)
    if status=="Optimal"

        println("The optimization problem to be solved is:")
        println("==================================")
        #
        r_value=getValue(r)
        print(r_value)
        return r_value
    else
        print("unsolved")
        return rate_f
    end
end


sigma=0.5
gamma=0.6
delta=100
new_rate=zeros(F,P)

for k=1:100
    #optimize for each flow
    temp=zeros(F,P)
    for i=1:F

        new_rate[i,:]=Optimize_flow_f(rate[i,:],i)
        #print(new_rate)
        for j=1:P
            temp[i,j]=rate[i,j]/(new_rate[i,j]+delta)
        end
        rate=new_rate

    end
    #dual update
    for i=1:F
        lambda[i]=lambda[i]+sigma*(sum(temp[i,p] for p=1:P)-1)
    end
    for p=1:P
        tmp_sum=0
        for i=1:F
            for j=1:P
                tmp_sum=tmp_sum+rate[i,j]
            end
        end
        beta[p]=beta[p]+gamma*(tmp_sum-path_bw[p])
    end
end
