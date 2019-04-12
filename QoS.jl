using JuMP
using GLPKMathProgInterface
using Ipopt
#using ECOS
#using SCS
#using PyPlot
#Number of nodes
N=5

#Number of flows
#F=3

#Number of flows
F=2
#link cost matrix
E=[[0,  1, 2,-1,-1],
   [1,  0, 2, 3,-1],
   [2,  2, 0,-1, 1],
   [-1, 3,-1, 0, 2],
   [-1,-1, 1, 2, 0]]



#finding paths from source to destination
paths=[]
#flag=0
label=zeros(N,N)
numberof_path=0
list_path=[]
iteration=10

function find_paths(source_node, destination_node)
    #flag=0
    visited_node=[]
    temp=[]
    k=1
    for i=1:N
        #print("\n--------------------------------------------------------------\n")
        #print("\nE[",source_node,"][",i,"]: ", E[source_node][i])
        if E[source_node][i]>0 && i==destination_node
            for j=1:length(path)
                append!(temp,path[j])
                append!(list_path,path[j])
            end
            append!(temp,i)
            append!(list_path,i)
            println("\npath found:",temp)
            #print("\nThis is list path after temp appending: ", list_path)
        elseif E[source_node][i]>0  &&in(i,path)==false
            append!(path,i)
            find_paths(i,destination_node)
            pop!(path)
        end
    end
end
source=1
destination=5
path=[]
append!(path,1)
find_paths(source,destination)
print("\nThis is list path: ", list_path)
print("This is list_path: ", list_path)
print("\n\nThis the element list_path[7]: ", list_path[7])
print("\nThis is the length: ", length(list_path))
# Finding all the links of E that are bloging to p in P
#           1  2  3   4  5
    #E =1[[0, 1, 2, -1,-1],
    #   2 [1, 0, 2,  3,-1],
    #   3 [2, 2, 0, -1, 1],
    #   4 [-1, 3,-1, 0, 2],
    #   5 [-1,-1, 1, 2, 0]]

#bandwidth matrix
B=[[0,  10, 20, 0,  0],
   [10, 0,  20, 10 ,0],
   [20, 20, 0,  0, 20],
   [0,  10, 0,  0, 20],
   [0,  0,  20, 20, 0]]

#if we kip this dataset, the results are different for lambdae
#B=[[0,  10, 20, 0,  0],
#   [10, 0,  20, 50 ,0],
#   [20, 20, 0,  0, 60],
#   [0,  50, 0,  0, 20],
#   [0,  0,  60, 20, 0]]
#assume that all the flows can go by 4 paths
P=4
path_cost  =[0.6,0.1,0.9,0.2]
path_delay =[0.5,0.1,0.2,0.3]
path_bw    =[0.1,0.1,0.2,0.1] #calculate based on the minimum BW of the link on that path
packet_lost=[0.5,0.4,0.8,0.3]

# 5 flows requirement
#Flow_bw    =[11,31,31,12,21]
#Flow_delay =[40,60,50,70,90]

# 7 flows requirement
#Flow_bw    =[11,31,31,12,21,22,13]
#Flow_delay =[40,60,50,70,90,75,80]

# 2 flows requirement
Flow_bw    =[0.07,0.05] #[11,5]
Flow_delay =[0.4,0.6]

# 3 flows requirement
#Flow_bw    =[4,8,2]
#Flow_delay =[60,20,40]

# FLOWS WITH DELAY AND BANDWIDTH CONSTRAINT ROUTING
###################################################



#
#Initialize previous variable with 2 flows
rate=zeros(F,P)

#Initialize previous variable with 3 flows
#rate=[10.0 0.0 0.0 0.0; 0.0 0.0 15.0 0.0; 0.0 0.0 0.0 25.0]

#Initialize previous variable with 7 flows
#rate=[10.0 0.0 0.0 0.0; 0.0 15.0 0.0 0.0; 10.0 0.0 0.0 0.0; 0.0 0.0 0.0 50.0; 0.0 10.0 0.0 0.0; 0.0 0.0 10.0 0.0; 10.0 0.0 0.0 0.0]

#dual variables
#lambda=0.1*ones(P)
#duallambdaprevious=zeros(P)
#duallambda=zeros(P)

#lambdae = 0.1
duallambdaprevious = 0
duallambda = 0
delta=0.01
cons=[[1,1,0.1],
        [1,2,0.1],
        [2,2,0.2],
        [2,1,0.1],
        [3,3,0.1],
        [3,4,0.2],
        [4,4,0.1],
        [1,1,0.1],
        [4,4,0.1],
        [1,4,0.1],
        [2,2,0.2],
        [2,3,0.2],
        [3,3,0.2]]

lambdae=zeros(length(cons))
scale=20
#delta=10

#objective_value=0.3
objsaving=[]
flow_objective_value=zeros(F, iteration) # 5 is the number of iteration
function Optimize_flow_f(rate_f,flow_index) #a vectoc 1xp
    println("lamdba",lambdae)
    #myModel = Model(solver=GLPKSolverLP())
    myModel = Model(solver=IpoptSolver())
    #myModel= Model(solver=ECOSSolver())
    #myModel= Model(solver=SCSSolver())
    #define activation interface of devices

    @variable(myModel,r[1:P]>=0)
    @variable(myModel,0<=temp_r[1:P]<=1)
    t=0
    for p=1:P

        @constraint(myModel, temp_r[p]<=r[p]/(rate_f[p]+delta))
        @constraint(myModel, r[p]<=path_bw[p])
        #relaxation variable should be less than 1

        #delay constraint
        @constraint(myModel, path_delay[p]*temp_r[p]<=temp_r[p]*Flow_delay[flow_index])

    end

    #requirement throughput constraint
    @constraint(myModel, sum(r[p] for p=1:P) >= Flow_bw[flow_index])

    #path routing constraint
    @constraint(myModel, sum(temp_r[p] for p=1:P) <= 1)

    #The primal
    term1=0
    for p=1:P
        term1=term1 + r[p]*path_cost[p]*r[p] + r[p]*packet_lost[p]
        #term1=term1 + r[p]*(path_cost[p] + packet_lost[p])
    end
    print("\nterm 1 ",term1)
    #The dual
    B=[[0,  10, 20, 0,  0],
   [10, 0,  20, 10 ,0],
   [20, 20, 0,  0, 20],
   [0,  10, 0,  0, 20],
   [0,  0,  20, 20, 0]]

    term2=0
    for i=1:length(cons)
        p1=Int(cons[i][1])
        p2=Int(cons[i][2])
        term2=term2+lambdae[i]*(-cons[i][3]+new_rate[2,p2]+r[p1])
    end
    print("\n term2 ",term2)
    #print("\n\nChecking term 2:\n")
    print("\n\n***********************************************************")
    #cost function
    obj=term1/scale+term2
    println("\nObjective function: ", obj)

    @objective(myModel, Min, obj)
    status = solve(myModel)
    r_value = getValue(r)
    #println("rate",r_value)
    println("\nThis is my model:\n", myModel)
    println("test status:",status)
    #println("test temp",getValue(temp_r))

    global objective_value = getObjectiveValue(myModel)
    print(objective_value)
    if(NaN in Set(r_value))==false
        print("solve")
        return r_value
    end
end
sigma=0.5
gamma=5.9
epsilon=0.1
#k=0
boolean="False"
new_rate=zeros(F,P)
value_rate=[]
for k=1:iteration
#while boolean=="False"
    #k=k+1
    # http://www.juliaopt.org/JuMP.jl/v0.13/refmodel.html
    print("\n\n==========================================================================================\n")
    print("\n**************************** This is the iteration number: ", k, " ****************************\n")
    print("\n==========================================================================================\n")

    #print("\n\n======== Checking: rate & new_rate values: \n", "rate: ", rate, " & new_rate: ", new_rate)

    #optimize for each flow
    temp=zeros(F,P)
    for i=1:F
        #print("\n\n\nThe optimization for flow: ", i, " starts...")

        #new_rate[i,:]=Optimize_flow_f(rate[i,:],i)
        rate[i,:]=Optimize_flow_f(new_rate[i,:],i)
        flow_objective_value[i,k] = objective_value
        println("flow ",i,"rate:",rate)
        for j=1:P
            if rate[i,j]<0.0001
                new_rate[i,j] =0
            else
                new_rate[i,j]=rate[i,j]
            end
            #print("\nThe calculated rate, of flow ", i, " in path ", j, " is equal to: ", rate[i,j])
            #temp[i,j]=rate[i,j]/(new_rate[i,j]+delta)
            #print("\nThe temp[i,j] is equal to:", temp[i,j])
        end
        #println("\nCheck1: new_rate[i,:] : ", new_rate[i,:])
        #println("\nCheck1: rate[i,:] : ", rate[i,:])
        #new_rate[i,:] =rate[i,:]

        #println("\nCheck2: new_rate[i,:] : ", new_rate[i,:])
        #println("\nCheck2: rate[i,:] : ", rate[i,:])
    end

    #dual update
    print("\n\n****************************************")
    print("\nDual update caculation for each link...")
    print("newrate:",new_rate)
    append!(value_rate,new_rate)
    for s=1:length(cons)
        index_pathflow1=Int(cons[s][1])
        index_pathflow2=Int(cons[s][2])

        #println("violate cap ",cons[s][3],":",(-new_rate[1,indexflow1]-new_rate[2,indexflow2]+cons[s][3]))
        temp=lambdae[s]
        lambdae[s] =lambdae[s]+gamma*(+new_rate[1,index_pathflow1]+new_rate[2,index_pathflow2]-cons[s][3])
        if(-new_rate[1,index_pathflow1]-new_rate[2,index_pathflow2]+cons[s][3])<0
            print("\nbw:",-new_rate[1,index_pathflow1]-new_rate[2,index_pathflow2]+cons[s][3])
            print("\nconstraint:",index_pathflow1,":",index_pathflow2,":",cons[s][3])
            print("\nnegative ",new_rate[1,index_pathflow1],new_rate[2,index_pathflow2])
                        #lambdae[s]=temp+100*abs((-new_rate[1,index_pathflow1]-new_rate[2,index_pathflow2]+cons[s][3]))
        end
        if lambdae[s]<0
            lambdae[s]=0
        end
        #println(lambdae[s])
    end
end
