"""
Subset quality metrics
"""
function coverage(fileCount::Int, subset::Array{Int}, totalN::Int)
  reachableNeighbours = @parallel (vcat) for i in subset
    tempNeighbours = Array(Array{Int},0)
    for j in 1:fileCount
      index = findfirst(fileData[j][1], i)
      if index != 0
        neighbours = fileData[j][3][:, index]
        found = findfirst(neighbours, -1)
        if found != 0
          deleteat!(neighbours, found:length(neighbours))
        end
        push!(tempNeighbours, neighbours)
      end
    end
    tempNeighbours = vcat(tempNeighbours...)
  end
  C = length(unique(reachableNeighbours))/totalN
  return C
end

"""
Clustering quality metrics
"""
function modularity(fileCount::Int, k::Int, weighted::Int)
  println("Modularity metric:")
  clusterResult = fill(0.0, k, 2)
  @inbounds for i in 1:k
    println("Cluster $(i)")
    tic()
    # nodes belonging to cluster i
    clusterNodes = qData[(qData[:,2] .== i), 1]
    println("Cluster size $(length(clusterNodes))")
    clusterValues = @parallel (.+) for node in clusterNodes
      clusterAssoc = 0.0
      clusterDeg = 0.0
      foundNeighbours = Int[]
      foundClusterNeighbours = Int[]
      @inbounds for j in 1:fileCount
        indexNode = findfirst(fileData[j][1], node)
        if indexNode != 0
          # degree of new node neighbours
          nodeNeighbours = find(x -> x ∉ foundNeighbours, fileData[j][3][:, indexNode])
          if length(nodeNeighbours) != 0
            push!(foundNeighbours, fileData[j][3][nodeNeighbours, indexNode]...)
            # if unweighted
            if weighted == 0
              clusterDeg += float(length(nodeNeighbours))
            else
              clusterDeg += sum(fileData[j][4][nodeNeighbours, indexNode])
            end
          end
          # degree of new node cluster neighbours
          nodeNeighbours = find(x -> x ∉ foundClusterNeighbours, fileData[j][3][:, indexNode])
          clusterNeighbours = find(x -> x in clusterNodes, fileData[j][3][nodeNeighbours, indexNode])
          if length(clusterNeighbours) != 0
            push!(foundClusterNeighbours, fileData[j][3][nodeNeighbours, indexNode][clusterNeighbours]...)
            # if unweighted
            if weighted == 0
              clusterAssoc += float(length(clusterNeighbours))
            else
              clusterAssoc += sum(fileData[j][4][nodeNeighbours, indexNode][clusterNeighbours])
            end
          end
        end
      end
      # end parallel loop
      [clusterAssoc, clusterDeg]
    end
    clusterResult[i, :] = clusterValues
    toc()
  end
  # total sum of weights of all edges in network
  edgesSum = sum(clusterResult[:, 2])
  Q = 0.0
  @inbounds for i in 1:k
    Q += clusterResult[i, 1]/edgesSum - (clusterResult[i, 2]/edgesSum)^2
  end
  return Q
end


function cut_conductance(fileCount::Int, qtest::Array{Int, 2})
  volumeG, volumeS = @parallel (+) for j in 1:fileCount
    #sumS = fileData[j][1]
    #sumG = sum(fileData[j][2])
    d = sum(fileData[j][2], 2)
    d, 0.0
  end
  return volumeG
end

#=
function modularity1(fileCount::Int, k::Int, weighted::Int)
  println("Modularity metric:")
  # total sum of weights of all edges in network
  edgesSum = @parallel (+) for j in 1:fileCount
    sum(fileData[j][2])
  end
  Xsum = fill(0.0, k, 2)
  @inbounds for j in 1:fileCount
    println("File $(j)")
    tic()
    @inbounds for i in 1:k
      println("Cluster $(i)")
      # node indices of file j belonging to cluster i
      clusterNodes = find(x -> x in qData[(qData[:,2] .== i), 1], fileData[j][1])
      X = @parallel (.+) for index in clusterNodes
        # find indices of neighbours of node also in cluster i
        clusterNeighbours = find(x -> x in clusterNodes, fileData[j][3][:, index])
        # cluster association and cluster degree value
        if length(clusterNeighbours) != 0
          if weighted == 0
            [float(length(clusterNeighbours)) fileData[j][2][index]]
          else
            [sum(fileData[j][4][clusterNeighbours, index]) sum(fileData[j][4][:, index])]
          end
        else
          [0.0 0.0]
        end
      end
      Xsum[i, :] .+= [X[1] X[2]]
    end
    toc()
  end
  Q = 0.0
  @inbounds for i in 1:k
    Q = Q + Xsum[i, 1]/edgesSum - (Xsum[i, 2]/edgesSum)^2
  end
  return Q
end

function modularity3(fileCount::Int, k::Int, weighted::Int)
  println("Modularity metric:")
  # total sum of weights of all edges in network
  edgesSum = @parallel (+) for j in 1:fileCount
    sum(fileData[j][2])
  end
  Xsum = @parallel (.+) for j in 1:fileCount
    tic()
    X = fill(0.0, k, 2)
    @inbounds for i in 1:k
      println("Cluster $(i) on File $(j)")
      clusterAssoc = 0.0
      clusterDeg = 0.0
      # node indices of file j belonging to cluster i
      clusterNodes = find(x -> x in qData[(qData[:,2] .== i), 1], fileData[j][1])
      println("$(length(clusterNodes)) Nodes on Processor $(myid())")
      @inbounds for index in clusterNodes
        # find indices of neighbours of node also in cluster i
        clusterNeighbours = find(x -> x in clusterNodes, fileData[j][3][:, index])
        # cluster association and cluster degree value
        if length(clusterNeighbours) != 0
          if weighted == 0
            clusterAssoc += float(length(clusterNeighbours))
            clusterDeg += fileData[j][2][index]
          else
            clusterAssoc += sum(fileData[j][4][clusterNeighbours, index])
            clusterDeg += sum(fileData[j][4][:, index])
          end
        end
      end
      X[i, :] = [clusterAssoc clusterDeg]
    end
    toc()
    X
  end
  Q = 0.0
  @inbounds for i in 1:k
    Q = Q + Xsum[i, 1]/edgesSum - (Xsum[i, 2]/edgesSum)^2
  end
  return Q
end
=#
