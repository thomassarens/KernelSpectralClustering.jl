"""
Subset quality metrics
"""
function coverage(fileCount::Int, subset::Array{Int}, totalN::Int)
  println("Coverage metric:")
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

function coverageApprox(fileCount::Int, subset::Array{Int}, minN::Int, maxN::Int)
  println("Coverage metric:")
  totalNodes = maxN-minN+1
  reachableNeighbours = 0
  # BitArray to keep track of unique set of neighbours
  foundNeighbours = falses(totalNodes)
  @inbounds for j in 1:fileCount
    println("File $(j)")
    subsetIndices = find(x -> x in subset, fileData[j][1])'
    @inbounds for nodeIndex in subsetIndices
      nodeNeighbours = fileData[j][3][:, nodeIndex]
      @inbounds for node in nodeNeighbours
        if node == -1
          break
        end
        if !foundNeighbours[node-minN+1]
          reachableNeighbours += 1
          foundNeighbours[node-minN+1] = true
        end
      end
    end
  end
  C = reachableNeighbours / totalNodes
  return C, totalNodes
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

function modularityApprox(fileCount::Int, k::Int, weighted::Int)
  println("Modularity metric:")
  clusterResult = fill(0.0, k, 2)
  @inbounds for i in 1:k
    println("Cluster $(i)")
    # nodes belonging to cluster i
    clusterNodes = qData[(qData[:,2] .== i), 1]
    println("Cluster size $(length(clusterNodes))")
    clusterDeg, clusterAssoc = @parallel (.+) for j in 1:fileCount
      tic()
      fileDeg = 0.0
      fileAssoc = 0.0
      fileNodes = find(x -> x in clusterNodes, fileData[j][1])
      fileNeighbours = fileData[j][3][:, fileNodes]
      # sum of degree of cluster nodes' neighbours
      if weighted == 0
        fileDeg = float(length(fileNeighbours[fileNeighbours .>= 0]))
      else
        fileDeg = sum(fileData[j][4][:, fileNodes])
      end
      # sum of degree of cluster nodes' in cluster neighbours
      clusterNeighbours = find(x -> x in clusterNodes, fileNeighbours)
      if length(clusterNeighbours) != 0
        if weighted == 0
          fileAssoc = float(length(clusterNeighbours))
        else
          fileAssoc = sum(fileData[j][4][:, fileNodes][clusterNeighbours])
        end
      end
      println("File $(j) done")
      toc()
      [fileDeg, fileAssoc]
    end
    clusterResult[i, :] = [clusterDeg, clusterAssoc]
  end
  # total sum of weights of all edges in network
  edgesSum = sum(clusterResult[:, 2])
  Q = 0.0
  @inbounds for i in 1:k
    Q += clusterResult[i, 1]/edgesSum - (clusterResult[i, 2]/edgesSum)^2
  end
  return Q
end

#=
function modularityApprox(fileCount::Int, k::Int, weighted::Int)
  println("Modularity metric:")
  clusterResult = fill(0.0, k, 2)
  @inbounds for i in 1:k
    println("Cluster $(i)")
    # nodes belonging to cluster i
    clusterNodes = qData[(qData[:,2] .== i), 1]
    println("Cluster size $(length(clusterNodes))")
    @inbounds for j in 1:fileCount
      println("File $(j)")
      tic()
      fileNodes = find(x -> x in clusterNodes, fileData[j][1])
      fileNeighbours = fileData[j][3][:, fileNodes]
      # sum of degree of cluster nodes' neighbours
      if weighted == 0
        clusterDeg = float(length(fileNeighbours[fileNeighbours .>= 0]))
      else
        clusterDeg = sum(fileData[j][4][:, fileNodes])
      end
      # sum of degree of cluster nodes' in cluster neighbours
      clusterNeighbours = find(x -> x in clusterNodes, fileNeighbours)
      if length(clusterNeighbours) != 0
        if weighted == 0
          clusterAssoc = float(length(clusterNeighbours))
        else
          clusterAssoc = sum(fileData[j][4][:, fileNodes][clusterNeighbours])
        end
      end
      clusterResult[i, :] .+= [clusterDeg clusterAssoc]
      toc()
    end
  end
  # total sum of weights of all edges in network
  edgesSum = sum(clusterResult[:, 2])
  Q = 0.0
  @inbounds for i in 1:k
    Q += clusterResult[i, 1]/edgesSum - (clusterResult[i, 2]/edgesSum)^2
  end
  return Q
end
=#

function cut_conductance(fileCount::Int)
end
