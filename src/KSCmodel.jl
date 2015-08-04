"""
Train model through eigendecomposition
"""
function modelTrain(trainSparse::SparseMatrixCSC{Float64, Int}, trainN::Int, maxK::Int, decomp::String)
  # kernel matrix train nodes
  trainDiag = spdiagm(vec(sqrt(sum(trainSparse.^2, 2))), 0, trainN, trainN)
  trainNorm = \(trainDiag, trainSparse)
  trainDiag = []
  kernelTrain = trainNorm * trainNorm'
  kernelTrain = 0.5*(kernelTrain + kernelTrain')
  inverseD = 1./(sum(kernelTrain, 1)')
  inverseDsparse = spdiagm(vec(inverseD), 0)
  # centering matrix
  centerMatrix = eye(trainN) - ones(trainN, 1) * inverseD' / sum(inverseD)
  # eigendecomposition
  if decomp == "svd"
    # svd support for sparse matrices in version 0.4
    # use svds() with economy size
    error("modelTrain: svd not yet available.")
  elseif decomp == "eigs"
    # Lanczos method
    λ,α = eigs(inverseDsparse * centerMatrix * kernelTrain, nev=maxK+4, which=:LM)
  elseif decomp == "eig"
    # full eigendecomposition
    λ,α = eig(inverseDsparse * centerMatrix * kernelTrain)
  else
    error("modelTrain: wrong decomp argument.")
  end
  λ = real(λ)
  λ[isnan(λ)] = 0.0
  α = real(α)
  if length(λ) < maxK
    println("Less eigenvalues than requested: k = $(maxK), λ = $(length(λ))")
  end
  # sort eigenvalues in descending order and only keep top k-1
  indexλ = resize!(sortperm(λ, rev=true), maxK-1)
  λ = λ[indexλ]
  α = α[:, indexλ]
  # bias terms
  β = (1 / (sum(inverseD))) * (inverseD' * kernelTrain * α);
  # score variables
  eTrain = kernelTrain * α + repeat(β, outer = [trainN, 1])
  CB = Array(Array{Int,2}, maxK-1)
  # generate codebooks for each possible number of clusters k
  for k in 2:maxK
    CB[k-1] = codebook(eTrain[:, 1:k-1], k)
  end
  return trainNorm, λ, α, β, CB
end

"""
Validate model through Balanced Angular Fit criterion (BAF)
"""
function modelValid(validSparse::SparseMatrixCSC{Float64, Int}, validN::Int, trainNorm::SparseMatrixCSC{Float64, Int}, maxK::Int, minK::Int, α::Array{Float64, 2}, β::Array{Float64, 2}, CB::Array{Array{Int, 2}})
  # kernel matrix valid nodes
  kernelValid = validSparse*trainNorm';
  eValid = kernelValid * α + repeat(β, outer = [validN, 1])
  # out-of-sample extension
  qValid = Array(Array{Int,2}, maxK-1)
  for numk in 2:maxK
    qValid[numk-1] = membership(eValid[:,1:numk-1], CB[numk-1])
  end
  # estimate optimal number of clusters using Balanced Angular Fit (BAF)
  baf = Array(Float64, maxK-minK+1)
  for k in 1:maxK-minK+1
    baf[k] = mean(balanced_angular_fit(eValid[:, 1:minK-2+k], qValid[minK-2+k]))
  end
  valueBAF,indexBAF = findmax(baf)
  kBAF = indexBAF+minK-1
  return baf, kBAF, indexBAF, valueBAF
end

"""
Assign cluster indices to unweighted network
Blocks of nodes in parallel
"""
function modelTestU(fileCount::Int, trainNorm::SparseMatrixCSC{Float64, Int}, minN::Int, maxN::Int, α::Array{Float64,2}, β::Array{Float64,2}, CB::Array{Array{Int,2}}, kBAF::Int, maxMB::Int; limitN=-1)
  lastNeighbour = maxN
  if minN == 0
    lastNeighbour += 1
  end
  if limitN >= 0 && limitN <= maxN
    endNode = limitN
  else
    endNode = maxN
  end
  avgLength = @parallel (+) for i in 1:fileCount
    mean(fileData[i][2])
  end
  # determine number of nodes tested simultaneously
  stepSize = min(int(floor(maxMB*1000000*0.5 / (16*nprocs()*avgLength))), int(ceil((endNode-minN+1)/nprocs())))
  println("Step size $(stepSize)")
  countN = @parallel (+) for n in colon(minN,stepSize,endNode+stepSize)
    if n > endNode
      return 0
    end
    if n+stepSize > endNode
      nodes = [n:endNode]
    else
      nodes = [n:n+stepSize-1]
    end
    testNodes = Array(Int,0)
    testNeighbours = Array(Int,0)
    @inbounds for fileIndex in 1:fileCount
      foundNodes = 0
      @inbounds for index in 1:length(fileData[fileIndex][1])
        if fileData[fileIndex][1][index] in nodes
          neighbours = fileData[fileIndex][3][:, index]
          foundEnd = findfirst(neighbours, -1)
          if foundEnd != 0
            deleteat!(neighbours, foundEnd:length(neighbours))
          end
          foundNodes += 1
          append!(testNodes, fill(fileData[fileIndex][1][index], length(neighbours)))
          append!(testNeighbours, neighbours)
        end
        if foundNodes == length(nodes)
          break
        end
      end
    end
    permIndices = unique(testNodes)
    testData = fill(-1, length(nodes), 2)
    if length(permIndices) != 0
      # avoid zero index in neighbours
      if minN == 0
        testNeighbours[testNeighbours .== 0] = lastNeighbour
      end
      @inbounds for i in 1:length(testNodes)
        testNodes[i] = findfirst(permIndices, testNodes[i])
      end
      # kernel matrix single test node
      testSparse = sparse(testNodes, testNeighbours, fill(1.0, length(testNodes)), length(permIndices), lastNeighbour, (w1,w2)->w1)
      kernelTest = testSparse*trainNorm';
      eTest = kernelTest * α + repeat(β, outer = [length(permIndices), 1])
      qTest = membership(eTest[:,1:kBAF-1], CB[kBAF-1])
      # store node and its cluster membership
      @inbounds for i in 1:length(permIndices)
        foundPerm = findfirst(nodes, permIndices[i])
        if foundPerm != 0
          testData[foundPerm, :] = [permIndices[i], qTest[i]]
        end
      end
    end
    qData[nodes .- (minN-1), :] = testData
    length(permIndices)
  end
  println("Nodes tested $(countN)")
end

"""
Assign cluster indices to weighted network
Blocks of nodes in parallel
"""
function modelTestW(fileCount::Int, trainNorm::SparseMatrixCSC{Float64, Int}, minN::Int, maxN::Int, α::Array{Float64,2}, β::Array{Float64,2}, CB::Array{Array{Int,2}}, kBAF::Int, maxMB::Int; limitN=-1)
  lastNeighbour = maxN
  if minN == 0
    lastNeighbour += 1
  end
  if limitN >= 0 && limitN <= maxN
    endNode = limitN
  else
    endNode = maxN
  end
  avgLength = @parallel (+) for i in 1:fileCount
    mean(fileData[i][2])
  end
  # determine number of nodes tested simultaneously
  stepSize = min(int(floor(maxMB*1000000*0.5 / (24*nprocs()*avgLength))), int(ceil((endNode-minN+1)/nprocs())))
  println("Step size $(stepSize)")
  countN = @parallel (+) for n in colon(minN,stepSize,endNode+stepSize)
    if n > endNode
      return 0
    end
    if n+stepSize > endNode
      nodes = [n:endNode]
    else
      nodes = [n:n+stepSize-1]
    end
    testNodes = Array(Int,0)
    testNeighbours = Array(Int,0)
    testWeights = Array(Float64,0)
    @inbounds for fileIndex in 1:fileCount
      foundNodes = 0
      @inbounds for index in 1:length(fileData[fileIndex][1])
        if fileData[fileIndex][1][index] in nodes
          neighbours = fileData[fileIndex][3][:, index]
          weights = fileData[fileIndex][4][:, index]
          foundEnd = findfirst(neighbours, -1)
          if foundEnd != 0
            deleteat!(neighbours, foundEnd:length(neighbours))
            deleteat!(weights, foundEnd:length(weights))
          end
          foundNodes += 1
          append!(testNodes, fill(fileData[fileIndex][1][index], length(neighbours)))
          append!(testNeighbours, neighbours)
          append!(testWeights, weights)
        end
        if foundNodes == length(nodes)
          break
        end
      end
    end
    permIndices = unique(testNodes)
    testData = fill(-1, length(nodes), 2)
    if length(permIndices) != 0
      # avoid zero index in neighbours
      if minN == 0
        testNeighbours[testNeighbours .== 0] = lastNeighbour
      end
      @inbounds for i in 1:length(testNodes)
        testNodes[i] = findfirst(permIndices, testNodes[i])
      end
      # kernel matrix single test node
      testSparse = sparse(testNodes, testNeighbours, testWeights, length(permIndices), lastNeighbour, (w1,w2)->w1)
      kernelTest = testSparse*trainNorm';
      eTest = kernelTest * α + repeat(β, outer = [length(permIndices), 1])
      qTest = membership(eTest[:,1:kBAF-1], CB[kBAF-1])
      # store node and its cluster membership
      @inbounds for i in 1:length(permIndices)
        foundPerm = findfirst(nodes, permIndices[i])
        if foundPerm != 0
          testData[foundPerm, :] = [permIndices[i], qTest[i]]
        end
      end
    end
    qData[nodes .- (minN-1), :] = testData
    length(permIndices)
  end
  println("Nodes tested $(countN)")
end

"""
Assign cluster indices to network
Single node in parallel -> inefficient
"""
function modelTestSingle(fileCount::Int, trainNorm::SparseMatrixCSC{Float64, Int}, minN::Int, maxN::Int, α::Array{Float64,2}, β::Array{Float64,2}, CB::Array{Array{Int,2}}, k::Int)
  lastNeighbour = maxN
  if minN == 0
    lastNeighbour += 1
  end
  countN = @parallel (+) for n in minN:maxN
    testNeighbours = Array(Array{Int},0)
    testWeights = Array(Array{Float64},0)
    for fileIndex in 1:fileCount
      index = findfirst(fileData[fileIndex][1], n)
      if index != 0
        neighbours = fileData[fileIndex][3][:, index]
        weights = fileData[fileIndex][4][:, index]
        foundIndex = findfirst(neighbours, -1)
        if foundIndex != 0
          deleteat!(neighbours, foundIndex:length(neighbours))
          deleteat!(weights, foundIndex:length(weights))
        end
        push!(testNeighbours, neighbours)
        push!(testWeights, weights)
      end
    end
    testNeighbours = vcat(testNeighbours...)
    testWeights = vcat(testWeights...)
    if length(testNeighbours) != 0
      # avoid zero index in neighbours
      if minN == 0
        testNeighbours[testNeighbours .== 0] = lastNeighbour
      end
      # kernel matrix single test node
      testSparse = sparse(fill(1, length(testNeighbours)), testNeighbours, testWeights, 1, lastNeighbour, (w1,w2)->w1)
      kernelTest = testSparse*trainNorm';
      eTest = kernelTest * α + β
      qTest = membership(eTest[:,1:k-1], CB[k-1])
      # store node and its cluster membership
      qData[n-minN+1, :] = [n qTest]
    else
      qData[n-minN+1, :] = [-1 -1]
    end
    1
  end
  println("Nodes tested $(countN)")
end

"""
Helper functions
"""
function codebook(e::Array{Float64,2}, numk::Int)
  eBin = int(sign(e))
  eUniBin = unique(eBin, 1)
  # TODO method unique does not return index vectors (JuliaLang issue #1845)
  # current solution = nested for-loop comparison
  eCW = Array(Int, size(eBin,1))
  @inbounds for i in 1:size(eBin, 1)
    @inbounds for j in 1:size(eUniBin,1)
      if eUniBin[j,:] == eBin[i,:]
        eCW[i] = j
        break;
      end
    end
  end
  # calculate occurence of each codeword
  m = size(eUniBin, 1)
  eSizeCW = Array(Int, m)
  @inbounds for i in 1:m
    eSizeCW[i] = sum(eCW .== i)
  end
  # sort codewords by descending occurence
  indexCW = sortperm(eSizeCW, rev=true)
  if m < numk
    println("Redundant codewords for $(numk) clusters")
    return eUniBin[indexCW[1:m], :]
  end
  return eUniBin[indexCW[1:numk], :]
end

function membership(e::Array{Float64,2}, CB::Array{Int,2})
  eBin = int(sign(e))
  # Hamming distance is equal to squared Euclidean when using binary vectors
  minValues, qVectors = findmin(dist_euclidean(eBin, CB), 2)
  return int(ceil(qVectors ./ size(qVectors,1)))
end

function balanced_angular_fit(e::Array{Float64,2}, q::Array{Int,2})
  eN, eK = size(e)
  # calculate mean of each cluster
  clustermean = zeros(eK+1, eK)
  for i in 1:eK+1
    # TODO mean function along dim, JuliaLang issue #2265
    clustermean[i, :] = [mean(e[q[:] .== i, :][:,j]) for j in 1:eK]
  end
  clustermeanNorm = sqrt(sum(clustermean.^2, 2))
  # angular similarity between eigenspace projections and each cluster mean
  clusterValues = Array(Float64, eN)
  clusterIndices = Array(Int, eN)
  for i in 1:eN
    clusterSim = clustermean * e[i, :]' ./ (clustermeanNorm * norm(e[i, :]));
    maxValue, maxIndex = findmax(clusterSim, 1)
    clusterValues[i] = maxValue[1]
    clusterIndices[i] = maxIndex[1]
  end
  # cosine similarity entire cluster
  clusterFit = Array(Float64, eK+1)
  for i in 1:eK+1
    indices = findin(clusterIndices, i)
    if isempty(indices)
      clusterFit[i] = 0.0
    else
      clusterFit[i] = mean(clusterValues[indices])
    end
  end
  return clusterFit
end

function dist_euclidean(x::Array{Int,2}, y::Array{Int,2})
  xDim = size(x, 1)
  yDim = size(y, 1)
  xSum = sum(x.*x, 2)
  ySum = sum(y'.*y', 1)
  return sqrt(repeat(xSum, outer=[1,yDim]) + repeat(ySum, outer=[xDim,1]) - 2*x*y')
end
