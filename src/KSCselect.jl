"""
Training and Validation subset selection through either:
- Fast and Unique Representative Subset selection
- Random Representative Subset selection
"""
function selectSubset(method::String, fileIndex::Int, trainSize::Int, validSize::Int, degreeThreshold::Float64)
  if method == "furs"
    sortedNodeIndices = sortIndicesDegree(fileData[fileIndex][2], degreeThreshold)
    trainSubset, sortedNodeIndices, trainComplete = furs(fileIndex, trainSize, sortedNodeIndices)
    # TODO: remove training subset from neighbours/degree of remaining nodes
    validSubset, sortedNodeIndices, validComplete = furs(fileIndex, validSize, sortedNodeIndices)
    return (trainSubset, trainComplete, validSubset, validComplete)
  elseif method == "randrs"
    # get all node indices above the degree threshold
    nodeIndices = find(x -> x > degreeThreshold, fileData[fileIndex][2])
    trainSubset, trainComplete = randrs(trainSize, nodeIndices)
    validSubset, validComplete = randrs(validSize, nodeIndices)
    return (trainSubset, trainComplete, validSubset, validComplete)
  else
    error("selectSubset: wrong method argument.")
  end
end

"""
Fast and Unique Representative Subset selection (FURS)
Return array of selected node indices and boolean value (true if equal to requested sampleSize)
"""
function furs(fileIndex::Int, sampleSize::Int, sortedNodeIndices::Vector{Int})
  selectedSize = 0
  selectedSubset = Array(Int, 0)
  tempIndices = Array(Int, 0)
  while selectedSize < sampleSize
    # check if still viable nodes left
    if isempty(sortedNodeIndices)
      if isempty(tempIndices)
        # FURS unable to finnish
        return (selectedSubset, sortedNodeIndices, false)
      else
        # add and sort the removed indices back into sortedNodeIndices
        sortedNodeIndices = tempIndices[sortIndicesDegree(fileData[fileIndex][2][tempIndices])]
        empty!(tempIndices)
      end
    else
      # get index of node with highest degree
      nodeIndex = splice!(sortedNodeIndices, 1)
      neighbourArray = fileData[fileIndex][3][:, nodeIndex]
      push!(selectedSubset, nodeIndex)
      selectedSize += 1
      # remove selected node's neighbours from sortedNodeIndices and store them in tempIndices
      for neighbour in neighbourArray
        # stop loop if no more candidate nodes
        if isempty(sortedNodeIndices)
          break
        end
        # stop loop when -1 padding is reached
        if neighbour == -1
          break
        end
        index = findfirst(fileData[fileIndex][1][sortedNodeIndices], neighbour)
        if index != 0
          push!(tempIndices, splice!(sortedNodeIndices, index))
        end
      end
    end
  end
  # add tempIndices back into sortedNodeIndices
  if !isempty(tempIndices)
    push!(sortedNodeIndices, tempIndices...)
    sortedNodeIndices = sortedNodeIndices[sortIndicesDegree(fileData[fileIndex][2][sortedNodeIndices])]
  end
  return (selectedSubset, sortedNodeIndices, true)
end

"""
Random representative subset selection
Return array of randomly selected node indices and boolean value (true if equal to requested sampleSize)
"""
function randrs(sampleSize::Int, nodeIndices::Vector{Int})
  selectedSize = 0
  selectedSubset = Array(Int, 0)
  while selectedSize < sampleSize
    if isempty(nodeIndices)
      return (selectedSubset, false)
    end
    # randomly sample nodeIndices
    sample = rand(1:length(nodeIndices))
    push!(selectedSubset, splice!(nodeIndices, sample))
    selectedSize += 1
  end
  return (selectedSubset, true)
end

"""
Node index sorting functions
"""
function sortIndicesDegree(degreeArray::Array{Float64, 2}, threshold::Float64)
  # remove all node indices with degree lower than threshold
  nodeIndices = find(x -> x >= threshold, degreeArray)
  # get sorted node indices based on descending degree
  return nodeIndices[sortperm(degreeArray[nodeIndices], rev=true)]
end

function sortIndicesDegree(degreeArray::Array{Float64, 1})
  # get sorted node indices based on descending degree
  return sortperm(degreeArray, rev=true)
end
