"""
Memory-Mapped implementation of Kernel Spectral Clustering
"""
function kscbignet(network_file::String, file_dir::String, file_name::String, delimiter::Char, maxMB::Int; fraction=0.15::Float64, method="furs"::String, minK=2::Int, maxK=100::Int, eigfull=false::Bool, convertF=true::Bool)
  if minK < 2 || minK > maxK
    error("minK must be at least 2 and smaller than maxK")
  end
  if !isdir("$(file_dir)/$(file_name)")
    mkdir("$(file_dir)/$(file_name)")
  end
  resultfile = open("$(file_dir)/$(file_name)/resultKSC.txt", "w")
  if convertF
    println("File Conversion")
    tic()
    fileCount, firstNode, lastNode, fileWeighted = convertDSVNet(network_file, file_dir, file_name, delimiter, maxMB)
    timeConversion = toq()
    println("File Conversion step finished, time elapsed $(timeConversion)s")
    write(resultfile, "File Conversion step: time $(timeConversion)s\n")
  else
    file = open("$(file_dir)/$(file_name)/info.bin")
    fileCount = read(file, Int)
    firstNode = read(file, Int)
    lastNode = read(file, Int)
    fileWeighted = read(file, Int)
    close(file)
  end
  println("Assigning global variables")
  # stores mmaps of nodeArray, degreeArray, neighbourMatrix, weightMatrix for each file on each processor
  for p in procs()
    if fileWeighted == 0
      @spawnat(p, netmapU(file_dir, file_name, fileCount))
    else
      @spawnat(p, netmapW(file_dir, file_name, fileCount))
    end
    @spawnat(p, qmap("$(file_dir)/$(file_name)/qtest.bin", firstNode, lastNode))
  end
  # TODO make max subset size dynamic based on maxMB?
  sampleSize = min(5000, int(ceil((lastNode-firstNode)*fraction)))
  sampleSizeBlock = int(floor(sampleSize/fileCount))
  println("Total sample size for $(fileCount) files = $(sampleSize)")
  println("Subset selection")
  tic()
  # perform the subset selection on each file in parallel (map) and concatenate them together (reduce)
  if fileWeighted == 0
    trainNodes, trainNeighbours, validNodes, validNeighbours = @parallel (subsetCombineU) for fileIndex = 1:fileCount
      # subset selection degree threshold
      degreeThreshold = mean(fileData[fileIndex][2])
      trainIndices, trainComplete, validIndices, validComplete = selectSubset(method, fileIndex, sampleSizeBlock, sampleSizeBlock, degreeThreshold)
      if !trainComplete || !validComplete
        println("Incomplete subset selection on file set $(fileIndex), retrying with median degree threshold")
        degreeThreshold = median(fileData[fileIndex][2])
        trainIndices, trainComplete, validIndices, validComplete = selectSubset(method, fileIndex, sampleSizeBlock, sampleSizeBlock, degreeThreshold)
        if !trainComplete || !validComplete
          println("Incomplete subset selection on file set $(fileIndex), continuing")
        end
      end
      trainNodeSeq = Array(Int, 0)
      trainNeighbourSeq = Array(Int, 0)
      validNodeSeq = Array(Int, 0)
      validNeighbourSeq = Array(Int, 0)
      generateSubsetU(fileIndex, trainIndices, trainNodeSeq, trainNeighbourSeq)
      generateSubsetU(fileIndex, validIndices, validNodeSeq, validNeighbourSeq)
      # return training and validation data as tuple
      (trainNodeSeq, trainNeighbourSeq, validNodeSeq, validNeighbourSeq)
    end
  else
    trainNodes, trainNeighbours, trainWeights, validNodes, validNeighbours, validWeights = @parallel (subsetCombineW) for fileIndex = 1:fileCount
      # subset selection degree threshold
      degreeThreshold = mean(fileData[fileIndex][2])
      trainIndices, trainComplete, validIndices, validComplete = selectSubset(method, fileIndex, sampleSizeBlock, sampleSizeBlock, degreeThreshold)
      if !trainComplete || !validComplete
        println("Incomplete subset selection on file set $(fileIndex), retrying with median degree threshold")
        degreeThreshold = median(fileData[fileIndex][2])
        trainIndices, trainComplete, validIndices, validComplete = selectSubset(method, fileIndex, sampleSizeBlock, sampleSizeBlock, degreeThreshold)
        if !trainComplete || !validComplete
          println("Incomplete subset selection on file set $(fileIndex), continuing")
        end
      end
      trainNodeSeq = Array(Int, 0)
      trainNeighbourSeq = Array(Int, 0)
      trainWeightSeq = Array(Float64, 0)
      validNodeSeq = Array(Int, 0)
      validNeighbourSeq = Array(Int, 0)
      validWeightSeq = Array(Float64, 0)
      generateSubsetW(fileIndex, trainIndices, trainNodeSeq, trainNeighbourSeq, trainWeightSeq)
      generateSubsetW(fileIndex, validIndices, validNodeSeq, validNeighbourSeq, validWeightSeq)
      # return training and validation data as tuple
      (trainNodeSeq, trainNeighbourSeq, trainWeightSeq, validNodeSeq, validNeighbourSeq, validWeightSeq)
    end
  end
  timeSubset = toq()
  println("Subset selection step finished, time elapsed $(timeSubset)s")
  write(resultfile, "Subset selection step: time $(timeSubset)s\n")
  # TODO: store to file?
  trainSubset = unique(trainNodes)
  validSubset = unique(validNodes)
  # change node indices to be sequential to avoid zeros in matrix representation
  permIndices = trainSubset
  trainN = length(permIndices)
  #[trainNodes[trainNodes .== permIndices[value]] = value for value in 1:length(permIndices)]
  @inbounds for i in 1:length(trainNodes)
    trainNodes[i] = findfirst(permIndices, trainNodes[i])
  end
  permIndices = validSubset
  validN = length(permIndices)
  #[validNodes[validNodes .== permIndices[value]] = value for value in 1:length(permIndices)]
  @inbounds for i in 1:length(validNodes)
    validNodes[i] = findfirst(permIndices, validNodes[i])
  end
  permIndices = []
  # avoid zero index in neighbours
  lastNeighbour = lastNode
  if firstNode == 0
    lastNeighbour += 1
    trainNeighbours[trainNeighbours .== 0] = lastNeighbour
    validNeighbours[validNeighbours .== 0] = lastNeighbour
  end
  # training phase
  println("Training model")
  tic()
  if fileWeighted == 0
    modelSparse = sparse(trainNodes, trainNeighbours, fill(1.0, length(trainNodes)), trainN, lastNeighbour)
  else
    modelSparse = sparse(trainNodes, trainNeighbours, trainWeights, trainN, lastNeighbour)
  end
  trainNorm, λ, α, β, CB = modelTrain(modelSparse, trainN, maxK, eigfull)
  timeTrain = toq()
  println("Training phase finished, time elapsed $(timeTrain)s")
  write(resultfile, "Training phase: time $(timeTrain)s\n")
  # validation phase
  if fileWeighted == 0
    modelSparse = sparse(validNodes, validNeighbours, fill(1.0, length(validNodes)), validN, lastNeighbour)
  else
    modelSparse = sparse(validNodes, validNeighbours, validWeights, validN, lastNeighbour)
  end
  println("Validating model")
  tic()
  baf, kBAF, indexBAF, valueBAF = modelValid(modelSparse, validN, trainNorm, maxK, minK, α, β, CB)
  timeValid = toq()
  println("Validation phase finished (k = $(kBAF)), time elapsed $(timeValid)s")
  write(resultfile, "Validation phase: BAF $(valueBAF), k $(kBAF), time $(timeValid)s\n")
  modelSparse = []
  # test phase
  println("Test phase")
  tic()
  if fileWeighted == 0
    modelTestU(fileCount, trainNorm, firstNode, lastNode, α, β, CB, kBAF, maxMB)
  else
    modelTestW(fileCount, trainNorm, firstNode, lastNode, α, β, CB, kBAF, maxMB)
  end
  timeTest = toq()
  println("Test phase finished, time elapsed $(timeTest)s")
  write(resultfile, "Test phase: time $(timeTest)s\n")
  close(resultfile)
  return qData, baf, kBAF, firstNode, lastNode, trainSubset, fileCount, fileWeighted
end

function subsetCombineU(x::Tuple, y::Tuple)
  return vcat(x[1], y[1]), vcat(x[2], y[2]), vcat(x[3], y[3]), vcat(x[4], y[4])
end

function subsetCombineW(x::Tuple, y::Tuple)
  return vcat(x[1], y[1]), vcat(x[2], y[2]), vcat(x[3], y[3]), vcat(x[4], y[4]), vcat(x[5], y[5]), vcat(x[6], y[6])
end

function generateSubsetU(fileIndex::Int, indicesS::Array{Int}, nodesS::Array{Int}, neighboursS::Array{Int})
  for i in indicesS
    tempNeighbourSeq = fileData[fileIndex][3][:, i]
    # remove the -1 and 0 padding
    foundIndex = findfirst(tempNeighbourSeq, -1)
    if foundIndex != 0
      lengthSeq = length(tempNeighbourSeq)
      splice!(tempNeighbourSeq, foundIndex:lengthSeq)
    end
    # remove links between training nodes
    for j in indicesS
      foundIndex = findfirst(tempNeighbourSeq, j)
      if foundIndex != 0
        splice!(tempNeighbourSeq, foundIndex)
      end
    end
    push!(nodesS, fill(fileData[fileIndex][1][i], length(tempNeighbourSeq))...)
    push!(neighboursS, tempNeighbourSeq...)
  end
end

function generateSubsetW(fileIndex::Int, indicesS::Array{Int}, nodesS::Array{Int}, neighboursS::Array{Int}, weightsS::Array{Float64})
  for i in indicesS
    tempNeighbourSeq = fileData[fileIndex][3][:, i]
    tempWeightSeq = fileData[fileIndex][4][:, i]
    # remove the -1 and 0 padding
    foundIndex = findfirst(tempNeighbourSeq, -1)
    if foundIndex != 0
      lengthSeq = length(tempNeighbourSeq)
      splice!(tempNeighbourSeq, foundIndex:lengthSeq)
      splice!(tempWeightSeq, foundIndex:lengthSeq)
    end
    # remove links between training nodes
    for j in indicesS
      foundIndex = findfirst(tempNeighbourSeq, j)
      if foundIndex != 0
        splice!(tempNeighbourSeq, foundIndex)
        splice!(tempWeightSeq, foundIndex)
      end
    end
    push!(nodesS, fill(fileData[fileIndex][1][i], length(tempNeighbourSeq))...)
    push!(neighboursS, tempNeighbourSeq...)
    push!(weightsS, tempWeightSeq...)
  end
end

function netmapU(file_dir::String, file_name, fileCount::Int)
  global fileData = Array((Array{Int, 2}, Array{Float64, 2}, Array{Int, 2}), fileCount)
  @inbounds for i = 1:fileCount
    fileData[i] = mmapNetU(file_dir, file_name, i)
  end
end

function netmapW(file_dir::String, file_name, fileCount::Int)
  global fileData = Array((Array{Int, 2}, Array{Float64, 2}, Array{Int, 2}, Array{Float64, 2}), fileCount)
  @inbounds for i = 1:fileCount
    fileData[i] = mmapNetW(file_dir, file_name, i)
  end
end

function qmap(qFile::String, minN::Int, maxN::Int)
  global qData
  file = open(qFile, "w+")
  qData = mmap_array(Int, (maxN-minN+1, 2), file)
  close(file)
end
