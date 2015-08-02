"""
Memory-Map a set of 3 binary files created by convertDSVNet
Return Vector of nodes, Vector of the corresponding node's degree, Matrix of neighbours of each node
"""
function mmapNetU(file_dir::String, file_name::String, index::Int)
  file = open("$(file_dir)/$(file_name)/nodes$(index).bin")
  columns = read(file, Int)
  nodeArray = mmap_array(Int, (columns, 1), file)
  close(file)
  file = open("$(file_dir)/$(file_name)/degrees$(index).bin")
  degreeArray = mmap_array(Float64, (columns, 1), file)
  close(file)
  file = open("$(file_dir)/$(file_name)/neighbours$(index).bin")
  rows = read(file, Int)
  neighbourMatrix = mmap_array(Int, (rows, columns), file)
  close(file)
  return (nodeArray, degreeArray, neighbourMatrix)
end

"""
Memory-Map a set of 4 binary files created by convertDSVNet
Return Vector of nodes, Vector of the corresponding node's degree, Matrix of neighbours of each node, Matrix of corresponding neighbour's weight
"""
function mmapNetW(file_dir::String, file_name::String, index::Int)
  file = open("$(file_dir)/$(file_name)/nodes$(index).bin")
  columns = read(file, Int)
  nodeArray = mmap_array(Int, (columns, 1), file)
  close(file)
  file = open("$(file_dir)/$(file_name)/degrees$(index).bin")
  degreeArray = mmap_array(Float64, (columns, 1), file)
  close(file)
  file = open("$(file_dir)/$(file_name)/neighbours$(index).bin")
  rows = read(file, Int)
  neighbourMatrix = mmap_array(Int, (rows, columns), file)
  close(file)
  file = open("$(file_dir)/$(file_name)/weights$(index).bin")
  weightMatrix = mmap_array(Float64, (rows, columns), file)
  close(file)
  return (nodeArray, degreeArray, neighbourMatrix, weightMatrix)
end

"""
Converts the network file into sets of binary files, one set each time the maxMB threshold is reached
- *nodes.bin: vector to map internal to external node index
- *degrees.bin: vector of the node degee values
- *neighbours.bin: square array of the neighbours of each node
- *weights.bin: square array of the weight of the corresponding neighbour
Return number of file sets created
"""
function convertDSVNet(dsv_file::String, file_dir::String, file_name::String, delim::Char, maxMB::Int)
  fileCount = 0
  firstNode = -1
  lastNode = -1
  fileWeighted = 1
  # read in file as stream
  fstream = open(dsv_file)
  positionConvert = position(fstream)
  startConvert = false
  lineValues = 0
  while !startConvert || !eof(fstream)
    line = strip(readline(fstream), ['\n', '\r'])
    sline = split(line, delim)
    pline = parse(sline[1]; raise=false)
    if !isa(pline, Int)
      positionConvert = position(fstream)
      println("Skipped line: $(line)")
    else
      startConvert = true
      fstream = seek(fstream, positionConvert)
      if length(sline) == 2
        fileWeighted = 0
        fileCount, firstNode, lastNode = convertDSVunweighted(fstream, file_dir, file_name, delim, maxMB)
      elseif length(sline) == 3
        fileCount, firstNode, lastNode = convertDSVweighted(fstream, file_dir, file_name, delim, maxMB)
      else
        error("Wrong file format.")
      end
    end
  end
  file = open("$(file_dir)/$(file_name)/info.bin", "w")
  write(file, fileCount, firstNode, lastNode, fileWeighted)
  close(file)
  return fileCount, firstNode, lastNode, fileWeighted
end

function convertDSVweighted(fstream::IOStream, file_dir::String, file_name::String, delim::Char, maxMB::Int)
  fileCount = 0
  firstNode = typemax(Int)
  lastNode = 0
  # hashtable: key = node, value = (degree, neighbours, weights)
  data = Dict{Int, (Float64, Vector{Int}, Vector{Float64})}()
  totalNodes = 0
  maxNeighbours = 0
  while !eof(fstream)
    # parse stream line per line
    line = split(strip(readline(fstream), ['\n', '\r']), delim)
    node = int(line[1])
    neighbour = int(line[2])
    weight = float64(line[3])
    # update the data
    firstNode = min(firstNode, node, neighbour)
    lastNode = max(lastNode, node, neighbour)
    totalNodes, maxNeighbours = updateData!(data, node, neighbour, weight, totalNodes, maxNeighbours)
    totalNodes, maxNeighbours = updateData!(data, neighbour, node, weight, totalNodes, maxNeighbours)
    # check if no more data or size is greater than bound
    if eof(fstream) || checkMemSize(data, maxMB, totalNodes, maxNeighbours)
      # array to map internal to external node index
      nodeArray = collect(keys(data))
      degreeArray = Array(Float64, length(nodeArray))
      # pad neighbour vector with -1's and weight vector with 0's as empty neighbours
      neighbourMatrix = fill(-1, (maxNeighbours, length(nodeArray)))
      weightMatrix = fill(0.0, (maxNeighbours, length(nodeArray)))
      for i in 1:length(nodeArray)
        degreeArray[i] = data[nodeArray[i]][1]
        for j in 1:length(data[nodeArray[i]][2])
          neighbourMatrix[j, i] = data[nodeArray[i]][2][j]
          weightMatrix[j, i] = data[nodeArray[i]][3][j]
        end
      end
      fileCount += 1
      println("Writing components to file set $(fileCount)")
      file = open("$(file_dir)/$(file_name)/nodes$(fileCount).bin", "w")
      write(file, totalNodes, nodeArray)
      close(file)
      file = open("$(file_dir)/$(file_name)/degrees$(fileCount).bin", "w")
      write(file, degreeArray)
      close(file)
      file = open("$(file_dir)/$(file_name)/neighbours$(fileCount).bin", "w")
      write(file, maxNeighbours, neighbourMatrix)
      close(file)
      file = open("$(file_dir)/$(file_name)/weights$(fileCount).bin", "w")
      write(file, weightMatrix)
      close(file)
      empty!(data)
      gc()
      totalNodes = 0
      maxNeighbours = 0
    end
  end
  return fileCount, firstNode, lastNode
end

function convertDSVunweighted(fstream::IOStream, file_dir::String, file_name::String, delim::Char, maxMB::Int)
  fileCount = 0
  firstNode = typemax(Int)
  lastNode = 0
  # hashtable: key = node, value = neighbours
  data = Dict{Int, Vector{Int}}()
  totalNodes = 0
  maxNeighbours = 0
  while !eof(fstream)
    # parse stream line per line
    line = split(strip(readline(fstream), ['\n', '\r']), delim)
    node = int(line[1])
    neighbour = int(line[2])
    # update the data
    firstNode = min(firstNode, node, neighbour)
    lastNode = max(lastNode, node, neighbour)
    totalNodes, maxNeighbours = updateData!(data, node, neighbour, totalNodes, maxNeighbours)
    totalNodes, maxNeighbours = updateData!(data, neighbour, node, totalNodes, maxNeighbours)
    # check if no more data or size is greater than bound
    if eof(fstream) || checkMemSize(data, maxMB, totalNodes, maxNeighbours)
      fileCount += 1
      println("Writing components to file set $(fileCount)")
      # array to map internal to external node index
      nodeArray = collect(keys(data))
      degreeArray = Array(Float64, length(nodeArray))
      # pad neighbour vector with -1's and weight vector with 0's as empty neighbours
      neighbourMatrix = fill(-1, (maxNeighbours, length(nodeArray)))
      for i in 1:length(nodeArray)
        degreeArray[i] = float(length(data[nodeArray[i]]))
        for j in 1:length(data[nodeArray[i]])
          neighbourMatrix[j, i] = data[nodeArray[i]][j]
        end
      end
      file = open("$(file_dir)/$(file_name)/nodes$(fileCount).bin", "w")
      write(file, totalNodes, nodeArray)
      close(file)
      file = open("$(file_dir)/$(file_name)/degrees$(fileCount).bin", "w")
      write(file, degreeArray)
      close(file)
      file = open("$(file_dir)/$(file_name)/neighbours$(fileCount).bin", "w")
      write(file, maxNeighbours, neighbourMatrix)
      close(file)
      empty!(data)
      totalNodes = 0
      maxNeighbours = 0
    end
  end
  return fileCount, firstNode, lastNode
end

"""
Update the data dictionary with the new node data
Return totalNodes and maxNeighbours
"""
function updateData!(data::Dict{Int, (Float64, Vector{Int}, Vector{Float64})}, node::Int, neighbour::Int, weight::Float64, totalNodes::Int, maxNeighbours::Int)
  # check if node already stored
  if haskey(data, node)
    # neighbour not yet stored
    if findfirst(data[node][2], neighbour) == 0
      # update the node degree and add the neighbour
      data[node] = (data[node][1]+weight, push!(data[node][2], neighbour), push!(data[node][3], weight))
      # check if maxNeighbours needs to be increased
      if length(data[node][2]) > maxNeighbours
        maxNeighbours = length(data[node][2])
      end
    end
  else
    data[node] = (weight, [neighbour], [weight])
    totalNodes += 1
  end
  return (totalNodes, maxNeighbours)
end

function updateData!(data::Dict{Int, Vector{Int}}, node::Int, neighbour::Int, totalNodes::Int, maxNeighbours::Int)
  # check if node already stored
  if haskey(data, node)
    # check if neighbour not yet stored
    if findfirst(data[node], neighbour) == 0
      # update the node degree and add the neighbour
      data[node] = push!(data[node], neighbour)
      # check if maxNeighbours needs to be increased
      if length(data[node]) > maxNeighbours
        maxNeighbours = length(data[node])
      end
    end
  else
    data[node] = [neighbour]
    totalNodes += 1
  end
  return (totalNodes, maxNeighbours)
end

"""
Check if the memory footprint is still smaller than the threshold
Return true if threshold reached
"""
function checkMemSize(data::Dict{Int, (Float64, Vector{Int}, Vector{Float64})}, maxMB::Int, totalNodes::Int, maxNeighbours::Int)
  # estimate current needed memory in Bytes: size of dictionary + arrays that need to be generated for storage
  size = sizeof(data) + (totalNodes + totalNodes*maxNeighbours)*sizeof(Int) + totalNodes*sizeof(Float64)
  # return true if size larger than maximum allowed MB threshold
  return (size < maxMB*1000000)?false:true
end

function checkMemSize(data::Dict{Int, Vector{Int}}, maxMB::Int, totalNodes::Int, maxNeighbours::Int)
  # estimate current needed memory in Bytes: size of dictionary + arrays that need to be generated for storage
  size = sizeof(data) + (totalNodes + totalNodes*maxNeighbours)*(sizeof(Int) + sizeof(Float64))
  # return true if size larger than maximum allowed MB threshold
  return (size < maxMB*1000000)?false:true
end

"""
Create an undirected network file with duplicate entries from the original network file
"""
function createMirrored(dsvFile::String, delim::Char, newFile::String)
  # read in dsvfile as stream
  fstream = open(dsvFile)
  # write to newfile as stream
  nstream = open(newFile, "w")
  while !eof(fstream)
    # parse stream line per line
    line = readline(fstream)
    sline = split(strip(line, ['\n', '\r']), delim)
    pline = parse(sline[1]; raise=false)
    if !isa(pline, Int)
      println("Skipped line: $(line)")
    else
      if length(sline) == 2
        node = int(sline[1])
        neighbour = int(sline[2])
        write(nstream, line)
        write(nstream, "$(neighbour)$(delim)$(node)\n")
      elseif length(sline) == 3
        node = int(sline[1])
        neighbour = int(sline[2])
        weight = float(sline[3])
        write(nstream, line)
        write(nstream, "$(neighbour)$(delim)$(node)$(delim)$(weight)\n")
      end
    end
  end
  close(nstream)
  close(fstream)
end
