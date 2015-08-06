module KernelSpectralClustering
# core files
include("KSCdata.jl")
include("KSCselect.jl")
include("KSCmodel.jl")
include("KSCevaluation.jl")
# algorithm implementation files
include("KSCnet.jl")
include("KSCbignet.jl")

export runKSCnetwork, rerunKSCnetwork
export kscnet, kscbignet
export coverage, modularity, modularityApprox
export convertDSVNet, createMirrored

function runKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int; sizeMB=4000::Int, undirected=true::Bool, eval=false::Bool)
  filedir = dirname(networkfile)
  filename = splitext(basename(networkfile))[1]
  if !undirected
    createMirrored(networkfile, delimiter, "$(filedir)/$(filename)_mirror.txt")
    networkfile = "$(filedir)/$(filename)_mirror.txt"
  end
  if sizeMB > 0
    # kscbignet algo
    qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted = @time kscbignet(networkfile, filedir, filename, delimiter, sizeMB; fraction=0.15, method="furs", decomp="eigs", minK=mink, maxK=maxk, convertF=true)
    println("$(numFiles) files with training subset size $(length(trainSubset)) and BAF value $(baf[k-mink+1]) for $(k) clusters")
    if eval
      println("Evaluation step")
      metricsfile = open("$(filedir)/$(filename)/metricsKSC.txt", "w")
      # Coverage
      tic()
      coverageM = coverage(numFiles, trainSubset, length(unique(qTest)))
      timeCoverage = toq()
      println("Coverage $(coverageM) finished, time elapsed $(timeCoverage)s")
      write(metricsfile, "Coverage: value $(coverageM), time $(timeCoverage)s\n")
      # Modularity
      tic()
      modularityM = modularityApprox(numFiles, k, fileWeighted)
      timeModularity = toq()
      println("Modularity $(modularityM) finished, time elapsed $(timeModularity)s")
      write(metricsfile, "Modularity: value $(modularityM), time $(timeModularity)s\n")
      close(metricsfile)
    end
    return qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted
  else
    # kscnet algo
    result = @time kscnet()
    return result
  end
end

function rerunKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int; sizeMB=4000::Int, undirected=true::Bool, eval=false::Bool)
  filedir = dirname(networkfile)
  filename = splitext(basename(networkfile))[1]
  if !undirected
    createMirrored(networkfile, delimiter, "$(filedir)$(filename)_mirror.txt")
    networkfile = "$(filedir)$(filename)_mirror.txt"
  end
  if sizeMB > 0
    # kscbignet algo
    qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted = @time kscbignet(networkfile, filedir, filename, delimiter, sizeMB; fraction=0.15, method="furs", decomp="eigs", minK=mink, maxK=maxk, convertF=false)
    println("KSCbignet: $(numFiles) files with training subset size $(length(trainSubset)) and BAF value $(baf[k-mink+1]) for $(k) clusters")
    if eval
      println("Evaluation step")
      metricsfile = open("$(filedir)/$(filename)/metricsKSC.txt", "w")
      # Coverage
      tic()
      coverageM = coverage(numFiles, trainSubset, length(unique(qTest)))
      timeCoverage = toq()
      println("Coverage $(coverageM) finished, time elapsed $(timeCoverage)s")
      write(metricsfile, "Coverage: value $(coverageM), time $(timeCoverage)s\n")
      # Modularity
      tic()
      modularityM = modularityApprox(numFiles, k, fileWeighted)
      timeModularity = toq()
      println("Modularity $(modularityM) finished, time elapsed $(timeModularity)s")
      write(metricsfile, "Modularity: value $(modularityM), time $(timeModularity)s\n")
      close(metricsfile)
    end
    return qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted
  else
    # kscnet algo
    result = @time kscnet()
    return result
  end
end
end # module
