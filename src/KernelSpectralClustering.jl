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
export coverage, modularity1, modularity3
export convertDSVNet, createMirrored

function runKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int; sizeMB=4000::Int, eval=false::Bool)
  filedir = dirname(networkfile)
  filename = splitext(basename(networkfile))[1]
  if sizeMB > 0
    # kscbignet algo
    qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted = @time kscbignet(networkfile, filedir, filename, delimiter, sizeMB; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false, convertF=true)
    println("$(numFiles) files with training subset size $(length(trainSubset)) and BAF value $(baf[k-mink+1]) for $(k) clusters")
    if eval
      println("Evaluation step")
      metricsfile = open("$(filedir)/$(filename)/metricsKSC.txt", "w")
      # Coverage
      tic()
      coverageM = coverage(numFiles, trainSubset, lastNode)
      timeCoverage = toq()
      println("Coverage $(coverageM) finished, time elapsed $(timeCoverage)s")
      write(metricsfile, "Coverage: value $(coverageM), time $(timeCoverage)s\n")
      # Modularity
      tic()
      modularityM = modularity3(numFiles, k, fileWeighted)
      timeModularity = toq()
      println("Modularity $(modularityM) finished, time elapsed $(timeModularity)s")
      write(metricsfile, "Modularity: value $(modularityM), time $(timeModularity)s\n")
      close(metricsfile)
    end
    return qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted
  else
    # kscnet algo
    result = @time kscnet(networkfile, filedir, filename, delimiter; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false)
    return result
  end
end

function rerunKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int; sizeMB=4000::Int, eval=false::Bool)
  filedir = dirname(networkfile)
  filename = splitext(basename(networkfile))[1]
  if sizeMB > 0
    # kscbignet algo
    qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted = @time kscbignet(networkfile, filedir, filename, delimiter, sizeMB; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false, convertF=false)
    println("KSCbignet: $(numFiles) files with training subset size $(length(trainSubset)) and BAF value $(baf[k-mink+1]) for $(k) clusters")
    if eval
      println("Evaluation step")
      metricsfile = open("$(filedir)/$(filename)/metricsKSC.txt", "w")
      # Coverage
      tic()
      coverageM = coverage(numFiles, trainSubset, lastNode)
      timeCoverage = toq()
      println("Coverage $(coverageM) finished, time elapsed $(timeCoverage)s")
      write(metricsfile, "Coverage: value $(coverageM), time $(timeCoverage)s\n")
      # Modularity
      tic()
      modularityM = modularity3(fileCount, kBAF, fileWeighted)
      timeModularity = toq()
      println("Modularity $(modularityM) finished, time elapsed $(timeModularity)s")
      write(metricsfile, "Modularity: value $(modularityM), time $(timeModularity)s\n")
      close(metricsfile)
    end
    return qTest, baf, k, firstNode, lastNode, trainSubset, numFiles, fileWeighted
  else
    # kscnet algo
    result = @time kscnet(networkfile, filedir, filename, delimiter; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false)
    return result
  end
end
end # module
