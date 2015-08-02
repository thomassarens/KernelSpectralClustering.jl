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
export createMirrored

function runKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int; sizeMB=4000::Int, eval=false::Bool)
  filedir = dirname(networkfile)
  filename = splitext(basename(networkfile))[1]
  if sizeMB > 0
    # kscbignet algo
    if eval
      qTest, baf, k, sizeSubset, coverageM, modularityM = @time kscbignet(networkfile, filedir, filename, delimiter, sizeMB; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false, metrics=true, convertF=true)
      println("Coverage $(coverageM) for $(sizeSubset) nodes, Modularity $(modularityM) and BAF $(baf[k-mink+1]) for $(k) clusters")
      return qTest, baf, k, coverageM, modularityM
    else
      qTest, baf, k, sizeSubset, networkData, numFiles = @time kscbignet(networkfile, filedir, filename, delimiter, sizeMB; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false, metrics=false, convertF=true)
      println("$(numFiles) files with training subset size $(sizeSubset) and BAF value $(baf[k-mink+1]) for $(k) clusters")
      return qTest, baf, k, sizeSubset, networkData, numFiles
    end
  else
    # kscnet algo
    result = @time kscnet(networkfile, filedir, filename, delimiter; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false, metrics=false)
    return result
  end
end

function rerunKSCnetwork(networkfile::String, delimiter::Char, mink::Int, maxk::Int; sizeMB=4000::Int, eval=false::Bool)
  filedir = dirname(networkfile)
  filename = splitext(basename(networkfile))[1]
  if sizeMB > 0
    # kscbignet algo
    if eval
      qTest, baf, k, sizeSubset, coverageM, modularityM = @time kscbignet(networkfile, filedir, filename, delimiter, sizeMB; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false, metrics=true, convertF=false)
      println("Coverage $(coverageM) for $(sizeSubset) nodes, Modularity $(modularityM) and BAF $(baf[k-mink+1]) for $(k) clusters")
      return qTest, baf, k, coverageM, modularityM
    else
      qTest, baf, k, sizeSubset, networkData, numFiles = @time kscbignet(networkfile, filedir, filename, delimiter, sizeMB; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false, metrics=false, convertF=false)
      println("$(numFiles) files with training subset size $(sizeSubset) and BAF value $(baf[k-mink+1]) for $(k) clusters")
      return qTest, baf, k, sizeSubset, networkData, numFiles
    end
  else
    # kscnet algo
    result = @time kscnet(networkfile, filedir, filename, delimiter; fraction=0.15, method="furs", minK=mink, maxK=maxk, eigfull=false, metrics=false)
    return result
  end
end
end # module
