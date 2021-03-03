## Network_Prune.jl


This repository contains two main elements:
 - prune_network.jl : a tool to reduce a network model by voltage level.
 - psse_to_spine.jl : a tool to parse PSSE raw files using Powersystems.jl and create a minimal SpineOpt model
  
 ## prune_network.jl 
 The algorithm calculates the power transfer distribution factors (PTDFs) of the network and uses tree traversal to find the nearest higher-voltage level node above the minimum voltage level specified. When only a single  higher level node is identified, generation and loads are simply moved to the nearest higher level node. Whem multiple paths exist to a higher voltage level, the network PTDFs are used to distribute the load and generation in the correct proportion. This approach will yield an identical power flow solution to the original network but will contain fewer nodes and connections - the trade-off is that the flows on the lower level network are unknown.

 ### usage
To run the tool on a Spine database located at `db_url`, simply call :

```julia
db_url="sqlite:///powersystems_test.sqlite"
include("prune_network.jl")
prune_network(db_url)
```

Requirements for the input Spine database
 - The network is defined by nodes and connections. It is assumed each connection has a single `connection__from_node` and a single `connection__to_node relationship`
 - To be included in the network, each node must be related to a single `commodity` with `commodity_physics(commodity=c) in(:commodity_physics_ptdf, :commodity_physics_lodf)`


## prune_network.jl 