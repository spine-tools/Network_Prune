## Network_Prune.jl


This repository contains two main elements:
 - psse_to_spine.jl : a tool to parse PSSE raw files using Powersystems.jl and create a minimal SpineOpt model
 - prune_network.jl : a tool to reduce a network model by voltage level. The algorithm calculates the power transfer distribution factors (PTDFs) of the network and uses tree traversal to find the nearest higher-voltage level node above the minimum voltage level specified. When only a single  higher level node is identified, generation and loads are simply moved to the nearest higher level node. Whem multiple paths exist to a higher voltage level, the network PTDFs are used to distribute the load and generation in the correct proportion. This approach will yield an identical power flow solution to the original network but will contain fewer nodes and connections - the trade-off is that the flows on the lower level network are unknown.

 