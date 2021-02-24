#############################################################################
# Copyright (C) 2017 - 2018  Spine Project
#
# This file is part of Spine Model.
#
# Spine Model is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Spine Model is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################
using SpineInterface
using InfrastructureSystems
using PyCall
using URIParser

"""
    psse_to_spine(psse_path, db_url)

Parse the psse raw file (`psse_path`) using PowerSystems.jl and convert it
to a Spine model at `db_url` using `nodes`, `units` and `conenctions`
"""


function psse_to_spine(psse_path, db_url::String)
    pm_data = PowerSystems.parse_psse(psse_path)
    write_powersystem(pm_data, db_url)
end


"""
    write_powersystem!(ps_system, db_url)

Given the PowerSystems.jl system dict `ps_system`, create a Spine Model db at `db_url`
"""

function import_data_structure()
    object_classes = []
    object_parameters = []
    relationship_classes = []
    relationship_parameters = []

    push!(object_parameters,("node", "minimum_voltage"))
    push!(object_parameters,("node", "voltage"))
    push!(object_parameters,("node", "demand"))

    object_classes=["node", "connection", "unit", "commodity"]

    push!(relationship_classes, ("node_group__node", ["node", "node"]))
    push!(relationship_classes, ("connection__to_node", ["connection", "node"]))
    push!(relationship_classes, ("connection__from_node", ["connection", "node"]))
    push!(relationship_classes, ("unit__to_node", ["unit", "node"]))
    push!(relationship_classes, ("node__commodity", ["node", "commodity"]))

    push!(object_parameters, ("commodity", "commodity_lodf_tolerance"))
    push!(object_parameters, ("commodity", "commodity_physics"))
    push!(object_parameters, ("commodity", "commodity_ptdf_flow_tolerance"))
    push!(object_parameters, ("commodity", "commodity_ptdf_threshold"))
    push!(object_parameters, ("commodity", "commodity_slack_penalty"))

    push!(object_parameters, ("connection", "connection_resistance"))
    push!(object_parameters, ("connection", "connection_reactance"))
    push!(object_parameters, ("connection", "connection_monitored"))
    push!(object_parameters, ("connection", "connection_contingency"))
    push!(object_parameters, ("connection", "connection_length"))
    push!(object_parameters, ("connection", "connection_availability_factor"))

    push!(object_parameters, ("unit", "fix_units_on"))
    push!(object_parameters, ("unit", "fom_cost"))
    push!(object_parameters, ("unit", "min_down_time"))
    push!(object_parameters, ("unit", "min_up_time"))
    push!(object_parameters, ("unit", "number_of_units"))
    push!(object_parameters, ("unit", "online_variable_type"))
    push!(object_parameters, ("unit", "shut_down_cost"))
    push!(object_parameters, ("unit", "start_up_cost"))
    push!(object_parameters, ("unit", "unit_availability_factor"))

    push!(relationship_parameters, ("unit__to_node", "unit_capacity"))
    push!(relationship_parameters, ("unit__to_node", "unit_conv_cap_to_flow"))
    push!(relationship_parameters, ("unit__to_node", "minimum_operating_point"))
    push!(relationship_parameters, ("unit__from_node", "unit_capacity"))
    push!(relationship_parameters, ("unit__from_node", "unit_conv_cap_to_flow"))
    push!(relationship_parameters, ("unit__from_node", "minimum_operating_point"))
    push!(relationship_parameters, ("connection__from_node", "connection_capacity"))
    push!(relationship_parameters, ("connection__to_node", "connection_capacity"))
    push!(relationship_parameters, ("connection__from_node", "connection_emergency_capacity"))
    push!(relationship_parameters, ("connection__to_node", "connection_emergency_capacity"))

    return  object_classes,
            object_parameters,
            relationship_classes,
            relationship_parameters
end

function write_powersystem(ps_system::Dict, db_url::String)

    object_classes,
    object_parameters,
    relationship_classes,
    relationship_parameters = import_data_structure()

    objects = []
    relationships = []
    object_parameter_values = []
    relationship_parameter_values = []

    commodity_name="elec"

    baseMVA=ps_system["baseMVA"]

    areas=[]
    zones=[]

    push!(objects, ["commodity", commodity_name])

    push!(object_parameter_values, ("commodity", commodity_name, "commodity_lodf_tolerance", 0.1))
    push!(object_parameter_values,("commodity", commodity_name, "commodity_physics", "commodity_physics_ptdf"))
    push!(object_parameter_values,("commodity", commodity_name, "commodity_ptdf_flow_tolerance",0.1))
    push!(object_parameter_values,("commodity", commodity_name, "commodity_ptdf_threshold",0.0001))
    push!(object_parameter_values,("commodity", commodity_name, "commodity_slack_penalty", 1000))

    node_demand = Dict{Int64,Float64}()
    node_name = Dict{Int64,String}()

    for b in ps_system["bus"]
        data=b[2]
        i=data["bus_i"]
        node_demand[i]=0
        name = string(i, "_" , data["name"][1:3], "_", Int(round(data["base_kv"], digits=0)))
        node_name[i] = name
        n = ("node", name)
        push!(objects, n)
        push!(object_parameter_values, ("node", name, "voltage", data["base_kv"]))

        #if data["bus_type"] == BusTypes.REF
        if data["bus_type"] == 3
            push!(object_parameter_values, ("node", name, "node_opf_type", "node_opf_type_reference"))
        end

        area_name = string("area_", data["area"])
        zone_name = string("zone_", data["zone"])
        if !(area_name in areas)
            push!(areas, area_name)
            push!(objects, ("node", area_name))
            push!(object_parameter_values, ("node", area_name, "minimum_voltage", 110.0))
        end
        obj_list=[]
        push!(obj_list, area_name)
        push!(obj_list, name)
        push!(relationships,("node_group__node",obj_list))

        if !(zone_name in zones)
            push!(zones, zone_name)
            push!(objects, ("node", zone_name))
        end
        obj_list=[]
        push!(obj_list, zone_name)
        push!(obj_list, name)
        push!(relationships,("node_group__node",obj_list))

        obj_list=[]
        push!(obj_list, name)
        push!(obj_list, "elec")
        push!(relationships,("node__commodity", obj_list))
    end

    for b in ps_system["branch"]
        data=b[2]
        if data["br_status"] in (0, 1)
            from_bus_name = node_name[data["f_bus"]]
            to_bus_name =  node_name[data["t_bus"]]
            ckt = rstrip(string(data["source_id"][4]))
            name = string(from_bus_name, "__", to_bus_name, "_", ckt)
            connection = ("connection", name)
            push!(objects, connection)
            push!(object_parameter_values, ("connection", name, "connection_resistance", data["br_r"]))
            push!(object_parameter_values, ("connection", name, "connection_reactance", data["br_x"]))
            push!(object_parameter_values, ("connection", name, "connection_monitored", 1))
            push!(object_parameter_values, ("connection", name, "connection_contingency", 1))
            #push!(object_parameter_values, ("connection", name, "connection_length", data["br_r"]))
            push!(object_parameter_values, ("connection", name, "connection_availability_factor", 1.0))
            obj_list=[]
            push!(obj_list,name)
            push!(obj_list,to_bus_name)
            push!(relationships,("connection__to_node",obj_list))
            rate_a = round(data["rate_a"] * baseMVA, digits=2)
            if haskey(data, "rate_b")
                rate_b = round(data["rate_b"] * baseMVA, digits=2)
            else
                rate_b = rate_a
            end

            push!(relationship_parameter_values,("connection__to_node", obj_list, "connection_capacity", rate_a))
            push!(relationship_parameter_values,("connection__to_node", obj_list, "connection_emergency_capacity", rate_b))

            obj_list=[]
            push!(obj_list,name)
            push!(obj_list,from_bus_name)
            push!(relationships,("connection__from_node",obj_list))
            push!(relationship_parameter_values,("connection__from_node", obj_list, "connection_capacity", rate_a))
            push!(relationship_parameter_values,("connection__from_node", obj_list, "connection_emergency_capacity", rate_b))
        end
    end

    for dc in ps_system["dcline"]
        data=dc[2]
        from_bus_name = node_name[data["source_id"][2]]
        to_bus_name = node_name[data["source_id"][3]]
        ckt = 1
        name = string(from_bus_name, "__", to_bus_name, "_", ckt)
        connection = ("connection", name)
        pmaxt=round(data["pmaxt"]  * baseMVA, digits=2)
        pmaxf=round(data["pmaxf"]  * baseMVA, digits=2)

        push!(objects, connection)
        push!(object_parameter_values, ("connection", name, "connection_resistance", 0))
        push!(object_parameter_values, ("connection", name, "connection_reactance", 0.0001))
        push!(object_parameter_values, ("connection", name, "connection_monitored", 0))
        push!(object_parameter_values, ("connection", name, "connection_contingency", 0))
        #push!(object_parameter_values, ("connection", name, "connection_length", data["br_r"]))
        push!(object_parameter_values, ("connection", name, "connection_availability_factor", 1.0))

        obj_list=[]
        push!(obj_list,name)
        push!(obj_list,to_bus_name)
        push!(relationships,("connection__to_node",obj_list))
        push!(relationship_parameter_values,("connection__to_node", obj_list, "connection_capacity", pmaxt))
        push!(relationship_parameter_values,("connection__to_node", obj_list, "connection_emergency_capacity", pmaxt))

        obj_list=[]
        push!(obj_list,name)
        push!(obj_list,from_bus_name)
        push!(relationships,("connection__from_node",obj_list))
        push!(relationship_parameter_values,("connection__from_node", obj_list, "connection_capacity", pmaxf))
        push!(relationship_parameter_values,("connection__from_node", obj_list, "connection_emergency_capacity", pmaxf))
    end

    gen_ids=[]
    for g in ps_system["gen"]
        data=g[2]
        bus_name=node_name[data["gen_bus"]]
        name = string(ps_system["bus"][string(data["gen_bus"])]["name"][1:3], "_", rstrip(data["source_id"][3]))
        if name in gen_ids
            sub_index = 1
            new_name = string(name, "_", sub_index)
            while new_name in gen_ids
                sub_index = sub_index + 1
                new_name = string(name, "_", sub_index)
            end
            name = new_name
        end
        push!(gen_ids, name)

        unit = ("unit", name)
        push!(objects, unit)
        push!(object_parameter_values, ("unit", name, "number_of_units", 1))
        push!(object_parameter_values, ("unit", name, "online_variable_type", "unit_online_variable_type_binary"))
        push!(object_parameter_values, ("unit", name, "unit_availability_factor", 1))

        obj_list=[]
        push!(obj_list,name)
        push!(obj_list,bus_name)
        push!(relationships,("unit__to_node",obj_list))

        push!(relationship_parameter_values,("unit__to_node", obj_list, "unit_capacity",round(data["pmax"]*baseMVA, digits=2)))
        pmin=round(data["pmin"]/data["pmax"], digits=4)
        if isa(pmin, Number)
            push!(relationship_parameter_values,("unit__to_node", obj_list, "minimum_operating_point", pmin))
        end
    end

    for l in ps_system["load"]
        data = l[2]
        if data["status"] == 1
            node_demand[data["load_bus"]] = node_demand[data["load_bus"]] + data["pd"]
        end
    end

    for b in ps_system["bus"]
        data = b[2]
        name = node_name[data["bus_i"]]
        load_contender = round(baseMVA * node_demand[data["bus_i"]], digits=4)
        if load_contender > 0
            push!(object_parameter_values, ("node", name, "demand", load_contender))
        end
    end
@info "writing PSSE data to $(db_url)"

db_map=db_api.DiffDatabaseMapping(db_url; upgrade=true)

added, err_log = db_api.import_data(
    db_map,
    object_classes=object_classes,
    relationship_classes=relationship_classes,
    object_parameters=object_parameters,
    relationship_parameters=relationship_parameters,
    objects=objects,
    relationships=relationships,
    object_parameter_values=object_parameter_values,
    relationship_parameter_values=relationship_parameter_values,
)
comment="powersystems to spine import"
db_map.commit_session(comment)

@info "data imported to $(db_url)"

end

function aggregate_network(db_url::String)
    con__mon = Tuple{Object,Object}[] # this is a set of monitored and contingent line tuples that must be considered as defined by the connection_monitored and connection_contingency parmaeters
    monitored_lines=[]
    ptdf_conn_n = Dict{Tuple{Object,Object},Float64}() #ptdfs returned by PowerSystems.jl
    lodf_con_mon = Dict{Tuple{Object,Object},Float64}() #lodfs calcuated based on ptdfs returned by PowerSystems.jl
    net_inj_nodes=[] # this is the set of nodes with demand or generation

    object_classes,
    object_parameters,
    relationship_classes,
    relationship_parameters = import_data_structure()

    objects = []
    object_parameter_values = []
    relationships = []
    relationship_parameter_values = []

    using_spinedb(db_url)

    @info "pruning system at $(db_url)"

    for c in commodity()
        if commodity_physics(commodity=c) in (:commodity_physics_ptdf, :commodity_physics_lodf)

            @info "Processing network for commodity $(c) with network_physics $(commodity_physics(commodity=c))"
            n_islands, island_node = islands()
            @info "Your network consists of $(n_islands) islands"
            if n_islands > 1
                @warn "Your network consists of multiple islands, this may end badly."
                print(island_node)
            end

            net_inj_nodes=get_net_inj_nodes() # returns list of nodes that have demand and/or generation

            @info "calculating ptdfs"
            ptdf_conn_n = calculate_ptdfs()

            if commodity_physics(commodity = c) == :commodity_physics_lodf
                con__mon = Tuple{Object,Object}[]
                @info "calculating lodfs"  lodf_con_mon = calculate_lodfs(ptdf_conn_n, con__mon)
            end
        end
    end

    inj_nodes = []
    comm_nodes = []

    traversed = Dict{Object,Bool}()
    node__new_nodes = Dict{Object,Array}()
    min_v = 0
    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                push!(comm_nodes, n)
                if !isempty(unit__to_node(node=n)) || demand(node=n) > 0
                    for ng in node_group__node(node2=n)
                        if ! (minimum_voltage(node=ng) == nothing)
                            min_v = minimum_voltage(node=ng)
                            break
                        end
                    end
                    if voltage(node=n) < min_v
                        push!(inj_nodes, n)
                    end
                end
            end
        end
    end

    #@info "Writing ptdf diagnostic file"
    #write_ptdfs(ptdf_conn_n, comm_nodes)

    @info "Traversing nodes"
    for n in inj_nodes
        node__new_nodes[n]=[]
    end

    for n in inj_nodes
        for n2 in comm_nodes
            traversed[n2] = false
        end
        min_v=0
        for ng in node_group__node(node2=n)
            if ! (minimum_voltage(node=ng) == nothing)
                min_v = minimum_voltage(node=ng)
                break
            end
        end
        traverse(n, n, traversed, node__new_nodes, min_v, ptdf_conn_n, 1)
    end

    to_prune_object_ids = []
    nodes_pruned = 0
    connections_pruned=0
    units_moved = 0
    units_distributed = 0
    demands_moved = 0
    demands_distributed = 0

    for n in comm_nodes
        for ng in node_group__node(node2=n)
            if ! (minimum_voltage(node=ng) == nothing)
                min_v = minimum_voltage(node=ng)
                if voltage(node=n) < min_v
                    push!(to_prune_object_ids, n.id)
                    nodes_pruned += 1
                    for conn in connection__to_node(node=n)
                        if !(conn.id in to_prune_object_ids)
                            push!(to_prune_object_ids, conn.id)
                            connections_pruned += 1
                        end
                    end
                    for conn in connection__from_node(node=n)
                        if !(conn.id in to_prune_object_ids)
                            push!(to_prune_object_ids, conn.id)
                            connections_pruned += 1
                        end
                    end
                end
                break
            end
        end
    end

    # line below generates the error `using_spinedb(db_url)`
    # notusing_spinedb(db_url)

    @info "Writing output"
    write_node__new_nodes(node__new_nodes)

    new_demand_dict = Dict{Object,Float64}()
    new_gen_dict = Dict{Object,Float64}()
    gens_to_move = Dict{Object,Object}()

    for (n, data) in node__new_nodes
        if demand(node=n) == nothing
            demand_to_shift = 0
        else
            demand_to_shift = demand(node=n)
        end
        if size(data, 1) == 1 # only one connected higher voltage node, move all the demand here
            if demand_to_shift > 0
                if haskey(new_demand_dict, data[1])
                    new_demand_dict[data[1][1]] = new_demand_dict[data[1][1]] + demand_to_shift
                else
                    new_demand_dict[data[1][1]] = demand_to_shift
                end
                demands_moved += 1
            end
            for u in unit__to_node(node=n)
                gens_to_move[u] = data[1][1]
                units_moved += 1
            end
        else
            gen_to_shift = 0
            for u in unit__to_node(node=n)
                gen_to_shift = gen_to_shift + unit_capacity(unit=u, node=n)
                push!(to_prune_object_ids, u.id)
                units_distributed += 1
            end
            for (n2, ptdf) in data
                ptdf = abs(ptdf)
                if demand_to_shift > 0
                    if haskey(new_demand_dict, n2)
                        new_demand_dict[n2] = new_demand_dict[n2] + demand_to_shift * ptdf
                    else
                        new_demand_dict[n2] = demand_to_shift * ptdf
                    end
                    demands_distributed += 1
                end
                if haskey(new_gen_dict, n2)
                    new_gen_dict[n2] = new_gen_dict[n2] + gen_to_shift * ptdf
                else
                    new_gen_dict[n2] = gen_to_shift * ptdf
                end
            end
        end
    end

    # Update demand parameter of higher voltage nodes to add demand of pruned nodes

    for (n, new_demand) in new_demand_dict
        if demand(node=n) == nothing
            updated_demand = new_demand
        else
            updated_demand = demand(node=n) + new_demand
        end
        push!(object_parameter_values, ("node", string(n), "demand", new_demand))
    end

    for (u, new_node) in gens_to_move
        obj_list=[]
        push!(obj_list,string(u))
        push!(obj_list,string(new_node))
        push!(relationships,("unit__to_node",obj_list))
        for old_node in unit__to_node(unit=u)
            push!(relationship_parameter_values,("unit__to_node", obj_list, "unit_capacity", unit_capacity(unit=u, node=old_node)))
            if !(minimum_operating_point(unit=u, node=old_node) == nothing)
                push!(relationship_parameter_values,("unit__to_node", obj_list, "minimum_operating_point", minimum_operating_point(unit=u, node=old_node)))
            end
        end
    end

    for (n, new_gen) in new_gen_dict
        if new_gen > 0
            unit_name = string("U_DIST_", n)
            push!(objects, ("unit", unit_name))
            push!(object_parameter_values, ("unit", unit_name, "number_of_units", 1))
            push!(object_parameter_values, ("unit", unit_name, "online_variable_type", "unit_online_variable_type_binary"))
            push!(object_parameter_values, ("unit", unit_name, "unit_availability_factor", 1))

            obj_list=[]
            push!(obj_list,unit_name)
            push!(obj_list,string(n))
            push!(relationships,("unit__to_node",obj_list))
            push!(relationship_parameter_values,("unit__to_node", obj_list, "unit_capacity", new_gen))
        end
    end

    db_uri = URI(db_url)
    db_path = db_uri.path[2:length(db_uri.path)]
    new_db_path_root = string(db_path[1:findlast(isequal('.'), db_path) - 1], "_pruned")
    path_ext = ".sqlite"
    if isfile(string(new_db_path_root, path_ext))
        i = 1
        while isfile(string(new_db_path_root, "_", i, path_ext))
            i = i + 1
        end
        new_db_path = string(new_db_path_root, "_", i, path_ext)
    else
        new_db_path = string(new_db_path_root, path_ext)
    end

    cp(db_path, new_db_path)
    new_db_url=string("sqlite:///", new_db_path)

    @info "new database copied to $(new_db_path)"

    db_map=db_api.DiffDatabaseMapping(new_db_url; upgrade=true)

    db_map.remove_items(object_ids=to_prune_object_ids)
    comment="network pruning demand and generation shifts"
    db_map.commit_session(comment)

    added, err_log = db_api.import_data(
        db_map,
        object_classes=object_classes,
        relationship_classes=relationship_classes,
        object_parameters=object_parameters,
        relationship_parameters=relationship_parameters,
        objects=objects,
        relationships=relationships,
        object_parameter_values=object_parameter_values,
        relationship_parameter_values=relationship_parameter_values,
    )
    @info "added $(added) items"
    for err in err_log
        @info "import error: " err.msg
    end

    comment="network pruning demand and generation shifts"
    db_map.commit_session(comment)

    @info "network pruned successfully" nodes_pruned connections_pruned units_moved units_distributed demands_moved demands_distributed
end

function traverse(n_t, n, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf)
    #@info "traversing $(n) for $(n_t) with voltage $(voltage(node=n)) and min voltage $(min_v)"
    traversed[n] = true
    if voltage(node=n) >= min_v
        push!(node__new_nodes[n_t], (n, ptdf) )
    else
        for conn in connection__from_node(node=n)
            for n2 in connection__to_node(connection=conn)
                if !traversed[n2]
                    traverse(n_t, n2, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf_conn_n[(conn,n_t)] )
                end
            end
        end
        for conn in connection__to_node(node=n)
            for n2 in connection__from_node(connection=conn)
                if !traversed[n2]
                    traverse(n_t, n2, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf_conn_n[(conn,n_t)])
                end
            end
        end
    end
end

function write_node__new_nodes(node__new_nodes)
    io = open("node_mapping.csv", "w")
    print(io, "node,V,total_gens,total_generation,total_demand,node1,V1,ptdf1,node2,V2,ptdf2,node3,V3,ptdf3,node4,V4,ptdf4,node5,V5,ptdf5\n")
    for n in node__new_nodes
        total_generators=0
        total_generation=0
        total_demand=0
        for u in unit__to_node(node=n[1])
            total_generators = total_generators + 1
            total_generation = total_generation + unit_capacity(unit=u, node=n[1])
        end
        if demand(node=n[1]) == nothing
            total_demand = 0
        else
            total_demand = demand(node=n[1])
        end
        print(io, string(n[1]), ",", voltage(node=n[1]), ",", total_generators, ",", total_generation, ",", total_demand)
        if size(n[2],1) == 1
            print(io, ",", string(n[2][1][1]), ",", voltage(node=n[2][1][1]), ",", 1)
        else
            for n2 in n[2]
                print(io, ",", string(n2[1]), ",", voltage(node=n2[1]), ",", string(abs(n2[2])))
            end
        end
        print(io, "\n")
    end
    close(io)
end


"""
    islands()

Determines the number of islands in a commodity network - used for diagnostic purposes

"""
function islands()
    visited_d = Dict{Object,Bool}()
    island_node = Dict{Int64,Array}()
    island = 0

    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                visited_d[n] = false
            end
        end
    end

    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                if !visited_d[n]
                    island = island + 1
                    island_node[island] = Object[]
                    visit(n, island, visited_d, island_node)
                end
            end
        end
    end
    return island, island_node
end

"""
    visit()

Function called recursively to visit nodes in the network to determine number of islands

"""
function visit(n, island, visited_d, island_node)
    visited_d[n] = true
    push!(island_node[island], n)
    for conn in connection__from_node(node=n)
        for n2 in connection__to_node(connection=conn)
            if !visited_d[n2]
                visit(n2, island, visited_d, island_node)
            end
        end
    end
    for conn in connection__to_node(node=n)
        for n2 in connection__from_node(connection=conn)
            if !visited_d[n2]
                visit(n2, island, visited_d, island_node)
            end
        end
    end
end

"""
    check_x()

Check for low reactance values

"""
function check_x()
    @info "Checking reactances"
    for conn in connection()
        if conn_reactance(connection=conn) < 0.0001
            @info "Low reactance may cause problems for line " conn
        end
    end
end

"""
    calculate_ptdfs()

Returns a dict indexed on tuples of (connection, node) containing the ptdfs of the system currently in memory.

"""
function calculate_ptdfs()
    ps_busses = Bus[]
    ps_lines = Line[]

    node_ps_bus = Dict{Object,Bus}()
    i = 1
    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                if node_opf_type(node=n) == :node_opf_type_reference
                    bustype = BusTypes.REF
                else
                    bustype = BusTypes.PV
                end
                ps_bus = Bus(
                    number = i,
                    name = string(n),
                    bustype = bustype,
                    angle = 0.0,
                    voltage = 0.0,
                    voltagelimits = (min = 0.0, max = 0.0),
                    basevoltage = nothing,
                    area = nothing,
                    load_zone = LoadZone(nothing),
                    ext = Dict{String, Any}()
                )

                push!(ps_busses,ps_bus)
                node_ps_bus[n] = ps_bus
                i = i + 1
            end

            PowerSystems.buscheck(ps_busses)
            PowerSystems.slack_bus_check(ps_busses)

            for conn in connection()
                for n_from in connection__from_node(connection=conn)
                    for n_to in connection__to_node(connection=conn)
                        for c in node__commodity(node=n_from)
                            if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
                                ps_arc = Arc(node_ps_bus[n_from], node_ps_bus[n_to])
                                new_line = Line(;
                                    name = string(conn),
                                    available = true,
                                    activepower_flow = 0.0,
                                    reactivepower_flow = 0.0,
                                    arc = ps_arc,
                                    r = connection_resistance(connection=conn),
                                    x = max(connection_reactance(connection=conn), 0.00001),
                                    b = (from=0.0, to=0.0),
                                    rate = 0.0,
                                    anglelimits = (min = 0.0, max = 0.0)
                                )
                                push!(ps_lines,new_line)
                            end  #in case there are somehow multiple commodities
                            break
                        end

                    end
                end
            end
        end
    end

    ps_ptdf=PowerSystems.PTDF(ps_lines,ps_busses)

    ptdf=Dict{Tuple{Object,Object},Float64}()

    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                for conn in connection()
                    ptdf[(conn,n)] = ps_ptdf[string(conn), node_ps_bus[n].number]
                end
            end
        end
    end

    return ptdf

    #buildlodf needs to be updated to account for cases
    #lodfs=PowerSystems.buildlodf(ps_lines,ps_busses)
end

#function (ptdf::PTDF)(conn::Object, n::Object)
#     # Do something with ptdf, conn, and n, and return the value
#
#     return ptdf[string(conn),string(n)]
#end



"""
    calculate_lodfs(ptdf_b_n, b_con__b_mon)

Returns lodfs for the system specified by ptdf_b_n ,b_con__b_mon as a dict of tuples: contingent_line, monitored_line

"""
# This function takes a long time. PowerSystems has a function that does it faster using linear algebra but doesn't handle the case of tails like
# I would like.

function calculate_lodfs(ptdf_conn_n, con__mon)
    lodf_con_mon = Dict{Tuple{Object,Object},Float64}()
    considered_contingencies = 0
    skipped = 0
    tolerance = 0
    for conn_con in connection()
        if connection_contingency(connection = conn_con) == 1
            for n_from in connection__from_node(connection=conn_con)
                for n_to in connection__to_node(connection=conn_con)
                    demoninator = 1 - (ptdf_conn_n[(conn_con, n_from)]-ptdf_conn_n[(conn_con, n_to)])
                    if abs(demoninator) < 0.001
                        demoninator = -1
                    end
                    for conn_mon in connection()
                        if connection_monitored(connection=conn_mon) == 1 && conn_con != conn_mon
                            if demoninator == -1
                                lodf_trial = (ptdf_conn_n[(conn_mon, n_from)] - ptdf_conn_n[(conn_mon, n_to)]) / demoninator
                            else
                                lodf_trial = -ptdf_conn_n[(conn_mon, n_from)]
                            end
                            c = first(indices(commodity_lodf_tolerance))
                            tolerance = commodity_lodf_tolerance(commodity=c)
                            if abs(lodf_trial) > tolerance
                                considered_contingencies = considered_contingencies + 1
                                push!(con__mon, (conn_con, conn_mon))
                                lodf_con_mon[(conn_con, conn_mon)] = lodf_trial
                            else
                                skipped = skipped + 1
                            end
                        end
                    end
                end
            end
        end
    end
    #@info "Contingencies summary " considered_contingencies skipped tolerance
    return lodf_con_mon
end

function get_net_inj_nodes()
    net_inj_nodes = []
    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                for u in unit__to_node(node=n)
                    if !(n in net_inj_nodes)
                        push!(net_inj_nodes, n)
                    end
                end
                for u in unit__from_node(node=n)
                    if !(n in net_inj_nodes)
                        push!(net_inj_nodes, n)
                    end
                end
                for ng in node_group__node(node2=n)
                    if fractional_demand(node1=ng, node2=n) > 0 || demand(node=n) > 0
                        if !(n in net_inj_nodes)
                            push!(net_inj_nodes, n)
                        end
                    end
                end
            end
        end
    end
    return net_inj_nodes
end


function write_ptdfs(ptdfs, net_inj_nodes)
    io = open("ptdfs.csv", "w")
    print(io, "connection,")
    for n in net_inj_nodes
        print(io, string(n), ",")
    end
    print(io, "\n")
    for conn in connection()
        print(io, string(conn), ",")
        for n in net_inj_nodes
            if haskey(ptdfs, (conn,n))
                print(io, ptdfs[(conn,n)], ",")
            else
                print(io, "_nan", ",")
            end
        end
        print(io, "\n")
    end
    close(io)
end
