using Assortativity, NetworkInference, Test

testdir = dirname(Base.source_path())
datadir = joinpath(joinpath(dirname(testdir), "examples"), "data")
data_filename = "data.csv"
test_net = infer_network(joinpath(datadir, data_filename), PIDCNetworkInference(), delim = ',')
labels_to_groups = get_labels_to_groups(test_net.nodes, joinpath(datadir, "groups.tsv"))
groups_to_indices = get_groups_to_indices(labels_to_groups)

@test typeof(test_net) == InferredNetwork
@test typeof(assortativity(test_net)) == AssortativityObject

@test assortativity(test_net, labels_to_groups, groups_to_indices, 10).value == assortativity(set_threshold(test_net, (1:10)), labels_to_groups, groups_to_indices).value
