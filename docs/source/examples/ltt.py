import dendropy

tree = dendropy.Tree.get(
            path="hiv1.nexus",
            schema="nexus")

# Returns distance of node furthest from root, i.e., maximum time available on
# tree
total_time = tree.max_distance_from_root()

# Divide time span into 10 steps
step = float(total_time) / 10

# To store tuples of (time, number of lineages)
ltt = []

# Start at first time step
current_time = step
while current_time <= total_time:
    # Get number of lineages at current time
    num_lineages = tree.num_lineages_at(current_time)
    # Store it
    ltt.append( (current_time, num_lineages) )
    # Move to next time step
    current_time += step

# Get the final number of lineages
# Note: may not be the same as the number of tips if the tree has extinct
# tips/taxa; though, if this were the case, we would not be dealing with an
# ultrametric tree.
if current_time < total_time:
    ltt.append( (total_time, tree.num_lineages_at(total_time)) )

# Print results
for t, num_lineages in ltt:
    print("{:12.8f}\t{}".format(t, num_lineages))
