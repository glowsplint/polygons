# Approach One: Publish to all relevant channels

Each vertex has its own buffered channel. When the buffered channel is full, the goroutine kicks off.
At the end of the goroutine, it writes to all buffered channels that requires its result.

The goroutine needs to take in:

1. Its input buffered channel where it will read from
2. The slice of output buffered channels it needs to write to

The benefit is that we don't have to topo-sort ahead of time, everything just works as-is.
You can delegate the responsibility of sending the message to the primary goroutine in a pub-sub model.
The pub-sub model can allow for symmetry to be handled inside the primary goroutine.

# Approach Two: Topological sort

1. Topological sort on the nodes
2. Create a queue with these nodes
3. Create a pool of worker nodes
4. When a task is done, the goroutine sends its output on a single channel.
5. The queue processor maintains a single channel and stores the output. It spins up new goroutines when they can be calculated.
6. An indegrees array is maintained to trigger the spinning up of new goroutines.

# Optimisations

1. We can discard the values of earlier nodes once they are no longer required
2. Function-level optimisation to ensure we allocate on the stack instead of the heap.
3. Using goroutines to speed up processing by saturating available cores
4. Using bottom-up dynamic programming and avoiding the recursion limit
5. Using symmetry to reduce duplicate calculations
6. Reduce type conversions

# Using float64 vs Decimal struct

If I use float64, I cannot compare values, because there will be infinitesimal differences between values like 1e-14 and -1e-14. I need to truncate first and then place it into the map/set data structure.

However, there are also problems in using Decimals. If I use Decimal, I will also truncate first and then place it into the set structure. However, Decimals will not equal to each other because the struct definition contains a pointer. Two Decimal structs can have the same value, but will not be equal to each other if their pointers point to different memory addresses

A possible workaround is to take .StringFixed(places int32) of the truncated Decimal struct and place it into a hashmap as the key, with the Decimal struct as the value.
I can then compare two maps for equivalence by only comparing their keys, by checking that all the keys in map1 are found in map2, and that their lengths are the same.
