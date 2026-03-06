# Deduplication of DG nodes

## Config file

See the attached config.json example.
Need to list:
* names of connectivity, coordinates and variables
* max_block_id=-1 means all blocks
* max_step=-1 means all time steps
* n_nodes_in_element=8 is hardcoded in hpMusic code
* convert_to_cgns_standard enables additional fields with types of elements as specified in CGNS standard
* verbose enables printing progress messages

## Build
```bash ./build.sh```

## Deduplicate DG data to make it continuous
```./bin/dedup config.json input.bp output.bp```
