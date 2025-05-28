# MGARD_unstructured

## Build
```bash ./build.sh```

## Config file

See the attached config.json example.
Need to list:
* names of connectivity, coordinates and variables
* path to meshfile to be used in plugin
* name of the ordering vector in the meshfile
* s, rel_tol are compression parameters
* max_block_id=-1 means all blocks
* max_step=-1 means all time steps
* n_nodes_in_element=8 is hardcoded


## Extract data variables only (remove mesh)
```./extract_data config.json input.bp output.bp```

## Extract mesh data only and add node ordering (remove variables) - create meshfile
```./make_meshfile_with_ordering config.json input.bp output.bp```

## Merge DG data to make it continuous
```./merge config.json input.bp output.bp```

## Compress data
```./compress_with_adios_plugin config.json input.bp output.bp```

## Decompress data
```./decompress_with_adios_plugin config.json input.bp output.bp```

## Calculate errors
```./evaluate_compression_error config.json input.bp output.bp```

## Full workflow
```./make_meshfile_with_ordering config.json input.bp meshfile.bp```

```./merge config.json input.bp merged.bp```

```./compress_with_adios_plugin config.json merged.bp compressed.bp```

```./decompress_with_adios_plugin config.json compressed.bp decompressed.bp```

```./evaluate_compression_error config.json merged.bp decompressed.bp```
