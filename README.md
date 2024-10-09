# MGARD_unstructured

## Build
```bash ./build.sh```

## Config file

Description goes here

## Merge data
```./merge config.json input.bp output.bp```

## Reorder data
```./reorder config.json input.bp output.bp```

## Compress data
```./compress_adios config.json input.bp output.bp```

```./compress_mgard config.json input.bp output.bp```

## Full workflow
```./merge config.json input.bp merged.bp```

```./reorder config.json merged.bp reindexed.bp```

```./compress_mgard config.json reindexed.bp output.bp```
