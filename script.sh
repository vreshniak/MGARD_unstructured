# ./merge   config.json /lustre/orion/cfd164/proj-shared/pugmire/cgns2adios/output_vki.bp  /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged.bp
# ./reorder config.json /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged.bp    /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_reordered.bp

./compress_adios config.json /lustre/orion/cfd164/proj-shared/pugmire/cgns2adios/output_vki.bp    /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_adios.bp
# ./compress_adios config.json /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged.bp    /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_adios_merged.bp
# ./compress_adios config.json /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_reordered.bp /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_adios_reordered.bp

# ./compress_mgard config.json /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged.bp    /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_mgard_merged.bp    > compressed_mgard_merged_out
# ./compress_mgard config.json /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_reordered.bp /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_mgard_reordered.bp > compressed_mgard_reordered_out

# ./extract_data config.json /lustre/orion/cfd164/proj-shared/pugmire/cgns2adios/output_vki.bp   /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_data.bp
# ./extract_data config.json /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged.bp     /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged_data.bp
# ./extract_data config.json /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_reordered.bp  /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_reordered_data.bp


du -h /lustre/orion/cfd164/proj-shared/pugmire/cgns2adios/output_vki.bp
du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_data.bp
# du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged.bp
du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged_data.bp
# du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_reordered.bp
du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_reordered_data.bp
# du -b /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_merged.bp
# du -b /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_reordered.bp

# du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_adios_merged.bp
# du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_adios_reordered.bp

du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_mgard_merged.bp
du -h /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_compressed_mgard_reordered.bp

# bpls /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_reordered_data.bp -D FlowSolution/RHE | head -n 10
# bpls /lustre/orion/cfd164/proj-shared/reshniakv/output_vki_merged_data.bp -D GridCoordinates/CoordinateX | head -n 10