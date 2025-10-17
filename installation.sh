First install dorado, minimap2, and modkit softwares.
^ Start by setting up apptainer environment within your scratch directory

##Add this into script to create the apptainer cache in your scratch directory.

mkdir -p /scratch/$USER/apptainer_images
apptainer pull /scratch/$USER/apptainer_images/ubuntu20.sif docker://ubuntu:20.04 
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp 
mkdir -p $APPTAINER_CACHEDIR $APPTAINER_TMPDIR

#dorado installation (just copy over for 1.0.1), any newer versions will need to be downloaded from 
https://github.com/nanoporetech/dorado
### Heres the path to my installation of the new version as of ~ August 2025 (unlocked permissions) 
"/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-1.0.1-linux-x64/"

##minimap 
curl -L https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.30_x64-linux/minimap2

#modkit 
#heres my path with unlocked permissions 
"/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/modkit_v0.5/dist_modkit_v0.5.0_5120ef7/modkit"
#otherwise install from here
https://github.com/nanoporetech/modkit/releases
