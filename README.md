Docker container for haddock-antibody protocol
-----------------------------------------------
Docker container for [haddock-antibody protocol](https://github.com/haddocking/HADDOCK-antibody-antigen) 
The container also contains haddock-tools, freesasa and additional scripts required to prepare intput for the local version of haddock

Usage
-----

For example, let's prepare data for docking 4G6K antibody (therapeutic antibody binding fragment of gevokizumab) to 4I1B.pdb as antigen.
```bash
docker run -v $(pwd)/examples:/data -i -t quay.io/antonkulaga/haddock-antibody:latest /bin/bash
```
and then inside of the container:
```bash
pdb_fetch  4G6K > 4G6K.pdb
pdb_fetch 4I1B > 4I1B.pdb
run.py run --antibody 4G6K.pdb --antigen 4I1B.pdb --output output
```
Building container
------------------

Just run:
```bash
./build-docker.sh
```