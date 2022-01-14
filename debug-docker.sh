#!/usr/bin/env bash
docker run -v $(pwd)/examples:/data -v $(pwd)/run.py:/opt/run.py -i -t quay.io/antonkulaga/haddock-antibody:latest /bin/bash
