#!/usr/bin/env bash
docker run -v $(pwd)/examples:/data -v $(pwd)/start.py:/opt/start.py -i -t quay.io/antonkulaga/haddock-antibody:latest /bin/bash
