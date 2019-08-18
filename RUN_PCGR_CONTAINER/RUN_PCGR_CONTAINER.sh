#!/bin/bash

if [ $# -lt 1 ]
then
        echo "ERROR: path to output directory argument missing"
        exit
fi

# pull latest version of pcgr container; will skip if up-to-date
docker pull tsharding/pcgr_v1.0

#supplies arg is path-to-your-output-directory
docker-compose run --rm  -v $1:/OUTPUT main bash

#stop and remove selenium server container
docker-compose down -v
