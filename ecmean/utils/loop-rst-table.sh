#dumb script to generate all the tables from the clim/reference

#!/bin/bash

for name in EC24 EC26-HIST EC26-CMIP; do
    echo "Generating table for $kind $name"
    python generate-rst-table.py ../climatology/$name/pi_climatology_$name.yml --output ../../docs/sphinx/source/tables/climatology_$name.rst
done

for name in EC23 EC26-HIST EC26-CMIP EC26-PDAY; do
    echo "Generating table for $kind $name"
    python generate-rst-table.py ../reference/gm_reference_$name.yml --output ../../docs/sphinx/source/tables/reference_$name.rst
done
