conda install --yes --force numpy scipy networkx decorator cython nose numba llvmlite funcsigs pandas enum34 dateutil pytz rtree jsonschema pyproj six && \
conda install --yes gdal

# for networkplanner integration
conda install --yes -c https://conda.binstar.org/sel networkplanner-metrics
