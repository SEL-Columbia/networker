conda install --yes --force numpy scipy networkx decorator cython nose numba llvmlite funcsigs pandas enum34 dateutil pytz && \
conda install --yes -c https://conda.binstar.org/osgeo gdal && \
conda install --yes -c https://conda.binstar.org/dougal libspatialindex && \
conda install --yes -c https://conda.binstar.org/sel rtree

# for networkplanner integration
# TODO: Rethink this
conda install --yes --force sqlalchemy==0.7.8 shapely && \
conda install --yes -c https://conda.binstar.org/rsignell geojson
