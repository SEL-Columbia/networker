## 0.2.7
- Fixes issue #80, where lon/lat to x/y/z conversions were incorrect
- Added SpatialRef checks (issue #81)

## 0.2.6
- updated library versions and handle osx install

## 0.2.5
- support networkx 1.10 connected_components function change that now returns sets

## 0.2.4

- Now handles MultiLineString geometry from shapefile
- Adds some more error handling/reporting when network is not as expected
- pep8 compliance changes (better)
- Issues addressed:
    #71, 72

## 0.2.3

- Handle network inputs that have self-referencing (i.e. zero-length) edges
  (via removal and warning messages)
- Add clean_network script to remove these edges from any shapefile
- pep8 compliance changes (still not perfect)
- Issues addressed:
    #67
