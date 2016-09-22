## 0.4.1
- Address coordinates as ndarray issue (#101)
- Allow label offset in draw_geograph (#102)

## 0.4.0
- Added GeoJSON support (see spec in [docs/geograph_geojson.md](docs/geograph_geojson.md))
- Added `compose` tool for combining and converting geographs (see demo in [docs/geograph_compose_demo.md](docs/geograph_compose_demo.md))
- Included changes from @invisibleroads fork in support of [infrastructure-planning](https://github.com/SEL-Columbia/infrastructure-planning)

## 0.3.0
- Added spherical_accuracy to network input configuration (defaults to false). 
  This addresses [#82](https://github.com/SEL-Columbia/networker/issues/82)
- Added more tools in scripts plus [documentation](https://github.com/SEL-Columbia/networker/tree/9e10ac319ef2d8531d951ad9eb9f3cbb524758da#other-tools) on them in README

## 0.2.7
- Fixes issue [#80](https://github.com/SEL-Columbia/networker/issues/80), where lon/lat to x/y/z conversions were incorrect
- Added SpatialRef checks, issue [#81](https://github.com/SEL-Columbia/networker/issues/81)

## 0.2.6
- updated library versions (issue [#76](https://github.com/SEL-Columbia/networker/issues/76)) and handle osx install (issue [#78](https://github.com/SEL-Columbia/networker/issues/78))

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
