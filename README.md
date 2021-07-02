# oharvester

**Goal**: A repository with a library of functions, scripts, and tutorials for harvesting oceanographic data.

* Each tool works basically like this:

```
toolX(time range, spatial range, data source)
  {
    gets [data] from specified source
    returns [data]
    (or possibly saves [data] to a file)
  }
```
    
Sometimes the needs change a little bit, so the arguments to the functions can change, too.  For example, if you need to retrieve point data the too might work as shown below. In this case, `points` might be a data frame with `lon`, `lat`, and `date` columns.

```
toolX(points, data source)
  {
    gets [data] from specified source
    returns [data]
    (or possibly saves [data] to a file)
  }
```


## Other resources

Cath Mitchell's [oc-satellite-misc](https://github.com/cathmmitchell/oc-satellite-misc) repos.