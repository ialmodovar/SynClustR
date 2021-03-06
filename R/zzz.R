### Lastest load into a package.

### Export Namespace does not use .First.lib() and .Last.lib(), but use
### .onLoad() and .onUnload().
# .First.lib <- function(lib, pkg){
# } # End of .First.lib().

# .Last.lib <- function(libpath){
# } # End of .Last.lib().

.onLoad <- function(lib, pkg) {
    library.dynam("SynClustR", pkg, lib)
} ## End of .onLoad()

.onUnload <- function(libpath){
  library.dynam.unload("SynClustR", libpath)
} ## End of .onUnload().
