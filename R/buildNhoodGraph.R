### Neighbourhood abstracted graph ###

#' Build an abstracted graph of neighbourhoods for visualization
#'
#'
#' @param x A \code{\linkS4class{Milo}} object with a non-empty \code{nhoods}
#' slot.
#' @param overlap A numeric scalar that thresholds graph edges based on  the number
#' of overlapping cells between neighbourhoods.
#'
#' @details
#' This constructs a weighted graph where nodes represent neighbourhoods and edges represent the number of overlapping
#' cells between two neighbourhoods.
#'
#' @return A \code{\linkS4class{Milo}} object containg an \code{igraph} graph in the \code{nhoodGraph} slot.
#'
#' @author
#' Emma Dann
#'
#' @examples
#'
#' NULL
#'
#' @importFrom igraph graph.adjacency set_vertex_attr vertex.attributes
#' @export
#' @rdname buildNhoodGraph
buildNhoodGraph <- function(x, overlap=1, graph_name='miloR', assay='RNA'){

    if(!(is(x, "Milo") || is(x, 'Seurat'))){
      stop("Not a valid object. Class should be Milo or Seurat")
    }

    # are neighbourhoods calculated?
    if (is(x, "Milo")) {
      if(ncol(nhoods(x)) == 1 & nrow(nhoods(x)) == 1){
        stop("No neighbourhoods found - run makeNhoods first")
      }
      ## Build adjacency matrix for nhoods
      nh_intersect_mat <- .build_nhood_adjacency(nhoods(x), overlap=overlap)
      nhood_sizes <- colSums(nhoods(x))
      # add to slot if empty
      nhoodAdjacency(x) <- nh_intersect_mat
    }
    
    if (is(x, "Seurat")) {
      if(ncol(x@neighbors[[graph_name]]$nhoods) == 1 & nrow(x@neighbors[[graph_name]]$nhoods) == 1){
        stop("No neighbourhoods found - run makeNhoods first")
      }
      ## Build adjacency matrix for nhoods
      nh_intersect_mat <- .build_nhood_adjacency(x@neighbors[[graph_name]]$nhoods, overlap=overlap)
      nhood_sizes <- colSums(x@neighbors[[graph_name]]$nhoods)
      # add to slot if empty
      x@neighbors[[graph_name]]$nhoodAdjacency <- nh_intersect_mat
    }

    ## Make igraph object
    ig <- graph.adjacency(nh_intersect_mat, mode="undirected", weighted=TRUE)
    ig <- set_vertex_attr(ig, name = 'size', value = nhood_sizes[vertex.attributes(ig)$name])
    if (is(x, 'Milo')) {
      ## Add to nhoodGraph slot in milo object
      nhoodGraph(x) <- ig
    }
    if (is(x, 'Seurat')) {
      x@neighbors[[graph_name]]$nhoodGraph <- ig
      command_meta_info <- new("SeuratCommand",
                               name="MiloR::buildNhoodGraph",
                               assay.used=assay,
                               call.string="MiloR::buildNhoodGraph(object)",
                               params=list(overlap=overlap,
                                           graph_name=graph_name),
                               time.stamp=as.POSIXct(Sys.time()))
      x@commands[['MiloR::buildNhoodGraph']] <- command_meta_info
    }
    return(x)
}


