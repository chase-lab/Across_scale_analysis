#' pivot_longer and filter model results in 3D arrays
#' @param x A 3D array object where dimension 1 is dataset_id, dimension 2 is
#' for model results and dimension 3 is for land uses.
#' @param filter a vector of Dataset_id_Land_use values that we need.
#' @returns a tibble with columns Dataset_id, Land_use, Estimate, EstError (without
#' point), lower and upper.
pivot_filter <- function(x, filter_vector) {
  dimnames(x)[[2]] <- c("Estimate","EstError","lower","upper")
  as_tibble(x, rownames = NA) |> 
    mutate(Dataset_id = rownames(x)) |> 
    pivot_longer(cols = !Dataset_id, names_sep = "\\.", names_to = c("variable", "Land_use")) |> 
    pivot_wider(names_from = variable) |> 
    mutate(Land_use = gsub("Land_use", replacement = "", x = Land_use, fixed = TRUE)) |> 
    # reduce to study-Land_use combinations in the data
    unite(col = "filter", c(Dataset_id, Land_use), remove = FALSE) |> 
    filter(filter %in% filter_vector)
}