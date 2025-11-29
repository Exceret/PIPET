PIPET.optimized <- function(
    sc_data,
    markers,
    group = NULL,
    rm_NA = TRUE,
    freq_counts = NULL,
    normalize = TRUE,
    scale = TRUE,
    nPerm = 1000,
    distance = "cosine",
    ...
) {
    # Handle group-wise analysis
    if (!is.null(group)) {
        return(PIPET_GroupAnalysis(
            Seurat_data = sc_data,
            markers = markers,
            group = group,
            rm_NA = rm_NA,
            freq_counts = freq_counts,
            normalize = normalize,
            scale = scale,
            nPerm = nPerm,
            distance = distance,
            nCores = nCores
        ))
    }
    PIPET_SingleAnalysis(
        sc_data = sc_data,
        markers = markers,
        rm_NA = rm_NA,
        freq_counts = freq_counts,
        normalize = normalize,
        scale = scale,
        nPerm = nPerm,
        distance = distance,
        nCores = nCores
    )
}
