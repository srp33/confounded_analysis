function (data, isProximityMatrix = FALSE, XdNum = NULL, XcNum = NULL, 
    YdNum = NULL, YcNum = NULL, row.name = NULL, Xd.name = NULL, 
    Xc.name = NULL, row.prox = NULL, col.prox = NULL, show.row.prox = TRUE, 
    show.col.prox = TRUE, row.order = NULL, col.order = NULL, 
    row.flip = NULL, col.flip = NULL, row.externalOrder = NULL, 
    col.externalOrder = NULL, original.color = NULL, row.color = NULL, 
    col.color = NULL, Xd.color = NULL, Xc.color = NULL, Yd.color = NULL, 
    Yc.color = NULL, row.label.size = NULL, col.label.size = NULL, 
    Xd.label.size = NULL, Xc.label.size = NULL, Yd.label.size = NULL, 
    Yc.label.size = NULL, colorbar.margin = 1.5, border = FALSE, 
    border.width = 1, isContainMissingValue = 0, MissingValue.color = "gray", 
    exp.row_order = FALSE, exp.column_order = FALSE, exp.row_names = FALSE, 
    exp.column_names = FALSE, exp.Xc = FALSE, exp.Yc = FALSE, 
    exp.Xd = FALSE, exp.Yd = FALSE, exp.Xd_codebook = FALSE, 
    exp.Yd_codebook = FALSE, exp.originalmatrix = FALSE, exp.row_prox = FALSE, 
    exp.col_prox = FALSE, PNGfilename = NULL, PNGwidth = 1800, 
    PNGheight = 1200, PNGres = 150, show.plot = FALSE) 
{
    is_binary_data <- function(mat, allow_na = TRUE) {
        if (isFALSE(isProximityMatrix)) {
            if (allow_na) {
                all(mat %in% c(0, 1, NA), na.rm = TRUE)
            }
            else {
                all(mat %in% c(0, 1))
            }
        }
    }
    colorbar_list <- list()
    GAP_Blue_White_Red_function <- colorRamp2(breaks = seq(-1, 
        1, length.out = length(GAP_Blue_White_Red)), colors = GAP_Blue_White_Red)
    colorbar_fun <- function(data, colorspectrum) {
        if (colorspectrum %in% c("GAP_Rainbow", "GAP_Blue_White_Red", 
            "GAP_d")) {
            num_colorspectrum <- length(get(colorspectrum, envir = .GlobalEnv))
            colorbar_colors <- get(colorspectrum, envir = .GlobalEnv)
        }
        else if (colorspectrum %in% rownames(brewer.pal.info)) {
            num_colorspectrum <- brewer.pal.info[colorspectrum, 
                "maxcolors"]
            colorbar_colors <- brewer.pal(num_colorspectrum, 
                colorspectrum)
        }
        else if (colorspectrum == "grayscale_palette") {
            num_colorspectrum <- 2
            colorbar_colors <- get(colorspectrum, envir = .GlobalEnv)
        }
        if (colorspectrum %in% c("BrBG", "PiYG", "PRGn", "PuOr", 
            "RdBu", "RdGy", "GAP_Blue_White_Red")) {
            colorbar_min <- -1
            colorbar_max <- 1
            colorbar_mid <- 0
            colorbar_breaks <- seq(-1, 1, length.out = num_colorspectrum)
            colorbar_function <- colorRamp2(breaks = colorbar_breaks, 
                colors = colorbar_colors)
        }
        else {
            colorbar_min <- round(min(data, na.rm = TRUE), digits = 2)
            colorbar_max <- round(max(data, na.rm = TRUE), digits = 2)
            colorbar_breaks <- round(seq(colorbar_min, colorbar_max, 
                length.out = num_colorspectrum), digits = 2)
            if (length(colorbar_breaks)%%2 == 1) {
                colorbar_mid <- colorbar_breaks[ceiling(length(colorbar_breaks)/2)]
            }
            else {
                middle_indices <- c(length(colorbar_breaks)/2, 
                  length(colorbar_breaks)/2 + 1)
                colorbar_mid <- mean(colorbar_breaks[middle_indices])
            }
            colorbar_function <- colorRamp2(breaks = colorbar_breaks, 
                colors = colorbar_colors)
        }
        return(list(colorbar_min = colorbar_min, colorbar_max = colorbar_max, 
            colorbar_mid = colorbar_mid, num_colorspectrum = num_colorspectrum, 
            colorbar_breaks = colorbar_breaks, colorbar_function = colorbar_function))
    }
    if (isTRUE(isProximityMatrix)) {
        XcNum <- NULL
        XdNum <- NULL
        YcNum <- NULL
        YdNum <- NULL
        row.prox <- NULL
        col.prox <- NULL
        col.order <- NULL
        col.flip <- NULL
        col.externalOrder <- NULL
    }
    rows_to_remove <- if (!is.null(c(XcNum, XdNum))) 
        sort(c(XcNum, XdNum))
    else NULL
    rowname_cols <- NULL
    if (!is.null(row.name) && length(row.name) != nrow(data)) {
        if (is.numeric(row.name)) {
            rowname_cols <- row.name
        }
    }
    cols_to_remove <- sort(unique(c(YcNum, YdNum, rowname_cols)))
    cols_to_remove <- if (length(cols_to_remove) == 0) 
        NULL
    else cols_to_remove
    if (!is.null(XcNum)) {
        if (is.null(cols_to_remove)) {
            Xc <- as.data.frame(data[XcNum, , drop = FALSE])
        }
        else {
            Xc <- as.data.frame(data[XcNum, -cols_to_remove, 
                drop = FALSE])
        }
    }
    else {
        Xc <- NULL
    }
    if (!is.null(XdNum)) {
        if (is.null(cols_to_remove)) {
            Xd <- as.data.frame(data[XdNum, , drop = FALSE])
        }
        else {
            Xd <- as.data.frame(data[XdNum, -cols_to_remove, 
                drop = FALSE])
        }
    }
    else {
        Xd <- NULL
    }
    if (!is.null(YcNum)) {
        if (is.null(rows_to_remove)) {
            Yc <- as.data.frame(data[, YcNum, drop = FALSE])
        }
        else {
            Yc <- as.data.frame(data[-rows_to_remove, YcNum, 
                drop = FALSE])
        }
    }
    else {
        Yc <- NULL
    }
    if (!is.null(YdNum)) {
        if (is.null(rows_to_remove)) {
            Yd <- as.data.frame(data[, YdNum, drop = FALSE])
        }
        else {
            Yd <- as.data.frame(data[-rows_to_remove, YdNum, 
                drop = FALSE])
        }
    }
    else {
        Yd <- NULL
    }
    if (isFALSE(isProximityMatrix)) {
        if (is.null(cols_to_remove) && is.null(rows_to_remove)) {
            data_matrix <- as.matrix(sapply(data, as.numeric))
            if (is.null(row.name)) {
                row_names <- 1:nrow(data)
            }
            else if (length(row.name) == nrow(data)) {
                row_names <- row.name
            }
            else if (length(row.name) == 1) {
                row_names <- data[[row.name]]
            }
            else {
                row_names <- do.call(paste0, as.list(data[, row.name]))
            }
        }
        else if (is.null(cols_to_remove) && !is.null(rows_to_remove)) {
            data_matrix <- as.matrix(sapply(data[-rows_to_remove, 
                ], as.numeric))
            if (is.null(row.name)) {
                row_names <- 1:nrow(data[-rows_to_remove, ])
            }
            else if (length(row.name) == nrow(data)) {
                row_names <- row.name
            }
            else if (length(row.name) == 1) {
                row_names <- data[[row.name]][-rows_to_remove]
            }
            else {
                row_names <- do.call(paste0, as.list(data[-rows_to_remove, 
                  row.name]))
            }
        }
        else if (!is.null(cols_to_remove) && is.null(rows_to_remove)) {
            data_matrix <- as.matrix(sapply(data[, -cols_to_remove], 
                as.numeric))
            if (is.null(row.name)) {
                row_names <- 1:nrow(data)
            }
            else if (length(row.name) == nrow(data)) {
                row_names <- row.name
            }
            else if (length(row.name) == 1) {
                row_names <- data[[row.name]]
            }
            else {
                row_names <- do.call(paste0, as.list(data[, row.name]))
            }
        }
        else {
            data_matrix <- as.matrix(sapply(data[-rows_to_remove, 
                -cols_to_remove], as.numeric))
            if (is.null(row.name)) {
                row_names <- 1:nrow(data[-rows_to_remove, ])
            }
            else if (length(row.name) == nrow(data)) {
                row_names <- row.name
            }
            else if (length(row.name) == 1) {
                row_names <- data[[row.name]][-rows_to_remove]
            }
            else {
                row_names <- do.call(paste0, as.list(data[-rows_to_remove, 
                  row.name]))
            }
        }
    }
    else {
        data_matrix <- data
        if (is.null(row.name)) {
            row_names <- 1:nrow(data)
        }
        else if (length(row.name) == nrow(data)) {
            row_names <- row.name
        }
        else if (length(row.name) == 1) {
            row_names <- data[[row.name]]
        }
        else {
            row_names <- do.call(paste0, as.list(data[, row.name]))
        }
    }
    if (row.order == "original" && length(row.name) != nrow(data)) {
        row_names <- rev(row_names)
    }
    column_names <- colnames(data_matrix)
    column_names <- tolower(column_names)
    if (is.null(Xd.name)) {
        Xd.name <- 1:length(XdNum)
    }
    if (is.null(Xc.name)) {
        Xc.name <- 1:length(XcNum)
    }
    if (isContainMissingValue == 1 && !(row.order %in% c("original", 
        "random", "reverse"))) {
        stop("Error: Missing values detected, but row.order must be one of 'original', 'random', or 'reverse'.")
    }
    if (isContainMissingValue == 1 && !(col.order %in% c("original", 
        "random", "reverse"))) {
        stop("Error: Missing values detected, but col.order must be one of 'original', 'random', or 'reverse'.")
    }
    if (is.null(row.prox)) {
        show.row.prox <- FALSE
    }
    if (is.null(col.prox)) {
        show.col.prox <- FALSE
    }
    if (isFALSE(isProximityMatrix)) {
        if ((is.null(row.prox) && !(row.order %in% c("original", 
            "random", "reverse"))) || (is.null(col.prox) && !(col.order %in% 
            c("original", "random", "reverse")))) {
            stop("Error: The proximity matrix must be computed prior to applying the specified ordering method.")
        }
    }
    if (isFALSE(isProximityMatrix)) {
        if (!isTRUE(is_binary_data(data_matrix, TRUE)) && ((!is.null(row.prox) && 
            row.prox %in% c("hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule")) || (!is.null(col.prox) && 
            col.prox %in% c("hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule")))) {
            problematic_prox <- if (row.prox %in% c("hamman", 
                "jaccard", "phi", "rao", "rogers", "simple", 
                "sneath", "yule")) {
                row.prox
            }
            else {
                col.prox
            }
            stop(sprintf("Error: The proximity method '%s' requires binary data. Please ensure your data matrix contains only 0 and 1 values.", 
                problematic_prox))
        }
    }
    if (isFALSE(isProximityMatrix)) {
        if (isTRUE(is_binary_data(data_matrix, TRUE)) && ((!is.null(row.prox) && 
            !row.prox %in% c("hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule")) || (!is.null(col.prox) && 
            !col.prox %in% c("hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule")))) {
            problematic_prox <- if (!row.prox %in% c("hamman", 
                "jaccard", "phi", "rao", "rogers", "simple", 
                "sneath", "yule")) {
                row.prox
            }
            else {
                col.prox
            }
            stop(sprintf("Error: The proximity method %s is not compatible with binary data. Please choose a binary-compatible method.", 
                problematic_prox))
        }
    }
    proxType_map <- list(euclidean = 0, pearson = 1, kendall = 2, 
        spearman = 3, atancorr = 4, `city-block` = 5, abs_pearson = 6, 
        uncenteredcorr = 7, abs_uncenteredcorr = 8, hamman = 20, 
        jaccard = 21, phi = 22, rao = 23, rogers = 24, simple = 25, 
        sneath = 26, yule = 27)
    if (!is.null(row.prox)) {
        if (row.prox %in% names(proxType_map)) {
            proxType <- proxType_map[[row.prox]]
            distance_matrix_for_row <- computeProximity(data_matrix, 
                proxType, side = 0, isContainMissingValue)
            if (row.prox %in% c("euclidean", "city-block")) {
                distance_matrix_for_row_hc <- distance_matrix_for_row
            }
            else {
                distance_matrix_for_row_hc <- as.dist(1 - distance_matrix_for_row)
            }
        }
        else if (row.prox %in% c("maximum", "canberra")) {
            distance_matrix_for_row <- dist(data_matrix, method = row.prox)
            distance_matrix_for_row_hc <- distance_matrix_for_row
        }
    }
    if (is.null(row.prox)) {
        if (isFALSE(isProximityMatrix)) {
            distance_matrix_for_row <- NULL
            distance_matrix_for_row_hc <- NULL
        }
    }
    else {
        distance_matrix_for_row <- as.matrix(distance_matrix_for_row)
        distance_matrix_for_row_hc <- as.matrix(distance_matrix_for_row_hc)
        distance_matrix_for_row_dist <- as.dist(distance_matrix_for_row_hc)
    }
    if (!is.null(col.prox)) {
        if (col.prox %in% names(proxType_map)) {
            proxType <- proxType_map[[col.prox]]
            distance_matrix_for_column <- computeProximity(data_matrix, 
                proxType, side = 1, isContainMissingValue)
            if (col.prox %in% c("euclidean", "city-block")) {
                distance_matrix_for_column_hc <- distance_matrix_for_column
            }
            else {
                distance_matrix_for_column_hc <- as.dist(1 - 
                  distance_matrix_for_column)
            }
        }
        else if (col.prox %in% c("maximum", "canberra")) {
            distance_matrix_for_column <- dist(t(data_matrix), 
                method = col.prox)
            distance_matrix_for_column_hc <- distance_matrix_for_column
        }
    }
    if (is.null(col.prox)) {
        if (isFALSE(isProximityMatrix)) {
            distance_matrix_for_column <- NULL
            distance_matrix_for_column_hc <- NULL
        }
    }
    else {
        distance_matrix_for_column <- as.matrix(distance_matrix_for_column)
        distance_matrix_for_column_hc <- as.matrix(distance_matrix_for_column_hc)
        distance_matrix_for_column_dist <- as.dist(distance_matrix_for_column_hc)
    }
    if (isTRUE(row.order == "r2e") || isTRUE(col.order == "r2e") || 
        isTRUE(row.flip == "r2e") || isTRUE(col.flip == "r2e")) {
        if (!is.null(row.prox)) {
            r2e_order_row <- ellipse_sort(distance_matrix_for_row)
        }
        else if (isTRUE(isProximityMatrix)) {
            r2e_order_row <- ellipse_sort(data_matrix)
        }
        if (!is.null(col.prox)) {
            r2e_order_col <- ellipse_sort(distance_matrix_for_column)
        }
    }
    seriation_dist <- c("ARSA", "BBURCG", "BBWRCG", "Enumerate", 
        "GSA", "GW", "GW_average", "GW_complete", "GW_single", 
        "GW_ward", "isomap", "isoMDS", "MDS", "MDS_angle", "metaMDS", 
        "monoMDS", "OLO", "OLO_average", "OLO_complete", "OLO_single", 
        "OLO_ward", "QAP_2SUM", "QAP_BAR", "QAP_Inertia", "QAP_LS", 
        "Sammon_mapping", "SGD", "Spectral", "Spectral_norm", 
        "SPIN_NH", "SPIN_STS", "TSP", "VAT")
    seriation_matrix <- c("AOE", "BEA", "BEA_TSP", "BK_unconstrained", 
        "CA", "Heatmap", "LLE", "Mean", "PCA", "PCA_angle")
    if (row.order == "random") {
        random_order_row <- sample(1:nrow(data_matrix))
        if (isFALSE(isProximityMatrix)) {
            data_matrix <- data_matrix[random_order_row, ]
            distance_matrix_for_row <- distance_matrix_for_row[random_order_row, 
                random_order_row]
        }
        else {
            data_matrix <- data_matrix[random_order_row, random_order_row]
        }
        row_names <- row_names[random_order_row]
        row_names <- rev(row_names)
        order_row <- random_order_row
    }
    else if (row.order == "reverse") {
        reverse_order_row <- nrow(data_matrix):1
        if (isFALSE(isProximityMatrix)) {
            data_matrix <- data_matrix[reverse_order_row, ]
            distance_matrix_for_row <- distance_matrix_for_row[reverse_order_row, 
                reverse_order_row]
        }
        else {
            data_matrix <- data_matrix[reverse_order_row, reverse_order_row]
        }
        row_names <- row_names[reverse_order_row]
        row_names <- rev(row_names)
        order_row <- reverse_order_row
    }
    else if (row.order == "r2e") {
        if (isFALSE(isProximityMatrix)) {
            data_matrix <- data_matrix[r2e_order_row, ]
            distance_matrix_for_row <- distance_matrix_for_row[r2e_order_row, 
                r2e_order_row]
        }
        else {
            data_matrix <- data_matrix[r2e_order_row, r2e_order_row]
        }
        row_names <- row_names[r2e_order_row]
        row_names <- rev(row_names)
        order_row <- r2e_order_row
    }
    else if (row.order %in% seriation_dist) {
        if (isFALSE(isProximityMatrix)) {
            seriation_order_row <- get_order(seriate(distance_matrix_for_row_dist, 
                method = row.order))
            data_matrix <- data_matrix[seriation_order_row, ]
            distance_matrix_for_row <- distance_matrix_for_row[seriation_order_row, 
                seriation_order_row]
        }
        else {
            seriation_order_row <- get_order(seriate(as.dist(data_matrix), 
                method = row.order))
            data_matrix <- data_matrix[seriation_order_row, seriation_order_row]
        }
        row_names <- row_names[seriation_order_row]
        row_names <- rev(row_names)
        order_row <- seriation_order_row
    }
    else if (row.order %in% seriation_matrix) {
        if (isFALSE(isProximityMatrix)) {
            seriation_order_row <- get_order(seriate(distance_matrix_for_row_hc, 
                method = row.order))
            data_matrix <- data_matrix[seriation_order_row, ]
            distance_matrix_for_row <- distance_matrix_for_row[seriation_order_row, 
                seriation_order_row]
        }
        else {
            seriation_order_row <- get_order(seriate(data_matrix, 
                method = row.order))
            data_matrix <- data_matrix[seriation_order_row, seriation_order_row]
        }
        row_names <- row_names[seriation_order_row]
        row_names <- rev(row_names)
        order_row <- seriation_order_row
    }
    if (isFALSE(isProximityMatrix)) {
        if (col.order == "random") {
            random_order_col <- sample(1:ncol(data_matrix))
            data_matrix <- data_matrix[, random_order_col]
            distance_matrix_for_column <- distance_matrix_for_column[random_order_col, 
                random_order_col]
            column_names <- column_names[random_order_col]
            order_col <- random_order_col
        }
        else if (col.order == "reverse") {
            reverse_order_col <- ncol(data_matrix):1
            data_matrix <- data_matrix[, reverse_order_col]
            distance_matrix_for_column <- distance_matrix_for_column[reverse_order_col, 
                reverse_order_col]
            column_names <- column_names[reverse_order_col]
            order_col <- reverse_order_col
        }
        else if (col.order == "r2e") {
            data_matrix <- data_matrix[, r2e_order_col]
            distance_matrix_for_column <- distance_matrix_for_column[r2e_order_col, 
                r2e_order_col]
            column_names <- column_names[r2e_order_col]
            order_col <- r2e_order_col
        }
        else if (col.order %in% seriation_dist) {
            seriation_order_col <- get_order(seriate(distance_matrix_for_column_dist, 
                method = col.order))
            data_matrix <- data_matrix[, seriation_order_col]
            distance_matrix_for_column <- distance_matrix_for_column[seriation_order_col, 
                seriation_order_col]
            column_names <- column_names[seriation_order_col]
            order_col <- seriation_order_col
        }
        else if (col.order %in% seriation_matrix) {
            seriation_order_col <- get_order(seriate(distance_matrix_for_column_hc, 
                method = col.order))
            data_matrix <- data_matrix[, seriation_order_col]
            distance_matrix_for_column <- distance_matrix_for_column[seriation_order_col, 
                seriation_order_col]
            column_names <- column_names[seriation_order_col]
            order_col <- seriation_order_col
        }
    }
    orderType_map <- list(single = 0, complete = 1, average = 2)
    flipType_map <- list(r2e = 1, uncle = 2, grandpa = 3)
    if (isTRUE(row.order == "r2e") && !is.null(row.flip) && !isTRUE(row.flip == 
        "r2e")) {
        stop("Error: 'r2e' (Rank Two Ellipse) is not a Hierarchical Clustering Tree (HCT) method.")
    }
    if (isFALSE(isProximityMatrix)) {
        if (isTRUE(col.order == "r2e") && !is.null(col.flip) && 
            !isTRUE(col.flip == "r2e")) {
            stop("Error: 'r2e' (Rank Two Ellipse) is not a Hierarchical Clustering Tree (HCT) method.")
        }
    }
    if (!isTRUE(row.order %in% c("r2e", "single", "complete", 
        "average")) && !is.null(row.flip)) {
        stop(sprintf("Error: '%s' is not a Hierarchical Clustering Tree (HCT) method.", 
            row.order))
    }
    if (isFALSE(isProximityMatrix)) {
        if (!isTRUE(col.order %in% c("r2e", "single", "complete", 
            "average")) && !is.null(col.flip)) {
            stop(sprintf("Error: '%s' is not a Hierarchical Clustering Tree (HCT) method.", 
                col.order))
        }
    }
    if (!is.null(row.externalOrder) && !is.null(row.flip)) {
        stop("Error: row.externalOrder and row.flip cannot both have values at the same time.")
    }
    if (!is.null(col.externalOrder) && !is.null(col.flip)) {
        stop("Error: col.externalOrder and col.flip cannot both have values at the same time.")
    }
    if (is.null(row.externalOrder) && isTRUE(row.flip == "r2e")) {
        default_rowexternalOrder <- r2e_order_row - 1
    }
    else if (!(is.null(row.externalOrder)) && is.null(row.flip)) {
        default_rowexternalOrder <- row.externalOrder
    }
    else {
        default_rowexternalOrder <- seq(0, nrow(data_matrix) - 
            1)
    }
    if (is.null(col.externalOrder) && isTRUE(col.flip == "r2e")) {
        default_colexternalOrder <- r2e_order_col - 1
    }
    else if (!(is.null(col.externalOrder)) && is.null(col.flip)) {
        default_colexternalOrder <- col.externalOrder
    }
    else {
        default_colexternalOrder <- seq(0, ncol(data_matrix) - 
            1)
    }
    if (row.order %in% c("single", "complete", "average")) {
        orderType <- orderType_map[[row.order]]
        flipType <- ifelse(is.null(row.flip), 0, flipType_map[[row.flip]])
        if (isFALSE(isProximityMatrix)) {
            res <- hctree_sort(distance_matrix_for_row_hc, default_rowexternalOrder, 
                orderType, flipType)
            hc <- list()
            hc$merge <- matrix(NA, nrow = nrow(distance_matrix_for_row_hc) - 
                1, ncol = 2)
            for (i in seq_along(res$left)) {
                hc$merge[i, 1] <- ifelse(res$left[i] < nrow(distance_matrix_for_row_hc), 
                  -res$left[i] - 1, res$left[i] - nrow(distance_matrix_for_row_hc) + 
                    1)
                hc$merge[i, 2] <- ifelse(res$right[i] < nrow(distance_matrix_for_row_hc), 
                  -res$right[i] - 1, res$right[i] - nrow(distance_matrix_for_row_hc) + 
                    1)
            }
        }
        else {
            res <- hctree_sort(data_matrix, default_rowexternalOrder, 
                orderType, flipType)
            hc <- list()
            hc$merge <- matrix(NA, nrow = nrow(data_matrix) - 
                1, ncol = 2)
            for (i in seq_along(res$left)) {
                hc$merge[i, 1] <- ifelse(res$left[i] < nrow(data_matrix), 
                  -res$left[i] - 1, res$left[i] - nrow(data_matrix) + 
                    1)
                hc$merge[i, 2] <- ifelse(res$right[i] < nrow(data_matrix), 
                  -res$right[i] - 1, res$right[i] - nrow(data_matrix) + 
                    1)
            }
        }
        hc$height <- res$height
        hc$order <- res$order + 1
        class(hc) <- "hclust"
        row_dend <- as.dendrogram(hc)
        row_names <- row_names[hc$order]
        row_names <- rev(row_names)
    }
    if (isFALSE(isProximityMatrix)) {
        if (col.order %in% c("single", "complete", "average")) {
            orderType <- orderType_map[[col.order]]
            flipType <- ifelse(is.null(col.flip), 0, flipType_map[[col.flip]])
            res <- hctree_sort(distance_matrix_for_column_hc, 
                default_colexternalOrder, orderType, flipType)
            hc <- list()
            hc$merge <- matrix(NA, nrow = nrow(distance_matrix_for_column_hc) - 
                1, ncol = 2)
            for (i in seq_along(res$left)) {
                hc$merge[i, 1] <- ifelse(res$left[i] < nrow(distance_matrix_for_column_hc), 
                  -res$left[i] - 1, res$left[i] - nrow(distance_matrix_for_column_hc) + 
                    1)
                hc$merge[i, 2] <- ifelse(res$right[i] < nrow(distance_matrix_for_column_hc), 
                  -res$right[i] - 1, res$right[i] - nrow(distance_matrix_for_column_hc) + 
                    1)
            }
            hc$height <- res$height
            hc$order <- res$order + 1
            class(hc) <- "hclust"
            column_dend <- as.dendrogram(hc)
            column_names <- column_names[hc$order]
        }
    }
    if (!is.null(Yc)) {
        for (i in seq_along(YcNum)) {
            if (row.order == "random") {
                assign(paste0("Yc_encoded_", i), as.numeric(Yc[, 
                  i])[random_order_row])
            }
            else if (row.order == "reverse") {
                assign(paste0("Yc_encoded_", i), as.numeric(Yc[, 
                  i])[reverse_order_row])
            }
            else if (row.order == "r2e") {
                assign(paste0("Yc_encoded_", i), as.numeric(Yc[, 
                  i])[r2e_order_row])
            }
            else if (row.order %in% c(seriation_dist, seriation_matrix)) {
                assign(paste0("Yc_encoded_", i), as.numeric(Yc[, 
                  i])[seriation_order_row])
            }
            else {
                assign(paste0("Yc_encoded_", i), as.numeric(Yc[, 
                  i]))
            }
        }
    }
    if (!is.null(Yd)) {
        Yd_encoding_maps <- list()
        list_Yd <- list()
        list_cb_Yd <- list()
        for (i in seq_along(YdNum)) {
            Yd_col <- Yd[, i]
            if (is.numeric(Yd_col) && all(sort(unique(Yd_col)) == 
                0:(length(unique(Yd_col)) - 1))) {
                Yd_numeric_encoding <- as.numeric(Yd_col)
                df_Yd <- data.frame(ordered_Yd = as.character(Yd_col))
                list_Yd[[i]] <- df_Yd
                Yd_encoding_map <- data.frame(Original = as.character(Yd_col), 
                  Encoding = Yd_numeric_encoding)
                Yd_encoding_map <- unique(Yd_encoding_map)
                rownames(Yd_encoding_map) <- NULL
                list_cb_Yd[[i]] <- Yd_encoding_map
                Yd_encoding_maps[[paste0("Yd_", i)]] <- Yd_encoding_map
            }
            else {
                Yd_factor_column <- factor(Yd_col, levels = unique(Yd_col))
                Yd_numeric_encoding <- as.numeric(Yd_factor_column) - 
                  1
                Yd_encoding_map_all <- data.frame(Original = as.character(Yd_factor_column), 
                  Encoding = Yd_numeric_encoding)
                df_Yd <- data.frame(ordered_Yd = Yd_encoding_map_all$Original)
                list_Yd[[i]] <- df_Yd
                Yd_encoding_map <- unique(Yd_encoding_map_all[complete.cases(Yd_encoding_map_all), 
                  ])
                rownames(Yd_encoding_map) <- NULL
                list_cb_Yd[[i]] <- Yd_encoding_map
                Yd_encoding_maps[[paste0("Yd_", i)]] <- Yd_encoding_map
            }
            if (row.order == "random") {
                assign(paste0("Yd_encoded_", i), Yd_numeric_encoding[random_order_row])
            }
            else if (row.order == "reverse") {
                assign(paste0("Yd_encoded_", i), Yd_numeric_encoding[reverse_order_row])
            }
            else if (row.order == "r2e") {
                assign(paste0("Yd_encoded_", i), Yd_numeric_encoding[r2e_order_row])
            }
            else if (row.order %in% c(seriation_dist, seriation_matrix)) {
                assign(paste0("Yd_encoded_", i), Yd_numeric_encoding[seriation_order_row])
            }
            else {
                assign(paste0("Yd_encoded_", i), Yd_numeric_encoding)
            }
        }
    }
    if (!is.null(Xc)) {
        for (i in seq_along(XcNum)) {
            if (col.order == "random") {
                assign(paste0("Xc_encoded_", i), as.numeric(Xc[i, 
                  ])[random_order_col])
            }
            else if (col.order == "reverse") {
                assign(paste0("Xc_encoded_", i), as.numeric(Xc[i, 
                  ])[reverse_order_col])
            }
            else if (col.order == "r2e") {
                assign(paste0("Xc_encoded_", i), as.numeric(Xc[i, 
                  ])[r2e_order_col])
            }
            else if (col.order %in% c(seriation_dist, seriation_matrix)) {
                assign(paste0("Xc_encoded_", i), as.numeric(Xc[i, 
                  ])[seriation_order_col])
            }
            else {
                assign(paste0("Xc_encoded_", i), as.numeric(Xc[i, 
                  ]))
            }
        }
    }
    if (!is.null(Xd)) {
        Xd_encoding_maps <- list()
        list_Xd <- list()
        list_cb_Xd <- list()
        for (i in seq_along(XdNum)) {
            Xd_row <- t(Xd[i, ])
            if (is.numeric(Xd_row) && all(sort(unique(Xd_row)) == 
                0:(length(unique(Xd_row)) - 1))) {
                Xd_numeric_encoding <- as.numeric(Xd_row)
                df_Xd <- data.frame(ordered_Xd = as.character(Xd_row))
                list_Xd[[i]] <- df_Xd
                Xd_encoding_map <- data.frame(Original = as.character(Xd_row), 
                  Encoding = Xd_numeric_encoding)
                Xd_encoding_map <- unique(Xd_encoding_map)
                rownames(Xd_encoding_map) <- NULL
                list_cb_Xd[[i]] <- Xd_encoding_map
                Xd_encoding_maps[[paste0("Xd_", i)]] <- Xd_encoding_map
            }
            else {
                Xd_factor_column <- factor(Xd_row, levels = unique(Xd_row))
                Xd_numeric_encoding <- as.numeric(Xd_factor_column) - 
                  1
                Xd_encoding_map_all <- data.frame(Original = as.character(Xd_factor_column), 
                  Encoding = Xd_numeric_encoding)
                df_Xd <- data.frame(ordered_Xd = Xd_encoding_map_all$Original)
                list_Xd[[i]] <- df_Xd
                Xd_encoding_map <- unique(Xd_encoding_map_all[complete.cases(Xd_encoding_map_all), 
                  ])
                rownames(Xd_encoding_map) <- NULL
                list_cb_Xd[[i]] <- Xd_encoding_map
                Xd_encoding_maps[[paste0("Xd_", i)]] <- Xd_encoding_map
            }
            if (col.order == "random") {
                assign(paste0("Xd_encoded_", i), Xd_numeric_encoding[random_order_col])
            }
            else if (col.order == "reverse") {
                assign(paste0("Xd_encoded_", i), Xd_numeric_encoding[reverse_order_col])
            }
            else if (col.order == "r2e") {
                assign(paste0("Xd_encoded_", i), Xd_numeric_encoding[r2e_order_col])
            }
            else if (col.order %in% c(seriation_dist, seriation_matrix)) {
                assign(paste0("Xd_encoded_", i), Xd_numeric_encoding[seriation_order_col])
            }
            else {
                assign(paste0("Xd_encoded_", i), Xd_numeric_encoding)
            }
        }
    }
    if (!(row.order %in% c("single", "complete", "average"))) {
        cluster_rows <- FALSE
    }
    else {
        cluster_rows <- row_dend
    }
    if (isFALSE(isProximityMatrix)) {
        if (!(col.order %in% c("single", "complete", "average"))) {
            cluster_columns <- FALSE
        }
        else {
            cluster_columns <- column_dend
        }
    }
    if (isTRUE(is_binary_data(data_matrix, TRUE))) {
        if (!is.null(original.color)) {
            default_originalcolor <- original.color
        }
        else {
            default_originalcolor <- grayscale_palette
        }
    }
    else {
        if (!is.null(original.color) && length(original.color) == 
            1 && original.color %in% rownames(brewer.pal.info)) {
            originalcolor_category <- brewer.pal.info[original.color, 
                "category"]
            if (originalcolor_category == "qual") {
                stop(sprintf("Error: '%s' is a qualitative color spectrum and cannot be used for continuous original data matrix.", 
                  original.color))
            }
            else {
                max_originalcolors <- brewer.pal.info[original.color, 
                  "maxcolors"]
                default_originalcolor <- brewer.pal(max_originalcolors, 
                  original.color)
            }
        }
        else if (!is.null(original.color) && length(original.color) == 
            1 && original.color %in% c("GAP_Rainbow", "GAP_Blue_White_Red")) {
            default_originalcolor <- get(original.color, envir = .GlobalEnv)
        }
        else if (!is.null(original.color) && !(original.color %in% 
            c(rownames(brewer.pal.info), "GAP_Rainbow", "GAP_Blue_White_Red"))) {
            default_originalcolor <- original.color
        }
        else {
            default_originalcolor <- GAP_Rainbow
        }
    }
    if (!is.null(row.prox)) {
        if (!is.null(row.color) && length(row.color) == 1 && 
            row.color %in% rownames(brewer.pal.info)) {
            rowcolor_category <- brewer.pal.info[row.color, "category"]
            if (rowcolor_category == "qual") {
                stop(sprintf("Error: '%s' is a qualitative color spectrum and cannot be used for proximity matrix.", 
                  row.color))
            }
            else if (!row.prox %in% c("euclidean", "city-block", 
                "maximum", "canberra", "hamman", "jaccard", "phi", 
                "rao", "rogers", "simple", "sneath", "yule") && 
                !row.color %in% c("BrBG", "PiYG", "PRGn", "PuOr", 
                  "RdBu", "RdGy")) {
                stop(sprintf("Error: '%s' cannot be used for proximity matrix.", 
                  row.color))
            }
            else {
                max_rowcolors <- brewer.pal.info[row.color, "maxcolors"]
                default_rowcolor <- brewer.pal(max_rowcolors, 
                  row.color)
            }
        }
        else if (!is.null(row.color) && length(row.color) == 
            1 && row.color %in% c("GAP_Rainbow", "GAP_Blue_White_Red")) {
            if (row.prox %in% c("euclidean", "city-block", "maximum", 
                "canberra", "hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule") && row.color == 
                "GAP_Blue_White_Red") {
                stop(sprintf("Error: The color spectrum 'GAP_Blue_White_Red' is not compatible with proximity methods such as '%s'. Please choose a different color spectrum or proximity measure.", 
                  row.prox))
            }
            else if (!row.prox %in% c("euclidean", "city-block", 
                "maximum", "canberra", "hamman", "jaccard", "phi", 
                "rao", "rogers", "simple", "sneath", "yule") && 
                row.color == "GAP_Rainbow") {
                stop(sprintf("Error: The color spectrum 'GAP_Rainbow' is not compatible with proximity methods such as '%s'. Please choose a different color spectrum or proximity measure.", 
                  row.prox))
            }
            else if (row.prox %in% c("euclidean", "city-block", 
                "maximum", "canberra", "hamman", "jaccard", "phi", 
                "rao", "rogers", "simple", "sneath", "yule") && 
                row.color == "GAP_Rainbow") {
                default_rowcolor <- GAP_Rainbow
            }
            else if (!row.prox %in% c("euclidean", "city-block", 
                "maximum", "canberra", "hamman", "jaccard", "phi", 
                "rao", "rogers", "simple", "sneath", "yule") && 
                row.color == "GAP_Blue_White_Red") {
                default_rowcolor <- GAP_Blue_White_Red_function
            }
        }
        else if (!is.null(row.color) && !(row.color %in% c(rownames(brewer.pal.info), 
            "GAP_Rainbow", "GAP_Blue_White_Red"))) {
            default_rowcolor <- row.color
        }
        else {
            if (row.prox %in% c("euclidean", "city-block", "maximum", 
                "canberra", "hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule")) {
                default_rowcolor <- GAP_Rainbow
            }
            else {
                default_rowcolor <- GAP_Blue_White_Red_function
            }
        }
    }
    if (!is.null(col.prox)) {
        if (!is.null(col.color) && length(col.color) == 1 && 
            col.color %in% rownames(brewer.pal.info)) {
            colcolor_category <- brewer.pal.info[col.color, "category"]
            if (colcolor_category == "qual") {
                stop(sprintf("Error: '%s' is a qualitative color spectrum and cannot be used for proximity matrix.", 
                  col.color))
            }
            else if (!col.prox %in% c("euclidean", "city-block", 
                "maximum", "canberra", "hamman", "jaccard", "phi", 
                "rao", "rogers", "simple", "sneath", "yule") && 
                !col.color %in% c("BrBG", "PiYG", "PRGn", "PuOr", 
                  "RdBu", "RdGy")) {
                stop(sprintf("Error: '%s' cannot be used for proximity matrix.", 
                  col.color))
            }
            else {
                max_colcolors <- brewer.pal.info[col.color, "maxcolors"]
                default_colcolor <- brewer.pal(max_colcolors, 
                  col.color)
            }
        }
        else if (!is.null(col.color) && length(col.color) == 
            1 && col.color %in% c("GAP_Rainbow", "GAP_Blue_White_Red")) {
            if (col.prox %in% c("euclidean", "city-block", "maximum", 
                "canberra", "hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule") && col.color == 
                "GAP_Blue_White_Red") {
                stop(sprintf("Error: The color spectrum 'GAP_Blue_White_Red' is not compatible with proximity methods such as '%s'. Please choose a different color spectrum or proximity measure.", 
                  col.prox))
            }
            else if (!col.prox %in% c("euclidean", "city-block", 
                "maximum", "canberra", "hamman", "jaccard", "phi", 
                "rao", "rogers", "simple", "sneath", "yule") && 
                col.color == "GAP_Rainbow") {
                stop(sprintf("Error: The color spectrum 'GAP_Rainbow' is not compatible with proximity methods such as '%s'. Please choose a different color spectrum or proximity measure.", 
                  col.prox))
            }
            else if (col.prox %in% c("euclidean", "city-block", 
                "maximum", "canberra", "hamman", "jaccard", "phi", 
                "rao", "rogers", "simple", "sneath", "yule") && 
                col.color == "GAP_Rainbow") {
                default_colcolor <- GAP_Rainbow
            }
            else if (!col.prox %in% c("euclidean", "city-block", 
                "maximum", "canberra", "hamman", "jaccard", "phi", 
                "rao", "rogers", "simple", "sneath", "yule") && 
                col.color == "GAP_Blue_White_Red") {
                default_colcolor <- GAP_Blue_White_Red_function
            }
        }
        else if (!is.null(col.color) && !(col.color %in% c(rownames(brewer.pal.info), 
            "GAP_Rainbow", "GAP_Blue_White_Red"))) {
            default_colcolor <- col.color
        }
        else {
            if (col.prox %in% c("euclidean", "city-block", "maximum", 
                "canberra", "hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule")) {
                default_colcolor <- GAP_Rainbow
            }
            else {
                default_colcolor <- GAP_Blue_White_Red_function
            }
        }
    }
    if (is.null(original.color)) {
        if (isTRUE(is_binary_data(data_matrix, TRUE))) {
            original.color <- "grayscale_palette"
        }
        else {
            original.color <- "GAP_Rainbow"
        }
    }
    original_colorbar_info <- colorbar_fun(data_matrix, original.color)
    if (isTRUE(is_binary_data(data_matrix, TRUE))) {
        original_colorbar <- Legend(title = "data matrix", at = c(0, 
            1), labels = c("0", "1"), legend_gp = gpar(fill = c("white", 
            "black")), direction = "horizontal", grid_height = unit(0.5, 
            "cm"), grid_width = unit(0.5, "cm"), title_gp = gpar(fontsize = 10, 
            fontface = "bold"), labels_gp = gpar(fontsize = 8), 
            tick_length = unit(0, "mm"), border = "#00000000")
    }
    else {
        original_colorbar <- Legend(title = if (isFALSE(isProximityMatrix)) 
            "data matrix"
        else "data matrix (proximity)", at = c(original_colorbar_info$colorbar_min, 
            original_colorbar_info$colorbar_mid, original_colorbar_info$colorbar_max), 
            labels = c(original_colorbar_info$colorbar_min, original_colorbar_info$colorbar_mid, 
                original_colorbar_info$colorbar_max), col_fun = original_colorbar_info$colorbar_function, 
            direction = "horizontal", legend_width = unit(4, 
                "cm"), legend_height = unit(0.8, "cm"), title_gp = gpar(fontsize = 10, 
                fontface = "bold"), labels_gp = gpar(fontsize = 8), 
            tick_length = unit(0, "mm"), border = "#00000000")
    }
    if (!is.null(row.prox)) {
        if (is.null(row.color)) {
            if (row.prox %in% c("euclidean", "city-block", "maximum", 
                "canberra", "hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule")) {
                row.color <- "GAP_Rainbow"
            }
            else {
                row.color <- "GAP_Blue_White_Red"
            }
        }
        row_colorbar_info <- colorbar_fun(distance_matrix_for_row_hc, 
            row.color)
        row_colorbar <- Legend(title = "row proximity matrix", 
            at = c(row_colorbar_info$colorbar_min, row_colorbar_info$colorbar_mid, 
                row_colorbar_info$colorbar_max), labels = c(row_colorbar_info$colorbar_min, 
                row_colorbar_info$colorbar_mid, row_colorbar_info$colorbar_max), 
            col_fun = row_colorbar_info$colorbar_function, direction = "horizontal", 
            legend_width = unit(4, "cm"), legend_height = unit(0.8, 
                "cm"), title_gp = gpar(fontsize = 10, fontface = "bold"), 
            labels_gp = gpar(fontsize = 8), tick_length = unit(0, 
                "mm"), border = "#00000000")
    }
    if (!is.null(col.prox)) {
        if (is.null(col.color)) {
            if (col.prox %in% c("euclidean", "city-block", "maximum", 
                "canberra", "hamman", "jaccard", "phi", "rao", 
                "rogers", "simple", "sneath", "yule")) {
                col.color <- "GAP_Rainbow"
            }
            else {
                col.color <- "GAP_Blue_White_Red"
            }
        }
        col_colorbar_info <- colorbar_fun(distance_matrix_for_column_hc, 
            col.color)
        col_colorbar <- Legend(title = "column proximity matrix", 
            at = c(col_colorbar_info$colorbar_min, col_colorbar_info$colorbar_mid, 
                col_colorbar_info$colorbar_max), labels = c(col_colorbar_info$colorbar_min, 
                col_colorbar_info$colorbar_mid, col_colorbar_info$colorbar_max), 
            col_fun = col_colorbar_info$colorbar_function, direction = "horizontal", 
            legend_width = unit(4, "cm"), legend_height = unit(0.8, 
                "cm"), title_gp = gpar(fontsize = 10, fontface = "bold"), 
            labels_gp = gpar(fontsize = 8), tick_length = unit(0, 
                "mm"), border = "#00000000")
    }
    colorbar_list <- append(colorbar_list, list(original_colorbar))
    if (isTRUE(show.row.prox)) {
        colorbar_list <- append(colorbar_list, list(row_colorbar))
    }
    if (isTRUE(show.col.prox)) {
        colorbar_list <- append(colorbar_list, list(col_colorbar))
    }
    ht_opt$message = FALSE
    original_plot <- Heatmap(data_matrix, show_heatmap_legend = FALSE, 
        cluster_rows = cluster_rows, cluster_columns = if (isFALSE(isProximityMatrix)) 
            cluster_columns
        else cluster_rows, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_row_names = FALSE, show_column_names = FALSE, col = default_originalcolor, 
        width = if (isFALSE(isProximityMatrix)) 
            unit(3, "cm")
        else unit(8, "cm"), height = unit(8, "cm"), na_col = MissingValue.color)
    if (!is.null(row.prox) || isTRUE(isProximityMatrix)) {
        row_prox <- Heatmap(if (isTRUE(isProximityMatrix)) 
            data_matrix
        else distance_matrix_for_row, show_heatmap_legend = FALSE, 
            cluster_rows = cluster_rows, cluster_columns = cluster_rows, 
            show_row_dend = TRUE, show_column_dend = FALSE, row_dend_side = "right", 
            show_row_names = FALSE, show_column_names = FALSE, 
            col = if (!isTRUE(isProximityMatrix) && isTRUE(show.row.prox)) 
                default_rowcolor
            else NULL, rect_gp = if (!isTRUE(isProximityMatrix) && 
                isTRUE(show.row.prox)) 
                gpar(col = NA)
            else gpar(type = "none"), width = if (!isTRUE(isProximityMatrix) && 
                isTRUE(show.row.prox)) 
                unit(8, "cm")
            else unit(0, "cm"), height = unit(8, "cm"), na_col = MissingValue.color)
    }
    if (!is.null(col.prox)) {
        col_prox <- Heatmap(distance_matrix_for_column, show_heatmap_legend = FALSE, 
            cluster_rows = cluster_columns, cluster_columns = cluster_columns, 
            show_row_dend = FALSE, show_column_dend = TRUE, column_dend_side = "top", 
            show_row_names = FALSE, show_column_names = FALSE, 
            col = if (isTRUE(show.col.prox)) 
                default_colcolor
            else NULL, rect_gp = if (isTRUE(show.col.prox)) 
                gpar(col = NA)
            else gpar(type = "none"), height = if (isTRUE(show.col.prox)) 
                unit(3, "cm")
            else unit(0, "cm"), width = unit(3, "cm"), column_dend_height = unit(1, 
                "cm"), na_col = MissingValue.color)
    }
    if (!is.null(XcNum)) {
        Xc_plots <- list()
        Xc_colorbars <- list()
        for (i in seq_along(XcNum)) {
            Xc_encoded <- get(paste0("Xc_encoded_", i))
            if (!is.null(Xc.color)) {
                if (length(Xc.color) == 1 && Xc.color %in% rownames(brewer.pal.info)) {
                  Xccolor_category <- brewer.pal.info[Xc.color, 
                    "category"]
                  if (Xccolor_category == "qual") {
                    stop(sprintf("Error: '%s' is a qualitative color spectrum and cannot be used for continuous data.", 
                      Xc.color))
                  }
                  else {
                    max_Xccolors <- brewer.pal.info[Xc.color, 
                      "maxcolors"]
                    default_Xccolor <- brewer.pal(max_Xccolors, 
                      Xc.color)
                  }
                }
                else if (length(Xc.color) == 1 && Xc.color %in% 
                  c("GAP_Rainbow", "GAP_Blue_White_Red")) {
                  default_Xccolor <- get(Xc.color, envir = .GlobalEnv)
                }
                else {
                  default_Xccolor <- Xc.color
                }
            }
            else {
                default_Xccolor <- GAP_Rainbow
                Xc.color <- "GAP_Rainbow"
            }
            Xc_colorbar_info <- colorbar_fun(Xc_encoded, Xc.color)
            Xc_colorbars[[i]] <- Legend(title = paste("Xc matrix of", 
                Xc.name[i]), at = c(Xc_colorbar_info$colorbar_min, 
                Xc_colorbar_info$colorbar_mid, Xc_colorbar_info$colorbar_max), 
                labels = c(Xc_colorbar_info$colorbar_min, Xc_colorbar_info$colorbar_mid, 
                  Xc_colorbar_info$colorbar_max), col_fun = Xc_colorbar_info$colorbar_function, 
                direction = "horizontal", legend_width = unit(4, 
                  "cm"), legend_height = unit(0.8, "cm"), title_gp = gpar(fontsize = 10, 
                  fontface = "bold"), labels_gp = gpar(fontsize = 8), 
                tick_length = unit(0, "mm"), border = "#00000000")
            Xc_plots[[i]] <- Heatmap(t(Xc_encoded), show_heatmap_legend = FALSE, 
                cluster_rows = FALSE, cluster_columns = cluster_columns, 
                show_row_dend = FALSE, show_column_dend = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE, 
                col = default_Xccolor, width = unit(3, "cm"), 
                height = unit(0.5, "cm"), na_col = MissingValue.color)
        }
        colorbar_list <- c(colorbar_list, Xc_colorbars)
    }
    if (!is.null(YcNum)) {
        Yc_plots <- list()
        Yc_colorbars <- list()
        for (i in seq_along(YcNum)) {
            Yc_encoded <- get(paste0("Yc_encoded_", i))
            if (!is.null(Yc.color)) {
                if (length(Yc.color) == 1 && Yc.color %in% rownames(brewer.pal.info)) {
                  Yccolor_category <- brewer.pal.info[Yc.color, 
                    "category"]
                  if (Yccolor_category == "qual") {
                    stop(sprintf("Error: '%s' is a qualitative color spectrum and cannot be used for continuous data.", 
                      Yc.color))
                  }
                  else {
                    max_Yccolors <- brewer.pal.info[Yc.color, 
                      "maxcolors"]
                    default_Yccolor <- brewer.pal(max_Yccolors, 
                      Yc.color)
                  }
                }
                else if (length(Yc.color) == 1 && Yc.color %in% 
                  c("GAP_Rainbow", "GAP_Blue_White_Red")) {
                  default_Yccolor <- get(Yc.color, envir = .GlobalEnv)
                }
                else {
                  default_Yccolor <- Yc.color
                }
            }
            else {
                default_Yccolor <- GAP_Rainbow
                Yc.color <- "GAP_Rainbow"
            }
            Yc_colorbar_info <- colorbar_fun(Yc_encoded, Yc.color)
            Yc_colorbars[[i]] <- Legend(title = paste("Yc matrix of", 
                colnames(data)[YcNum[i]]), at = c(Yc_colorbar_info$colorbar_min, 
                Yc_colorbar_info$colorbar_mid, Yc_colorbar_info$colorbar_max), 
                labels = c(Yc_colorbar_info$colorbar_min, Yc_colorbar_info$colorbar_mid, 
                  Yc_colorbar_info$colorbar_max), col_fun = Yc_colorbar_info$colorbar_function, 
                direction = "horizontal", legend_width = unit(4, 
                  "cm"), legend_height = unit(0.8, "cm"), title_gp = gpar(fontsize = 10, 
                  fontface = "bold"), labels_gp = gpar(fontsize = 8), 
                tick_length = unit(0, "mm"), border = "#00000000")
            Yc_plots[[i]] <- Heatmap(Yc_encoded, show_heatmap_legend = FALSE, 
                cluster_rows = cluster_rows, cluster_columns = FALSE, 
                show_row_dend = FALSE, show_column_dend = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE, 
                col = default_Yccolor, width = unit(0.5, "cm"), 
                height = unit(8, "cm"), na_col = MissingValue.color)
        }
        colorbar_list <- c(colorbar_list, Yc_colorbars)
    }
    if (!is.null(XdNum)) {
        Xd_plots <- list()
        Xd_colorbar <- list()
        for (i in seq_along(XdNum)) {
            Xd_encoded <- get(paste0("Xd_encoded_", i))
            if (!is.null(Xd.color)) {
                if (length(Xd.color) == 1 && Xd.color %in% c("GAP_Rainbow", 
                  "GAP_Blue_White_Red")) {
                  stop(sprintf("Error: '%s' cannot be used for categorical data. Please use a qualitative color spectrum instead.", 
                    Xd.color))
                }
                else if (length(Xd.color) == 1 && Xd.color %in% 
                  rownames(brewer.pal.info)) {
                  Xdcolor_category <- brewer.pal.info[Xd.color, 
                    "category"]
                  if (Xdcolor_category == "qual") {
                    max_Xdcolors <- brewer.pal.info[Xd.color, 
                      "maxcolors"]
                    if (length(unique(Xd_encoded)) > max_Xdcolors) {
                      warning(sprintf("Warning: Number of unique categories (%d) exceeds the maximum colors in '%s' spectrum (%d). Colors will be recycled.", 
                        length(unique(Xd_encoded)), Xd.color, 
                        max_Xdcolors))
                      default_Xdcolor <- rep(brewer.pal(max_Xdcolors, 
                        Xd.color), length.out = length(unique(Xd_encoded)))
                    }
                    else {
                      default_Xdcolor <- brewer.pal(max_Xdcolors, 
                        Xd.color)[1:length(unique(Xd_encoded))]
                    }
                  }
                  else {
                    stop(sprintf("Error: '%s' is a %s color spectrum and cannot be used for categorical data. Please use a qualitative color spectrum instead.", 
                      Xd.color, ifelse(Xdcolor_category == "seq", 
                        "Sequential", "Diverging")))
                  }
                }
                else if (length(Xd.color) == 1 && Xd.color == 
                  "GAP_d") {
                  default_Xdcolor <- GAP_d[1:length(unique(Xd_encoded))]
                }
                else {
                  default_Xdcolor <- Xd.color
                }
            }
            else {
                default_Xdcolor <- GAP_d[1:length(unique(Xd_encoded))]
                if (length(unique(Xd_encoded)) > 16) {
                  warning(sprintf("Warning: Number of unique categories (%d) exceeds the maximum colors in default color spectrum (%d).", 
                    length(unique(Xd_encoded)), 16))
                }
            }
            Xd_label <- Xd_encoding_maps[[paste0("Xd_", i)]]
            Xd_colorbar[[i]] <- Legend(title = paste("Xd categories of", 
                Xd.name[i]), at = Xd_label$Encoding, labels = Xd_label$Original, 
                legend_gp = gpar(fill = default_Xdcolor), direction = "horizontal", 
                grid_height = unit(0.5, "cm"), grid_width = unit(0.5, 
                  "cm"), title_gp = gpar(fontsize = 10, fontface = "bold"), 
                labels_gp = gpar(fontsize = 8), tick_length = unit(0, 
                  "mm"), border = "#00000000")
            Xd_plots[[i]] <- Heatmap(t(Xd_encoded), show_heatmap_legend = FALSE, 
                cluster_rows = FALSE, cluster_columns = cluster_columns, 
                show_row_dend = FALSE, show_column_dend = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE, 
                col = default_Xdcolor, width = unit(3, "cm"), 
                height = unit(0.5, "cm"), row_dend_width = unit(1, 
                  "cm"), na_col = MissingValue.color)
        }
        colorbar_list <- c(colorbar_list, Xd_colorbar)
    }
    if (!is.null(YdNum)) {
        Yd_plots <- list()
        Yd_colorbar <- list()
        for (i in seq_along(YdNum)) {
            Yd_encoded <- get(paste0("Yd_encoded_", i))
            if (!is.null(Yd.color)) {
                if (length(Xd.color) == 1 && Xd.color %in% c("GAP_Rainbow", 
                  "GAP_Blue_White_Red")) {
                  stop(sprintf("Error: '%s' cannot be used for categorical data. Please use a qualitative color spectrum instead.", 
                    Yd.color))
                }
                else if (length(Yd.color) == 1 && Yd.color %in% 
                  rownames(brewer.pal.info)) {
                  Ydcolor_category <- brewer.pal.info[Yd.color, 
                    "category"]
                  if (Ydcolor_category == "qual") {
                    max_Ydcolors <- brewer.pal.info[Yd.color, 
                      "maxcolors"]
                    if (length(unique(Yd_encoded)) > max_Ydcolors) {
                      warning(sprintf("Warning: Number of unique categories (%d) exceeds the maximum colors in '%s' spectrum (%d). Colors will be recycled.", 
                        length(unique(Yd_encoded)), Yd.color, 
                        max_Ydcolors))
                      default_Ydcolor <- rep(brewer.pal(max_Ydcolors, 
                        Yd.color), length.out = length(unique(Yd_encoded)))
                    }
                    else {
                      default_Ydcolor <- brewer.pal(max_Ydcolors, 
                        Yd.color)[1:length(unique(Yd_encoded))]
                    }
                  }
                  else {
                    stop(sprintf("Error: '%s' is a %s color spectrum and cannot be used for categorical data. Please use a qualitative color spectrum instead.", 
                      Yd.color, ifelse(Ydcolor_category == "seq", 
                        "Sequential", "Diverging")))
                  }
                }
                else if (length(Yd.color) == 1 && Yd.color == 
                  "GAP_d") {
                  default_Ydcolor <- GAP_d[1:length(unique(Yd_encoded))]
                }
                else {
                  default_Ydcolor <- Yd.color
                }
            }
            else {
                default_Ydcolor <- GAP_d[1:length(unique(Yd_encoded))]
                if (length(unique(Yd_encoded)) > 16) {
                  warning(sprintf("Warning: Number of unique categories (%d) exceeds the maximum colors in default color spectrum (%d).", 
                    length(unique(Yd_encoded)), 16))
                }
            }
            Yd_label <- Yd_encoding_maps[[paste0("Yd_", i)]]
            Yd_colorbar[[i]] <- Legend(title = paste("Yd categories of", 
                colnames(data)[YdNum[i]]), at = Yd_label$Encoding, 
                labels = Yd_label$Original, legend_gp = gpar(fill = default_Ydcolor), 
                direction = "horizontal", grid_height = unit(0.5, 
                  "cm"), grid_width = unit(0.5, "cm"), title_gp = gpar(fontsize = 10, 
                  fontface = "bold"), labels_gp = gpar(fontsize = 8), 
                tick_length = unit(0, "mm"), border = "#00000000")
            Yd_plots[[i]] <- Heatmap(Yd_encoded, show_heatmap_legend = FALSE, 
                cluster_rows = cluster_rows, cluster_columns = FALSE, 
                show_row_dend = FALSE, show_column_dend = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE, 
                col = default_Ydcolor, width = unit(0.5, "cm"), 
                height = unit(8, "cm"), na_col = MissingValue.color)
        }
        colorbar_list <- c(colorbar_list, Yd_colorbar)
    }
    grid.newpage()
    pushViewport(viewport())
    has_value <- function(...) {
        any(!sapply(list(...), is.null))
    }
    all_have_value <- function(...) {
        all(!sapply(list(...), is.null))
    }
    if (has_value(YdNum, YcNum) && !has_value(XdNum, XcNum)) {
        if (is.null(row.prox) && is.null(col.prox)) {
            layout_row <- 1
            layout_col <- if (all_have_value(YdNum, YcNum)) 
                3
            else 2
        }
        else if (!is.null(row.prox) && is.null(col.prox)) {
            layout_row <- 1
            layout_col <- if (all_have_value(YdNum, YcNum)) 
                4
            else 3
        }
        else if (is.null(row.prox) && !is.null(col.prox)) {
            layout_row <- 2
            layout_col <- if (all_have_value(YdNum, YcNum)) 
                3
            else 2
        }
        else {
            layout_row <- 2
            layout_col <- if (all_have_value(YdNum, YcNum)) 
                4
            else 3
        }
    }
    else if (!has_value(YdNum, YcNum) && has_value(XdNum, XcNum)) {
        if (is.null(row.prox) && is.null(col.prox)) {
            layout_row <- if (all_have_value(XdNum, XcNum)) 
                3
            else 2
            layout_col <- 1
        }
        else if (!is.null(row.prox) && is.null(col.prox)) {
            layout_row <- if (all_have_value(XdNum, XcNum)) 
                3
            else 2
            layout_col <- 2
        }
        else if (is.null(row.prox) && !is.null(col.prox)) {
            layout_row <- if (all_have_value(XdNum, XcNum)) 
                4
            else 3
            layout_col <- 1
        }
        else {
            layout_row <- if (all_have_value(XdNum, XcNum)) 
                4
            else 3
            layout_col <- 2
        }
    }
    else if (has_value(YdNum, YcNum) && has_value(XdNum, XcNum)) {
        if (is.null(row.prox) && is.null(col.prox)) {
            layout_row <- if (all_have_value(XdNum, XcNum)) 
                3
            else 2
            layout_col <- if (all_have_value(YdNum, YcNum)) 
                3
            else 2
        }
        else if (!is.null(row.prox) && is.null(col.prox)) {
            layout_row <- if (all_have_value(XdNum, XcNum)) 
                3
            else 2
            layout_col <- if (all_have_value(YdNum, YcNum)) 
                4
            else 3
        }
        else if (is.null(row.prox) && !is.null(col.prox)) {
            layout_row <- if (all_have_value(XdNum, XcNum)) 
                4
            else 3
            layout_col <- if (all_have_value(YdNum, YcNum)) 
                3
            else 2
        }
        else {
            layout_row <- if (all_have_value(XdNum, XcNum)) 
                4
            else 3
            layout_col <- if (all_have_value(YdNum, YcNum)) 
                4
            else 3
        }
    }
    else if (!has_value(YdNum, YcNum) && !has_value(XdNum, XcNum)) {
        if (is.null(row.prox) && is.null(col.prox)) {
            if (isFALSE(isProximityMatrix)) {
                layout_row <- 1
                layout_col <- 1
            }
            else {
                layout_row <- 1
                layout_col <- 2
            }
        }
        else if (!is.null(row.prox) && is.null(col.prox)) {
            layout_row <- 1
            layout_col <- 2
        }
        else if (is.null(row.prox) && !is.null(col.prox)) {
            layout_row <- 2
            layout_col <- 1
        }
        else {
            layout_row <- 2
            layout_col <- 2
        }
    }
    else {
        stop("Unexpected input")
    }
    if (!is.null(col.prox) && isTRUE(col.order %in% c("single", 
        "complete", "average")) && isFALSE(show.col.prox)) {
        layout_row <- layout_row + 1
    }
    if (has_value(YdNum) && has_value(YcNum)) {
        if (is.null(row.prox)) {
            widths <- unit(c(length(Yd_plots) * 0.5 + 0.5, length(Yc_plots) * 
                0.5 + 0.5, 3.5), "cm")
        }
        else {
            if (isFALSE(show.row.prox)) {
                widths <- unit(c(length(Yd_plots) * 0.5 + 0.5, 
                  length(Yc_plots) * 0.5 + 0.5, 3.5, 1), "cm")
            }
            else {
                widths <- unit(c(length(Yd_plots) * 0.5 + 0.5, 
                  length(Yc_plots) * 0.5 + 0.5, 3.5, 9.5), "cm")
            }
        }
    }
    else if (has_value(YdNum)) {
        if (is.null(row.prox)) {
            widths <- unit(c(length(Yd_plots) * 0.5 + 0.5, 3.5), 
                "cm")
        }
        else {
            if (isFALSE(show.row.prox)) {
                widths <- unit(c(length(Yd_plots) * 0.5 + 0.5, 
                  3.5, 1), "cm")
            }
            else {
                widths <- unit(c(length(Yd_plots) * 0.5 + 0.5, 
                  3.5, 9.5), "cm")
            }
        }
    }
    else if (has_value(YcNum)) {
        if (is.null(row.prox)) {
            widths <- unit(c(length(Yc_plots) * 0.5 + 0.5, 3.5), 
                "cm")
        }
        else {
            if (isFALSE(show.row.prox)) {
                widths <- unit(c(length(Yc_plots) * 0.5 + 0.5, 
                  3.5, 1), "cm")
            }
            else {
                widths <- unit(c(length(Yc_plots) * 0.5 + 0.5, 
                  3.5, 9.5), "cm")
            }
        }
    }
    else {
        if (is.null(row.prox)) {
            if (isFALSE(isProximityMatrix)) {
                widths <- unit(c(3.5), "cm")
            }
            else {
                widths <- unit(c(8.5, 1), "cm")
            }
        }
        else {
            if (isFALSE(show.row.prox)) {
                widths <- unit(c(3.5, 1), "cm")
            }
            else {
                widths <- unit(c(3.5, 9.5), "cm")
            }
        }
    }
    pushViewport(viewport(gp = gpar(fontsize = ifelse(is.null(col.label.size), 
        8, col.label.size))))
    column_names_length <- sapply(column_names, function(label) {
        convertUnit(grobWidth(textGrob(label)), "cm", valueOnly = TRUE)
    })
    max_column_names <- which.max(column_names_length)
    max_column_names_cm <- column_names_length[max_column_names]
    if (has_value(XdNum) && has_value(XcNum)) {
        if (is.null(col.prox)) {
            heights <- unit(c(length(Xc_plots) * 0.5 + 0.5, length(Xd_plots) * 
                0.5 + 0.5, 8.5), "cm")
        }
        else {
            if (isFALSE(show.col.prox)) {
                if (!is.null(col.prox) && isTRUE(col.order %in% 
                  c("single", "complete", "average"))) {
                  heights <- unit(c(1.2, max_column_names_cm, 
                    length(Xc_plots) * 0.5 + 0.5, length(Xd_plots) * 
                      0.5 + 0.5, 8.5), "cm")
                }
                else {
                  heights <- unit(c(1, length(Xc_plots) * 0.5 + 
                    0.5, length(Xd_plots) * 0.5 + 0.5, 8.5), 
                    "cm")
                }
            }
            else {
                heights <- unit(c(4.5, length(Xc_plots) * 0.5 + 
                  0.5, length(Xd_plots) * 0.5 + 0.5, 8.5), "cm")
            }
        }
    }
    else if (has_value(XdNum)) {
        if (is.null(col.prox)) {
            heights <- unit(c(length(Xd_plots) * 0.5 + 0.5, 8.5), 
                "cm")
        }
        else {
            if (isFALSE(show.col.prox)) {
                if (!is.null(col.prox) && isTRUE(col.order %in% 
                  c("single", "complete", "average"))) {
                  heights <- unit(c(1.2, max_column_names_cm, 
                    length(Xd_plots) * 0.5 + 0.5, 8.5), "cm")
                }
                else {
                  heights <- unit(c(1, length(Xd_plots) * 0.5 + 
                    0.5, 8.5), "cm")
                }
            }
            else {
                heights <- unit(c(4.5, length(Xd_plots) * 0.5 + 
                  0.5, 8.5), "cm")
            }
        }
    }
    else if (has_value(XcNum)) {
        if (is.null(col.prox)) {
            heights <- unit(c(length(Xc_plots) * 0.5 + 0.5, 8.5), 
                "cm")
        }
        else {
            if (isFALSE(show.col.prox)) {
                if (!is.null(col.prox) && isTRUE(col.order %in% 
                  c("single", "complete", "average"))) {
                  heights <- unit(c(1.2, max_column_names_cm, 
                    length(Xc_plots) * 0.5 + 0.5, 8.5), "cm")
                }
                else {
                  heights <- unit(c(1, length(Xc_plots) * 0.5 + 
                    0.5, 8.5), "cm")
                }
            }
            else {
                heights <- unit(c(4.5, length(Xc_plots) * 0.5 + 
                  0.5, 8.5), "cm")
            }
        }
    }
    else {
        if (is.null(col.prox)) {
            if (isFALSE(isProximityMatrix)) {
                heights <- unit(c(8.5), "cm")
            }
            else {
                heights <- unit(c(8.5), "cm")
            }
        }
        else {
            if (isFALSE(show.col.prox)) {
                if (!is.null(col.prox) && isTRUE(col.order %in% 
                  c("single", "complete", "average"))) {
                  heights <- unit(c(1.2, max_column_names_cm, 
                    8.5), "cm")
                }
                else {
                  heights <- unit(c(1, 8.5), "cm")
                }
            }
            else {
                heights <- unit(c(4.5, 8.5), "cm")
            }
        }
    }
    while (dev.cur() > 1) dev.off()
    if (is.null(PNGfilename)) {
        PNGfilename <- file.path(tempdir(), "output_plot.png")
    }
    png(filename = PNGfilename, width = PNGwidth, height = PNGheight, 
        res = PNGres, units = "px")
    pushViewport(viewport(layout = grid.layout(layout_row, layout_col, 
        widths = widths, heights = heights)))
    add_border <- function(plot_ht) {
        heatmap_name <- plot_ht@ht_list[[1]]@name
        decorate_heatmap_body(heatmap_name, {
            grid.rect(gp = gpar(col = "black", fill = NA, lwd = border.width))
        })
    }
    pushViewport(viewport(layout.pos.row = layout_row, layout.pos.col = if (is.null(row.prox) && 
        isFALSE(isProximityMatrix)) 
        layout_col
    else layout_col - 1))
    original_plot_ht <- draw(original_plot, newpage = FALSE)
    if (isTRUE(border)) {
        add_border(original_plot_ht)
    }
    popViewport()
    if (!is.null(row.prox) || isTRUE(isProximityMatrix)) {
        pushViewport(viewport(layout.pos.row = layout_row, layout.pos.col = layout_col))
        row_prox_ht <- draw(row_prox, newpage = FALSE)
        if (isTRUE(border) && isTRUE(show.row.prox)) {
            add_border(row_prox_ht)
        }
        popViewport()
    }
    pushViewport(viewport(layout.pos.row = layout_row, layout.pos.col = 1))
    for (i in seq_along(row_names)) {
        grid.text(row_names[i], x = unit(0, "cm"), y = unit(0.3 + 
            (8/length(row_names)) * (i - 0.5), "cm"), just = "right", 
            gp = gpar(fontsize = ifelse(is.null(row.label.size), 
                2, row.label.size)))
    }
    popViewport()
    if (isFALSE(isProximityMatrix)) {
        if (!is.null(col.prox)) {
            pushViewport(viewport(layout.pos.row = 1, layout.pos.col = if (is.null(row.prox) && 
                isFALSE(isProximityMatrix)) 
                layout_col
            else layout_col - 1))
            col_prox_ht <- draw(col_prox, newpage = FALSE)
            if (isTRUE(border) && isTRUE(show.col.prox)) {
                add_border(col_prox_ht)
            }
            popViewport()
            if (isTRUE(show.col.prox)) {
                pushViewport(viewport(layout.pos.row = 1, layout.pos.col = layout_col))
                for (i in seq_along(column_names)) {
                  grid.text(column_names[i], x = if (!is.null(row.prox)) 
                    unit(0, "cm")
                  else unit(1, "npc"), y = unit(if (col.order %in% 
                    c("original", "random", "reverse", "r2e", 
                      seriation_dist, seriation_matrix)) {
                    3.75 - (i - 0.5) * (3/length(column_names))
                  }
                  else {
                    3.25 - (i - 0.5) * (3/length(column_names))
                  }, "cm"), just = "left", gp = gpar(fontsize = ifelse(is.null(col.label.size), 
                    8, col.label.size)))
                }
                popViewport()
            }
            else {
                pushViewport(viewport(layout.pos.row = if (!is.null(col.prox) && 
                  isTRUE(col.order %in% c("single", "complete", 
                    "average"))) 
                  2
                else 1, layout.pos.col = if (is.null(row.prox) && 
                  isFALSE(isProximityMatrix)) 
                  layout_col
                else layout_col - 1))
                for (i in seq_along(column_names)) {
                  grid.text(column_names[i], x = unit(0.2 + ((i - 
                    0.5) * (3/length(column_names))), "cm"), 
                    y = unit(0, "cm"), just = "left", gp = gpar(fontsize = ifelse(is.null(col.label.size), 
                      8, col.label.size)), rot = 90)
                }
                popViewport()
            }
        }
        else {
            pushViewport(viewport(layout.pos.row = 1, layout.pos.col = if (is.null(row.prox) && 
                isFALSE(isProximityMatrix)) 
                layout_col
            else layout_col - 1))
            for (i in seq_along(column_names)) {
                grid.text(column_names[i], x = unit(0.2 + ((i - 
                  0.5) * (3/length(column_names))), "cm"), y = unit(1, 
                  "npc"), just = "left", gp = gpar(fontsize = ifelse(is.null(col.label.size), 
                  8, col.label.size)), rot = 90)
            }
            popViewport()
        }
    }
    if (!is.null(YcNum)) {
        list_Yc <- list()
        if (length(Yc_plots) == 1) {
            pushViewport(viewport(layout.pos.row = layout_row, 
                layout.pos.col = if (!is.null(YdNum)) 
                  2
                else 1))
            Yc_plots_ht <- draw(Yc_plots[[1]], newpage = FALSE)
            if (isTRUE(border)) {
                add_border(Yc_plots_ht)
            }
            popViewport()
        }
        else {
            pushViewport(viewport(layout.pos.row = layout_row, 
                layout.pos.col = if (!is.null(YdNum)) 
                  2
                else 1))
            Yc_plots_ht <- draw(Reduce(`+`, Yc_plots), newpage = FALSE, 
                gap = unit(0, "cm"))
            if (isTRUE(border)) {
                for (ht in Yc_plots_ht@ht_list) {
                  decorate_heatmap_body(ht@name, {
                    grid.rect(gp = gpar(col = "black", fill = NA, 
                      lwd = border.width))
                  })
                }
            }
            popViewport()
        }
        pushViewport(viewport(layout.pos.row = if (!is.null(col.prox) || 
            isTRUE(isProximityMatrix)) 
            layout_row - 1
        else layout_row, layout.pos.col = if (!is.null(YdNum)) 
            2
        else 1))
        for (i in seq_along(YcNum)) {
            grid.text(colnames(data)[YcNum[i]], x = unit(0.5 * 
                i - 0.05, "cm"), y = if (!is.null(col.prox) || 
                isTRUE(isProximityMatrix)) 
                unit(0, "cm")
            else unit(1, "npc"), just = "left", gp = gpar(fontsize = ifelse(is.null(Yc.label.size), 
                8, Yc.label.size)), rot = 90)
        }
        popViewport()
        for (i in seq_along(YcNum)) {
            list_Yc[[i]] <- get(paste0("Yc_encoded_", i))[if (is.list(row_order(Yc_plots_ht))) 
                row_order(Yc_plots_ht)[[1]]
            else row_order(Yc_plots_ht)]
        }
    }
    if (!is.null(YdNum)) {
        if (length(Yd_plots) == 1) {
            pushViewport(viewport(layout.pos.row = layout_row, 
                layout.pos.col = 1))
            Yd_plots_ht <- draw(Yd_plots[[1]], newpage = FALSE)
            if (isTRUE(border)) {
                add_border(Yd_plots_ht)
            }
            popViewport()
        }
        else {
            pushViewport(viewport(layout.pos.row = layout_row, 
                layout.pos.col = 1))
            Yd_plots_ht <- draw(Reduce(`+`, Yd_plots), newpage = FALSE, 
                gap = unit(0, "cm"))
            if (isTRUE(border)) {
                for (ht in Yd_plots_ht@ht_list) {
                  decorate_heatmap_body(ht@name, {
                    grid.rect(gp = gpar(col = "black", fill = NA, 
                      lwd = border.width))
                  })
                }
            }
            popViewport()
        }
        pushViewport(viewport(layout.pos.row = if (!is.null(col.prox) || 
            isTRUE(isProximityMatrix)) 
            layout_row - 1
        else layout_row, layout.pos.col = 1))
        for (i in seq_along(YdNum)) {
            grid.text(colnames(data)[YdNum[i]], x = unit(0.5 * 
                i - 0.05, "cm"), y = if (!is.null(col.prox) || 
                isTRUE(isProximityMatrix)) 
                unit(0, "cm")
            else unit(1, "npc"), just = "left", gp = gpar(fontsize = ifelse(is.null(Yd.label.size), 
                8, Yd.label.size)), rot = 90)
        }
        popViewport()
        for (i in seq_along(YdNum)) {
            list_Yd[[i]] <- list_Yd[[i]][row_order(Yd_plots_ht), 
                ]
        }
    }
    if (!is.null(XcNum)) {
        list_Xc <- list()
        if (length(Xc_plots) == 1) {
            pushViewport(viewport(layout.pos.row = if (!is.null(col.prox) && 
                isTRUE(col.order %in% c("single", "complete", 
                  "average")) && isFALSE(show.col.prox)) {
                3
            }
            else if (!is.null(col.prox) || isTRUE(isProximityMatrix)) {
                2
            }
            else {
                1
            }, layout.pos.col = if (is.null(row.prox) && isFALSE(isProximityMatrix)) 
                layout_col
            else layout_col - 1))
            Xc_plots_ht <- draw(Xc_plots[[1]], newpage = FALSE)
            if (isTRUE(border)) {
                add_border(Xc_plots_ht)
            }
            popViewport()
        }
        else {
            pushViewport(viewport(layout.pos.row = if (!is.null(col.prox) && 
                isTRUE(col.order %in% c("single", "complete", 
                  "average")) && isFALSE(show.col.prox)) {
                3
            }
            else if (!is.null(col.prox) || isTRUE(isProximityMatrix)) {
                2
            }
            else {
                1
            }, layout.pos.col = if (is.null(row.prox) && isFALSE(isProximityMatrix)) 
                layout_col
            else layout_col - 1))
            Xc_plots_ht <- draw(Reduce(`%v%`, Xc_plots), newpage = FALSE, 
                gap = unit(0, "cm"))
            if (isTRUE(border)) {
                for (ht in Xc_plots_ht@ht_list) {
                  decorate_heatmap_body(ht@name, {
                    grid.rect(gp = gpar(col = "black", fill = NA, 
                      lwd = border.width))
                  })
                }
            }
            popViewport()
        }
        pushViewport(viewport(layout.pos.row = if (!is.null(col.prox) && 
            isTRUE(col.order %in% c("single", "complete", "average")) && 
            isFALSE(show.col.prox)) {
            3
        }
        else if (!is.null(col.prox) || isTRUE(isProximityMatrix)) {
            2
        }
        else {
            1
        }, layout.pos.col = layout_col))
        for (i in seq_along(XcNum)) {
            grid.text(Xc.name[i], x = if (!is.null(row.prox) || 
                isTRUE(isProximityMatrix)) 
                unit(0, "cm")
            else unit(1, "npc"), y = unit(0.5 * i + 0.05, "cm"), 
                just = "left", gp = gpar(fontsize = ifelse(is.null(Xc.label.size), 
                  8, Xc.label.size)))
        }
        popViewport()
        for (i in seq_along(XcNum)) {
            list_Xc[[i]] <- get(paste0("Xc_encoded_", i))[column_order(Xc_plots_ht)]
        }
    }
    if (!is.null(XdNum)) {
        if (length(Xd_plots) == 1) {
            pushViewport(viewport(layout.pos.row = layout_row - 
                1, layout.pos.col = if (is.null(row.prox) && 
                isFALSE(isProximityMatrix)) 
                layout_col
            else layout_col - 1))
            Xd_plots_ht <- draw(Xd_plots[[1]], newpage = FALSE)
            if (isTRUE(border)) {
                add_border(Xd_plots_ht)
            }
            popViewport()
        }
        else {
            pushViewport(viewport(layout.pos.row = layout_row - 
                1, layout.pos.col = if (is.null(row.prox) && 
                isFALSE(isProximityMatrix)) 
                layout_col
            else layout_col - 1))
            Xd_plots_ht <- draw(Reduce(`%v%`, Xd_plots), newpage = FALSE, 
                gap = unit(0, "cm"))
            if (isTRUE(border)) {
                for (ht in Xd_plots_ht@ht_list) {
                  decorate_heatmap_body(ht@name, {
                    grid.rect(gp = gpar(col = "black", fill = NA, 
                      lwd = border.width))
                  })
                }
            }
            popViewport()
        }
        pushViewport(viewport(layout.pos.row = layout_row - 1, 
            layout.pos.col = layout_col))
        for (i in seq_along(XdNum)) {
            grid.text(Xd.name[i], x = if (!is.null(row.prox) || 
                isTRUE(isProximityMatrix)) 
                unit(0, "cm")
            else unit(1, "npc"), y = unit(0.5 * i + 0.05, "cm"), 
                just = "left", gp = gpar(fontsize = ifelse(is.null(Xd.label.size), 
                  8, Xd.label.size)))
        }
        popViewport()
        for (i in seq_along(XdNum)) {
            list_Xd[[i]] <- list_Xd[[i]][column_order(Xd_plots_ht), 
                ]
        }
    }
    if ((!is.null(XdNum) || !is.null(YdNum)) && (length(XcNum) > 
        1 || length(YcNum) > 1)) {
        group2_list <- list()
        if (!is.null(XdNum)) 
            group2_list <- c(group2_list, Xd_colorbar)
        if (!is.null(YdNum)) 
            group2_list <- c(group2_list, Yd_colorbar)
        col_group2 <- group2_list
        col_group1 <- setdiff(colorbar_list, col_group2)
        packed_col1 <- do.call(packLegend, c(col_group1, list(direction = "vertical", 
            gap = unit(0.5, "cm"))))
        packed_col2 <- do.call(packLegend, c(col_group2, list(direction = "vertical", 
            gap = unit(0.5, "cm"))))
        pushViewport(viewport(layout.pos.col = layout_col, ))
        draw(packed_col1, x = unit(1, "npc") + unit(colorbar.margin, 
            "cm"), y = unit(1, "npc"), just = c("left", "top"))
        draw(packed_col2, x = unit(1, "npc") + unit(colorbar.margin + 
            5, "cm"), y = unit(1, "npc"), just = c("left", "top"))
        popViewport()
    }
    else {
        packed_colorbar <- do.call(packLegend, c(colorbar_list, 
            list(direction = "vertical", gap = unit(0.5, "cm"))))
        pushViewport(viewport(layout.pos.col = layout_col, ))
        draw(packed_colorbar, x = unit(1, "npc") + unit(colorbar.margin, 
            "cm"), y = unit(1, "npc"), just = c("left", "top"))
        popViewport()
    }
    output_plot <- grid.grab()
    dev.off()
    if (isTRUE(show.plot)) {
        grid.newpage()
        grid.draw(output_plot)
    }
    result <- list()
    if ((row.order %in% c("single", "complete", "average"))) {
        if (isFALSE(isProximityMatrix)) {
            order_row <- row_order(row_prox_ht)
        }
        else {
            order_row <- row_order(original_plot_ht)
        }
    }
    else if (row.order == "original") {
        order_row <- 1:length(row_names)
    }
    else {
        order_row <- order_row
    }
    if (isTRUE(exp.row_order)) {
        result$row_order <- order_row
    }
    if (isFALSE(isProximityMatrix)) {
        if ((col.order %in% c("single", "complete", "average"))) {
            order_col <- column_order(col_prox_ht)
        }
        else if (col.order == "original") {
            order_col <- 1:length(column_names)
        }
        else {
            order_col <- order_col
        }
    }
    if (isTRUE(exp.column_order) && isFALSE(isProximityMatrix)) {
        result$column_order <- order_col
    }
    if (isTRUE(exp.row_names)) {
        result$row_names <- row_names
    }
    if (isTRUE(exp.column_names) && isFALSE(isProximityMatrix)) {
        result$column_names <- column_names
    }
    if (isTRUE(exp.Xc) && isFALSE(isProximityMatrix)) {
        result$Xc <- list_Xc
    }
    if (isTRUE(exp.Yc) && isFALSE(isProximityMatrix)) {
        result$Yc <- list_Yc
    }
    if (isTRUE(exp.Xd) && isFALSE(isProximityMatrix)) {
        result$Xd <- list_Xd
    }
    if (isTRUE(exp.Yd) && isFALSE(isProximityMatrix)) {
        result$Yd <- list_Yd
    }
    if (isTRUE(exp.Xd_codebook) && isFALSE(isProximityMatrix)) {
        result$Xd_codebook <- list_cb_Xd
    }
    if (isTRUE(exp.Yd_codebook) && isFALSE(isProximityMatrix)) {
        result$Yd_codebook <- list_cb_Yd
    }
    original_plot_data <- data_matrix[row_order(original_plot_ht), 
        if (isFALSE(isProximityMatrix)) 
            column_order(original_plot_ht)
        else row_order(original_plot_ht)]
    if (isTRUE(exp.originalmatrix)) {
        result$originalmatrix <- original_plot_data
    }
    if (isFALSE(isProximityMatrix)) {
        row_proximity_matrix <- distance_matrix_for_row[row_order(row_prox_ht), 
            row_order(row_prox_ht)]
    }
    if (isTRUE(exp.row_prox) && isFALSE(isProximityMatrix)) {
        result$row_prox <- row_proximity_matrix
    }
    if (isFALSE(isProximityMatrix)) {
        col_proximity_matrix <- distance_matrix_for_column[column_order(col_prox_ht), 
            column_order(col_prox_ht)]
    }
    if (isTRUE(exp.col_prox) && isFALSE(isProximityMatrix)) {
        result$col_prox <- col_proximity_matrix
    }
    if (length(result) > 0) {
        return(result)
    }
}
<bytecode: 0x561a3e223700>
<environment: namespace:GAPR>
