# ─── Data Loading Reactives ────────────────────────────────────────────────────

# Determine which file parser to use based on extension
read_delim_auto <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("tsv", "txt")) {
        fread(path, sep="\t")
    } else {
        fread(path, sep=",")
    }
}

# ── msigdbr version-compatibility shim ─────────────────────────────────────────
# msigdbr >= 10 renamed category/subcategory -> collection/subcollection and split
# KEGG into CP:KEGG_LEGACY/CP:KEGG_MEDICUS. Take the legacy argument names, try the
# new API first, fall back to the old one. gs_name / gene_symbol are stable.
msigdbr_compat <- function(species, category, subcategory = NULL) {
    sub_new <- if (identical(subcategory, "CP:KEGG")) "CP:KEGG_LEGACY" else subcategory
    tryCatch(
        do.call(msigdbr, c(list(species = species, collection = category),
                           if (!is.null(sub_new)) list(subcollection = sub_new))),
        error = function(e)
            do.call(msigdbr, c(list(species = species, category = category),
                               if (!is.null(subcategory)) list(subcategory = subcategory)))
    )
}

# Reactive: expression matrix (rows=genes, cols=samples; first col = gene ID)
MatrixReactive <- reactive({
    if (input$DemoData == FALSE) {
        tryCatch(
            fread('www/demo_matrix.csv'),
            error = function(e) { showNotification(paste("Matrix load error:", e), type='error', duration=NULL); NULL }
        )
    } else {
        shiny::validate(need(!is.null(input$matrix_file), "Please upload an expression matrix file."))
        tryCatch(
            read_delim_auto(input$matrix_file$datapath),
            error = function(e) { showNotification(paste("Matrix parse error:", e), type='error', duration=NULL); NULL }
        )
    }
})

# Reactive: sample metadata (rows=samples, first col = sample name)
MetadataReactive <- reactive({
    if (input$DemoData == FALSE) {
        tryCatch(
            fread('www/demo_metadata.csv'),
            error = function(e) { showNotification(paste("Metadata load error:", e), type='error', duration=NULL); NULL }
        )
    } else {
        shiny::validate(need(!is.null(input$metadata_file), "Please upload a metadata file."))
        tryCatch(
            read_delim_auto(input$metadata_file$datapath),
            error = function(e) { showNotification(paste("Metadata parse error:", e), type='error', duration=NULL); NULL }
        )
    }
})

# ─── Gene Column Selector (for custom uploads) ─────────────────────────────────
output$gene_col_selector <- renderUI({
    req(input$matrix_file)
    mat <- MatrixReactive()
    req(mat)
    tagList(
        selectInput("gene_col", "Which column contains gene names?",
            choices = colnames(mat), selected = colnames(mat)[1]),
        hr()
    )
})

# ─── Preview Tables ────────────────────────────────────────────────────────────
output$matrix_table <- DT::renderDataTable({
    mat <- MatrixReactive()
    req(mat)
    DT::datatable(mat, style='bootstrap', options=list(pageLength=15, scrollX=TRUE))
})

output$metadata_table <- DT::renderDataTable({
    meta <- MetadataReactive()
    req(meta)
    DT::datatable(meta, style='bootstrap', options=list(pageLength=15, scrollX=TRUE))
})

# ─── Processed matrix for heatmap (numeric matrix, rownames = genes) ──────────
ProcessedMatrix <- reactive({
    mat  <- MatrixReactive()
    req(mat)

    # Determine gene column: user-selected (custom upload) or column 1 (demo)
    gene_col <- if (!is.null(input$gene_col) && input$DemoData) input$gene_col else colnames(mat)[1]

    gene_names <- as.character(mat[[gene_col]])
    num_cols   <- setdiff(colnames(mat), gene_col)

    m <- as.matrix(mat[, ..num_cols])
    rownames(m) <- make.unique(gene_names)
    storage.mode(m) <- "numeric"

    # Restrict to samples present in BOTH the matrix and the metadata (shared set,
    # in matrix-column order). Matches the qc/deg idiom and keeps NA-metadata
    # samples out of annotations/grouping instead of showing them as an "NA" group.
    meta_raw <- MetadataReactive()
    if (!is.null(meta_raw) && ncol(meta_raw) >= 1) {
        sample_ids <- as.character(meta_raw[[colnames(meta_raw)[1]]])
        shared <- intersect(colnames(m), sample_ids)
        shiny::validate(need(length(shared) >= 1,
            "No samples are shared between the matrix columns and the metadata sample names."))
        if (length(shared) < ncol(m)) {
            showNotification(sprintf("Using %d of %d matrix samples that are present in the metadata.",
                length(shared), ncol(m)), type='warning', duration=8)
        }
        m <- m[, shared, drop=FALSE]
    }
    m
})

# ─── Aligned metadata (rows = samples, ordered to match matrix columns) ────────
ProcessedMeta <- reactive({
    meta  <- MetadataReactive()
    mat_m <- ProcessedMatrix()
    req(meta, mat_m)

    sample_col   <- colnames(meta)[1]
    sample_names <- as.character(meta[[sample_col]])

    # Match matrix column order
    col_order <- colnames(mat_m)
    idx       <- match(col_order, sample_names)

    if (any(is.na(idx))) {
        showNotification(
            paste("Warning: Some matrix sample names were not found in metadata:",
                  paste(col_order[is.na(idx)], collapse=", ")),
            type='warning', duration=8
        )
    }

    df <- as.data.frame(meta)[idx, , drop=FALSE]
    rownames(df) <- col_order
    df
})

# ─── Annotation column selector UI (rendered on Heatmap tab sidebar) ──────────
output$annotation_col_selector <- renderUI({
    meta <- MetadataReactive()
    req(meta)
    meta_cols <- colnames(meta)[-1]   # exclude sample-name column
    tagList(
        checkboxGroupInput("anno_cols", "Metadata columns to annotate:",
            choices  = meta_cols,
            selected = meta_cols[1]
        )
    )
})

# ─── Split column selector (metadata-based splitting) ─────────────────────────
output$split_col_selector <- renderUI({
    meta <- MetadataReactive()
    req(meta)
    meta_cols <- colnames(meta)[-1]
    selectInput("split_col", "Split columns by:", choices=meta_cols, selected=meta_cols[1])
})

# ─── Dynamic annotation color pickers ─────────────────────────────────────────
output$anno_color_ui <- renderUI({
    req(input$anno_cols)
    meta    <- ProcessedMeta()
    req(meta)

    ui_elements <- lapply(input$anno_cols, function(col) {
        vals <- unique(as.character(meta[[col]]))
        vals <- vals[!is.na(vals)]
        if (length(vals) <= 12) {
            color_inputs <- lapply(seq_along(vals), function(i) {
                default_colors <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                                    "#FF7F00","#FFFF33","#A65628","#F781BF",
                                    "#999999","#66C2A5","#FC8D62","#8DA0CB")
                column(4,
                    colourInput(
                        inputId  = paste0("anno_color_", col, "_", gsub("[^A-Za-z0-9]", "_", vals[i])),
                        label    = vals[i],
                        value    = default_colors[((i - 1) %% 12) + 1]
                    )
                )
            })
            tagList(
                strong(paste("Colors for:", col)),
                fluidRow(color_inputs)
            )
        } else {
            p(paste(col, "has too many levels for manual color assignment (> 12)."))
        }
    })

    do.call(tagList, ui_elements)
})

# ─── MSigDB gene set selector ─────────────────────────────────────────────────
observe({
    req(input$gene_source == "database")
    tryCatch({
        species_sel <- input$db_species %||% "Homo sapiens"
        src         <- input$db_source  %||% "H"

        if (src == "H") {
            gs_df <- msigdbr_compat(species_sel, "H")
        } else if (src == "C2_KEGG") {
            gs_df <- msigdbr_compat(species_sel, "C2", "CP:KEGG")
        } else if (src == "C5_BP") {
            gs_df <- msigdbr_compat(species_sel, "C5", "GO:BP")
        } else if (src == "C5_MF") {
            gs_df <- msigdbr_compat(species_sel, "C5", "GO:MF")
        }

        gene_sets <- sort(unique(gs_df$gs_name))
        updateSelectizeInput(session, "db_geneset", choices=gene_sets, server=TRUE,
            options=list(placeholder='Type to search gene sets...'))
    }, error=function(e) {
        showNotification(paste("Database query error:", e$message), type='error', duration=8)
    })
})

# Null-coalescing operator (not in base R < 4.4)
`%||%` <- function(a, b) {
    if (is.null(a)) return(b)
    if (length(a) == 0) return(b)
    if (length(a) == 1 && is.character(a) && !nzchar(a)) return(b)
    a
}
