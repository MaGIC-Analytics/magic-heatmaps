# ─── Active Gene List Reactive ─────────────────────────────────────────────────

ActiveGenes <- reactiveVal(NULL)

observeEvent(input$apply_genes, {
    mat <- ProcessedMatrix()
    req(mat)
    all_genes <- rownames(mat)

    if (input$gene_source == "all") {
        ActiveGenes(all_genes)

    } else if (input$gene_source == "custom") {
        raw <- input$custom_genes
        shiny::validate(need(nchar(trimws(raw)) > 0, "Please enter at least one gene name."))
        # Support both newline and comma separation
        genes <- trimws(unlist(strsplit(raw, "[,\n]+")))
        genes <- genes[genes != ""]
        found <- genes[genes %in% all_genes]
        missing <- genes[!genes %in% all_genes]
        if (length(missing) > 0) {
            showNotification(
                paste("Genes not found in matrix:", paste(missing, collapse=", ")),
                type='warning', duration=8
            )
        }
        shiny::validate(need(length(found) > 0, "None of the entered genes were found in the matrix."))
        ActiveGenes(found)

    } else if (input$gene_source == "database") {
        req(input$db_geneset)
        tryCatch({
            species_sel <- input$db_species %||% "Homo sapiens"
            src         <- input$db_source  %||% "H"

            if (src == "H") {
                gs_df <- msigdbr(species=species_sel, category="H")
            } else if (src == "C2_KEGG") {
                gs_df <- msigdbr(species=species_sel, category="C2", subcategory="CP:KEGG")
            } else if (src == "C5_BP") {
                gs_df <- msigdbr(species=species_sel, category="C5", subcategory="GO:BP")
            } else if (src == "C5_MF") {
                gs_df <- msigdbr(species=species_sel, category="C5", subcategory="GO:MF")
            }

            gs_genes <- gs_df$gene_symbol[gs_df$gs_name == input$db_geneset]
            found    <- gs_genes[gs_genes %in% all_genes]
            showNotification(
                paste0("Gene set '", input$db_geneset, "': ", length(gs_genes),
                       " genes in set, ", length(found), " found in matrix."),
                type='message', duration=6
            )
            shiny::validate(need(length(found) > 0, "No genes from the selected gene set were found in the matrix."))
            ActiveGenes(found)
        }, error=function(e) {
            showNotification(paste("Database error:", e$message), type='error', duration=NULL)
        })
    }
}, ignoreNULL=FALSE)

# Initialize with all genes when data first loads (so heatmap renders immediately)
observeEvent(list(input$submit, input$demo_submit), {
    mat <- ProcessedMatrix()
    req(mat)
    ActiveGenes(rownames(mat))
    showNotification("Data loaded. Generating heatmap...", type='message', duration=3)
}, ignoreNULL=TRUE)

# ─── Color Function Builder ─────────────────────────────────────────────────────

build_color_fun <- function(mat_vals, color_mode_two, color_high, color_low,
                            color_palette, color_reverse) {
    rng <- range(mat_vals, na.rm=TRUE)
    mid <- 0  # midpoint for diverging (z-score scaled data centers on 0)

    if (color_mode_two) {
        colorRamp2(c(rng[1], mid, rng[2]),
                   c(color_low, "white", color_high))
    } else {
        pal_colors <- tryCatch(
            brewer.pal(11, color_palette),
            error=function(e) {
                if (color_palette == "viridis") {
                    c("#440154","#482878","#3E4989","#31688E","#26828E",
                      "#1F9E89","#35B779","#6DCD59","#B4DE2C","#FDE725")
                } else if (color_palette == "plasma") {
                    c("#0D0887","#41049D","#6A00A8","#8F0DA4","#B12A90",
                      "#CC4678","#E16462","#F2844B","#FCA636","#FCCE25","#F0F921")
                } else {
                    colorRampPalette(c("#4575b4","white","#d73027"))(11)
                }
            }
        )
        if (color_reverse) pal_colors <- rev(pal_colors)
        colorRamp2(seq(rng[1], rng[2], length.out=length(pal_colors)), pal_colors)
    }
}

# ─── Annotation Builder ────────────────────────────────────────────────────────

build_column_annotation <- function(meta_df, anno_cols, anno_bar_size,
                                    anno_font_size, anno_position, input,
                                    legend_title_size, legend_font_size,
                                    legend_grid_height, legend_grid_width) {
    if (is.null(anno_cols) || length(anno_cols) == 0) return(NULL)

    anno_list   <- list()
    color_list  <- list()

    for (col in anno_cols) {
        vals <- as.character(meta_df[[col]])
        uniq <- unique(vals[!is.na(vals)])

        # Build color mapping for this column from colourInput widgets
        col_colors <- sapply(uniq, function(v) {
            id <- paste0("anno_color_", col, "_", gsub("[^A-Za-z0-9]", "_", v))
            clr <- input[[id]]
            if (is.null(clr)) "#999999" else clr
        })
        names(col_colors) <- uniq

        anno_list[[col]]  <- vals
        color_list[[col]] <- col_colors
    }

    # Mirror the same legend params across all annotation columns
    anno_legend_params <- lapply(anno_cols, function(col) {
        list(
            title_gp    = gpar(fontsize = legend_title_size, fontface = "bold"),
            labels_gp   = gpar(fontsize = legend_font_size),
            grid_height = unit(legend_grid_height, "mm"),
            grid_width  = unit(legend_grid_width,  "mm")
        )
    })
    names(anno_legend_params) <- anno_cols

    HeatmapAnnotation(
        df                    = as.data.frame(anno_list),
        col                   = color_list,
        annotation_height     = unit(anno_bar_size, "mm"),
        annotation_name_gp    = gpar(fontsize = anno_font_size),
        annotation_legend_param = anno_legend_params,
        which                 = "column"
    )
}

# ─── Main Heatmap Reactive ─────────────────────────────────────────────────────

HeatmapPlotter <- reactive({
    mat  <- ProcessedMatrix()
    meta <- ProcessedMeta()
    req(mat, meta)

    genes <- ActiveGenes()
    if (is.null(genes) || length(genes) == 0) genes <- rownames(mat)

    # Subset to active genes (keep only those present)
    genes <- genes[genes %in% rownames(mat)]
    shiny::validate(need(length(genes) >= 2, "At least 2 genes are required to draw a heatmap."))
    mat_sub <- mat[genes, , drop=FALSE]

    # Row scaling (z-score)
    if (isTRUE(input$scale_rows)) {
        mat_sub <- t(scale(t(mat_sub)))
        mat_sub[is.nan(mat_sub)] <- 0
    }

    # Color function
    col_fun <- build_color_fun(
        mat_vals       = mat_sub,
        color_mode_two = identical(input$color_mode, "two"),
        color_high     = input$color_high    %||% "#d73027",
        color_low      = input$color_low     %||% "#4575b4",
        color_palette  = input$color_palette %||% "RdBu",
        color_reverse  = isTRUE(input$color_reverse)
    )

    # Column annotation
    anno_cols <- input$anno_cols
    top_anno  <- bottom_anno <- NULL
    if (!is.null(anno_cols) && length(anno_cols) > 0) {
        col_anno <- build_column_annotation(
            meta_df            = meta,
            anno_cols          = anno_cols,
            anno_bar_size      = input$anno_bar_size      %||% 5,
            anno_font_size     = input$anno_font_size     %||% 10,
            anno_position      = input$anno_position      %||% "top",
            input              = input,
            legend_title_size  = input$legend_title_size  %||% 12,
            legend_font_size   = input$legend_font_size   %||% 10,
            legend_grid_height = input$legend_grid_height %||% 4,
            legend_grid_width  = input$legend_grid_width  %||% 4
        )
        if (identical(input$anno_position, "top")) {
            top_anno <- col_anno
        } else {
            bottom_anno <- col_anno
        }
    }

    # Row / column splitting
    row_split <- col_split <- NULL
    if (!is.null(input$split_mode)) {
        if (input$split_mode == "kmeans") {
            row_split <- input$row_km %||% 2
            col_split <- input$col_km %||% 1
            if (col_split == 1) col_split <- NULL
        } else if (input$split_mode == "metadata") {
            req(input$split_col)
            col_split <- factor(as.character(meta[[input$split_col]]))
        }
    }

    # Legend params
    heatmap_legend_param <- list(
        title             = if (isTRUE(input$scale_rows)) "z-score" else "Expression",
        title_gp          = gpar(fontsize = input$legend_title_size %||% 12, fontface="bold"),
        labels_gp         = gpar(fontsize = input$legend_font_size  %||% 10),
        legend_height     = unit(input$legend_bar_height %||% 40, "mm"),
        grid_height       = unit(input$legend_grid_height %||% 4, "mm"),
        grid_width        = unit(input$legend_grid_width  %||% 4, "mm")
    )

    # Build heatmap
    Heatmap(
        mat_sub,
        col                 = col_fun,
        top_annotation      = top_anno,
        bottom_annotation   = bottom_anno,

        # Row settings
        show_row_names      = isTRUE(input$show_row_names),
        row_names_gp        = gpar(fontsize = input$row_font_size %||% 8),
        row_names_rot       = input$row_font_angle %||% 0,
        cluster_rows        = isTRUE(input$cluster_rows),
        clustering_distance_rows = "pearson",
        row_dend_reorder    = isTRUE(input$reorder_rows),
        row_dend_width      = unit(input$row_dend_size %||% 10, "mm"),
        row_split           = row_split,

        # Column settings
        show_column_names   = isTRUE(input$show_col_names),
        column_names_gp     = gpar(fontsize = input$col_font_size %||% 8),
        column_names_rot    = input$col_font_angle %||% 45,
        cluster_columns     = isTRUE(input$cluster_cols),
        clustering_distance_columns = "pearson",
        column_dend_reorder = isTRUE(input$reorder_cols),
        column_dend_height  = unit(input$col_dend_size %||% 10, "mm"),
        column_split        = col_split,

        # Legend
        heatmap_legend_param = heatmap_legend_param,

        # Aesthetics
        border              = TRUE,
        rect_gp             = gpar(col="white", lwd=0.5)
    )
})

# ─── Render Heatmap Output ─────────────────────────────────────────────────────

output$heatmap_out <- renderPlot({
    ht <- HeatmapPlotter()
    draw(ht, merge_legend=TRUE)
}, height = function() input$hm_height %||% 700,
   width  = function() input$hm_width  %||% 900)

# Keep rendering even when the Heatmap tab is hidden so the plot is ready
# the moment the tab is revealed (avoids the blank-until-tab-switch bug).
outputOptions(output, "heatmap_out", suspendWhenHidden=FALSE)

# ─── Reproducible Code Generator ─────────────────────────────────────────────
# Plain function — no reactive wrapper, called directly from the observer.
# Reads ProcessedMatrix/Meta via isolate() to pull current cached values
# without creating reactive dependencies or triggering req() conditions.

build_code_string <- function(mat, meta, genes, inp) {

    genes <- genes[genes %in% rownames(mat)]

    # Gene list block
    if (length(genes) == nrow(mat)) {
        gene_block <- "genes <- rownames(mat)  # all genes\n"
    } else if (length(genes) <= 30) {
        gene_block <- paste0('genes <- c(\n  "', paste(genes, collapse = '",\n  "'), '"\n)\n')
    } else {
        gene_block <- paste0(
            sprintf("# %d genes selected — first 5 shown; paste your full list here\n", length(genes)),
            'genes <- c(\n  "', paste(head(genes, 5), collapse = '",\n  "'), '",\n  ...\n)\n'
        )
    }

    # Color block
    color_mode <- inp$color_mode %||% "two"
    if (color_mode == "two") {
        color_block <- sprintf(
            '# Two-color balance\ncol_fun <- colorRamp2(c(-2, 0, 2), c("%s", "white", "%s"))\n',
            inp$color_low %||% "#4575b4", inp$color_high %||% "#d73027")
    } else {
        pal      <- inp$color_palette %||% "RdBu"
        rev_line <- if (isTRUE(inp$color_reverse)) "pal_colors <- rev(pal_colors)\n" else ""
        color_block <- sprintf(
            'pal_colors <- RColorBrewer::brewer.pal(11, "%s")\n%scol_fun <- colorRamp2(seq(-2, 2, length.out = 11), pal_colors)\n',
            pal, rev_line)
    }

    # Scale block
    scale_block <- if (isTRUE(inp$scale_rows))
        "mat_sub <- t(scale(t(mat_sub)))\nmat_sub[is.nan(mat_sub)] <- 0\n" else ""

    # Annotation block
    anno_cols <- inp$anno_cols
    if (!is.null(anno_cols) && length(anno_cols) > 0) {
        col_defs <- lapply(anno_cols, function(col) {
            vals   <- unique(as.character(meta[[col]]))
            vals   <- vals[!is.na(vals)]
            colors <- sapply(vals, function(v) {
                id  <- paste0("anno_color_", col, "_", gsub("[^A-Za-z0-9]", "_", v))
                clr <- inp[[id]]
                if (is.null(clr)) "#999999" else clr
            })
            color_str <- paste0('c(', paste0('"', vals, '" = "', colors, '"', collapse = ", "), ')')
            sprintf('    %s = anno_simple(meta$%s, col = %s)', col, col, color_str)
        })
        anno_block <- paste0("top_anno <- HeatmapAnnotation(\n", paste(col_defs, collapse = ",\n"), "\n)\n")
        anno_arg   <- "    top_annotation = top_anno,\n"
    } else {
        anno_block <- anno_arg <- ""
    }

    # Split block
    split_block <- split_args <- ""
    split_mode  <- inp$split_mode %||% "none"
    if (split_mode == "kmeans") {
        split_args <- sprintf("    row_split = %d,\n    column_split = %d,\n",
                              inp$row_km %||% 2, inp$col_km %||% 1)
    } else if (split_mode == "metadata" && !is.null(inp$split_col)) {
        split_block <- sprintf('col_split <- factor(meta$%s)\n', inp$split_col)
        split_args  <- "    column_split = col_split,\n"
    }

    paste0(
        "library(ComplexHeatmap)\nlibrary(circlize)\nlibrary(RColorBrewer)\n\n",
        "# ── Gene selection ──────────────────────────────────────────────────────────\n",
        gene_block,
        "mat_sub <- mat[genes, , drop = FALSE]\n\n",
        if (nchar(scale_block) > 0) paste0("# ── Row scaling (z-score) ──────────────────────────────────────────────────\n", scale_block, "\n") else "",
        "# ── Color scale ─────────────────────────────────────────────────────────────\n",
        color_block, "\n",
        if (nchar(anno_block) > 0) paste0("# ── Column annotation ──────────────────────────────────────────────────────\n", anno_block, "\n") else "",
        if (nchar(split_block) > 0) paste0("# ── Column splitting ───────────────────────────────────────────────────────\n", split_block, "\n") else "",
        "# ── Heatmap ─────────────────────────────────────────────────────────────────\n",
        "ht <- Heatmap(\n",
        "    mat_sub,\n",
        "    col                         = col_fun,\n",
        anno_arg, split_args,
        sprintf("    cluster_rows                = %s,\n",  tolower(as.character(isTRUE(inp$cluster_rows)))),
        sprintf("    cluster_columns             = %s,\n",  tolower(as.character(isTRUE(inp$cluster_cols)))),
        "    clustering_distance_rows    = \"pearson\",\n",
        "    clustering_distance_columns = \"pearson\",\n",
        sprintf("    show_row_names              = %s,\n",  tolower(as.character(isTRUE(inp$show_row_names)))),
        sprintf("    show_column_names           = %s,\n",  tolower(as.character(isTRUE(inp$show_col_names)))),
        sprintf("    row_names_rot               = %s,\n",  as.character(inp$row_font_angle  %||% 0)),
        sprintf("    column_names_rot            = %s,\n",  as.character(inp$col_font_angle  %||% 45)),
        sprintf("    row_names_gp                = gpar(fontsize = %s),\n", as.character(inp$row_font_size %||% 8)),
        sprintf("    column_names_gp             = gpar(fontsize = %s),\n", as.character(inp$col_font_size %||% 8)),
        "    border                      = TRUE,\n",
        "    rect_gp                     = gpar(col = \"white\", lwd = 0.5),\n",
        "    heatmap_legend_param        = list(\n",
        sprintf("        title     = \"%s\",\n", if (isTRUE(inp$scale_rows)) "z-score" else "Expression"),
        sprintf("        title_gp  = gpar(fontsize = %s, fontface = \"bold\"),\n", as.character(inp$legend_title_size %||% 12)),
        sprintf("        labels_gp = gpar(fontsize = %s)\n",  as.character(inp$legend_font_size  %||% 10)),
        "    )\n)\n\n",
        "draw(ht, merge_legend = TRUE)\n"
    )
}

# ─── Code Modal ───────────────────────────────────────────────────────────────

observeEvent(input$show_code_modal, {
    # Read current data via isolate — avoids re-triggering reactive graph
    mat   <- isolate(tryCatch(ProcessedMatrix(), error = function(e) NULL))
    meta  <- isolate(tryCatch(ProcessedMeta(),   error = function(e) NULL))
    genes <- isolate(ActiveGenes())

    if (is.null(mat) || is.null(meta)) {
        code <- "# Generate a heatmap first, then click here to get the reproducible code."
    } else {
        if (is.null(genes) || length(genes) == 0) genes <- rownames(mat)
        code <- tryCatch(
            build_code_string(mat, meta, genes, input),
            error = function(e) paste0("# Code generation error: ", conditionMessage(e))
        )
    }

    showModal(modalDialog(
        title     = tagList(icon("file-code"), " Reproducible R Code"),
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        p("Copy this code to reproduce your current heatmap in an offline R session.",
          style = "color:#555; margin-bottom:12px;"),
        tags$pre(
            style = paste(
                "background:#1e1e1e; color:#d4d4d4; border-radius:6px;",
                "padding:16px; font-size:12px; max-height:520px; overflow-y:auto;",
                "white-space:pre; font-family:'Courier New', monospace;"
            ),
            code
        )
    ))
})

# ─── Help Modal ───────────────────────────────────────────────────────────────

show_help_modal_ui <- function() {
    showModal(modalDialog(
        title     = tagList(icon("circle-question"), " Heatmap Tool Help"),
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        tabsetPanel(
            tabPanel("Overview",
                br(),
                h4("What is a Heatmap?"),
                p("A heatmap encodes a numeric matrix as a colour grid. Each cell represents one gene (row) in one sample (column). Colour intensity reflects the expression level, making it easy to spot patterns across genes and samples simultaneously."),
                h4("Reading the Colour Scale"),
                tags$ul(
                    tags$li(strong("Two-colour balance:"), " a diverging scale anchored at white (midpoint). Blue → low expression, Red → high expression (or whichever colours you pick)."),
                    tags$li(strong("Colour palette:"), " a sequential or diverging palette applied across the full data range.")
                ),
                h4("Z-score Scaling"),
                p("When 'Scale rows (z-score)' is ON, each gene is centred and scaled independently: z = (x − mean) / sd. This removes between-gene magnitude differences and highlights relative variation across samples. The legend then shows z-score units, not raw expression."),
                h4("Workflow"),
                tags$ol(
                    tags$li("Upload your matrix (genes × samples) and metadata, then click Submit."),
                    tags$li("Choose which genes to display using the Gene Selection dropdown."),
                    tags$li("Toggle options in the sidebar — the heatmap updates immediately."),
                    tags$li("Click ", icon("file-code"), " to get reproducible R code for your current view."),
                    tags$li("Expand Resize, choose a format, and download.")
                )
            ),
            tabPanel("Gene Selection",
                br(),
                h4("All Genes"), p("Displays every gene present in your uploaded matrix."),
                h4("Custom Gene List"),
                p("Paste gene symbols, one per line or comma-separated. Genes not found in the matrix are ignored with a warning notification."),
                h4("Database Gene Set"),
                tags$ul(
                    tags$li(strong("MSigDB Hallmark (H):"), " 50 well-defined biological states and processes."),
                    tags$li(strong("MSigDB C2 KEGG:"), " Pathway gene sets from the KEGG database."),
                    tags$li(strong("MSigDB C5 GO BP:"), " Gene Ontology Biological Process terms."),
                    tags$li(strong("MSigDB C5 GO MF:"), " Gene Ontology Molecular Function terms.")
                ),
                p("Type in the ", strong("Gene Set"), " field to search. Only genes present in your matrix are used; a notification shows how many matched.")
            ),
            tabPanel("Controls",
                br(),
                h4("Color Options"),
                tags$ul(
                    tags$li(strong("Two-color balance:"), " Pick high and low colours. White anchors the midpoint."),
                    tags$li(strong("Color palette:"), " Select a named RColorBrewer palette and optionally reverse it."),
                    tags$li(strong("Scale rows:"), " Toggle z-score normalisation per gene.")
                ),
                h4("Annotation"),
                p("Select one or more metadata columns to draw as colour bars above (or below) the heatmap. Expand to customise bar thickness, label size, and the colour assigned to each group level."),
                h4("Fonts & Labels"),
                p("Toggle row (gene) and column (sample) name visibility independently. Adjust font size and rotation angle for each axis."),
                h4("Clustering"),
                p("Enable hierarchical clustering for rows and/or columns using Pearson distance. Adjust the dendrogram arm width and toggle reordering."),
                h4("Splitting"),
                tags$ul(
                    tags$li(strong("K-means:"), " Split rows and/or columns into k clusters automatically."),
                    tags$li(strong("Metadata Column:"), " Split columns by a categorical metadata variable (e.g. Group).")
                ),
                h4("Legend"),
                p("Sliders control font sizes for the legend title and labels, the bar height (continuous legend), and the grid cell dimensions (discrete annotation legends). All legend types are linked to the same set of sliders."),
                h4("Resize"),
                p("Drag the height and width sliders to scale the rendered plot canvas. The same dimensions are used when downloading.")
            )
        )
    ))
}

# Both the Heatmap-tab button and the global floating button share the same modal
observeEvent(input$show_help_modal,  { show_help_modal_ui() })
observeEvent(input$show_help_float,  { show_help_modal_ui() })

# ─── Active Gene List Table ────────────────────────────────────────────────────

output$gene_list_table <- DT::renderDataTable({
    genes <- ActiveGenes()
    req(genes)
    df <- data.frame(Gene=genes, stringsAsFactors=FALSE)
    DT::datatable(df, style='bootstrap', options=list(pageLength=20, scrollX=TRUE))
})

# ─── Download Handler ─────────────────────────────────────────────────────────

output$download_heatmap <- downloadHandler(
    filename = function() {
        paste0("heatmap.", input$hm_download_format)
    },
    content = function(file) {
        ht   <- HeatmapPlotter()
        h_in <- (input$hm_height %||% 700) / 96
        w_in <- (input$hm_width  %||% 900) / 96
        fmt  <- input$hm_download_format

        if (fmt == "png") {
            png(file, height=input$hm_height, width=input$hm_width, res=96)
        } else if (fmt == "jpeg") {
            jpeg(file, height=input$hm_height, width=input$hm_width, res=96)
        } else if (fmt == "tiff") {
            tiff(file, height=input$hm_height, width=input$hm_width, res=96)
        } else if (fmt == "pdf") {
            pdf(file, height=h_in, width=w_in)
        } else if (fmt == "svg") {
            svg(file, height=h_in, width=w_in)
        } else if (fmt == "eps") {
            setEPS()
            postscript(file, height=h_in, width=w_in)
        }
        draw(ht, merge_legend=TRUE)
        dev.off()
    }
)
