# ─── Tab Visibility Management ─────────────────────────────────────────────────
# The "Heatmap" tab is hidden by default and only shown after data is submitted.

# Hide on initial load / when no data is present
observe({
    hideTab(inputId="NAVTABS", target="Heatmap")
})

# Show tab after custom data is submitted
observeEvent(input$submit, {
    mat  <- MatrixReactive()
    meta <- MetadataReactive()
    if (!is.null(mat) && !is.null(meta) && nrow(mat) > 0 && nrow(meta) > 0) {
        showTab(inputId="NAVTABS", target="Heatmap")
        updateTabsetPanel(session, inputId="NAVTABS", selected="Heatmap")
        shinyjs::delay(300, shinyjs::runjs("$(window).trigger('resize');"))
    }
})

# Show tab after demo data is loaded
observeEvent(input$demo_submit, {
    showTab(inputId="NAVTABS", target="Heatmap")
    updateTabsetPanel(session, inputId="NAVTABS", selected="Heatmap")
    shinyjs::delay(300, shinyjs::runjs("$(window).trigger('resize');"))
})
