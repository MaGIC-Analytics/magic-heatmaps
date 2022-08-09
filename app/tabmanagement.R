# Hiding and showing tabs on command!
############################################################################
observe({
    if(is.null(input$counts_file) || is.null(input$metadata_file)){
        hideTab(inputId = "NAVTABS", target = "Heatmaps")
    }
})

observeEvent(input$submit,{
    if( (!is.null(input$counts_file)) & (!is.null(input$metadata_file))){
        showTab(inputId = "NAVTABS", target = "Heatmaps")
        updateTabsetPanel(session, inputId="NAVTABS", selected="Heatmaps")
    }      
})

observeEvent(input$demo_submit,{
    showTab(inputId = "NAVTABS", target = "Heatmaps")
    updateTabsetPanel(session, inputId="NAVTABS", selected="Heatmaps")
})
