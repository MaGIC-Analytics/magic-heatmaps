# Heatmap
############################################################################
observe({
    updateSelectInput(session, "HMGroup", choices=c("None",colnames(InputReactive()$meta_data)), selected='None')
})
observe({
    updateSelectInput(session, "CPal", choices=c(rownames(brewer.pal.info)), selected='RdBu')
})

heatmapper <- reactive({
    if(input$HMGroup=='None'){
        annotation_col=NULL
    } else {
        mat_col <- data.frame(InputReactive()$meta_data[[paste(input$HMGroup)]])
        names(mat_col) <- paste(input$HMGroup)
        rownames(mat_col) <- colnames(InputReactive()$count_data)
        annotation_col=mat_col
    }

    my_colors = brewer.pal(n = 11, name = input$CPal)
    my_colors = colorRampPalette(my_colors)(50)
    my_colors = rev(my_colors)

    hm_plot <- pheatmap(InputReactive()$count_data, cluster_rows=input$RClust, cluster_cols=input$CClust, 
        scale=input$Scale, annotation_col=annotation_col,, color=my_colors,
        show_rownames=input$RName, show_colnames=input$CName, fontsize=input$HXsize, angle_col=input$Hang)
    return(hm_plot)

})

observe({
    output$heatmap_out <- renderPlot({
        heatmapper()
    }, height=input$HHeight, width=input$HWidth)
})