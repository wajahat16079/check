
## Libraries in use
#I don't think all of these are needed- try to reomve some at somepoint? 
#NOTE: installation of packages may use BiocManager::install('package_name')



# BiocManager::install('reactome.db')
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(org.Hs.eg.db)
# library(ReactomePA)
library(pathview)
library(tidyverse)

library(shiny)
library("shinydashboard")
library(plotly)
library(dplyr)
library(gprofiler2)

unlink("www",recursive = TRUE)
dir.create("./www")


ui <- dashboardPage(
  dashboardHeader(title = "Simulation Tool"),
  dashboardSidebar(disable = TRUE
                   # sidebarMenu(
                   #   menuItem("Numerical vs Numerical", tabName = "g1" ),
                   #   menuItem("Categorical vs Categorical", tabName = "g2" ),
                   #   menuItem("Categorical vs Numerical ", tabName = "g3" ),
                   #   menuItem("Categorical Only ", tabName = "g4" )
                   # )
                   # downloadButton("downloadPlot", "Download")
  ),
  dashboardBody(
    tabsetPanel(type = "tabs",
                tabPanel("Settings",
                         uiOutput("Settings")),
                tabPanel("Uploads",
                         uiOutput("Uploads")),
                tabPanel("Simulation",
                         uiOutput("Simulation")),
                tabPanel("Results",
                         uiOutput("Results")),
                tabPanel("Downloads",
                         uiOutput("Downloads"))
    ),
    
    fluidPage( 
      # uiOutput("text"),
      # # plotOutput("plot1"),
      # uiOutput("plot1"),
      # uiOutput("image1")
    )
    
  )
)


server <- function(input, output, session) {
  
  ## ------------ Reactive Values Initialization -------------------  
  
  RV = reactiveValues(
    
    ## Varaiables
    pvalue_cut_off= NULL,
    qvalue_cut_off= NULL,
    ensembl_db_ver= NULL,
    
    ## upload files 
    candidate_file_name_1 = NULL,
    universe_file_name_1 = NULL,
    
    ## download files 
    CP_ego_BP_candidates_symbol_1 = NULL,
    CP_ekegg_candidates_symbol_1 = NULL,
    CP_ereact_candidates_symbol_1 = NULL,
    filtered_pw = NULL,
    
    ## Plot Objects
    Plot1filtered_pw = NULL
    
  )
  
  ## ------------ Settings Page -------------------  
  
  output$Settings = renderUI({
    UI = list()
    
    UI = append(UI, list( box(width = 2, numericInput("pvalue_cut_off", "Enter P-Value Cut-Off", value = 1, min = -100 , max = +100))))
    UI = append(UI, list( box(width = 2, numericInput("qvalue_cut_off", "Enter Q-Value Cut-Off", value = 1, min = -100 , max = +100))))
    UI = append(UI, list( box(width = 2, numericInput("ensembl_db_ver", "Enter ensembl db ver", value = 103, min = NA , max = NA))))
    
    UI
  }) 
  
  
  ## updating the reactive values based on the input values
  observe({
    req(input$pvalue_cut_off)
    req(input$qvalue_cut_off)
    req(input$ensembl_db_ver)
    
    RV$pvalue_cut_off= input$pvalue_cut_off
    RV$qvalue_cut_off= input$qvalue_cut_off
    RV$ensembl_db_ver= input$ensembl_db_ver
    
  }) 
  
  
  ## ------------ Uploads Page -------------------  
  
  output$Uploads = renderUI({
    UI = list()
    
    UI = append(UI, list( box(width = 2,fileInput("candidate_file_name_1", "Upload Candidate File", accept = ".txt"), )))
    UI = append(UI, list( box(width = 2,fileInput("universe_file_name_1", "Upload Universe File", accept = ".txt"), )))
    
    UI
  })
  
  ## file upload handling
  observe({
    req(input$candidate_file_name_1)
    file <- input$candidate_file_name_1
    RV$candidate_file_name_1 = read_delim(file = file$datapath, delim="\n", col_names=FALSE)
  }) 
  
  ## file upload handling
  observe({
    req(input$universe_file_name_1)
    file <- input$universe_file_name_1
    RV$universe_file_name_1 = read_delim(file = file$datapath, delim="\n", col_names=FALSE)
  }) 
  
  ## ------------ Simulation Page -------------------  
  
  output$Simulation = renderUI({
    req(input$candidate_file_name_1)
    req(input$universe_file_name_1)
    
    UI = list()
    UI = append(UI, list( box(width = 2, actionButton("Run", "Run Simulation", width = '100%') )))
    UI
    
  })
  
  ## functionality is in here
  observeEvent(input$Run,{
    
    ## wait message: will be on screen until simulation finishes 
    
    showModal(modalDialog(
      title = "Simulation Running",
      "Please Wait!",
      footer = NULL
    ))
    # function here
    
    print("in func")
    
    ## --------------Global variables-------------
    
    ## input files-these are the file names of what needs to analysed and can be changed, ensure the file with that exacts file name is also in the folder. 
    
    # candidate_file_name_1 <- "candidate_genes_TGFBUP.txt"
    # universe_file_name_1 <- "universe_genes.txt"
    
    ## pval cutoffs for the pathways returned in .csv files
    
    # -------- >>>>>>>>>> three value below updating from inputs on settings page
    pvalue_cut_off <- RV$pvalue_cut_off
    qvalue_cut_off <- RV$qvalue_cut_off
    ensembl_db_ver <- RV$ensembl_db_ver # where possible the version of the Esembl database should coincide with the gene annotation used, e.g. Gencode v37 matches Ensembl v103 
    
    ##additional filters for control of GOBP in plot
    GOBP_specific_padj<-0.01 #to filter more GOBP than KEGG or REACT
    GOBP_slim_factor<-0.7 #removes redundant term
    
    ##control plots
    EMAP_padj<-0.05 #filter for significant terms to go into emap plot
    cex_category_size <- 20 #changes size of the node
    min_edge_value <- 0.25 #between 0 and 1, change how many nodes are connected. 
    cex_label_category_size <- 1 #controls the size of text label
    node_label_list <-'none' #label the notes on 'list' emap
    node_label_DB <- 'category' #label the notes on 'DB' emap
    node_alpha_value<- 0.7
    cex_line_thickness<- 1
    layout_function<-'kk'#use function 'kk' to keep plot static
    
    color_list1<-"#D55E00"
    color_GOBP<-"#D55E00"
    color_KEGG<-"#0072B2"
    color_REACT<-"#CC79A7"
    
    
    # ------------- shares select with clusterProfiler, so put last ---------- 
    
    # Make sure graphics device is enabled
    options(bitmapType="cairo")
    
    ## Record package usage  
    #Of particular use is the version of R (at the top), and the version of DESeq2 (buried within 'other attached packages').
    writeLines(capture.output(sessionInfo()), "RsessionInfo.txt")
    
    ensemblVer = useEnsembl(biomart = "ensembl", version = ensembl_db_ver)
    
    # see the available marts
    #listMarts(ensemblVer)
    
    # see all genomes in ensembl
    #View(listDatasets(ensemblVer))
    
    ensemblVer = useDataset("hsapiens_gene_ensembl", 
                            mart = ensemblVer)
    
    
    ### Setup output for ORA output
    
    # dir.create("./ORA")
    
    ### clusterProfiler ORA (over-representation analysis/enrichment) for different gene groupings Uses Ensembl gene IDs - use of Entrez Gene ids (from biomaRt) results in multiple matches between nomenclatures.  
    
    
    
    
    print("## candidate list 1 analysis")
    # Set up the 'Universe' for list 1
    
    ## Universe list of ensembl genes (e.g. ENSG...) will be input as a plain-text file containing ENSG or ENSG.ver
    # Import gene list as a data.frame 
    
    # universe <- read_delim(file = universe_file_name_1, delim="\n", col_names=FALSE)
    
    ## getting the data from the uploaded file
    universe <- RV$universe_file_name_1
    
    
    # Extract column of genes from data.frame to make a vector
    universe <- universe$X1  
    
    ##Universe set up from versioned ENSG - could go directly from entrez_gene IDs...?
    
    # available keytypes 'keytypes(org.Hs.eg.db)'
    
    # Remove versions from gene names to ensure matching
    universe_nover <- str_replace(universe, pattern = "\\.[0-9]+$", replacement = "")
    universe_nover <- as.character(universe_nover)
    
    print("# Get entrezgene ids")# Get entrezgene ids
    
    universe_entrezgene <- getBM(attributes=c('entrezgene_id'),
                                 filters = c('ensembl_gene_id'),
                                 values = universe_nover,
                                 mart = ensemblVer)
    
    universe_entrezgene_unique <- as.data.frame(universe_entrezgene) %>% drop_na() %>% distinct()
    
    universe_entrezgene_unique <- as.character(universe_entrezgene_unique$entrezgene_id)
    
    
    
    
    print("getting list of candidate genes for list 1")
    # getting list of candidate genes for list 1
    
    #input .txt file with list of candidate ESNG.... here its set to a global variable
    # candidates_ver <- read_delim(file = candidate_file_name_1, delim="\n", col_names=FALSE)
    
    
    ## getting the data from the uploaded file
    candidates_ver <- RV$candidate_file_name_1
    
    candidates_ver <- candidates_ver$X1  
    
    # Remove Ensembl gene version if present
    candidates <- str_replace(candidates_ver, pattern = "\\.[0-9]+$", replacement = "")
    
    # Significant gene list
    #candidates <- candidates_mutated %>% filter(test == TRUE) %>% pull(entrezgene_id)
    #candidates <- candidates_mutated %>% filter(test == TRUE) %>% pull(ensgene_nover)
    candidates <- as.character(candidates)
    
    # Convert to entrezgene
    candidates_df <- data.frame(ensgene = candidates)
    
    candidates_entrezgene <- getBM(attributes=c('entrezgene_id'),
                                   filters = c('ensembl_gene_id'),
                                   values = candidates_df$ensgene,
                                   mart = ensemblVer)
    
    candidates_entrezgene_unique <- candidates_entrezgene %>% 
      drop_na() %>% 
      distinct()
    
    candidates_entrezgene_unique <- as.character(candidates_entrezgene$entrezgene_id)
    
    # print(candidates_entrezgene_unique)
    
    print("# Run all clusterProfiler for candidate gene list 1")
    
    # Run all clusterProfiler for candidate gene list 1
    
    # clusterProfiler Biological Property
    CP_ego_BP_candidates_1 = enrichGO(gene = candidates_entrezgene_unique,
                                      OrgDb = org.Hs.eg.db,
                                      pvalueCutoff = pvalue_cut_off,
                                      pAdjustMethod = "BH",
                                      ont = "BP",
                                      keyType = "ENTREZID",
                                      universe = universe_entrezgene_unique, 
                                      qvalueCutoff = qvalue_cut_off)
    
    
    CP_ego_BP_candidates_symbol_1 <- setReadable(CP_ego_BP_candidates_1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    
    
    # clusterProfiler KEGG
    # Not with Ensembl IDs - requires Entrez Gene IDs  
    CP_ekegg_candidates_1 = enrichKEGG(gene = candidates_entrezgene_unique,
                                       organism = "hsa",
                                       pvalueCutoff = pvalue_cut_off,
                                       pAdjustMethod = "BH",
                                       universe = universe_entrezgene_unique, 
                                       qvalueCutoff = qvalue_cut_off)
    
    CP_ekegg_candidates_symbol_1 <- setReadable(CP_ekegg_candidates_1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    
    # clusterProfiler Reactome
    # Not with Ensembl IDs - requires Entrez Gene IDs  
    CP_ereact_candidates_1 = enrichPathway(gene = candidates_entrezgene_unique,
                                           organism = "human",
                                           pvalueCutoff = pvalue_cut_off,
                                           pAdjustMethod = "BH",
                                           universe = universe_entrezgene_unique, 
                                           qvalueCutoff = qvalue_cut_off)
    
    CP_ereact_candidates_symbol_1 <- setReadable(CP_ereact_candidates_1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    
    
    
    ## save write xls here
    
    
    
    print("# Generate excel outputs for candidate gene list 1")
    # Generate excel outputs for candidate gene list 1
    
    # Output to plain text and tab-delimited file (.xls helps it to be opened in a spreadsheet)
    
    # write_tsv(as.data.frame(CP_ego_BP_candidates_symbol_1),"ORA/CP_ego_BP_candidates_1.xls")
    # write_tsv(as.data.frame(CP_ekegg_candidates_symbol_1),"ORA/CP_ekegg_candidates_1.xls")
    # write_tsv(as.data.frame(CP_ereact_candidates_symbol_1),"ORA/CP_ereact_candidates_1.xls")
    
    ## data saved for downloading: will be downloaded when clicked on button in Downloads
    RV$CP_ego_BP_candidates_symbol_1 =  as.data.frame(CP_ego_BP_candidates_symbol_1)
    RV$CP_ekegg_candidates_symbol_1 = as.data.frame(CP_ekegg_candidates_symbol_1)
    RV$CP_ereact_candidates_symbol_1 = as.data.frame(CP_ereact_candidates_symbol_1)
    
    
    
    
    
    
    print("## simplifyGOBP")
    ## simplifyGOBP
    
    
    CP_ego_BP_candidates_symbol_1_filtered <- CP_ego_BP_candidates_symbol_1 %>% filter(p.adjust < 0.01)
    CP_ego_BP_candidates_symbol_1_filtered_pw <- pairwise_termsim(CP_ego_BP_candidates_symbol_1_filtered)
    GOSIM1 <- clusterProfiler::simplify(CP_ego_BP_candidates_symbol_1_filtered_pw,cutoff=0.7, by="p.adjust",select_fun=min)
    
    
    
    
    ## Merge all results
    
    
    MERGEALL<-merge_result(list("GOBP_1"=GOSIM1,"KEGG_1"=CP_ekegg_candidates_symbol_1,"REACT_1"=CP_ereact_candidates_symbol_1)) 
    
    #make a combined emap with labels FOR MERGEMERGE
    filtered <- MERGEALL %>% filter(p.adjust < 0.05)
    filtered_pw <- pairwise_termsim(filtered)
    
    # write_tsv(as.data.frame(filtered_pw),"ORA/MERGEALL_filtered.xls")
    ## data saved for downloading: will be downloaded when clicked on button in Downloads
    RV$filtered_pw =  as.data.frame(filtered_pw)
    
    ## saving the object,filtered_pw, for plotting in the Results
  
    RV$Plot1filtered_pw = filtered_pw
    
    # MERGEALL_COLORBYLIST_emapplot <- emapplot(filtered_pw,showCategory=200,layout=layout_function,cex_category=cex_category_size,min_edge=min_edge_value,cex_label_category = cex_label_category_size,repel = TRUE,cex_line=cex_line_thickness,node_label=node_label_list)+scale_fill_manual(values=alpha(c(color_list1,color_list1,color_list1),node_alpha_value))
    
    # png(file = "ORA/MERGEALL_COLORBYLIST_emapplot.png", width = 1684, height = 1190, unit = "px")
    # MERGEALL_COLORBYLIST_emapplot
    # dev.off()
    # RStudioGD 
    # 2 
    # 
    # 
    # 
  
    
    # MERGEALL_COLORBYDB_emapplot2 <- emapplot(filtered_pw,showCategory=200,layout=layout_function,cex_category=cex_category_size,min_edge=min_edge_value,cex_label_category = cex_label_category_size,repel = TRUE,cex_line=cex_line_thickness,node_label=node_label_DB)+scale_fill_manual(values=alpha(c(color_GOBP,color_KEGG,color_REACT),node_alpha_value))
    
    # png(file = "ORA/MERGEALL_COLORBYDB_emapplot.png", width = 1684, height = 1190, unit = "px")
    # MERGEALL_COLORBYDB_emapplot
    # dev.off()
    # 
    # 
    # dev.off()
    # RStudioGD 
    # 2 
    
    
    
    
  })
  
  
  
  ## ------------ Results Page -------------------
  
  ## using ggsave to save the plots in www folder then loading the images
  observe({
    req(RV$Plot1filtered_pw)
    
    print("Plot Ob")
    
    filtered_pw = RV$Plot1filtered_pw
    
    MERGEALL_COLORBYLIST_emapplot <- emapplot(filtered_pw,showCategory=200,layout=layout_function,cex_category=cex_category_size,min_edge=min_edge_value,cex_label_category = cex_label_category_size,repel = TRUE,cex_line=cex_line_thickness,node_label=node_label_list)+scale_fill_manual(values=alpha(c(color_list1,color_list1,color_list1),node_alpha_value))
    ggsave("www/MERGEALL_COLORBYLIST_emapplot.png",width = 3000, height = 3000, unit = "px")
    
    MERGEALL_COLORBYDB_emapplot <- emapplot(filtered_pw,showCategory=200,layout=layout_function,cex_category=cex_category_size,min_edge=min_edge_value,cex_label_category = cex_label_category_size,repel = TRUE,cex_line=cex_line_thickness,node_label=node_label_DB)+scale_fill_manual(values=alpha(c(color_GOBP,color_KEGG,color_REACT),node_alpha_value))
    ggsave("www/MERGEALL_COLORBYDB_emapplot.png",width = 3000, height = 3000, unit = "px")
    
    ## removes the waiting message after simulation is complete and plots are saved
    removeModal()
    
  })
  
  
  output$Results = renderUI({
    req(input$Run)
    UI = list()
    UI = append(UI, list(box(width = 12,h2("Plot 1"))))
    UI = append(UI, list(box(width = 12,img(src='MERGEALL_COLORBYLIST_emapplot.png', align = "center",height="100%", width="100%"))))
    UI = append(UI, list(box(width = 12,h2("Plot 2"))))
    UI = append(UI, list(box(width = 12,img(src='MERGEALL_COLORBYDB_emapplot.png', align = "center",height="100%", width="80%"))))
    UI
  })  
  
  ## ------------ Download Page -------------------
  
  output$Downloads = renderUI({
    req(input$Run)
    
    UI = list()
    UI = append(UI, list(box( width = 3,downloadButton("CP_ego_BP_candidates_symbol_1", "Download CP_ego_BP_candidates_symbol_1"))))
    UI = append(UI, list(box( width = 3,downloadButton("CP_ekegg_candidates_symbol_1", "Download CP_ekegg_candidates_symbol_1"))))
    UI = append(UI, list(box( width = 3,downloadButton("CP_ereact_candidates_symbol_1", "Download CP_ereact_candidates_symbol_1"))))
    UI = append(UI, list(box( width = 3,downloadButton("filtered_pw", "Download filtered_pw"))))
    
    UI = append(UI, list(box( width = 3,downloadButton("MERGEALL_COLORBYLIST_emapplot", "Download MERGEALL_COLORBYLIST_emapplot"))))
    UI = append(UI, list(box( width = 3,downloadButton("MERGEALL_COLORBYDB_emapplot", "Download MERGEALL_COLORBYDB_emapplot"))))
    
    UI
  })
  
  

  output$CP_ego_BP_candidates_symbol_1 <- downloadHandler(
    filename = function() {
      "CP_ego_BP_candidates_symbol_1.xls"
    },
    content = function(con) {
      write_tsv(RV$CP_ego_BP_candidates_symbol_1,con)
    }
  )
  
  output$CP_ekegg_candidates_symbol_1 <- downloadHandler(
    filename = function() {
      "CP_ekegg_candidates_symbol_1.xls"
    },
    content = function(con) {
      write_tsv(RV$CP_ekegg_candidates_symbol_1,con)
    }
  )
  
  output$CP_ereact_candidates_symbol_1 <- downloadHandler(
    filename = function() {
      "CP_ereact_candidates_symbol_1.xls"
    },
    content = function(con) {
      write_tsv(RV$CP_ereact_candidates_symbol_1,con)
    }
  )
  
  output$filtered_pw <- downloadHandler(
    filename = function() {
      "filtered_pw.xls"
    },
    content = function(con) {
      write_tsv(RV$filtered_pw,con)
    }
  )  
  
  output$MERGEALL_COLORBYLIST_emapplot <- downloadHandler(
    filename = function() {
      "www/MERGEALL_COLORBYLIST_emapplot.png"
    },
    content = function(con) {
      filtered_pw = RV$Plot1filtered_pw
      MERGEALL_COLORBYLIST_emapplot <- emapplot(filtered_pw,showCategory=200,layout=layout_function,cex_category=cex_category_size,min_edge=min_edge_value,cex_label_category = cex_label_category_size,repel = TRUE,cex_line=cex_line_thickness,node_label=node_label_list)+scale_fill_manual(values=alpha(c(color_list1,color_list1,color_list1),node_alpha_value))
      ggsave(con,width = 3000, height = 3000, unit = "px")
    }
  )
  
  output$MERGEALL_COLORBYDB_emapplot <- downloadHandler(
    filename = function() {
      "www/MERGEALL_COLORBYDB_emapplot.png"
    },
    content = function(con) {
      filtered_pw = RV$Plot1filtered_pw
      MERGEALL_COLORBYDB_emapplot <- emapplot(filtered_pw,showCategory=200,layout=layout_function,cex_category=cex_category_size,min_edge=min_edge_value,cex_label_category = cex_label_category_size,repel = TRUE,cex_line=cex_line_thickness,node_label=node_label_DB)+scale_fill_manual(values=alpha(c(color_GOBP,color_KEGG,color_REACT),node_alpha_value))
      ggsave(con,width = 3000, height = 3000, unit = "px")
    }
  )
  
  
  
  
}

shinyApp(ui, server)

