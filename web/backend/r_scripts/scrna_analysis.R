#!/usr/bin/env Rscript
# EasyMultiProfiler - 单细胞RNA-seq分析脚本 (调用Seurat等R包)
# 完整的单细胞分析流程

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

# 命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="基因表达矩阵文件路径"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="样本元数据文件"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="输出目录"),
  make_option(c("-p", "--params"), type="character", default="{}", help="分析参数"),
  make_option(c("-t", "--task-id"), type="character", default=NULL, help="任务ID")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# 解析参数
params <- fromJSON(args$params)

cat(sprintf("开始单细胞RNA-seq分析 - 任务ID: %s\n", args$task_id))

# 检查并安装依赖
check_dependencies <- function() {
  packages <- c("Seurat", "dplyr", "ggplot2")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("📦 安装 %s...\n", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  }
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  
  # 可选依赖
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    cat("⚠️  SingleR 未安装，细胞注释功能将不可用\n")
  }
  cat("✅ 依赖包检查完成\n")
}

tryCatch({
  check_dependencies()
}, error = function(e) {
  cat(sprintf("❌ 依赖安装失败: %s\n", e$message))
  quit(status=1)
})

try {
  # 读取数据
  cat("📁 读取表达矩阵...\n")
  counts <- read.csv(args$input, row.names=1, check.names=FALSE)
  cat(sprintf("   维度: %d 基因 x %d 细胞\n", nrow(counts), ncol(counts)))
  
  # 读取或创建元数据
  if (!is.null(args$metadata) && file.exists(args$metadata)) {
    metadata <- read.csv(args$metadata, row.names=1)
    cat(sprintf("   读取元数据: %d 样本\n", nrow(metadata)))
  } else {
    metadata <- NULL
  }
  
  # 获取参数
  min_genes <- params$qc$min_genes %||% 200
  min_cells <- params$qc$min_cells %||% 3
  max_mt <- params$qc$max_mt_percent %||% 5
  resolution <- params$cluster$resolution %||% 0.8
  dims <- 1:(params$cluster$dims %||% 30)
  min_pct <- params$markers$min_pct %||% 0.25
  logfc <- params$markers$logfc %||% 0.25
  
  # Step 1: 创建Seurat对象
  cat("🧫 Step 1: 创建Seurat对象...\n")
  seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, project = "EasyMultiProfiler")
  cat(sprintf("   初始细胞数: %d\n", ncol(seurat_obj)))
  
  # Step 2: 质控
  cat("🔍 Step 2: 质量控制...\n")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # QC前统计
  cat(sprintf("   QC前: %d 细胞\n", ncol(seurat_obj)))
  
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > min_genes & 
                         percent.mt < max_mt)
  
  cat(sprintf("   QC后: %d 细胞 (过滤掉 %d 低质量细胞)\n", 
              ncol(seurat_obj), ncol(counts) - ncol(seurat_obj)))
  
  # Step 3: 标准化
  cat("⚗️  Step 3: 标准化...\n")
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  cat(sprintf("   高变基因数: %d\n", length(VariableFeatures(object = seurat_obj))))
  
  # Step 4: 缩放和PCA
  cat("📊 Step 4: 降维 (PCA)...\n")
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  # 保存PCA图
  png(file.path(args$output, "pca_plot.png"), width=800, height=600)
  print(DimPlot(seurat_obj, reduction = "pca"))
  dev.off()
  
  # Step 5: 聚类
  cat("🔗 Step 5: 聚类分析...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  cat(sprintf("   发现 %d 个聚类\n", length(unique(Idents(seurat_obj)))))
  
  # Step 6: UMAP
  cat("🗺️  Step 6: UMAP降维...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  
  # 保存UMAP图
  png(file.path(args$output, "umap_clusters.png"), width=800, height=600)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE))
  dev.off()
  
  # 如果有样本信息，绘制分样本UMAP
  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    png(file.path(args$output, "umap_by_sample.png"), width=1000, height=800)
    print(DimPlot(seurat_obj, reduction = "umap", group.by = "sample"))
    dev.off()
  }
  
  # Step 7: 标记基因
  cat("🧬 Step 7: 寻找标记基因...\n")
  markers <- FindAllMarkers(seurat_obj, 
                              only.pos = TRUE, 
                              min.pct = min_pct, 
                              logfc.threshold = logfc)
  
  # 每个cluster的top标记基因
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
  
  write.csv(top_markers, file.path(args$output, "cluster_markers.csv"), row.names=FALSE)
  cat(sprintf("   找到 %d 个标记基因\n", nrow(markers)))
  
  # Step 8: 可视化
  cat("📈 Step 8: 生成可视化...\n")
  
  # 标记基因热图
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  png(file.path(args$output, "marker_heatmap.png"), width=1000, height=1200)
  print(DoHeatmap(seurat_obj, features = top10$gene) + NoLegend())
  dev.off()
  
  # 小提琴图 - nCount 和 nFeature
  png(file.path(args$output, "qc_violin.png"), width=1000, height=600)
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p1)
  dev.off()
  
  # FeaturePlot - 标记基因
  top_genes <- markers %>% group_by(cluster) %>% slice_max(n = 1, order_by = avg_log2FC) %>% pull(gene)
  if (length(top_genes) > 0) {
    png(file.path(args$output, "feature_plot.png"), width=1200, height=1000)
    print(FeaturePlot(seurat_obj, features = head(top_genes, 6)))
    dev.off()
  }
  
  # Step 9: 细胞注释 (如果有SingleR)
  if (requireNamespace("SingleR", quietly = TRUE) && 
      requireNamespace("celldex", quietly = TRUE)) {
    cat("🏷️  Step 9: 细胞类型注释...\n")
    tryCatch({
      ref <- celldex::HumanPrimaryCellAtlasData()
      predictions <- SingleR(test = GetAssayData(seurat_obj, slot = "data"),
                              ref = ref,
                              labels = ref$label.main)
      seurat_obj$cell_type <- predictions$labels
      
      # 带注释的UMAP
      png(file.path(args$output, "umap_annotated.png"), width=1000, height=800)
      print(DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE))
      dev.off()
      
      # 保存注释结果
      write.csv(predictions, file.path(args$output, "cell_annotation.csv"))
      cat("   细胞注释完成\n")
    }, error = function(e) {
      cat(sprintf("   细胞注释失败: %s\n", e$message))
    })
  }
  
  # Step 10: 保存统计信息
  stats <- list(
    module = "scrna",
    initial_cells = ncol(counts),
    filtered_cells = ncol(seurat_obj),
    n_clusters = length(unique(Idents(seurat_obj))),
    n_markers = nrow(markers),
    resolution = resolution,
    min_genes = min_genes,
    max_mt = max_mt,
    task_id = args$task_id
  )
  write_json(stats, file.path(args$output, "stats.json"))
  
  # 保存Seurat对象 (RDS格式)
  saveRDS(seurat_obj, file.path(args$output, "seurat_object.rds"))
  
  # 生成PDF报告
  pdf(file.path(args$output, "scRNAseq_report.pdf"), width=12, height=10)
  
  # 1. QC
  p_qc <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p_qc)
  
  # 2. UMAP聚类
  p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
  print(p_umap)
  
  # 3. 标记基因热图
  if (nrow(top10) > 0) {
    p_heatmap <- DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
    print(p_heatmap)
  }
  
  # 4. 特征图
  if (length(top_genes) >= 4) {
    p_feature <- FeaturePlot(seurat_obj, features = head(top_genes, 4))
    print(p_feature)
  }
  
  dev.off()
  
  cat("✅ 单细胞RNA-seq分析完成！\n")
  cat(sprintf("   结果保存在: %s\n", args$output))
  cat(sprintf("   聚类数: %d, 标记基因: %d\n", stats$n_clusters, stats$n_markers))
  
} catch (e) {
  cat(sprintf("❌ 分析失败: %s\n", e$message))
  writeLines(as.character(e), file.path(args$output, "error.log"))
  quit(status=1)
}

# 辅助函数
`%||%` <- function(x, y) if (is.null(x)) y else x

cat("单细胞RNA-seq分析脚本执行成功\n")
