#!/usr/bin/env Rscript
# EasyMultiProfiler - RNA-seq 分析脚本 (基于验证的DESeq2标准流程)
# 转录组数据分析：差异表达、富集分析、可视化

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# 命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="基因表达矩阵文件路径"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="样本分组信息文件"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="输出目录"),
  make_option(c("-p", "--params"), type="character", default="{}", help="分析参数"),
  make_option(c("-t", "--task-id"), type="character", default=NULL, help="任务ID")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# 解析参数
params <- fromJSON(args$params)

cat(sprintf("开始 RNA-seq 分析 - 任务ID: %s\n", args$task_id))
cat("使用验证的标准DESeq2流程\n\n")

# 确保输出目录存在
if (!is.null(args$output)) {
  dir.create(args$output, recursive=TRUE, showWarnings=FALSE)
}

# 主分析流程
tryCatch({
  # 1. 读取表达矩阵
  cat("步骤1: 读取表达矩阵...\n")
  count_data <- read.csv(args$input, row.names=1, check.names=FALSE)
  cat(sprintf("✅ 表达矩阵: %d 基因 x %d 样本\n", nrow(count_data), ncol(count_data)))
  
  # 2. 读取分组信息
  cat("\n步骤2: 读取分组信息...\n")
  if (!is.null(args$metadata) && file.exists(args$metadata)) {
    group_info <- read.csv(args$metadata, row.names=1)
    cat(sprintf("✅ 分组信息: %d 样本\n", nrow(group_info)))
  } else {
    # 创建默认分组
    n_samples <- ncol(count_data)
    group_info <- data.frame(
      Group = factor(c(rep("Control", floor(n_samples/2)), rep("Treatment", ceiling(n_samples/2)))),
      row.names = colnames(count_data)
    )
    cat("✅ 使用默认分组: Control vs Treatment\n")
  }
  
  # 验证Mapping
  if (!all(colnames(count_data) %in% rownames(group_info))) {
    stop("样本名不匹配！请检查数据文件和元数据文件")
  }
  cat("✅ Mapping验证通过\n")
  
  # 3. 创建DESeq2对象
  cat("\n步骤3: 创建DESeq2对象...\n")
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = group_info,
    design = ~ Group
  )
  
  # 4. 过滤低表达基因
  cat("\n步骤4: 过滤低表达基因...\n")
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  cat(sprintf("✅ 过滤后: %d 基因\n", nrow(dds)))
  
  # 5. 运行DESeq2差异分析
  cat("\n步骤5: 运行DESeq2差异分析...\n")
  cat("   (这可能需要几分钟)...\n")
  dds <- DESeq(dds, quiet=TRUE)
  
  # 6. 获取结果
  cat("\n步骤6: 提取差异分析结果...\n")
  group_levels <- levels(dds$Group)
  cat(sprintf("对比组: %s vs %s\n", group_levels[2], group_levels[1]))
  res <- results(dds, contrast=c("Group", group_levels[2], group_levels[1]))
  
  # 转换为数据框并添加基因名
  de_table <- as.data.frame(res)
  de_table$GeneID <- rownames(de_table)
  de_table <- de_table[, c("GeneID", setdiff(colnames(de_table), "GeneID"))]
  
  # 计算显著性
  sig_table <- de_table[!is.na(de_table$padj) & de_table$padj < 0.05, ]
  up_genes <- sum(sig_table$log2FoldChange > 0, na.rm=TRUE)
  down_genes <- sum(sig_table$log2FoldChange < 0, na.rm=TRUE)
  
  cat(sprintf("✅ 差异分析完成\n"))
  cat(sprintf("   总基因: %d\n", nrow(de_table)))
  cat(sprintf("   显著差异基因 (p < 0.05): %d\n", nrow(sig_table)))
  cat(sprintf("   上调: %d | 下调: %d\n", up_genes, down_genes))
  
  # 7. 生成可视化
  cat("\n步骤7: 生成可视化...\n")
  output_dir <- args$output
  
  # 火山图
  tryCatch({
    de_table$significant <- ifelse(!is.na(de_table$padj) & de_table$padj < 0.05,
                                     ifelse(de_table$log2FoldChange > 0, "Up", "Down"),
                                     "Not Sig")
    
    p_volcano <- ggplot(de_table, aes(x=log2FoldChange, y=-log10(pvalue), color=significant)) +
      geom_point(alpha=0.6, size=1.5) +
      scale_color_manual(values=c("Up"="red", "Down"="blue", "Not Sig"="grey")) +
      theme_bw() +
      labs(title="Volcano Plot",
           x="Log2 Fold Change",
           y="-Log10 P-value",
           color="Significance") +
      theme(legend.position="top")
    
    ggsave(file.path(output_dir, "volcano_plot.png"), p_volcano, width=10, height=8, dpi=150)
    cat("✅ 火山图已保存: volcano_plot.png\n")
  }, error=function(e) {
    cat(sprintf("⚠️ 火山图生成失败: %s\n", e$message))
  })
  
  # MA图
  tryCatch({
    p_ma <- ggplot(de_table, aes(x=log10(baseMean + 1), y=log2FoldChange, color=significant)) +
      geom_point(alpha=0.6, size=1.5) +
      scale_color_manual(values=c("Up"="red", "Down"="blue", "Not Sig"="grey")) +
      theme_bw() +
      labs(title="MA Plot",
           x="Log10 Mean Expression",
           y="Log2 Fold Change") +
      theme(legend.position="none")
    
    ggsave(file.path(output_dir, "ma_plot.png"), p_ma, width=10, height=8, dpi=150)
    cat("✅ MA图已保存: ma_plot.png\n")
  }, error=function(e) {
    cat(sprintf("⚠️ MA图生成失败: %s\n", e$message))
  })
  
  # 热图 (Top 50差异基因)
  tryCatch({
    if (nrow(sig_table) >= 10) {
      top_genes <- head(order(sig_table$padj), min(50, nrow(sig_table)))
      top_gene_names <- sig_table$GeneID[top_genes]
      
      # 标准化计数
      vst_counts <- vst(dds, blind=FALSE)
      mat <- assay(vst_counts)[top_gene_names, ]
      mat <- mat - rowMeans(mat)
      
      # 保存热图
      png(file.path(output_dir, "heatmap.png"), width=800, height=1000)
      pheatmap(mat, 
               cluster_cols=TRUE, 
               cluster_rows=TRUE,
               annotation_col=as.data.frame(colData(dds)),
               show_rownames=FALSE,
               main="Top Differential Genes Heatmap")
      dev.off()
      cat("✅ 热图已保存: heatmap.png\n")
    }
  }, error=function(e) {
    cat(sprintf("⚠️ 热图生成失败: %s\n", e$message))
  })
  
  # 8. 保存结果
  cat("\n步骤8: 保存结果...\n")
  write.csv(de_table, file.path(output_dir, "differential_expression.csv"), row.names=FALSE)
  write.csv(sig_table, file.path(output_dir, "significant_genes.csv"), row.names=FALSE)
  cat("✅ 结果表格已保存\n")
  
  # 保存统计数据
  stats <- list(
    module = "rnaseq",
    samples = ncol(count_data),
    genes_total = nrow(count_data),
    genes_filtered = nrow(dds),
    differential_genes = nrow(sig_table),
    up_regulated = up_genes,
    down_regulated = down_genes,
    de_method = "DESeq2",
    group_comparison = paste(group_levels[2], "vs", group_levels[1]),
    task_id = args$task_id,
    status = "completed"
  )
  write_json(stats, file.path(output_dir, "stats.json"), pretty=TRUE)
  cat("✅ 统计信息已保存\n")
  
  # 生成简单的HTML报告
  html_report <- sprintf("
  <html>
  <head><title>RNA-seq分析报告</title></head>
  <body>
    <h1>EasyMultiProfiler RNA-seq 分析报告</h1>
    <p>任务ID: %s</p>
    <p>分析时间: %s</p>
    <h2>统计摘要</h2>
    <ul>
      <li>总样本数: %d</li>
      <li>总基因数: %d</li>
      <li>显著差异基因: %d</li>
      <li>上调基因: %d</li>
      <li>下调基因: %d</li>
    </ul>
    <h2>可视化结果</h2>
    <img src='volcano_plot.png' width='600'/><br/>
    <img src='ma_plot.png' width='600'/><br/>
    <img src='heatmap.png' width='600'/><br/>
  </body>
  </html>
  ", 
  args$task_id, format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stats$samples, stats$genes_total, stats$differential_genes,
  stats$up_regulated, stats$down_regulated)
  
  writeLines(html_report, file.path(output_dir, "report.html"))
  cat("✅ HTML报告已生成: report.html\n")
  
  cat("\n✅ RNA-seq 分析全部完成！\n")
  
}, error = function(e) {
  cat(sprintf("\n❌ 分析失败: %s\n", e$message))
  
  # 保存错误信息
  error_info <- list(
    task_id = args$task_id,
    status = "failed",
    error = e$message,
    timestamp = format(Sys.time())
  )
  write_json(error_info, file.path(args$output, "error.json"))
  
  quit(status=1)
})
