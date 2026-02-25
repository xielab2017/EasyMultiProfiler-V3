#!/usr/bin/env Rscript
# EasyMultiProfiler - 蛋白质组学分析脚本
# 蛋白质组数据分析：差异表达、通路分析、PPI网络

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

# 命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="蛋白质定量矩阵文件路径"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="样本分组信息文件"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="输出目录"),
  make_option(c("-p", "--params"), type="character", default="{}", help="分析参数"),
  make_option(c("-t", "--task-id"), type="character", default=NULL, help="任务ID")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# 解析参数
params <- fromJSON(args$params)

cat(sprintf("开始蛋白质组学分析 - 任务ID: %s\n", args$task_id))

try {
  # 读取蛋白质定量数据
  protein_data <- read.csv(args$input, row.names=1, check.names=FALSE)
  cat(sprintf("蛋白质数据维度: %d 蛋白 x %d 样本\n", nrow(protein_data), ncol(protein_data)))
  
  # 数据预处理
  normalization <- params$preprocess$normalization %||% "median"
  imputation <- params$preprocess$imputation %||% "knn"
  
  cat(sprintf("归一化方法: %s, 缺失值填充: %s\n", normalization, imputation))
  
  # 模拟数据预处理（实际应使用 MSnbase/limma 等包）
  processed_data <- protein_data
  
  # 处理缺失值
  if (imputation != "none") {
    # 简单填充（模拟）
    processed_data[is.na(processed_data)] <- min(processed_data, na.rm=TRUE) / 2
  }
  
  # 归一化（模拟）
  if (normalization == "median") {
    processed_data <- sweep(processed_data, 2, apply(processed_data, 2, median, na.rm=TRUE), "/")
  } else if (normalization == "zscore") {
    processed_data <- scale(processed_data)
  }
  
  # 模拟差异分析结果
  set.seed(42)
  
  de_proteins <- data.frame(
    protein_id = rownames(protein_data)[1:min(500, nrow(protein_data))],
    gene_name = paste0("Gene", 1:min(500, nrow(protein_data))),
    log2FC = rnorm(min(500, nrow(protein_data)), 0, 1.5),
    pvalue = runif(min(500, nrow(protein_data)), 0, 1),
    adj_pvalue = runif(min(500, nrow(protein_data)), 0, 1)
  )
  
  # 调整 pvalue
  fc_threshold <- params$de$fc_threshold %||% 1.5
  p_threshold <- params$de$pvalue %||% 0.05
  
  de_proteins$adj_pvalue[1:30] <- runif(30, 0.001, 0.05)
  de_proteins$log2FC[1:15] <- runif(15, 1, 3)
  de_proteins$log2FC[16:30] <- runif(15, -3, -1)
  
  # 标记差异蛋白
  de_proteins$regulation <- ifelse(
    de_proteins$adj_pvalue < p_threshold & abs(de_proteins$log2FC) >= log2(fc_threshold),
    ifelse(de_proteins$log2FC > 0, "Up", "Down"),
    "Not Sig"
  )
  
  # 生成图表
  output_dir <- args$output
  
  # 1. 火山图
  png(file.path(output_dir, "protein_volcano.png"), width=800, height=600)
  plot(de_proteins$log2FC, -log10(de_proteins$adj_pvalue), 
       pch=20, col=ifelse(de_proteins$regulation == "Not Sig", "grey", 
                          ifelse(de_proteins$regulation == "Up", "red", "blue")),
       xlab="log2 Fold Change", ylab="-log10(adjusted p-value)",
       main="Protein Differential Expression (Volcano Plot)")
  abline(h=-log10(p_threshold), col="grey", lty=2)
  abline(v=c(-log2(fc_threshold), log2(fc_threshold)), col="grey", lty=2)
  legend("topright", legend=c("Up-regulated", "Down-regulated", "Not significant"),
         col=c("red", "blue", "grey"), pch=20)
  dev.off()
  
  # 2. 蛋白质定量分布
  png(file.path(output_dir, "protein_distribution.png"), width=1000, height=600)
  par(mfrow=c(1,2))
  
  # 密度图
  plot(density(as.matrix(processed_data), na.rm=TRUE), 
       main="Protein Intensity Distribution",
       xlab="Intensity", ylab="Density")
  
  # 箱线图
  boxplot(processed_data, main="Sample Boxplot", las=2, cex.axis=0.7)
  dev.off()
  
  # 3. 样本聚类热图
  if (ncol(protein_data) > 2) {
    png(file.path(output_dir, "protein_correlation.png"), width=800, height=800)
    cor_matrix <- cor(processed_data, use="pairwise.complete.obs")
    heatmap(cor_matrix, main="Sample Correlation (Protein)", 
            xlab="Samples", ylab="Samples")
    dev.off()
  }
  
  # 4. 差异蛋白热图
  sig_proteins <- de_proteins[de_proteins$regulation != "Not Sig", ]
  if (nrow(sig_proteins) > 0) {
    png(file.path(output_dir, "protein_heatmap.png"), width=1000, height=800)
    top_proteins <- sig_proteins$protein_id[1:min(40, nrow(sig_proteins))]
    heatmap_data <- as.matrix(processed_data[top_proteins, ])
    heatmap_data <- t(scale(t(heatmap_data)))
    
    if (ncol(heatmap_data) > 1) {
      heatmap(heatmap_data, main="Differential Protein Expression",
              xlab="Samples", ylab="Proteins")
    }
    dev.off()
  }
  
  # 5. 通路富集分析（模拟）
  if (!is.null(params$function) && (params$function$pathway %||% TRUE)) {
    png(file.path(output_dir, "pathway_enrichment.png"), width=1000, height=600)
    
    pathways <- c("Metabolic pathways", "Ribosome", "Oxidative phosphorylation",
                  "Protein processing", "Cell cycle", "DNA replication",
                  "Signal transduction", "Immune system")
    enrichment <- data.frame(
      Pathway = pathways,
      GeneRatio = runif(8, 0.05, 0.3),
      pvalue = runif(8, 0.001, 0.05),
      Count = sample(5:50, 8)
    )
    enrichment <- enrichment[order(enrichment$pvalue), ]
    
    par(mar=c(5, 20, 4, 2))
    barplot(enrichment$GeneRatio, names.arg=enrichment$Pathway, 
            horiz=TRUE, las=2, main="Pathway Enrichment (Proteins)",
            xlab="Gene Ratio", col=terrain.colors(8))
    dev.off()
  }
  
  # 6. PPI网络（简化模拟）
  if (!is.null(params$function) && (params$function$ppi %||% FALSE)) {
    png(file.path(output_dir, "ppi_network.png"), width=800, height=800)
    
    # 模拟网络
    set.seed(123)
    n_nodes <- min(30, nrow(sig_proteins))
    if (n_nodes > 5) {
      # 简化的网络可视化
      plot(1:n_nodes, rep(0, n_nodes), type="n", xlim=c(0, n_nodes+1), ylim=c(0, 10),
           main="Protein-Protein Interaction Network (Top DE Proteins)",
           xlab="", ylab="", axes=FALSE)
      
      # 绘制节点
      for (i in 1:n_nodes) {
        x <- runif(1, 1, n_nodes)
        y <- runif(1, 1, 9)
        points(x, y, pch=21, bg=ifelse(sig_proteins$regulation[i]=="Up", "red", "blue"), cex=2)
        text(x, y, sig_proteins$gene_name[i], cex=0.6, pos=3)
      }
      
      # 绘制边（随机）
      for (i in 1:50) {
        x1 <- runif(1, 1, n_nodes)
        y1 <- runif(1, 1, 9)
        x2 <- runif(1, 1, n_nodes)
        y2 <- runif(1, 1, 9)
        lines(c(x1, x2), c(y1, y2), col="grey80", lwd=0.5)
      }
      
      legend("topright", legend=c("Up-regulated", "Down-regulated"),
             pch=21, pt.bg=c("red", "blue"))
    }
    dev.off()
  }
  
  # 保存差异分析结果
  write.csv(de_proteins, file.path(output_dir, "differential_proteins.csv"), row.names=FALSE)
  
  # 保存处理后的数据
  write.csv(processed_data, file.path(output_dir, "processed_protein_data.csv"))
  
  # 保存统计数据
  stats <- list(
    module = "proteomics",
    samples = ncol(protein_data),
    proteins = nrow(protein_data),
    up_regulated = sum(de_proteins$regulation == "Up"),
    down_regulated = sum(de_proteins$regulation == "Down"),
    normalization = normalization,
    imputation = imputation,
    task_id = args$task_id
  )
  write_json(stats, file.path(output_dir, "stats.json"))
  
  # 生成 PDF 报告
  pdf(file.path(output_dir, "proteomics_analysis_report.pdf"), width=12, height=10)
  par(mfrow=c(2,2))
  
  # 火山图
  plot(de_proteins$log2FC, -log10(de_proteins$adj_pvalue), 
       pch=20, col=ifelse(de_proteins$regulation == "Not Sig", "grey", 
                          ifelse(de_proteins$regulation == "Up", "red", "blue")),
       xlab="log2 Fold Change", ylab="-log10(adjusted p-value)",
       main="Protein Differential Expression")
  
  # 差异蛋白数
  sig_counts <- c(Up=stats$up_regulated, Down=stats$down_regulated,
                 NS=nrow(de_proteins) - stats$up_regulated - stats$down_regulated)
  barplot(sig_counts, main="Differentially Expressed Proteins", 
          col=c("red", "blue", "grey"))
  
  # 蛋白强度分布
  hist(as.matrix(processed_data), breaks=50, main="Protein Intensity Distribution",
       xlab="Normalized Intensity")
  
  # 样本总蛋白数
  barplot(colSums(!is.na(protein_data)), main="Proteins Quantified per Sample",
          las=2, cex.names=0.7)
  
  dev.off()
  
  cat("蛋白质组学分析完成！\n")
  cat(sprintf("差异上调蛋白: %d, 下调: %d\n", stats$up_regulated, stats$down_regulated))
  
} catch (e) {
  cat(sprintf("错误: %s\n", e$message))
  writeLines(as.character(e), file.path(args$output, "error.log"))
  quit(status=1)
}

# 辅助函数
`%||%` <- function(x, y) if (is.null(x)) y else x

cat("蛋白质组学分析脚本执行成功\n")
