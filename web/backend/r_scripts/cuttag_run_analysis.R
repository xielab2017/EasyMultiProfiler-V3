#!/usr/bin/env Rscript
# EasyMultiProfiler - CUT&Tag/CUT&RUN分析脚本
# 靶向切割分析：SEACR Peak calling, 富集分析, 可视化

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

# 命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Peak文件或计数矩阵"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="样本信息"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="输出目录"),
  make_option(c("-p", "--params"), type="character", default="{}", help="分析参数"),
  make_option(c("-t", "--task-id"), type="character", default=NULL, help="任务ID"),
  make_option(c("--mode"), type="character", default="cutntag", help="模式: cutntag 或 cutnrun")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# 解析参数
params <- fromJSON(args$params)
mode <- args$mode  # cutntag 或 cutnrun

cat(sprintf("开始%s分析 - 任务ID: %s\n", toupper(mode), args$task_id))

# 检查依赖
check_dependencies <- function() {
  bioc_packages <- c("ChIPseeker", "ChIPpeakAnno", "clusterProfiler", 
                     "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
  
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("📦 安装 %s...\n", pkg))
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE)
    }
  }
  
  library(ChIPseeker)
  library(ChIPpeakAnno)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  cat("✅ 依赖包检查完成\n")
}

tryCatch({
  check_dependencies()
}, error = function(e) {
  cat(sprintf("⚠️ 依赖问题: %s\n", e$message))
})

try {
  output_dir <- args$output
  
  # 读取数据
  cat("📁 读取数据...\n")
  
  is_peak_file <- grepl("\\.(bed|narrowPeak|broadPeak|seacr)$", args$input, ignore.case = TRUE)
  
  if (is_peak_file) {
    # 读取SEACR/MACS2输出的Peak文件
    cat(sprintf("   读取Peak文件: %s\n", args$input))
    
    peaks <- tryCatch({
      ChIPseeker::readPeakFile(args$input)
    }, error = function(e) {
      # 手动读取BED
      bed <- read.table(args$input, sep="\t", header=FALSE)[, 1:3]
      colnames(bed) <- c("chr", "start", "end")
      GenomicRanges::GRanges(seqnames=bed$chr, 
                             ranges=IRanges::IRanges(bed$start, bed$end))
    })
    
    cat(sprintf("   读取 %d peaks\n", length(peaks)))
    
    # Step 1: Peak质量评估
    cat("📊 Step 1: Peak质量评估...\n")
    
    # Peak长度分布
    peak_widths <- width(peaks)
    
    png(file.path(output_dir, "peak_width_distribution.png"), width=800, height=600)
    hist(peak_widths, breaks=50, main="Peak Width Distribution",
         xlab="Peak Width (bp)", ylab="Frequency")
    abline(v=median(peak_widths), col="red", lty=2, lwd=2)
    legend("topright", legend=paste("Median:", round(median(peak_widths)), "bp"),
           col="red", lty=2)
    dev.off()
    
    cat(sprintf("   Peak长度: median=%d bp, mean=%d bp\n", 
                median(peak_widths), round(mean(peak_widths))))
    
    # Step 2: Peak注释
    cat("📝 Step 2: Peak注释...\n")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    peak_anno <- annotatePeak(peaks, 
                               tssRegion=c(-3000, 3000),
                               TxDb=txdb,
                               annoDb="org.Hs.eg.db")
    
    write.csv(as.data.frame(peak_anno), file.path(output_dir, "peak_annotation.csv"))
    
    # 注释可视化
    png(file.path(output_dir, "annotation_pie.png"), width=800, height=600)
    print(plotAnnoPie(peak_anno))
    dev.off()
    
    png(file.path(output_dir, "annotation_bar.png"), width=800, height=600)
    print(plotAnnoBar(peak_anno))
    dev.off()
    
    png(file.path(output_dir, "tss_distribution.png"), width=800, height=600)
    print(plotDistToTSS(peak_anno))
    dev.off()
    
    # Step 3: 基因组分布
    cat("🗺️  Step 3: 基因组分布...\n")
    
    png(file.path(output_dir, "genome_coverage.png"), width=1200, height=800)
    tryCatch({
      print(covplot(peaks, weightCol="score"))
    }, error = function(e) {
      print(covplot(peaks, title=paste(toupper(mode), "Peak Coverage")))
    })
    dev.off()
    
    # Step 4: 转录因子/靶基因富集分析
    cat("🧬 Step 4: 富集分析...\n")
    
    genes <- unique(peak_anno@anno$geneId)
    genes <- genes[!is.na(genes)]
    
    if (length(genes) > 10) {
      # GO富集
      go_enrich <- enrichGO(gene = genes,
                             OrgDb = org.Hs.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             readable = TRUE)
      
      if (!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
        write.csv(go_enrich@result, file.path(output_dir, "go_enrichment.csv"))
        
        png(file.path(output_dir, "go_barplot.png"), width=1000, height=800)
        print(barplot(go_enrich, showCategory=20))
        dev.off()
        
        cat(sprintf("   GO富集: %d terms\n", nrow(go_enrich@result)))
      }
      
      # KEGG富集
      kegg_enrich <- tryCatch({
        enrichKEGG(gene = genes, organism = 'hsa', pvalueCutoff = 0.05)
      }, error = function(e) NULL)
      
      if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
        write.csv(kegg_enrich@result, file.path(output_dir, "kegg_enrichment.csv"))
        
        png(file.path(output_dir, "kegg_dotplot.png"), width=1000, height=800)
        print(dotplot(kegg_enrich, showCategory=20))
        dev.off()
        
        cat(sprintf("   KEGG: %d pathways\n", nrow(kegg_enrich@result)))
      }
    }
    
    # Step 5: CUT&Tag特异性分析
    if (mode == "cutntag") {
      cat("✂️  Step 5: CUT&Tag特异性分析...\n")
      
      # Peak在启动子区的比例 (CUT&Tag特点)
      promoter_peaks <- sum(peak_anno@anno$annotation == "Promoter (<=1kb)") +
        sum(peak_anno@anno$annotation == "Promoter (1-2kb)") +
        sum(peak_anno@anno$annotation == "Promoter (2-3kb)")
      
      promoter_ratio <- promoter_peaks / length(peaks) * 100
      
      cat(sprintf("   启动子区Peak比例: %.1f%% (CUT&Tag特异性指标)\n", promoter_ratio))
      
      # 保存QC指标
      qc_metrics <- data.frame(
        Metric = c("Total_Peaks", "Median_Peak_Width", "Promoter_Peaks", 
                   "Promoter_Ratio", "Gene_Body_Peaks", "Intergenic_Peaks"),
        Value = c(length(peaks), median(peak_widths), promoter_peaks,
                  sprintf("%.2f%%", promoter_ratio),
                  sum(grepl("Intron|Exon|UTR", peak_anno@anno$annotation)),
                  sum(peak_anno@anno$annotation == "Distal Intergenic"))
      )
      write.csv(qc_metrics, file.path(output_dir, "qc_metrics.csv"), row.names=FALSE)
    }
    
    # 保存统计
    stats <- list(
      module = mode,
      n_peaks = length(peaks),
      median_peak_width = median(peak_widths),
      annotated_genes = length(unique(genes)),
      promoter_ratio = ifelse(exists("promoter_ratio"), promoter_ratio, NA),
      task_id = args$task_id
    )
    
  } else {
    # 计数矩阵模式
    cat("   读取计数矩阵...\n")
    count_matrix <- read.csv(args$input, row.names=1, check.names=FALSE)
    
    png(file.path(output_dir, "count_distribution.png"), width=800, height=600)
    hist(as.matrix(count_matrix), breaks=50, 
         main=paste(toupper(mode), "Peak Count Distribution"))
    dev.off()
    
    stats <- list(
      module = mode,
      n_peaks = nrow(count_matrix),
      n_samples = ncol(count_matrix),
      task_id = args$task_id
    )
  }
  
  write_json(stats, file.path(output_dir, "stats.json"))
  
  # 生成PDF报告
  pdf(file.path(output_dir, paste0(mode, "_report.pdf")), width=12, height=10)
  
  if (exists("peak_anno")) {
    print(plotAnnoPie(peak_anno))
    print(plotAnnoBar(peak_anno))
    print(plotDistToTSS(peak_anno))
  }
  
  if (exists("go_enrich") && !is.null(go_enrich) && nrow(go_enrich@result) > 0) {
    print(barplot(go_enrich, showCategory=10))
  }
  
  dev.off()
  
  cat("✅ 分析完成！\n")
  cat(sprintf("   结果保存在: %s\n", output_dir))
  
} catch (e) {
  cat(sprintf("❌ 分析失败: %s\n", e$message))
  writeLines(as.character(e), file.path(args$output, "error.log"))
  quit(status=1)
}

# 辅助函数
`%||%` <- function(x, y) if (is.null(x)) y else x

cat(sprintf("%s分析脚本执行成功\n", toupper(mode)))
