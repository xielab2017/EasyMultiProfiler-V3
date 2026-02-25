#!/usr/bin/env Rscript
# EasyMultiProfiler - RNA-seq 分析脚本 (调用 EasyMultiProfiler R包)
# 转录组数据分析：差异表达、富集分析、可视化

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
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

# 检查 EasyMultiProfiler 是否安装
check_emp <- function() {
  if (!requireNamespace("EasyMultiProfiler", quietly = TRUE)) {
    cat("⚠️  EasyMultiProfiler 包未安装，尝试安装...\n")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools", repos = "https://cloud.r-project.org/")
    }
    devtools::install_github("xielab2017/EasyMultiProfiler", upgrade = "never")
  }
  library(EasyMultiProfiler)
}

tryCatch({
  check_emp()
  cat("✅ EasyMultiProfiler 加载成功\n")
}, error = function(e) {
  cat(sprintf("❌ EasyMultiProfiler 加载失败: %s\n", e$message))
  cat("将使用基础分析流程...\n")
})

# 主分析流程
try {
  # 读取表达矩阵
  count_data <- read.csv(args$input, row.names=1, check.names=FALSE)
  cat(sprintf("表达矩阵维度: %d 基因 x %d 样本\n", nrow(count_data), ncol(count_data)))
  
  # 准备分组信息（从文件名或元数据文件）
  group_info <- NULL
  if (!is.null(args$metadata) && file.exists(args$metadata)) {
    group_info <- read.csv(args$metadata, row.names=1)
    cat(sprintf("读取分组信息: %d 样本\n", nrow(group_info)))
  } else {
    # 创建默认分组（假设前一半是A组，后一半是B组）
    n_samples <- ncol(count_data)
    group_info <- data.frame(
      Group = c(rep("Control", floor(n_samples/2)), rep("Treatment", ceiling(n_samples/2))),
      row.names = colnames(count_data)
    )
    cat("使用默认分组: Control vs Treatment\n")
  }
  
  # 使用 EasyMultiProfiler 进行分析
  if (exists("EMP_easy_normal_import") && exists("EMP_diff_analysis")) {
    
    cat("🔄 使用 EasyMultiProfiler 进行分析...\n")
    
    # 1. 导入数据创建 EMPT 对象
    cat("步骤1: 导入数据...\n")
    MAE <- EMP_easy_normal_import(
      data = count_data,
      assay = "rnaseq",
      assay_name = "counts",
      coldata = group_info,
      output = "MAE"
    )
    
    cat(sprintf("✅ 创建 MAE 对象: %d 特征 x %d 样本\n", 
                nrow(SummarizedExperiment::assay(MAE[[1]])),
                ncol(SummarizedExperiment::assay(MAE[[1]]))))
    
    # 2. 差异分析
    cat("步骤2: 差异表达分析...\n")
    de_method <- params$de$method %||% "DESeq2"
    fc_threshold <- params$de$fc_threshold %||% 2
    p_threshold <- params$de$pvalue %||% 0.05
    
    # 构建 formula
    group_col <- colnames(group_info)[1]
    formula_str <- paste0("~ 0 + ", group_col)
    
    # 获取组水平用于比较
    group_levels <- unique(group_info[[group_col]])
    if (length(group_levels) >= 2) {
      group_level_vec <- c(group_levels[1], group_levels[2])
    } else {
      group_level_vec <- NULL
    }
    
    # 执行差异分析
    diff_result <- tryCatch({
      MAE |\u003e 
        EMP_diff_analysis(
          experiment = "rnaseq",
          .formula = as.formula(formula_str),
          method = de_method,
          p.adjust = "fdr",
          group_level = group_level_vec,
          action = "add"
        )
    }, error = function(e) {
      cat(sprintf("差异分析警告: %s\n", e$message))
      # 降级使用 t.test
      cat("降级使用 t.test...\n")
      MAE |\u003e 
        EMP_diff_analysis(
          experiment = "rnaseq",
          method = "t.test",
          estimate_group = group_col,
          p.adjust = "fdr",
          action = "add"
        )
    })
    
    # 3. 获取差异分析结果
    cat("步骤3: 提取差异分析结果...\n")
    de_table <- tryCatch({
      diff_result |\u003e EMP_filter(.data$pvalue < p_threshold) |\u003e .get.result.EMPT()
    }, error = function(e) {
      # 直接从对象提取
      rowData(diff_result[[1]]) |\u003e as.data.frame()
    })
    
    cat(sprintf("✅ 差异分析完成: %d 差异基因\n", nrow(de_table)))
    
    # 4. 富集分析
    enrichment_enabled <- !is.null(params$enrichment) && (params$enrichment$database %||% "") != ""
    
    if (enrichment_enabled && exists("EMP_enrich")) {
      cat("步骤4: 富集分析...\n")
      
      db_type <- switch(params$enrichment$database %||% "go_kegg",
                         "go_kegg" = c("GO", "KEGG"),
                         "go_only" = "GO",
                         "kegg_only" = "KEGG",
                         "reactome" = "Reactome",
                         c("GO", "KEGG"))
      
      enrich_result <- tryCatch({
        diff_result |\u003e
          EMP_enrich(
            experiment = "rnaseq",
            OrgDb = "org.Hs.eg.db",  # 人类数据库
            type = db_type,
            pvalueCutoff = p_threshold,
            action = "add"
          )
      }, error = function(e) {
        cat(sprintf("富集分析警告: %s\n", e$message))
        NULL
      })
      
      if (!is.null(enrich_result)) {
        cat("✅ 富集分析完成\n")
      }
    }
    
    # 5. 生成可视化
    cat("步骤5: 生成可视化...\n")
    output_dir <- args$output
    
    # 火山图
    if (exists("EMP_volcano_plot")) {
      tryCatch({
        p_volcano <- diff_result |\u003e EMP_volcano_plot()
        ggplot2::ggsave(file.path(output_dir, "volcano_plot.png"), p_volcano, 
                       width = 10, height = 8, dpi = 300)
        cat("✅ 火山图已保存\n")
      }, error = function(e) {
        cat(sprintf("火山图生成失败: %s\n", e$message))
      })
    }
    
    # 热图
    if (exists("EMP_heatmap_plot")) {
      tryCatch({
        p_heatmap <- diff_result |\u003e EMP_heatmap_plot()
        ggplot2::ggsave(file.path(output_dir, "heatmap.png"), p_heatmap,
                       width = 12, height = 10, dpi = 300)
        cat("✅ 热图已保存\n")
      }, error = function(e) {
        cat(sprintf("热图生成失败: %s\n", e$message))
      })
    }
    
    # MA Plot
    tryCatch({
      png(file.path(output_dir, "ma_plot.png"), width=800, height=600)
      # 使用基础图形
      if ("log2FoldChange" %in% colnames(de_table) && "baseMean" %in% colnames(de_table)) {
        plot(log2(de_table$baseMean + 1), de_table$log2FoldChange,
             pch=20, col=ifelse(de_table$pvalue < p_threshold, "red", "grey"),
             xlab="log2 Mean Expression", ylab="log2 Fold Change",
             main="MA Plot")
        abline(h=c(-1, 0, 1), col=c("blue", "black", "blue"), lty=c(2,1,2))
      }
      dev.off()
      cat("✅ MA Plot 已保存\n")
    }, error = function(e) {
      cat(sprintf("MA Plot 生成失败: %s\n", e$message))
    })
    
    # 6. 保存结果
    cat("步骤6: 保存结果...\n")
    write.csv(de_table, file.path(output_dir, "differential_expression.csv"), row.names=FALSE)
    
    # 保存统计数据
    stats <- list(
      module = "rnaseq",
      samples = ncol(count_data),
      genes = nrow(count_data),
      differential_genes = nrow(de_table),
      up_regulated = sum(de_table$log2FoldChange > 0, na.rm=TRUE),
      down_regulated = sum(de_table$log2FoldChange < 0, na.rm=TRUE),
      de_method = de_method,
      fc_threshold = fc_threshold,
      pvalue_threshold = p_threshold,
      task_id = args$task_id
    )
    write_json(stats, file.path(output_dir, "stats.json"))
    
    # 生成报告
    if (exists("EMP_report")) {
      tryCatch({
        diff_result |\u003e EMP_report(output = file.path(output_dir, "report.html"))
        cat("✅ HTML 报告已生成\n")
      }, error = function(e) {
        cat(sprintf("报告生成失败: %s\n", e$message))
      })
    }
    
    cat("✅ RNA-seq 分析完成！\n")
    
  } else {
    # 降级到基础分析
    cat("⚠️  EasyMultiProfiler 函数不可用，使用基础分析...\n")
    source(file.path(dirname(getScriptPath()), "generic_analysis.R"))
  }
  
} catch (e) {
  cat(sprintf("❌ 分析失败: %s\n", e$message))
  writeLines(as.character(e), file.path(args$output, "error.log"))
  quit(status=1)
}

# 辅助函数
`%||%` <- function(x, y) if (is.null(x)) y else x

getScriptPath <- function() {
  cmd_args <- commandArgs(trailingOnly=FALSE)
  needle <- "--file="
  match <- grep(needle, cmd_args)
  if (length(match) > 0) {
    return(normalizePath(sub(needle, "", cmd_args[match])))
  }
  return(normalizePath(sys.frames()[[1]]$ofile))
}

cat("RNA-seq 分析脚本执行成功\n")
