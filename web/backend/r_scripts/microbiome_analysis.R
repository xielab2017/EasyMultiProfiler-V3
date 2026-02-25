#!/usr/bin/env Rscript
# EasyMultiProfiler - 微生物组分析脚本 (调用 EasyMultiProfiler R包)
# 16S/宏基因组数据分析

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

# 命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="OTU/ASV表路径"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="样本元数据路径"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="输出目录"),
  make_option(c("-p", "--params"), type="character", default="{}", help="分析参数"),
  make_option(c("-t", "--task-id"), type="character", default=NULL, help="任务ID")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# 解析参数
params <- fromJSON(args$params)

cat(sprintf("开始微生物组分析 - 任务ID: %s\n", args$task_id))

# 检查 EasyMultiProfiler
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
  cat(sprintf("⚠️ EasyMultiProfiler 加载失败: %s\n", e$message))
})

try {
  # 读取数据
  otu_table <- read.csv(args$input, row.names=1, check.names=FALSE)
  cat(sprintf("OTU表维度: %d 特征 x %d 样本\n", nrow(otu_table), ncol(otu_table)))
  
  # 准备分组信息
  group_info <- NULL
  if (!is.null(args$metadata) && file.exists(args$metadata)) {
    group_info <- read.csv(args$metadata, row.names=1)
  } else {
    n_samples <- ncol(otu_table)
    group_info <- data.frame(
      Group = c(rep("Group_A", floor(n_samples/2)), rep("Group_B", ceiling(n_samples/2))),
      row.names = colnames(otu_table)
    )
  }
  
  output_dir <- args$output
  
  # 使用 EasyMultiProfiler
  if (exists("EMP_easy_taxonomy_import") && exists("EMP_diff_analysis")) {
    
    cat("🔄 使用 EasyMultiProfiler 进行微生物组分析...\n")
    
    # 1. 导入数据
    cat("步骤1: 导入微生物组数据...\n")
    
    # 假设行名包含分类学信息，用分号分隔
    has_taxonomy <- any(grepl(";", rownames(otu_table)))
    
    if (has_taxonomy) {
      # 有分类学信息，使用 taxonomy_import
      MAE <- EMP_easy_taxonomy_import(
        data = otu_table,
        assay = "microbiome",
        assay_name = "counts",
        coldata = group_info,
        type = "tax",  # 16S 数据
        output = "MAE"
      )
    } else {
      # 无分类学信息，使用 normal_import
      MAE <- EMP_easy_normal_import(
        data = otu_table,
        assay = "microbiome",
        assay_name = "counts",
        coldata = group_info,
        output = "MAE"
      )
    }
    
    cat("✅ 微生物组数据导入完成\n")
    
    # 2. Alpha 多样性分析
    alpha_metric <- params$alpha$metric %||% "shannon"
    cat(sprintf("步骤2: Alpha多样性分析 (%s)...\n", alpha_metric))
    
    if (exists("EMP_alpha_analysis")) {
      tryCatch({
        MAE_alpha <- MAE |\u003e 
          EMP_alpha_analysis(
            experiment = "microbiome",
            method = alpha_metric,
            action = "add"
          )
        cat("✅ Alpha多样性分析完成\n")
        
        # 绘制Alpha多样性箱线图
        if (exists("EMP_boxplot_alpha")) {
          p_alpha <- MAE_alpha |\u003e EMP_boxplot_alpha()
          ggplot2::ggsave(file.path(output_dir, "alpha_diversity.png"), p_alpha,
                         width = 8, height = 6, dpi = 300)
        }
      }, error = function(e) {
        cat(sprintf("Alpha多样性分析警告: %s\n", e$message))
      })
    }
    
    # 3. Beta 多样性分析
    beta_method <- params$beta$method %||% "bray"
    cat(sprintf("步骤3: Beta多样性分析 (%s)...\n", beta_method))
    
    if (exists("EMP_dimension_analysis")) {
      tryCatch({
        MAE_beta <- MAE |\u003e
          EMP_dimension_analysis(
            experiment = "microbiome",
            method = beta_method,
            dimension = "pcoa",  # 或 NMDS
            action = "add"
          )
        cat("✅ Beta多样性分析完成\n")
        
        # 绘制PCoA图
        if (exists("EMP_scatterplot_reduce_dimension")) {
          p_pcoa <- MAE_beta |\u003e EMP_scatterplot_reduce_dimension()
          ggplot2::ggsave(file.path(output_dir, "beta_pcoa.png"), p_pcoa,
                         width = 8, height = 6, dpi = 300)
        }
      }, error = function(e) {
        cat(sprintf("Beta多样性分析警告: %s\n", e$message))
      })
    }
    
    # 4. 差异分析
    cat("步骤4: 差异分析...\n")
    diff_method <- params$diff$method %||% "wilcox"
    group_col <- colnames(group_info)[1]
    
    diff_result <- tryCatch({
      MAE |\u003e
        EMP_diff_analysis(
          experiment = "microbiome",
          method = diff_method,
          estimate_group = group_col,
          p.adjust = "fdr",
          action = "add"
        )
    }, error = function(e) {
      cat(sprintf("差异分析警告: %s\n", e$message))
      MAE
    })
    
    # 5. 保存结果
    cat("步骤5: 保存结果...\n")
    
    # 获取差异分析结果表
    de_table <- tryCatch({
      diff_result |\u003e .get.result.EMPT()
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(de_table) && nrow(de_table) > 0) {
      write.csv(de_table, file.path(output_dir, "differential_features.csv"), row.names=FALSE)
    }
    
    # 统计数据
    stats <- list(
      module = "microbiome",
      samples = ncol(otu_table),
      features = nrow(otu_table),
      alpha_metric = alpha_metric,
      beta_method = beta_method,
      diff_method = diff_method,
      task_id = args$task_id
    )
    write_json(stats, file.path(output_dir, "stats.json"))
    
    cat("✅ 微生物组分析完成！\n")
    
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

cat("微生物组分析脚本执行成功\n")
