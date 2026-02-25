#!/usr/bin/env Rscript
# EasyMultiProfiler Web - 统一的R包调用包装器
# 所有分析模块都通过此脚本调用R包中的EMP_xxx_analysis函数

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

# 命令行参数
option_list <- list(
  make_option(c("-f", "--function"), type="character", default=NULL, 
              help="R包函数名, 如EMP_scrnaseq_analysis"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="输入数据文件路径"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="元数据文件路径"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="输出目录"),
  make_option(c("-p", "--params"), type="character", default="{}",
              help="分析参数(JSON格式)"),
  make_option(c("-t", "--task-id"), type="character", default=NULL,
              help="任务ID")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

cat(sprintf("╔══════════════════════════════════════════════════╗\n"))
cat(sprintf("║   EasyMultiProfiler Web - R Package Wrapper      ║\n"))
cat(sprintf("╚══════════════════════════════════════════════════╝\n\n"))

# 验证参数
if (is.null(args$function)) {
  cat("❌ 错误: 必须指定函数名 (--function)\n")
  quit(status=1)
}

if (is.null(args$input)) {
  cat("❌ 错误: 必须指定输入文件 (--input)\n")
  quit(status=1)
}

cat(sprintf("🎯 任务ID: %s\n", args$task_id %||% "N/A"))
cat(sprintf("🔧 函数: %s\n", args$function))
cat(sprintf("📁 输入: %s\n", args$input))
cat(sprintf("📂 输出: %s\n", args$output %||% "N/A"))
cat("\n")

# 加载EasyMultiProfiler
cat("📦 加载 EasyMultiProfiler...\n")
tryCatch({
  if (!requireNamespace("EasyMultiProfiler", quietly = TRUE)) {
    cat("⚠️  EasyMultiProfiler 未安装，尝试安装...\n")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools", repos = "https://cloud.r-project.org/")
    }
    devtools::install_github("xielab2017/EasyMultiProfiler", upgrade = "never")
  }
  library(EasyMultiProfiler)
  cat(sprintf("   版本: %s\n", as.character(packageVersion("EasyMultiProfiler"))))
}, error = function(e) {
  cat(sprintf("❌ 加载失败: %s\n", e$message))
  quit(status=1)
})

# 解析参数
params <- tryCatch({
  fromJSON(args$params)
}, error = function(e) {
  cat(sprintf("⚠️ 参数解析失败，使用默认参数: %s\n", e$message))
  list()
})

# 确保输出目录存在
if (!is.null(args$output)) {
  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
}

# 执行对应函数
cat(sprintf("\n🚀 执行 %s...\n\n", args$function))

tryCatch({
  
  result <- switch(args$function,
                   
    # ====== 单细胞RNA-seq ======
    "EMP_scrnaseq_analysis" = {
      cat("🧫 单细胞RNA-seq分析流程\n")
      
      # 读取数据
      counts <- read.csv(args$input, row.names=1, check.names=FALSE)
      cat(sprintf("   表达矩阵: %d 基因 x %d 细胞\n", nrow(counts), ncol(counts)))
      
      # 读取metadata（可选）
      metadata <- NULL
      if (!is.null(args$metadata) && file.exists(args$metadata)) {
        metadata <- read.csv(args$metadata, row.names=1)
        cat(sprintf("   元数据: %d 样本\n", nrow(metadata)))
      }
      
      # 调用R包函数
      EMP_scrnaseq_analysis(
        counts = counts,
        metadata = metadata,
        params = params,
        output_dir = args$output
      )
    },
    
    # ====== ChIP-seq ======
    "EMP_chipseq_analysis" = {
      cat("🧬 ChIP-seq分析流程\n")
      
      cat(sprintf("   Peak文件: %s\n", args$input))
      
      # 读取metadata（可选）
      metadata <- NULL
      if (!is.null(args$metadata) && file.exists(args$metadata)) {
        metadata <- read.csv(args$metadata, row.names=1)
      }
      
      EMP_chipseq_analysis(
        peak_file = args$input,
        metadata = metadata,
        params = params,
        output_dir = args$output
      )
    },
    
    # ====== CUT&Tag ======
    "EMP_cutntag_analysis" = {
      cat("✂️  CUT&Tag分析流程\n")
      
      metadata <- NULL
      if (!is.null(args$metadata) && file.exists(args$metadata)) {
        metadata <- read.csv(args$metadata, row.names=1)
      }
      
      EMP_cutntag_analysis(
        peak_file = args$input,
        metadata = metadata,
        params = params,
        output_dir = args$output
      )
    },
    
    # ====== CUT&RUN ======
    "EMP_cutnrun_analysis" = {
      cat("🔬 CUT&RUN分析流程\n")
      
      metadata <- NULL
      if (!is.null(args$metadata) && file.exists(args$metadata)) {
        metadata <- read.csv(args$metadata, row.names=1)
      }
      
      EMP_cutnrun_analysis(
        peak_file = args$input,
        metadata = metadata,
        params = params,
        output_dir = args$output
      )
    },
    
    # ====== RNA-seq (原有功能) ======
    "EMP_rnaseq_analysis" = {
      cat("📊 RNA-seq分析流程\n")
      
      counts <- read.csv(args$input, row.names=1, check.names=FALSE)
      
      metadata <- NULL
      if (!is.null(args$metadata) && file.exists(args$metadata)) {
        metadata <- read.csv(args$metadata, row.names=1)
      }
      
      # 使用现有的EMP分析流程
      # 这里可以调用EMP的核心差异分析函数
      EMP_diff_analysis(...)  # 根据实际参数调整
    },
    
    # ====== 微生物组 (原有功能) ======
    "EMP_microbiome_analysis" = {
      cat("🦠 微生物组分析流程\n")
      
      otu_table <- read.csv(args$input, row.names=1, check.names=FALSE)
      
      metadata <- NULL
      if (!is.null(args$metadata) && file.exists(args$metadata)) {
        metadata <- read.csv(args$metadata, row.names=1)
      }
      
      # 调用微生物组分析
      # 使用EMP_easy_taxonomy_import等函数
      EMP_microbiome_analysis(...)  # 根据实际参数调整
    },
    
    # ====== 多组学整合 ======
    "EMP_multiomics_integration" = {
      cat("🔗 多组学整合分析\n")
      
      # 多组学数据需要特殊处理
      # 输入目录包含多个组学数据文件
      data_files <- list.files(args$input, pattern = "\\.(csv|txt|tsv)$", full.names = TRUE)
      
      data_list <- lapply(data_files, function(f) {
        read.csv(f, row.names=1, check.names=FALSE)
      })
      names(data_list) <- tools::file_path_sans_ext(basename(data_files))
      
      cat(sprintf("   组学类型: %s\n", paste(names(data_list), collapse = " + ")))
      
      EMP_multiomics_integration(
        data_list = data_list,
        method = params$method %||% "MOFA2",
        params = params,
        output_dir = args$output
      )
    },
    
    # ====== 未知函数 ======
    {
      stop(sprintf("未知的函数: %s", args$function))
    }
  )
  
  cat("\n✅ 分析完成!\n")
  
  # 保存结果摘要
  if (!is.null(args$output)) {
    summary_file <- file.path(args$output, "EMP_summary.txt")
    sink(summary_file)
    print(result)
    sink()
    cat(sprintf("📄 结果摘要: %s\n", summary_file))
  }
  
}, error = function(e) {
  cat(sprintf("\n❌ 分析失败: %s\n", e$message))
  
  # 保存错误日志
  if (!is.null(args$output)) {
    error_file <- file.path(args$output, "error.log")
    writeLines(as.character(e), error_file)
  }
  
  quit(status=1)
})

cat("\n══════════════════════════════════════════════════\n")
cat("R包调用成功完成\n")

# 辅助函数
`%||%` <- function(x, y) if (is.null(x)) y else x
