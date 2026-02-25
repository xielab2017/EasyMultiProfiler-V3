#!/bin/bash
# EasyMultiProfiler V3.0 一键安装脚本

set -e

echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║     EasyMultiProfiler V3.0 - 一键安装程序                    ║"
echo "║     R包 + Web + Docker 完整环境                              ║"
echo "╚═══════════════════════════════════════════════════════════════╝"
echo ""

INSTALL_DIR="${HOME}/EasyMultiProfiler-V3"

echo "📁 安装目录: $INSTALL_DIR"
echo ""

# 创建目录
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# 克隆仓库
echo "📥 下载 EasyMultiProfiler V3..."
if [ -d ".git" ]; then
    git pull origin main
else
    git clone https://github.com/xielab2017/EasyMultiProfiler-V3.git .
fi

echo ""
echo "请选择安装方式:"
echo "  1) Docker安装 (推荐，最简单)"
echo "  2) 本地安装 (需要R和Python环境)"
echo ""
read -p "请输入选项 [1-2]: " choice

case $choice in
    1)
        echo ""
        echo "🐳 Docker安装..."
        
        # 检查Docker是否安装
        if ! command -v docker >/dev/null 2>&1; then
            echo "❌ Docker未安装"
            echo ""
            echo "请安装Docker:"
            echo "  Mac: https://docs.docker.com/desktop/install/mac-install/"
            echo "  Linux: https://docs.docker.com/engine/install/"
            echo ""
            echo "安装后请启动Docker Desktop或docker服务"
            exit 1
        fi
        
        # 检查Docker守护进程是否运行
        if ! docker info >/dev/null 2>&1; then
            echo "❌ Docker守护进程未运行"
            echo ""
            echo "请启动Docker:"
            echo "  Mac: 打开Docker Desktop应用"
            echo "  Linux: sudo systemctl start docker"
            exit 1
        fi
        
        echo "✅ Docker已安装并运行"
        docker-compose up -d
        echo ""
        echo "✅ 安装完成！"
        echo "🌐 访问: http://localhost:8080"
        ;;
    
    2)
        echo ""
        echo "🔧 本地安装..."
        
        # 安装R包
        echo "📦 安装R包..."
        echo "   这可能需要几分钟，请耐心等待..."
        cd r-package
        
        # 检查并修复 DESCRIPTION 文件格式
        # 修复1: 删除注释行
        if [ -f "DESCRIPTION" ]; then
            if head -1 DESCRIPTION | grep -q "^#"; then
                echo "   修复 DESCRIPTION 格式 (移除注释行)..."
                sed -i '' '1d' DESCRIPTION
            fi
        fi
        
        # 修复2: 修复 Authors@R 字段的缩进（确保 ) 行有前导空格）
        if [ -f "DESCRIPTION" ]; then
            # 检查是否有单独的 ) 行（没有前导空格）
            if grep -n "^)$" DESCRIPTION > /dev/null 2>&1; then
                echo "   修复 DESCRIPTION 格式 (Authors@R 缩进)..."
                sed -i '' 's/^)$/    )/' DESCRIPTION
            fi
        fi
        
        # 创建临时R脚本文件
        R_SCRIPT_FILE=$(mktemp)
        cat > "$R_SCRIPT_FILE" << 'RSCRIPT_EOF'
options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(timeout = 300)

# 安装 devtools（如果不存在）
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org/")
}

# 配置 Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# 安装 Bioconductor 依赖
bioc_deps <- c("MOFA2", "ChIPseeker", "clusterProfiler")
for (dep in bioc_deps) {
    if (!requireNamespace(dep, quietly = TRUE)) {
        message(sprintf("从 Bioconductor 安装: %s", dep))
        tryCatch({
            BiocManager::install(dep, ask = FALSE)
        }, error = function(e) {
            message(sprintf("警告: %s 安装失败 - %s", dep, conditionMessage(e)))
        })
    }
}

# 安装 CRAN 核心依赖
cran_deps <- c("Seurat", "optparse", "dplyr", "ggplot2", "jsonlite")
for (dep in cran_deps) {
    if (!requireNamespace(dep, quietly = TRUE)) {
        message(sprintf("从 CRAN 安装: %s", dep))
        tryCatch({
            install.packages(dep, repos = "https://cloud.r-project.org/")
        }, error = function(e) {
            message(sprintf("警告: %s 安装失败 - %s", dep, conditionMessage(e)))
        })
    }
}

# 安装主包
devtools::install(".", dependencies = TRUE, upgrade = "never")
RSCRIPT_EOF
        
        # 执行R脚本
        Rscript "$R_SCRIPT_FILE" 2>&1 | tee r_install.log
        R_EXIT_CODE=${PIPESTATUS[0]}
        rm -f "$R_SCRIPT_FILE"
        
        if [ $R_EXIT_CODE -eq 0 ]; then
            echo "✅ R包安装成功"
        else
            echo "❌ R包安装失败，查看日志: r-package/r_install.log"
            echo ""
            echo "常见问题:"
            echo "  1. 缺少系统依赖 (如libcurl, libssl等)"
            echo "  2. 网络连接问题"
            echo "  3. 某些Bioconductor包需要单独安装"
            exit 1
        fi
        cd ..
        
        # 安装Python依赖
        echo "📦 安装Python依赖..."
        cd web/backend
        python3 -m venv venv
        source venv/bin/activate
        pip install -r requirements.txt
        cd ../..
        
        # 检查前端
        if [ -f "web/frontend/package.json" ]; then
            echo "📦 安装前端依赖..."
            cd web/frontend
            npm install
            npm run build
            # 复制构建文件到后端
            mkdir -p ../backend/static
            cp -r build/* ../backend/static/ 2>/dev/null || true
            cd ../..
        else
            echo "⚠️  前端代码未找到，跳过前端构建"
        fi
        
        echo ""
        echo "✅ 安装完成！"
        echo "🚀 启动命令: cd $INSTALL_DIR && ./scripts/start.sh"
        ;;
    *)
        echo "❌ 无效选项"
        exit 1
        ;;
esac
