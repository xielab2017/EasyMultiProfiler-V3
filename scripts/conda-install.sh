#!/bin/bash
# EasyMultiProfiler V3.0 - Conda 简化安装脚本

set -e

echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║     EasyMultiProfiler V3.0 - Conda 简化安装                 ║"
echo "║     使用 Conda/Mamba 加速依赖安装                          ║"
echo "╚═══════════════════════════════════════════════════════════════╝"
echo ""

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 检查 conda/mamba
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo -e "${GREEN}✅ 检测到 Mamba (推荐)${NC}"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo -e "${GREEN}✅ 检测到 Conda${NC}"
else
    echo -e "${RED}❌ 未检测到 Conda 或 Mamba${NC}"
    echo ""
    echo "请安装 Miniforge 或 Miniconda:"
    echo "  - Miniforge (推荐): https://github.com/conda-forge/miniforge"
    echo "  - Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

INSTALL_DIR="${PWD}"
ENV_NAME="easymultiprofiler"

echo "📁 安装目录: $INSTALL_DIR"
echo "🐍 环境名称: $ENV_NAME"
echo ""

# 克隆仓库（如果不在仓库目录中）
if [ ! -f "environment.yml" ]; then
    echo "📥 下载 EasyMultiProfiler V3..."
    git clone https://github.com/xielab2017/EasyMultiProfiler-V3.git .
fi

# 创建/更新 conda 环境
echo ""
echo "🔧 创建 Conda 环境 (这可能需要 10-20 分钟)..."
echo "   使用 $CONDA_CMD 安装依赖..."

if $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
    echo -e "${YELLOW}⚠️  环境 '$ENV_NAME' 已存在，将更新...${NC}"
    $CONDA_CMD env update -f environment.yml -n $ENV_NAME
else
    $CONDA_CMD env create -f environment.yml -n $ENV_NAME
fi

# 激活环境
echo ""
echo "🚀 激活环境..."
eval "$($CONDA_CMD shell.bash hook)"
$CONDA_CMD activate $ENV_NAME

# 安装额外的 pip 依赖（如果存在）
if [ -f "web/backend/requirements.txt" ]; then
    echo ""
    echo "📦 安装 Python pip 依赖..."
    cd web/backend
    pip install -r requirements.txt
    cd ../..
fi

# 修复 DESCRIPTION 文件格式
echo ""
echo "🔧 修复 R 包描述文件格式..."
cd r-package

# 修复1: 删除注释行
if head -1 DESCRIPTION | grep -q "^#"; then
    echo "   修复 DESCRIPTION 格式 (移除注释行)..."
    sed -i.bak '1d' DESCRIPTION
fi

# 修复2: 修复 Authors@R 字段的缩进
if grep -n "^)$" DESCRIPTION > /dev/null 2>&1; then
    echo "   修复 DESCRIPTION 格式 (Authors@R 缩进)..."
    sed -i.bak 's/^)$/    )/' DESCRIPTION
fi

rm -f DESCRIPTION.bak
cd ..

# 构建并安装 R 包
echo ""
echo "📦 构建并安装 EasyMultiProfiler R 包..."
cd r-package

# 使用 Rscript 安装
Rscript -e "
if (!requireNamespace('devtools', quietly = TRUE)) {
    install.packages('devtools', repos = 'https://cloud.r-project.org/')
}
devtools::install('.', dependencies = FALSE, upgrade = 'never')
"

cd ..

# 构建前端
echo ""
echo "🎨 构建前端..."
if [ -f "web/frontend/package.json" ]; then
    cd web/frontend
    
    # 检查 node_modules 是否存在
    if [ ! -d "node_modules" ]; then
        echo "   安装 npm 依赖..."
        npm install
    fi
    
    echo "   构建生产版本..."
    npm run build
    
    # 复制到后端静态目录
    echo "   复制构建文件到后端..."
    mkdir -p ../backend/static
    cp -r build/* ../backend/static/ 2>/dev/null || true
    
    cd ../..
else
    echo -e "${YELLOW}⚠️  前端代码未找到${NC}"
fi

# 创建启动脚本
echo ""
echo "📝 创建启动脚本..."
cat > start.sh <> 'EOF'
#!/bin/bash
# EasyMultiProfiler V3.0 启动脚本

echo "🚀 启动 EasyMultiProfiler V3.0..."

# 激活 conda 环境
eval "$(conda shell.bash hook)"
conda activate easymultiprofiler

# 启动后端
cd web/backend
python app.py &
BACKEND_PID=$!

echo "✅ 后端服务已启动 (PID: $BACKEND_PID)"
echo "🌐 访问: http://localhost:5000"
echo ""
echo "按 Ctrl+C 停止服务"

# 等待中断
trap "kill $BACKEND_PID; exit" INT
wait $BACKEND_PID
EOF

chmod +x start.sh

echo ""
echo -e "${GREEN}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║                    ✅ 安装完成！                              ║${NC}"
echo -e "${GREEN}╚═══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "📋 总结:"
echo "  - Conda 环境: $ENV_NAME"
echo "  - R 包: 已安装"
echo "  - Python: 已配置"
echo "  - 前端: 已构建"
echo ""
echo "🚀 启动命令:"
echo "  conda activate $ENV_NAME"
echo "  ./start.sh"
echo ""
echo "🌐 访问地址: http://localhost:5000"
echo ""
