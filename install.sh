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
        if command -v docker >/dev/null; then
            docker-compose up -d
            echo ""
            echo "✅ 安装完成！"
            echo "🌐 访问: http://localhost:8080"
        else
            echo "❌ 请先安装Docker: https://docs.docker.com/get-docker/"
            exit 1
        fi
        ;;
    2)
        echo ""
        echo "🔧 本地安装..."
        
        # 安装R包
        echo "📦 安装R包..."
        R -e "devtools::install('.', dependencies=TRUE)"
        
        # 安装Python依赖
        echo "📦 安装Python依赖..."
        cd web/backend
        python3 -m venv venv
        source venv/bin/activate
        pip install -r requirements.txt
        
        # 安装前端依赖
        echo "📦 安装前端依赖..."
        cd ../frontend
        npm install
        npm run build
        
        echo ""
        echo "✅ 安装完成！"
        echo "🚀 启动命令: cd $INSTALL_DIR && ./start.sh"
        ;;
    *)
        echo "❌ 无效选项"
        exit 1
        ;;
esac
