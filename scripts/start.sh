#!/bin/bash
# EasyMultiProfiler V3.0 启动脚本

# 获取脚本所在目录的父目录（项目根目录）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_DIR="$(dirname "$SCRIPT_DIR")"

echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║     EasyMultiProfiler V3.0 - 启动程序                        ║"
echo "╚═══════════════════════════════════════════════════════════════╝"
echo ""
echo "📁 项目目录: $INSTALL_DIR"
echo ""

# 检查Docker模式
if [ -f "$INSTALL_DIR/docker-compose.yml" ] && command -v docker-compose >/dev/null; then
    echo "🐳 Docker模式启动..."
    cd "$INSTALL_DIR"
    docker-compose up -d
    echo ""
    echo "✅ 服务已启动"
    echo "🌐 访问: http://localhost:8080"
    echo ""
    echo "查看日志: docker-compose logs -f"
    echo "停止服务: docker-compose down"
else
    echo "🔧 本地模式启动..."
    
    # 启动后端
    echo "🚀 启动后端..."
    cd "$INSTALL_DIR/web/backend"
    
    if [ -f "venv/bin/activate" ]; then
        source venv/bin/activate
    fi
    
    # 检查静态文件
    if [ ! -d "static" ] || [ -z "$(ls -A static 2>/dev/null)" ]; then
        echo "⚠️  警告: 静态文件不存在，Web界面可能无法正常显示"
        echo "   请先运行: cd web/frontend && npm install && npm run build"
    fi
    
    python app.py &
    BACKEND_PID=$!
    echo $BACKEND_PID > "$INSTALL_DIR/backend.pid"
    
    echo ""
    echo "✅ 后端服务已启动"
    echo "🌐 访问地址: http://localhost:5000"
    echo ""
    echo "🛑 停止服务: kill $BACKEND_PID"
fi
