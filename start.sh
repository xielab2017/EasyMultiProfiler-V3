#!/bin/bash
# EasyMultiProfiler V3.0 启动脚本

INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║     EasyMultiProfiler V3.0 - 启动程序                        ║"
echo "╚═══════════════════════════════════════════════════════════════╝"
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
    source venv/bin/activate
    python app.py &
echo $! > "$INSTALL_DIR/backend.pid"
    
    # 启动前端
    echo "🌐 启动前端..."
    cd "$INSTALL_DIR/web/frontend"
    npm start &
echo $! > "$INSTALL_DIR/frontend.pid"
    
    echo ""
    echo "✅ 服务已启动"
    echo "🌐 Web界面: http://localhost:3000"
    echo "🔌 API地址: http://localhost:5000"
fi
