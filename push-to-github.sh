#!/bin/bash
# EasyMultiProfiler V3 GitHub推送脚本
# 使用方法: ./push-to-github.sh YOUR_GITHUB_TOKEN

TOKEN=$1

if [ -z "$TOKEN" ]; then
    echo "使用方法: ./push-to-github.sh YOUR_GITHUB_TOKEN"
    echo ""
    echo "获取GitHub Token:"
    echo "1. 访问 https://github.com/settings/tokens"
    echo "2. 点击 'Generate new token (classic)'"
    echo "3. 勾选 'repo' 权限"
    echo "4. 生成token并复制"
    exit 1
fi

echo "🚀 推送到GitHub..."

# 配置远程使用token
git remote set-url origin https://${TOKEN}@github.com/xielab2017/EasyMultiProfiler-V3.git

# 推送
git push -u origin main

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ 推送成功！"
    echo "🌐 访问: https://github.com/xielab2017/EasyMultiProfiler-V3"
else
    echo ""
    echo "❌ 推送失败"
fi

# 重置远程URL（删除token）
git remote set-url origin https://github.com/xielab2017/EasyMultiProfiler-V3.git
