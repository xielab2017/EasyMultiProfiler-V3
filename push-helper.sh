#!/bin/bash
# EasyMultiProfiler V3 完整推送方案
# 此脚本提供多种推送方式，选择最适合的一种

REPO_DIR="/Users/liweixie/.openclaw/workspaces/Tele_bot/Feishu_Bot_001/github-projects/EasyMultiProfiler-V3"

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║     EasyMultiProfiler V3 - GitHub推送方案                     ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

cd "$REPO_DIR"

echo "请选择推送方式:"
echo ""
echo "[方式1] 使用GitHub Token (推荐，最简单)"
echo "  - 访问 https://github.com/settings/tokens"
echo "  - 生成token后运行: ./push-with-token.sh YOUR_TOKEN"
echo ""
echo "[方式2] 使用SSH密钥 (如果已配置)"
echo "  运行: git push git@github.com:xielab2017/EasyMultiProfiler-V3.git main"
echo ""
echo "[方式3] 先创建仓库，再推送"
echo "  1. 访问 https://github.com/new"
echo "  2. 输入: EasyMultiProfiler-V3"
echo "  3. 不勾选README，点击Create"
echo "  4. 运行: git push -u origin main"
echo ""

# 检查是否有可用的认证方式
if git ls-remote https://github.com/xielab2017/EasyMultiProfiler-V3.git HEAD 2>/dev/null | grep -q "refs/heads"; then
    echo "✅ 检测到仓库已存在且有访问权限"
    echo ""
    read -p "是否直接推送? [Y/n]: " confirm
    if [[ $confirm =~ ^[Yy]$ ]] || [ -z "$confirm" ]; then
        git push -u origin main
        exit 0
    fi
else
    echo "⚠️  仓库不存在或无访问权限"
fi

echo ""
echo "推荐使用方式1 (Token):"
echo ""
echo "步骤:"
echo "1. 打开: https://github.com/settings/tokens/new"
echo "2. Note填写: EasyMultiProfiler-V3 Deploy"
echo "3. 勾选 'repo' 权限"
echo "4. 点击 Generate token"
echo "5. 复制token"
echo "6. 运行: ./push-with-token.sh ghp_xxxxxxxx"
echo ""
echo "或者直接运行:"
echo "  ./push-with-token.sh \$(pbpaste)  # Mac"
echo "  ./push-with-token.sh \$(xclip -o)  # Linux"
